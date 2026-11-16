"""Motor de la simulación: corte transversal de un conducto mamario.

Geometría (de dentro a fuera): luz del conducto · epitelio luminal (donde nace
casi todo el cáncer de mama) · capa mioepitelial + membrana basal (la barrera
que separa in situ de invasivo) · estroma con fibroblastos y vasos sanguíneos.
"""

from __future__ import annotations

import random
from collections import deque
from dataclasses import replace

import numpy as np

from breast_cancer_abm.model.agents import (
    Agent,
    Analyze,
    Attack,
    Move,
    Observation,
    Perception,
    Repair,
    macrophage_policy,
    nk_policy,
    t_cell_policy,
)
from breast_cancer_abm.model.biology import (
    GENES,
    IMMUNE_GENES,
    SCENARIOS,
    STEPS_PER_DAY,
    TARGET_PREVALENCE,
    Parameters,
)
from breast_cancer_abm.model.policies import THERAPIES

LUMEN, EPITHELIUM, MYO, STROMA = range(4)
ZONE_NAME = ("lumen", "epithelium", "myoepithelium", "stroma")

LUMINAL, MYOEPITHELIAL, FIBROBLAST = range(3)
TISSUE_NAME = ("luminal", "myoepithelial", "fibroblast")

EMPTY, FIBRO, MYOEP, NORMAL, REPAIRED, INITIATED, HYPERPLASIA, IN_SITU, INVASIVE, NECROSIS, QUEEN = range(11)
STATE_NAME = {
    EMPTY: "vacío",
    FIBRO: "fibroblasto",
    MYOEP: "mioepitelial",
    NORMAL: "normal",
    REPAIRED: "reparada",
    INITIATED: "iniciada",
    HYPERPLASIA: "hiperplasia/atipia",
    IN_SITU: "carcinoma in situ",
    INVASIVE: "carcinoma invasivo",
    NECROSIS: "necrosis",
    QUEEN: "célula madre tumoral (reina)",
}
# Etiquetas en español para lo que ve el usuario (las claves internas van en inglés).
ZONE_LABEL = ("luz", "epitelio", "mioepitelio", "estroma")
TISSUE_LABEL = ("luminal", "mioepitelial", "fibroblasto")
TUMOR_STATES = (INITIATED, HYPERPLASIA, IN_SITU, INVASIVE, QUEEN)
O2_CONSUMPTION = np.array([0.0005, 0.002, 0.004, 0.004, 0.004, 0.005, 0.010, 0.016, 0.016, 0.0, 0.016])
CHEMOKINE_SOURCE = np.array([0, 0, 0, 0, 0, 0, 0.3, 1.0, 1.0, 1.0, 1.0])

NEIGHBORS = tuple((dy, dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx)
MAX_AGENTS = 500
NECROSIS_LIFETIME = 40


def _sign(v: float) -> int:
    return (v > 0) - (v < 0)


class Cell:
    __slots__ = (
        "id",
        "kind",
        "genome",
        "er",
        "mutations",
        "repaired",
        "necrosis",
        "last_division",
        "altered",
        "n_drivers",
        "prolif_signal",
        "p53_ok",
        "hrd",
        "cdh1_lost",
        "mhc",
        "pdl1",
        "her2",
        "anoikis_res",
        "invasive",
        "factor",
        "stem",
        "divisions",
        "guardian",
        "neo",
        "immortal",
        "telomerase_ok",
        "senescent",
        "bait",
        "aneuploid",
        "cd47",
        "caf",
        "corrector",
        "pdl1_until",
    )

    def __init__(
        self,
        id_: int,
        kind: int,
        genome: dict | None = None,
        er: bool = False,
        mutations: int = 0,
        stem: bool = False,
    ):
        self.id = id_
        self.stem = stem  # célula madre: se divide sin límite (la "reina")
        self.divisions = None  # divisiones que le quedan a una obrera tumoral (None = sin límite)
        self.guardian = None  # reina curada y reprogramada: "simple" | "doble"
        self.neo = 0  # neoantígenos que su MHC-I puede llegar a presentar
        self.telomerase_ok = True  # False si un inhibidor de telomerasa la ha alcanzado
        self.senescent = False  # jubilada: no se divide, pero sigue viva y segrega SASP
        self.bait = False  # cargada con p53 + 6-tio-dG: envenena a quien se la coma
        self.corrector = False  # lleva la maquinaria de corrección y la pasa a sus vecinas
        self.pdl1_until = -1  # PD-L1 ADAPTATIVO (inducido por IFN-γ) hasta este paso
        self.caf = False  # fibroblasto activado (myCAF): el que fabrica la malla
        self.cd47 = True  # "no me comas": lo lleva de serie TODA célula propia
        self.aneuploid = (
            False  # cromosomas descuadrados: casi exclusivo del tumor, pero INVISIBLE desde fuera
        )
        self.kind = kind
        self.genome = genome or {}
        self.er = er
        self.mutations = mutations  # pasajeras: alimentan los neoantígenos
        self.repaired = False
        self.necrosis = 0
        self.last_division = -999
        self.update_phenotype()

    def update_phenotype(self):
        """Traduce el genotipo (golpes por gen) a rasgos de comportamiento."""
        alt = set()
        for name, gene in GENES.items():
            h = self.genome.get(name, 0)
            if (gene.kind == "oncogene" and h >= 1) or h >= 2:
                alt.add(name)
        self.altered = frozenset(alt)
        self.n_drivers = len(alt - IMMUNE_GENES)
        self.prolif_signal = bool(alt & {"PIK3CA", "ERBB2", "MYC", "PTEN"})
        self.p53_ok = "TP53" not in alt
        self.hrd = bool(alt & {"BRCA1", "BRCA2"})
        self.cdh1_lost = "CDH1" in alt
        self.mhc = "HLA" not in alt
        self.pdl1 = "CD274" in alt
        self.her2 = "ERBB2" in alt
        # Sobrevivir sin anclaje (llenar la luz) exige señal de crecimiento + otro fallo.
        self.immortal = "TERT" in alt
        self.anoikis_res = self.prolif_signal and self.n_drivers >= 2
        self.invasive = self.anoikis_res and (self.cdh1_lost or self.n_drivers >= 4)
        f = 1.0
        if "PIK3CA" in alt:
            f *= 2.5
        if "ERBB2" in alt:
            f *= 3.0
        if "MYC" in alt:
            f *= 2.0
        if "PTEN" in alt:
            f *= 1.5
        if not self.p53_ok:
            f *= 1.3
        self.factor = min(f, 5.0)  # el ciclo celular tiene un mínimo físico (~1 día)


class Simulation:
    def __init__(self, scenario: str = "luminal", seed: int | None = None, params: Parameters | None = None):
        scen = SCENARIOS[scenario]
        self.scenario = scenario
        self.p = replace(params or Parameters(), estrogen=scen.estrogen)
        self.rng = random.Random(seed)
        self.t = 0
        self._next_id = 0
        self.cells: dict[tuple, Cell] = {}
        self.agents: list[Agent] = []
        self.events: deque[str] = deque(maxlen=12)
        self._milestones: set[str] = set()
        self.clones: dict[frozenset, int] = {}
        self.stats = dict(
            apoptosis_p53=0,
            anoikis=0,
            necrosis=0,
            immune_kills=0,
            therapy_kills=0,
            collateral=0,
            repairs=0,
            repairs_useless=0,
            membrane_breaches=0,
            new_vessels=0,
            dedifferentiations=0,
            guardians=0,
            phagocytosed=0,
            phagocytosed_healthy=0,
            dose=0,
            senescent_cells=0,
            healthy_stems_touched=0,
            propagations=0,
            baits=0,
            entosis=0,
            bait_kills=0,
            steps_blocked=0,
            cd47_blocked=0,
            cd47_healthy=0,
            macrophage_eaten=0,
            macrophage_eaten_healthy=0,
            cafs_activated=0,
            invasions_blocked=0,
            trojan_deliveries=0,
            sterilized=0,
            sterilized_healthy=0,
            guardians_reverted=0,
            propagated_corrections=0,
            fully_corrected=0,
            editing_errors=0,
            excised=0,
            healthy_excised=0,
            # Desglose de en qué acaba cada reparación génica (E0).
            # "reparaciones" cuenta solo las de _repair; el fármaco p53 y el
            # cebo llevan contador propio para no mezclar tres mecanismos.
            repairs_that_cure=0,
            repairs_p53_lethal=0,
            repairs_sterile=0,
            repairs_drug=0,
            repairs_bait=0,
            # ADC: molécula, no agente. Ver release_adc().
            adc_reached=0,
            adc_killed=0,
            adc_healthy=0,
            adc_quiescent_spared=0,
        )
        self.history: dict[str, list] = {
            k: [] for k in ("day", *TUMOR_STATES, REPAIRED, NECROSIS, "T", "NK", "AI")
        }
        self._event(f"Escenario '{scenario}': {scen.description}")
        self._build_tissue(scen)
        self._update_map()
        for _ in range(60):
            self._diffuse()
        self._record()

    # ------------------------------------------------------------------ montaje
    def _new_id(self) -> int:
        self._next_id += 1
        return self._next_id

    def _build_tissue(self, scen):
        n, r = self.p.grid_size, self.p.duct_radius
        self.center = (n - 1) / 2
        yy, xx = np.mgrid[0:n, 0:n]
        self.dist = np.hypot(yy - self.center, xx - self.center)
        self.zone = np.full((n, n), STROMA, np.uint8)
        self.zone[self.dist < r] = MYO
        self.zone[self.dist < r - 1.5] = EPITHELIUM
        self.zone[self.dist < r - 3.2] = LUMEN

        self.vessels = np.zeros((n, n), bool)
        periductal = (self.dist >= r + 2) & (self.dist < r + 3)
        far = self.dist >= r + 3
        rand = np.array([[self.rng.random() for _ in range(n)] for _ in range(n)])
        self.vessels[(periductal & (rand < 0.15)) | (far & (rand < 0.025))] = True

        for y in range(n):
            for x in range(n):
                z = self.zone[y, x]
                if z == EPITHELIUM:
                    self.cells[(y, x)] = Cell(
                        self._new_id(),
                        LUMINAL,
                        dict(scen.germline),
                        self.rng.random() < scen.er_fraction,
                        stem=self.rng.random() < self.p.normal_stem_fraction,
                    )
                elif z == MYO:
                    self.cells[(y, x)] = Cell(self._new_id(), MYOEPITHELIAL)
                elif z == STROMA and not self.vessels[y, x] and self.rng.random() < 0.08:
                    self.cells[(y, x)] = Cell(self._new_id(), FIBROBLAST)

        if scen.seed:
            pos = self.rng.choice([q for q, c in self.cells.items() if c.kind == LUMINAL and c.stem])
            c = self.cells[pos]
            c.genome.update(scen.seed)
            c.neo = sum(self.rng.random() < self.p.prob_neoantigen for _ in scen.seed)
            c.update_phenotype()
            self._record_clone(c, pos)
        # La madre SANA solo es madre en su nicho basal (ver _healthy_stem_off_niche): las que
        # el sorteo inicial dejó en filas suprabasales nacen ya diferenciadas.
        for pos, c in self.cells.items():
            if self._healthy_stem_off_niche(pos, c):
                c.stem = False

        self.tumor_vessels = np.zeros((n, n), bool)
        self.vessel_mode = "normal"
        self.O2 = np.full((n, n), 0.8)
        self.Q = np.zeros((n, n))  # quimiocinas (CXCL9/10...): atraen a las células inmunes
        self.M = np.zeros((n, n))  # densidad de malla: colágeno I + ácido hialurónico
        self.Ms = np.zeros((n, n))  # la misma, extendida: la malla no acaba en la celda del CAF
        self.P = np.zeros((n, n))  # presión intersticial: empuja hacia fuera de la masa
        self.pressure_factor = 1.0  # lo baja degradar la matriz o normalizar los vasos
        self.cd47_until = -1  # hasta qué paso sigue circulando el anti-CD47
        self.pd1_until = -1  # hasta qué paso sigue circulando el anti-PD-1 (ver block_pd1)

    # ------------------------------------------------------------------ bucle
    def step(self):
        self.t += 1
        order = list(self.cells.items())
        self.rng.shuffle(order)
        for pos, c in order:
            if self.cells.get(pos) is c:
                self._live(pos, c)
        self._update_map()
        self._diffuse()
        self._move_agents()
        self._recruit()
        self._record()

    @property
    def day(self) -> float:
        return self.t / STEPS_PER_DAY

    def _neighbors(self, pos):
        n = self.p.grid_size
        y, x = pos
        for dy, dx in NEIGHBORS:
            ny, nx = y + dy, x + dx
            if 0 <= ny < n and 0 <= nx < n:
                yield ny, nx

    def _live(self, pos, c: Cell):
        rng, p = self.rng, self.p
        if c.kind == FIBROBLAST:
            if not c.caf:
                # Dos vías de activación: el TGF-β del tumor vecino, y la propia rigidez
                # de la malla (YAP/TAZ). La segunda es la que la mantiene sola.
                tumor_near = any(
                    (v := self.cells.get(q)) is not None and v.n_drivers for q in self._neighbors(pos)
                )
                risk = (p.caf_activation_tgfb if tumor_near else 0.0) + p.caf_activation_stiffness * float(
                    self.M[pos]
                ) * p.lox
                if risk and rng.random() < risk:
                    c.caf = True
                    self.stats["cafs_activated"] += 1
                    self._milestone("caf", "Un fibroblasto se activa (myCAF) y empieza a tejer malla")
            elif self.M[pos] < 1.0:
                self.M[pos] = min(1.0, self.M[pos] + p.matrix_production)
            return
        if c.kind == MYOEPITHELIAL:
            if rng.random() < 0.003:  # la capa mioepitelial se autorrepara si puede
                gaps = [q for q in self._neighbors(pos) if self.zone[q] == MYO and q not in self.cells]
                if gaps:
                    self.cells[rng.choice(gaps)] = Cell(self._new_id(), MYOEPITHELIAL)
            return

        if c.necrosis:
            c.necrosis += 1
            if c.necrosis > NECROSIS_LIFETIME:  # los macrófagos retiran los restos
                del self.cells[pos]
            return
        o2 = self.O2[pos]
        # MUERTE POR O2 GRADUAL (variante B, elegida el 2026-10-28). Antes la única muerte por
        # O2 era O2 < o2_necrosis con p = 0,08/paso (vida 3,1 d), y entre o2_necrosis y
        # o2_hypoxia la tasa era ~0,0001/paso (vida ~1500-2500 d): las hipóxicas eran casi
        # inmortales y acumulaban el 28 % de las célula-pasos tumorales. Datos:
        #  - las hipóxicas humanas viven 4-10 d en xenoinjerto y 3-5 d en esferoide (Durand y
        #    Sham, IJROBP 1998;42:711): la hipoxia crónica SÍ mata, en días;
        #  - con hipoxia moderada (~0,5 % O2) no hay efecto apreciable sobre la supervivencia
        #    (Papandreou et al., Cancer Res 2005;65:3171): en el borde de la hipoxia ~0.
        # Por eso la probabilidad por paso crece LINEALMENTE de 0 en o2_hypoxia a 0,08 en
        # o2_necrosis (y 0,08 por debajo). La rampa es la interpolación más simple entre esos
        # dos extremos, no un ajuste a una curva medida. Vida en la banda: ~1/p (O2 = 0,095
        # -> p = 0,04 -> 6,25 d; ver biologia.STEPS_PER_DAY).
        if o2 < p.o2_necrosis:
            p_necrosis = 0.08
        elif o2 < p.o2_hypoxia:
            p_necrosis = 0.08 * (p.o2_hypoxia - o2) / (p.o2_hypoxia - p.o2_necrosis)
        else:
            p_necrosis = 0.0
        if p_necrosis and rng.random() < p_necrosis:
            c.necrosis = 1
            self.stats["necrosis"] += 1
            self._milestone(
                "necrosis", "Primera necrosis: el centro del tumor se queda sin oxígeno (comedonecrosis)"
            )
            return
        if c.n_drivers == 0 and rng.random() < p.prob_normal_death * (0.1 if c.stem else 1):
            del self.cells[pos]
            return
        if c.n_drivers and not c.stem:
            if c.divisions == 0 and rng.random() < p.prob_differentiated_death:
                del self.cells[pos]  # obrera agotada: envejece y muere sola
                return
            sasp = any((v := self.cells.get(q)) is not None and v.senescent for q in self._neighbors(pos))
            if rng.random() < (
                p.plasticity * (3 if o2 < p.o2_hypoxia else 1)
                + (p.prob_sasp_plasticity * (1 - p.drug_anti_sasp) if sasp else 0)
            ):
                # Plasticidad: una obrera se desdiferencia y se convierte en nueva reina.
                c.stem, c.divisions = True, None
                self.stats["dedifferentiations"] += 1
                self._milestone("plasticity", "Plasticidad: una obrera tumoral se convierte en nueva reina")

        if self._healthy_stem_off_niche(pos, c):
            c.stem = False  # ha perdido el nicho (o dejó de ser guardiana fuera de él)
        z = self.zone[pos]
        out_of_place = (z == LUMEN and not c.anoikis_res) or (z in (MYO, STROMA) and not c.invasive)
        if out_of_place and rng.random() < p.prob_anoikis:
            del self.cells[pos]
            self.stats["anoikis"] += 1
            return

        if c.invasive and z != STROMA and rng.random() < p.prob_invasion * 8:
            q = rng.choice(NEIGHBORS)
            q = (pos[0] + q[0], pos[1] + q[1])
            v = self.cells.get(q)
            if v is not None and v.kind == MYOEPITHELIAL:
                del self.cells[q]
                self.stats["membrane_breaches"] += 1
                self._milestone(
                    "breach", "Una célula invasiva rompe la capa mioepitelial (MMPs): empieza la invasión"
                )

        if c.corrector and p.correction_spread and p.correction_reach:
            self._correct_neighbors(pos, c)
        if c.guardian and p.guardian_reversion and rng.random() < p.guardian_reversion:
            # El estado impuesto se apaga: la célula sigue sana, pero deja de vigilar.
            c.guardian = None
            self.stats["guardians_reverted"] += 1
            self._milestone("reversion", "Una guardiana pierde el estado reprogramado y deja de vigilar")
            return
        if c.guardian and rng.random() < p.prob_phagocytosis * (
            p.internal_read_cost if c.guardian in self.INTERNAL_GATES else 1.0
        ):
            self._phagocytose(pos, c, c.guardian == "bait")
            if self.cells.get(pos) is not c:
                return
        if c.n_drivers and not c.bait and not c.senescent and rng.random() < p.prob_entosis:
            self._entosis(pos, c)
            if self.cells.get(pos) is not c:
                return
        if o2 < p.o2_hypoxia and c.n_drivers >= 2 and rng.random() < p.prob_angiogenesis:
            self._angiogenesis(pos)

        pdiv = p.prob_division * c.factor * ((1 + p.estrogen) if c.er else 1)
        if p.context_strength and c.n_drivers:
            # El tejido recupera la autoridad: la célula sigue mutada, pero obedece.
            pdiv *= 1 - p.context_strength * self._context(pos)
        if o2 < p.o2_hypoxia:
            pdiv *= 0.15
        if c.stem and c.n_drivers:
            pdiv *= p.stem_quiescence_factor  # las reinas se dividen poco: la quimio no las ve
        if rng.random() < pdiv:
            self._divide(pos, c)

    def _context(self, pos) -> float:
        """Cuánta autoridad tiene el tejido sobre la célula que ocupa `pos` (0 a 1).

        No es una señal que fabrique nadie: es la ARQUITECTURA. Una célula obedece
        cuando está en su sitio (epitelio), rodeada de vecinas sanas que la sujetan, y
        anclada a una membrana basal en condiciones. Es el contexto que Weaver y Bissell
        restauraron con el anticuerpo anti-β1-integrina para revertir el fenotipo
        maligno sin corregir un solo gen.

        Y por eso el tumor grande se defiende solo: a medida que crece destruye la
        arquitectura, y con ella la autoridad que podría frenarlo."""
        if self.zone[pos] not in (EPITHELIUM, MYO):
            return 0.0  # fuera de sitio no hay arquitectura que valga
        neighbors = [self.cells.get(q) for q in self._neighbors(pos)]
        present = [v for v in neighbors if v is not None and not v.necrosis]
        if not present:
            return 0.0  # sin vecinas no hay tejido: no hay quien mande
        healthy = sum(1 for v in present if v.kind != LUMINAL or v.n_drivers == 0)
        order = healthy / len(present)
        anchorage = 1.0 if any(self.zone[q] == MYO for q in self._neighbors(pos)) else 0.4
        if self.p.context_needs_matrix:
            # Una malla rígida y rota no sostiene la arquitectura: la destruye. Por eso
            # la reversión funcionó en matriz de laboratorio y no se ha logrado in vivo.
            anchorage *= max(0.0, 1.0 - float(self.Ms[pos]) * self.p.lox)
        return order * anchorage

    def _allowed(self, c: Cell, origin, dest) -> bool:
        zd = self.zone[dest]
        if self.vessels[dest]:
            return False
        if zd == EPITHELIUM:
            return True
        if zd == LUMEN:
            return c.anoikis_res
        if zd == MYO:
            return c.invasive and not self._blocked_by_mesh(dest)
        if not (c.invasive and self.zone[origin] in (MYO, STROMA)):
            return False
        return not self._blocked_by_mesh(dest)

    def _blocked_by_mesh(self, dest) -> bool:
        """La malla no solo estorba a la terapia: CONTIENE al tumor. Cuanto más densa y
        más entrecruzada (LOX), más cuesta atravesarla. Esta es la pieza sin la cual el
        modelo no podría reproducir por qué eliminar los CAFs empeoró los tumores."""
        brake = self.p.matrix_containment * float(self.Ms[dest]) * self.p.lox
        if brake and self.rng.random() < min(1.0, brake):
            self.stats["invasions_blocked"] += 1
            return True
        return False

    def _divide(self, pos, c: Cell):
        if c.divisions == 0:
            return  # obrera sin divisiones: senescencia
        gaps = [q for q in self._neighbors(pos) if q not in self.cells and self._allowed(c, pos, q)]
        if not gaps and c.prolif_signal:
            # Competición celular: el clon hiperproliferativo expulsa a una vecina sana.
            gaps = [
                q
                for q in self._neighbors(pos)
                if (v := self.cells.get(q)) is not None
                and v.kind == LUMINAL
                and v.n_drivers == 0
                and not v.necrosis
                and self._allowed(c, pos, q)
            ]
            if gaps and self.rng.random() < 0.3:
                del self.cells[gaps[0]]
                gaps = gaps[:1]
            else:
                gaps = []
        if not gaps:
            return  # inhibición por contacto: sin sitio no hay división
        q = self.rng.choice(gaps)
        daughter = Cell(self._new_id(), LUMINAL, dict(c.genome), c.er, c.mutations)
        daughter.repaired = c.repaired
        daughter.neo = c.neo
        # Mientras el anticuerpo circula, también alcanza a las células recién nacidas.
        daughter.cd47 = c.cd47 if self.t < self.cd47_until else True
        daughter.corrector = c.corrector
        # La aneuploidía se hereda, y cada mitosis sin control (p53 roto) puede generarla.
        pa = self.p.prob_aneuploidy if not c.p53_ok else self.p.prob_aneuploidy_p53ok
        daughter.aneuploid = c.aneuploid or self.rng.random() < pa
        if not c.aneuploid and daughter.aneuploid:
            self._milestone("aneuploidy", "Una mitosis sin control reparte mal los cromosomas: aneuploidía")
        immortal = (c.immortal or c.stem) and c.telomerase_ok
        drug_hit = (
            getattr(self, "drug_active", False)
            and self.p.drug_persistent
            and immortal
            and self.rng.random() < self.p.drug_persistent
        )
        if c.stem and c.n_drivers:
            # División asimétrica: casi siempre reina + obrera con divisiones contadas.
            daughter.stem = self.rng.random() < self.p.prob_symmetric_division
            daughter.telomerase_ok = c.telomerase_ok
            if daughter.stem and immortal:
                daughter.divisions = None
            else:
                daughter.divisions = (
                    self.p.worker_divisions if not daughter.immortal or not c.telomerase_ok else None
                )
            if not immortal and c.divisions is not None:
                c.divisions -= 1
        elif c.stem:
            # División orientada: la hija que sigue siendo madre es la que queda pegada a la
            # membrana basal; la que cae hacia la luz sale diferenciada. La guardiana es un
            # estado IMPUESTO y se hereda aunque la hija caiga fuera del nicho.
            daughter.stem = self.rng.random() < self.p.normal_stem_fraction and (
                bool(c.guardian and self.p.guardian_inheritance) or self._in_niche(q)
            )
            daughter.guardian = c.guardian if (daughter.stem and self.p.guardian_inheritance) else None
        else:
            daughter.telomerase_ok = c.telomerase_ok
            if immortal:
                daughter.divisions = c.divisions = None
            else:
                if c.divisions is None:
                    c.divisions = self.p.worker_divisions
                c.divisions -= 1
                daughter.divisions = c.divisions
        c.last_division = daughter.last_division = self.t
        if drug_hit:
            for cel in (c, daughter):
                if not cel.p53_ok and cel.genome.get("TP53") and self.rng.random() < self.p.drug_p53:
                    cel.genome["TP53"] = 0
                    cel.update_phenotype()
                if cel.p53_ok:
                    cel.senescent = True  # se marca; la apoptosis la resuelve el paso siguiente
                    cel.divisions = 0
                else:
                    cel.senescent, cel.divisions = True, 0
                    if self.p.drug_bait:
                        cel.bait = True
                        self.stats["baits"] += 1
        for cel, site in ((c, pos), (daughter, q)):
            hits = self._mutate(cel)
            damage = hits + (1 if cel.hrd and self.rng.random() < self.p.prob_hrd_damage else 0)
            if damage and cel.p53_ok and self.rng.random() < 1 - (1 - self.p.prob_apoptosis_p53) ** damage:
                self.stats["apoptosis_p53"] += 1  # p53 detecta el daño y la célula se suicida
                if cel is c:
                    del self.cells[pos]
                continue
            if cel is daughter:
                self.cells[q] = daughter
            if hits:
                self._record_clone(cel, site)

    def _mutate(self, c: Cell) -> int:
        mult = self.p.mutation_rate * (self.p.hrd_factor if c.hrd else 1) * (1 if c.p53_ok else 3)
        hits = 0
        for name, gene in GENES.items():
            if self.rng.random() < mult * gene.weight:
                h = c.genome.get(name, 0)
                if h < (1 if gene.kind == "oncogene" else 2):
                    c.genome[name] = h + 1
                    hits += 1
                    c.neo += self.rng.random() < self.p.prob_neoantigen
        lam = 0.3 * (10 if c.hrd else 1)
        fresh = int(lam) + (self.rng.random() < lam - int(lam))
        c.mutations += fresh
        c.neo += sum(self.rng.random() < self.p.prob_neoantigen / 15 for _ in range(fresh))
        if hits:
            c.update_phenotype()
            if c.n_drivers and not c.stem and c.divisions is None and not c.immortal:
                c.divisions = self.p.worker_divisions  # mutación en célula de vida corta
        return hits

    # Puertas lógicas de la guardiana. Cada una es un IDENTIFICADOR distinto del tumor;
    # ver INFORME_IDENTIFICADORES.md: ninguna señal tiene a la vez cobertura y exclusividad,
    # así que la especificidad tiene que venir de exigir DOS condiciones a la vez.
    INTERNAL_GATES = ("tert", "tert+niche")

    def _gate(self, mode: str, v: Cell, o) -> bool:
        # Señales de SUPERFICIE (lo único que un receptor tipo synNotch puede leer de verdad).
        tumor_signal = o.neoantigens > 0 or not o.mhc_i or o.senescent
        activity_signal = o.proliferating or o.hypoxic or o.stress_ligands
        if mode == "simple":
            return activity_signal  # una sola señal: se equivoca
        if mode in ("double", "bait"):
            return tumor_signal and activity_signal  # neoantígeno/MHC/senescencia + actividad
        if mode == "niche":
            # Puerta Y posicional: marcador de célula madre Y fuera de la capa basal.
            # La madre SANA lleva el mismo marcador, pero vive pegada a la membrana basal.
            return o.stem_marker and o.off_niche
        if mode in self.INTERNAL_GATES:
            # Puerta Y interna: telomerasa activa Y cromosomas descuadrados. La madre sana
            # tiene telomerasa pero NO es aneuploide; el resto del tejido no tiene telomerasa.
            # Exige leer dentro del núcleo: más lenta (internal_read_cost) y con falsos
            # positivos (internal_read_error).
            real = (v.immortal or v.stem) and v.telomerase_ok and v.aneuploid
            if mode == "tert+niche":
                real = real or (o.stem_marker and o.off_niche)
            return real or self.rng.random() < self.p.internal_read_error
        return tumor_signal and activity_signal

    def _correct_neighbors(self, pos, c: Cell):
        """Efecto espectador: la célula corregida pasa la corrección a sus vecinas.

        Corrige UN gen por intento, como haría un vector que transfiere una copia sana
        por las uniones comunicantes. Cuando una vecina queda limpia, se vuelve corrector
        a su vez: de ahí el crecimiento exponencial, y de ahí también el peligro."""
        p, rng = self.p, self.rng
        neighbors = [q for q in self._neighbors(pos)][: p.correction_reach]
        for q in neighbors:
            v = self.cells.get(q)
            if v is None or v.kind != LUMINAL or v.necrosis or v.corrector:
                continue
            if rng.random() >= p.correction_spread:
                continue
            broken = [g for g in v.genome if v.genome[g]]
            if not broken:
                v.corrector = True  # ya estaba limpia: ahora también propaga
                continue
            gene = min(broken, key=lambda g: ("TP53" != g, g))  # p53 primero
            if rng.random() < p.repair_efficacy:
                v.genome[gene] -= 1
                self.stats["propagated_corrections"] += 1
                if rng.random() < p.editing_error:
                    # Error fuera de diana: la edición rompe otro gen. El precio de que
                    # algo se copie solo es que los errores también se copian.
                    other = rng.choice([g for g in GENES if g != gene])
                    v.genome[other] = v.genome.get(other, 0) + 1
                    self.stats["editing_errors"] += 1
                v.update_phenotype()
                if not any(v.genome.values()):
                    v.corrector = True
                    v.repaired = True
                    self.stats["fully_corrected"] += 1
                    self._milestone(
                        "propagation", "Una célula corregida transmite la corrección a una vecina"
                    )

    def _phagocytose(self, pos, c: Cell, ag_baits: bool = False):
        """La guardiana se come a una vecina que pasa su puerta lógica y ocupa el hueco
        con una hija sana (competición celular / extrusión epitelial). Si ag_baits, a las
        senescentes las deja cargadas como cebo en lugar de retirarlas."""
        for q in self._neighbors(pos):
            v = self.cells.get(q)
            if v is None or v.kind != LUMINAL or v.necrosis or v.guardian:
                continue
            o = self._observe(v, q, 0, 0)
            target = self._gate(c.guardian, v, o)
            if not target:
                continue
            if ag_baits and v.senescent and not v.bait:
                v.bait = True  # trampa cargada: el SASP deja de ser abono y pasa a ser anzuelo
                self.stats["baits"] += 1
                self._milestone("bait", "Una célula senescente queda CARGADA como cebo (p53 + 6-tio-dG)")
                return
            del self.cells[q]
            self.stats["phagocytosed"] += 1
            if v.n_drivers == 0:
                self.stats["phagocytosed_healthy"] += 1
                self.stats["collateral"] += 1
                self._milestone("guardian_error", "ERROR: una guardiana se come a una célula SANA")
            else:
                self._milestone("guardian", "Una guardiana se come a una obrera corrupta y ocupa su sitio")
            if self.rng.random() < 0.5:
                self._divide(pos, c)
            return

    def resect(self, coverage: float = 0.95, margin: int = 2) -> dict:
        """Cirugía: retira un BLOQUE de tejido que cubre la mayor parte de la masa.

        No es una terapia elegante y no pretende serlo. Es la única forma conocida de
        convertir un problema de VOLUMEN en uno de superficie, que es lo que este
        simulador ha acabado señalando como el muro de los tumores grandes.

        Se modela como es: un corte contiguo con margen, no una retirada selectiva.
        Se lleva por delante tejido sano (eso es el precio) y deja fuera lo que ya había
        escapado del bloque: eso es la ENFERMEDAD RESIDUAL, y es el problema real de
        después, no el tumor que se quitó.

        CUIDADO AL INTERPRETARLA (medido 2026-10-29, ver T18 en validation.py). Esto NO
        es una tumorectomía: el tumor acaba rodeando la luz (ocupa ~90 % del anillo en
        ~5 focos), así que el centroide cae en el centro del conducto y el disco engloba
        el conducto ENTERO. Con cobertura 0,95 y margen 2 se extirpa el 93-96 % del
        conducto, sea a 400 células o el día 120. En términos clínicos es la resección
        completa del único órgano del campo, no una cirugía conservadora, y no se puede
        comparar con las series de conservadora. Bajar la cobertura NO lo arregla: un
        margen quirúrgico real (1-2 mm = 67-133 celdas) no cabe en un campo de 110.
        """
        tumor = [q for q, c in self.cells.items() if c.kind == LUMINAL and c.n_drivers and not c.necrosis]
        if not tumor:
            return dict(excised=0, healthy=0, residual=0)
        cy = sum(q[0] for q in tumor) / len(tumor)
        cx = sum(q[1] for q in tumor) / len(tumor)
        dist = sorted(((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5 for q in tumor)
        radius = dist[min(len(dist) - 1, int(len(dist) * coverage))] + margin
        excised = healthy = 0
        for q, c in list(self.cells.items()):
            if ((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5 > radius:
                continue
            if c.kind == LUMINAL and c.n_drivers:
                excised += 1
            else:
                healthy += 1
            del self.cells[q]
        self._update_map()
        residual = sum(1 for c in self.cells.values() if c.kind == LUMINAL and c.n_drivers and not c.necrosis)
        self.stats["excised"] = excised
        self.stats["healthy_excised"] = healthy
        self._event(
            f"CIRUGÍA: se extirpan {excised} células tumorales y {healthy} sanas; "
            f"queda un residuo de {residual}"
        )
        return dict(excised=excised, healthy=healthy, residual=residual)

    def ablate(
        self,
        coverage: float = 0.95,
        margin: int = 2,
        rim: int = 2,
        rim_survival: float = 0.0,
        guided: bool = False,
        visible_threshold: int = 1,
    ) -> dict:
        """Ablación (crio o radioterapia ablativa): destruir el tumor SIN sacarlo.

        Inspirada en el ensayo FIRST (crioablación, 48 pacientes, 2026: ablación completa
        89,6 %, CDIS residual en el BORDE en el 12,5 %) y en la radioterapia ablativa de
        MD Anderson (ESTRO 2026). Misma geometría que resect() para que sean comparables,
        con dos diferencias que son las reales:
          · lo destruido NO se retira: queda como necrosis hasta que los macrófagos lo
            limpian (NECROSIS_LIFETIME), y esa necrosis ya cuenta como inflamación en
            _recruit() — el efecto "vacuna in situ" sale del modelo, no se le añade;
          · el filo de la bola de hielo no es letal del todo: en la corona exterior de
            `borde` celdas cada célula sobrevive con probabilidad `rim_survival`.
            SUPUESTO: no hay dato que traduzca isotermas a celdas de 15 µm; se barre.
        Igual que resect(), a esta escala el volumen engloba el conducto casi entero:
        destruye tejido sano y lo dice.

        `guiada=True` (2026-11-05): ablación guiada por imagen. Apunta SOLO al carcinoma
        (in situ, invasivo, reina), que es lo que se ve; las células iniciadas o
        hiperplásicas son invisibles y no arrastran el disco. Destruye todo lo que esté a
        ≤ `margen` celdas de algún foco, con varias sondas si es multifocal.
        SUPUESTO OPTIMISTA, declarado: se asume que la imagen localiza focos de muy pocas
        células; la resonancia real no ve por debajo de ~mm."""
        if guided:
            return self._guided_ablation(margin, rim, rim_survival, visible_threshold)
        tumor = [q for q, c in self.cells.items() if c.kind == LUMINAL and c.n_drivers and not c.necrosis]
        if not tumor:
            return dict(destroyed=0, healthy=0, residual=0)
        cy = sum(q[0] for q in tumor) / len(tumor)
        cx = sum(q[1] for q in tumor) / len(tumor)
        self._lesion_center = (cy, cx)
        dist = sorted(((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5 for q in tumor)
        radius = dist[min(len(dist) - 1, int(len(dist) * coverage))] + margin
        self._lesion_radius = radius
        destroyed = healthy = 0
        for q, c in self.cells.items():
            d = ((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5
            if d > radius or c.necrosis:
                continue
            if d > radius - rim and self.rng.random() < rim_survival:
                continue
            if c.kind == LUMINAL and c.n_drivers:
                destroyed += 1
            else:
                healthy += 1
            c.necrosis = 1
        self._update_map()
        residual = sum(1 for c in self.cells.values() if c.kind == LUMINAL and c.n_drivers and not c.necrosis)
        self.stats["ablated"] = destroyed
        self.stats["healthy_ablated"] = healthy
        self._event(
            f"ABLACIÓN: se destruyen in situ {destroyed} células tumorales y "
            f"{healthy} sanas; queda un residuo de {residual}"
        )
        return dict(destroyed=destroyed, healthy=healthy, residual=residual)

    def carcinoma_foci(self) -> list[list]:
        """Focos de carcinoma = componentes conexas (8 vecinas) de células de carcinoma."""
        carcinoma = (IN_SITU, INVASIVE, QUEEN)
        sites = {q for q, c in self.cells.items() if c.kind == LUMINAL and self.state(c, q) in carcinoma}
        foci, seen = [], set()
        for q in sites:
            if q in seen:
                continue
            focus, stack_ = [], [q]
            seen.add(q)
            while stack_:
                y, x = stack_.pop()
                focus.append((y, x))
                for dy, dx in NEIGHBORS:
                    v = (y + dy, x + dx)
                    if v in sites and v not in seen:
                        seen.add(v)
                        stack_.append(v)
            foci.append(focus)
        return foci

    def _guided_ablation(
        self, margin: int, rim: int, rim_survival: float, visible_threshold: int = 1
    ) -> dict:
        """`visible_threshold`: tamaño mínimo (células) de un foco para que la imagen lo vea.
        A la escala del modelo (15 µm/celda) una resonancia real (≥2 mm) no vería NINGÚN
        foco: el valor realista es 'infinito'. Se barre para ver cuánto depende de ello."""
        carcinoma = (IN_SITU, INVASIVE, QUEEN)
        foci = self.carcinoma_foci()
        visible = [f for f in foci if len(f) >= visible_threshold]
        targets = [q for f in visible for q in f]
        if not targets:
            self._event(
                f"IMAGEN: {len(foci)} focos de carcinoma, ninguno visible "
                f"(umbral {visible_threshold} células): no hay nada que ablacionar"
            )
            return dict(destroyed=0, healthy=0, residual=sum(map(len, foci)), foci_seen=0)
        self._lesion_center = (
            sum(q[0] for q in targets) / len(targets),
            sum(q[1] for q in targets) / len(targets),
        )
        self._lesion_radius = max(
            3,
            max(
                ((q[0] - self._lesion_center[0]) ** 2 + (q[1] - self._lesion_center[1]) ** 2) ** 0.5
                for q in targets
            ),
        )
        dmin: dict = {}
        for y, x in targets:
            for dy in range(-margin, margin + 1):
                for dx in range(-margin, margin + 1):
                    d = (dy * dy + dx * dx) ** 0.5
                    if d <= margin:
                        q = (y + dy, x + dx)
                        if d < dmin.get(q, 1e9):
                            dmin[q] = d
        destroyed = healthy = 0
        for q, d in dmin.items():
            c = self.cells.get(q)
            if c is None or c.necrosis:
                continue
            if d > margin - rim and self.rng.random() < rim_survival:
                continue
            if c.kind == LUMINAL and c.n_drivers:
                destroyed += 1
            else:
                healthy += 1
            c.necrosis = 1
        self._update_map()
        residual = sum(
            1 for q, c in self.cells.items() if c.kind == LUMINAL and self.state(c, q) in carcinoma
        )
        self.stats["ablated"] = self.stats.get("ablated", 0) + destroyed
        self.stats["healthy_ablated"] = self.stats.get("healthy_ablated", 0) + healthy
        self._event(
            f"ABLACIÓN GUIADA: {len(visible)}/{len(foci)} focos visibles "
            f"({len(targets)} células de carcinoma); se destruyen {destroyed} células "
            f"con driver y {healthy} sanas; carcinoma residual {residual}"
        )
        return dict(destroyed=destroyed, healthy=healthy, residual=residual, foci_seen=len(visible))

    def biopsy(self, cores: int = 6, cells_per_core: int = 12) -> dict:
        """Biopsia por vacío del lecho tumoral (la llave para NO operar, Kuerer 2022/2025).

        Toma `cilindros` muestras de `cells_per_core` celdas contiguas alrededor de
        donde estaba la lesión (el clip). Positiva si sale UNA sola célula tumoral viva.
        No se le programa ninguna tasa de falsos negativos: sale del muestreo, y se compara
        después con la de los ensayos (RESPONDER 17,8 %)."""
        if getattr(self, "_lesion_center", None) is None:
            tumor = [q for q, c in self.cells.items() if c.kind == LUMINAL and c.n_drivers and not c.necrosis]
            if not tumor:
                return dict(positive=False, sampled=0, real=0)
            self._lesion_center = (
                sum(q[0] for q in tumor) / len(tumor),
                sum(q[1] for q in tumor) / len(tumor),
            )
            dist = sorted(
                ((q[0] - self._lesion_center[0]) ** 2 + (q[1] - self._lesion_center[1]) ** 2) ** 0.5
                for q in tumor
            )
            self._lesion_radius = dist[int(len(dist) * 0.95)] if dist else 5
        cy, cx = self._lesion_center
        r = max(self._lesion_radius, 3)
        seen_set = set()
        for _ in range(cores):
            ang = self.rng.uniform(0, 2 * np.pi)
            rad = r * self.rng.random() ** 0.5
            oy, ox = cy + rad * np.sin(ang), cx + rad * np.cos(ang)
            for k in range(cells_per_core):  # el cilindro: una fila de celdas
                seen_set.add(
                    (
                        int(round(oy + (k - cells_per_core / 2) * np.sin(ang + 1.57))),
                        int(round(ox + (k - cells_per_core / 2) * np.cos(ang + 1.57))),
                    )
                )
        # Lo que ve un patólogo: CARCINOMA (in situ, invasivo o célula madre tumoral). Una
        # célula "iniciada" (un driver, morfología normal) es invisible al microscopio.
        # Corregido 2026-11-05: la primera versión contaba cualquier célula con driver
        # como positiva y como "tumor real", e inflaba los falsos negativos.
        carcinoma = (IN_SITU, INVASIVE, QUEEN)
        positive = any(
            (c := self.cells.get(q)) is not None and c.kind == LUMINAL and self.state(c, q) in carcinoma
            for q in seen_set
        )
        real = sum(1 for q, c in self.cells.items() if c.kind == LUMINAL and self.state(c, q) in carcinoma)
        self._event(
            f"BIOPSIA: {len(seen_set)} celdas muestreadas → "
            f"{'POSITIVA' if positive else 'negative'} (carcinoma real: {real})"
        )
        return dict(positive=positive, sampled=len(seen_set), real=real)

    def ctdna(self, fragments_per_cell: float, necrosis_weight: float = 3.0) -> dict:
        """Biopsia líquida: ADN tumoral circulante en una extracción de sangre.

        Cada célula de carcinoma viva vierte ADN al morir y renovarse; las necróticas
        tumorales lo vierten de golpe (`necrosis_weight`). Los fragmentos mutantes que caen
        en el tubo siguen una Poisson; positivo si aparece ≥1 (ensayo tumor-informado).
        `fragments_per_cell` es la única incógnita y NO se inventa: se calibra fuera
        contra la sensibilidad clínica antes del tratamiento (ver residual_test.py).
        La célula 'iniciada' también vierte ADN, pero el ensayo busca las mutaciones del
        TUMOR ORIGINAL: cuentan las de carcinoma, el mismo criterio que biopsy()."""
        carcinoma = (IN_SITU, INVASIVE, QUEEN)
        alive = sum(1 for q, c in self.cells.items() if c.kind == LUMINAL and self.state(c, q) in carcinoma)
        dead = sum(1 for c in self.cells.values() if c.kind == LUMINAL and c.n_drivers and c.necrosis)
        lam = fragments_per_cell * (alive + necrosis_weight * dead)
        n = int(np.random.default_rng(self.rng.getrandbits(32)).poisson(lam))
        return dict(positive=n > 0, fragments=n, expected=lam, alive=alive)

    def restore_architecture(self, strength: float = 0.8):
        """Anti-β1-integrina: devuelve al tejido la autoridad sobre sus células.

        No corrige genes, no mata, no entra en ninguna célula. Es una señal que baña el
        tejido entero, así que su coste NO crece con el tamaño del tumor: es la única
        intervención de este simulador que trata igual de barato 10 células que 10.000.

        Es reversible: `fuerza=0` la retira y las células vuelven a mandar sobre el tejido."""
        self.p.context_strength = strength
        self._event(
            f"ARQUITECTURA: el tejido recupera autoridad sobre sus células "
            f"(fuerza {strength:.0%}); las mutaciones siguen ahí"
        )

    def block_pd1(self, days: float | None = None):
        """Anti-PD-1 / anti-PD-L1 (dostarlimab, pembrolizumab...): quita el freno.

        Intervención MÍNIMA y aislada (añadida 2026-10-28 para el test T15 de
        validation.py). Lo único que hace es anular, mientras el anticuerpo circula, el
        efecto que ya existía en el modelo de que una célula con PD-L1 (CD274 alterado)
        agote al linfocito T que la reconoce por TCR (o CAR). No toca el reclutamiento,
        ni los neoantígenos, ni el MHC-I, ni las NK, ni ningún parámetro: si el linfocito
        no veía a la célula, sigue sin verla. `dias=None` = bloqueo hasta el final.
        Por defecto (sin llamarla) el modelo se comporta exactamente igual que antes."""
        self.pd1_until = 10**12 if days is None else self.t + int(days * STEPS_PER_DAY)
        self._event("TERAPIA: anti-PD-1 — el PD-L1 deja de agotar a los linfocitos T")

    def block_cd47(self, efficacy: float | None = None, targeted: bool = False):
        """Anti-CD47 (como el magrolimab): le quita a las células el cartel de 'no me comas'.

        No es selectivo y no puede serlo: CD47 lo lleva todo el cuerpo. Lo que impide que
        los macrófagos se coman el tejido sano NO es este fármaco, es la segunda señal:
        la célula sana no expone calreticulina. En la clínica el precio de esto es la
        anemia, porque los glóbulos rojos sí llevan mucho CD47 y son frágiles."""
        ef = self.p.anti_cd47_efficacy if efficacy is None else efficacy
        self.cd47_until = self.t + self.p.anti_cd47_duration * STEPS_PER_DAY
        touched = healthy = 0
        for pos, c in self.cells.items():
            if c.kind != LUMINAL or self.rng.random() >= ef:
                continue
            if c.n_drivers == 0 and not self.p.cd47_reaches_healthy:
                continue
            if targeted:
                # Segunda diana del biespecífico: algo que la célula sana no muestra.
                o = self._observe(c, pos, 0, 0)
                if not (o.her2_high or o.stress_ligands or o.neoantigens > 0 or not o.mhc_i):
                    continue
            c.cd47 = False
            touched += 1
            healthy += c.n_drivers == 0
        self.stats["cd47_blocked"] = touched
        self.stats["cd47_healthy"] = healthy
        self._event(
            f"ANTI-CD47: se retira el 'no me comas' a {touched} células "
            f"({healthy} de ellas SANAS: el fármaco no distingue)"
        )

    def damage_telomeres(self):
        """Rompe el telómero AL INSTANTE en las células con telomerasa activa.

        Si p53 funciona, la célula se autodestruye. Si no, solo se jubila y se queda
        gritando (SASP). Las cualidades de Parameters.farmaco_* modifican el resultado:
        restaurar p53, dejar la carga dentro (cebo), silenciar el SASP o persistir."""
        p = self.p
        self.drug_active = True
        dead = retired = healthy = 0
        for pos, c in list(self.cells.items()):
            telomerase = c.immortal or c.stem
            if c.kind != LUMINAL or not telomerase or not c.telomerase_ok:
                continue
            if self.rng.random() >= self.p.efficacy_6tiodg:
                continue
            if c.n_drivers == 0 and c.stem:
                healthy += 1
            if not c.p53_ok and c.genome.get("TP53") and self.rng.random() < p.drug_p53:
                c.genome["TP53"] = 0  # el propio fármaco restaura p53 (vector tipo Gendicine)
                c.update_phenotype()
                self.stats["repairs_drug"] += 1
            if c.p53_ok:
                del self.cells[pos]  # daño irreparable + p53 = apoptosis
                self.stats["apoptosis_p53"] += 1
                dead += 1
            else:
                c.senescent = True  # sin p53 solo se jubila: hay que retirarla aparte
                c.divisions = 0
                if p.drug_bait:
                    c.bait = True  # la carga se queda dentro: la jubilada es una mina
                    self.stats["baits"] += 1
                retired += 1
        self.stats["healthy_stems_touched"] += healthy
        self._event(
            f"6-TIO-dG: {dead} células se autodestruyen (p53 intacto), "
            f"{retired} solo se jubilan (p53 roto), {healthy} madres sanas afectadas"
        )

    def inhibit_telomerase(self):
        """La 'castración': bloquea la telomerasa. Toda célula alcanzada recupera un contador
        de divisiones finito y se jubila al agotarlo. Daño colateral: las células madre SANAS
        del tejido también usan telomerasa."""
        touched = healthy = 0
        for c in self.cells.values():
            if c.kind != LUMINAL or self.rng.random() >= self.p.antitelomerase_efficacy:
                continue
            c.telomerase_ok = False
            touched += 1
            if c.divisions is None:
                c.divisions = self.p.worker_divisions
            if c.n_drivers == 0 and c.stem:
                healthy += 1
        self.stats["healthy_stems_touched"] += healthy
        self._event(
            f"INHIBIDOR DE TELOMERASA: {touched} células pierden el rebobinado "
            f"({healthy} madres sanas afectadas)"
        )

    def _entosis(self, pos, c: Cell):
        """La célula tumoral se traga a una vecina (entosis, descrita en cáncer de mama).
        Si era un cebo, recibe p53 reparado + telómero roto: el banquete la mata."""
        baits = [q for q in self._neighbors(pos) if (v := self.cells.get(q)) is not None and v.bait]
        if not baits:
            return
        q = self.rng.choice(baits)
        del self.cells[q]  # la digiere en su lisosoma
        self.stats["entosis"] += 1
        if self.rng.random() >= self.p.bait_efficacy:
            return
        if c.genome.get("TP53"):
            c.genome["TP53"] = 0  # carga 1: p53 restaurado
            c.update_phenotype()
            self.stats["repairs_bait"] += 1
        if c.p53_ok and (c.immortal or c.stem):
            del self.cells[pos]  # carga 2: telómero roto + p53 = autodestrucción
            self.stats["bait_kills"] += 1
            self.stats["apoptosis_p53"] += 1
            self._milestone("bait_kills_cell", "Una célula tumoral se come un cebo y se autodestruye")
            return
        c.senescent, c.divisions = True, 0
        if self.rng.random() < self.p.prob_bait_chain:
            c.bait = True  # la comedora se vuelve cebo: la trampa se propaga
            self.stats["baits"] += 1

    def release_adc(
        self,
        target: str = "trop2",
        affinity: str = "high",
        prevalence: float | None = None,
        dose: float = 0.35,
        bystander_radius: int = 8,
        cycle_dependent: bool = True,
    ) -> dict:
        """Anticuerpo conjugado (ADC). NO es un agente que camina: es una molécula.

        Tres diferencias con las terapias celulares del modelo, y las tres importan:

        1. PENETRA POCO Y DESDE EL VASO. Medido en humanos con variantes de afinidad
           de un anti-HER2: alta afinidad se queda a ~20 µm del vaso, baja afinidad
           llega a ~85 µm. Es la "barrera del sitio de unión": agarrarse fuerte a la
           primera célula impide llegar al fondo. A 15 µm/celda → 1,3 y 5,7 celdas.
        2. EL VENENO SE DERRAMA 50-200 µm, así que mata también a la vecina que NO
           expresa la diana. Es lo que resolvió la heterogeneidad antigénica, y la
           razón de que la FDA apruebe estos fármacos SIN exigir expresión de diana.
        3. EL PAYLOAD DEPENDE DE LA DIVISIÓN. Auristatinas, maitansinoides e
           inhibidores de topoisomerasa I necesitan una célula en ciclo. Nunca se han
           probado contra células dormidas. Por eso `cycle_dependent`.

        `prevalencia` es la fracción de células tumorales que expresan la diana EN
        CRIBADO PROSPECTIVO, que no es la de la literatura: con folato-α se pasó de un
        80 % publicado a un 10 % real medido en 80 pacientes consecutivas, y la tasa de
        respuesta fue 0. Por defecto se usa la del cribado, no la publicada.
        """
        p, rng = self.p, self.rng
        reach = {"high": 1, "mean": 3, "low": 6}.get(affinity, 3)  # celdas desde el vaso
        prev = TARGET_PREVALENCE.get(target, 0.56) if prevalence is None else prevalence

        # Distancia al vaso más cercano (la molécula sale de ahí y difunde)
        dist = np.full((p.grid_size, p.grid_size), 1e9)
        front = [tuple(q) for q in np.argwhere(self.vessels)]
        for q in front:
            dist[q] = 0
        d = 0
        while front and d <= reach + bystander_radius:
            d += 1
            new_front = []
            for y, x in front:
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    q = (y + dy, x + dx)
                    if 0 <= q[0] < p.grid_size and 0 <= q[1] < p.grid_size and dist[q] > d:
                        dist[q] = d
                        new_front.append(q)
            front = new_front

        # 1) Células ALCANZADAS: dentro del alcance Y que expresan la diana
        reached = []
        for pos, c in self.cells.items():
            if c.kind != LUMINAL or not c.n_drivers or c.necrosis:
                continue
            if dist[pos] > reach:
                continue
            if rng.random() < prev:  # ¿expresa la diana?
                reached.append(pos)

        # 2) El veneno se derrama: zona letal alrededor de cada célula alcanzada
        zone = set()
        for pos in reached:
            for dy in range(-bystander_radius, bystander_radius + 1):
                for dx in range(-bystander_radius, bystander_radius + 1):
                    if dy * dy + dx * dx <= bystander_radius * bystander_radius:
                        zone.add((pos[0] + dy, pos[1] + dx))

        dead = healthy = quiescent_spared = 0
        for q in list(zone):
            c = self.cells.get(q)
            if c is None or c.kind != LUMINAL or c.necrosis:
                continue
            if cycle_dependent and (c.divisions == 0 or c.senescent or self.t - c.last_division >= 40):
                quiescent_spared += c.n_drivers > 0
                continue  # el payload no la toca: no se divide
            if rng.random() >= dose:
                continue
            if c.n_drivers == 0:
                healthy += 1
                self.stats["collateral"] += 1
            del self.cells[q]
            dead += 1
        self._update_map()
        self.stats["adc_reached"] += len(reached)
        self.stats["adc_killed"] += dead
        self.stats["adc_healthy"] += healthy
        self.stats["adc_quiescent_spared"] += quiescent_spared
        self._event(
            f"ADC anti-{target} (afinidad {affinity}): alcanza {len(reached)} "
            f"células con diana, mata {dead} ({healthy} sanas); "
            f"{quiescent_spared} tumorales quiescentes SOBREVIVEN"
        )
        return dict(reached=len(reached), dead=dead, healthy=healthy, quiescent=quiescent_spared)

    def treat_matrix(self, mode: str, intensity: float = 0.6):
        """Las cinco maneras de atacar la malla. NO son equivalentes, y ahí está el asunto:
        en la clínica todo lo que la DESTRUYE ha fracasado y lo que la ABLANDA funciona.

        'hialuronidasa'  degrada el ácido hialurónico (PEGPH20). Baja la presión de golpe,
                         pero los CAFs siguen ahí y vuelven a tejer. Fase III falló (2019).
        'anti_lox'       impide el entrecruzamiento (simtuzumab). No retira material:
                         lo destensa. Baja presión Y contención a la vez.
        'losartan'       reduce la PRODUCCIÓN de colágeno y HA. Barato, positivo en
                         páncreas con quimio (Murphy 2019).
        'anti_tgfb'      corta la orden de activación: los fibroblastos no se vuelven CAF.
        'matar_caf'      elimina al fabricante. Es lo que empeoró los tumores en ratones.
        """
        if mode == "hyaluronidase":
            self.M *= 1 - intensity
            message = f"se degrada el {intensity:.0%} de la malla; los CAFs siguen ahí"
        elif mode == "anti_lox":
            self.p.lox *= 1 - intensity
            message = f"entrecruzamiento al {self.p.lox:.0%}: la malla se destensa sin retirarse"
        elif mode == "losartan":
            self.p.matrix_production *= 1 - intensity
            message = f"producción de malla al {1 - intensity:.0%}"
        elif mode == "anti_tgfb":
            self.p.caf_activation_tgfb *= 1 - intensity
            self.p.caf_activation_stiffness *= 1 - intensity
            message = "se corta la activación de nuevos CAFs"
        elif mode == "kill_caf":
            return self.degrade_matrix(intensity)
        else:
            raise ValueError(f"modo de matriz desconocido: {mode}")
        self._event(f"MALLA ({mode}): {message}")

    def degrade_matrix(self, fraction: float = 0.25):
        """Enzima tipo PEGPH20/hialuronidasa: deshace la matriz densa, baja la presión y
        abre huecos por los que puede circular una célula terapéutica. No mata nada:
        retira fibroblastos y descomprime el estroma."""
        removed = 0
        for pos, c in list(self.cells.items()):
            if c.kind == FIBROBLAST and self.rng.random() < fraction:
                del self.cells[pos]
                removed += 1
        self.pressure_factor *= 1 - self.p.pressure_drop_matrix * fraction
        self._event(
            f"MATRIZ: se degradan {removed} fibroblastos; la presión baja al "
            f"{self.pressure_factor:.0%} y se abren huecos de paso"
        )

    def treat_vessels(self, mode: str):
        """normal · antiangiogenico (bloquea VEGF: retira vasos tumorales y no deja crecer más)
        · normalizar (los vasos tumorales riegan bien) · asfixiar (disruptor vascular cerca del tumor)
        · reconducir (desvía oxígeno del tejido sano lejano a la zona tumoral: menos hipoxia)."""
        self.vessel_mode = mode
        if mode == "normalize":
            # Normalizar la vasculatura cierra las fugas: baja la presión Y mejora el riego.
            self.pressure_factor *= 1 - self.p.pressure_drop_normalize
        if mode == "antiangiogenic":
            self.vessels &= ~self.tumor_vessels
            self.tumor_vessels[:] = False
        self._event(f"VASOS: modo {mode}")

    def _angiogenesis(self, pos):
        if self.vessel_mode == "antiangiogenic":
            return
        y, x = pos
        ny, nx = y + self.rng.randint(-5, 5), x + self.rng.randint(-5, 5)
        n = self.p.grid_size
        if 0 <= ny < n and 0 <= nx < n and self.zone[ny, nx] == STROMA and not self.vessels[ny, nx]:
            v = self.cells.get((ny, nx))
            if v is None or v.kind == FIBROBLAST:
                self.cells.pop((ny, nx), None)
                self.vessels[ny, nx] = True
                self.tumor_vessels[ny, nx] = True
                self.stats["new_vessels"] += 1
                self._milestone("vegf", "El tumor hipóxico segrega VEGF: crece un vaso nuevo (angiogénesis)")

    # ------------------------------------------------------------------ campos
    def state(self, c: Cell, pos) -> int:
        if c.kind == FIBROBLAST:
            return FIBRO
        if c.kind == MYOEPITHELIAL:
            return MYOEP
        if c.necrosis:
            return NECROSIS
        if c.n_drivers == 0:
            return REPAIRED if c.repaired else NORMAL
        if c.stem:
            return QUEEN
        if self.zone[pos] in (MYO, STROMA):
            return INVASIVE
        if c.anoikis_res:
            return IN_SITU
        if c.prolif_signal:
            return HYPERPLASIA
        return INITIATED

    def _update_map(self):
        self.code = np.zeros(self.zone.shape, np.uint8)
        for pos, c in self.cells.items():
            self.code[pos] = self.state(c, pos)
        self.counts = np.bincount(self.code.ravel(), minlength=11)

    def _diffuse(self, iterations: int = 6):
        k = O2_CONSUMPTION[self.code]
        source = CHEMOKINE_SOURCE[self.code]
        perfusion = self.vessels.astype(float)
        if self.vessel_mode not in ("normalize", "reroute"):
            perfusion[self.tumor_vessels] = self.p.o2_tumor_vessel
        targeted = None
        if self.vessel_mode in ("starve", "reroute"):
            near = np.isin(self.code, TUMOR_STATES)
            for _ in range(4):
                near = (
                    near
                    | np.roll(near, 1, 0)
                    | np.roll(near, -1, 0)
                    | np.roll(near, 1, 1)
                    | np.roll(near, -1, 1)
                )
            if self.vessel_mode == "starve":
                perfusion[near & self.vessels] = self.p.o2_starved_vessel
            else:
                # Reconducir: el flujo se desvía del tejido sano lejano hacia la zona tumoral.
                perfusion[~near & self.vessels] = 0.85
                targeted = near.astype(float)
        self.M *= 1 - self.p.matrix_turnover
        # La malla que teje un CAF no se queda en su celda: impregna el entorno, incluida
        # la membrana basal que el tumor tiene que atravesar para invadir.
        self.Ms = self.M.copy()
        for axis_mid in (0, 1):
            self.Ms = np.maximum(self.Ms, 0.7 * np.roll(self.M, 1, axis_mid))
            self.Ms = np.maximum(self.Ms, 0.7 * np.roll(self.M, -1, axis_mid))
        self._diffuse_pressure(iterations)
        for cfield, sink, supply in ((self.O2, k, None), (self.Q, 0.08, source)):
            for _ in range(iterations):
                pad = np.pad(cfield, 1, mode="edge")
                lap = pad[:-2, 1:-1] + pad[2:, 1:-1] + pad[1:-1, :-2] + pad[1:-1, 2:] - 4 * cfield
                cfield += 0.2 * lap - sink * cfield
                if supply is not None:
                    cfield += 0.05 * supply
                else:
                    cfield[self.vessels] = perfusion[self.vessels]
                    if targeted is not None:
                        cfield += 0.03 * targeted * (1 - cfield)
                np.clip(cfield, 0.0, 1.0, out=cfield)

    def _diffuse_pressure(self, iterations: int = 6):
        """Campo de presión intersticial. Lo generan las células tumorales (vasos porosos
        que gotean plasma) y los fibroblastos (matriz densa que impide el drenaje); fuera
        de la masa los linfáticos sanos lo drenan deprisa. El resultado es una meseta de
        presión alta dentro y una caída brusca en el borde: justo el perfil que se mide
        en tumores reales, y la razón de que el borde se trate solo y el centro no."""
        # El ácido hialurónico de la malla es un gel que retiene agua: la presión sale
        # de ahí, no de que haya fibroblastos. Por eso se puede ABLANDAR sin destruir.
        source = np.isin(self.code, TUMOR_STATES).astype(float)
        source += 1.2 * self.Ms * self.p.lox
        source *= self.p.interstitial_pressure * self.pressure_factor
        drainage = np.where(source > 0, 0.02, 0.25)  # linfáticos: colapsados dentro, bien fuera
        for _ in range(iterations):
            pad = np.pad(self.P, 1, mode="constant")
            lap = pad[:-2, 1:-1] + pad[2:, 1:-1] + pad[1:-1, :-2] + pad[1:-1, 2:] - 4 * self.P
            self.P += 0.2 * lap + 0.08 * source - drainage * self.P
            np.clip(self.P, 0.0, 1.0, out=self.P)

    def _blocked_by_pressure(self, origin, dest, trojan: float = 0.0) -> bool:
        """¿La corriente de salida impide este paso? Solo frena a quien intenta ENTRAR
        (moverse hacia donde la presión sube); salir es cuesta abajo y es gratis.

        A un vehículo troyano casi no le afecta: no nada contra la corriente, lo lleva
        dentro una célula que el tumor ha reclutado y a la que le han abierto el paso."""
        if not self.p.interstitial_pressure or trojan >= 1.0:
            return False
        rise = float(self.P[dest] - self.P[origin])
        if rise <= 0:
            return False
        brake = self.p.pressure_effect * rise * (1 - trojan)
        return self.rng.random() < min(1.0, brake)

    # ------------------------------------------------------------------ agentes
    def _create_agent(self, kind: str, pos) -> Agent:
        if kind == "T":
            return Agent(
                "T", "linfocito T CD8+", t_cell_policy, frozenset({"TCR"}), pos=pos, expansion=0.2, life=160
            )
        if kind == "MAC":
            # El tumor polariza a la mayoría hacia M2: presentes, pero no fagocitan.
            m2 = self.rng.random() < self.p.m2_macrophage_fraction
            return Agent(
                "MAC",
                "macrófago M2 (no come)" if m2 else "macrófago",
                macrophage_policy,
                frozenset() if m2 else frozenset({"FAGOCITO"}),
                pos=pos,
                life=240,
            )
        return Agent("NK", "célula NK", nk_policy, frozenset({"NK"}), pos=pos, life=80)

    def _recruit(self):
        natural = sum(1 for a in self.agents if a.kind in ("T", "NK") and not a.clonal)
        inflammation = (self.counts[IN_SITU] + self.counts[INVASIVE] + self.counts[NECROSIS]) / 60
        if natural < self.p.max_natural_immune and self.rng.random() < self.p.immune_recruitment * (
            1 + inflammation
        ):
            kind = "T" if self.rng.random() < 0.6 else "NK"
            self.agents.append(self._create_agent(kind, self._random_vessel()))
        # Los macrófagos están siempre ahí: no hay que infundirlos, hay que desbloquearlos.
        n_macrophages = sum(1 for a in self.agents if a.kind == "MAC")
        if n_macrophages < self.p.max_macrophages and self.rng.random() < self.p.macrophage_recruitment * (
            1 + inflammation
        ):
            self.agents.append(self._create_agent("MAC", self._random_vessel()))

    def _random_vessel(self, near: bool = False):
        ys, xs = np.nonzero(self.vessels)
        if near:
            d = self.dist[ys, xs]
            sel = np.argsort(d)[:40]
            ys, xs = ys[sel], xs[sel]
        i = self.rng.randrange(len(ys))
        return int(ys[i]), int(xs[i])

    def _delivery_positions(self, n: int, route: str):
        """De dónde salen las células terapéuticas.

        'vena'          por la sangre: aparecen en un vaso y tienen que abrirse paso (lo de siempre)
        'intratumoral'  aguja guiada por ecografía: aparecen DENTRO de la masa, repartidas
        'aguja'         una sola punción: todas en un punto del tumor (lo realista con una aguja)
        """
        tumor = [q for q, c in self.cells.items() if c.kind == LUMINAL and c.n_drivers and not c.necrosis]
        if route == "vein" or not tumor:
            return [self._random_vessel() for _ in range(n)]
        if route == "needle":
            center = self.rng.choice(tumor)
            near = sorted(tumor, key=lambda q: sum((a - b) ** 2 for a, b in zip(q, center)))[:n]
            return [self.rng.choice(near) for _ in range(n)]
        return [self.rng.choice(tumor) for _ in range(n)]

    def release_therapy(self, name: str, amount: int | None = None, route: str = "vein"):
        """Admite combinaciones: "guardiana+memoria". `via`: vena | intratumoral | aguja."""
        if "+" in name:
            for part in name.split("+"):
                self.release_therapy(part, amount, route)
            return
        cfg = THERAPIES[name]
        n = amount or cfg["n"]
        self.stats["dose"] += n
        for pos in self._delivery_positions(n, route):
            ag = cfg["create"]()
            ag.pos = pos
            self.agents.append(ag)
        self._event(f"TERAPIA: se infunden {amount or cfg['n']} × {name}")

    def _observe(self, c: Cell, pos, dy, dx) -> Observation:
        hypoxic = self.O2[pos] < self.p.o2_hypoxia
        neo = c.neo if (c.mhc and c.n_drivers > 0) else 0
        if c.stem:
            neo //= 3  # las madres tumorales presentan menos antígeno (MHC-I bajo)
        return Observation(
            dy,
            dx,
            c.id,
            TISSUE_NAME[c.kind],
            c.mhc,
            neo,
            self._pdl1(c),
            c.her2,
            c.senescent or c.n_drivers >= 2 or (hypoxic and c.n_drivers > 0),
            self.t - c.last_division < 8,
            hypoxic,
            c.stem,
            c.senescent,
            c.cd47,
            self._exposes_eat_me(c, hypoxic),
            c.stem and not self._in_niche(pos),
        )

    @staticmethod
    def _exposes_eat_me(c: Cell, hypoxic: bool) -> bool:
        """Calreticulina en la superficie: la señal de 'cómeme'. La expone la célula
        jubilada (senescente), la dañada y la muy estresada. Una célula sana que sigue
        su ciclo normal NO la expone: esa es la mitad de la puerta que la protege
        aunque le hayan bloqueado el CD47."""
        if c.n_drivers == 0 and not c.senescent:
            return False
        return c.senescent or c.divisions == 0 or c.bait or (c.n_drivers >= 2 and hypoxic)

    def _healthy_stem_off_niche(self, pos, c: Cell) -> bool:
        """La identidad de célula madre SANA depende del nicho: contacto con el mioepitelio y
        la lámina basal (Notch/Wnt). Fuera de él se diferencia y pierde CD44+/CD24-/ALDH1.
        La reina TUMORAL no depende del nicho (independencia de nicho = rasgo del tumor), y la
        guardiana tampoco: su estado es impuesto. Corrige un defecto (2026-10-28): el modelo
        dejaba madres sanas suprabasales indefinidamente con el marcador de madre, que es
        exactamente lo que la puerta "nicho" busca -> falso positivo (T13 con el flag
        adaptativo: semilla 1, día 232; gates.py sin tratar: 6/301 sanas tocadas el día 120).
        No da a la guardiana ninguna percepción nueva ni toca la puerta.
        MEDIDO: daño a sanas en la terapia de T13, semillas 1-24: antes 3/24 semillas sin flag
        (7, 16, 18) y 2/24 con flag (1, 16); después 0/24 en ambos. Control del tumor igual
        (semillas 9-24: 12/16 → 11/16 sin flag, 12/16 → 12/16 con flag; a igual tamaño el
        día 120 se controla lo mismo). gates.py día 120: nicho toca 0 sanas (antes 6)."""
        return (
            c.stem and c.n_drivers == 0 and not c.guardian and c.kind == LUMINAL and not self._in_niche(pos)
        )

    def _in_niche(self, pos) -> bool:
        """El nicho de la célula madre mamaria es la capa basal: epitelio EN CONTACTO con
        el mioepitelio. Una reina tumoral no respeta esa arquitectura y acaba fuera."""
        if self.zone[pos] != EPITHELIUM:
            return False
        return any(self.zone[q] == MYO for q in self._neighbors(pos))

    def _perceive(self, ag: Agent) -> Perception:
        y, x = ag.pos
        n = self.p.grid_size
        neighbors = []
        for dy in range(-2, 3):
            for dx in range(-2, 3):
                c = self.cells.get((y + dy, x + dx))
                if c is not None and not c.necrosis:
                    neighbors.append(self._observe(c, (y + dy, x + dx), dy, dx))
        best, chemo = self.Q[y, x] + 1e-4, None
        for dy, dx in NEIGHBORS:
            ny, nx = y + dy, x + dx
            if 0 <= ny < n and 0 <= nx < n and self.Q[ny, nx] > best:
                best, chemo = self.Q[ny, nx], (dy, dx)
        return Perception(
            neighbors,
            chemo,
            (_sign(self.center - y), _sign(self.center - x)),
            ZONE_NAME[self.zone[y, x]],
            float(self.O2[y, x]),
            self.day,
            self.rng,
        )

    def _move_agents(self):
        n = self.p.grid_size
        living = []
        adaptive = self.p.adaptive_immunity
        if adaptive:
            self._n_clonal = sum(1 for a in self.agents if a.clonal and a.life > 1)
        for ag in self.agents:
            ag.life -= 1
            if ag.life <= 0 or ag.exhaustion >= 1:
                continue
            living.append(ag)
            if adaptive and ag.kind == "T" and ag.t_state:
                self._t_cycle(ag, living)
            if ag.delivers and self._deliver_trojan(ag, living):
                continue
            act = ag.policy(ag, self._perceive(ag))
            if isinstance(act, Move):
                dest = (
                    min(n - 1, max(0, ag.pos[0] + _sign(act.dy))),
                    min(n - 1, max(0, ag.pos[1] + _sign(act.dx))),
                )
                if not self._blocked_by_pressure(ag.pos, dest, ag.trojan):
                    ag.pos = dest
                    if ag.oxygen:
                        # La bombona va soltando O2 por donde pasa.
                        self.O2[ag.pos] = min(1.0, self.O2[ag.pos] + ag.oxygen)
                else:
                    self.stats["steps_blocked"] += 1
                continue
            if not isinstance(act, (Attack, Analyze, Repair)) or max(abs(act.dy), abs(act.dx)) > 1:
                continue
            q = (ag.pos[0] + act.dy, ag.pos[1] + act.dx)
            c = self.cells.get(q)
            if c is None or c.necrosis:
                continue
            if isinstance(act, Attack):
                self._attack(ag, c, q, living)
            elif isinstance(act, Analyze):
                # Secuenciar: devuelve TODOS los golpes, incluidos alelos heterocigotos.
                ag.memory.setdefault("genomes", {})[c.id] = {g: h for g, h in c.genome.items() if h}
            else:
                self._repair(ag, c, q, act.gene)
        self.agents = living

    def _deliver_trojan(self, ag: Agent, living: list) -> bool:
        """El vehículo ha llegado: suelta la carga y desaparece. La terapia real nace
        AQUÍ DENTRO, que es la gracia — ya no tiene que atravesar nada."""
        surrounded = sum(
            1 for q in self._neighbors(ag.pos) if (c := self.cells.get(q)) is not None and c.n_drivers
        )
        if surrounded < 5:
            return False
        cfg = THERAPIES[ag.delivers]
        for _ in range(3):
            if len(living) >= MAX_AGENTS:
                break
            new = cfg["create"]()
            new.pos = ag.pos
            living.append(new)
        living.remove(ag)
        self.stats["trojan_deliveries"] += 1
        self._milestone("trojan", f"Un vehículo llega al interior del tumor y suelta {ag.delivers}")
        return True

    def _recognizes(self, ag: Agent, c: Cell, q) -> str | None:
        r = ag.receptors
        if "NONSPECIFIC" in r:
            return "nonspecific"
        if c.kind != LUMINAL:
            return None
        obs = self._observe(c, q, 0, 0)
        if "FAGOCITO" in r:
            # Puerta Y natural: hace falta que falte el freno Y que esté el acelerador.
            return "phagocyte" if (obs.eat_me and not obs.cd47) else None
        if "HER2" in r and c.her2:
            return "car"
        if "TCR" in r and obs.neoantigens > 0:
            return "tcr"
        if "NK" in r and (not c.mhc or obs.stress_ligands):
            return "nk"
        return None

    def _attack(self, ag: Agent, c: Cell, q, living: list):
        route = self._recognizes(ag, c, q)
        if route is None:
            return
        adaptive = self.p.adaptive_immunity
        t_natural = adaptive and ag.kind == "T"
        if t_natural and route == "tcr":
            self._stimulate_t(ag)  # cada encuentro con su antígeno cuenta (Wherry 2011)
        pdl1 = self._pdl1(c)
        if route in ("tcr", "car") and pdl1 and not ag.pdl1_resistant and self.t >= self.pd1_until:
            self.stats["pdl1_brakes"] = self.stats.get("pdl1_brakes", 0) + 1
            if c.pdl1:
                self.stats["pdl1_brakes_genetic"] = self.stats.get("pdl1_brakes_genetic", 0) + 1
            if not t_natural:
                ag.exhaustion += 0.2  # PD-L1 frena al linfocito: se agota
            # Con inmunidad adaptativa el T natural NO desaparece por el freno: queda
            # inhibido este turno. Lo que lo agota es la estimulación crónica (_stimulate_t);
            # sin PD-1 se agota igual o más (Odorizzi, J Exp Med 2015;212:1125).
            # DESENGANCHE (corrección 2026-10-28; generalizada el 2026-10-29 a CUALQUIER
            # agente frenado, no sólo a los T endógenos: los de terapia son de tipo "IA" y
            # nunca lo recibían. Censo de politicas.THERAPIES: el único agente de terapia
            # frenable por PD-L1 es car_her2, el único con receptor y sin pdl1_resistant.
            # Medido antes de corregirlo: 1.556 frenos en 12 semillas, el 45 % en pasos
            # CONSECUTIVOS sobre la misma célula, rachas de hasta 25 pasos = 6,25 días.
            # Antes, la política del T no tenía
            # memoria y elegía siempre la diana visible más cercana: tras un freno volvía
            # a la MISMA célula PD-L1+ en el paso siguiente. Medido (luminal, neo 0,1,
            # semillas 1-3, 250 días): sin el flag, cada T frenado recibía sus 5 frenos
            # en 5 pasos seguidos sobre la misma célula y moría; con el flag adaptativo,
            # el 94 % de los ~114.000 frenos eran repeticiones sobre la misma célula, en
            # rachas de hasta 165 pasos (41 días, la vida entera del T).
            # POR QUÉ ES HONESTO: la señal PD-1 inhibe al CTL pero no lo ancla. Al
            # contrario, PD-1/PD-L1 BLOQUEA la señal de parada del TCR y el T sigue
            # moviéndose (Fife, Nat Immunol 2009;10:1185); los contactos T-diana duran
            # minutos-horas, no semanas. Se conserva TODO lo demás: el freno sigue
            # costando agotamiento (sin flag) o contando como estímulo (con flag), y el
            # T puede volver a esa célula pasado el desenganche si sigue ahí.
            # EXPECTATIVA, escrita ANTES de ejecutar: menos frenos totales (sobre todo
            # con el flag), T frenados que sobreviven más y matan algo más de células
            # PD-L1−; efecto pequeño en el tumor porque las alternativas PD-L1− al
            # alcance eran escasas (adyacentes en 3-8 % de los frenos sin flag; <2 % con
            # flag por el halo de IFN-γ). Veredictos de la batería sin cambios; riesgo
            # en T12 (con flag ya iba justo, recaída 1/3). T15 sigue fallando.
            # MEDIDO DESPUÉS: rachas sobre la misma célula 0 (antes hasta 165 pasos).
            # Sin flag: T1-T13 pasan igual; T15 falla igual. Con flag: T1-T12 pasan y
            # T13 FALLA por 1 célula sana (semilla 1, día 232) comida por una GUARDIANA
            # de puerta "nicho", no por un T: la trayectoria cambia y destapa un falso
            # positivo de esa puerta. No se ha tocado nada aquí para arreglarlo: la
            # causa era una madre SANA suprabasal y se corrigió en el tejido
            # (_healthy_stem_off_niche, 2026-10-28), no en el T ni en la puerta.
            ag.memory.setdefault("blocked_by", {})[c.id] = self.day + self.p.pd1_disengage_days
            self._milestone("pdl1", "Un clon con PD-L1 alto agota a los linfocitos T (escape inmune)")
            return
        if route in ("tcr", "car") and pdl1 and not ag.pdl1_resistant:
            # Solo se llega aquí con el anti-PD-1 circulando: el freno estaba y se ha quitado.
            self.stats["pdl1_brakes_blocked"] = self.stats.get("pdl1_brakes_blocked", 0) + 1
        prob = self.p.prob_attack_kill
        if t_natural and ag.t_state == "tex_prog":
            prob *= self.p.tex_prog_efficacy
        elif t_natural and ag.t_state == "tex_term":
            prob *= self.p.tex_term_efficacy
        kill_now = self.rng.random() < prob
        if adaptive and route in ("tcr", "car"):
            # El linfocito que ha reconocido suelta IFN-γ: PD-L1 en la diana y alrededor.
            self._induce_pdl1(q)
        if not kill_now:
            return
        if ag.sterilizes:
            # No la mata: la jubila. Telómero roto, sin más divisiones, y marcada con el
            # "cómeme" para que los macrófagos la retiren.
            if c.senescent:
                return
            c.senescent = True
            c.divisions = 0
            c.telomerase_ok = False
            ag.kills += 1
            self.stats["sterilized"] += 1
            if c.n_drivers == 0:
                self.stats["sterilized_healthy"] += 1
                self.stats["collateral"] += 1
            self._milestone("sterilized_one", "Una célula tumoral queda estéril: ya no puede dividirse")
            if ag.expansion and len(living) + 1 < MAX_AGENTS and self.rng.random() < ag.expansion:
                living.append(replace(ag, memory={}, exhaustion=0.0, kills=0, repairs=0))
            return
        del self.cells[q]
        ag.kills += 1
        if ag.kind == "MAC":
            self.stats["macrophage_eaten"] += 1
            self.stats["macrophage_eaten_healthy"] += c.n_drivers == 0
            self._milestone("macrophage", "Un macrófago se come una célula sin 'no me comas' y con 'cómeme'")
        self.stats["immune_kills" if ag.kind in ("T", "NK", "MAC") else "therapy_kills"] += 1
        if c.kind != LUMINAL or c.n_drivers == 0:
            self.stats["collateral"] += 1
            self._milestone("collateral", f"Daño colateral: {ag.name} mata una célula sana")
        if ag.spreads:
            # Replicación lítica: la progenie salta SOLO a vecinas tumorales en contacto.
            neighbors = [
                w
                for w in self._neighbors(q)
                if (u := self.cells.get(w)) is not None
                and u.kind == LUMINAL
                and u.n_drivers > 0
                and not u.necrosis
            ]
            self.rng.shuffle(neighbors)
            for w in neighbors[: int(ag.spreads)]:
                if len(living) + 1 >= MAX_AGENTS:
                    break
                living.append(replace(ag, pos=w, memory={}, exhaustion=0.0, kills=0, repairs=0, life=ag.life))
                self.stats["propagations"] += 1
        if (
            ag.expansion
            and not t_natural
            and len(living) + 1 < MAX_AGENTS
            and self.rng.random() < ag.expansion
        ):
            living.append(replace(ag, memory={}, exhaustion=0.0, kills=0, repairs=0))

    # ------------------------------------------------- inmunidad adaptativa (T15)
    # Ver biologia.Parameters.adaptive_immunity para fuentes y SUPUESTOS.
    def _pdl1(self, c: Cell) -> bool:
        """PD-L1 en superficie: por mutación (CD274) o inducido por IFN-γ (adaptativo)."""
        return c.pdl1 or self.t < c.pdl1_until

    def _induce_pdl1(self, q):
        """Resistencia adaptativa (Taube 2012; Ribas 2015): el IFN-γ del linfocito que
        ataca enciende PD-L1 en la diana y en las células epiteliales de alrededor."""
        r = self.p.ifng_radius
        until = self.t + int(self.p.induced_pdl1_days * STEPS_PER_DAY)
        y, x = q
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                c = self.cells.get((y + dy, x + dx))
                if c is not None and c.kind == LUMINAL and c.pdl1_until < until:
                    if c.pdl1_until < self.t:
                        self.stats["pdl1_induced"] = self.stats.get("pdl1_induced", 0) + 1
                    c.pdl1_until = until
        self._milestone(
            "adaptive_pdl1", "El IFN-γ del ataque enciende PD-L1 alrededor (resistencia adaptativa)"
        )

    def _stimulate_t(self, ag: Agent):
        """Un encuentro con el antígeno. El primero arranca la expansión clonal; los
        siguientes, acumulados, agotan (estimulación crónica, Wherry 2011)."""
        p = self.p
        ag.stimuli += 1
        if not ag.t_state:
            ag.t_state, ag.program = "effector", p.program_divisions
            self._milestone(
                "expansion", "Un linfocito T reconoce un neoantígeno y arranca su expansión clonal"
            )
        elif ag.t_state == "effector" and ag.stimuli >= p.exhaustion_threshold:
            ag.program = 0
            if ag.derived:
                self._to_terminal(ag)
            else:
                ag.t_state = "tex_prog"
                ag.life = max(ag.life, 160)
                self.stats["tex_progenitors"] = self.stats.get("tex_progenitors", 0) + 1
                self._milestone(
                    "tex", "Estimulación crónica: aparecen linfocitos T agotados progenitores (TCF1+)"
                )
        elif ag.t_state == "tex_prog" and self.rng.random() < p.prob_tex_terminal:
            self._to_terminal(ag)

    def _to_terminal(self, ag: Agent):
        ag.t_state, ag.program = "tex_term", 0
        ag.life = min(ag.life, int(self.p.tex_term_lifespan_days * STEPS_PER_DAY))
        self.stats["tex_terminals"] = self.stats.get("tex_terminals", 0) + 1

    def _clonal_room(self, living: list) -> bool:
        return self._n_clonal < self.p.max_clonal_t and len(living) + 1 < MAX_AGENTS

    def _spawn_clone(self, ag: Agent, living: list, **changes) -> Agent:
        daughter = replace(ag, memory={}, exhaustion=0.0, kills=0, repairs=0, clonal=True, **changes)
        living.append(daughter)
        self._n_clonal += 1
        self.stats["clonal_t"] = self.stats.get("clonal_t", 0) + 1
        return daughter

    def _contract(self, ag: Agent):
        """Fin del programa: el 90-95 % de los efectores muere (Kaech 2002)."""
        if self.rng.random() < self.p.t_contraction:
            ag.life = min(ag.life, self.rng.randint(1, int(self.p.contraction_days * STEPS_PER_DAY)))

    def _t_cycle(self, ag: Agent, living: list):
        """Un paso de vida de un linfocito T natural ya activado."""
        p, rng = self.p, self.rng
        if ag.t_state == "effector" and ag.program > 0:
            if rng.random() < p.prob_division_t:
                ag.program -= 1
                daughter = None
                if self._clonal_room(living):
                    daughter = self._spawn_clone(ag, living, life=160)
                if ag.program == 0:
                    self._contract(ag)
                    if daughter is not None:
                        self._contract(daughter)
        elif ag.t_state == "tex_prog":
            if self.t < self.pd1_until and rng.random() < p.prob_reinvigoration:
                # Anti-PD-1: el progenitor se autorrenueva y suelta progenie efectora
                # que vuelve a matar y acaba terminal (Siddiqui 2019; Miller 2019).
                if self._clonal_room(living):
                    self._spawn_clone(
                        ag,
                        living,
                        t_state="effector",
                        program=p.reinvigoration_divisions,
                        stimuli=p.exhaustion_threshold // 2,
                        derived=True,
                        life=160,
                    )
                    self.stats["reinvigorations"] = self.stats.get("reinvigorations", 0) + 1
                    self._milestone(
                        "reinvig",
                        "Anti-PD-1: los T agotados progenitores vuelven a dividirse (revigorización)",
                    )
            elif rng.random() < p.tex_self_renewal and self._clonal_room(living):
                self._spawn_clone(ag, living)

    def _repair(self, ag: Agent, c: Cell, q, gene: str):
        if gene not in ag.vectors or ag.payload <= 0 or c.kind != LUMINAL:
            return
        ag.payload -= 1
        if not c.genome.get(gene):
            self.stats["repairs_useless"] += 1
            return
        if self.rng.random() >= self.p.repair_efficacy:
            return
        had_broken_p53 = not c.p53_ok
        c.genome[gene] = 0
        c.update_phenotype()
        ag.repairs += 1
        self.stats["repairs"] += 1
        if c.n_drivers == 0:
            c.repaired = True
            self.stats["repairs_that_cure"] += 1
            # Chispa del efecto espectador: la célula corregida pasa a corregir vecinas.
            if self.p.correction_reach:
                c.corrector = True
            self._milestone("repaired", f"{ag.name} devuelve una célula a la normalidad ({gene} corregido)")
            if c.stem and ag.reprograms:
                c.guardian = ag.reprograms
                self.stats["guardians"] += 1
                self._milestone("rerouted", f"Reina curada y reprogramada como guardiana ({ag.reprograms})")
        elif had_broken_p53 and c.p53_ok and self.rng.random() < 0.5:
            # p53 restaurado detecta el resto de daños y ordena la apoptosis.
            del self.cells[q]
            self.stats["apoptosis_p53"] += 1
            self.stats["repairs_p53_lethal"] += 1
            self._milestone(
                "p53", "p53 restaurado: la célula tumoral reconoce su propio daño y se autodestruye"
            )
        else:
            # Gene corregido pero la célula sigue siendo tumoral: la reparación no cambió nada.
            self.stats["repairs_sterile"] += 1

    # ------------------------------------------------------------------ registro
    def _event(self, message: str):
        self.events.append(f"d{self.day:6.1f}  {message}")

    def _milestone(self, key: str, message: str):
        if key not in self._milestones:
            self._milestones.add(key)
            self._event(message)

    def _record_clone(self, c: Cell, pos):
        signature = c.altered
        if signature and signature not in self.clones:
            self.clones[signature] = len(self.clones) + 1
            self._event(
                f"clon #{self.clones[signature]}: {' + '.join(sorted(signature))} "
                f"→ {STATE_NAME[self.state(c, pos)]}"
            )

    def _record(self):
        h = self.history
        h["day"].append(self.day)
        for e in (*TUMOR_STATES, REPAIRED, NECROSIS):
            h[e].append(int(self.counts[e]))
        for kind in ("T", "NK", "AI"):
            h[kind].append(sum(1 for a in self.agents if a.kind == kind))

    def tumor_frequencies(self) -> dict[str, float]:
        tumor = [c for c in self.cells.values() if c.kind == LUMINAL and c.n_drivers > 0 and not c.necrosis]
        if not tumor:
            return {}
        return {g: sum(g in c.altered for c in tumor) / len(tumor) for g in GENES}

    def describe(self, pos) -> list[str]:
        c = self.cells.get(pos)
        if c is None:
            return [
                f"({pos[0]},{pos[1]}) {ZONE_LABEL[self.zone[pos]]}: vacío"
                + (" · vaso sanguíneo" if self.vessels[pos] else ""),
                f"O2 {self.O2[pos]:.2f}",
            ]
        hits = ", ".join(f"{g}×{h}" for g, h in c.genome.items() if h) or "ninguno"
        return [
            f"Célula #{c.id} · {TISSUE_LABEL[c.kind]} · {STATE_NAME[self.state(c, pos)]}",
            f"Zona {ZONE_LABEL[self.zone[pos]]} · O2 {self.O2[pos]:.2f} · ER{'+' if c.er else '-'}",
            f"Golpes: {hits}",
            f"Alterados: {', '.join(sorted(c.altered)) or '-'} · pasajeras {c.mutations}",
            f"p53 {'OK' if c.p53_ok else 'PERDIDO'} · HRD {'sí' if c.hrd else 'no'} · "
            f"MHC-I {'sí' if c.mhc else 'NO'} · PD-L1 {'alto' if c.pdl1 else 'bajo'} · "
            f"prolif ×{c.factor:.1f}",
            f"GUARDIANA (puerta {c.guardian})"
            if c.guardian
            else "Célula madre (reina)"
            if c.stem
            else f"Obrera · divisiones restantes {'∞' if c.divisions is None else c.divisions}"
            + (" · TELOMERASA BLOQUEADA" if not c.telomerase_ok else "")
            + (" · SENESCENTE (jubilada)" if c.senescent else "")
            + (" · CEBO CARGADO" if c.bait else ""),
        ]

    def summary(self) -> dict:
        tumor = int(sum(self.counts[e] for e in TUMOR_STATES))
        self.stats["senescent_cells"] = sum(
            1
            for c in self.cells.values()
            if c.kind == LUMINAL and c.n_drivers and (c.senescent or c.divisions == 0)
        )
        self.stats["aneuploid_cells"] = sum(1 for c in self.cells.values() if c.aneuploid)
        self.stats["aneuploid_healthy"] = sum(
            1 for c in self.cells.values() if c.aneuploid and c.n_drivers == 0
        )
        return dict(
            day=round(self.day, 1),
            tumor_cells=tumor,
            in_situ=int(self.counts[IN_SITU]),
            dangerous=int(self.counts[IN_SITU] + self.counts[INVASIVE] + self.counts[QUEEN]),
            queens=int(self.counts[QUEEN]),
            invasive_cells=int(self.counts[INVASIVE]),
            repaired_cells=int(self.counts[REPAIRED]),
            normal_cells=int(self.counts[NORMAL]),
            **self.stats,
        )
