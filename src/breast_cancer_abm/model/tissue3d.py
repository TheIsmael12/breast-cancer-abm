"""Paso A5: versión 3D del conducto mamario (cilindro, no círculo).

Por qué importa el 3D y no es un capricho estético: en 2D toda célula del tumor está
cerca del borde, así que cualquier terapia la alcanza. En 3D una masa tumoral es una
ESFERA, y su interior queda lejos de todo. La relación superficie/volumen cae con el
radio (∝ 1/r), de modo que cuanto más crece, menor fracción de sus células es accesible.
Esa es la razón principal por la que terapias que funcionan en placa fracasan en un tumor.

Reutiliza la misma biología y las mismas políticas que la versión 2D: aquí solo cambia
la geometría, la difusión y el vecindario (26 vecinos en vez de 8).

    bcabm study tissue3d --days 400 --therapy guardian+memory --therapy-day 200
"""

from __future__ import annotations

import argparse
import random
import time
from dataclasses import replace

import numpy as np

from breast_cancer_abm.model.agents import Agent, Analyze, Attack, Move, Observation, Perception, Repair
from breast_cancer_abm.model.biology import GENES, SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.policies import THERAPIES
from breast_cancer_abm.model.tissue import FIBROBLAST, LUMINAL, MYOEPITHELIAL, Cell, Simulation

LUMEN, EPITHELIUM, MYO, STROMA = range(4)
ZONE_NAME = ("lumen", "epithelium", "myoepithelium", "stroma")
NEIGHBORS = tuple(
    (dz, dy, dx) for dz in (-1, 0, 1) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dz or dy or dx
)
MAX_AGENTS = 1200
NECROSIS_LIFETIME = 40


def _sign(v: float) -> int:
    return (v > 0) - (v < 0)


class Simulation3D:
    """Conducto mamario como cilindro hueco a lo largo del eje Z."""

    def __init__(
        self,
        scenario: str = "luminal",
        seed: int | None = None,
        params: Parameters | None = None,
        side: int = 54,
        z_len: int = 40,
        radius: float = 13.0,
    ):
        scen = SCENARIOS[scenario]
        self.scenario = scenario
        self.p = replace(params or Parameters(), estrogen=scen.estrogen)
        self.rng = random.Random(seed)
        self.t = 0
        self._next_id = 0
        self.grid_shape = (z_len, side, side)
        self.cells: dict[tuple, Cell] = {}
        self.agents: list[Agent] = []
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
            healthy_stems_touched=0,
            propagations=0,
            steps_blocked=0,
        )
        self._build(scen, radius)
        for _ in range(40):
            self._diffuse()

    # ------------------------------------------------------------------ montaje
    def _new_id(self) -> int:
        self._next_id += 1
        return self._next_id

    def _build(self, scen, radius: float):
        z_len, side, _ = self.grid_shape
        self.center = (side - 1) / 2
        zz, yy, xx = np.mgrid[0:z_len, 0:side, 0:side]
        self.dist = np.hypot(yy - self.center, xx - self.center)  # distancia al eje del cilindro
        self.zone = np.full(self.grid_shape, STROMA, np.uint8)
        self.zone[self.dist < radius] = MYO
        self.zone[self.dist < radius - 1.5] = EPITHELIUM
        self.zone[self.dist < radius - 3.2] = LUMEN

        rand = self.rng
        self.vessels = np.zeros(self.grid_shape, bool)
        periductal = (self.dist >= radius + 2) & (self.dist < radius + 3.5)
        far = self.dist >= radius + 3.5
        noise = np.array(rand.random()) * 0 + np.random.default_rng(rand.randrange(1 << 30)).random(
            self.grid_shape
        )
        self.vessels[(periductal & (noise < 0.05)) | (far & (noise < 0.012))] = True

        for pos in zip(*np.nonzero(self.zone == EPITHELIUM)):
            self.cells[tuple(map(int, pos))] = Cell(
                self._new_id(),
                LUMINAL,
                dict(scen.germline),
                rand.random() < scen.er_fraction,
                stem=rand.random() < self.p.normal_stem_fraction,
            )
        for pos in zip(*np.nonzero(self.zone == MYO)):
            self.cells[tuple(map(int, pos))] = Cell(self._new_id(), MYOEPITHELIAL)
        for pos in zip(*np.nonzero((self.zone == STROMA) & ~self.vessels & (noise < 0.04))):
            self.cells[tuple(map(int, pos))] = Cell(self._new_id(), FIBROBLAST)

        if scen.seed:
            stems = [q for q, c in self.cells.items() if c.kind == LUMINAL and c.stem]
            c = self.cells[rand.choice(stems)]
            c.genome.update(scen.seed)
            c.neo = sum(rand.random() < self.p.prob_neoantigen for _ in scen.seed)
            c.update_phenotype()

        self.O2 = np.full(self.grid_shape, 0.8)
        self.Q = np.zeros(self.grid_shape)
        self.P = np.zeros(self.grid_shape)  # presión intersticial: empuja hacia fuera de la masa
        self.pressure_factor = 1.0
        self.cd47_until = -1
        self.tumor_vessels = np.zeros(self.grid_shape, bool)
        self.vessel_mode = "normal"

    # ------------------------------------------------------------------ bucle
    @property
    def day(self) -> float:
        return self.t / STEPS_PER_DAY

    def _neighbors(self, pos):
        z_len, side, _ = self.grid_shape
        z, y, x = pos
        for dz, dy, dx in NEIGHBORS:
            nz, ny, nx = z + dz, y + dy, x + dx
            if 0 <= nz < z_len and 0 <= ny < side and 0 <= nx < side:
                yield nz, ny, nx

    def step(self):
        self.t += 1
        order = list(self.cells.items())
        self.rng.shuffle(order)
        for pos, c in order:
            if self.cells.get(pos) is c:
                self._live(pos, c)
        self._diffuse()
        self._move_agents()

    def _live(self, pos, c: Cell):
        rng, p = self.rng, self.p
        if c.kind in (FIBROBLAST, MYOEPITHELIAL):
            return
        if c.necrosis:
            c.necrosis += 1
            if c.necrosis > NECROSIS_LIFETIME:
                del self.cells[pos]
            return
        o2 = self.O2[pos]
        # MUERTE POR O2 GRADUAL (variante B). Misma regla que tissue.py (ver allí la
        # justificación completa: Durand y Sham IJROBP 1998;42:711, y Papandreou et al.
        # Cancer Res 2005;65:3171). Hasta el 2026-10-29 este módulo seguía con la regla
        # vieja —sólo morían las de O2 < o2_necrosis, con p = 0,08— que hacía a las
        # hipóxicas casi inmortales. Las dos mallas tienen que envejecer igual o ninguna
        # medida 3D es comparable con la 2D.
        if o2 < p.o2_necrosis:
            p_necrosis = 0.08
        elif o2 < p.o2_hypoxia:
            p_necrosis = 0.08 * (p.o2_hypoxia - o2) / (p.o2_hypoxia - p.o2_necrosis)
        else:
            p_necrosis = 0.0
        if p_necrosis and rng.random() < p_necrosis:
            c.necrosis = 1
            self.stats["necrosis"] += 1
            return
        if c.n_drivers == 0 and rng.random() < p.prob_normal_death * (0.1 if c.stem else 1):
            del self.cells[pos]
            return
        if c.n_drivers and not c.stem:
            if c.divisions == 0 and rng.random() < p.prob_differentiated_death:
                del self.cells[pos]
                return
            sasp = any((v := self.cells.get(q)) is not None and v.senescent for q in self._neighbors(pos))
            if rng.random() < (
                p.plasticity * (3 if o2 < p.o2_hypoxia else 1) + (p.prob_sasp_plasticity if sasp else 0)
            ):
                c.stem, c.divisions = True, None
                self.stats["dedifferentiations"] += 1

        z = self.zone[pos]
        outside = (z == LUMEN and not c.anoikis_res) or (z in (MYO, STROMA) and not c.invasive)
        if outside and rng.random() < p.prob_anoikis:
            del self.cells[pos]
            self.stats["anoikis"] += 1
            return
        if c.invasive and z != STROMA and rng.random() < p.prob_invasion * 8:
            q = next(iter(self._neighbors(pos)))
            v = self.cells.get(q)
            if v is not None and v.kind == MYOEPITHELIAL:
                del self.cells[q]
                self.stats["membrane_breaches"] += 1
        if c.guardian and rng.random() < p.prob_phagocytosis * (
            p.internal_read_cost if c.guardian in Simulation.INTERNAL_GATES else 1.0
        ):
            self._phagocytose(pos, c)
            if self.cells.get(pos) is not c:
                return
        if c.senescent:
            return
        pdiv = p.prob_division * c.factor * ((1 + p.estrogen) if c.er else 1)
        if o2 < p.o2_hypoxia:
            pdiv *= 0.15
        if c.stem and c.n_drivers:
            pdiv *= p.stem_quiescence_factor
        if rng.random() < pdiv:
            self._divide(pos, c)

    def _allowed(self, c: Cell, origin, dest) -> bool:
        zd = self.zone[dest]
        if self.vessels[dest]:
            return False
        if zd == EPITHELIUM:
            return True
        if zd == LUMEN:
            return c.anoikis_res
        if zd == MYO:
            return c.invasive
        return c.invasive and self.zone[origin] in (MYO, STROMA)

    def _divide(self, pos, c: Cell):
        if c.divisions == 0:
            return
        gaps = [q for q in self._neighbors(pos) if q not in self.cells and self._allowed(c, pos, q)]
        if not gaps and c.prolif_signal:
            rivals = [
                q
                for q in self._neighbors(pos)
                if (v := self.cells.get(q)) is not None
                and v.kind == LUMINAL
                and v.n_drivers == 0
                and not v.necrosis
                and self._allowed(c, pos, q)
            ]
            if rivals and self.rng.random() < 0.3:
                del self.cells[rivals[0]]
                gaps = rivals[:1]
        if not gaps:
            return
        q = self.rng.choice(gaps)
        daughter = Cell(self._new_id(), LUMINAL, dict(c.genome), c.er, c.mutations)
        daughter.repaired, daughter.neo = c.repaired, c.neo
        daughter.cd47 = c.cd47 if self.t < self.cd47_until else True
        pa = self.p.prob_aneuploidy if not c.p53_ok else self.p.prob_aneuploidy_p53ok
        daughter.aneuploid = c.aneuploid or self.rng.random() < pa
        daughter.telomerase_ok = c.telomerase_ok
        immortal = (c.immortal or c.stem) and c.telomerase_ok
        if c.stem and c.n_drivers:
            daughter.stem = self.rng.random() < self.p.prob_symmetric_division
            daughter.divisions = None if (daughter.stem and immortal) else self.p.worker_divisions
            daughter.guardian = c.guardian if daughter.stem else None
        elif c.stem:
            daughter.stem = self.rng.random() < self.p.normal_stem_fraction
            daughter.guardian = c.guardian if daughter.stem else None
        elif immortal:
            daughter.divisions = c.divisions = None
        else:
            if c.divisions is None:
                c.divisions = self.p.worker_divisions
            c.divisions -= 1
            daughter.divisions = c.divisions
        for cel, site in ((c, pos), (daughter, q)):
            hits = self._mutate(cel)
            damage = hits + (1 if cel.hrd and self.rng.random() < self.p.prob_hrd_damage else 0)
            if damage and cel.p53_ok and self.rng.random() < 1 - (1 - self.p.prob_apoptosis_p53) ** damage:
                self.stats["apoptosis_p53"] += 1
                if cel is c:
                    del self.cells[pos]
                continue
            if cel is daughter:
                self.cells[q] = daughter

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
                c.divisions = self.p.worker_divisions
        return hits

    # ------------------------------------------------------------------ campos
    def _diffuse(self, iterations: int = 4):
        consumption = np.zeros(self.grid_shape)
        source = np.zeros(self.grid_shape)
        for pos, c in self.cells.items():
            if c.necrosis:
                continue
            consumption[pos] = 0.016 if c.n_drivers else 0.004
            if c.n_drivers:
                source[pos] = 1.0
        # Presión intersticial (ver tissue.py): la generan el tumor y los fibroblastos,
        # los linfáticos sanos la drenan fuera. En 3D es donde de verdad importa: la
        # meseta de presión cubre todo el interior de la esfera.
        if self.p.interstitial_pressure:
            pressure_src = source.copy()
            for pos, c in self.cells.items():
                if c.kind == FIBROBLAST:
                    pressure_src[pos] = 0.4
            pressure_src *= self.p.interstitial_pressure * self.pressure_factor
            drainage = np.where(pressure_src > 0, 0.02, 0.25)
            for _ in range(iterations):
                lap = -6 * self.P
                for axis_mid in (0, 1, 2):
                    lap = lap + np.roll(self.P, 1, axis_mid) + np.roll(self.P, -1, axis_mid)
                self.P += 0.14 * lap + 0.08 * pressure_src - drainage * self.P
                np.clip(self.P, 0.0, 1.0, out=self.P)
        for cfield, sink, supply in ((self.O2, consumption, None), (self.Q, 0.08, source)):
            for _ in range(iterations):
                lap = -6 * cfield
                for axis_mid in (0, 1, 2):
                    lap = lap + np.roll(cfield, 1, axis_mid) + np.roll(cfield, -1, axis_mid)
                cfield += 0.14 * lap - sink * cfield
                if supply is None:
                    cfield[self.vessels] = 1.0
                else:
                    cfield += 0.05 * supply
                np.clip(cfield, 0.0, 1.0, out=cfield)

    # ------------------------------------------------------------------ agentes
    def _observe(self, c: Cell, pos, dz, dy, dx) -> Observation:
        hypoxic = self.O2[pos] < self.p.o2_hypoxia
        neo = c.neo if (c.mhc and c.n_drivers > 0) else 0
        tissue = ("luminal", "myoepithelial", "fibroblast")[c.kind]
        return Observation(
            dy,
            dx,
            c.id,
            tissue,
            c.mhc,
            neo,
            c.pdl1,
            c.her2,
            c.senescent or c.n_drivers >= 2 or (hypoxic and c.n_drivers > 0),
            self.t - c.last_division < 8,
            hypoxic,
            c.stem,
            c.senescent,
            c.cd47,
            Simulation._exposes_eat_me(c, hypoxic),
            c.stem and not self._in_niche(pos),
        )

    def _perceive(self, ag: Agent) -> Perception:
        z, y, x = ag.pos
        self._dz = {}
        neighbors = []
        for dz in range(-1, 2):
            for dy in range(-2, 3):
                for dx in range(-2, 3):
                    c = self.cells.get((z + dz, y + dy, x + dx))
                    if c is not None and not c.necrosis:
                        o = self._observe(c, (z + dz, y + dy, x + dx), dz, dy, dx)
                        self._dz[(o.dy, o.dx)] = dz
                        neighbors.append(o)
        best, chemo = self.Q[z, y, x] + 1e-4, None
        for q in self._neighbors(ag.pos):
            if self.Q[q] > best:
                best, chemo = self.Q[q], (q[1] - y, q[2] - x)
        return Perception(
            neighbors,
            chemo,
            (_sign(self.center - y), _sign(self.center - x)),
            ZONE_NAME[self.zone[z, y, x]],
            float(self.O2[z, y, x]),
            self.day,
            self.rng,
        )

    def _move_agents(self):
        z_len, side, _ = self.grid_shape
        living = []
        for ag in self.agents:
            ag.life -= 1
            if ag.life <= 0 or ag.exhaustion >= 1:
                continue
            living.append(ag)
            p = self._perceive(ag)
            act = ag.policy(ag, p)
            if isinstance(act, Move):
                dz = self.rng.randint(-1, 1)
                dest = (
                    min(z_len - 1, max(0, ag.pos[0] + dz)),
                    min(side - 1, max(0, ag.pos[1] + _sign(act.dy))),
                    min(side - 1, max(0, ag.pos[2] + _sign(act.dx))),
                )
                if Simulation._blocked_by_pressure(self, ag.pos, dest):
                    self.stats["steps_blocked"] += 1
                else:
                    ag.pos = dest
                continue
            if not isinstance(act, (Attack, Analyze, Repair)):
                continue
            dz = self._dz.get((act.dy, act.dx), 0)
            q = (ag.pos[0] + dz, ag.pos[1] + act.dy, ag.pos[2] + act.dx)
            c = self.cells.get(q)
            if c is None or c.necrosis:
                continue
            if isinstance(act, Attack):
                self._attack(ag, c, q, living)
            elif isinstance(act, Analyze):
                ag.memory.setdefault("genomes", {})[c.id] = {g: h for g, h in c.genome.items() if h}
            else:
                self._repair(ag, c, q, act.gene)
        self.agents = living

    def _recognizes(self, ag: Agent, c: Cell, q) -> str | None:
        r = ag.receptors
        if "NONSPECIFIC" in r:
            return "nonspecific"
        if c.kind != LUMINAL:
            return None
        o = self._observe(c, q, 0, 0, 0)
        if "HER2" in r and c.her2:
            return "car"
        if "TCR" in r and o.neoantigens > 0:
            return "tcr"
        if "NK" in r and (not c.mhc or o.stress_ligands):
            return "nk"
        return None

    def _attack(self, ag: Agent, c: Cell, q, living: list):
        route = self._recognizes(ag, c, q)
        if route is None:
            return
        if route in ("tcr", "car") and c.pdl1 and not ag.pdl1_resistant:
            ag.exhaustion += 0.2
            self.stats["pdl1_brakes"] = self.stats.get("pdl1_brakes", 0) + 1
            # DESENGANCHE: el mismo artefacto que se corrigió en tissue.py el
            # 2026-10-28 y que aquí seguía sin corregir hasta el 2026-10-29.
            # agentes.t_cell_policy YA consultaba ag.memoria["frenado_por"],
            # pero este módulo no lo escribía nunca: sin memoria, hunt() vuelve en el
            # paso siguiente a la MISMA célula PD-L1+ y el T se queda anclado a ella
            # hasta agotarse. La señal PD-1 inhibe al CTL pero no lo ancla — al
            # contrario, bloquea la señal de parada del TCR y el T sigue moviéndose
            # (Fife, Nat Immunol 2009;10:1185). El freno sigue costando agotamiento y
            # el T puede volver a esa célula pasado el desenganche si sigue ahí.
            ag.memory.setdefault("blocked_by", {})[c.id] = self.day + self.p.pd1_disengage_days
            return
        if self.rng.random() >= self.p.prob_attack_kill:
            return
        del self.cells[q]
        ag.kills += 1
        self.stats["immune_kills" if ag.kind in ("T", "NK") else "therapy_kills"] += 1
        if c.kind != LUMINAL or c.n_drivers == 0:
            self.stats["collateral"] += 1
        if ag.spreads:
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
                living.append(replace(ag, pos=w, memory={}, exhaustion=0.0, kills=0, repairs=0))
                self.stats["propagations"] += 1
        if ag.expansion and len(living) + 1 < MAX_AGENTS and self.rng.random() < ag.expansion:
            living.append(replace(ag, memory={}, exhaustion=0.0, kills=0, repairs=0))

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
            if c.stem and ag.reprograms:
                c.guardian = ag.reprograms
                self.stats["guardians"] += 1
        elif had_broken_p53 and c.p53_ok and self.rng.random() < 0.5:
            del self.cells[q]
            self.stats["apoptosis_p53"] += 1

    def _in_niche(self, pos) -> bool:
        """Capa basal: epitelio en contacto con el mioepitelio (ver tissue.py)."""
        if self.zone[pos] != EPITHELIUM:
            return False
        return any(self.zone[q] == MYO for q in self._neighbors(pos))

    def _phagocytose(self, pos, c: Cell):
        for q in self._neighbors(pos):
            v = self.cells.get(q)
            if v is None or v.kind != LUMINAL or v.necrosis or v.guardian:
                continue
            o = self._observe(v, q, 0, 0, 0)
            target = Simulation._gate(self, c.guardian, v, o)
            if not target:
                continue
            del self.cells[q]
            self.stats["phagocytosed"] += 1
            if v.n_drivers == 0:
                self.stats["phagocytosed_healthy"] += 1
                self.stats["collateral"] += 1
            if self.rng.random() < 0.5:
                self._divide(pos, c)
            return

    # ------------------------------------------------------------------ terapias
    def _random_vessel(self):
        zs, ys, xs = np.nonzero(self.vessels)
        i = self.rng.randrange(len(zs))
        return int(zs[i]), int(ys[i]), int(xs[i])

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

    def treat_vessels(self, mode: str):
        """Mismos modos que en 2D. En 3D solo cambia el efecto sobre el riego; la caída
        de presión al normalizar es la misma (cerrar las fugas baja la PIF)."""
        self.vessel_mode = mode
        if mode == "normalize":
            self.pressure_factor *= 1 - self.p.pressure_drop_normalize
        # OJO: el 3D todavía NO modela la malla (CAFs, colágeno, LOX) ni, por tanto,
        # `treat_matrix`. Solo tiene presión y geometría. degrade_matrix aquí se limita
        # a retirar fibroblastos.

    # ------------------------------------------------------------------ sin bisturí
    # Mismas reglas que tissue.py (ablacionar, ctdna, carcinoma_foci), en 3D.
    # Estudio: docs/avances/2026-11-05-tratar-sin-operar.md
    def is_carcinoma(self, pos, c: Cell) -> bool:
        """El mismo criterio que 'peligrosas' en summary(): in situ, invasivo o reina."""
        return (
            c.kind == LUMINAL
            and not c.necrosis
            and c.n_drivers > 0
            and (c.stem or c.anoikis_res or self.zone[pos] in (MYO, STROMA))
        )

    def carcinoma_foci(self) -> list[list]:
        """Componentes conexas (26 vecinas) de células de carcinoma."""
        sites = {q for q, c in self.cells.items() if self.is_carcinoma(q, c)}
        foci, seen = [], set()
        for q in sites:
            if q in seen:
                continue
            focus, stack_ = [], [q]
            seen.add(q)
            while stack_:
                z, y, x = stack_.pop()
                focus.append((z, y, x))
                for dz, dy, dx in NEIGHBORS:
                    v = (z + dz, y + dy, x + dx)
                    if v in sites and v not in seen:
                        seen.add(v)
                        stack_.append(v)
            foci.append(focus)
        return foci

    def ablate(
        self,
        coverage: float = 0.95,
        margin: int = 2,
        rim: int = 1,
        rim_survival: float = 0.3,
        guided: bool = True,
        visible_threshold: int = 1,
    ) -> dict:
        """Ablación en 3D. guiada=True: bola de radio `margen` alrededor de cada célula de
        carcinoma de los focos visibles (≥ visible_threshold células). guiada=False: una
        esfera alrededor del centroide de todas las células con driver (el «disco»).
        Las luminales destruidas quedan como necrosis; las mioepiteliales y los
        fibroblastos se retiran (en 3D _live no limpia su necrosis)."""
        if guided:
            foci = self.carcinoma_foci()
            targets = [q for f in foci if len(f) >= visible_threshold for q in f]
            seen = sum(len(f) >= visible_threshold for f in foci)
            ball_radius = margin
        else:
            targets = [
                q for q, c in self.cells.items() if c.kind == LUMINAL and c.n_drivers and not c.necrosis
            ]
            seen = 1 if targets else 0
        if not targets:
            return dict(destroyed=0, healthy=0, residual=0, foci_seen=0)
        dmin: dict = {}
        if guided:
            ofs = [
                (dz, dy, dx)
                for dz in range(-margin, margin + 1)
                for dy in range(-margin, margin + 1)
                for dx in range(-margin, margin + 1)
                if (dz * dz + dy * dy + dx * dx) ** 0.5 <= ball_radius
            ]
            for z, y, x in targets:
                for dz, dy, dx in ofs:
                    q = (z + dz, y + dy, x + dx)
                    d = (dz * dz + dy * dy + dx * dx) ** 0.5
                    if d < dmin.get(q, 1e9):
                        dmin[q] = d
            limit = margin
        else:
            cz, cy, cx = (sum(q[i] for q in targets) / len(targets) for i in range(3))
            dist = sorted(((q[0] - cz) ** 2 + (q[1] - cy) ** 2 + (q[2] - cx) ** 2) ** 0.5 for q in targets)
            limit = dist[min(len(dist) - 1, int(len(dist) * coverage))] + margin
            for q in self.cells:
                d = ((q[0] - cz) ** 2 + (q[1] - cy) ** 2 + (q[2] - cx) ** 2) ** 0.5
                if d <= limit:
                    dmin[q] = d
        destroyed = healthy = 0
        for q, d in dmin.items():
            c = self.cells.get(q)
            if c is None or c.necrosis:
                continue
            if d > limit - rim and self.rng.random() < rim_survival:
                continue
            if c.kind == LUMINAL and c.n_drivers:
                destroyed += 1
            else:
                healthy += 1
            if c.kind == LUMINAL:
                c.necrosis = 1
            else:
                del self.cells[q]
        residual = sum(1 for q, c in self.cells.items() if self.is_carcinoma(q, c))
        return dict(destroyed=destroyed, healthy=healthy, residual=residual, foci_seen=seen)

    def ctdna(self, fragments_per_cell: float, necrosis_weight: float = 3.0) -> dict:
        alive = sum(1 for q, c in self.cells.items() if self.is_carcinoma(q, c))
        dead = sum(1 for c in self.cells.values() if c.kind == LUMINAL and c.n_drivers and c.necrosis)
        lam = fragments_per_cell * (alive + necrosis_weight * dead)
        n = int(np.random.default_rng(self.rng.getrandbits(32)).poisson(lam))
        return dict(positive=n > 0, fragments=n, expected=lam, alive=alive)

    def degrade_matrix(self, fraction: float = 0.25):
        removed = 0
        for pos, c in list(self.cells.items()):
            if c.kind == FIBROBLAST and self.rng.random() < fraction:
                del self.cells[pos]
                removed += 1
        self.pressure_factor *= 1 - self.p.pressure_drop_matrix * fraction
        return removed

    def damage_telomeres(self):
        dead = retired = healthy = 0
        for pos, c in list(self.cells.items()):
            if c.kind != LUMINAL or not (c.immortal or c.stem) or not c.telomerase_ok:
                continue
            if self.rng.random() >= self.p.efficacy_6tiodg:
                continue
            if c.n_drivers == 0 and c.stem:
                healthy += 1
            if c.p53_ok:
                del self.cells[pos]
                self.stats["apoptosis_p53"] += 1
                dead += 1
            else:
                c.senescent, c.divisions = True, 0
                retired += 1
        self.stats["healthy_stems_touched"] += healthy
        return dead, retired, healthy

    # ------------------------------------------------------------------ informe
    def summary(self) -> dict:
        tumor = dangerous = queens = senescent_cells = normal_cells = 0
        for pos, c in self.cells.items():
            if c.kind != LUMINAL or c.necrosis:
                continue
            if c.n_drivers == 0:
                normal_cells += 1
                continue
            tumor += 1
            if c.stem:
                queens += 1
            if c.senescent or c.divisions == 0:
                senescent_cells += 1
            if c.stem or c.anoikis_res or self.zone[pos] in (MYO, STROMA):
                dangerous += 1
        return dict(
            day=round(self.day, 1),
            tumor_cells=tumor,
            dangerous=dangerous,
            queens=queens,
            senescent_cells=senescent_cells,
            normal_cells=normal_cells,
            agents=len(self.agents),
            **self.stats,
        )

    def accessible_fraction(self) -> float:
        """Fracción de células tumorales con al menos un hueco adyacente: por ahí, y solo
        por ahí, puede llegar una célula terapéutica. Es la medida clave del 3D."""
        tumor = [q for q, c in self.cells.items() if c.kind == LUMINAL and c.n_drivers and not c.necrosis]
        if not tumor:
            return 1.0
        free = sum(any(w not in self.cells for w in self._neighbors(q)) for q in tumor)
        return free / len(tumor)


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--days", type=float, default=400)
    ap.add_argument("--therapy-day", type=float, default=200)
    ap.add_argument("--therapy", default=None)
    ap.add_argument("--tiodg", action="store_true")
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--route", default="vein", choices=("vein", "intratumoral", "needle"))
    ap.add_argument("--dose", type=int, default=None, help="células por terapia (por defecto, la suya)")
    ap.add_argument("--matrix", type=float, default=0.0, help="fracción de matriz degradada (0-1)")
    ap.add_argument("--neoantigens", type=float, default=0.1)
    a = ap.parse_args(argv)

    print(
        f"3D · {a.scenario} · {a.days:g} días · terapia {a.therapy or 'none'} "
        f"el día {a.therapy_day:g} · vía {a.route}"
        f"{f' · dosis {a.dose}' if a.dose else ''}"
        f"{f' · matriz -{a.matrix:.0%}' if a.matrix else ''}"
        f"{' + 6-tio-dG' if a.tiodg else ''}"
    )
    for seed in range(a.seeds):
        sim = Simulation3D(a.scenario, seed, Parameters(prob_neoantigen=a.neoantigens))
        t0 = time.time()
        therapy_step = int(a.therapy_day * STEPS_PER_DAY)
        for step in range(int(a.days * STEPS_PER_DAY)):
            if step == therapy_step:
                if a.matrix:
                    sim.degrade_matrix(a.matrix)
                if a.therapy:
                    sim.release_therapy(a.therapy, a.dose, a.route)
                if a.tiodg:
                    sim.damage_telomeres()
            sim.step()
        r = sim.summary()
        print(
            f"  semilla {seed}: tumorales {r['tumor_cells']:6} · peligrosas {r['dangerous']:6} "
            f"· reinas {r['queens']:5} · senescentes {r['senescent_cells']:5} "
            f"· colateral {r['collateral']:4} · accesible {sim.accessible_fraction():.0%} "
            f"({time.time() - t0:.0f}s)"
        )


if __name__ == "__main__":
    main()
