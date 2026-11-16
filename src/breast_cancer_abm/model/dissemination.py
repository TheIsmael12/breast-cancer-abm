"""Diseminación: las células que se escapan antes de operar.

    bcabm study dissemination                 # valida la fase 1 (D1-D4)
    bcabm study dissemination --seeds 10
    bcabm study dissemination --d4            # D4: varía los SUPUESTO (32 semillas)

QUÉ AÑADE AL MODELO. El simulador base es un conducto mamario: un compartimento.
Pero lo que mata al sujeto no es el tumor que se extirpa, es la célula que ya se
había ido — y está en otro órgano, no se divide, y nadie la ve. Este módulo añade
los dos compartimentos que faltan:

    mama (tissue.py)  ──intravasación──>  circulación  ──siembra──>  médula ósea
                                          (casi todas              (dormancia,
                                           mueren)                  años)

NADA DE ESTO SE LE PROGRAMA AL MODELO como objetivo. Se deja que el tumor crezca,
invada y siembre por su cuenta, y después se mide contra los criterios de
`DISEMINACION.md`. Si no los pasa, la fase 1 no está lista.

Mecanismo de la dormancia, que NO es inventado (Ghajar, Nat Cell Biol 2013):
  · el microvaso ESTABLE mantiene dormida a la célula (trombospondina-1)
  · el BROTE vascular la despierta (TGF-β1, periostina)
Desde la v2 (2026-10-28) esto es literal en el código: la médula es un conjunto de
microvasos que brotan y se estabilizan; la DTC se aloja junto a uno y su estado
(dormida o no) es el de su vaso, y cambia cuando el vaso cambia. Ya no hay una
probabilidad de «entrar dormida». Diseño y supuestos: DISEMINACION.md, «Fase 1 · v2».
Y lo que NO hay que modelar como protección: la quiescencia no protege del ataque
de células T (93 % de eliminación, Goddard, Cancer Cell 2024). Lo que protege a la
DTC es la ESCASEZ RELATIVA — son poquísimas y raramente se encuentran con un
linfocito que las reconozca. Es una propiedad espacial, no bioquímica.
"""

from __future__ import annotations

import argparse
import math
import random
import statistics
from dataclasses import dataclass, field

from breast_cancer_abm.model.biology import STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import LUMINAL, Simulation


# --------------------------------------------------------------------- parámetros
@dataclass
class DisseminationParams:
    """Todo lo que aquí no tiene fuente va marcado como SUPUESTO y entra en D4."""

    # DOS VÍAS DE ESCAPE, no una. El modelo original solo tenía la segunda, y por eso
    # D1 daba un interruptor (0 % con tumor pequeño, 100 % con tumor grande) en vez del
    # 30-40 % real: suponía, sin decirlo, que solo disemina el tumor invasivo.
    #
    # (1) DISEMINACIÓN TEMPRANA. Células que escapan desde lesiones AÚN NO INVASIVAS.
    #     Demostrado en modelos murinos (Hosseini et al. y Harper et al., Nature 2016):
    #     cambió cómo se entiende la cascada metastásica. Explica que la fracción de
    #     pacientes con DTC sea intermedia y relativamente independiente del tamaño:
    #     unas siembran pronto y otras no, en vez de sembrar todas al crecer.
    prob_early_dissemination: float = 5e-6  # CALIBRADO (ver nota al final)
    early_competence: float = 0.25  # SUPUESTO: fracción de tumores capaces de ello
    # (2) Intravasación clásica: una célula invasiva que toca un vaso se suelta.
    prob_intravasation: float = 3e-4  # CALIBRADO (ver nota al final)
    # Tránsito: la cascada metastásica es MUY ineficiente (criterio D3).
    prob_survives_transit: float = 0.02  # SUPUESTO dentro de un rango muy aceptado
    prob_dies_in_niche: float = 0.0015  # SUPUESTO, por paso: el nicho no es gratis
    # NICHO MEDULAR (v2). La dormancia ya NO es un parámetro (`prob_entra_dormida` y
    # `prob_brote_vascular` eliminados el 2026-10-28): sale del vaso junto al que se
    # aloja la DTC. Vaso estable -> dormida (trombospondina-1); vaso en brote ->
    # despierta (TGF-β1, periostina). Ghajar et al., Nat Cell Biol 2013.
    # MEDIDO: volumen de vaso tipo H (angiogénico) / volumen vascular total en fémur
    # humano, 18,9 % a 45-59 años y 13,9 % a 60-75 (Wang et al., Cell Death Dis 2017).
    sprouting_vessel_fraction: float = 0.16
    sprout_duration_days: float = 14.0  # SUPUESTO: vida media de un brote
    # Peso del vaso estable al alojarse (1 = proporcional al volumen vascular). Price
    # 2016 (nicho perisinusoidal E-selectina/SDF-1) sugiere > 1; magnitud SUPUESTA.
    stable_homing_bias: float = 1.0
    n_vessels: int = 500  # SUPUESTO de discretización
    # La DTC vuelve a dormirse si su vaso se estabiliza (Sosa et al., Nat Rev Cancer
    # 2014: dormancia reversible dictada por el entorno). Respaldo cualitativo.
    reversible_dormancy: bool = True
    # Detección: sensibilidad imperfecta, y esto SÍ está medido.
    ihc_sensitivity: float = 0.24  # biopsia de médula por IHC
    cytometry_sensitivity: float = 0.47  # la misma población, otro método


@dataclass
class DTC:
    """Una célula tumoral diseminada, viviendo en la médula."""

    id: int
    n_drivers: int
    stem: bool
    vessel: int = 0
    dormant: bool = True
    arrival_day: float = 0.0
    wake_day: float | None = None


@dataclass
class Organism:
    """El sujeto completo: la mama más los compartimentos que faltaban."""

    sim: Simulation
    p: DisseminationParams = field(default_factory=DisseminationParams)
    marrow: list[DTC] = field(default_factory=list)
    in_transit: int = 0
    stats: dict = field(
        default_factory=lambda: dict(
            intravasations=0,
            dead_in_transit=0,
            seeded=0,
            awakened=0,
            redormant=0,
            dead_in_niche=0,
            early_route=0,
            classic_route=0,
        )
    )
    # La médula tiene su propio generador: variar el nicho no cambia tumor ni siembra.
    rng_marrow: random.Random | None = None
    sprouting_vessels: list[bool] = field(default_factory=list)
    _next_id: int = 0
    _competent: bool | None = None

    def __post_init__(self):
        if self.rng_marrow is None:
            self.rng_marrow = random.Random(0)
        f = self.p.sprouting_vessel_fraction
        # La médula adulta ya está en equilibrio: estado inicial estacionario.
        self.sprouting_vessels = [self.rng_marrow.random() < f for _ in range(self.p.n_vessels)]
        dt = 1 / STEPS_PER_DAY
        k_minus = 1 / self.p.sprout_duration_days
        k_plus = f * k_minus / (1 - f) if f < 1 else float("inf")
        self._p_sprout = 1 - math.exp(-k_plus * dt)
        self._p_stabilize = 1 - math.exp(-k_minus * dt)

    # ---------------------------------------------------------------- un paso
    def step(self):
        self.sim.step()
        self._intravasate()
        self._niche()

    def _competent_early(self) -> bool:
        """¿Este sujeto tiene un tumor capaz de diseminar antes de ser invasivo?

        Se decide una vez por sujeto, no por célula: es una propiedad del clon, y es
        lo que produce una fracción intermedia de pacientes con DTC en vez de un
        todo-o-nada ligado al tamaño."""
        if self._competent is None:
            self._competent = self.sim.rng.random() < self.p.early_competence
        return self._competent

    def _intravasate(self):
        """Dos vías: la temprana (no invasiva) y la clásica (invasiva junto a vaso)."""
        sim, p = self.sim, self.p
        # (1) VÍA TEMPRANA: no necesita invasión ni contacto con vaso
        if p.prob_early_dissemination and self._competent_early():
            for pos, c in list(sim.cells.items()):
                if c.kind != LUMINAL or not c.n_drivers or c.necrosis:
                    continue
                if sim.rng.random() >= p.prob_early_dissemination:
                    continue
                self.stats["intravasations"] += 1
                self.stats["early_route"] += 1
                if sim.rng.random() >= p.prob_survives_transit:
                    self.stats["dead_in_transit"] += 1
                    continue
                self._seeding(c)
        # (2) VÍA CLÁSICA
        for pos, c in list(sim.cells.items()):
            if c.kind != LUMINAL or not c.n_drivers or c.necrosis or not c.invasive:
                continue
            if not any(sim.vessels[q] for q in sim._neighbors(pos)):
                continue
            if sim.rng.random() >= p.prob_intravasation:
                continue
            self.stats["intravasations"] += 1
            # Tránsito: la inmensa mayoría no sobrevive (criterio D3)
            if sim.rng.random() >= p.prob_survives_transit:
                self.stats["dead_in_transit"] += 1
                continue
            self.stats["classic_route"] += 1
            self._seeding(c)

    def _seeding(self, c):
        """La célula sobrevivió al tránsito y se aloja junto a un microvaso.

        Su estado NO se sortea: es el del vaso donde cae (estable -> dormida)."""
        p, rng = self.p, self.rng_marrow
        weight_list = [1.0 if b else p.stable_homing_bias for b in self.sprouting_vessels]
        v = rng.choices(range(len(weight_list)), weights=weight_list)[0]
        self._next_id += 1
        sprout = self.sprouting_vessels[v]
        self.marrow.append(
            DTC(
                id=self._next_id,
                n_drivers=c.n_drivers,
                stem=c.stem,
                vessel=v,
                dormant=not sprout,
                arrival_day=self.sim.day,
                wake_day=self.sim.day if sprout else None,
            )
        )
        self.stats["seeded"] += 1

    def _niche(self):
        """La médula: los vasos brotan y se estabilizan, y las DTC siguen a su vaso."""
        p, rng = self.p, self.rng_marrow
        vb = self.sprouting_vessels
        for i, b in enumerate(vb):
            if rng.random() < (self._p_stabilize if b else self._p_sprout):
                vb[i] = not b
        alive = []
        for d in self.marrow:
            if rng.random() < p.prob_dies_in_niche:
                self.stats["dead_in_niche"] += 1
                continue
            sprout = vb[d.vessel]
            if d.dormant and sprout:
                # El brote vascular rompe la dormancia (TGF-β1, periostina)
                d.dormant = False
                d.wake_day = self.sim.day
                self.stats["awakened"] += 1
            elif not d.dormant and not sprout and p.reversible_dormancy:
                # El vaso se estabiliza: vuelve la trombospondina-1
                d.dormant = True
                self.stats["redormant"] += 1
            alive.append(d)
        self.marrow = alive

    # ------------------------------------------------------------- observación
    def biopsy(self, method_name: str = "ihc") -> int:
        """Mirar la médula con sensibilidad IMPERFECTA, que es lo que pasa de verdad.

        Una terapia perfecta sobre una población que no sabes detectar no sirve."""
        s = self.p.ihc_sensitivity if method_name == "ihc" else self.p.cytometry_sensitivity
        return sum(1 for _ in self.marrow if self.sim.rng.random() < s)

    def summary(self) -> dict:
        dormant_count = sum(1 for d in self.marrow if d.dormant)
        return dict(
            dtc=len(self.marrow),
            dormant_count=dormant_count,
            awake=len(self.marrow) - dormant_count,
            dormant_frac=dormant_count / len(self.marrow) if self.marrow else 0.0,
            tumor=sum(
                1 for c in self.sim.cells.values() if c.kind == LUMINAL and c.n_drivers and not c.necrosis
            ),
            **self.stats,
        )


# ------------------------------------------------------------------ validación
def run(seed: int, days: float, p: DisseminationParams | None = None) -> dict:
    sim = Simulation("luminal", seed, Parameters(prob_neoantigen=0.1))
    org = Organism(sim, p or DisseminationParams(), rng_marrow=random.Random(1_000_003 + seed))
    for _ in range(int(days * STEPS_PER_DAY)):
        org.step()
    r = org.summary()
    r["biopsy_ihc"] = org.biopsy("ihc")
    r["biopsy_cytometry"] = org.biopsy("cytometry")
    return r


def _run_task(args):
    return run(*args)


def run_batch(seeds, days, p=None, processes=None) -> list[dict]:
    """Las mismas corridas que `correr`, repartidas por semillas entre procesos."""
    tasks = [(s, days, p) for s in seeds]
    if processes == 1 or len(tasks) == 1:
        return [_run_task(t) for t in tasks]
    import multiprocessing as mp

    with mp.Pool(processes or max(1, (mp.cpu_count() or 2) - 4)) as pool:
        return pool.map(_run_task, tasks)


def _d2_pooled(res: list[dict]) -> float:
    tot = sum(r["dtc"] for r in res)
    return sum(r["dormant_count"] for r in res) / tot if tot else 0.0


def se_d2_bootstrap(res: list[dict], n: int = 1000) -> float:
    """Error estándar de la D2 agregada remuestreando PACIENTES (generador fijo)."""
    rng = random.Random(0)
    vals = [_d2_pooled(rng.choices(res, k=len(res))) for _ in range(n)]
    return statistics.stdev(vals)


def metrics(res: list[dict]) -> dict:
    """D1, D2 y D3 tal como los define DISEMINACION.md (v2).

    D2 es la fracción AGREGADA (dormidas de todos / DTC de todos), anunciada antes de
    correr la v2; la mediana por paciente (v1) se devuelve como `d2_median`, sin
    veredicto: con 1-2 DTC por paciente solo vale 0, ½ o 1 y se satura."""
    with_dtc = sum(r["dtc"] > 0 for r in res)
    fracs = [r["dormant_frac"] for r in res if r["dtc"] > 0]
    total_intra = sum(r["intravasations"] for r in res)
    total_dead = sum(r["dead_in_transit"] for r in res)
    d1 = with_dtc / len(res)
    d2 = _d2_pooled(res)
    d3 = total_dead / total_intra if total_intra else 0.0
    return dict(
        d1=d1,
        d2=d2,
        d3=d3,
        with_dtc=with_dtc,
        n=len(res),
        intra=total_intra,
        dtc=sum(r["dtc"] for r in res),
        d2_median=statistics.median(fracs) if fracs else 0.0,
        ok1=0.25 <= d1 <= 0.50,
        ok2=d2 > 0.80,
        ok3=d3 > 0.90,
    )


# Supuestos que entran en D4, y cómo se varía cada uno (fijado antes de correr).
D4_ASSUMPTIONS = (
    "early_competence",
    "prob_survives_transit",
    "prob_dies_in_niche",
    "sprout_duration_days",
    "stable_homing_bias",
    "n_vessels",
    "sprouting_vessel_fraction",
)
# Los parámetros MEDIDOS no se varían ×0,5/×2 sino a los extremos de su rango medido.
MEASURED_RANGES = {"sprouting_vessel_fraction": [("60-75 años", 0.139), ("45-59 años", 0.189)]}


def _variants(name: str, value: float) -> list[tuple[str, float]]:
    if name in MEASURED_RANGES:
        return MEASURED_RANGES[name]
    if isinstance(value, int):
        return [("×0,5", max(1, round(value * 0.5))), ("×2", value * 2)]
    cap = 1.0 if name.startswith("prob_") or name == "early_competence" else None
    double = value * 2 if cap is None else min(cap, value * 2)
    return [("×0,5", value * 0.5), ("×2", double)]


def control_d4(seeds, days, processes=None):
    """D4: ¿depende el resultado de un parámetro inventado sin respaldo?

    CRITERIO DE ACEPTACIÓN v2 (escrito antes de ejecutar; copia en DISEMINACION.md,
    apartado «D4 v2»):
      · Cada supuesto de D4_ASSUMPTIONS se corre a ×0,5 y ×2 (probabilidades recortadas
        a [0, 1]); los MEDIDOS (sprouting_vessel_fraction) a los extremos de su rango
        medido (0,139 y 0,189).
      · Mismas semillas en la base y en todas las variantes.
      · Para cada supuesto y cada Dk (D1-D3) el veredicto PASA/FALLA «se mantiene» si
        es igual que en la base en ambas variantes; «cambia» si no, y la diferencia
        supera el ruido (D1: 2·√(2·p(1−p)/n); D2 agregada: 2·√2·EE bootstrap por
        pacientes en la base; D3 tal cual); «no concluyente» si cambia dentro del ruido.
      · D4 PASA si ningún supuesto cambia ningún veredicto; FALLA si alguno lo cambia
        por encima del ruido; NO CONCLUYENTE si solo hay cambios dentro del ruido.
      · Una dependencia por construcción que cambie el veredicto cuenta como FALLO.
      · Control estructural INFORMATIVO (no cuenta para D4): despertar irreversible.
    """
    base_p = DisseminationParams()
    n = len(seeds)
    print(f"D4 v2 · control de dependencia · {n} semillas · {days:g} días\n")
    res_baseline = run_batch(seeds, days, base_p, processes)
    base = metrics(res_baseline)

    def row(tg, m):
        mark = lambda ok: "ok" if ok else "XX"
        print(
            f"  {tg:42s} D1 {m['d1']:5.0%} [{mark(m['ok1'])}]  "
            f"D2 {m['d2']:5.1%} [{mark(m['ok2'])}]  "
            f"D3 {m['d3']:6.1%} [{mark(m['ok3'])}]  "
            f"(DTC {m['dtc']}, mediana D2 {m['d2_median']:.0%})",
            flush=True,
        )

    row("BASE", base)
    noise_d1 = 2 * (2 * max(base["d1"] * (1 - base["d1"]), 1 / n) / n) ** 0.5
    noise_d2 = 2 * 2**0.5 * se_d2_bootstrap(res_baseline)
    noise = {"d1": noise_d1, "d2": noise_d2}
    worst = "PASA"
    for name in D4_ASSUMPTIONS:
        v0 = getattr(base_p, name)
        print(f"\n  {name} (base {v0:g})")
        states = {}
        for tg, v in _variants(name, v0):
            m = metrics(run_batch(seeds, days, DisseminationParams(**{name: v}), processes))
            row(f"  {tg:10s} = {v:g}", m)
            for k, dk in (("ok1", "d1"), ("ok2", "d2"), ("ok3", "d3")):
                if m[k] == base[k]:
                    continue
                inside = dk in noise and abs(m[dk] - base[dk]) <= noise[dk]
                status = "no concluyente" if inside else "CAMBIA"
                if states.get(dk) != "CAMBIA":
                    states[dk] = status
        for dk in ("d1", "d2", "d3"):
            status = states.get(dk, "se mantiene")
            if status != "se mantiene":
                print(f"    -> {dk.upper()}: veredicto {status}")
            if status == "CAMBIA":
                worst = "FALLA"
            elif status == "no concluyente" and worst == "PASA":
                worst = "NO CONCLUYENTE"
    print(f"\n  (ruido tolerado: D1 ±{noise_d1:.0%}, D2 ±{noise_d2:.1%})")
    print(f"\n  D4 {worst}")
    print("\n  Control estructural (INFORMATIVO, no cuenta para D4):")
    m = metrics(run_batch(seeds, days, DisseminationParams(reversible_dormancy=False), processes))
    row("  despertar irreversible", m)
    return worst


def validate(seeds, days, processes=None):
    """Criterios D1-D4 de DISEMINACION.md. Se pasan o no se pasan."""
    print(f"FASE 1 · {len(seeds)} semillas · {days:g} días (momento del diagnóstico)\n")
    res = run_batch(seeds, days, None, processes)
    for s, r in zip(seeds, res):
        print(
            f"  s{s:2d}  tumor {r['tumor']:6d} · DTC {r['dtc']:4d} "
            f"({r['dormant_frac']:.0%} dormidas) · intravasan {r['intravasations']:5d} "
            f"· mueren en tránsito {r['dead_in_transit']:5d} · "
            f"biopsia IHC {r['biopsy_ihc']:3d}",
            flush=True,
        )

    m = metrics(res)
    with_dtc, d1, d2, d3, total_intra = m["with_dtc"], m["d1"], m["d2"], m["d3"], m["intra"]

    print(f"\n{'=' * 70}")

    def line(key, message, value, ok, target):
        print(f"  [{'PASA ' if ok else 'FALLA'}] {key} {message}")
        print(f"          medido {value} · diana {target}")
        return ok

    ok1 = line(
        "D1",
        "pacientes con DTC al diagnóstico",
        f"{d1:.0%} ({with_dtc}/{len(res)})",
        0.25 <= d1 <= 0.50,
        "30-40 % (Braun 2005)",
    )
    ok2 = line(
        "D2",
        "fracción de DTC quiescentes (agregada)",
        f"{d2:.1%} de {m['dtc']} DTC (mediana por paciente {m['d2_median']:.0%}, sin veredicto)",
        d2 > 0.80,
        ">80 % (Pantel 2004)",
    )
    ok3 = line(
        "D3",
        "mueren en el trayecto",
        f"{d3:.1%} de {total_intra} intravasaciones",
        d3 > 0.90,
        "la gran mayoría",
    )
    print("\n  [ -- ] D4 control de dependencia: ejecutar aparte (--d4) antes de dar la fase por cerrada")
    print(
        f"\n  Detección: IHC ve {statistics.median(r['biopsy_ihc'] for r in res):.0f} "
        f"y citometría {statistics.median(r['biopsy_cytometry'] for r in res):.0f} "
        f"de {statistics.median(r['dtc'] for r in res):.0f} reales "
        f"(el 24 % y el 47 % medidos en clínica)"
    )
    print(f"\n{'FASE 1 LISTA' if (ok1 and ok2 and ok3) else 'FASE 1 NO LISTA: no se construye encima'}")


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=None, help="por defecto 8 (validación) o 32 (--d4)")
    ap.add_argument("--days", type=float, default=250, help="momento del 'diagnóstico' en el que se mide")
    ap.add_argument("--d4", action="store_true", help="control de dependencia de los parámetros SUPUESTO")
    ap.add_argument("--processes", type=int, default=None, help="procesos en paralelo (1 = en serie)")
    a = ap.parse_args(argv)
    n = a.seeds or (32 if a.d4 else 8)
    if a.d4:
        control_d4(tuple(range(1, n + 1)), a.days, a.processes)
    else:
        validate(tuple(range(1, n + 1)), a.days, a.processes)


if __name__ == "__main__":
    main()


# --------------------------------------------------------------- NOTA DE CALIBRACIÓN
# [HISTÓRICO, v1] prob_brote_vascular ya no existe (2026-10-28): el despertar sale
# del recambio vascular del nicho. prob_intravasation y prob_early_dissemination
# siguen siendo los elegidos aquí; NO se han recalibrado para la v2.
#
# Medido el 2026-10-23, 10 semillas, día 250. Barrido de prob_intravasation,
# prob_brote_vascular y prob_early_dissemination contra D1 y D2:
#
#   intravasación   brote   temprana | D1 % con DTC | D2 % dormidas
#       4e-03       8e-04    2e-05   |     80 %     |     76 %
#       1e-03       8e-04    2e-05   |     70 %     |    100 %
#    >> 3e-04       8e-04    5e-06   |     40 %     |     83 %   <- ELEGIDA
#       1e-04       3e-04    2e-06   |     30 %     |    100 %
#       3e-05       3e-04    1e-06   |     10 %     |    100 %
#       1e-05       1e-04    5e-07   |      0 %     |      0 %
#
# POR QUÉ LA ELEGIDA Y NO LA DE 1e-04, que también pasa los dos criterios: aquella
# da D2 = 100 %, es decir NINGUNA célula despierta nunca. Un modelo así no puede
# producir recurrencias, y la fase 2 consiste exactamente en reproducir las recaídas
# repartidas entre los años 5 y 20. El ajuste elegido deja un 17 % de células
# despiertas, que es el material con el que la fase 2 tiene que trabajar.
#
# ESTO ES CALIBRACIÓN, NO VALIDACIÓN, Y LA DIFERENCIA IMPORTA. Un parámetro ajustado
# contra un dato no predice ese dato: lo reproduce por construcción. Que el modelo dé
# ahora un 40 % de pacientes con DTC NO es un acierto — es el ajuste funcionando.
# La validación, si llega, vendrá de lo que el modelo haga DESPUÉS con estos números
# fijos: la curva de recurrencias (fase 2) y la asimetría célula T / ADC contra células
# dormidas (fase 3). Compárese con T1 en validation.py, donde la geometría del oxígeno
# ACERTÓ los 150-180 µm de Thomlinson & Gray sin que nadie la ajustara: eso sí es
# validación.
