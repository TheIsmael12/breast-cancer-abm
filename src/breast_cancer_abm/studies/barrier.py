"""¿Se puede llegar al interior de un tumor sólido? Se mide ANTES de probar ninguna terapia.

Un tumor no está envuelto en una membrana. Lo que impide entrar son tres cosas distintas,
y conviene no confundirlas porque cada una se levanta de una manera:

  1. GEOMETRÍA      una esfera esconde su interior: la fracción accesible cae como 1/r.
                    No se arregla con fármacos; se arregla entrando por dentro (aguja).
  2. PRESIÓN        vasos porosos + linfáticos colapsados = el líquido empuja hacia fuera.
                    Es corriente en contra, no una pared. Baja al degradar la matriz o al
                    normalizar los vasos (Jain 1987-2014).
  3. OXÍGENO        el O2 difunde ~100-150 µm desde un capilar: el centro está lejos, no
                    tapiado. Y la hipoxia no solo apaga a la terapia: activa la
                    desdiferenciación del tumor (HIF-1α), que es lo que fabrica reinas.

Este banco mide las tres por separado, con una terapia inerte que solo camina, para que
lo que se vea sea el ACCESO y no la eficacia.

    bcabm study barrier --days 300 --seeds 2
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from breast_cancer_abm.model.biology import SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import LUMINAL
from breast_cancer_abm.model.tissue3d import Simulation3D

# (etiqueta, vía de entrega, fracción de matriz degradada, modo de vasos, sin presión)
CONDITIONS = (
    ("nada (vena, tal cual)", "vein", 0.0, "normal", False),
    ("degradar matriz 50 %", "vein", 0.5, "normal", False),
    ("normalizar vasos", "vein", 0.0, "normalize", False),
    ("matriz + normalizar", "vein", 0.5, "normalize", False),
    ("reconducir O2 al tumor", "vein", 0.0, "reroute", False),
    ("aguja (una punción)", "needle", 0.0, "normal", False),
    ("intratumoral repartido", "intratumoral", 0.0, "normal", False),
    ("intratumoral + matriz + vasos", "intratumoral", 0.5, "normalize", False),
    ("SIN presión (control irreal)", "vein", 0.0, "normal", True),
)


def _tumor(sim: Simulation3D) -> list:
    return [q for q, c in sim.cells.items() if c.kind == LUMINAL and c.n_drivers and not c.necrosis]


def _radii(sim: Simulation3D, tumor: list):
    """Distancia de cada célula tumoral al centroide de la masa."""
    a = np.array(tumor, float)
    return np.linalg.norm(a - a.mean(0), axis=1)


def measure(sim: Simulation3D, reached: set) -> dict:
    tumor = _tumor(sim)
    if not tumor:
        return dict(tumor_cells=0, o2_rim=0.0, o2_center=0.0, accessible=0.0, was_reached=0.0, is_hypoxic=0.0)
    r = _radii(sim, tumor)
    order = np.argsort(r)
    n = max(1, len(tumor) // 5)
    o2 = np.array([sim.O2[q] for q in tumor])
    core = {tumor[i] for i in order[:n]}
    return dict(
        tumor_cells=len(tumor),
        o2_rim=float(o2[order[-n:]].mean()),
        o2_center=float(o2[order[:n]].mean()),
        accessible=sim.accessible_fraction(),
        # de las células del NÚCLEO (el quinto más interior), ¿a cuántas llegó la terapia?
        was_reached=len(core & reached) / len(core),
        is_hypoxic=float((o2 < sim.p.o2_hypoxia).mean()),
    )


def run(
    scenario: str,
    seed: int,
    days: float,
    therapy_day: float,
    route: str,
    matrix: float,
    vessels: str,
    no_pressure: bool,
    neoantigens: float,
) -> dict:
    p = Parameters(prob_neoantigen=neoantigens, interstitial_pressure=0.0 if no_pressure else 1.0)
    sim = Simulation3D(scenario, seed, p)
    therapy_step = int(therapy_day * STEPS_PER_DAY)
    reached: set = set()
    for step in range(int(days * STEPS_PER_DAY)):
        if step == therapy_step:
            if matrix:
                sim.degrade_matrix(matrix)
            if vessels != "normal":
                sim.treat_vessels(vessels)
            # Sonda inerte: camina y nada más. Mide acceso, no eficacia.
            sim.release_therapy("repairer", 40, route)
        sim.step()
        if step > therapy_step and step % 4 == 0:
            for ag in sim.agents:
                if ag.kind != "AI":
                    continue
                z, y, x = ag.pos
                for dz in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dx in (-1, 0, 1):
                            q = (z + dz, y + dy, x + dx)
                            c = sim.cells.get(q)
                            if c is not None and c.kind == LUMINAL and c.n_drivers:
                                reached.add(q)
    r = measure(sim, reached)
    r["blocked_count"] = sim.stats["steps_blocked"]
    return r


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--days", type=float, default=300)
    ap.add_argument("--therapy-day", type=float, default=200)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--neoantigens", type=float, default=0.1)
    a = ap.parse_args(argv)

    print(
        f"3D · {a.scenario} · {a.days:g} días · acceso medido desde el día {a.therapy_day:g}"
        f" · {a.seeds} semillas"
    )
    print(
        "núcleo alcanzado = del quinto MÁS INTERIOR del tumor, qué fracción llegó a tener "
        "una célula terapéutica al lado\n"
    )
    header = (
        f"{'condición':<32}{'tumorales':>10}{'O2 borde':>10}{'O2 centro':>11}"
        f"{'accesible':>11}{'núcleo alcanz.':>16}{'bloqueados':>12}"
    )
    print(header)
    for tag, route, matrix, vessels, no_pressure in CONDITIONS:
        t0 = time.time()
        res = [
            run(a.scenario, s, a.days, a.therapy_day, route, matrix, vessels, no_pressure, a.neoantigens)
            for s in range(a.seeds)
        ]
        m = {k: sum(r[k] for r in res) / len(res) for k in res[0]}
        print(
            f"{tag:<32}{m['tumor_cells']:>10.0f}{m['o2_rim']:>10.2f}"
            f"{m['o2_center']:>11.2f}{m['accessible']:>10.0%}{m['was_reached']:>16.0%}"
            f"{m['blocked_count']:>12.0f}   ({time.time() - t0:.0f}s)"
        )


if __name__ == "__main__":
    main()
