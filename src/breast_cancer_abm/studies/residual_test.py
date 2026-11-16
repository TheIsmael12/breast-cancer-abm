"""¿Quién tiene todavía tumor? La prueba que decide si se puede NO tocar a la paciente.

no_surgery.py (2026-11-05) mostró que la vía Kuerer funciona en el modelo y que su
punto débil es la biopsia de comprobación. Aquí se miden las dos llaves posibles:

  A) biopsia por vacío con más o menos cilindros (lo que hacen OPTIMIST, Kuerer)
  B) biopsia líquida: ADN tumoral circulante (ctDNA) con ensayo tumor-informado

FASE 'foto': todas las pacientes reciben el arsenal el día 120 y se miran el día 180,
60 días después. Para cada una se repite cada prueba muchas veces (sin avanzar el
tiempo) y se obtiene la probabilidad de detectar el carcinoma que de verdad queda.

CALIBRACIÓN DEL ctDNA (lo único que no sale del modelo): los fragmentos por célula se
fijan para que, ANTES del tratamiento (día 120), la sensibilidad sea la clínica de un
ensayo tumor-informado en estadio precoz (~75 %; menor en luminal pequeño). Se dan
también 50 % y 90 % para ver si la conclusión depende de ese número.

    bcabm study residual-test --seeds 24 --processes 6
"""

from __future__ import annotations

import argparse
import math
import statistics
from multiprocessing import Pool

from breast_cancer_abm.model.biology import STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation

ARSENAL = "trojan+niche_guardian+memory"
CORES = (6, 12, 24, 48, 96)
REPEATS = 30


def snapshot(seed: int) -> dict:
    sim = Simulation("luminal", seed, Parameters(prob_neoantigen=0.1, editing_error=0.01))
    for _ in range(int(120 * STEPS_PER_DAY)):
        sim.step()
    base_pre = sim.ctdna(1.0)["expected"]  # λ por unidad de fragmentos/célula
    sim.release_therapy(ARSENAL, 300)
    for _ in range(int(60 * STEPS_PER_DAY)):
        sim.step()
    out = dict(seed=seed, base_pre=base_pre, base_post=sim.ctdna(1.0)["expected"])
    out["real"] = sim.biopsy(1)["real"]
    for k in CORES:
        out[f"b{k}"] = sum(sim.biopsy(k)["positive"] for _ in range(REPEATS)) / REPEATS
    return out


def calibrate(bases: list[float], goal: float) -> float:
    """Fragmentos/célula tal que la sensibilidad media antes de tratar sea `objetivo`."""
    lo, hi = 1e-7, 10.0
    for _ in range(80):
        f = math.sqrt(lo * hi)
        s = statistics.mean(1 - math.exp(-f * b) for b in bases)
        lo, hi = (f, hi) if s < goal else (lo, f)
    return f


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=24)
    ap.add_argument("--processes", type=int, default=6)
    a = ap.parse_args(argv)
    with Pool(a.processes, maxtasksperchild=1) as pool:
        fs = pool.map(snapshot, range(1, a.seeds + 1), chunksize=1)

    with_ = [f for f in fs if f["real"] > 0]
    without = [f for f in fs if f["real"] == 0]
    print(f"FOTO DEL DÍA 180 (60 d tras el arsenal) · {len(fs)} pacientes")
    print(
        f"  con carcinoma residual: {len(with_)}  ·  sin carcinoma (respuesta completa real): {len(without)}"
    )
    print(f"  carcinoma residual en las que tienen: {sorted(f['real'] for f in with_)}\n")

    print("BIOPSIA POR VACÍO — falsos negativos (sobre las que SÍ tienen carcinoma)")
    for k in CORES:
        fnr = statistics.mean(1 - f[f"b{k}"] for f in with_) if with_ else float("nan")
        print(f"  {k:>3} cilindros: falsos negativos {fnr:6.1%}")

    print("\nBIOPSIA LÍQUIDA (ctDNA) — una sola extracción el día 180")
    bases = [f["base_pre"] for f in fs]
    for goal in (0.50, 0.75, 0.90):
        fr = calibrate(bases, goal)
        fnr = statistics.mean(math.exp(-fr * f["base_post"]) for f in with_) if with_ else float("nan")
        fp = statistics.mean(1 - math.exp(-fr * f["base_post"]) for f in without) if without else float("nan")
        print(
            f"  sensibilidad pre-tratamiento {goal:.0%} → {fr:.2e} fragmentos/célula · "
            f"falsos negativos el día 180: {fnr:6.1%} · positivos sin carcinoma: {fp:5.1%}"
        )
    print("\n  por paciente (carcinoma residual · ADN que vierte, relativo al de antes de tratar):")
    print("  " + "  ".join(f"{f['real']}·{f['base_post'] / max(f['base_pre'], 1):.2%}" for f in with_))


if __name__ == "__main__":
    main()
