"""¿Quién tiene todavía tumor? La prueba que decide si se puede NO tocar a la paciente.

sin_cirugia.py (2026-11-05) mostró que la vía Kuerer funciona en el modelo y que su
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

    python prueba_residuo.py --semillas 24 --procesos 6
"""
from __future__ import annotations

import argparse
import math
import statistics
from multiprocessing import Pool

from biologia import PASOS_POR_DIA, Parametros
from tejido import Simulacion

ARSENAL = "troyano+guardiana_nicho+memoria"
CILINDROS = (6, 12, 24, 48, 96)
REPETICIONES = 30


def foto(semilla: int) -> dict:
    sim = Simulacion("luminal", semilla, Parametros(prob_neoantigeno=0.1, error_edicion=0.01))
    for _ in range(int(120 * PASOS_POR_DIA)):
        sim.paso()
    base_pre = sim.ctdna(1.0)["esperado"]          # λ por unidad de fragmentos/célula
    sim.liberar_terapia(ARSENAL, 300)
    for _ in range(int(60 * PASOS_POR_DIA)):
        sim.paso()
    out = dict(semilla=semilla, base_pre=base_pre, base_post=sim.ctdna(1.0)["esperado"])
    out["real"] = sim.biopsia(1)["real"]
    for k in CILINDROS:
        out[f"b{k}"] = sum(sim.biopsia(k)["positiva"] for _ in range(REPETICIONES)) / REPETICIONES
    return out


def calibrar(bases: list[float], objetivo: float) -> float:
    """Fragmentos/célula tal que la sensibilidad media antes de tratar sea `objetivo`."""
    lo, hi = 1e-7, 10.0
    for _ in range(80):
        f = math.sqrt(lo * hi)
        s = statistics.mean(1 - math.exp(-f * b) for b in bases)
        lo, hi = (f, hi) if s < objetivo else (lo, f)
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--semillas", type=int, default=24)
    ap.add_argument("--procesos", type=int, default=6)
    a = ap.parse_args()
    with Pool(a.procesos, maxtasksperchild=1) as pool:
        fs = pool.map(foto, range(1, a.semillas + 1), chunksize=1)

    con = [f for f in fs if f["real"] > 0]
    sin = [f for f in fs if f["real"] == 0]
    print(f"FOTO DEL DÍA 180 (60 d tras el arsenal) · {len(fs)} pacientes")
    print(f"  con carcinoma residual: {len(con)}  ·  sin carcinoma (respuesta completa real): "
          f"{len(sin)}")
    print(f"  carcinoma residual en las que tienen: {sorted(f['real'] for f in con)}\n")

    print("BIOPSIA POR VACÍO — falsos negativos (sobre las que SÍ tienen carcinoma)")
    for k in CILINDROS:
        fnr = statistics.mean(1 - f[f"b{k}"] for f in con) if con else float("nan")
        print(f"  {k:>3} cilindros: falsos negativos {fnr:6.1%}")

    print("\nBIOPSIA LÍQUIDA (ctDNA) — una sola extracción el día 180")
    bases = [f["base_pre"] for f in fs]
    for objetivo in (0.50, 0.75, 0.90):
        fr = calibrar(bases, objetivo)
        fnr = statistics.mean(math.exp(-fr * f["base_post"]) for f in con) if con else float("nan")
        fp = statistics.mean(1 - math.exp(-fr * f["base_post"]) for f in sin) if sin else float("nan")
        print(f"  sensibilidad pre-tratamiento {objetivo:.0%} → {fr:.2e} fragmentos/célula · "
              f"falsos negativos el día 180: {fnr:6.1%} · positivos sin carcinoma: {fp:5.1%}")
    print("\n  por paciente (carcinoma residual · ADN que vierte, relativo al de antes de tratar):")
    print("  " + "  ".join(f"{f['real']}·{f['base_post'] / max(f['base_pre'], 1):.2%}" for f in con))


if __name__ == "__main__":
    main()
