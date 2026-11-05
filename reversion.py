"""¿Cuánto tiene que durar la reprogramación para que la estrategia siga valiendo?

La guardiana —reina curada y reprogramada que vigila el tejido— es la pieza INVENTADA
de este proyecto: no existe hoy en ningún laboratorio (ver HOJA_DE_RUTA.md, fase B). Y
en el modelo tenía dos propiedades que nadie ha justificado: dura para siempre y se
hereda al dividirse.

La literatura de terapia de diferenciación dice lo contrario. En la leucemia
promielocítica (APL) el ATRA impone un estado que se mantiene y cura a >90 %. Pero en
tumores SÓLIDOS la diferenciación inducida REVIERTE: al retirar el inductor las células
recuperan su fenotipo proliferativo, y esa es una de las razones por las que ningún
agente de diferenciación ha curado un tumor sólido.

Así que la pregunta no es "¿funciona la guardiana?", sino "¿cuánto tiene que durar?".
Eso sí es contestable en un laboratorio, y es barato comparado con todo lo demás.

    python reversion.py --semillas 3
"""
from __future__ import annotations

import argparse
import time

from biologia import ESCENARIOS, PASOS_POR_DIA, Parametros
from tejido import Simulacion

# (etiqueta, probabilidad de revertir por paso) — la vida media sale de 1/p pasos.
NIVELES = (
    ("permanente (lo que asumíamos)", 0.0),
    ("vida media ~500 días", 0.0005),
    ("vida media ~125 días", 0.002),
    ("vida media ~25 días", 0.01),
    ("vida media ~5 días", 0.05),
)


def correr(terapia: str, dosis: int | None, escenario: str, semilla: int, dias: float,
           dia_terapia: float, reversion: float, herencia: bool, neoantigenos: float) -> dict:
    p = Parametros(prob_neoantigeno=neoantigenos, reversion_guardiana=reversion,
                   herencia_guardiana=herencia)
    sim = Simulacion(escenario, semilla, p)
    for paso in range(int(dias * PASOS_POR_DIA)):
        if paso == int(dia_terapia * PASOS_POR_DIA):
            sim.liberar_terapia(terapia, dosis)
        sim.paso()
    return sim.resumen()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", type=float, default=450)
    ap.add_argument("--dia-terapia", type=float, default=250)
    ap.add_argument("--semillas", type=int, default=3)
    ap.add_argument("--neoantigenos", type=float, default=0.1)
    ap.add_argument("--terapia", default="troyano+guardiana_nicho+memoria")
    ap.add_argument("--dosis", type=int, default=300)
    a = ap.parse_args()

    print(f"REVERSIÓN · {a.terapia} (dosis {a.dosis}) · {a.escenario} · {a.dias:g} días · "
          f"{a.semillas} semillas")
    print("La guardiana es la pieza inventada del proyecto. Esto mide cuánto tiene que durar.\n")
    print(f"{'duración del estado':<32}{'herencia':>10}{'peligrosas':>12}{'reinas':>8}"
          f"{'revertidas':>12}   por semilla")
    for herencia in (True, False):
        for etiqueta, rev in NIVELES:
            t0 = time.time()
            res = [correr(a.terapia, a.dosis, a.escenario, s, a.dias, a.dia_terapia,
                          rev, herencia, a.neoantigenos) for s in range(a.semillas)]
            m = {k: sum(r[k] for r in res) / len(res)
                 for k in ("peligrosas", "reinas", "guardianas_revertidas")}
            pel = sorted(r["peligrosas"] for r in res)
            print(f"{etiqueta:<32}{'sí' if herencia else 'NO':>10}{m['peligrosas']:>12.0f}"
                  f"{m['reinas']:>8.0f}{m['guardianas_revertidas']:>12.0f}"
                  f"   {pel}   ({time.time() - t0:.0f}s)", flush=True)
        print()


if __name__ == "__main__":
    main()
