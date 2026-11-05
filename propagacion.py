"""Corrección que se propaga: ¿cuánta fidelidad necesita para merecer la pena?

El cuello de botella medido en este proyecto es que cada célula corregida cuesta un
agente que la encuentre, la secuencie y acierte con todos sus genes rotos. En 450 días
solo se consiguen ~23. Es coste lineal.

La idea: que una célula ya corregida corrija a sus vecinas, y que las que quedan limpias
corrijan a su vez. Pasa a crecimiento exponencial, y basta conseguirlo una vez.

Equivalentes reales a medias:
  - EFECTO ESPECTADOR: existe en terapia génica. El gen suicida HSV-TK pasa el tóxico a
    las vecinas por las uniones comunicantes, sin necesidad de transducirlas todas.
  - GENE DRIVE: un elemento de CRISPR que se copia a sí mismo en la población. Funciona
    en insectos; por eso mismo está sujeto a moratorias.

Y el peligro, que es la razón de ser de este banco: si algo se copia solo, **los errores
también se copian**. La edición génica tiene errores fuera de diana. Una corrección que
se auto-propaga con un 1 % de error va sembrando mutaciones nuevas justo donde más
células hay. Por eso aquí no se pregunta "¿funciona?", sino:

    ¿qué fidelidad MÍNIMA necesita la edición para que propagarla sea mejor que no hacerlo?

Eso sí es una especificación que un laboratorio puede medir.

    python propagacion.py --semillas 3
"""
from __future__ import annotations

import argparse
import time

from biologia import ESCENARIOS, PASOS_POR_DIA, Parametros
from tejido import Simulacion

# (propagación por paso, error fuera de diana, alcance en vecinas)
CASOS = (
    (0.00, 0.00, 0),   # control: sin propagación, como hasta ahora
    (0.02, 0.00, 8),
    (0.05, 0.00, 8),
    (0.15, 0.00, 8),
    (0.05, 0.01, 8),   # 1 % de error fuera de diana
    (0.05, 0.05, 8),   # 5 %
    (0.15, 0.05, 8),   # propagación alta Y error alto: el peor caso
)


def correr(prop: float, err: float, alcance: int, escenario: str, semilla: int,
           dias: float, dia_terapia: float, terapia: str, dosis: int,
           neoantigenos: float) -> dict:
    p = Parametros(prob_neoantigeno=neoantigenos, alcance_correccion=alcance,
                   propagacion_correccion=prop, error_edicion=err)
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

    print(f"PROPAGACIÓN · {a.terapia} (dosis {a.dosis}) · {a.escenario} · {a.dias:g} días · "
          f"{a.semillas} semillas")
    print("La pregunta no es si funciona, sino con cuánto error fuera de diana deja de "
          "compensar.\n")
    print(f"{'propag.':>9}{'error':>7}{'alcance':>9}{'peligrosas':>12}{'reinas':>8}"
          f"{'propagadas':>12}{'limpias':>9}{'errores':>9}   por semilla")
    for prop, err, alcance in CASOS:
        t0 = time.time()
        res = [correr(prop, err, alcance, a.escenario, s, a.dias, a.dia_terapia,
                      a.terapia, a.dosis, a.neoantigenos) for s in range(a.semillas)]
        m = {k: sum(r[k] for r in res) / len(res) for k in
             ("peligrosas", "reinas", "correcciones_propagadas", "corregidas_completas",
              "errores_edicion")}
        pel = sorted(r["peligrosas"] for r in res)
        print(f"{prop:>9.2f}{err:>7.2f}{alcance:>9}{m['peligrosas']:>12.0f}{m['reinas']:>8.0f}"
              f"{m['correcciones_propagadas']:>12.0f}{m['corregidas_completas']:>9.0f}"
              f"{m['errores_edicion']:>9.0f}   {pel}   ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
