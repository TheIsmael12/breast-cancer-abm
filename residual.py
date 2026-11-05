"""El escenario real: quitar la masa primero y pelear por el residuo.

Este banco existe por un resultado, no por una idea. Todo lo construido en este proyecto
—el caballo de Troya, la corrección que se propaga, la vigilancia residente, la autoridad
del tejido— actúa sobre la SUPERFICIE del tumor:

  acceso en 3D ................. 26 % accesible, 1 % del núcleo alcanzado
  macrófagos ................... no atraviesan la corriente de salida
  autoridad del tejido ......... 148 células con contexto, tenga el tumor 1.089 o 8.832

La superficie crece con r², el volumen con r³. Por eso existe un umbral (~4.000 células
en este modelo) a partir del cual el tumor gana por geometría, sin necesidad de evolucionar.

La cirugía no es elegante, pero es la única forma conocida de convertir un problema de
volumen en uno de superficie. La pregunta que contesta este banco no es si la cirugía
funciona —eso ya se sabe—, sino si el arsenal de este proyecto sirve para lo que de verdad
importa después: la ENFERMEDAD RESIDUAL y la recaída.

    python residual.py --semillas 3
"""
from __future__ import annotations

import argparse
import time

from biologia import ESCENARIOS, PASOS_POR_DIA, Parametros
from tejido import Simulacion

# (etiqueta, ¿opera?, terapia tras la cirugía, propagación de la corrección)
BRAZOS = (
    ("nada", False, None, 0.0),
    ("solo el arsenal (sin operar)", False, "troyano+guardiana_nicho+memoria", 0.0),
    ("solo cirugía", True, None, 0.0),
    ("cirugía + vigilancia", True, "memoria", 0.0),
    ("cirugía + arsenal", True, "troyano+guardiana_nicho+memoria", 0.0),
    ("cirugía + arsenal + propagación", True, "troyano+guardiana_nicho+memoria", 0.15),
)


def correr(brazo, escenario: str, semilla: int, dias: float, dia: float,
           dosis: int, cobertura: float, neoantigenos: float) -> dict:
    _, opera, terapia, prop = brazo
    p = Parametros(prob_neoantigeno=neoantigenos,
                   alcance_correccion=8 if prop else 0,
                   propagacion_correccion=prop, error_edicion=0.01)
    sim = Simulacion(escenario, semilla, p)
    paso_0 = int(dia * PASOS_POR_DIA)
    for paso in range(int(dias * PASOS_POR_DIA)):
        if paso == paso_0 and opera:
            sim.resecar(cobertura)
        # La terapia va DESPUÉS de operar, como la adyuvancia real.
        if paso == paso_0 + (4 * PASOS_POR_DIA if opera else 0) and terapia:
            sim.liberar_terapia(terapia, dosis)
        sim.paso()
    return sim.resumen()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", type=float, default=500)
    ap.add_argument("--dia", type=float, default=250, help="día de la cirugía")
    ap.add_argument("--semillas", type=int, default=3)
    ap.add_argument("--neoantigenos", type=float, default=0.1)
    ap.add_argument("--dosis", type=int, default=300)
    ap.add_argument("--cobertura", type=float, default=0.95,
                    help="fracción de la masa que se lleva el bloque extirpado")
    a = ap.parse_args()

    print(f"RESIDUAL · {a.escenario} · {a.dias:g} días · cirugía el día {a.dia:g} "
          f"(cobertura {a.cobertura:.0%}) · terapia 4 días después · {a.semillas} semillas")
    print("La pregunta: ¿sirve nuestro arsenal para la recaída, que es el problema real?\n")
    print(f"{'brazo':<34}{'peligrosas':>12}{'reinas':>8}{'extirp.':>9}"
          f"{'sanas ext.':>12}   por semilla")
    for brazo in BRAZOS:
        t0 = time.time()
        res = [correr(brazo, a.escenario, s, a.dias, a.dia, a.dosis, a.cobertura,
                      a.neoantigenos) for s in range(a.semillas)]
        m = {k: sum(r[k] for r in res) / len(res)
             for k in ("peligrosas", "reinas", "extirpadas", "sanas_extirpadas")}
        pel = sorted(r["peligrosas"] for r in res)
        print(f"{brazo[0]:<34}{m['peligrosas']:>12.0f}{m['reinas']:>8.0f}"
              f"{m['extirpadas']:>9.0f}{m['sanas_extirpadas']:>12.0f}"
              f"   {pel}   ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
