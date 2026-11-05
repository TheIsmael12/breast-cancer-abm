"""Simulador de cáncer de mama y células terapéuticas programables.

    python main.py                                   # escenario luminal ER+
    python main.py --escenario brca1                 # portadora de BRCA1
    python main.py --escenario her2 --terapia hibrido --dia-terapia 150
"""
from __future__ import annotations

import argparse

from biologia import ESCENARIOS
from politicas import TERAPIAS
from visor import Visor


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--semilla", type=int, default=None, help="fija el azar para repetir una partida")
    ap.add_argument("--terapia", default=None,
                    help=f"se libera sola en --dia-terapia; admite '+': {', '.join(TERAPIAS)}")
    ap.add_argument("--dia-terapia", type=float, default=150)
    a = ap.parse_args()
    if a.terapia and not all(t in TERAPIAS for t in a.terapia.split("+")):
        ap.error(f"terapia desconocida: {a.terapia}")
    Visor(a.escenario, a.semilla, a.terapia, a.dia_terapia).mostrar()


if __name__ == "__main__":
    main()
