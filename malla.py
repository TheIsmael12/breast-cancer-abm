"""Las cinco maneras de atacar la malla: ¿cuál abre el paso sin soltar al tumor?

La malla desmoplásica hace DOS cosas a la vez, y son opuestas para nosotros:

  ESTORBA   genera la presión intersticial que impide entrar a cualquier terapia.
  CONTIENE  frena la invasión. Es, en parte, lo que mantiene al tumor donde está.

Por eso en la clínica todo lo que la DESTRUYE ha fracasado (PEGPH20 en fase III, 2019;
eliminar CAFs empeoró los tumores en ratones) y lo que la ABLANDA o normaliza funciona
(losartán con quimio, Murphy 2019; normalización vascular, Jain).

Este banco mide las dos caras por separado, para cada intervención:

  acceso     cuántas fagocitosis consiguen los macrófagos (que ya están ahí de serie)
             y cuántos pasos se bloquean contra la corriente de salida
  contención cuántas células invasivas hay al final y cuántas invasiones frenó la malla

Una intervención solo sirve si sube el acceso SIN disparar las invasivas.

    python malla.py --semillas 3 --dias 450
"""
from __future__ import annotations

import argparse
import time

from biologia import ESCENARIOS, PASOS_POR_DIA, Parametros
from tejido import Simulacion

MODOS = ("ninguno", "hialuronidasa", "anti_lox", "losartan", "anti_tgfb", "matar_caf")

EXPLICACION = {
    "ninguno": "sin tocar la malla",
    "hialuronidasa": "degrada el HA (PEGPH20): quita material, los CAFs siguen",
    "anti_lox": "destensa sin retirar (simtuzumab)",
    "losartan": "reduce la producción (Murphy 2019)",
    "anti_tgfb": "corta la activación de nuevos CAFs",
    "matar_caf": "elimina al fabricante (lo que empeoró en ratones)",
}


def correr(modo: str, escenario: str, semilla: int, dias: float, dia: float,
           params: Parametros, intensidad: float, cadena: bool) -> dict:
    sim = Simulacion(escenario, semilla, params)
    paso_0 = int(dia * PASOS_POR_DIA)
    for paso in range(int(dias * PASOS_POR_DIA)):
        if paso == paso_0 and modo != "ninguno":
            sim.tratar_matriz(modo, intensidad)
        if cadena:
            if paso == paso_0 + int(5 * PASOS_POR_DIA):
                sim.bloquear_cd47(dirigido=True)
            elif paso == paso_0 + int(15 * PASOS_POR_DIA):
                sim.danar_telomeros()
            elif paso == paso_0 + int(20 * PASOS_POR_DIA):
                sim.liberar_terapia("guardiana_nicho+memoria")
        sim.paso()
    r = sim.resumen()
    r["malla"] = float(sim.Ms.max())
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", type=float, default=450)
    ap.add_argument("--dia", type=float, default=250, help="día en que se trata la malla")
    ap.add_argument("--semillas", type=int, default=3)
    ap.add_argument("--neoantigenos", type=float, default=0.1)
    ap.add_argument("--intensidad", type=float, default=0.6)
    ap.add_argument("--cadena", action="store_true",
                    help="además de la malla, aplica la cadena completa")
    ap.add_argument("--modos", nargs="*", default=list(MODOS))
    a = ap.parse_args()
    params = Parametros(prob_neoantigeno=a.neoantigenos)

    print(f"MALLA · {a.escenario} · {a.dias:g} días · intervención el día {a.dia:g} al "
          f"{a.intensidad:.0%} · {a.semillas} semillas{' · + cadena completa' if a.cadena else ''}")
    print("sirve solo si sube ACCESO sin disparar INVASIVAS\n")
    print(f"{'intervención':<16}{'peligros':>10}{'invasivas':>11}{'comidas':>9}"
          f"{'bloqueados':>12}{'frenadas':>10}{'malla':>8}   qué es")
    for modo in a.modos:
        t0 = time.time()
        res = [correr(modo, a.escenario, s, a.dias, a.dia, params, a.intensidad, a.cadena)
               for s in range(a.semillas)]
        m = {k: sum(r.get(k, 0) for r in res) / len(res) for k in
             ("peligrosas", "invasivas", "comidas_macrofago", "pasos_bloqueados",
              "invasiones_frenadas", "malla")}
        print(f"{modo:<16}{m['peligrosas']:>10.0f}{m['invasivas']:>11.0f}"
              f"{m['comidas_macrofago']:>9.0f}{m['pasos_bloqueados']:>12.0f}"
              f"{m['invasiones_frenadas']:>10.0f}{m['malla']:>8.2f}   {EXPLICACION[modo]}")
        peligro = sorted(r["peligrosas"] for r in res)
        print(f"{'':<16}mediana {peligro[len(peligro) // 2]:>7} · {peligro}"
              f"   ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
