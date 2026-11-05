"""La cadena completa: abrir · desenmascarar · castrar · que se lo coman · cerrar la fuente.

La idea, en los términos en que la planteó el usuario: si al tumor le quitas la careta
(CD47) y además no puede duplicarse (telómeros), se convierte en objetivo único y el
cuerpo se lo come solo. Lo natural, sin cirugía y sin quimio.

El modelo añade dos pasos que la versión original no tenía, y los dos salieron de medir:

  PASO 0 (abrir)   sin esto no llega nadie. Con presión intersticial los macrófagos
                   fagocitaron 3 células; sin ella, 920. Es el cuello de botella real.
  PASO 4 (cerrar)  comerse 920 células bajó las peligrosas de 957 a 932. Limpiar la masa
                   no cura si las células madre tumorales siguen reponiéndola.

Los brazos de control no están de adorno: sirven para saber si la cadena funciona por lo
que creemos que funciona, o si en realidad la hace sola una de sus piezas.

    python cadena.py --semillas 3 --dias 500
"""
from __future__ import annotations

import argparse
import time

from biologia import ESCENARIOS, PASOS_POR_DIA, Parametros
from tejido import Simulacion

# Cada paso es (días después del inicio del tratamiento, qué se hace).
ABRIR = "abrir"          # hialuronidasa + normalizar vasos: baja la presión intersticial
CARETA = "cd47"          # biespecífico anti-CD47: quita el "no me comas" solo donde hay sospecha
CASTRAR = "6tiodg"       # 6-tio-dG: rompe el telómero YA en las células con telomerasa
FUENTE = "guardiana"     # cierra la fuente: cura las reinas y deja memoria residente

BRAZOS = {
    "nada": (),
    "solo abrir": ((0, ABRIR),),
    "solo guardiana+memoria": ((20, FUENTE),),
    "cadena SIN abrir": ((5, CARETA), (15, CASTRAR), (20, FUENTE)),
    "cadena SIN cerrar": ((0, ABRIR), (5, CARETA), (15, CASTRAR)),
    "CADENA COMPLETA": ((0, ABRIR), (5, CARETA), (15, CASTRAR), (20, FUENTE)),
}

COLUMNAS = (("peligrosas", "peligrosas"), ("reinas", "reinas"),
            ("comidas_macrofago", "comidas"), ("comidas_macrofago_sanas", "comidas SANAS"),
            ("colateral", "colateral"), ("pasos_bloqueados", "bloqueados"))


def correr(brazo, escenario: str, semilla: int, dias: float, dia_terapia: float,
           params: Parametros, puerta: str) -> dict:
    sim = Simulacion(escenario, semilla, params)
    agenda = {int((dia_terapia + d) * PASOS_POR_DIA): paso for d, paso in brazo}
    for paso in range(int(dias * PASOS_POR_DIA)):
        accion = agenda.get(paso)
        if accion == ABRIR:
            sim.degradar_matriz(0.5)
            sim.tratar_vasos("normalizar")
        elif accion == CARETA:
            sim.bloquear_cd47(dirigido=True)
        elif accion == CASTRAR:
            sim.danar_telomeros()
        elif accion == FUENTE:
            sim.liberar_terapia(puerta)
        sim.paso()
    return sim.resumen()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", type=float, default=500)
    ap.add_argument("--dia-terapia", type=float, default=250)
    ap.add_argument("--semillas", type=int, default=3)
    ap.add_argument("--neoantigenos", type=float, default=0.1)
    ap.add_argument("--puerta", default="guardiana_nicho+memoria",
                    help="terapia del paso 'cerrar la fuente'")
    ap.add_argument("--brazos", nargs="*", default=list(BRAZOS))
    a = ap.parse_args()
    params = Parametros(prob_neoantigeno=a.neoantigenos)

    print(f"CADENA · {a.escenario} · {a.dias:g} días · inicio el día {a.dia_terapia:g} · "
          f"{a.semillas} semillas · neoantígenos {a.neoantigenos:g} · fuente: {a.puerta}")
    print("abrir(0) - careta(+5) - castrar(+15) - cerrar la fuente(+20)")
    print("peligrosas = in situ + invasivas + reinas · entre paréntesis, el rango de semillas\n")
    print(f"{'brazo':<24}" + "".join(f"{e:>16}" for _, e in COLUMNAS))
    for nombre in a.brazos:
        t0 = time.time()
        res = [correr(BRAZOS[nombre], a.escenario, s, a.dias, a.dia_terapia, params, a.puerta)
               for s in range(a.semillas)]
        fila = f"{nombre:<24}"
        for clave, _ in COLUMNAS:
            v = [r.get(clave, 0) for r in res]
            fila += f"{sum(v) / len(v):>16.0f}"
        print(fila + f"   ({time.time() - t0:.0f}s)")
        peligro = sorted(r["peligrosas"] for r in res)
        mediana = peligro[len(peligro) // 2]
        print(f"{'':<24}mediana {mediana:>8} · por semilla {peligro}")


if __name__ == "__main__":
    main()
