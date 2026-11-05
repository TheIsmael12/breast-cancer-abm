"""Banco de pruebas sin ventana: compara terapias con varias semillas.

    python experimento.py --escenario luminal --dias 200 --dia-terapia 120 --semillas 3
    python experimento.py --escenario her2 --plasticidad 0   # colmena sin plasticidad
    python experimento.py --terapias guardiana+memoria --vigilancia --hormonal
"""
from __future__ import annotations

import argparse
import time

from biologia import CALIBRACIONES, ESCENARIOS, PASOS_POR_DIA, Parametros
from politicas import TERAPIAS
from tejido import Simulacion

COLUMNAS = ("tumorales", "peligrosas", "reinas", "senescentes",
            "comidas_macrofago", "comidas_macrofago_sanas", "cd47_sanas", "colateral")


def correr(escenario: str, terapia: str | None, semilla: int, dias: float, dia_terapia: float,
           params: Parametros, vasos: str = "normal", vigilancia: bool = False,
           hormonal: bool = False, castrar: bool = False, tiodg: bool = False,
           dia_tiodg: float | None = None, anti_cd47: bool = False,
           dia_cd47: float | None = None, dia_castrar: float | None = None) -> dict:
    sim = Simulacion(escenario, semilla, params)
    paso_terapia = int(dia_terapia * PASOS_POR_DIA)
    for paso in range(int(dias * PASOS_POR_DIA)):
        if terapia and paso == paso_terapia:
            sim.liberar_terapia(terapia)
        if vasos != "normal" and paso == paso_terapia:
            sim.tratar_vasos(vasos)
        if anti_cd47 and paso == int((dia_cd47 if dia_cd47 is not None else dia_terapia)
                                     * PASOS_POR_DIA):
            sim.bloquear_cd47()
        if castrar and paso == int((dia_castrar if dia_castrar is not None else dia_terapia)
                                   * PASOS_POR_DIA):
            sim.inhibir_telomerasa()
        if tiodg and paso == int((dia_tiodg if dia_tiodg is not None else dia_terapia)
                                 * PASOS_POR_DIA):
            sim.danar_telomeros()
        if hormonal and paso == paso_terapia:
            sim.p.estrogeno = 0.1   # tamoxifeno / inhibidor de aromatasa
            sim._evento("PREVENCIÓN HORMONAL: estrógenos bloqueados")
        if vigilancia and paso > paso_terapia and paso % (30 * PASOS_POR_DIA) == 0:
            # Biopsia líquida mensual: si detecta ADN tumoral (>= 5 células peligrosas), dosis pequeña.
            if sim.resumen()["peligrosas"] >= 5:
                sim.liberar_terapia("guardiana", 10)
        sim.paso()
    return sim.resumen()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", type=float, default=200)
    ap.add_argument("--dia-terapia", type=float, default=120)
    ap.add_argument("--semillas", type=int, default=3)
    ap.add_argument("--terapias", nargs="*", default=["ninguna", *TERAPIAS],
                    help="admite combinaciones con '+', p. ej. guardiana+memoria")
    ap.add_argument("--plasticidad", type=float, default=Parametros.plasticidad,
                    help="probabilidad por paso de que una obrera tumoral se vuelva reina")
    ap.add_argument("--vasos", default="normal",
                    choices=("normal", "antiangiogenico", "normalizar", "asfixiar", "reconducir"))
    ap.add_argument("--vigilancia", action="store_true", help="biopsia líquida mensual + dosis si positiva")
    ap.add_argument("--hormonal", action="store_true", help="bloqueo de estrógenos desde el día de terapia")
    ap.add_argument("--castrar", action="store_true",
                    help="inhibidor de telomerasa el día de la terapia (bloquea el contador de divisiones)")
    ap.add_argument("--anti-cd47", action="store_true",
                    help="bloquea CD47 ('no me comas'); afecta también a las células sanas")
    ap.add_argument("--dia-cd47", type=float, default=None,
                    help="día del anti-CD47; por defecto, el mismo que la terapia")
    ap.add_argument("--dia-castrar", type=float, default=None,
                    help="día del inhibidor de telomerasa; por defecto, el mismo que la terapia")
    ap.add_argument("--tiodg", action="store_true",
                    help="6-tio-dG: rompe el telómero al instante en células con telomerasa")
    ap.add_argument("--farmaco-p53", type=float, default=0.0,
                    help="cualidad: el fármaco también restaura p53 (0-1)")
    ap.add_argument("--farmaco-cebo", action="store_true",
                    help="cualidad: la carga se queda dentro; quien se coma la célula se envenena")
    ap.add_argument("--farmaco-anti-sasp", type=float, default=0.0,
                    help="cualidad: silencia el SASP (0-1), como un inhibidor de JAK")
    ap.add_argument("--farmaco-persistente", type=float, default=0.0,
                    help="cualidad: sigue alcanzando a las células que nacen después (0-1)")
    ap.add_argument("--dia-tiodg", type=float, default=None,
                    help="día del 6-tio-dG; por defecto, el mismo que la terapia celular")
    ap.add_argument("--neoantigenos", type=float, default=Parametros.prob_neoantigeno,
                    help="1 = optimista; 0.1 = calibrado con pacientes reales (TCGA)")
    ap.add_argument("--calibracion", choices=list(CALIBRACIONES), default="tratamiento",
                    help="tasa de mutación. 'tratamiento' (por defecto) da tumores siempre "
                         "y rápido; las otras bajan la incidencia y hacen estudiable la "
                         "prevención, a costa de necesitar más semillas y más días. "
                         "Ver la tabla en biologia.py")
    a = ap.parse_args()
    params = Parametros(plasticidad=a.plasticidad, prob_neoantigeno=a.neoantigenos,
                        tasa_mutacion=CALIBRACIONES[a.calibracion],
                        farmaco_p53=a.farmaco_p53, farmaco_cebo=a.farmaco_cebo,
                        farmaco_anti_sasp=a.farmaco_anti_sasp,
                        farmaco_persistente=a.farmaco_persistente)

    print(f"Escenario {a.escenario} · {a.dias:g} días · terapia el día {a.dia_terapia:g} · "
          f"{a.semillas} semillas (media) · plasticidad {a.plasticidad:g} · neoantígenos {a.neoantigenos:g} · "
          f"vasos {a.vasos}"
          f"{' · vigilancia' if a.vigilancia else ''}{' · hormonal' if a.hormonal else ''}"
          f"{' · ANTI-CD47' if a.anti_cd47 else ''}"
          f"{' · ANTITELOMERASA' if a.castrar else ''}{f' · 6-TIO-dG dia {a.dia_tiodg or a.dia_terapia:g}' if a.tiodg else ''}")
    print("peligrosas = in situ + invasivas + reinas (lo que de verdad es cáncer)\n")
    print(f"{'terapia':<20}" + "".join(f"{c[:12]:>13}" for c in COLUMNAS))
    for terapia in a.terapias:
        t0 = time.time()
        res = [correr(a.escenario, None if terapia == "ninguna" else terapia, s, a.dias,
                      a.dia_terapia, params, a.vasos, a.vigilancia, a.hormonal, a.castrar,
                      a.tiodg, a.dia_tiodg, a.anti_cd47, a.dia_cd47, a.dia_castrar)
               for s in range(a.semillas)]
        medias = {c: sum(r.get(c, 0) for r in res) / len(res) for c in COLUMNAS}
        print(f"{terapia:<20}" + "".join(f"{medias[c]:>13.1f}" for c in COLUMNAS)
              + f"   ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
