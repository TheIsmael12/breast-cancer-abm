"""¿Qué puerta lógica identifica mejor a la célula tumoral y SOLO a ella?

INFORME_IDENTIFICADORES.md concluyó, sobre 1.084 pacientes reales, que ningún marcador
tiene a la vez cobertura y exclusividad, y que la salida son las puertas Y. Esto lo
comprueba DENTRO del simulador: deja crecer un tumor sin tratar y, en ese tejido, mide
para cada puerta dos números que son los de un test diagnóstico:

  SENSIBILIDAD  de las células tumorales, ¿a cuántas ve?        (si falla, el tumor rebrota)
  ESPECIFICIDAD de las células SANAS, ¿a cuántas respeta?       (si falla, hay daño colateral)

No simula tratamiento: solo pregunta a cada puerta qué haría con cada célula del tejido.

    python puertas.py --dias 400 --semillas 3
"""
from __future__ import annotations

import argparse

from biologia import ESCENARIOS, PASOS_POR_DIA, Parametros
from tejido import LUMINAL, Simulacion

PUERTAS = ("simple", "doble", "nicho", "tert", "tert+nicho")
ETIQUETA = {
    "simple": "1 señal (actividad)",
    "doble": "Y: tumoral + actividad",
    "nicho": "Y: madre + fuera del nicho",
    "tert": "Y: telomerasa + aneuploidía",
    "tert+nicho": "las dos Y en paralelo",
}


def evaluar(sim: Simulacion) -> dict[str, tuple[int, int, int, int]]:
    """Por puerta: (tumorales vistas, tumorales totales, sanas tocadas, sanas totales)."""
    r = {p: [0, 0, 0, 0] for p in PUERTAS}
    sim.p = sim.p.__class__(**{**sim.p.__dict__, "error_lectura_interna": 0.0})
    for pos, c in sim.celulas.items():
        if c.tipo != LUMINAL or c.necrosis or c.guardiana:
            continue
        o = sim._observar(c, pos, 0, 0)
        tumoral = c.n_drivers > 0
        for p in PUERTAS:
            d = r[p]
            d[1 if tumoral else 3] += 1
            if sim._puerta(p, c, o):
                d[0 if tumoral else 2] += 1
    return {p: tuple(v) for p, v in r.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", nargs="*", type=float, default=[120, 200, 300],
                    help="momentos en los que se mide (el tejido sigue creciendo entre ellos)")
    ap.add_argument("--semillas", type=int, default=3)
    ap.add_argument("--neoantigenos", type=float, default=0.1)
    a = ap.parse_args()

    dias = sorted(a.dias)
    total = {d: {p: [0, 0, 0, 0] for p in PUERTAS} for d in dias}
    for semilla in range(a.semillas):
        sim = Simulacion(a.escenario, semilla, Parametros(prob_neoantigeno=a.neoantigenos))
        paso = 0
        for d in dias:
            while paso < int(d * PASOS_POR_DIA):
                sim.paso()
                paso += 1
            for p, v in evaluar(sim).items():
                for i in range(4):
                    total[d][p][i] += v[i]

    print(f"Escenario {a.escenario} · sin tratar · {a.semillas} semillas · "
          f"neoantígenos {a.neoantigenos:g}")
    print("sensibilidad = tumorales que ve · especificidad = sanas que respeta")
    for d in dias:
        t = total[d][PUERTAS[0]]
        print(f"\n— día {d:g}: {t[1]} tumorales · {t[3]} luminales sanas")
        print(f"  {'puerta lógica':<30}{'sensib.':>10}{'especif.':>11}{'sanas tocadas':>15}")
        for p in PUERTAS:
            vistas, tum, tocadas, sanas = total[d][p]
            sens = vistas / tum if tum else float("nan")
            espe = 1 - tocadas / sanas if sanas else float("nan")
            print(f"  {ETIQUETA[p]:<30}{sens:>9.1%}{espe:>11.2%}{tocadas:>15}")


if __name__ == "__main__":
    main()
