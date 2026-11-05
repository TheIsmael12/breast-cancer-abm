"""¿Se puede llegar al interior de un tumor sólido? Se mide ANTES de probar ninguna terapia.

Un tumor no está envuelto en una membrana. Lo que impide entrar son tres cosas distintas,
y conviene no confundirlas porque cada una se levanta de una manera:

  1. GEOMETRÍA      una esfera esconde su interior: la fracción accesible cae como 1/r.
                    No se arregla con fármacos; se arregla entrando por dentro (aguja).
  2. PRESIÓN        vasos porosos + linfáticos colapsados = el líquido empuja hacia fuera.
                    Es corriente en contra, no una pared. Baja al degradar la matriz o al
                    normalizar los vasos (Jain 1987-2014).
  3. OXÍGENO        el O2 difunde ~100-150 µm desde un capilar: el centro está lejos, no
                    tapiado. Y la hipoxia no solo apaga a la terapia: activa la
                    desdiferenciación del tumor (HIF-1α), que es lo que fabrica reinas.

Este banco mide las tres por separado, con una terapia inerte que solo camina, para que
lo que se vea sea el ACCESO y no la eficacia.

    python barrera.py --dias 300 --semillas 2
"""
from __future__ import annotations

import argparse
import time

import numpy as np

from biologia import ESCENARIOS, PASOS_POR_DIA, Parametros
from tejido import LUMINAL
from tejido3d import Simulacion3D

# (etiqueta, vía de entrega, fracción de matriz degradada, modo de vasos, sin presión)
CONDICIONES = (
    ("nada (vena, tal cual)",          "vena",         0.0,  "normal",     False),
    ("degradar matriz 50 %",           "vena",         0.5,  "normal",     False),
    ("normalizar vasos",               "vena",         0.0,  "normalizar", False),
    ("matriz + normalizar",            "vena",         0.5,  "normalizar", False),
    ("reconducir O2 al tumor",         "vena",         0.0,  "reconducir", False),
    ("aguja (una punción)",            "aguja",        0.0,  "normal",     False),
    ("intratumoral repartido",         "intratumoral", 0.0,  "normal",     False),
    ("intratumoral + matriz + vasos",  "intratumoral", 0.5,  "normalizar", False),
    ("SIN presión (control irreal)",   "vena",         0.0,  "normal",     True),
)


def _tumor(sim: Simulacion3D) -> list:
    return [q for q, c in sim.celulas.items()
            if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]


def _radios(sim: Simulacion3D, tumor: list):
    """Distancia de cada célula tumoral al centroide de la masa."""
    a = np.array(tumor, float)
    return np.linalg.norm(a - a.mean(0), axis=1)


def medir(sim: Simulacion3D, alcanzadas: set) -> dict:
    tumor = _tumor(sim)
    if not tumor:
        return dict(tumorales=0, o2_borde=0.0, o2_centro=0.0, accesible=0.0,
                    alcanzado=0.0, hipoxico=0.0)
    r = _radios(sim, tumor)
    orden = np.argsort(r)
    n = max(1, len(tumor) // 5)
    o2 = np.array([sim.O2[q] for q in tumor])
    nucleo = {tumor[i] for i in orden[:n]}
    return dict(
        tumorales=len(tumor),
        o2_borde=float(o2[orden[-n:]].mean()),
        o2_centro=float(o2[orden[:n]].mean()),
        accesible=sim.fraccion_accesible(),
        # de las células del NÚCLEO (el quinto más interior), ¿a cuántas llegó la terapia?
        alcanzado=len(nucleo & alcanzadas) / len(nucleo),
        hipoxico=float((o2 < sim.p.o2_hipoxia).mean()),
    )


def correr(escenario: str, semilla: int, dias: float, dia_terapia: float, via: str,
           matriz: float, vasos: str, sin_presion: bool, neoantigenos: float) -> dict:
    p = Parametros(prob_neoantigeno=neoantigenos,
                   presion_intersticial=0.0 if sin_presion else 1.0)
    sim = Simulacion3D(escenario, semilla, p)
    paso_terapia = int(dia_terapia * PASOS_POR_DIA)
    alcanzadas: set = set()
    for paso in range(int(dias * PASOS_POR_DIA)):
        if paso == paso_terapia:
            if matriz:
                sim.degradar_matriz(matriz)
            if vasos != "normal":
                sim.tratar_vasos(vasos)
            # Sonda inerte: camina y nada más. Mide acceso, no eficacia.
            sim.liberar_terapia("reparador", 40, via)
        sim.paso()
        if paso > paso_terapia and paso % 4 == 0:
            for ag in sim.agentes:
                if ag.tipo != "IA":
                    continue
                z, y, x = ag.pos
                for dz in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dx in (-1, 0, 1):
                            q = (z + dz, y + dy, x + dx)
                            c = sim.celulas.get(q)
                            if c is not None and c.tipo == LUMINAL and c.n_drivers:
                                alcanzadas.add(q)
    r = medir(sim, alcanzadas)
    r["bloqueados"] = sim.stats["pasos_bloqueados"]
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--escenario", choices=ESCENARIOS, default="luminal")
    ap.add_argument("--dias", type=float, default=300)
    ap.add_argument("--dia-terapia", type=float, default=200)
    ap.add_argument("--semillas", type=int, default=2)
    ap.add_argument("--neoantigenos", type=float, default=0.1)
    a = ap.parse_args()

    print(f"3D · {a.escenario} · {a.dias:g} días · acceso medido desde el día {a.dia_terapia:g}"
          f" · {a.semillas} semillas")
    print("núcleo alcanzado = del quinto MÁS INTERIOR del tumor, qué fracción llegó a tener "
          "una célula terapéutica al lado\n")
    cab = f"{'condición':<32}{'tumorales':>10}{'O2 borde':>10}{'O2 centro':>11}" \
          f"{'accesible':>11}{'núcleo alcanz.':>16}{'bloqueados':>12}"
    print(cab)
    for etiqueta, via, matriz, vasos, sin_presion in CONDICIONES:
        t0 = time.time()
        res = [correr(a.escenario, s, a.dias, a.dia_terapia, via, matriz, vasos,
                      sin_presion, a.neoantigenos) for s in range(a.semillas)]
        m = {k: sum(r[k] for r in res) / len(res) for k in res[0]}
        print(f"{etiqueta:<32}{m['tumorales']:>10.0f}{m['o2_borde']:>10.2f}"
              f"{m['o2_centro']:>11.2f}{m['accesible']:>10.0%}{m['alcanzado']:>16.0%}"
              f"{m['bloqueados']:>12.0f}   ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
