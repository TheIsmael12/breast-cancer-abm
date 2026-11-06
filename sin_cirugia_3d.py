"""¿Aguanta en 3D el protocolo sin bisturí?

En 2D (sin_cirugia.py, 2026-11-05) ctDNA mensual + ablación guiada controló 24/24
pacientes destruyendo 8 células sanas. Pero en 2D todo el tumor está cerca de la
superficie; en 3D las terapias no llegan al núcleo (capturas 05-06: la fracción
accesible cae al 33 %). Este banco repite los brazos clave en el cilindro 3D.

Diseño (CRITERIO escrito antes de ejecutar, 2026-11-05):
  · tratamiento el día 60 (en 3D el tumor crece mucho más deprisa que en 2D),
    seguimiento hasta el día 240; tasa de mutación realista (2e-4) desde el día 60.
  · ctDNA calibrado con estas mismas pacientes 3D: sensibilidad del 75 % el día 60.
  · El protocolo AGUANTA si ctDNA + guiada controla (<5 % de las peligrosas del brazo
    "nada" en la mediana) al menos tantas pacientes como el disco + arsenal (el
    equivalente a la cirugía) destruyendo menos de la mitad de tejido sano.
    FALLA si controla menos de la mitad que el disco.

    python sin_cirugia_3d.py --semillas 8 --procesos 6

Umbral de imagen en 3D (--umbrales; CRITERIO escrito antes de ejecutar, 2026-11-06):
  · la imagen solo ve focos de ≥N células (N = 1, 5, 25, 100 o ciega). Si el ctDNA da
    positivo y no se ve nada, no se ablaciona: solo se repite el arsenal (como en 2D).
  · Predicción: en 3D los focos ya miden cientos o miles de células el día 90, así que
    ≥25 y ≥100 deberían controlar casi igual que ≥1; la ciega, como el arsenal solo
    (0/8), porque en 3D los fármacos no llegan al núcleo.
  · El control NO depende de ver focos diminutos si ≥100 controla al menos 7/8.
    SÍ depende si ≥25 baja de 6/8. Entre medias, resultado ambiguo.

    python sin_cirugia_3d.py --umbrales --semillas 8 --procesos 3
"""
from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import time
from multiprocessing import Pool

from biologia import CALIBRACIONES, PASOS_POR_DIA, Parametros
from tejido3d import Simulacion3D

ARSENAL = "troyano+guardiana_nicho+memoria"
CIEGA = 10 ** 9
BRAZOS_UMBRAL = (
    ("nada", None, False),
    ("ctDNA → guiada, imagen ve todo (≥1)", "ctdna_guiada", True),
    ("ctDNA → guiada, imagen ve ≥5", "ctdna_guiada", True, 5),
    ("ctDNA → guiada, imagen ve ≥25", "ctdna_guiada", True, 25),
    ("ctDNA → guiada, imagen ve ≥100", "ctdna_guiada", True, 100),
    ("ctDNA → imagen ciega (solo arsenal)", "ctdna_guiada", True, CIEGA),
)

BRAZOS = (
    ("nada", None, False),
    ("solo arsenal", None, True),
    ("ablación en disco + arsenal (≈ cirugía)", "disco", True),
    ("ctDNA mensual → ablación guiada", "ctdna_guiada", True),
    ("ctDNA mensual → ablación en disco", "ctdna_disco", True),
)


def _nueva(semilla: int) -> Simulacion3D:
    return Simulacion3D("luminal", semilla, Parametros(prob_neoantigeno=0.1, error_edicion=0.01))


CACHE = os.path.join("resultados", "cache_3d")


def _guardado(nombre: str, calcular):
    """Lee el resultado de disco si ya existe; si no, lo calcula y lo guarda.
    Así un banco cortado (falta de memoria, apagón) se retoma donde se quedó."""
    ruta = os.path.join(CACHE, nombre + ".json")
    if os.path.exists(ruta):
        with open(ruta, encoding="utf-8") as fh:
            return json.load(fh)
    valor = calcular()
    os.makedirs(CACHE, exist_ok=True)
    with open(ruta + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(valor, fh)
    os.replace(ruta + ".tmp", ruta)
    return valor


def base_pre(args) -> float:
    semilla, dia = args
    return _guardado(f"base_s{semilla}_d{dia:g}", lambda: _base_pre(semilla, dia))


def _base_pre(semilla: int, dia: float) -> float:
    sim = _nueva(semilla)
    for _ in range(int(dia * PASOS_POR_DIA)):
        sim.paso()
    return sim.ctdna(1.0)["esperado"]


def correr(args) -> dict:
    brazo, semilla, dia, dias, dosis, fragmentos = args
    umbral = brazo[3] if len(brazo) > 3 else 1
    clave = (f"{brazo[1] or ('arsenal' if brazo[2] else 'nada')}"
             f"{'' if umbral == 1 else f'_u{umbral}'}"
             f"_s{semilla}_d{dia:g}-{dias:g}_x{dosis}_f{fragmentos:.4e}")
    return _guardado(clave, lambda: _correr(brazo, semilla, dia, dias, dosis, fragmentos))


def _correr(brazo, semilla, dia, dias, dosis, fragmentos) -> dict:
    _, local, arsenal = brazo[:3]
    umbral = brazo[3] if len(brazo) > 3 else 1
    sim = _nueva(semilla)
    p0 = int(dia * PASOS_POR_DIA)
    controles = {int(d * PASOS_POR_DIA) for d in range(int(dia) + 30, int(dias), 30)}
    out = dict(sanas=0, rescate=None, ablacion=0)
    t0 = time.time()
    for paso in range(int(dias * PASOS_POR_DIA)):
        if paso == p0:
            sim.p.tasa_mutacion = CALIBRACIONES["realista"]
            if local == "disco":
                r = sim.ablacionar(guiada=False, borde=2)
                out["sanas"], out["ablacion"] = r["sanas"], 1
            if arsenal:
                sim.liberar_terapia(ARSENAL, dosis)
        if local and local.startswith("ctdna") and paso in controles and out["rescate"] is None:
            if sim.ctdna(fragmentos)["positiva"]:
                r = sim.ablacionar(guiada=local == "ctdna_guiada",
                                   borde=1 if local == "ctdna_guiada" else 2,
                                   umbral_visible=umbral)
                out["sanas"], out["ablacion"] = r["sanas"], int(r["focos_vistos"] > 0)
                out["rescate"] = round(sim.dia)
                sim.liberar_terapia(ARSENAL, dosis)
        sim.paso()
    r = sim.resumen()
    out.update(peligrosas=r["peligrosas"], tumorales=r["tumorales"],
               accesible=round(sim.fraccion_accesible(), 2), segundos=round(time.time() - t0))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--semillas", type=int, default=8)
    ap.add_argument("--dia", type=float, default=60)
    ap.add_argument("--dias", type=float, default=240)
    ap.add_argument("--dosis", type=int, default=300)
    ap.add_argument("--procesos", type=int, default=6)
    ap.add_argument("--umbrales", action="store_true",
                    help="umbral de imagen: la guiada solo ve focos de ≥N células")
    a = ap.parse_args()
    brazos = BRAZOS_UMBRAL if a.umbrales else BRAZOS
    semillas = list(range(1, a.semillas + 1))
    t0 = time.time()
    with Pool(a.procesos, maxtasksperchild=1) as pool:
        bases = pool.map(base_pre, [(s, a.dia) for s in semillas], chunksize=1)
        lo, hi = 1e-8, 10.0
        for _ in range(80):
            f = math.sqrt(lo * hi)
            sens = statistics.mean(1 - math.exp(-f * b) for b in bases)
            lo, hi = (f, hi) if sens < 0.75 else (lo, f)
        print(f"SIN CIRUGÍA EN 3D · luminal · tratamiento día {a.dia:g} · seguimiento a "
              f"{a.dias:g} d · {len(semillas)} semillas · dosis {a.dosis}")
        print(f"ctDNA calibrado en 3D: {f:.2e} fragmentos/célula (sensibilidad 75 % el día "
              f"{a.dia:g}; carcinoma-equivalente antes de tratar: "
              f"{sorted(round(b) for b in bases)})\n", flush=True)
        tareas = [(b, s, a.dia, a.dias, a.dosis, f) for b in brazos for s in semillas]
        res = pool.map(correr, tareas, chunksize=1)
    n = len(semillas)
    nada = statistics.median(r["peligrosas"] for r in res[:n])
    umbral = 0.05 * nada
    print(f"control = <5 % de la mediana sin tratar ({umbral:.0f} peligrosas)\n")
    print(f"{'brazo':<42}{'mediana':>9}{'control':>9}{'a 0':>5}{'sanas':>8}   por semilla")
    for i, b in enumerate(brazos):
        rs = res[i * n:(i + 1) * n]
        pel = [r["peligrosas"] for r in rs]
        print(f"{b[0]:<42}{statistics.median(pel):>9.0f}"
              f"{sum(p < umbral for p in pel):>6}/{n:<2}{sum(p == 0 for p in pel):>5}"
              f"{statistics.median(r['sanas'] for r in rs):>8.0f}   {pel}")
        if b[1] and b[1].startswith("ctdna"):
            print(f"{'':<42}rescates (día): {[r['rescate'] for r in rs]} · "
                  f"sin tocar: {sum(not r['ablacion'] for r in rs)}/{n} · "
                  f"sanas: {[r['sanas'] for r in rs]}")
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
