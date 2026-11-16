"""¿Aguanta en 3D el protocolo sin bisturí?

En 2D (no_surgery.py, 2026-11-05) ctDNA mensual + ablación guiada controló 24/24
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

    bcabm study no-surgery-3d --seeds 8 --processes 6

Umbral de imagen en 3D (--umbrales; CRITERIO escrito antes de ejecutar, 2026-11-06):
  · la imagen solo ve focos de ≥N células (N = 1, 5, 25, 100 o ciega). Si el ctDNA da
    positivo y no se ve nada, no se ablaciona: solo se repite el arsenal (como en 2D).
  · Predicción: en 3D los focos ya miden cientos o miles de células el día 90, así que
    ≥25 y ≥100 deberían controlar casi igual que ≥1; la ciega, como el arsenal solo
    (0/8), porque en 3D los fármacos no llegan al núcleo.
  · El control NO depende de ver focos diminutos si ≥100 controla al menos 7/8.
    SÍ depende si ≥25 baja de 6/8. Entre medias, resultado ambiguo.

    bcabm study no-surgery-3d --thresholds --seeds 8 --processes 3
"""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import time
from multiprocessing import Pool

from breast_cancer_abm.model.biology import CALIBRATIONS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue3d import Simulation3D
from breast_cancer_abm.paths import CACHE_3D_DIR

ARSENAL = "trojan+niche_guardian+memory"
BLIND = 10**9
THRESHOLD_ARMS = (
    ("nada", None, False),
    ("ctDNA → guiada, imagen ve todo (≥1)", "ctdna_guided", True),
    ("ctDNA → guiada, imagen ve ≥5", "ctdna_guided", True, 5),
    ("ctDNA → guiada, imagen ve ≥25", "ctdna_guided", True, 25),
    ("ctDNA → guiada, imagen ve ≥100", "ctdna_guided", True, 100),
    ("ctDNA → imagen ciega (solo arsenal)", "ctdna_guided", True, BLIND),
)

ARMS = (
    ("nada", None, False),
    ("solo arsenal", None, True),
    ("ablación en disco + arsenal (≈ cirugía)", "disc", True),
    ("ctDNA mensual → ablación guiada", "ctdna_guided", True),
    ("ctDNA mensual → ablación en disco", "ctdna_disc", True),
)


def _new_sim(seed: int) -> Simulation3D:
    return Simulation3D("luminal", seed, Parameters(prob_neoantigen=0.1, editing_error=0.01))


CACHE = str(CACHE_3D_DIR)


def _cached(name: str, compute):
    """Lee el resultado de disco si ya existe; si no, lo calcula y lo guarda.
    Así un banco cortado (falta de memoria, apagón) se retoma donde se quedó."""
    file_path = os.path.join(CACHE, name + ".json")
    if os.path.exists(file_path):
        with open(file_path, encoding="utf-8") as fh:
            return json.load(fh)
    value = compute()
    os.makedirs(CACHE, exist_ok=True)
    with open(file_path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(value, fh)
    os.replace(file_path + ".tmp", file_path)
    return value


def base_pre(args) -> float:
    seed, day = args
    return _cached(f"base_s{seed}_d{day:g}", lambda: _baseline_pre(seed, day))


def _baseline_pre(seed: int, day: float) -> float:
    sim = _new_sim(seed)
    for _ in range(int(day * STEPS_PER_DAY)):
        sim.step()
    return sim.ctdna(1.0)["expected"]


def run(args) -> dict:
    arm, seed, day, days, dose, fragments = args
    threshold = arm[3] if len(arm) > 3 else 1
    key = (
        f"{arm[1] or ('arsenal' if arm[2] else 'nothing')}"
        f"{'' if threshold == 1 else f'_u{threshold}'}"
        f"_s{seed}_d{day:g}-{days:g}_x{dose}_f{fragments:.4e}"
    )
    return _cached(key, lambda: _run(arm, seed, day, days, dose, fragments))


def _run(arm, seed, day, days, dose, fragments) -> dict:
    _, local, arsenal = arm[:3]
    threshold = arm[3] if len(arm) > 3 else 1
    sim = _new_sim(seed)
    p0 = int(day * STEPS_PER_DAY)
    checkpoints = {int(d * STEPS_PER_DAY) for d in range(int(day) + 30, int(days), 30)}
    out = dict(healthy=0, rescue=None, ablation=0)
    t0 = time.time()
    for step in range(int(days * STEPS_PER_DAY)):
        if step == p0:
            sim.p.mutation_rate = CALIBRATIONS["realistic"]
            if local == "disc":
                r = sim.ablate(guided=False, rim=2)
                out["healthy"], out["ablation"] = r["healthy"], 1
            if arsenal:
                sim.release_therapy(ARSENAL, dose)
        if local and local.startswith("ctdna") and step in checkpoints and out["rescue"] is None:
            if sim.ctdna(fragments)["positive"]:
                r = sim.ablate(
                    guided=local == "ctdna_guided",
                    rim=1 if local == "ctdna_guided" else 2,
                    visible_threshold=threshold,
                )
                out["healthy"], out["ablation"] = r["healthy"], int(r["foci_seen"] > 0)
                out["rescue"] = round(sim.day)
                sim.release_therapy(ARSENAL, dose)
        sim.step()
    r = sim.summary()
    out.update(
        dangerous=r["dangerous"],
        tumor_cells=r["tumor_cells"],
        accessible=round(sim.accessible_fraction(), 2),
        seconds=round(time.time() - t0),
    )
    return out


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--day", type=float, default=60)
    ap.add_argument("--days", type=float, default=240)
    ap.add_argument("--dose", type=int, default=300)
    ap.add_argument("--processes", type=int, default=6)
    ap.add_argument(
        "--thresholds", action="store_true", help="umbral de imagen: la guiada solo ve focos de ≥N células"
    )
    a = ap.parse_args(argv)
    arms = THRESHOLD_ARMS if a.thresholds else ARMS
    seeds = list(range(1, a.seeds + 1))
    t0 = time.time()
    with Pool(a.processes, maxtasksperchild=1) as pool:
        bases = pool.map(base_pre, [(s, a.day) for s in seeds], chunksize=1)
        lo, hi = 1e-8, 10.0
        for _ in range(80):
            f = math.sqrt(lo * hi)
            sens = statistics.mean(1 - math.exp(-f * b) for b in bases)
            lo, hi = (f, hi) if sens < 0.75 else (lo, f)
        print(
            f"SIN CIRUGÍA EN 3D · luminal · tratamiento día {a.day:g} · seguimiento a "
            f"{a.days:g} d · {len(seeds)} semillas · dosis {a.dose}"
        )
        print(
            f"ctDNA calibrado en 3D: {f:.2e} fragmentos/célula (sensibilidad 75 % el día "
            f"{a.day:g}; carcinoma-equivalente antes de tratar: "
            f"{sorted(round(b) for b in bases)})\n",
            flush=True,
        )
        tasks = [(b, s, a.day, a.days, a.dose, f) for b in arms for s in seeds]
        res = pool.map(run, tasks, chunksize=1)
    n = len(seeds)
    nothing = statistics.median(r["dangerous"] for r in res[:n])
    threshold = 0.05 * nothing
    print(f"control = <5 % de la mediana sin tratar ({threshold:.0f} peligrosas)\n")
    print(f"{'brazo':<42}{'mediana':>9}{'control':>9}{'a 0':>5}{'sanas':>8}   por semilla")
    for i, b in enumerate(arms):
        rs = res[i * n : (i + 1) * n]
        dang = [r["dangerous"] for r in rs]
        print(
            f"{b[0]:<42}{statistics.median(dang):>9.0f}"
            f"{sum(p < threshold for p in dang):>6}/{n:<2}{sum(p == 0 for p in dang):>5}"
            f"{statistics.median(r['healthy'] for r in rs):>8.0f}   {dang}"
        )
        if b[1] and b[1].startswith("ctdna"):
            print(
                f"{'':<42}rescates (día): {[r['rescue'] for r in rs]} · "
                f"sin tocar: {sum(not r['ablation'] for r in rs)}/{n} · "
                f"sanas: {[r['healthy'] for r in rs]}"
            )
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
