"""¿Se puede llegar al mismo sitio que la cirugía SIN operar?

Lo que ha aprendido el mundo entre 2022 y 2026, y que este banco traslada al modelo:

  · Kuerer (MD Anderson; Lancet Oncol 2022, JAMA Oncol 2025): tratamiento sistémico
    PRIMERO, biopsia por vacío del lecho, y si sale limpia NO se opera. 31 pacientes
    HER2+ o triple negativo sin cirugía: 0 % de recaída en la mama a ~5 años.
  · FIRST (crioablación, 2026): destruir el tumor in situ con frío. Ablación completa
    89,6 % (96,9 % si ≤2 cm), pero CDIS residual en el BORDE en el 12,5 %.
  · MD Anderson, ESTRO 2026: hormonoterapia + radioterapia ablativa en luminales
    pequeños → 53 % de respuestas completas en biopsia, sin operar.
  · El punto débil común: la BIOPSIA de comprobación (RESPONDER: falsos negativos 17,8 %).

Brazos (todo el día 120, seguimiento hasta el día 500, misma semilla = misma paciente):
  nada · solo arsenal · cirugía + arsenal (la referencia del 22-10) · ablación + arsenal
  (borde letal y borde con supervivientes) · vía Kuerer (arsenal → biopsia → ablación
  solo si la biopsia sale positiva).

    bcabm study no-surgery --seeds 12 --processes 6
"""

from __future__ import annotations

import argparse
import statistics
import time
from multiprocessing import Pool

from breast_cancer_abm.model.biology import CALIBRATIONS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation

ARSENAL = "trojan+niche_guardian+memory"
CONTROL_THRESHOLD = 100  # suelo absoluto de la métrica (corrección del 2026-10-22 (7))

# (etiqueta, local: None | "cirugia" | "ablacion" | "kuerer", rim_survival, arsenal)
ARMS = (
    ("nada", None, 0.0, False),
    ("solo arsenal (sin tocar)", None, 0.0, True),
    ("cirugía + arsenal", "surgery", 0.0, True),
    ("ablación + arsenal (borde letal)", "ablation", 0.0, True),
    ("ablación + arsenal (borde 30 % vivo)", "ablation", 0.3, True),
    ("Kuerer: arsenal → biopsia → ablación", "kuerer", 0.3, True),
    # Puertas de residual_test.py (2026-11-05): el ctDNA seriado, no de una vez.
    ("Kuerer + ctDNA mensual", "kuerer+ctdna", 0.3, True),
    ("solo ctDNA mensual (sin pinchar nunca)", "ctdna", 0.3, True),
    # 2026-11-05 (3): ablación guiada por imagen y tasa de mutación REALISTA tras tratar
    # (acelerada hasta el día del tratamiento para que haya tumor; 2e-4 después, para que
    # el campo no fabrique lesiones nuevas a ritmo de laboratorio).
    ("ctDNA + ablación GUIADA", "ctdna", 0.3, True, dict(guided=True)),
    ("REAL: solo arsenal", None, 0.0, True, dict(rate_after="realistic")),
    ("REAL: ctDNA + ablación disco", "ctdna", 0.3, True, dict(rate_after="realistic")),
    ("REAL: ctDNA + ablación GUIADA", "ctdna", 0.3, True, dict(guided=True, rate_after="realistic")),
    # 2026-11-05 (4): la imagen solo ve focos de >= N células. Si el ctDNA da positivo y
    # no se ve nada, solo se repite el arsenal (tratamiento sistémico, sin tocar la mama).
    (
        "IMAGEN≥5: ctDNA + guiada",
        "ctdna",
        0.3,
        True,
        dict(guided=True, rate_after="realistic", visible_threshold=5),
    ),
    (
        "IMAGEN≥25: ctDNA + guiada",
        "ctdna",
        0.3,
        True,
        dict(guided=True, rate_after="realistic", visible_threshold=25),
    ),
    (
        "IMAGEN≥100: ctDNA + guiada",
        "ctdna",
        0.3,
        True,
        dict(guided=True, rate_after="realistic", visible_threshold=100),
    ),
    (
        "IMAGEN ciega: ctDNA → solo arsenal",
        "ctdna",
        0.3,
        True,
        dict(guided=True, rate_after="realistic", visible_threshold=10**9),
    ),
)
# Fragmentos/célula calibrados en residual_test.py para sensibilidad pre-tratamiento 75 %.
CTDNA_FRAGMENTS = 9.91e-2
CTDNA_EVERY = 30  # días entre extracciones
CTDNA_UNTIL = 480


def run(args) -> dict:
    arm, seed, days, day, dose, biopsy_wait = args
    _, local, rim_surv, arsenal = arm[:4]
    rate_after = (arm[4] if len(arm) > 4 else {}).get("rate_after")
    if local in ("kuerer+ctdna", "ctdna"):
        return run_ctdna(arm, seed, days, day, dose, biopsy_wait)
    sim = Simulation("luminal", seed, Parameters(prob_neoantigen=0.1, editing_error=0.01))
    p0 = int(day * STEPS_PER_DAY)
    p_therapy = p0 + (4 * STEPS_PER_DAY if local in ("surgery", "ablation") else 0)
    p_biopsy = p0 + int(biopsy_wait * STEPS_PER_DAY)
    out = dict(healthy_destroyed=0, biopsy=None, real_in_biopsy=None, operated=0)
    for step in range(int(days * STEPS_PER_DAY)):
        if step == p0 and rate_after:
            sim.p.mutation_rate = CALIBRATIONS[rate_after]
        if step == p0 and local == "surgery":
            out["healthy_destroyed"] = sim.resect()["healthy"]
        elif step == p0 and local == "ablation":
            out["healthy_destroyed"] = sim.ablate(rim_survival=rim_surv)["healthy"]
        if step == p_therapy and arsenal:
            sim.release_therapy(ARSENAL, dose)
        if local == "kuerer" and step == p_biopsy:
            b = sim.biopsy()
            out["biopsy"], out["real_in_biopsy"] = b["positive"], b["real"]
            if b["positive"]:  # rescate local SIN bisturí, y segunda dosis del arsenal
                out["healthy_destroyed"] = sim.ablate(rim_survival=rim_surv)["healthy"]
                out["operated"] = 1
        if local == "kuerer" and out["operated"] and step == p_biopsy + 4 * STEPS_PER_DAY:
            sim.release_therapy(ARSENAL, dose)
        sim.step()
    r = sim.summary()
    out.update(dangerous=r["dangerous"], queens=r["queens"], invasive_cells=r["invasive_cells"])
    return out


def run_ctdna(arm, seed, days, day, dose, biopsy_wait) -> dict:
    """Arsenal el día `dia`; a los `biopsy_wait` días primera comprobación (biopsia o
    ctDNA según el brazo) y después ctDNA cada CTDNA_EVERY días. El primer positivo
    dispara ablación + arsenal UNA vez; después se sigue vigilando sin más rescates."""
    _, local, rim_surv, _ = arm[:4]
    op = arm[4] if len(arm) > 4 else {}
    guided, rate_after = op.get("guided", False), op.get("rate_after")
    threshold = op.get("visible_threshold", 1)
    sim = Simulation("luminal", seed, Parameters(prob_neoantigen=0.1, editing_error=0.01))
    p0 = int(day * STEPS_PER_DAY)
    first = day + biopsy_wait
    checkpoints = {int(d * STEPS_PER_DAY) for d in range(int(first), CTDNA_UNTIL + 1, CTDNA_EVERY)}
    out = dict(
        healthy_destroyed=0,
        biopsy=None,
        real_in_biopsy=None,
        operated=0,
        rescue_day=None,
        draws=0,
        ablation=0,
    )
    p_rearsenal = None
    for step in range(int(days * STEPS_PER_DAY)):
        if step == p0:
            sim.release_therapy(ARSENAL, dose)
            if rate_after:
                sim.p.mutation_rate = CALIBRATIONS[rate_after]
        if step in checkpoints and not out["operated"]:
            if local == "kuerer+ctdna" and step == int(first * STEPS_PER_DAY):
                b = sim.biopsy()
                out["biopsy"], out["real_in_biopsy"] = b["positive"], b["real"]
                is_positive = b["positive"] or sim.ctdna(CTDNA_FRAGMENTS)["positive"]
            else:
                is_positive = sim.ctdna(CTDNA_FRAGMENTS)["positive"]
            out["draws"] += 1
            if is_positive:
                ab = sim.ablate(
                    rim_survival=rim_surv, guided=guided, rim=1 if guided else 2, visible_threshold=threshold
                )
                out["healthy_destroyed"] = ab["healthy"]
                out["ablation"] = int(ab.get("foci_seen", 1) > 0)
                out["operated"], out["rescue_day"] = 1, round(sim.day)
                p_rearsenal = step + 4 * STEPS_PER_DAY
        if step == p_rearsenal:
            sim.release_therapy(ARSENAL, dose)
        sim.step()
    r = sim.summary()
    out.update(dangerous=r["dangerous"], queens=r["queens"], invasive_cells=r["invasive_cells"])
    return out


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=12)
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--days", type=float, default=500)
    ap.add_argument("--day", type=float, default=120)
    ap.add_argument("--dose", type=int, default=300)
    ap.add_argument(
        "--biopsy-wait", type=float, default=60, help="días entre el arsenal y la biopsia en la vía Kuerer"
    )
    ap.add_argument("--processes", type=int, default=6, help="procesos en paralelo (deja núcleos libres)")
    ap.add_argument("--only", default="", help="ejecutar solo los brazos cuyo nombre contenga este texto")
    a = ap.parse_args(argv)
    seeds = range(a.first, a.first + a.seeds)
    global ARMS
    ARMS = tuple(b for b in ARMS if a.only.lower() in b[0].lower())

    print(
        f"SIN CIRUGÍA · luminal · tratamiento día {a.day:g} · seguimiento a {a.days:g} d · "
        f"{a.seeds} semillas · control = <{CONTROL_THRESHOLD} peligrosas\n"
    )
    print(f"{'brazo':<40}{'mediana':>9}{'control':>9}{'sanas dest.':>13}   por semilla")
    tasks = [(b, s, a.days, a.day, a.dose, a.biopsy_wait) for b in ARMS for s in seeds]
    t0 = time.time()
    with Pool(a.processes, maxtasksperchild=1) as pool:
        res = pool.map(run, tasks, chunksize=1)
    n = len(seeds)
    for i, b in enumerate(ARMS):
        rs = res[i * n : (i + 1) * n]
        dang = [r["dangerous"] for r in rs]
        ctrl = sum(p < CONTROL_THRESHOLD for p in dang)
        healthy = statistics.median(r["healthy_destroyed"] for r in rs)
        print(f"{b[0]:<40}{statistics.median(dang):>9.0f}{ctrl:>6}/{n:<2}{healthy:>13.0f}   {dang}")
        if b[1] == "kuerer":
            neg = [r for r in rs if r["biopsy"] is False]
            fn = sum(r["real_in_biopsy"] > 0 for r in neg)
            living = sum(r["real_in_biopsy"] > 0 for r in rs)
            print(
                f"{'':<40}biopsia negativa (no se toca) en {len(neg)}/{n}; "
                f"de ellas con tumor vivo en el tejido: {fn} "
                f"(falsos negativos {fn}/{max(living, 1)} de las que tenían tumor)"
            )
            print(f"{'':<40}desenlace de las NO tocadas: {[r['dangerous'] for r in neg]}")
        if b[1] in ("kuerer+ctdna", "ctdna"):
            touched = [r for r in rs if r["operated"]]
            free = [r for r in rs if not r["ablation"]]
            print(
                f"{'':<40}sin tocar nunca la mama: {len(free)}/{n} → desenlace "
                f"{[r['dangerous'] for r in free]}"
            )
            print(
                f"{'':<40}ctDNA positivo pero nada visible (solo arsenal): "
                f"{sum(r['operated'] and not r['ablation'] for r in rs)}/{n}"
            )
            print(f"{'':<40}rescatadas: días {sorted(r['rescue_day'] for r in touched)}")
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
