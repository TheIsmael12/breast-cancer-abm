"""Experimentos diseñados sin programar: pacientes, brazos de tratamiento e intervenciones.

Un experimento se describe con datos (se guarda como JSON) y se ejecuta en paralelo,
una tarea por paciente y brazo. Las mismas semillas se usan en todos los brazos
(diseño emparejado): la única diferencia entre brazos es el tratamiento.
"""

from __future__ import annotations

import json
import os
import statistics
import time
from collections.abc import Callable
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass, field, replace

from breast_cancer_abm.catalog import VESSEL_MODES, therapy_label
from breast_cancer_abm.model.biology import CALIBRATIONS, STEPS_PER_DAY, Parameters

# tipo -> (nombre, disponible en 2D, disponible en 3D)
INTERVENTIONS = {
    "therapy": ("Liberar terapia celular", True, True),
    "ctdna_monitor": ("Vigilar con ctDNA y rescatar si da positivo", True, True),
    "ablate": ("Ablación (destruir el volumen sin cirugía)", True, True),
    "resect": ("Cirugía (extirpar)", True, False),
    "block_pd1": ("Inmunoterapia anti-PD-1", True, False),
    "block_cd47": ("Anticuerpo anti-CD47", True, False),
    "vessels": ("Tratar los vasos sanguíneos", True, True),
    "inhibit_telomerase": ("Inhibidor de telomerasa", True, False),
    "damage_telomeres": ("Dañar telómeros (6-tio-dG)", True, True),
    "realistic_rate": ("Pasar a la tasa de mutación realista", True, True),
}


@dataclass
class Intervention:
    kind: str
    day: float
    therapy: str = ""  # terapia o combinación "a+b"
    dose: int = 300
    guided: bool = True  # ablación guiada por imagen o en disco
    visible_threshold: int = 1  # la imagen solo ve focos de >= N células
    rim: int = 1
    duration: float = 60.0  # días del anticuerpo anti-PD-1
    mode: str = "reroute"  # modo de los vasos
    every: float = 30.0  # ctDNA: cada cuántos días
    until: float = 0.0  # ctDNA: hasta qué día (0 = hasta el final)
    fragments: float = 0.05  # ctDNA: fragmentos por célula (sensibilidad)
    rescue: str = "guided"  # ctDNA positivo -> "guided" | "disc" | "therapy"

    def describe(self) -> str:
        k = self.kind
        if k == "therapy":
            return f"Día {self.day:g} · {therapy_label(self.therapy)} ({self.dose} células)"
        if k == "ablate":
            how = "guiada por imagen" if self.guided else "en disco"
            seen = (
                f", la imagen ve focos ≥{self.visible_threshold}"
                if self.guided and self.visible_threshold > 1
                else ""
            )
            return f"Día {self.day:g} · ablación {how}{seen}"
        if k == "ctdna_monitor":
            rescue = {"guided": "ablación guiada", "disc": "ablación en disco"}.get(
                self.rescue, therapy_label(self.therapy) if self.therapy else "terapia"
            )
            end = f" hasta el día {self.until:g}" if self.until else ""
            return (
                f"Desde el día {self.day:g}, ctDNA cada {self.every:g} días{end}; si da positivo: {rescue}"
                + (f" + {therapy_label(self.therapy)}" if self.therapy and self.rescue != "therapy" else "")
            )
        if k == "block_pd1":
            return f"Día {self.day:g} · anti-PD-1 durante {self.duration:g} días"
        if k == "vessels":
            return f"Día {self.day:g} · vasos: {VESSEL_MODES.get(self.mode, self.mode).lower()}"
        return f"Día {self.day:g} · {INTERVENTIONS[k][0].lower()}"


@dataclass
class Arm:
    name: str
    interventions: list[Intervention] = field(default_factory=list)


@dataclass
class Experiment:
    name: str = "Nuevo experimento"
    question: str = ""
    model: str = "2d"  # "2d" | "3d"
    scenario: str = "luminal"
    patients: int = 6
    first_seed: int = 1
    days: float = 300
    size3d: str = "fast"  # "fast" | "standard"
    params: dict = field(default_factory=lambda: {"prob_neoantigen": 0.1, "editing_error": 0.01})
    arms: list[Arm] = field(default_factory=list)
    control_cells: int = 100  # controlada = menos de N células peligrosas al final

    # ------------------------------------------------------------------ persistencia
    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> Experiment:
        d = dict(d)
        arms = [
            Arm(a["name"], [Intervention(**i) for i in a.get("interventions", [])]) for a in d.pop("arms", [])
        ]
        known = {f for f in cls.__dataclass_fields__}
        return cls(**{k: v for k, v in d.items() if k in known}, arms=arms)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=1)

    def estimate_seconds(self, processes: int) -> float:
        per_step = 0.06 if self.model == "2d" else (0.15 if self.size3d == "fast" else 0.8)
        tasks = self.patients * max(1, len(self.arms))
        return tasks * self.days * STEPS_PER_DAY * per_step / max(1, min(processes, tasks))


SIZES_3D = {"fast": dict(side=34, z_len=16, radius=9.0), "standard": dict(side=54, z_len=40, radius=13.0)}


def _new_sim(exp: Experiment, seed: int):
    params = replace(Parameters(), **exp.params)
    if exp.model == "3d":
        from breast_cancer_abm.model.tissue3d import Simulation3D

        return Simulation3D(exp.scenario, seed, params, **SIZES_3D[exp.size3d])
    from breast_cancer_abm.model.tissue import Simulation

    return Simulation(exp.scenario, seed, params)


def _apply(sim, iv: Intervention, log: list, healthy: list) -> None:
    day = round(sim.day, 1)
    k = iv.kind
    if k == "therapy" and iv.therapy:
        sim.release_therapy(iv.therapy, iv.dose)
        log.append((day, iv.describe()))
    elif k == "ablate":
        r = sim.ablate(
            guided=iv.guided,
            rim=iv.rim if iv.guided else max(iv.rim, 2),
            visible_threshold=iv.visible_threshold,
        )
        healthy.append(r["healthy"])
        log.append((day, f"Ablación: {r['destroyed']} tumorales y {r['healthy']} sanas destruidas"))
    elif k == "resect":
        r = sim.resect()
        healthy.append(r["healthy"])
        log.append((day, f"Cirugía: {r['excised']} células extirpadas, {r['healthy']} sanas"))
    elif k == "block_pd1":
        sim.block_pd1(iv.duration)
        log.append((day, "Anti-PD-1"))
    elif k == "block_cd47":
        sim.block_cd47()
        log.append((day, "Anti-CD47"))
    elif k == "vessels":
        sim.treat_vessels(iv.mode)
        log.append((day, f"Vasos: {iv.mode}"))
    elif k == "inhibit_telomerase":
        sim.inhibit_telomerase()
        log.append((day, "Inhibidor de telomerasa"))
    elif k == "damage_telomeres":
        sim.damage_telomeres()
        log.append((day, "6-tio-dG"))
    elif k == "realistic_rate":
        sim.p.mutation_rate = CALIBRATIONS["realistic"]
        log.append((day, "Tasa de mutación realista"))


def run_patient(exp_dict: dict, arm_index: int, seed: int) -> dict:
    """Una paciente en un brazo. Función de nivel superior para poder ir a otro proceso."""
    exp = Experiment.from_dict(exp_dict)
    arm = exp.arms[arm_index]
    sim = _new_sim(exp, seed)
    total = int(exp.days * STEPS_PER_DAY)
    pending = sorted((iv for iv in arm.interventions if iv.kind != "ctdna_monitor"), key=lambda i: i.day)
    monitors = [iv for iv in arm.interventions if iv.kind == "ctdna_monitor"]
    rescued = {id(m): False for m in monitors}
    log: list = []
    healthy: list[int] = []
    series = {"day": [], "dangerous": [], "tumor_cells": [], "queens": []}
    rescue_day = None
    t0 = time.time()
    for step in range(total + 1):
        day = step / STEPS_PER_DAY
        while pending and pending[0].day <= day:
            _apply(sim, pending.pop(0), log, healthy)
        for m in monitors:
            until = m.until or exp.days
            if rescued[id(m)] or day < m.day or day > until:
                continue
            if step % max(1, int(m.every * STEPS_PER_DAY)) == int(m.day * STEPS_PER_DAY) % max(
                1, int(m.every * STEPS_PER_DAY)
            ):
                test = sim.ctdna(m.fragments)
                if test["positive"]:
                    rescued[id(m)] = True
                    rescue_day = round(day, 1)
                    log.append((rescue_day, f"ctDNA positivo ({test['fragments']} fragmentos) → rescate"))
                    if m.rescue in ("guided", "disc"):
                        _apply(
                            sim,
                            Intervention(
                                "ablate",
                                day,
                                guided=m.rescue == "guided",
                                visible_threshold=1,
                                rim=1 if m.rescue == "guided" else 2,
                            ),
                            log,
                            healthy,
                        )
                    if m.therapy:
                        _apply(
                            sim, Intervention("therapy", day, therapy=m.therapy, dose=m.dose), log, healthy
                        )
        if step % STEPS_PER_DAY == 0:
            r = sim.summary()
            series["day"].append(round(day, 2))
            for key in ("dangerous", "tumor_cells", "queens"):
                series[key].append(r[key])
        if step < total:
            sim.step()
    final = sim.summary()
    return {
        "arm": arm_index,
        "seed": seed,
        "dangerous": final["dangerous"],
        "tumor_cells": final["tumor_cells"],
        "queens": final["queens"],
        "healthy_destroyed": sum(healthy),
        "rescue_day": rescue_day,
        "untouched": not healthy,
        "series": series,
        "log": log,
        "seconds": round(time.time() - t0, 1),
    }


def default_processes() -> int:
    """Deja núcleos libres: nunca más de la mitad menos uno (el ordenador debe seguir usable)."""
    return max(1, min(6, (os.cpu_count() or 2) // 2 - 1))


def run_experiment(
    exp: Experiment,
    processes: int,
    on_result: Callable[[dict, int, int], None] | None = None,
    cancelled: Callable[[], bool] = lambda: False,
) -> list[dict]:
    """Ejecuta todas las pacientes de todos los brazos. Llama a on_result(resultado, hechas, total)."""
    seeds = list(range(exp.first_seed, exp.first_seed + exp.patients))
    tasks = [(i, s) for i in range(len(exp.arms)) for s in seeds]
    data = exp.to_dict()
    results: list[dict] = []
    if processes <= 1:
        for i, s in tasks:
            if cancelled():
                break
            results.append(run_patient(data, i, s))
            if on_result:
                on_result(results[-1], len(results), len(tasks))
        return results
    with ProcessPoolExecutor(max_workers=processes) as pool:
        futures = {pool.submit(run_patient, data, i, s) for i, s in tasks}
        while futures:
            done, futures = wait(futures, timeout=0.5, return_when=FIRST_COMPLETED)
            for f in done:
                results.append(f.result())
                if on_result:
                    on_result(results[-1], len(results), len(tasks))
            if cancelled():
                for f in futures:
                    f.cancel()
                break
    return results


def summarize(exp: Experiment, results: list[dict]) -> list[dict]:
    """Una fila por brazo con las cifras que importan."""
    rows = []
    for i, arm in enumerate(exp.arms):
        rs = sorted((r for r in results if r["arm"] == i), key=lambda r: r["seed"])
        if not rs:
            continue
        dang = [r["dangerous"] for r in rs]
        n = len(rs)
        rows.append(
            {
                "arm": arm.name,
                "n": n,
                "controlled": sum(d < exp.control_cells for d in dang),
                "eradicated": sum(d == 0 for d in dang),
                "median_dangerous": statistics.median(dang),
                "median_healthy_destroyed": statistics.median(r["healthy_destroyed"] for r in rs),
                "untouched": sum(r["untouched"] for r in rs),
                "per_patient": dang,
            }
        )
    return rows


def median_curves(
    results: list[dict], arm_index: int
) -> tuple[list[float], list[float], list[float], list[float]]:
    """Mediana, mínimo y máximo de células peligrosas por día en un brazo."""
    rs = [r for r in results if r["arm"] == arm_index]
    if not rs:
        return [], [], [], []
    days = rs[0]["series"]["day"]
    cols = list(zip(*(r["series"]["dangerous"] for r in rs)))
    return days, [statistics.median(c) for c in cols], [min(c) for c in cols], [max(c) for c in cols]
