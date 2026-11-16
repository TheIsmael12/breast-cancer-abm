"""Plantillas: experimentos completos y brazos listos para añadir con un clic."""

from __future__ import annotations

from collections.abc import Callable

from breast_cancer_abm.catalog import ARSENAL
from breast_cancer_abm.lab.experiment import Arm, Experiment, Intervention

I = Intervention  # noqa: E741 - abreviatura legible en las tablas de abajo


def arm_templates(day: float) -> dict[str, Callable[[], Arm]]:
    """Brazos típicos, todos tratando el día ``day``."""
    return {
        "Sin tratar (control)": lambda: Arm("Sin tratar"),
        "Arsenal celular (troyano + guardiana + memoria)": lambda: Arm(
            "Arsenal celular", [I("therapy", day, therapy=ARSENAL, dose=300)]
        ),
        "Cirugía + arsenal": lambda: Arm(
            "Cirugía + arsenal", [I("resect", day), I("therapy", day, therapy=ARSENAL, dose=300)]
        ),
        "Ablación guiada + arsenal": lambda: Arm(
            "Ablación guiada + arsenal", [I("ablate", day, guided=True), I("therapy", day, therapy=ARSENAL)]
        ),
        "Ablación en disco + arsenal": lambda: Arm(
            "Ablación en disco + arsenal",
            [I("ablate", day, guided=False, rim=2), I("therapy", day, therapy=ARSENAL)],
        ),
        "ctDNA mensual → ablación guiada (sin operar)": lambda: Arm(
            "ctDNA → ablación guiada",
            [
                I("therapy", day, therapy=ARSENAL),
                I("realistic_rate", day),
                I("ctdna_monitor", day + 30, every=30, rescue="guided", therapy=ARSENAL),
            ],
        ),
        "Quimioterapia": lambda: Arm("Quimioterapia", [I("therapy", day, therapy="chemo", dose=80)]),
        "Inmunoterapia anti-PD-1": lambda: Arm("Anti-PD-1", [I("block_pd1", day, duration=60)]),
        "Linfocito cazador (CAR-T)": lambda: Arm(
            "CAR-T cazador", [I("therapy", day, therapy="hunter", dose=40)]
        ),
        "Brazo vacío (lo diseño yo)": lambda: Arm("Mi brazo"),
    }


def _no_surgery() -> Experiment:
    t = arm_templates(120)
    return Experiment(
        name="¿Se puede tratar sin operar?",
        question="Comparar la cirugía con la ablación guiada y la vigilancia por ctDNA (avance del 05-11).",
        days=300,
        patients=6,
        arms=[
            t["Sin tratar (control)"](),
            t["Arsenal celular (troyano + guardiana + memoria)"](),
            t["Cirugía + arsenal"](),
            t["ctDNA mensual → ablación guiada (sin operar)"](),
        ],
    )


def _immunotherapy() -> Experiment:
    t = arm_templates(120)
    return Experiment(
        name="¿Funciona el anti-PD-1?",
        question="Con la inmunidad adaptativa activada, ¿mejora el anti-PD-1 a las pacientes?",
        days=260,
        patients=6,
        params={"prob_neoantigen": 0.1, "adaptive_immunity": True},
        arms=[t["Sin tratar (control)"](), t["Inmunoterapia anti-PD-1"]()],
    )


def _cell_therapies() -> Experiment:
    day = 120
    arms = [Arm("Sin tratar")] + [
        Arm(name, [I("therapy", day, therapy=key, dose=dose)])
        for name, key, dose in (
            ("CAR-T cazador", "hunter", 40),
            ("Híbrido repara o elimina", "hybrid", 40),
            ("Arreglar la reina", "queen", 40),
            ("Guardiana + memoria", "niche_guardian+memory", 40),
            ("Arsenal con troyano", ARSENAL, 300),
        )
    ]
    return Experiment(
        name="Comparar terapias celulares",
        question="¿Qué célula programada controla mejor el tumor sin dañar tejido sano?",
        days=260,
        patients=4,
        arms=arms,
    )


def _3d() -> Experiment:
    t = arm_templates(60)
    return Experiment(
        name="3D: ¿ablación guiada o en disco?",
        model="3d",
        size3d="fast",
        question="En el cilindro 3D los fármacos no llegan al núcleo: ¿qué tratamiento local funciona?",
        days=140,
        patients=4,
        arms=[
            t["Sin tratar (control)"](),
            t["Arsenal celular (troyano + guardiana + memoria)"](),
            t["Ablación en disco + arsenal"](),
            t["ctDNA mensual → ablación guiada (sin operar)"](),
        ],
    )


EXPERIMENT_TEMPLATES: dict[str, tuple[str, Callable[[], Experiment]]] = {
    "no_surgery": ("¿Se puede tratar sin operar?", _no_surgery),
    "immunotherapy": ("¿Funciona el anti-PD-1?", _immunotherapy),
    "cell_therapies": ("Comparar terapias celulares", _cell_therapies),
    "3d": ("3D: ablación guiada o en disco", _3d),
    "blank": ("Empezar de cero", lambda: Experiment(arms=[Arm("Sin tratar")])),
}
