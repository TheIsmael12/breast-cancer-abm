"""Breast Cancer ABM: simulador por agentes de cáncer de mama y células terapéuticas.

Modelo educativo y de investigación desarrollado por Claude (Anthropic) y dirigido por
Ismael Ben. Las cifras comparan brazos de tratamiento dentro del modelo; no son
predicciones clínicas.

Uso rápido::

    from breast_cancer_abm import Simulation, Parameters

    sim = Simulation("luminal", seed=1, params=Parameters(prob_neoantigen=0.1))
    for _ in range(120 * 4):          # 120 días, 4 pasos de 6 h por día
        sim.step()
    print(sim.summary())
"""

from __future__ import annotations

from breast_cancer_abm.model.biology import CALIBRATIONS, GENES, SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.policies import THERAPIES
from breast_cancer_abm.model.tissue import Simulation
from breast_cancer_abm.model.tissue3d import Simulation3D

__version__ = "1.0.0"

__all__ = [
    "CALIBRATIONS",
    "GENES",
    "SCENARIOS",
    "STEPS_PER_DAY",
    "THERAPIES",
    "Parameters",
    "Simulation",
    "Simulation3D",
    "__version__",
]
