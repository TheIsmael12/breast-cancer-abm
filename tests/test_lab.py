"""El laboratorio: plantillas válidas, experimentos que se guardan y se ejecutan."""

from __future__ import annotations

import dataclasses

from breast_cancer_abm.catalog import PARAMS
from breast_cancer_abm.lab.experiment import Arm, Experiment, Intervention, run_experiment, summarize
from breast_cancer_abm.lab.templates import EXPERIMENT_TEMPLATES, arm_templates
from breast_cancer_abm.model.biology import Parameters


def test_catalog_covers_every_parameter():
    fields = {f.name for f in dataclasses.fields(Parameters)} - {"grid_size", "duct_radius"}
    assert fields == set(PARAMS)


def test_templates_round_trip():
    for _, make in EXPERIMENT_TEMPLATES.values():
        exp = make()
        assert Experiment.from_dict(exp.to_dict()).to_dict() == exp.to_dict()
    assert all(callable(f) for f in arm_templates(100).values())


def test_small_experiment_runs():
    exp = Experiment(
        days=40,
        patients=1,
        arms=[
            Arm("control"),
            Arm(
                "tratado",
                [Intervention("therapy", 20, therapy="hunter", dose=20), Intervention("ablate", 30)],
            ),
            Arm("ctdna", [Intervention("ctdna_monitor", 10, every=10, rescue="guided", fragments=1.0)]),
        ],
    )
    res = run_experiment(exp, processes=1)
    rows = summarize(exp, res)
    assert [r["arm"] for r in rows] == ["control", "tratado", "ctdna"]
    assert all(r["n"] == 1 for r in rows)
    assert len(res[0]["series"]["day"]) == 41
    assert any("Ablación" in text for _, text in res[1]["log"])
