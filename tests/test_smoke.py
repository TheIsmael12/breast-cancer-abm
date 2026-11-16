"""Prueba de humo: que el motor arranca y que cada intervención devuelve lo que promete.

No mide resultados científicos (eso es ``bcabm study validation``, que tarda): solo que
nada se rompe. La huella exacta del modelo está en test_fingerprint.py.

    pytest
"""

from __future__ import annotations

from breast_cancer_abm.model.biology import STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation


def _sim(days: int = 60, seed: int = 6) -> Simulation:
    sim = Simulation("luminal", seed, Parameters(prob_neoantigen=0.1))
    for _ in range(days * STEPS_PER_DAY):
        sim.step()
    return sim


def test_reproducible():
    assert _sim(30, 3).summary() == _sim(30, 3).summary()


def test_tumor_appears():
    assert _sim(90).summary()["tumor_cells"] > 0


def test_local_interventions():
    base = _sim(90)
    for act in ("resect", "disc", "guided"):
        import copy

        sim = copy.deepcopy(base)
        if act == "resect":
            r = sim.resect()
        else:
            r = sim.ablate(guided=act == "guided", rim=1)
        assert set(("healthy", "residual")) <= set(r)
        assert r["healthy"] >= 0 and r["residual"] >= 0
        for _ in range(2 * STEPS_PER_DAY):
            sim.step()


def test_residual_disease_tests():
    sim = _sim(90)
    b = sim.biopsy(6)
    assert b["sampled"] > 0 and b["real"] >= 0
    c = sim.ctdna(0.1)
    assert c["fragments"] >= 0 and c["expected"] >= 0


def test_therapy_and_3d():
    sim = _sim(40)
    sim.release_therapy("trojan+niche_guardian+memory", 30)
    for _ in range(STEPS_PER_DAY):
        sim.step()
    from breast_cancer_abm.model.tissue3d import Simulation3D

    s3 = Simulation3D("luminal", 1, side=30, z_len=12, radius=8.0)
    for _ in range(STEPS_PER_DAY):
        s3.step()
    assert s3.summary()["normal_cells"] > 0
