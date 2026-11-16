"""Huella reducida del modelo: simulaciones deterministas cortas y lo que devuelven.

La usa ``tests/test_fingerprint.py``. Si un cambio la mueve, cambia la ciencia: o es un
error, o es un cambio de modelo que hay que declarar (y regenerar la huella con
``python tests/fingerprint.py``, explicándolo en la bitácora).
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from breast_cancer_abm.model.biology import STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation
from breast_cancer_abm.model.tissue3d import Simulation3D

FILE = Path(__file__).parent / "data" / "fingerprint.json"


def _advance(sim, days: float) -> None:
    for _ in range(int(days * STEPS_PER_DAY)):
        sim.step()


def _print2d(sim) -> dict:
    return dict(
        summary=sim.summary(),
        t=sim.t,
        code=hashlib.sha256(sim.code.tobytes()).hexdigest()[:16],
        agents=[(a.kind, list(a.pos), a.kills) for a in sim.agents][:30],
    )


def _print3d(sim) -> dict:
    state = sorted((pos, c.n_drivers, c.necrosis, bool(c.stem)) for pos, c in sim.cells.items())
    return dict(summary=sim.summary(), t=sim.t, cells=hashlib.sha256(repr(state).encode()).hexdigest()[:16])


def compute() -> dict:
    out = {}
    sim = Simulation("luminal", 6, Parameters(prob_neoantigen=0.1, editing_error=0.01))
    _advance(sim, 60)
    out["2d_day60"] = _print2d(sim)
    for therapy in ("hunter", "trojan+niche_guardian+memory"):
        t = copy.deepcopy(sim)
        t.release_therapy(therapy, 80)
        _advance(t, 6)
        out[f"therapy_{therapy}"] = _print2d(t)
    for guided in (False, True):
        t = copy.deepcopy(sim)
        r = t.ablate(guided=guided, rim=1 if guided else 2)
        out[f"ablation_guided_{guided}"] = dict(r=r, ctdna=t.ctdna(0.05), biopsy=t.biopsy(12))
    s3 = Simulation3D("luminal", 1, Parameters(prob_neoantigen=0.1), side=30, z_len=12, radius=8.0)
    _advance(s3, 20)
    s3.release_therapy("trojan+niche_guardian+memory", 60)
    _advance(s3, 3)
    out["3d"] = _print3d(s3)
    return json.loads(json.dumps(out, default=str))


if __name__ == "__main__":
    FILE.parent.mkdir(exist_ok=True)
    FILE.write_text(json.dumps(compute(), indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print("huella regenerada en", FILE)
