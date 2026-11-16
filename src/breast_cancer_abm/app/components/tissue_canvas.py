"""Componente del lienzo 2D: convierte una ``Simulation`` en fotogramas compactos."""

from __future__ import annotations

import base64
from collections.abc import Callable

import numpy as np
from nicegui.element import Element
from nicegui.events import GenericEventArguments

from breast_cancer_abm.model.tissue import LUMEN, STROMA, Simulation
from breast_cancer_abm.palette import AGENT_STYLE, VESSEL_COLOR, ZONE_BACKGROUND, state_palette

_SHAPES = {"T": "triangle", "NK": "diamond", "AI": "star", "MAC": "circle"}


def _b64(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a, dtype=np.uint8).tobytes()).decode("ascii")


class TissueCanvas(Element, component="tissue_canvas.js", default_classes="tissue-canvas"):
    def __init__(self, on_cell: Callable[[int, int], None] | None = None) -> None:
        super().__init__()
        self._props["palette"] = state_palette()
        self._props["background"] = {"lumen": ZONE_BACKGROUND[LUMEN], "stroma": ZONE_BACKGROUND[STROMA]}
        self._props["vessel-color"] = VESSEL_COLOR
        self._props["agent-style"] = {k: {"shape": _SHAPES[k], "color": v[1]} for k, v in AGENT_STYLE.items()}
        if on_cell:
            self.on("cell", lambda e: on_cell(*_cell_args(e)))

    def set_tissue(self, sim: Simulation) -> None:
        """Fija la geometría (zonas) de una simulación nueva."""
        self.run_method("setZone", int(sim.p.grid_size), _b64(sim.zone))

    def draw(self, sim: Simulation, *, oxygen: bool, agents: bool, selection: tuple[int, int] | None) -> None:
        frame = {
            "states": _b64(sim.code),
            "vessels": _b64(sim.vessels),
            "o2": _b64(np.clip(sim.O2 * 255, 0, 255)) if oxygen else None,
            "agents": [[a.kind, int(a.pos[0]), int(a.pos[1])] for a in sim.agents] if agents else [],
            "selection": list(selection) if selection else None,
        }
        self.run_method("draw", frame)


def _cell_args(e: GenericEventArguments) -> tuple[int, int]:
    return int(e.args["y"]), int(e.args["x"])
