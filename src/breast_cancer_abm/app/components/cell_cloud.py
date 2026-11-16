"""Nube de células 3D para ``ui.scene`` y su empaquetado desde una ``Simulation3D``."""

from __future__ import annotations

import base64

import numpy as np
from nicegui.elements.scene.scene_object3d import Object3D

from breast_cancer_abm.model.tissue import (
    FIBROBLAST,
    HYPERPLASIA,
    IN_SITU,
    INITIATED,
    INVASIVE,
    MYOEPITHELIAL,
    NECROSIS,
    NORMAL,
    QUEEN,
    TUMOR_STATES,
)
from breast_cancer_abm.model.tissue3d import MYO, STROMA, Simulation3D
from breast_cancer_abm.palette import AGENT_STYLE, STATE_COLORS

# Códigos de color propios de la vista 3D (índices de PALETTE_3D).
MYOEP_3D, AGENT_3D = 11, 12
PALETTE_3D = [STATE_COLORS.get(i, "#000000") for i in range(11)] + ["#8e7cc3", AGENT_STYLE["AI"][1]]


def state3d(sim: Simulation3D, pos: tuple, c) -> int:
    """Mismo código de colores que el visor 2D (-1 = fibroblasto, que no se pinta)."""
    if c.kind == FIBROBLAST:
        return -1
    if c.kind == MYOEPITHELIAL:
        return MYOEP_3D
    if c.necrosis:
        return NECROSIS
    if c.n_drivers == 0:
        return NORMAL
    if c.stem:
        return QUEEN
    if sim.zone[pos] in (MYO, STROMA):
        return INVASIVE
    if c.anoikis_res:
        return IN_SITU
    return HYPERPLASIA if c.prolif_signal else INITIATED


def pack(sim: Simulation3D, *, cutaway: bool, tumor_only: bool) -> tuple[dict, dict, dict]:
    """Devuelve (células, agentes, recuento por estado) listos para el cliente."""
    axis_mid = sim.center
    positions, states = [], []
    counts: dict[int, int] = {}
    for pos, c in sim.cells.items():
        z, y, x = pos
        if cutaway and x > axis_mid + 0.5:
            continue
        s = state3d(sim, pos, c)
        if s < 0 or (tumor_only and s not in (*TUMOR_STATES, NECROSIS)):
            continue
        positions.append(pos)
        states.append(s)
        counts[s] = counts.get(s, 0) + 1
    agents = [a.pos for a in sim.agents if a.kind == "AI" and (not cutaway or a.pos[2] <= axis_mid + 0.5)]
    zl, side, _ = sim.grid_shape
    center = [zl / 2, (side - 1) / 2, (side - 1) / 2]
    return (
        _payload(positions, states, center),
        _payload(agents, [AGENT_3D] * len(agents), center),
        counts,
    )


def _payload(positions: list, states: list, center: list) -> dict:
    pos = np.asarray(positions, dtype=np.uint8).reshape(-1, 3)
    return {
        "positions": base64.b64encode(pos.tobytes()).decode("ascii"),
        "states": base64.b64encode(np.asarray(states, dtype=np.uint8).tobytes()).decode("ascii"),
        "center": center,
    }


class CellCloud(Object3D, component="cell_cloud.js"):
    def __init__(self, payload: dict, size: float = 1.0) -> None:
        super().__init__(payload, PALETTE_3D, size)

    def set_cells(self, payload: dict) -> None:
        self.args[0] = payload  # para reenviarlo si el cliente se reconecta
        self.run_method("set_cells", payload)

    def set_size(self, size: float) -> None:
        self.args[2] = size
        self.run_method("set_size", size)
