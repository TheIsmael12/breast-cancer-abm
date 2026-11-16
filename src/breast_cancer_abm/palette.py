"""Colores del modelo, compartidos por el visor de matplotlib, las figuras y la aplicación."""

from __future__ import annotations

from breast_cancer_abm.model.tissue import (
    EMPTY,
    FIBRO,
    HYPERPLASIA,
    IN_SITU,
    INITIATED,
    INVASIVE,
    LUMEN,
    MYOEP,
    NECROSIS,
    NORMAL,
    QUEEN,
    REPAIRED,
    STATE_NAME,
    STROMA,
)

STATE_COLORS: dict[int, str] = {
    FIBRO: "#6b5a4a",
    MYOEP: "#8e7cc3",
    NORMAL: "#f4b6c8",
    REPAIRED: "#3ddc84",
    INITIATED: "#ffe066",
    HYPERPLASIA: "#ff9f1c",
    IN_SITU: "#e63946",
    INVASIVE: "#9d0208",
    NECROSIS: "#4a4a4a",
    QUEEN: "#ffffff",
}
ZONE_BACKGROUND = {LUMEN: "#120c16", STROMA: "#231b17"}
VESSEL_COLOR = "#ff2d55"
# tipo de agente -> (marcador de matplotlib, color, nombre)
AGENT_STYLE = {
    "T": ("^", "#4cc9f0", "linfocito T"),
    "NK": ("D", "#f9f871", "célula NK"),
    "AI": ("*", "#f72585", "célula terapéutica"),
    "MAC": ("o", "#a7f3d0", "macrófago"),
}
# orden de la leyenda (de sano a peor)
LEGEND_ORDER = (NORMAL, INITIATED, HYPERPLASIA, IN_SITU, INVASIVE, QUEEN, NECROSIS, REPAIRED, MYOEP, FIBRO)


def state_palette() -> list[str]:
    """Lista indexada por código de estado (0..10); el vacío no se pinta."""
    return [STATE_COLORS.get(code, "#000000") if code != EMPTY else "" for code in range(11)]


def legend() -> list[tuple[str, str]]:
    return [(STATE_COLORS[c], STATE_NAME[c]) for c in LEGEND_ORDER]
