"""Simulación 2D en vivo: el corte de un conducto mamario, célula a célula."""

from __future__ import annotations

import asyncio

from nicegui import run, ui

from breast_cancer_abm.app.components.tissue_canvas import TissueCanvas
from breast_cancer_abm.app.theme import card, fmt, ghost, kpi, legend_item, primary, secondary, shell
from breast_cancer_abm.catalog import THERAPIES as THERAPY_NAMES
from breast_cancer_abm.catalog import therapy_label
from breast_cancer_abm.model.biology import STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.policies import THERAPIES
from breast_cancer_abm.model.tissue import (
    HYPERPLASIA,
    IN_SITU,
    INITIATED,
    INVASIVE,
    NECROSIS,
    QUEEN,
    REPAIRED,
    STATE_NAME,
    Simulation,
)
from breast_cancer_abm.palette import AGENT_STYLE, STATE_COLORS, VESSEL_COLOR, legend

SCENARIO_LABELS = {
    "luminal": "Luminal ER+ (~70 %)",
    "her2": "HER2+",
    "basal": "Triple negativo / basal",
    "brca1": "Portadora de BRCA1",
}
VESSEL_MODES = {
    "normal": "Normales",
    "reroute": "Reconducir el O₂ hacia el tumor",
    "normalize": "Normalizar",
    "antiangiogenic": "Antiangiogénico",
    "starve": "Asfixiar",
}
THERAPY_OPTIONS = {k: v[0] for k, v in THERAPY_NAMES.items()}
CURVES = (INITIATED, HYPERPLASIA, IN_SITU, INVASIVE, QUEEN, NECROSIS, REPAIRED)
CHART_COLOR = {**STATE_COLORS, QUEEN: "#111827"}  # la reina es blanca en el lienzo, no en la gráfica
RESULT_LABELS = {
    "excised": "células extirpadas",
    "destroyed": "tumorales destruidas",
    "healthy": "sanas destruidas",
    "residual": "carcinoma residual",
    "foci_seen": "focos vistos",
    "positive": "positiva",
    "sampled": "células muestreadas",
    "real": "carcinoma real en el tejido",
    "fragments": "fragmentos detectados",
    "expected": "fragmentos esperados",
    "alive": "células tumorales vivas",
}


class Session:
    """Estado de una pestaña: la simulación y lo que se está mostrando."""

    def __init__(self) -> None:
        self.scenario = "luminal"
        self.seed: int | None = 6
        self.sim = self._new()
        self.playing = False
        self.speed = 8
        self.oxygen = False
        self.agents = True
        self.selection: tuple[int, int] | None = None
        self.lock = asyncio.Lock()
        self.ticks = 0

    def _new(self) -> Simulation:
        return Simulation(self.scenario, self.seed, Parameters(prob_neoantigen=0.1))

    def reset(self) -> None:
        self.sim = self._new()
        self.selection = None

    def advance(self, steps: int) -> None:
        for _ in range(steps):
            self.sim.step()


def page() -> None:
    s = Session()

    with shell(
        "/simulation",
        "Simulación 2D",
        "Corte de un conducto mamario: cada píxel es una célula (~15 µm) y cada paso son 6 horas.",
    ) as pg:
        with ui.dialog() as help_dialog, ui.card().classes("q-pa-lg").style("max-width: 560px"):
            ui.label("Cómo se usa").classes("text-lg font-semibold")
            ui.markdown(
                "1. Pulsa **Reproducir** y mira cómo nace el tumor (suele aparecer hacia el día 60-120).\n"
                "2. Cuando quieras, **libera una terapia** o haz una **intervención** desde el panel derecho.\n"
                "3. **Haz clic en una célula** del tejido para ver su genoma.\n\n"
                "Para comparar tratamientos con muchas pacientes a la vez, usa el **Laboratorio**."
            )
            with ui.row().classes("w-full justify-end"):
                primary("Entendido", help_dialog.close)
        with pg.actions:
            ghost("Cómo se usa", lambda: help_dialog.open(), icon="help_outline")

        # -------------------------------------------------------------- barra
        with ui.element("div").classes("toolbar"):
            with ui.element("div").classes("toolbar-group"):
                scenario = (
                    ui.select(SCENARIO_LABELS, value=s.scenario, label="Tipo de cáncer")
                    .props("outlined dense options-dense")
                    .style("width: 230px")
                )
                seed = (
                    ui.number("Paciente (semilla)", value=s.seed, min=0, step=1, format="%d")
                    .props("outlined dense")
                    .style("width: 150px")
                    .tooltip("Cada semilla es una paciente distinta; vacía = al azar")
                )
                secondary("Nueva", lambda: new_simulation(), icon="restart_alt")
            ui.element("div").classes("toolbar-sep")
            with ui.element("div").classes("toolbar-group"):
                play = primary("Reproducir", lambda: toggle_play(), icon="play_arrow").style("width: 140px")
                for label, steps, tip in (
                    ("6 h", 1, "Un paso"),
                    ("+1 día", STEPS_PER_DAY, "Avanzar un día"),
                    ("+10 días", 10 * STEPS_PER_DAY, "Avanzar diez días"),
                ):
                    secondary(label, lambda _, n=steps: step(n)).tooltip(tip)
            ui.element("div").classes("toolbar-sep")
            with ui.element("div").classes("toolbar-group"):
                ui.label("Velocidad").classes("field-label")
                ui.slider(min=1, max=40, step=1, value=s.speed).bind_value(s, "speed").style("width: 130px")

        with ui.element("div").classes("sim-grid"):
            # ---------------------------------------------------------- tejido
            with ui.column().classes("gap-0 w-full"):
                layers = []

                def layer_switches() -> None:
                    layers.append(ui.switch("Oxígeno").bind_value(s, "oxygen").props("dense"))
                    layers.append(ui.switch("Inmunes").bind_value(s, "agents").props("dense"))

                with card("Tejido", icon="blur_on", actions=layer_switches):
                    with ui.row().classes("w-full items-baseline justify-between no-wrap"):
                        title = ui.label().classes("text-base font-semibold")
                        status = ui.label().classes("hint")
                    with ui.element("div").classes("sim-stage"):
                        canvas = TissueCanvas(on_cell=lambda y, x: select_cell(y, x))
                    with ui.element("div").classes("legend"):
                        for color, label in legend():
                            legend_item(color, label)
                        legend_item(VESSEL_COLOR, "vaso sanguíneo")
                        for _, color, name in AGENT_STYLE.values():
                            legend_item(color, name)

            # ---------------------------------------------------------- panel
            with ui.column().classes("gap-4 w-full"):
                with ui.grid(columns=3).classes("w-full gap-2 kpi-compact"):
                    k_day, k_tumor, k_danger, k_queen, k_immune, k_ai = (
                        kpi(label, accent=accent)
                        for label, accent in (
                            ("Día", False),
                            ("Tumorales", False),
                            ("Peligrosas", True),
                            ("Reinas", False),
                            ("Inmunes", False),
                            ("Terapéuticas", False),
                        )
                    )
                with card(
                    "Tratamiento", "Libera células terapéuticas o actúa sobre el tejido.", icon="vaccines"
                ):
                    therapy = (
                        ui.select(
                            THERAPY_OPTIONS,
                            value=["trojan", "niche_guardian", "memory"],
                            label="Terapia celular (una o varias, se combinan)",
                            multiple=True,
                            with_input=True,
                        )
                        .props("outlined dense use-chips options-dense")
                        .classes("w-full")
                    )
                    description = ui.label().classes("hint")
                    with ui.row().classes("w-full no-wrap gap-2 items-center"):
                        dose = (
                            ui.number("Dosis (células)", value=300, min=1, step=10, format="%d")
                            .props("outlined dense")
                            .classes("col")
                        )
                        primary("Liberar terapia", lambda: release(), icon="vaccines")
                    ui.element("div").classes("divider")
                    ui.label("Intervenciones").classes("field-label")
                    with ui.grid(columns=2).classes("w-full gap-2"):
                        for label, icon, tip, fn in (
                            (
                                "Ablación guiada",
                                "my_location",
                                "Destruye solo los focos visibles y un margen",
                                lambda x: x.ablate(guided=True, rim=1),
                            ),
                            (
                                "Ablación disco",
                                "radio_button_checked",
                                "Destruye un disco alrededor del tumor",
                                lambda x: x.ablate(guided=False, rim=2),
                            ),
                            ("Cirugía", "content_cut", "Extirpa la masa con margen", lambda x: x.resect()),
                            (
                                "Biopsia",
                                "colorize",
                                "12 cilindros: ¿queda carcinoma?",
                                lambda x: x.biopsy(12),
                            ),
                            (
                                "ctDNA",
                                "bloodtype",
                                "Análisis de sangre: ADN tumoral",
                                lambda x: x.ctdna(0.05),
                            ),
                            (
                                "Anti-PD-1",
                                "shield",
                                "Inmunoterapia durante 60 días",
                                lambda x: x.block_pd1(60),
                            ),
                        ):
                            secondary(label, lambda _, f=fn, n=label: intervene(n, f), icon=icon).props(
                                "dense"
                            ).classes("w-full").tooltip(tip)
                    ui.select(
                        VESSEL_MODES,
                        value="normal",
                        label="Vasos sanguíneos",
                        on_change=lambda e: vessels(e.value),
                    ).props("outlined dense").classes("w-full")

                with card(
                    "Evolución",
                    "Células por estado a lo largo del tiempo (escala logarítmica).",
                    icon="show_chart",
                ):
                    chart = ui.echart(_chart_options()).classes("w-full").style("height: 260px")

                with card("Célula seleccionada", icon="ads_click"):
                    inspector = ui.column().classes("gap-0 w-full")
                with card("Genes alterados en el tumor", icon="science"):
                    genes = ui.column().classes("gap-1 w-full")

                with card("Eventos", icon="history"):
                    events = ui.column().classes("gap-0 w-full")

    # ---------------------------------------------------------------- lógica
    def describe_therapy() -> None:
        names = therapy.value or []
        description.set_text(
            " · ".join(THERAPIES[n]["description"].split("\n")[0] for n in names if n in THERAPIES)
        )

    therapy.on_value_change(lambda _: describe_therapy())
    describe_therapy()

    def redraw(full: bool = True) -> None:
        sim = s.sim
        canvas.draw(sim, oxygen=s.oxygen, agents=s.agents, selection=s.selection)
        r = sim.summary()
        kinds = [a.kind for a in sim.agents]
        k_day.set_text(fmt(round(sim.day, 1)))
        k_tumor.set_text(fmt(r["tumor_cells"]))
        k_danger.set_text(fmt(r["dangerous"]))
        k_queen.set_text(fmt(r["queens"]))
        k_immune.set_text(fmt(kinds.count("T") + kinds.count("NK")))
        k_ai.set_text(fmt(kinds.count("AI")))
        hours = (sim.t % STEPS_PER_DAY) * 24 // STEPS_PER_DAY
        title.set_text(f"Día {int(sim.day)} · {hours:02d}:00 h")
        status.set_text(
            f"{SCENARIO_LABELS[sim.scenario]} · semilla {s.seed if s.seed is not None else 'al azar'}"
            f" · vasos: {VESSEL_MODES.get(sim.vessel_mode, sim.vessel_mode).lower()}"
        )
        if full:
            update_chart()
            update_panels()

    def update_chart() -> None:
        h = s.sim.history
        days = h["day"]
        k = max(1, len(days) // 400)
        idx = range(0, len(days), k)
        series = chart.options["series"]
        for i, code in enumerate(CURVES):
            series[i]["data"] = [[days[j], h[code][j] or None] for j in idx]
        series[len(CURVES)]["data"] = [[days[j], (h["T"][j] + h["NK"][j]) or None] for j in idx]
        series[len(CURVES) + 1]["data"] = [[days[j], h["AI"][j] or None] for j in idx]
        chart.update()

    def update_panels() -> None:
        sim = s.sim
        inspector.clear()
        with inspector:
            if s.selection:
                for i, line in enumerate(sim.describe(s.selection)):
                    ui.label(line).classes("text-sm" + (" font-semibold" if i == 0 else " muted"))
            else:
                ui.label("Haz clic en una célula del tejido para ver su genoma.").classes("text-sm muted")
        genes.clear()
        with genes:
            freqs = sorted(((v, g) for g, v in sim.tumor_frequencies().items() if v > 0), reverse=True)
            if not freqs:
                ui.label("Todavía no hay tumor.").classes("text-sm muted")
            for v, g in freqs[:8]:
                with ui.row().classes("w-full items-center no-wrap gap-2"):
                    ui.label(g).classes("mono").style("width: 64px")
                    ui.linear_progress(value=v, show_value=False, size="8px").props(
                        "rounded color=primary track-color=grey-3"
                    ).classes("col")
                    ui.label(f"{v:.0%}").classes("text-xs muted").style("width: 34px; text-align: right")
        events.clear()
        with events:
            for e in list(sim.events)[-8:][::-1]:
                ui.label(e).classes("event-line")

    def select_cell(y: int, x: int) -> None:
        s.selection = (y, x)
        redraw()

    async def step(steps: int) -> None:
        async with s.lock:
            await run.io_bound(s.advance, steps)
        redraw()

    async def tick() -> None:
        if not s.playing or s.lock.locked():
            return
        async with s.lock:
            await run.io_bound(s.advance, int(s.speed))
        s.ticks += 1
        redraw(full=s.ticks % 4 == 0)

    def toggle_play() -> None:
        s.playing = not s.playing
        play.set_text("Pausa" if s.playing else "Reproducir")
        play.props(f"icon={'pause' if s.playing else 'play_arrow'}")
        if not s.playing:
            redraw()

    async def new_simulation() -> None:
        async with s.lock:
            s.scenario = scenario.value
            s.seed = None if seed.value is None else int(seed.value)
            await run.io_bound(s.reset)
        canvas.set_tissue(s.sim)
        redraw()

    async def release() -> None:
        names = therapy.value or []
        if not names:
            ui.notify("Elige al menos una terapia.", type="warning")
            return
        async with s.lock:
            s.sim.release_therapy("+".join(names), int(dose.value) if dose.value else None)
        ui.notify(
            f"Liberada: {therapy_label('+'.join(names))} ({int(dose.value or 0)} células)", type="positive"
        )
        redraw()

    async def intervene(name: str, fn) -> None:
        async with s.lock:
            result = fn(s.sim)
        if isinstance(result, dict):
            parts = [
                f"{RESULT_LABELS.get(k, k)}: {('sí' if v else 'no') if isinstance(v, bool) else fmt(v)}"
                for k, v in result.items()
            ]
            ui.notify(f"{name} · " + " · ".join(parts), multi_line=True, timeout=8000)
        else:
            ui.notify(f"{name} aplicado.", type="positive")
        redraw()

    async def vessels(mode: str) -> None:
        async with s.lock:
            s.sim.treat_vessels(mode)
        redraw()

    for switch in layers:
        switch.on_value_change(lambda _: redraw())
    ui.timer(0.06, tick)

    async def first_paint() -> None:
        canvas.set_tissue(s.sim)
        redraw()

    ui.timer(0.2, first_paint, once=True)


def _chart_options() -> dict:
    series = [
        {
            "name": STATE_NAME[c],
            "type": "line",
            "showSymbol": False,
            "data": [],
            "lineStyle": {"width": 1.8},
            "color": CHART_COLOR[c],
            "connectNulls": False,
        }
        for c in CURVES
    ]
    series += [
        {
            "name": "inmunes naturales",
            "type": "line",
            "showSymbol": False,
            "data": [],
            "lineStyle": {"width": 1.2, "type": "dashed"},
            "color": "#0ea5e9",
        },
        {
            "name": "terapéuticas",
            "type": "line",
            "showSymbol": False,
            "data": [],
            "lineStyle": {"width": 1.4, "type": "dashed"},
            "color": AGENT_STYLE["AI"][1],
        },
    ]
    axis = {
        "axisLine": {"lineStyle": {"color": "#d1d5db"}},
        "axisLabel": {"color": "#6b7280", "fontSize": 10},
        "splitLine": {"lineStyle": {"color": "#f3f4f6"}},
    }
    return {
        "animation": False,
        "grid": {"left": 40, "right": 10, "top": 56, "bottom": 26},
        "legend": {
            "type": "scroll",
            "top": 0,
            "icon": "roundRect",
            "itemWidth": 10,
            "itemHeight": 10,
            "textStyle": {"fontSize": 10, "color": "#374151"},
        },
        "tooltip": {"trigger": "axis"},
        "xAxis": {"type": "value", "name": "días", "nameGap": 4, "nameTextStyle": {"fontSize": 10}, **axis},
        "yAxis": {"type": "log", "min": 1, "logBase": 10, **axis},
        "series": series,
    }
