"""Simulación 3D: el conducto como cilindro, con la cámara en la mano."""

from __future__ import annotations

import asyncio
import time

from nicegui import run, ui

from breast_cancer_abm.app.components.cell_cloud import AGENT_3D, MYOEP_3D, PALETTE_3D, CellCloud, pack
from breast_cancer_abm.app.pages.simulation2d import SCENARIO_LABELS, THERAPY_OPTIONS
from breast_cancer_abm.app.theme import (
    CANVAS_BG,
    card,
    fmt,
    ghost,
    kpi,
    legend_item,
    primary,
    secondary,
    shell,
)
from breast_cancer_abm.catalog import therapy_label
from breast_cancer_abm.model.biology import CALIBRATIONS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import (
    HYPERPLASIA,
    IN_SITU,
    INITIATED,
    INVASIVE,
    NECROSIS,
    NORMAL,
    QUEEN,
    STATE_NAME,
)
from breast_cancer_abm.model.tissue3d import Simulation3D

SIZES = {
    "fast": ("Rápido (34 × 34 × 16)", dict(side=34, z_len=16, radius=9.0)),
    "standard": ("Estándar (54 × 54 × 40, el de los estudios)", dict(side=54, z_len=40, radius=13.0)),
}
LEGEND_3D = (NORMAL, INITIATED, HYPERPLASIA, IN_SITU, INVASIVE, QUEEN, NECROSIS, MYOEP_3D, AGENT_3D)
LABEL_3D = {**STATE_NAME, MYOEP_3D: "mioepitelial", AGENT_3D: "célula terapéutica"}


class Session:
    def __init__(self) -> None:
        self.scenario = "luminal"
        self.seed: int | None = 1
        self.size = "fast"
        self.sim = self._new()
        self.cutaway = True
        self.tumor_only = False
        self.point_size = 1.3
        self.lock = asyncio.Lock()
        self.cancel = False
        self.busy = False

    def _new(self) -> Simulation3D:
        return Simulation3D(
            self.scenario,
            self.seed,
            Parameters(prob_neoantigen=0.1, editing_error=0.01),
            **SIZES[self.size][1],
        )

    def reset(self) -> None:
        self.sim = self._new()


def page() -> None:
    s = Session()

    with shell(
        "/simulation-3d",
        "Simulación 3D",
        "El conducto como un cilindro hueco (26 vecinos por célula). Aquí se ve el muro que frena a las "
        "terapias: solo llegan a las células que tienen un hueco al lado.",
    ):
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
                    .style("width: 140px")
                )
                size = (
                    ui.select({k: v[0] for k, v in SIZES.items()}, value=s.size, label="Tamaño")
                    .props("outlined dense")
                    .style("width: 200px")
                )
                secondary("Nueva", lambda: new_simulation(), icon="restart_alt")
            ui.element("div").classes("toolbar-sep")
            with ui.element("div").classes("toolbar-group"):
                for days in (1, 10, 30):
                    secondary(f"+{days} {'día' if days == 1 else 'días'}", lambda _, d=days: advance(d))
                stop = ghost("Detener", lambda: setattr(s, "cancel", True), icon="stop")
                stop.visible = False
                progress = (
                    ui.linear_progress(value=0, show_value=False, size="6px")
                    .props("rounded track-color=grey-3")
                    .style("width: 120px")
                )
                progress.visible = False

        with ui.element("div").classes("sim-grid"):
            # ---------------------------------------------------------- escena
            with ui.column().classes("gap-0 w-full"):

                def camera_buttons() -> None:
                    ghost("Lateral", lambda: camera("side"))
                    ghost("Frontal", lambda: camera("front"))
                    ghost("Perspectiva", lambda: camera("iso"))

                with card(
                    "Conducto",
                    "Arrastra para girar y usa la rueda para acercar.",
                    icon="view_in_ar",
                    actions=camera_buttons,
                ):
                    with ui.row().classes("w-full items-baseline justify-between no-wrap"):
                        title = ui.label().classes("text-base font-semibold")
                        status = ui.label().classes("hint")
                    with ui.element("div").classes("scene-wrap"):
                        scene = ui.scene(
                            width=980, height=600, grid=False, background_color=CANVAS_BG
                        ).classes("w-full")
                        with scene:
                            cells_payload, agents_payload, _ = pack(s.sim, cutaway=True, tumor_only=False)
                            cloud = CellCloud(cells_payload, s.point_size)
                            agents_cloud = CellCloud(agents_payload, 2.2)
                    with ui.row().classes("w-full items-center no-wrap gap-3"):
                        view_switches = [
                            ui.switch("Corte longitudinal").bind_value(s, "cutaway").props("dense"),
                            ui.switch("Solo el tumor").bind_value(s, "tumor_only").props("dense"),
                        ]
                        ui.space()
                        ui.label("Tamaño del punto").classes("field-label")
                        point = ui.slider(min=0.4, max=2.4, step=0.1, value=s.point_size).style(
                            "width: 140px"
                        )

            # ---------------------------------------------------------- panel
            with ui.column().classes("gap-4 w-full"):
                with ui.grid(columns=3).classes("w-full gap-2 kpi-compact"):
                    k_day, k_tumor, k_danger, k_queen, k_access, k_ai = (
                        kpi(label, accent=accent)
                        for label, accent in (
                            ("Día", False),
                            ("Tumorales", False),
                            ("Peligrosas", True),
                            ("Reinas", False),
                            ("Accesible", False),
                            ("Terapéuticas", False),
                        )
                    )
                with card(
                    "Tratamiento",
                    "300 células por dosis; desde ese momento, tasa de mutación realista.",
                    icon="vaccines",
                ):
                    therapy = (
                        ui.select(
                            THERAPY_OPTIONS,
                            value=["trojan", "niche_guardian", "memory"],
                            label="Terapia celular",
                            multiple=True,
                            with_input=True,
                        )
                        .props("outlined dense use-chips options-dense")
                        .classes("w-full")
                    )
                    primary("Liberar terapia", lambda: release(), icon="vaccines").classes("w-full")
                    ui.element("div").classes("divider")
                    ui.label("Intervenciones").classes("field-label")
                    with ui.grid(columns=2).classes("w-full gap-2"):
                        secondary("Ablación guiada", lambda: ablate(True), icon="my_location").props(
                            "dense"
                        ).classes("w-full")
                        secondary("Ablación disco", lambda: ablate(False), icon="radio_button_checked").props(
                            "dense"
                        ).classes("w-full")
                    secondary("Análisis de ctDNA", lambda: ctdna(), icon="bloodtype").props("dense").classes(
                        "w-full"
                    )

                with card("Células en la vista", icon="bar_chart"):
                    counts_box = ui.column().classes("w-full gap-1")

                with card("Qué mirar", icon="lightbulb"):
                    ui.markdown(
                        "La **fracción accesible** es la parte del tumor con una célula vacía al lado: por "
                        "ahí, y solo por ahí, puede entrar una terapia. En los estudios cae del 82 % al 33 % al "
                        "crecer. Los **fármacos solos no llegan al núcleo**; la **ablación guiada** va a cada "
                        "foco visible."
                    ).classes("text-sm")

    # ---------------------------------------------------------------- lógica
    def refresh(with_access: bool = True) -> None:
        sim = s.sim
        cells, agents, counts = pack(sim, cutaway=s.cutaway, tumor_only=s.tumor_only)
        cloud.set_cells(cells)
        agents_cloud.set_cells(agents)
        r = sim.summary()
        k_day.set_text(fmt(round(sim.day, 1)))
        k_tumor.set_text(fmt(r["tumor_cells"]))
        k_danger.set_text(fmt(r["dangerous"]))
        k_queen.set_text(fmt(r["queens"]))
        k_ai.set_text(fmt(sum(a.kind == "AI" for a in sim.agents)))
        if with_access:
            k_access.set_text(f"{sim.accessible_fraction():.0%}")
        title.set_text(f"Día {int(sim.day)}")
        zl, side, _ = sim.grid_shape
        status.set_text(
            f"{SCENARIO_LABELS[sim.scenario]} · {side}×{side}×{zl} · {fmt(len(sim.cells))} células"
        )
        counts_box.clear()
        with counts_box:
            for code in LEGEND_3D:
                n = counts.get(code, 0) if code != AGENT_3D else sum(a.kind == "AI" for a in sim.agents)
                with ui.row().classes("w-full items-center justify-between no-wrap"):
                    legend_item(PALETTE_3D[code] if code != NORMAL else "#f4b6c8", LABEL_3D[code])
                    ui.label(fmt(n)).classes("text-sm font-medium").style(
                        "font-variant-numeric: tabular-nums"
                    )

    def camera(view: str) -> None:
        zl, side, _ = s.sim.grid_shape
        d = max(zl, side) * 1.05
        if view == "side":
            scene.move_camera(x=0, y=d, z=0.2 * d, look_at_x=0, look_at_y=0, look_at_z=0)
        elif view == "front":
            scene.move_camera(x=-d, y=0, z=0, look_at_x=0, look_at_y=0, look_at_z=0)
        else:
            scene.move_camera(x=-0.45 * d, y=0.8 * d, z=0.5 * d, look_at_x=0, look_at_y=0, look_at_z=0)

    async def advance(days: int) -> None:
        if s.lock.locked():
            return
        steps = days * STEPS_PER_DAY
        s.cancel = False
        progress.visible = stop.visible = True
        async with s.lock:
            done, last = 0, time.time()
            while done < steps and not s.cancel:
                chunk = min(STEPS_PER_DAY, steps - done)
                await run.io_bound(lambda: [s.sim.step() for _ in range(chunk)])
                done += chunk
                progress.set_value(done / steps)
                if time.time() - last > 1.5:
                    refresh(with_access=False)
                    last = time.time()
        progress.visible = stop.visible = False
        refresh()

    async def new_simulation() -> None:
        async with s.lock:
            s.scenario, s.size = scenario.value, size.value
            s.seed = None if seed.value is None else int(seed.value)
            await run.io_bound(s.reset)
        refresh()
        camera("iso")

    async def release() -> None:
        names = therapy.value or []
        if not names:
            ui.notify("Elige al menos una terapia.", type="warning")
            return
        async with s.lock:
            s.sim.release_therapy("+".join(names), 300)
            s.sim.p.mutation_rate = CALIBRATIONS["realistic"]
        ui.notify(
            f"Liberada: {therapy_label('+'.join(names))} (300 células). Tasa de mutación realista desde ahora.",
            type="positive",
        )
        refresh()

    async def ablate(guided: bool) -> None:
        async with s.lock:
            r = await run.io_bound(s.sim.ablate, guided=guided, rim=1 if guided else 2)
        ui.notify(
            f"Ablación {'guiada' if guided else 'en disco'}: {fmt(r['destroyed'])} tumorales y "
            f"{fmt(r['healthy'])} sanas destruidas · residuo {fmt(r['residual'])}",
            multi_line=True,
            timeout=8000,
        )
        refresh()

    async def ctdna() -> None:
        async with s.lock:
            r = s.sim.ctdna(2.25e-3)
        ui.notify(
            f"ctDNA {'POSITIVO' if r['positive'] else 'negativo'} · {fmt(r['fragments'])} fragmentos "
            f"(esperados {r['expected']:.1f}). Calibrado a una sensibilidad del 75 % el día 60.",
            type="warning" if r["positive"] else "positive",
            multi_line=True,
            timeout=8000,
        )

    async def first_paint() -> None:
        refresh()
        camera("iso")

    for switch in view_switches:
        switch.on_value_change(lambda _: refresh())
    point.on_value_change(lambda e: cloud.set_size(e.value))
    ui.timer(0.3, first_paint, once=True)
