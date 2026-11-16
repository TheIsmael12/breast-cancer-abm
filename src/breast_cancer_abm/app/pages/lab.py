"""Laboratorio: diseñar un experimento paso a paso, ejecutarlo y comparar los brazos.

Pensado para dos públicos a la vez: quien no sabe programar (plantillas, nombres en
español, valores por defecto razonables) y quien investiga (los 80 parámetros del modelo,
cualquier combinación de intervenciones, exportación de los datos brutos).
"""

from __future__ import annotations

import csv
import io
import json
import os
import shutil
import time
from datetime import datetime

from nicegui import run, ui

from breast_cancer_abm.app.theme import ACCENT, card, choice, fmt, ghost, primary, secondary, shell
from breast_cancer_abm.catalog import PARAM_GROUPS, SCENARIOS, THERAPIES, VESSEL_MODES
from breast_cancer_abm.lab.experiment import (
    INTERVENTIONS,
    Arm,
    Experiment,
    Intervention,
    default_processes,
    median_curves,
    run_experiment,
    summarize,
)
from breast_cancer_abm.lab.templates import EXPERIMENT_TEMPLATES, arm_templates
from breast_cancer_abm.model.biology import Parameters
from breast_cancer_abm.paths import LAB_DIR

TEMPLATE_ICONS = {
    "no_surgery": "content_cut",
    "immunotherapy": "shield",
    "cell_therapies": "vaccines",
    "3d": "view_in_ar",
    "blank": "add_circle_outline",
}
ARM_COLORS = ["#6b7280", ACCENT, "#2563eb", "#0d9488", "#d97706", "#7c3aed", "#dc2626", "#0891b2", "#65a30d"]
THERAPY_OPTIONS = {k: f"{v[0]} · {v[1]}" for k, v in THERAPIES.items()}
DEFAULTS = Parameters()


class Lab:
    """Estado de una pestaña del laboratorio."""

    def __init__(self) -> None:
        self.exp = EXPERIMENT_TEMPLATES["no_surgery"][1]()
        self.results: list[dict] = []
        self.running = False
        self.cancel = False
        self.done = 0
        self.total = 0
        self.started = 0.0
        self.processes = default_processes()
        self.saved_as: str | None = None


def page() -> None:
    lab = Lab()

    with shell(
        "/lab",
        "Laboratorio",
        "Diseña un experimento sin escribir código: elige las pacientes, los tratamientos a comparar y, si "
        "quieres, la biología. Las mismas pacientes pasan por todos los brazos, así que la única diferencia "
        "entre brazos es el tratamiento.",
    ) as pg:
        with pg.actions:
            secondary("Experimentos guardados", lambda: open_saved(), icon="folder_open")
            secondary("Importar", lambda: import_dialog(), icon="upload_file")

        # ------------------------------------------------------------ plantillas
        with card(
            "Empieza por una plantilla",
            "Carga un experimento completo y cámbialo a tu gusto.",
            icon="auto_awesome",
        ):
            with ui.element("div").classes("choice-grid"):
                for key, (title, make) in EXPERIMENT_TEMPLATES.items():
                    choice(
                        title,
                        make().question or "Un brazo de control y lo que tú añadas.",
                        False,
                        lambda _, k=key: load_template(k),
                        icon=TEMPLATE_ICONS.get(key, "science"),
                    )

        with card("Diseño del experimento", icon="tune"):
            with ui.stepper().props("flat animated header-nav").classes("w-full lab-stepper") as stepper:
                with ui.step("Pacientes", icon="person").props("done-color=primary active-color=primary"):
                    patients_box = ui.column().classes("w-full gap-5")
                    with ui.stepper_navigation():
                        primary("Siguiente: tratamientos", stepper.next, icon="arrow_forward")
                with ui.step("Tratamientos", icon="medication").props(
                    "done-color=primary active-color=primary"
                ):
                    arms_box = ui.column().classes("w-full gap-4")
                    with ui.stepper_navigation():
                        primary("Siguiente: biología", stepper.next, icon="arrow_forward")
                        ghost("Atrás", stepper.previous)
                with ui.step("Biología (opcional)", icon="biotech").props(
                    "done-color=primary active-color=primary"
                ):
                    biology_box = ui.column().classes("w-full gap-3")
                    with ui.stepper_navigation():
                        primary("Siguiente: ejecutar", stepper.next, icon="arrow_forward")
                        ghost("Atrás", stepper.previous)
                with ui.step("Ejecutar", icon="play_circle").props("done-color=primary active-color=primary"):
                    run_box = ui.column().classes("w-full gap-4")
                    with ui.stepper_navigation():
                        ghost("Atrás", stepper.previous)

        results_box = ui.column().classes("w-full gap-6")

    # =================================================================== paso 1
    @ui.refreshable
    def patients() -> None:
        e = lab.exp
        with ui.row().classes("w-full gap-4 items-start"):
            ui.input("Nombre del experimento").bind_value(e, "name").props("outlined dense").classes("col")
        ui.textarea("Pregunta que quieres responder (opcional)").bind_value(e, "question").props(
            "outlined dense autogrow"
        ).classes("w-full")
        ui.label("Modelo").classes("section-title")
        with ui.element("div").classes("choice-grid"):
            for model, title, text in (
                (
                    "2d",
                    "2D · corte de un conducto",
                    "Rápido. Ideal para comparar muchas estrategias con muchas pacientes.",
                ),
                (
                    "3d",
                    "3D · el conducto como cilindro",
                    "Más realista para el acceso al tumor, mucho más lento. "
                    "Úsalo para confirmar lo que veas en 2D.",
                ),
            ):
                choice(title, text, e.model == model, lambda m=model: (setattr(e, "model", m), refresh_all()))
        if e.model == "3d":
            ui.select(
                {"fast": "Rápido (34 × 34 × 16)", "standard": "Estándar (54 × 54 × 40, el de los estudios)"},
                label="Tamaño 3D",
            ).bind_value(e, "size3d").props("outlined dense").style("width: 360px")
        ui.label("Tipo de cáncer").classes("section-title")
        with ui.element("div").classes("choice-grid"):
            for key, (title, text) in SCENARIOS.items():
                choice(
                    title,
                    text,
                    e.scenario == key,
                    lambda k=key: (setattr(e, "scenario", k), patients.refresh()),
                )
        ui.label("Tamaño del estudio").classes("section-title")
        with ui.grid(columns=3).classes("w-full gap-4"):
            _slider_field(
                "Pacientes por brazo",
                "Más pacientes, más fiable (la guía pide ≥ 12 para publicar).",
                e,
                "patients",
                1,
                48,
                1,
            )
            _slider_field(
                "Duración (días)",
                "Un paso son 6 horas. El tumor suele aparecer hacia el día 60-120.",
                e,
                "days",
                30,
                800,
                10,
            )
            _slider_field(
                "«Controlada» si quedan menos de",
                "células peligrosas al final (carcinoma + reinas).",
                e,
                "control_cells",
                0,
                2000,
                10,
            )
        with ui.expansion("Semillas (reproducibilidad)", icon="tune").classes("w-full"):
            ui.number("Primera semilla", min=0, step=1, format="%d").bind_value(
                e, "first_seed", forward=lambda v: int(v or 0)
            ).props("outlined dense").style("width: 200px")
            ui.label(
                "Cada paciente es una semilla del azar. Con las mismas semillas, el experimento da exactamente "
                "los mismos números."
            ).classes("text-xs muted")

    # =================================================================== paso 2
    @ui.refreshable
    def arms() -> None:
        e = lab.exp
        with ui.row().classes("w-full items-center gap-3"):
            day = (
                ui.number("Día del tratamiento (plantillas)", value=120, min=0, step=5, format="%d")
                .props("outlined dense")
                .style("width: 230px")
            )
            with ui.button("Añadir brazo", icon="add").props("unelevated"):
                with ui.menu():
                    for title in arm_templates(0):
                        ui.menu_item(title, on_click=lambda _, t=title: add_arm(t, day.value or 0))
            ui.label(f"{len(e.arms)} brazos · {len(e.arms) * e.patients} simulaciones").classes(
                "text-sm muted"
            )
        for i, arm in enumerate(e.arms):
            color = ARM_COLORS[i % len(ARM_COLORS)]
            with ui.element("div").classes("panel w-full column gap-2"):
                with ui.row().classes("w-full items-center no-wrap gap-3"):
                    ui.element("span").classes("swatch").style(f"background:{color}; width:14px; height:14px")
                    ui.input("Nombre del brazo").bind_value(arm, "name").props("outlined dense").classes(
                        "col"
                    )
                    ui.button(icon="add", on_click=lambda _, a=arm: edit_intervention(a, None)).props(
                        "flat dense"
                    ).tooltip("Añadir intervención")
                    ui.button(icon="content_copy", on_click=lambda _, a=arm: duplicate_arm(a)).props(
                        "flat dense"
                    ).tooltip("Duplicar brazo")
                    ui.button(icon="delete_outline", on_click=lambda _, k=i: remove_arm(k)).props(
                        "flat dense color=negative"
                    ).tooltip("Quitar brazo")
                if not arm.interventions:
                    ui.label("Sin intervenciones: es el grupo de control (la evolución natural).").classes(
                        "text-sm muted"
                    )
                for iv in sorted(arm.interventions, key=lambda x: x.day):
                    with ui.row().classes("w-full items-center no-wrap gap-2 lab-iv"):
                        ui.icon("schedule", size="18px").classes("muted")
                        ui.label(iv.describe()).classes("text-sm col")
                        ui.button(icon="edit", on_click=lambda _, a=arm, v=iv: edit_intervention(a, v)).props(
                            "flat dense round size=sm"
                        )
                        ui.button(icon="close", on_click=lambda _, a=arm, v=iv: remove_iv(a, v)).props(
                            "flat dense round size=sm"
                        )

    # =================================================================== paso 3
    @ui.refreshable
    def biology() -> None:
        e = lab.exp
        changed = {k: v for k, v in e.params.items() if getattr(DEFAULTS, k, None) != v}
        with ui.row().classes("w-full items-center justify-between"):
            ui.label(
                f"{len(changed)} parámetros cambiados respecto al modelo por defecto. Si no sabes qué tocar, "
                "déjalo así: los valores por defecto son los de los estudios publicados."
            ).classes("text-sm muted col")
            ui.button(
                "Restablecer todo", icon="restart_alt", on_click=lambda: (e.params.clear(), biology.refresh())
            ).props("flat dense no-caps color=grey-8")
        for title, text, params in PARAM_GROUPS:
            n = sum(p.key in changed for p in params)
            with ui.expansion(f"{title}" + (f"  ·  {n} cambiados" if n else ""), icon="science").classes(
                "w-full panel q-pa-none"
            ):
                ui.label(text).classes("text-sm muted q-mb-sm")
                with (
                    ui.grid()
                    .classes("w-full gap-x-6 gap-y-3")
                    .style("grid-template-columns: repeat(auto-fill, minmax(340px, 1fr))")
                ):
                    for p in params:
                        _param_field(e, p)

    def _param_field(e: Experiment, p) -> None:
        default = getattr(DEFAULTS, p.key)
        value = e.params.get(p.key, default)
        modified = value != default

        def set_value(v) -> None:
            if v is None:
                return
            if isinstance(default, bool):
                v = bool(v)
            elif isinstance(default, int):
                v = int(v)
            else:
                v = float(v)
            if v == default:
                e.params.pop(p.key, None)
            else:
                e.params[p.key] = v

        with ui.column().classes("gap-0"):
            if isinstance(default, bool):
                ui.switch(p.label, value=value, on_change=lambda ev: set_value(ev.value))
            else:
                with ui.row().classes("w-full items-center no-wrap gap-2"):
                    ui.number(
                        p.label,
                        value=value,
                        min=p.lo,
                        max=p.hi,
                        step=p.step,
                        on_change=lambda ev: set_value(ev.value),
                    ).props("outlined dense" + (" color=primary" if modified else "")).classes("col")
            hint = f"Por defecto: {default}" + (" · CAMBIADO" if modified else "")
            ui.label(f"{p.help} {hint}".strip()).classes("text-xs muted").style(
                f"color:{ACCENT}" if modified else ""
            )

    # =================================================================== paso 4
    @ui.refreshable
    def run_panel() -> None:
        e = lab.exp
        cpus = os.cpu_count() or 2
        with ui.element("div").classes("panel w-full column gap-2"):
            ui.label("Resumen del diseño").classes("section-title")
            ui.label(
                f"{e.name} · modelo {e.model.upper()} · {SCENARIOS[e.scenario][0]} · {e.patients} pacientes "
                f"por brazo · {int(e.days)} días · {len(e.arms)} brazos"
            ).classes("text-sm")
            for i, arm in enumerate(e.arms):
                ui.label(
                    f"{arm.name}: " + ("; ".join(iv.describe() for iv in arm.interventions) or "sin tratar")
                ).classes("text-sm muted")
        with ui.row().classes("w-full items-center gap-4"):
            ui.label("Procesos en paralelo").classes("text-sm")
            ui.slider(min=1, max=max(1, cpus - 2), step=1).bind_value(lab, "processes").props("label").style(
                "width: 260px"
            )
            eta = ui.label().classes("text-sm muted")
            eta.bind_text_from(
                lab,
                "processes",
                backward=lambda p: (
                    f"≈ {_duration(e.estimate_seconds(p))} con {p} de {cpus} núcleos. Deja núcleos libres para "
                    "seguir usando el ordenador."
                ),
            )
        with ui.row().classes("gap-2"):
            ui.button("Ejecutar experimento", icon="play_arrow", on_click=lambda: start()).props(
                "unelevated size=lg"
            ).bind_enabled_from(lab, "running", backward=lambda r: not r)
            ui.button("Detener", icon="stop", on_click=lambda: setattr(lab, "cancel", True)).props(
                "outline color=negative"
            ).bind_visibility_from(lab, "running")
        bar = ui.linear_progress(value=0, show_value=False, size="10px").props("rounded track-color=grey-3")
        status = ui.label().classes("text-sm muted")

        def update() -> None:
            bar.visible = lab.running or lab.done > 0
            if lab.total:
                bar.set_value(lab.done / lab.total)
                left = (time.time() - lab.started) / max(lab.done, 1) * (lab.total - lab.done)
                status.set_text(
                    f"{lab.done} de {lab.total} simulaciones"
                    + (f" · quedan ≈ {_duration(left)}" if lab.running and lab.done else "")
                )

        ui.timer(0.5, update)

    # =================================================================== acciones
    def refresh_all() -> None:
        patients.refresh()
        arms.refresh()
        biology.refresh()
        run_panel.refresh()

    def load_template(key: str) -> None:
        lab.exp = EXPERIMENT_TEMPLATES[key][1]()
        lab.results = []
        refresh_all()
        results.refresh()
        stepper.set_value("Pacientes")
        ui.notify(f"Plantilla cargada: {EXPERIMENT_TEMPLATES[key][0]}")

    def add_arm(title: str, day: float) -> None:
        lab.exp.arms.append(arm_templates(day)[title]())
        arms.refresh()
        run_panel.refresh()

    def duplicate_arm(arm: Arm) -> None:
        copy = Experiment.from_dict(
            {
                "arms": [
                    {"name": arm.name + " (copia)", "interventions": [vars(iv) for iv in arm.interventions]}
                ]
            }
        ).arms[0]
        lab.exp.arms.append(copy)
        arms.refresh()

    def remove_arm(index: int) -> None:
        lab.exp.arms.pop(index)
        arms.refresh()
        run_panel.refresh()

    def remove_iv(arm: Arm, iv: Intervention) -> None:
        arm.interventions.remove(iv)
        arms.refresh()

    def edit_intervention(arm: Arm, iv: Intervention | None) -> None:
        model = lab.exp.model
        kinds = {k: v[0] for k, v in INTERVENTIONS.items() if (v[1] if model == "2d" else v[2])}
        draft = (
            Intervention(**vars(iv))
            if iv
            else Intervention("therapy", 120, therapy="trojan+niche_guardian+memory")
        )
        if draft.kind not in kinds:
            draft.kind = "therapy"
        with ui.dialog() as dialog, ui.card().classes("q-pa-lg").style("width: 620px; max-width: 94vw"):
            ui.label("Intervención" if iv else "Nueva intervención").classes("text-lg font-semibold")
            with ui.row().classes("w-full gap-3 no-wrap"):
                kind_select = (
                    ui.select(kinds, label="Tipo")
                    .bind_value(draft, "kind")
                    .props("outlined dense")
                    .classes("col")
                )
                ui.number("Día", min=0, step=1, format="%g").bind_value(draft, "day").props(
                    "outlined dense"
                ).style("width: 120px")
            box = ui.column().classes("w-full gap-3")

            @ui.refreshable
            def fields() -> None:
                k = draft.kind
                if k in ("therapy", "ctdna_monitor"):
                    selected = [t for t in draft.therapy.split("+") if t]
                    label = (
                        "Terapia de rescate (opcional)" if k == "ctdna_monitor" else "Terapia o combinación"
                    )
                    sel = (
                        ui.select(
                            THERAPY_OPTIONS, value=selected, multiple=True, label=label, with_input=True
                        )
                        .props("outlined dense use-chips options-dense")
                        .classes("w-full")
                    )
                    sel.on_value_change(lambda ev: setattr(draft, "therapy", "+".join(ev.value or [])))
                    ui.number("Dosis (células)", min=1, step=10, format="%d").bind_value(
                        draft, "dose", forward=lambda v: int(v or 1)
                    ).props("outlined dense").style("width: 200px")
                if k == "ablate":
                    ui.switch("Guiada por imagen (solo los focos visibles); si no, en disco").bind_value(
                        draft, "guided"
                    )
                    ui.number(
                        "La imagen ve focos de al menos (células)", min=1, step=1, format="%d"
                    ).bind_value(draft, "visible_threshold", forward=lambda v: int(v or 1)).props(
                        "outlined dense"
                    )
                if k == "ctdna_monitor":
                    with ui.row().classes("w-full gap-3 no-wrap"):
                        ui.number("Cada (días)", min=1, step=1, format="%g").bind_value(draft, "every").props(
                            "outlined dense"
                        ).classes("col")
                        ui.number("Hasta el día (0 = final)", min=0, step=10, format="%g").bind_value(
                            draft, "until"
                        ).props("outlined dense").classes("col")
                        ui.number(
                            "Sensibilidad (fragmentos/célula)", min=0.0001, step=0.01, format="%g"
                        ).bind_value(draft, "fragments").props("outlined dense").classes("col")
                    ui.select(
                        {
                            "guided": "Ablación guiada",
                            "disc": "Ablación en disco",
                            "therapy": "Solo la terapia",
                        },
                        label="Si da positivo",
                    ).bind_value(draft, "rescue").props("outlined dense")
                    ui.label(
                        "El ctDNA es ADN tumoral en sangre: un análisis por extracción. El rescate se hace una "
                        "sola vez, en el primer positivo."
                    ).classes("text-xs muted")
                if k == "block_pd1":
                    ui.number("Duración del anticuerpo (días)", min=1, step=5, format="%g").bind_value(
                        draft, "duration"
                    ).props("outlined dense")
                if k == "vessels":
                    ui.select(VESSEL_MODES, label="Tratamiento de los vasos").bind_value(draft, "mode").props(
                        "outlined dense"
                    )

            with box:
                fields()
            kind_select.on_value_change(lambda _: fields.refresh())

            def save() -> None:
                if draft.kind == "therapy" and not draft.therapy:
                    ui.notify("Elige al menos una terapia.", type="warning")
                    return
                if iv:
                    arm.interventions[arm.interventions.index(iv)] = draft
                else:
                    arm.interventions.append(draft)
                dialog.close()
                arms.refresh()
                run_panel.refresh()

            with ui.row().classes("w-full justify-end gap-2 q-mt-md"):
                ui.button("Cancelar", on_click=dialog.close).props("flat no-caps color=grey-8")
                ui.button("Guardar", on_click=save).props("unelevated")
        dialog.open()

    async def start() -> None:
        e = lab.exp
        if not e.arms:
            ui.notify("Añade al menos un brazo de tratamiento.", type="warning")
            return
        lab.running, lab.cancel, lab.results = True, False, []
        lab.done, lab.total, lab.started = 0, e.patients * len(e.arms), time.time()
        lab.saved_as = None
        results.refresh()

        def on_result(r: dict, done: int, total: int) -> None:
            lab.results.append(r)
            lab.done = done

        try:
            await run.io_bound(run_experiment, e, int(lab.processes), on_result, lambda: lab.cancel)
        finally:
            lab.running = False
        if lab.results:
            lab.saved_as = str(_write_run(e, lab.results))  # se guarda siempre, aunque se cierre la ventana
        ui.notify(
            ("Experimento terminado" if not lab.cancel else "Experimento detenido")
            + " y guardado en «Experimentos guardados».",
            type="positive" if not lab.cancel else "warning",
        )
        results.refresh()

    def open_saved() -> None:
        runs = (
            sorted((p for p in LAB_DIR.glob("*") if (p / "experiment.json").exists()), reverse=True)
            if LAB_DIR.exists()
            else []
        )
        with ui.dialog() as dialog, ui.card().classes("q-pa-lg").style("width: 640px; max-width: 94vw"):
            ui.label("Experimentos guardados").classes("text-lg font-semibold")
            ui.label(str(LAB_DIR)).classes("text-xs muted mono")
            if not runs:
                ui.label("Todavía no has guardado ningún experimento.").classes("text-sm muted")
            for p in runs[:50]:
                exp = json.loads((p / "experiment.json").read_text(encoding="utf-8"))
                with ui.row().classes("w-full items-center no-wrap gap-2"):
                    with ui.column().classes("col gap-0"):
                        ui.label(exp.get("name", p.name)).classes("text-sm font-semibold")
                        ui.label(p.name[:17].replace("_", " ")).classes("text-xs muted")
                    ui.button("Abrir", on_click=lambda _, f=p: (load_run(f), dialog.close())).props(
                        "flat dense"
                    )
                    ui.button(
                        icon="delete_outline",
                        on_click=lambda _, f=p: (shutil.rmtree(f, True), dialog.close(), open_saved()),
                    ).props("flat dense color=negative")
            with ui.row().classes("w-full justify-end"):
                ui.button("Cerrar", on_click=dialog.close).props("flat no-caps color=grey-8")
        dialog.open()

    def load_run(folder) -> None:
        lab.exp = Experiment.from_dict(json.loads((folder / "experiment.json").read_text(encoding="utf-8")))
        rpath = folder / "results.json"
        lab.results = json.loads(rpath.read_text(encoding="utf-8")) if rpath.exists() else []
        lab.saved_as = str(folder)
        refresh_all()
        results.refresh()

    def import_dialog() -> None:
        with ui.dialog() as dialog, ui.card().classes("q-pa-lg").style("width: 520px; max-width: 94vw"):
            ui.label("Importar un experimento").classes("text-lg font-semibold")
            ui.label("Un fichero .json descargado desde este laboratorio (con o sin resultados).").classes(
                "text-sm muted"
            )

            async def handle(e) -> None:
                await import_json(e)
                dialog.close()

            ui.upload(auto_upload=True, max_files=1, on_upload=handle).props(
                "accept=.json flat bordered"
            ).classes("w-full")
            with ui.row().classes("w-full justify-end"):
                ui.button("Cerrar", on_click=dialog.close).props("flat no-caps color=grey-8")
        dialog.open()

    async def import_json(e) -> None:
        try:
            data = await e.file.json()
            lab.exp = Experiment.from_dict(data.get("experiment", data))
            lab.results = data.get("results", [])
            refresh_all()
            results.refresh()
            ui.notify("Experimento importado.", type="positive")
        except Exception as exc:  # noqa: BLE001 - el usuario puede subir cualquier cosa
            ui.notify(f"No es un experimento válido: {exc}", type="negative")

    # =================================================================== resultados
    @ui.refreshable
    def results() -> None:
        e, res = lab.exp, lab.results
        if not res or lab.running:
            return
        rows = summarize(e, res)
        ui.label("Resultados").classes("page-title q-mt-md")
        if e.question:
            ui.label(e.question).classes("page-sub")
        with ui.row().classes("gap-2"):
            if lab.saved_as:
                ui.label(f"Guardado en {lab.saved_as}").classes("text-xs muted self-center")
            ui.button(
                "Descargar CSV",
                icon="table_view",
                on_click=lambda: ui.download.content(_summary_csv(e, res), f"{e.name}.csv"),
            ).props("outline no-caps color=grey-9").classes("btn-secondary")
            ui.button(
                "Descargar JSON completo",
                icon="data_object",
                on_click=lambda: ui.download.content(
                    json.dumps({"experiment": e.to_dict(), "results": res}, ensure_ascii=False),
                    f"{e.name}.json",
                ),
            ).props("outline no-caps color=grey-9").classes("btn-secondary")
        names = [r["arm"] for r in rows]
        with ui.row().classes("w-full gap-4 no-wrap items-stretch"):
            ui.echart(
                _grouped_bars(
                    names,
                    [
                        ("Controladas", [round(100 * r["controlled"] / r["n"]) for r in rows], ACCENT),
                        (
                            "Erradicadas (0 peligrosas)",
                            [round(100 * r["eradicated"] / r["n"]) for r in rows],
                            "#111827",
                        ),
                        (
                            "Sin tocar la mama",
                            [round(100 * r["untouched"] / r["n"]) for r in rows],
                            "#9ca3af",
                        ),
                    ],
                    "% de pacientes",
                    "{c} %",
                )
            ).classes("col panel").style("height: 300px")
            ui.echart(
                _grouped_bars(
                    names,
                    [
                        ("Peligrosas al final (mediana)", [r["median_dangerous"] for r in rows], "#b91c1c"),
                        (
                            "Sanas destruidas (mediana)",
                            [r["median_healthy_destroyed"] for r in rows],
                            "#4b5563",
                        ),
                    ],
                    "células",
                    "{c}",
                )
            ).classes("col panel").style("height: 300px")
        with ui.element("div").classes("panel w-full"):
            ui.label(
                "Células peligrosas a lo largo del tiempo (mediana de las pacientes, escala logarítmica)"
            ).classes("section-title")
            ui.echart(_curves(e, res)).classes("w-full").style("height: 340px")
        columns = [
            {"name": "arm", "label": "Brazo", "field": "arm", "align": "left"},
            {"name": "controlled", "label": f"Controladas (<{e.control_cells})", "field": "controlled"},
            {"name": "eradicated", "label": "Erradicadas", "field": "eradicated"},
            {"name": "median_dangerous", "label": "Peligrosas (mediana)", "field": "median_dangerous"},
            {
                "name": "median_healthy_destroyed",
                "label": "Sanas destruidas (mediana)",
                "field": "median_healthy_destroyed",
            },
            {"name": "untouched", "label": "Sin tocar la mama", "field": "untouched"},
            {"name": "per_patient", "label": "Por paciente", "field": "per_patient", "align": "left"},
        ]
        table_rows = [
            {
                **r,
                "controlled": f"{r['controlled']}/{r['n']}",
                "eradicated": f"{r['eradicated']}/{r['n']}",
                "untouched": f"{r['untouched']}/{r['n']}",
                "median_dangerous": fmt(r["median_dangerous"]),
                "median_healthy_destroyed": fmt(r["median_healthy_destroyed"]),
                "per_patient": " · ".join(fmt(x) for x in r["per_patient"]),
            }
            for r in rows
        ]
        ui.table(columns=columns, rows=table_rows, row_key="arm", pagination=0).props(
            "flat bordered wrap-cells"
        ).classes("w-full")
        with ui.expansion("Diario de cada paciente (qué se hizo y cuándo)", icon="list_alt").classes(
            "w-full panel"
        ):
            for r in sorted(res, key=lambda x: (x["arm"], x["seed"])):
                arm_name = e.arms[r["arm"]].name if r["arm"] < len(e.arms) else r["arm"]
                ui.label(
                    f"{arm_name} · paciente {r['seed']} · {fmt(r['dangerous'])} peligrosas al final · "
                    f"{fmt(r['healthy_destroyed'])} sanas destruidas"
                ).classes("text-sm font-semibold q-mt-sm")
                for day, text in r["log"] or [(None, "Sin intervenciones")]:
                    ui.label(f"   {'' if day is None else f'día {day:g} · '}{text}").classes("text-xs muted")
        ui.label(
            "Recuerda: las cifras comparan brazos dentro del modelo (tasas aceleradas, tiempo comprimido). "
            "No son predicciones clínicas."
        ).classes("text-xs muted")

    # =================================================================== montaje
    with patients_box:
        patients()
    with arms_box:
        arms()
    with biology_box:
        biology()
    with run_box:
        run_panel()
    with results_box:
        results()
    stepper.on_value_change(lambda _: run_panel.refresh())


# ======================================================================= ayudas
def _write_run(e: Experiment, res: list[dict]):
    stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    slug = "".join(c if c.isalnum() else "-" for c in e.name.lower()).strip("-")[:40] or "experimento"
    folder = LAB_DIR / f"{stamp}_{slug}"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "experiment.json").write_text(e.to_json(), encoding="utf-8")
    (folder / "results.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
    (folder / "summary.csv").write_text(_summary_csv(e, res), encoding="utf-8")
    return folder


def _slider_field(label: str, help_: str, obj, attr: str, lo: float, hi: float, step: float) -> None:
    with ui.column().classes("gap-0"):
        with ui.row().classes("w-full items-baseline justify-between"):
            ui.label(label).classes("text-sm font-medium")
            ui.label().bind_text_from(obj, attr, backward=lambda v: fmt(v)).classes(
                "text-base font-semibold"
            ).style(f"color:{ACCENT}")
        ui.slider(min=lo, max=hi, step=step).bind_value(obj, attr, forward=lambda v: int(v))
        ui.label(help_).classes("text-xs muted")


def _duration(seconds: float) -> str:
    if seconds < 90:
        return f"{max(1, round(seconds))} s"
    if seconds < 5400:
        return f"{round(seconds / 60)} min"
    return f"{seconds / 3600:.1f} h"


def _summary_csv(e: Experiment, res: list[dict]) -> str:
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\n")
    w.writerow(
        ["arm", "seed", "dangerous", "tumor_cells", "queens", "healthy_destroyed", "rescue_day", "untouched"]
    )
    for r in sorted(res, key=lambda x: (x["arm"], x["seed"])):
        w.writerow(
            [
                e.arms[r["arm"]].name,
                r["seed"],
                r["dangerous"],
                r["tumor_cells"],
                r["queens"],
                r["healthy_destroyed"],
                r["rescue_day"] if r["rescue_day"] is not None else "",
                r["untouched"],
            ]
        )
    return out.getvalue()


def _grouped_bars(names: list[str], series: list[tuple[str, list, str]], unit: str, label: str) -> dict:
    return {
        "animation": False,
        "legend": {
            "top": 0,
            "icon": "roundRect",
            "itemWidth": 10,
            "itemHeight": 10,
            "textStyle": {"fontSize": 11},
        },
        "grid": {"left": 50, "right": 10, "top": 36, "bottom": 60},
        "tooltip": {"trigger": "axis"},
        "xAxis": {
            "type": "category",
            "data": names,
            "axisLabel": {"fontSize": 11, "interval": 0, "width": 120, "overflow": "break"},
        },
        "yAxis": {
            "type": "value",
            "name": unit,
            "nameTextStyle": {"fontSize": 10},
            "splitLine": {"lineStyle": {"color": "#f3f4f6"}},
        },
        "series": [
            {
                "name": n,
                "type": "bar",
                "data": d,
                "color": c,
                "barMaxWidth": 28,
                "label": {"show": True, "position": "top", "formatter": label, "fontSize": 10},
                "itemStyle": {"borderRadius": [3, 3, 0, 0]},
            }
            for n, d, c in series
        ],
    }


def _curves(e: Experiment, res: list[dict]) -> dict:
    series = []
    for i, arm in enumerate(e.arms):
        days, med, _, _ = median_curves(res, i)
        if days:
            series.append(
                {
                    "name": arm.name,
                    "type": "line",
                    "showSymbol": False,
                    "color": ARM_COLORS[i % len(ARM_COLORS)],
                    "lineStyle": {"width": 2},
                    "data": [[d, v or None] for d, v in zip(days, med)],
                }
            )
    marks = sorted({iv.day for arm in e.arms for iv in arm.interventions})
    if series and marks:
        series[0]["markLine"] = {
            "symbol": "none",
            "silent": True,
            "label": {"formatter": "día {c}", "fontSize": 10},
            "lineStyle": {"color": "#d1d5db", "type": "dashed"},
            "data": [{"xAxis": m} for m in marks],
        }
    return {
        "animation": False,
        "legend": {
            "top": 0,
            "icon": "roundRect",
            "itemWidth": 10,
            "itemHeight": 10,
            "textStyle": {"fontSize": 11},
        },
        "grid": {"left": 50, "right": 20, "top": 40, "bottom": 30},
        "tooltip": {"trigger": "axis"},
        "xAxis": {"type": "value", "name": "días", "max": e.days, "splitLine": {"show": False}},
        "yAxis": {"type": "log", "min": 1, "splitLine": {"lineStyle": {"color": "#f3f4f6"}}},
        "series": series,
    }
