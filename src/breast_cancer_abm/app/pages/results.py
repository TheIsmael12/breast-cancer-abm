"""Resultados: el resumen de cada banco de pruebas, sus gráficas y las salidas brutas."""

from __future__ import annotations

from nicegui import ui

from breast_cancer_abm.app import content
from breast_cancer_abm.app.pages.home import _bars
from breast_cancer_abm.app.theme import ACCENT, NavLink, card, fmt, ghost, shell, sidepanel


def page(set: str = "") -> None:  # noqa: A002 - el parámetro de la URL se llama así
    sets = content.result_sets()
    with shell(
        "/results",
        "Resultados",
        "Cada carpeta de results/ guarda la salida bruta de un día de experimentos y un resumen. Las cifras "
        "comparan brazos dentro del modelo; no son predicciones clínicas.",
    ):
        if not sets:
            ui.label("Todavía no hay resultados.").classes("muted")
            return
        current = next((s for s in sets if s.date == set), sets[0])
        with ui.row().classes("w-full no-wrap items-start gap-6"):
            sidepanel(
                [
                    (
                        "Fecha",
                        [
                            NavLink(f"/results?set={s.date}", s.date, s is current, icon="event_note")
                            for s in sets
                        ],
                    )
                ]
            )
            with ui.column().classes("col gap-6").style("min-width: 0"):
                _summary(current)
                _outputs(current)


def _summary(s: content.ResultSet) -> None:
    if not s.summary:
        return
    rows = s.summary
    names = [r["arm"] for r in rows][::-1]
    control = [round(100 * int(r["controlled"]) / int(r["n"])) for r in rows][::-1]
    healthy = [int(float(r["healthy_destroyed"])) for r in rows][::-1]
    height = f"height: {max(220, 26 * len(rows) + 50)}px"
    with card(f"Resumen del {s.date}", "Una fila por brazo de tratamiento.", icon="insights"):
        with ui.row().classes("w-full gap-6 no-wrap"):
            ui.echart(_bars(names, control, "% de pacientes controladas", ACCENT, "{c} %")).classes(
                "col"
            ).style(height)
            ui.echart(_bars(names, healthy, "Células sanas destruidas (mediana)", "#4b5563", "{c}")).classes(
                "col"
            ).style(height)
        columns = [
            {
                "name": k,
                "label": v,
                "field": k,
                "align": "left" if k in ("arm", "file") else "right",
                "sortable": k != "file",
            }
            for k, v in content.SUMMARY_HEADERS.items()
            if k in rows[0]
        ]
        table_rows = [
            {
                **r,
                "healthy_destroyed": fmt(float(r["healthy_destroyed"])),
                "rate_after_treatment": {"accelerated": "acelerada", "realistic": "realista"}.get(
                    r["rate_after_treatment"], r["rate_after_treatment"]
                ),
            }
            for r in rows
        ]
        ui.table(columns=columns, rows=table_rows, row_key="arm", pagination=0).props(
            "flat bordered dense wrap-cells"
        ).classes("w-full")


def _outputs(s: content.ResultSet) -> None:
    if not s.outputs:
        return
    with card("Salidas brutas", "Lo que imprimió cada estudio, tal cual.", icon="terminal", flush=True):
        with (
            ui.tabs()
            .props("dense align=left no-caps active-color=primary indicator-color=primary")
            .classes("q-px-md") as tabs
        ):
            for p in s.outputs:
                ui.tab(p.name, label=p.stem)
        with ui.tab_panels(tabs, value=s.outputs[0].name).classes("w-full"):
            for p in s.outputs:
                with ui.tab_panel(p.name).classes("q-pa-none"):
                    text = p.read_text(encoding="utf-8")
                    with ui.row().classes("w-full justify-end q-px-md q-pt-sm"):
                        ghost(
                            "Copiar",
                            lambda _, t=text: (ui.clipboard.write(t), ui.notify("Copiado")),
                            icon="content_copy",
                        )
                    with ui.scroll_area().classes("w-full").style("height: 420px"):
                        ui.label(text).classes("raw-output")
