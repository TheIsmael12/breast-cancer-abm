"""Inicio: qué es el proyecto, por dónde empezar, el último avance y las figuras."""

from __future__ import annotations

from nicegui import ui

from breast_cancer_abm.app import content
from breast_cancer_abm.app.theme import ACCENT, REPO, card, choice, primary, secondary, shell

START = (
    (
        "/simulation",
        "blur_on",
        "Ver crecer un tumor",
        "Simulación 2D en vivo: reproduce, trata y haz clic en las células.",
    ),
    (
        "/simulation-3d",
        "view_in_ar",
        "Explorar en 3D",
        "El conducto como cilindro, con la cámara en la mano.",
    ),
    (
        "/lab",
        "biotech",
        "Diseñar un experimento",
        "Compara tratamientos con muchas pacientes, sin programar.",
    ),
    ("/results", "insights", "Ver los resultados", "Las salidas de cada estudio publicado, con gráficas."),
)


def page() -> None:
    with shell(
        "/",
        "Breast Cancer ABM",
        "Simulador por agentes de cáncer de mama y células terapéuticas programables. Cada idea se mide con "
        "varias pacientes, con el criterio escrito antes de ejecutar y con los errores propios documentados.",
    ) as pg:
        with pg.actions:
            secondary("GitHub", lambda: ui.navigate.to(REPO, new_tab=True), icon="code")
            primary("Abrir el laboratorio", lambda: ui.navigate.to("/lab"), icon="biotech")

        with ui.element("div").classes("choice-grid").style("grid-template-columns: repeat(4, 1fr)"):
            for path, icon, title, text in START:
                choice(title, text, False, lambda _, p=path: ui.navigate.to(p), icon=icon)

        with ui.row().classes("w-full no-wrap gap-6 items-stretch"):
            with ui.column().classes("col gap-0"):
                _key_results()
            latest = content.progress_docs()[:1]
            with ui.column().classes("gap-0").style("width: 420px; flex: none"):
                with card("Último avance", icon="campaign", classes="h-full"):
                    if latest:
                        d = latest[0]
                        ui.label(d.date).classes("hint")
                        ui.label(d.title).classes("text-base font-semibold")
                        ui.label(d.summary).classes("text-sm").style("line-height: 1.55")
                    ui.space()
                    secondary(
                        "Leer el avance",
                        lambda: ui.navigate.to(f"/docs?doc={latest[0].key}") if latest else None,
                        icon="article",
                    )

        with card("Figuras", "Haz clic en una para verla en grande.", icon="photo_library"):
            with (
                ui.grid()
                .classes("w-full gap-4")
                .style("grid-template-columns: repeat(auto-fill, minmax(280px, 1fr))")
            ):
                for url, caption in content.figures():
                    if "/app_" in url:
                        continue
                    with (
                        ui.column()
                        .classes("gap-2 cursor-pointer")
                        .on("click", lambda _, u=url, c=caption: _zoom(u, c))
                    ):
                        ui.image(url).classes("w-full").props("fit=cover ratio=1.6").style(
                            "border-radius: 6px; border: 1px solid var(--border)"
                        )
                        ui.label(caption).classes("hint")

        with card("Hecho con IA, y dicho claramente", icon="smart_toy"):
            ui.markdown(
                "Este proyecto lo desarrolla **Claude** (Anthropic): el código, la lectura de literatura, el diseño "
                "y la ejecución de los experimentos, el análisis y la documentación. Lo **dirige Ismael Ben**, que "
                "plantea las preguntas y decide en cada bifurcación. Es **educativo**: las tasas están aceleradas y "
                "el tejido es un corte de 1,65 mm, así que no sirve para decisiones médicas. Si ves un error "
                f"biológico, [dilo en GitHub]({REPO}/issues)."
            ).classes("text-sm")


def _key_results() -> None:
    """Las dos cifras que resumen el proyecto hoy: control y daño al tejido sano en 3D."""
    rows = [r for s in content.result_sets() for r in s.summary if r["arm"].startswith("3D:")]
    if not rows:
        return
    names = [r["arm"].removeprefix("3D: ") for r in rows][::-1]
    control = [round(100 * int(r["controlled"]) / int(r["n"])) for r in rows][::-1]
    healthy = [int(float(r["healthy_destroyed"])) for r in rows][::-1]
    with card(
        "Resultado clave · tratar sin operar en 3D",
        f"{rows[0]['n']} pacientes por brazo. Con ctDNA mensual y ablación guiada se controlan 8 de 8 destruyendo "
        "menos de la mitad de tejido sano que el equivalente a la cirugía.",
        icon="insights",
    ):
        with ui.row().classes("w-full no-wrap gap-4"):
            ui.echart(_bars(names, control, "% de pacientes controladas", ACCENT, "{c} %")).classes(
                "col"
            ).style("height: 220px")
            ui.echart(_bars(names, healthy, "Células sanas destruidas (mediana)", "#4b5563", "{c}")).classes(
                "col"
            ).style("height: 220px")


def _bars(names: list[str], values: list[int], title: str, color: str, label: str) -> dict:
    return {
        "animation": False,
        "title": {"text": title, "textStyle": {"fontSize": 12, "fontWeight": 600, "color": "#374151"}},
        "grid": {"left": 190, "right": 50, "top": 30, "bottom": 6},
        "xAxis": {"type": "value", "show": False},
        "yAxis": {
            "type": "category",
            "data": names,
            "axisTick": {"show": False},
            "axisLine": {"show": False},
            "axisLabel": {"fontSize": 11, "color": "#374151", "width": 180, "overflow": "truncate"},
        },
        "series": [
            {
                "type": "bar",
                "data": values,
                "color": color,
                "barWidth": 14,
                "label": {"show": True, "position": "right", "formatter": label, "fontSize": 11},
                "itemStyle": {"borderRadius": [0, 3, 3, 0]},
            }
        ],
    }


def _zoom(url: str, caption: str) -> None:
    with ui.dialog() as dialog, ui.card().classes("q-pa-md").style("max-width: 92vw"):
        ui.image(url).props("fit=contain").style("width: 86vw; max-height: 80vh")
        with ui.row().classes("w-full items-center justify-between"):
            ui.label(caption).classes("text-sm")
            secondary("Cerrar", dialog.close)
    dialog.open()
