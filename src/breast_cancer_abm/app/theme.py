"""Sistema de diseño de la aplicación: estructura de página, tarjetas, botones y cifras.

Todas las páginas se construyen con las mismas piezas para que se vean coherentes:

- ``shell``: navegación lateral + cabecera de página (título, subtítulo y acciones).
- ``card``: el único contenedor de contenido (borde fino, título claro, cuerpo con aire).
- ``primary`` / ``secondary`` / ``ghost``: tres niveles de botón, siempre iguales.
- ``kpi_strip``: la tira de cifras de una simulación.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass

from nicegui import app, ui

from breast_cancer_abm import __version__

ACCENT = "#be185d"  # rosa del lazo, en tono profundo para que sea legible
ACCENT_SOFT = "#fdf2f8"
INK = "#111827"
MUTED = "#6b7280"
BORDER = "#e5e7eb"
CANVAS_BG = "#0f0b12"
REPO = "https://github.com/TheIsmael12/breast-cancer-abm"

# Marca: corte de un conducto (capa mioepitelial, epitelio y un foco tumoral en la luz).
LOGO_SVG = """<svg viewBox="0 0 64 64" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Breast Cancer ABM">
<circle cx="32" cy="32" r="29" fill="#231b17"/>
<circle cx="32" cy="32" r="22" fill="none" stroke="#8e7cc3" stroke-width="4"/>
<circle cx="32" cy="32" r="16.5" fill="none" stroke="#f4b6c8" stroke-width="5" stroke-dasharray="4.3 2"/>
<circle cx="32" cy="32" r="11" fill="#120c16"/>
<circle cx="27" cy="29" r="4.2" fill="#e63946"/><circle cx="33.5" cy="34" r="3.2" fill="#9d0208"/>
<circle cx="31" cy="26.5" r="2.4" fill="#ff9f1c"/><circle cx="29.5" cy="33.5" r="1.6" fill="#ffffff"/>
</svg>"""

NAV = (
    ("", (("/", "home", "Inicio"),)),
    (
        "Explorar",
        (("/simulation", "blur_on", "Simulación 2D"), ("/simulation-3d", "view_in_ar", "Simulación 3D")),
    ),
    (
        "Investigar",
        (
            ("/lab", "biotech", "Laboratorio"),
            ("/results", "insights", "Resultados"),
            ("/studies", "science", "Estudios"),
        ),
    ),
    ("Proyecto", (("/docs", "menu_book", "Documentación"),)),
)

CSS = f"""
:root {{ --accent: {ACCENT}; --accent-soft: {ACCENT_SOFT}; --ink: {INK}; --muted: {MUTED};
         --border: {BORDER}; --canvas: {CANVAS_BG}; }}
html, body {{ background: #ffffff; color: var(--ink); font-size: 14px;
        font-family: Inter, "Segoe UI", system-ui, -apple-system, Roboto, sans-serif;
        -webkit-font-smoothing: antialiased; }}
* {{ scrollbar-width: thin; scrollbar-color: var(--accent) transparent; }}
::-webkit-scrollbar {{ width: 6px; height: 6px; }}
::-webkit-scrollbar-thumb {{ background: var(--accent); border-radius: 3px; }}
::-webkit-scrollbar-track {{ background: transparent; }}

/* Quasar a la medida: plano, radios de 8 px como máximo, sin sombras */
.q-card, .q-btn, .q-field__control, .q-menu, .q-dialog__inner > div, .q-table__container,
.q-tab, .q-item, .q-expansion-item {{ border-radius: 8px !important; }}
.q-card {{ box-shadow: none !important; border: 1px solid var(--border); }}
.q-btn {{ text-transform: none; font-weight: 600; letter-spacing: 0; min-height: 36px; }}
.q-btn .q-btn__content {{ flex-wrap: nowrap; white-space: nowrap; }}
.q-btn--dense {{ min-height: 32px; }}
.q-btn--unelevated, .q-btn--outline, .q-btn--flat {{ box-shadow: none; }}
.btn-secondary.q-btn--outline:before {{ border-color: #d1d5db; }}
.q-field--outlined .q-field__control:before {{ border-color: #d1d5db; }}
.q-field--dense .q-field__control {{ min-height: 38px; }}
.q-toggle__track, .q-badge, .q-chip {{ border-radius: 999px !important; }}
.q-table th {{ font-weight: 600; color: var(--muted); font-size: 12px; }}
.q-table td {{ font-size: 13px; }}
.q-drawer {{ border-right: 1px solid var(--border) !important; }}
.q-expansion-item__container > .q-item {{ padding: 8px 4px; min-height: 40px; }}
.q-stepper {{ box-shadow: none; background: transparent; }}
.q-stepper__header {{ border-bottom: 1px solid var(--border); }}
.q-stepper__tab {{ padding: 12px 8px; }}
.q-stepper__step-inner {{ padding: 24px 0 8px !important; }}
.q-stepper__nav {{ padding: 16px 0 0 !important; }}

/* Estructura de página */
.page {{ width: 100%; max-width: 1440px; margin: 0 auto; padding: 28px 32px 48px; gap: 24px; }}
.page-head {{ width: 100%; display: flex; align-items: flex-end; justify-content: space-between; gap: 16px;
              padding-bottom: 20px; border-bottom: 1px solid var(--border); }}
.page-title {{ font-size: 1.6rem; font-weight: 700; letter-spacing: -0.015em; line-height: 1.2; }}
.page-sub {{ color: var(--muted); font-size: 0.95rem; max-width: 820px; margin-top: 6px; line-height: 1.5; }}
.muted {{ color: var(--muted); }}

/* Tarjeta: el único contenedor */
.card {{ border: 1px solid var(--border); border-radius: 8px; background: #fff; width: 100%; min-width: 0; }}
.card-head {{ display: flex; align-items: center; justify-content: space-between; gap: 12px;
              padding: 14px 18px 0; }}
.card-title {{ font-size: 0.98rem; font-weight: 650; display: flex; align-items: center; gap: 8px; }}
.card-title .q-icon {{ color: var(--accent); }}
.card-sub {{ color: var(--muted); font-size: 0.82rem; padding: 2px 18px 0; line-height: 1.45; }}
.card-body {{ padding: 14px 18px 18px; gap: 12px; width: 100%; }}
.card-body.flush {{ padding: 0; }}
.field-label {{ font-size: 0.8rem; font-weight: 600; color: #374151; }}
.hint {{ font-size: 0.78rem; color: var(--muted); line-height: 1.45; }}
.divider {{ height: 1px; background: var(--border); width: 100%; }}

/* Alias de compatibilidad: mismo aspecto que .card */
.panel {{ border: 1px solid var(--border); border-radius: 8px; background: #fff; padding: 16px 18px; }}
.section-title {{ font-size: 0.92rem; font-weight: 650; color: var(--ink); }}

/* Barra de herramientas de las simulaciones */
.toolbar {{ display: flex; flex-wrap: wrap; align-items: center; gap: 12px 16px; width: 100%;
            border: 1px solid var(--border); border-radius: 8px; padding: 12px 16px; }}
.toolbar-group {{ display: flex; align-items: center; gap: 8px; }}
.toolbar-sep {{ width: 1px; align-self: stretch; background: var(--border); }}

/* Navegación */
.nav-section {{ font-size: 0.7rem; font-weight: 700; color: #9ca3af; text-transform: uppercase;
                letter-spacing: 0.07em; padding: 16px 12px 4px; }}
.nav-item {{ display: flex; align-items: center; gap: 12px; width: 100%; box-sizing: border-box;
             padding: 9px 12px; border-radius: 8px; color: #374151; font-weight: 500; font-size: 0.92rem;
             text-decoration: none; line-height: 1.3; }}
.nav-item:hover {{ background: #f9fafb; color: var(--ink); }}
.nav-item.active {{ background: var(--accent-soft); color: var(--accent); font-weight: 600; }}
.nav-item .q-icon {{ flex: none; }}
.nav-list {{ width: 100%; display: flex; flex-direction: column; gap: 2px; }}
/* Navegación de sección (documentación, resultados) */
.sidepanel {{ width: 280px; flex: none; position: sticky; top: 16px; border: 1px solid var(--border);
              border-radius: 8px; padding: 6px 8px 10px; background: #fff; max-height: calc(100vh - 32px);
              overflow-y: auto; }}
.sidepanel .nav-section {{ padding: 12px 10px 4px; }}
.sidepanel .nav-item {{ padding: 8px 10px; font-size: 0.86rem; align-items: flex-start; gap: 10px; }}
.sidepanel .nav-date {{ font-size: 0.72rem; color: #9ca3af; font-weight: 500; }}
.nav-item.active .nav-date {{ color: var(--accent); opacity: .75; }}
/* Simulaciones: lienzo a la izquierda, panel fijo a la derecha */
.sim-grid {{ display: grid; grid-template-columns: minmax(420px, 1fr) 420px; gap: 24px; width: 100%;
             align-items: start; }}
.sim-stage {{ width: 100%; max-width: 820px; margin: 0 auto; }}
.kpi-compact .kpi {{ padding: 8px 12px; }}
.kpi-compact .kpi-value {{ font-size: 1.2rem; }}
.brand {{ display: flex; align-items: center; gap: 10px; padding: 4px 6px 8px; }}
.brand svg {{ width: 34px; height: 34px; }}
.brand-name {{ font-weight: 700; font-size: 1.02rem; line-height: 1.1; }}
.brand-sub {{ color: var(--muted); font-size: 0.75rem; }}

/* Cifras */
.kpi-strip {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; width: 100%; }}
.kpi {{ border: 1px solid var(--border); border-radius: 8px; padding: 12px 16px; background: #fff; min-width: 0; }}
.kpi-label {{ color: var(--muted); font-size: 0.78rem; font-weight: 500; white-space: nowrap;
              overflow: hidden; text-overflow: ellipsis; }}
.kpi-value {{ font-size: 1.5rem; font-weight: 700; font-variant-numeric: tabular-nums; letter-spacing: -0.01em;
              line-height: 1.25; }}
.kpi-value.accent {{ color: var(--accent); }}

/* Leyendas */
.swatch {{ width: 10px; height: 10px; border-radius: 3px; flex: none; box-shadow: inset 0 0 0 1px rgba(0,0,0,.15); }}
.legend {{ display: flex; flex-wrap: wrap; gap: 6px 14px; }}
.legend-item {{ display: flex; align-items: center; gap: 6px; font-size: 0.78rem; color: #374151; }}

/* Lienzos */
.tissue-canvas {{ width: 100%; }}
.tissue-canvas-wrap {{ position: relative; width: 100%; aspect-ratio: 1 / 1; background: var(--canvas);
                       border-radius: 6px; overflow: hidden; cursor: crosshair; }}
.tissue-canvas-wrap canvas {{ display: block; }}
.tissue-canvas-hover {{ position: absolute; right: 8px; bottom: 8px; color: #e5e7eb; font-size: 0.75rem;
                        font-variant-numeric: tabular-nums; background: rgba(0,0,0,.55);
                        padding: 2px 6px; border-radius: 4px; }}
.scene-wrap {{ background: var(--canvas); border-radius: 6px; overflow: hidden; width: 100%; }}

/* Tarjetas de elección (plantillas, modelo, tipo de cáncer) */
.choice {{ border: 1px solid var(--border); border-radius: 8px; padding: 14px 16px; cursor: pointer;
           background: #fff; transition: border-color .15s, background .15s; min-width: 0; }}
.choice:hover {{ border-color: #f9a8d4; }}
.choice.active {{ border-color: var(--accent); background: var(--accent-soft); }}
.choice-title {{ font-weight: 600; font-size: 0.9rem; display: flex; align-items: center; gap: 8px; }}
.choice-title .q-icon {{ color: var(--accent); }}
.choice-text {{ color: var(--muted); font-size: 0.8rem; line-height: 1.45; margin-top: 4px; }}
.choice-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(230px, 1fr)); gap: 12px; width: 100%; }}

.mono {{ font-family: "Cascadia Code", "JetBrains Mono", Consolas, monospace; font-size: 0.8rem; }}
.raw-output {{ white-space: pre; font-family: "Cascadia Code", Consolas, monospace; font-size: 0.78rem;
               line-height: 1.45; color: #1f2937; padding: 4px 18px 18px; }}
.event-line {{ font-size: 0.8rem; color: #374151; font-variant-numeric: tabular-nums; line-height: 1.5; }}
.command {{ border: 1px solid var(--border); border-radius: 6px; padding: 6px 10px; }}

/* Documentos */
.doc {{ max-width: 860px; line-height: 1.7; font-size: 0.95rem; }}
.doc h1 {{ font-size: 1.6rem; font-weight: 700; margin: 0 0 14px; line-height: 1.25; }}
.doc h2 {{ font-size: 1.2rem; font-weight: 700; margin: 30px 0 10px; }}
.doc h3 {{ font-size: 1.02rem; font-weight: 700; margin: 22px 0 8px; }}
.doc p, .doc ul, .doc ol {{ margin: 0 0 12px; }}
.doc ul {{ list-style: disc; padding-left: 22px; }}
.doc ol {{ list-style: decimal; padding-left: 22px; }}
.doc a {{ color: var(--accent); text-decoration: none; font-weight: 500; }}
.doc a:hover {{ text-decoration: underline; }}
.doc img {{ max-width: 100%; border-radius: 8px; border: 1px solid var(--border); margin: 8px 0; }}
.doc code {{ background: #fdf2f8; color: #9d174d; padding: 1px 5px; border-radius: 4px; font-size: 0.86em; }}
.doc pre {{ border: 1px solid var(--border); border-radius: 8px; padding: 12px 14px; overflow-x: auto; }}
.doc pre code {{ background: none; color: var(--ink); padding: 0; }}
.doc table {{ border-collapse: collapse; margin: 10px 0 16px; font-size: 0.86rem; display: block; overflow-x: auto; }}
.doc th, .doc td {{ border: 1px solid var(--border); padding: 6px 10px; text-align: left; vertical-align: top; }}
.doc th {{ font-weight: 600; }}
.doc blockquote {{ border: 1px solid var(--border); border-radius: 8px; padding: 10px 14px; color: #374151;
                   margin: 0 0 12px; }}
"""


def setup() -> None:
    """Configuración global (una vez por proceso)."""
    app.colors(
        primary=ACCENT,
        secondary="#4b5563",
        accent=ACCENT,
        positive="#15803d",
        negative="#b91c1c",
        warning="#b45309",
        info="#1d4ed8",
    )


@dataclass
class Page:
    body: ui.column
    actions: ui.row


@contextmanager
def shell(active: str, title: str, subtitle: str = "") -> Iterator[Page]:
    """Estructura de cada página: navegación lateral y cabecera con título y acciones."""
    ui.add_css(CSS)
    ui.page_title(f"{title} · Breast Cancer ABM")
    with ui.left_drawer(value=True, fixed=True).props("width=232 bordered=false").classes("bg-white q-pa-md"):
        with ui.element("div").classes("brand"):
            ui.html(LOGO_SVG, sanitize=False)
            with ui.column().classes("gap-0"):
                ui.label("Breast Cancer ABM").classes("brand-name")
                ui.label("Simulador por agentes").classes("brand-sub")
        with ui.element("nav").classes("nav-list"):
            for section, items in NAV:
                if section:
                    ui.label(section).classes("nav-section")
                for path, icon, label in items:
                    with ui.link(target=path).classes("nav-item" + (" active" if path == active else "")):
                        ui.icon(icon, size="20px")
                        ui.label(label)
        ui.space()
        with ui.column().classes("w-full gap-1 q-mt-lg q-px-sm"):
            ui.label(
                "Desarrollado por Claude (Anthropic), dirigido por Ismael Ben. "
                "Educativo: no sirve para decisiones médicas."
            ).classes("hint")
            with ui.row().classes("items-center gap-2 q-mt-xs"):
                ui.link("GitHub", REPO, new_tab=True).classes("text-xs").style(f"color:{ACCENT}")
                ui.label(f"v{__version__}").classes("text-xs muted")
    with ui.column().classes("page") as page:
        with ui.element("div").classes("page-head"):
            with ui.column().classes("gap-0"):
                ui.label(title).classes("page-title")
                if subtitle:
                    ui.label(subtitle).classes("page-sub")
            actions = ui.row().classes("items-center gap-2 no-wrap")
        yield Page(page, actions)


@contextmanager
def card(
    title: str = "",
    subtitle: str = "",
    icon: str = "",
    *,
    flush: bool = False,
    classes: str = "",
    actions: Callable[[], None] | None = None,
) -> Iterator[ui.column]:
    """Tarjeta estándar: título (con icono y acciones opcionales), subtítulo y cuerpo."""
    with ui.element("div").classes(f"card column {classes}"):
        if title:
            with ui.element("div").classes("card-head"):
                with ui.element("div").classes("card-title"):
                    if icon:
                        ui.icon(icon, size="18px")
                    ui.label(title)
                if actions:
                    with ui.row().classes("items-center gap-1 no-wrap"):
                        actions()
            if subtitle:
                ui.label(subtitle).classes("card-sub")
        with ui.column().classes("card-body" + (" flush" if flush else "")) as body:
            yield body


def primary(label: str, on_click: Callable | None = None, icon: str | None = None) -> ui.button:
    return ui.button(label, icon=icon, on_click=on_click).props("unelevated no-caps")


def secondary(label: str = "", on_click: Callable | None = None, icon: str | None = None) -> ui.button:
    return (
        ui.button(label, icon=icon, on_click=on_click)
        .props("outline no-caps color=grey-9")
        .classes("btn-secondary")
    )


def ghost(label: str = "", on_click: Callable | None = None, icon: str | None = None) -> ui.button:
    return ui.button(label, icon=icon, on_click=on_click).props("flat no-caps color=grey-8")


def kpi(label: str, value: str = "–", accent: bool = False) -> ui.label:
    """Tarjeta de cifra; devuelve la etiqueta del valor para actualizarla en sitio."""
    with ui.element("div").classes("kpi"):
        ui.label(label).classes("kpi-label")
        return ui.label(value).classes("kpi-value" + (" accent" if accent else ""))


def kpi_strip(labels: list[tuple[str, bool]]) -> list[ui.label]:
    """Tira de cifras del mismo tamaño; devuelve las etiquetas de valor en orden."""
    with ui.element("div").classes("kpi-strip"):
        return [kpi(label, accent=accent) for label, accent in labels]


def legend_item(color: str, label: str) -> None:
    with ui.element("div").classes("legend-item"):
        ui.element("span").classes("swatch").style(f"background:{color}")
        ui.label(label)


def choice(title: str, text: str, active: bool, on_click: Callable, icon: str = "") -> None:
    """Tarjeta seleccionable (plantilla, modelo, tipo de cáncer)."""
    with ui.element("div").classes("choice" + (" active" if active else "")).on("click", on_click):
        with ui.element("div").classes("choice-title"):
            ui.icon(icon or ("radio_button_checked" if active else "radio_button_unchecked"), size="18px")
            ui.label(title)
        if text:
            ui.label(text).classes("choice-text")


def fmt(n: float | int) -> str:
    """Número con separador de miles español."""
    if isinstance(n, float) and not n.is_integer():
        return f"{n:,.1f}".replace(",", " ").replace(".", ",").replace(" ", ".")
    return f"{int(n):,}".replace(",", ".")


@dataclass
class NavLink:
    target: str
    title: str
    active: bool = False
    date: str = ""
    icon: str = ""


def sidepanel(groups: list[tuple[str, list[NavLink]]]) -> None:
    """Panel lateral de sección: una sola tarjeta con encabezados y enlaces a todo el ancho."""
    with ui.element("aside").classes("sidepanel"):
        for title, links in groups:
            if title:
                ui.label(title).classes("nav-section")
            with ui.element("div").classes("nav-list"):
                for link in links:
                    with ui.link(target=link.target).classes("nav-item" + (" active" if link.active else "")):
                        if link.icon:
                            ui.icon(link.icon, size="18px")
                        with ui.column().classes("gap-0").style("min-width: 0"):
                            if link.date:
                                ui.label(link.date).classes("nav-date")
                            ui.label(link.title)
