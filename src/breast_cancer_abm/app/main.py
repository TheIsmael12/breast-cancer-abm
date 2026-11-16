"""Arranque de la aplicación: ``bcabm app`` (ventana propia) o ``bcabm app --browser``."""

from __future__ import annotations

import argparse
import multiprocessing

from nicegui import app, ui

from breast_cancer_abm.app import theme
from breast_cancer_abm.app.pages import docs, home, lab, results, simulation2d, simulation3d, studies
from breast_cancer_abm.paths import DOCS_DIR


def _register() -> None:
    theme.setup()
    app.add_static_files("/files/docs", DOCS_DIR)
    ui.page("/")(home.page)
    ui.page("/simulation")(simulation2d.page)
    ui.page("/simulation-3d")(simulation3d.page)
    ui.page("/lab")(lab.page)
    ui.page("/results")(results.page)
    ui.page("/studies")(studies.page)
    ui.page("/docs")(docs.page)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="bcabm app", description=__doc__)
    ap.add_argument("--browser", action="store_true", help="abre el navegador en vez de una ventana propia")
    ap.add_argument("--port", type=int, default=None, help="puerto local (por defecto, uno libre)")
    ap.add_argument("--no-show", action="store_true", help="no abre nada (útil para pruebas)")
    a = ap.parse_args(argv)
    _register()
    native = not a.browser and not a.no_show
    ui.run(
        title="Breast Cancer ABM",
        favicon=theme.LOGO_SVG,
        native=native,
        window_size=(1500, 940) if native else None,
        port=a.port or (None if native else 8765),
        show=not a.no_show,
        reload=False,
        language="es",
        show_welcome_message=False,
        reconnect_timeout=60.0,
    )
    return 0


if __name__ in {"__main__", "__mp_main__"}:
    multiprocessing.freeze_support()
    main()
