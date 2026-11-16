"""Línea de comandos ``bcabm``: un único punto de entrada para todo el proyecto.

bcabm app                              # la aplicación (ventana propia o navegador)
bcabm viewer --scenario her2           # el visor clásico de matplotlib
bcabm study --list                     # los estudios disponibles
bcabm study no-surgery-3d --seeds 8 --processes 2
bcabm figures --only 3d                # regenera docs/figures
"""

from __future__ import annotations

import argparse
import sys

from breast_cancer_abm import __version__


def _viewer(argv: list[str]) -> int:
    from breast_cancer_abm.model.biology import SCENARIOS
    from breast_cancer_abm.model.policies import THERAPIES

    ap = argparse.ArgumentParser(
        prog="bcabm viewer",
        description="Visor clásico de matplotlib: Play, paso a paso, "
        "terapias y clic para inspeccionar células.",
    )
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--seed", type=int, default=None, help="fija el azar para repetir una partida")
    ap.add_argument(
        "--therapy", default=None, help=f"se libera sola en --therapy-day; admite '+': {', '.join(THERAPIES)}"
    )
    ap.add_argument("--therapy-day", type=float, default=150)
    a = ap.parse_args(argv)
    if a.therapy and not all(t in THERAPIES for t in a.therapy.split("+")):
        ap.error(f"terapia desconocida: {a.therapy}")
    from breast_cancer_abm.viz.viewer import Viewer

    Viewer(a.scenario, a.seed, a.therapy, a.therapy_day).show()
    return 0


def _study(argv: list[str]) -> int:
    from breast_cancer_abm.studies import STUDIES

    if not argv or argv[0] in ("-l", "--list", "-h", "--help"):
        print("Estudios disponibles (bcabm study <nombre> --help para sus opciones):\n")
        topic = None
        for s in sorted(STUDIES.values(), key=lambda s: (s.topic, s.name)):
            if s.topic != topic:
                topic = s.topic
                print(f"  {topic}")
            print(f"    {s.name:<16} {s.title}")
        return 0
    name, *rest = argv
    if name not in STUDIES:
        print(f"Estudio desconocido: {name}. Lista: bcabm study --list", file=sys.stderr)
        return 2
    STUDIES[name].run(rest)
    return 0


def _figures(argv: list[str]) -> int:
    from breast_cancer_abm.viz.figures import main as figures_main

    figures_main(argv)
    return 0


def _app(argv: list[str]) -> int:
    from breast_cancer_abm.app.main import main as app_main

    return app_main(argv)


COMMANDS = {
    "app": (_app, "abre la aplicación (simulación 2D/3D, resultados, avances y bitácora)"),
    "viewer": (_viewer, "visor clásico de matplotlib"),
    "study": (_study, "ejecuta un estudio (bcabm study --list)"),
    "figures": (_figures, "regenera las figuras de docs/figures"),
}


def _utf8_output() -> None:
    """En Windows la consola usa cp1252 y rompe al redirigir a fichero (flechas, ·, µ)."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, OSError, ValueError):
            pass


def main(argv: list[str] | None = None) -> int:
    _utf8_output()
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__.strip(), "\n\nComandos:")
        for name, (_, help_) in COMMANDS.items():
            print(f"  {name:<9} {help_}")
        return 0
    if argv[0] in ("-V", "--version"):
        print(f"breast-cancer-abm {__version__}")
        return 0
    command, *rest = argv
    if command not in COMMANDS:
        print(f"Comando desconocido: {command}. Ayuda: bcabm --help", file=sys.stderr)
        return 2
    return COMMANDS[command][0](rest)


if __name__ == "__main__":
    raise SystemExit(main())
