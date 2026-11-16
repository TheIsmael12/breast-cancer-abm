"""Rutas del proyecto, iguales desde el repositorio, desde una instalación y desde el ejecutable.

El código vive en ``src/breast_cancer_abm``; los datos, resultados y documentos en la raíz
del repositorio. En el ejecutable (PyInstaller) todo va empaquetado en ``sys._MEIPASS``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _repo_root() -> Path:
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
    env = os.environ.get("BCABM_ROOT")
    if env:
        return Path(env).resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pyproject.toml").exists():
            return parent
    return Path.cwd()


ROOT = _repo_root()
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
CACHE_3D_DIR = RESULTS_DIR / "cache_3d"
DOCS_DIR = ROOT / "docs"
FIGURES_DIR = DOCS_DIR / "figures"
REPORTS_DIR = DOCS_DIR / "reports"
PROGRESS_DIR = DOCS_DIR / "progress"
LOGBOOK = DOCS_DIR / "LOGBOOK.md"

# Experimentos del laboratorio: en el ejecutable, en Documentos (la carpeta interna es temporal).
LAB_DIR = (
    (Path.home() / "Documents" / "BreastCancerABM" / "lab")
    if getattr(sys, "frozen", False)
    else RESULTS_DIR / "lab"
)
