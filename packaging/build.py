"""Construye el ejecutable de la aplicación con PyInstaller.

    pip install -e ".[app,build]"
    python packaging/build.py            # dist/BreastCancerABM.exe (un solo fichero)
    python packaging/build.py --onedir   # dist/BreastCancerABM/ (arranca más rápido)

Lleva dentro el modelo, la aplicación, los documentos, las figuras, los datos de TCGA y
los resultados publicados (no la caché 3D local).
"""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

import PyInstaller.__main__

ROOT = Path(__file__).resolve().parent.parent
NAME = "BreastCancerABM"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--onedir", action="store_true", help="carpeta en vez de un único .exe")
    a = ap.parse_args()

    sep = os.pathsep
    data = [(ROOT / "docs", "docs"), (ROOT / "data", "data")]
    data += [(d, f"results/{d.name}") for d in sorted((ROOT / "results").glob("20*")) if d.is_dir()]
    args = [
        str(ROOT / "packaging" / "launcher.py"),
        "--name", NAME,
        "--windowed",
        "--noconfirm",
        "--clean",
        "--onedir" if a.onedir else "--onefile",
        "--icon", str(ROOT / "packaging" / "icon.ico"),
        "--paths", str(ROOT / "src"),
        "--distpath", str(ROOT / "dist"),
        "--workpath", str(ROOT / "build"),
        "--specpath", str(ROOT / "build"),
        "--collect-all", "nicegui",
        "--collect-all", "webview",
        "--collect-submodules", "breast_cancer_abm",
        "--collect-data", "breast_cancer_abm",
        # el visor de matplotlib y las figuras no se usan en la aplicación
        "--exclude-module", "tkinter",
        "--exclude-module", "breast_cancer_abm.viz",
    ]
    for src, dest in data:
        args += ["--add-data", f"{src}{sep}{dest}"]
    PyInstaller.__main__.run(args)
    out = ROOT / "dist" / (NAME if a.onedir else f"{NAME}.exe")
    shutil.rmtree(ROOT / "build", ignore_errors=True)
    print(f"\nListo: {out}")


if __name__ == "__main__":
    main()
