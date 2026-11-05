"""Prueba de humo: que el motor arranca y que cada intervención devuelve lo que promete.

No mide resultados científicos (eso es validacion.py, que tarda): solo que nada se rompe.
    python -m pytest -q tests
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from biologia import PASOS_POR_DIA, Parametros  # noqa: E402
from tejido import Simulacion  # noqa: E402


def _sim(dias: int = 60, semilla: int = 6) -> Simulacion:
    sim = Simulacion("luminal", semilla, Parametros(prob_neoantigeno=0.1))
    for _ in range(dias * PASOS_POR_DIA):
        sim.paso()
    return sim


def test_reproducible():
    assert _sim(30, 3).resumen() == _sim(30, 3).resumen()


def test_tumor_aparece():
    assert _sim(90).resumen()["tumorales"] > 0


def test_intervenciones_locales():
    base = _sim(90)
    for accion in ("resecar", "disco", "guiada"):
        import copy
        sim = copy.deepcopy(base)
        if accion == "resecar":
            r = sim.resecar()
        else:
            r = sim.ablacionar(guiada=accion == "guiada", borde=1)
        assert set(("sanas", "residuo")) <= set(r)
        assert r["sanas"] >= 0 and r["residuo"] >= 0
        for _ in range(2 * PASOS_POR_DIA):
            sim.paso()


def test_pruebas_de_residuo():
    sim = _sim(90)
    b = sim.biopsia(6)
    assert b["muestreadas"] > 0 and b["real"] >= 0
    c = sim.ctdna(0.1)
    assert c["fragmentos"] >= 0 and c["esperado"] >= 0


def test_terapia_y_3d():
    sim = _sim(40)
    sim.liberar_terapia("troyano+guardiana_nicho+memoria", 30)
    for _ in range(PASOS_POR_DIA):
        sim.paso()
    from tejido3d import Simulacion3D
    s3 = Simulacion3D("luminal", 1, lado=30, largo=12, radio=8.0)
    for _ in range(PASOS_POR_DIA):
        s3.paso()
    assert s3.resumen()["normales"] > 0
