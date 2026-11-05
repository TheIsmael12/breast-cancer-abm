"""Regenera las capturas de docs/capturas/ sin abrir ninguna ventana.

    python capturas.py

Todas salen del visor real (visor.py) con backend Agg, misma paciente (semilla 6, que
el día 120 tiene un carcinoma de ~900 células en 7 focos):
  01 · día 60, lesión inicial            02 · día 120, el tumor que se va a tratar
  03 · comparativa del tratamiento local el día después: tumor · cirugía · ablación en
       disco · ablación guiada (la misma paciente copiada cuatro veces)
  04 · gráfica de resultados de sin_cirugia.py (resultados/2026-11-05/resumen.csv)
"""
from __future__ import annotations

import copy
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from biologia import PASOS_POR_DIA  # noqa: E402
from visor import Visor  # noqa: E402

RAIZ = Path(__file__).parent
SALIDA = RAIZ / "docs" / "capturas"
SEMILLA = 6


def avanzar(v: Visor, dias: float):
    for _ in range(int(dias * PASOS_POR_DIA)):
        v.sim.paso()


def foto(v: Visor, nombre: str, titulo: str | None = None):
    v._dibujar()
    if titulo:
        v.ax.set_title(titulo, fontsize=10)
    v.fig.savefig(SALIDA / nombre, dpi=90)
    print("  ·", nombre)


def capturas_visor():
    v = Visor("luminal", SEMILLA, None, None)
    avanzar(v, 60)
    foto(v, "01_dia_060_lesion_inicial.png")
    avanzar(v, 60)
    foto(v, "02_dia_120_tumor.png")
    base = v.sim
    casos = (
        ("Tumor el día 120", None),
        ("Cirugía (resecar)", lambda s: s.resecar()),
        ("Ablación en disco", lambda s: s.ablacionar(supervivencia_borde=0.3)),
        ("Ablación guiada por imagen",
         lambda s: s.ablacionar(supervivencia_borde=0.3, guiada=True, borde=1)),
    )
    plt.style.use("dark_background")
    fig, ejes = plt.subplots(1, 4, figsize=(16, 4.9))
    for ax, (titulo, accion) in zip(ejes, casos):
        v.sim = copy.deepcopy(base)
        sub = ""
        if accion:
            r = accion(v.sim)
            avanzar(v, 1)
            sub = f" · sanas destruidas: {r['sanas']}"
        ax.imshow(v._rgb(), interpolation="nearest")
        ax.set_title(titulo + sub, fontsize=11)
        ax.set_xticks([]), ax.set_yticks([])
    fig.suptitle("Misma paciente, el día después del tratamiento local (rojo: carcinoma · "
                 "gris oscuro: tejido destruido · hueco: tejido extirpado)", fontsize=11)
    fig.tight_layout()
    fig.savefig(SALIDA / "03_comparativa_tratamiento_local.png", dpi=90)
    plt.close(fig)
    print("  · 03_comparativa_tratamiento_local.png")
    plt.close(v.fig)


def grafica_resultados():
    """Dos paneles (dos escalas distintas → dos gráficas, nunca doble eje)."""
    filas = list(csv.DictReader(open(RAIZ / "resultados" / "2026-11-05" / "resumen.csv",
                                     encoding="utf-8")))
    filas.reverse()   # la referencia arriba
    etiquetas = [f["brazo"] for f in filas]
    control = [int(f["controladas"]) for f in filas]
    erradica = [int(f["erradicadas"]) for f in filas]
    sanas = [float(f["sanas_destruidas"]) for f in filas]
    n = int(filas[0]["n"])

    plt.style.use("default")
    tinta, tinta2, rejilla, serie = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 0.5 * len(filas) + 2.2), sharey=True,
                                 gridspec_kw=dict(wspace=0.08))
    fig.patch.set_facecolor("#fcfcfb")
    y = range(len(filas))
    a1.barh(y, control, color="#b7d3f6", height=0.62, label="controladas (<100 peligrosas)")
    a1.barh(y, erradica, color=serie, height=0.62, label="erradicadas (0)")
    for i, (c, e) in enumerate(zip(control, erradica)):
        a1.text(c + 0.3, i, f"{c}/{n} · {e} a 0", va="center", fontsize=8.5, color=tinta2)
    a1.set_xlim(0, n + 7)
    a1.set_title("Pacientes controladas y erradicadas (de 24)", loc="left", fontsize=11,
                 color=tinta, pad=10)
    a1.legend(loc="upper left", bbox_to_anchor=(0, -0.06), ncol=2, frameon=False, fontsize=9)
    a2.barh(y, sanas, color=serie, height=0.62)
    for i, s in enumerate(sanas):
        a2.text(s + 4, i, f"{s:.0f}", va="center", fontsize=8.5, color=tinta2)
    a2.set_xlim(0, max(sanas) * 1.15)
    a2.set_title("Células sanas destruidas (mediana)", loc="left", fontsize=11, color=tinta,
                 pad=10)
    a1.set_yticks(list(y), etiquetas, fontsize=9, color=tinta)
    for a in (a1, a2):
        a.set_facecolor("#fcfcfb")
        a.grid(axis="x", color=rejilla, linewidth=0.8)
        a.set_axisbelow(True)
        for lado in ("top", "right", "left"):
            a.spines[lado].set_visible(False)
        a.spines["bottom"].set_color(rejilla)
        a.tick_params(axis="x", colors=tinta2, labelsize=8.5)
        a.tick_params(axis="y", length=0)
    fig.suptitle("¿Se puede tratar sin operar? · luminal, tratamiento día 120, seguimiento "
                 "500 días, 24 pacientes simuladas", x=0.01, ha="left", fontsize=12,
                 color=tinta)
    fig.subplots_adjust(left=0.25, right=0.98, top=0.86, bottom=0.12, wspace=0.08)
    fig.savefig(SALIDA / "04_resultados_sin_cirugia.png", dpi=110)
    plt.close(fig)
    print("  · 04_resultados_sin_cirugia.png")


if __name__ == "__main__":
    SALIDA.mkdir(parents=True, exist_ok=True)
    print("Capturas en", SALIDA)
    capturas_visor()
    grafica_resultados()
