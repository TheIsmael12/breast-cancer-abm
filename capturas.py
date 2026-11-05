"""Regenera las capturas de docs/capturas/ sin abrir ninguna ventana.

    python capturas.py

Todas salen del visor real (visor.py) con backend Agg, misma paciente (semilla 6, que
el día 120 tiene un carcinoma de ~900 células en 7 focos):
  01 · día 60, lesión inicial            02 · día 120, el tumor que se va a tratar
  03 · comparativa del tratamiento local el día después: tumor · cirugía · ablación en
       disco · ablación guiada (la misma paciente copiada cuatro veces)
  04 · gráfica de resultados de sin_cirugia.py (resultados/2026-11-05/resumen.csv)
  05 · 3D (tejido3d.py): el conducto en corte el día 45, el 80 y el 100 (arsenal el 80)
  06 · 3D: secciones transversales del día 100, de la superficie al núcleo

    python capturas.py --solo 3d      # solo las 3D (~2 min)
"""
from __future__ import annotations

import argparse
import copy
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from biologia import PASOS_POR_DIA  # noqa: E402
from visor import COLORES, Visor  # noqa: E402
from tejido import (FIBROBLASTO, HIPERPLASIA, IN_SITU, INICIADA, INVASIVA,  # noqa: E402
                    MIOEPITELIAL, NECROSIS, NORMAL, REINA)

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


def _estado3d(sim, pos, c) -> int:
    """El mismo código de colores que el visor 2D."""
    from tejido3d import ESTROMA, MIO
    if c.tipo == FIBROBLASTO:
        return -1
    if c.tipo == MIOEPITELIAL:
        return -2
    if c.necrosis:
        return NECROSIS
    if c.n_drivers == 0:
        return NORMAL
    if c.madre:
        return REINA
    if sim.zona[pos] in (MIO, ESTROMA):
        return INVASIVA
    if c.anoikis_res:
        return IN_SITU
    return HIPERPLASIA if c.senal_prolif else INICIADA


def _corte3d(ax, sim, titulo: str):
    """Medio cilindro (x <= eje) para ver dentro: lo que una célula terapéutica ve."""
    eje = sim.centro
    grupos: dict[int, list] = {}
    for pos, c in sim.celulas.items():
        z, y, x = pos
        if x > eje + 0.5:
            continue
        e = _estado3d(sim, pos, c)
        if e == -1:
            continue
        grupos.setdefault(e, []).append((x, y, z))
    estilo = {   # (color, tamaño, opacidad, etiqueta)
        -2: ("#8e7cc3", 3, 0.10, "mioepitelial"), NORMAL: ("#f4b6c8", 3, 0.18, "normal"),
        INICIADA: (COLORES[INICIADA], 5, 0.6, "iniciada"),
        HIPERPLASIA: (COLORES[HIPERPLASIA], 5, 0.6, "hiperplasia"),
        IN_SITU: (COLORES[IN_SITU], 6, 0.55, "carcinoma in situ"),
        INVASIVA: (COLORES[INVASIVA], 7, 0.8, "carcinoma invasivo"),
        REINA: ("#ffffff", 10, 1.0, "célula madre tumoral"),
        NECROSIS: (COLORES[NECROSIS], 6, 0.7, "necrosis"),
    }
    for e, (color, tam, alfa, nombre) in estilo.items():
        if e in grupos:
            xs, ys, zs = zip(*grupos[e])
            ax.scatter(xs, ys, zs, s=tam, c=color, alpha=alfa, linewidths=0,
                       label=f"{nombre} ({len(grupos[e])})", depthshade=False)
    ia = [a.pos for a in sim.agentes if a.tipo == "IA" and a.pos[2] <= eje + 0.5]
    if ia:
        zs, ys, xs = zip(*ia)
        ax.scatter(xs, ys, zs, s=26, c="#f72585", marker="*", linewidths=0,
                   label=f"célula terapéutica ({len(ia)})", depthshade=False)
    largo, lado, _ = sim.forma
    ax.set_xlim(0, lado), ax.set_ylim(0, lado), ax.set_zlim(0, largo)
    ax.set_box_aspect((lado / 2, lado, largo))
    ax.view_init(elev=18, azim=-28)
    ax.set_axis_off()
    r = sim.resumen()
    ax.set_title(f"{titulo}\ntumor {r['tumorales']:,} · peligrosas {r['peligrosas']:,} · "
                 f"accesible {sim.fraccion_accesible():.0%}".replace(",", "."), fontsize=10)
    ax.legend(loc="upper left", fontsize=7, frameon=False, markerscale=1.6)


def capturas_3d():
    from tejido3d import Simulacion3D

    def avanzar3d(sim, dias):
        for _ in range(int(dias * PASOS_POR_DIA)):
            sim.paso()

    plt.style.use("dark_background")
    sim = Simulacion3D("luminal", 1)
    fig = plt.figure(figsize=(18, 7.2))
    avanzar3d(sim, 45)
    _corte3d(fig.add_subplot(1, 3, 1, projection="3d"), sim, "Día 45 · lesión inicial")
    avanzar3d(sim, 35)
    _corte3d(fig.add_subplot(1, 3, 2, projection="3d"), sim,
             "Día 80 · la masa llena el conducto")
    sim.liberar_terapia("troyano+guardiana_nicho+memoria", 300)
    avanzar3d(sim, 20)
    _corte3d(fig.add_subplot(1, 3, 3, projection="3d"), sim,
             "Día 100 · 20 días después del arsenal")
    fig.suptitle("Modelo 3D (tejido3d.py): medio conducto mamario cortado a lo largo. "
                 "«Accesible» = células tumorales con un hueco al lado, por donde puede "
                 "entrar una terapia", fontsize=11)
    fig.tight_layout()
    fig.savefig(SALIDA / "05_3d_conducto_en_corte.png", dpi=85)
    plt.close(fig)
    print("  · 05_3d_conducto_en_corte.png")

    # Secciones transversales del mismo día 100: dónde están las terapéuticas.
    largo, lado, _ = sim.forma
    import numpy as np
    from matplotlib.colors import to_rgb
    fig, ejes = plt.subplots(1, 4, figsize=(16, 4.6))
    for ax, z in zip(ejes, (largo // 8, largo // 3, largo // 2, 3 * largo // 4)):
        img = np.full((lado, lado, 3), to_rgb("#231b17"))
        img[sim.zona[z] == 0] = to_rgb("#120c16")
        img[sim.vasos[z]] = to_rgb("#ff2d55")
        for (zz, y, x), c in sim.celulas.items():
            if zz != z:
                continue
            e = _estado3d(sim, (zz, y, x), c)
            img[y, x] = to_rgb({-1: COLORES[1], -2: "#8e7cc3"}.get(e, COLORES.get(e, "#ffffff")))
        ax.imshow(img, interpolation="nearest")
        pts = [(a.pos[2], a.pos[1]) for a in sim.agentes if a.tipo == "IA" and a.pos[0] == z]
        if pts:
            ax.scatter(*zip(*pts), s=40, marker="*", c="#f72585")
        ax.set_title(f"sección z = {z} · terapéuticas en el plano: {len(pts)}", fontsize=10)
        ax.set_xticks([]), ax.set_yticks([])
    fig.suptitle("Día 100 en 3D, secciones transversales del conducto (rojo: carcinoma · "
                 "estrellas: células terapéuticas)", fontsize=11)
    fig.tight_layout()
    fig.savefig(SALIDA / "06_3d_secciones.png", dpi=85)
    plt.close(fig)
    print("  · 06_3d_secciones.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solo", choices=("2d", "3d", "grafica"), default=None)
    a = ap.parse_args()
    SALIDA.mkdir(parents=True, exist_ok=True)
    print("Capturas en", SALIDA)
    if a.solo in (None, "2d"):
        capturas_visor()
    if a.solo in (None, "grafica"):
        grafica_resultados()
    if a.solo in (None, "3d"):
        capturas_3d()
