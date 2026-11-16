"""Regenera las capturas de docs/capturas/ sin abrir ninguna ventana.

    bcabm figures

Todas salen del visor real (viewer.py) con backend Agg, misma paciente (semilla 6, que
el día 120 tiene un carcinoma de ~900 células en 7 focos):
  01 · día 60, lesión inicial            02 · día 120, el tumor que se va a tratar
  03 · comparativa del tratamiento local el día después: tumor · cirugía · ablación en
       disco · ablación guiada (la misma paciente copiada cuatro veces)
  04 · gráfica de resultados de no_surgery.py (results/2026-11-05/summary.csv)
  05 · 3D (tissue3d.py): el conducto en corte el día 45, el 80 y el 100 (arsenal el 80)
  06 · 3D: secciones transversales del día 100, de la superficie al núcleo
  07 · 3D, no_surgery_3d.py: la paciente 4 el día 90 (arsenal el 60, ctDNA positivo)
       copiada tres veces: sin ablación · ablación en disco · ablación guiada

    bcabm figures --only ablation3d   # solo la 07 (~5 min)

    bcabm figures --only 3d      # solo las 3D (~2 min)
"""

from __future__ import annotations

import argparse
import copy
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from breast_cancer_abm.model.biology import STEPS_PER_DAY  # noqa: E402
from breast_cancer_abm.model.tissue import (  # noqa: E402
    FIBROBLAST,
    HYPERPLASIA,
    IN_SITU,
    INITIATED,
    INVASIVE,
    MYOEPITHELIAL,
    NECROSIS,
    NORMAL,
    QUEEN,
)
from breast_cancer_abm.paths import FIGURES_DIR, RESULTS_DIR  # noqa: E402
from breast_cancer_abm.viz.viewer import COLORS, Viewer  # noqa: E402

OUTPUT = FIGURES_DIR
SEED = 6


def advance(v: Viewer, days: float):
    for _ in range(int(days * STEPS_PER_DAY)):
        v.sim.step()


def snapshot(v: Viewer, name: str, title: str | None = None):
    v._draw()
    if title:
        v.ax.set_title(title, fontsize=10)
    v.fig.savefig(OUTPUT / name, dpi=90)
    print("  ·", name)


def viewer_figures():
    v = Viewer("luminal", SEED, None, None)
    advance(v, 60)
    snapshot(v, "01_day_060_initial_lesion.png")
    advance(v, 60)
    snapshot(v, "02_day_120_tumor.png")
    base = v.sim
    cases = (
        ("Tumor el día 120", None),
        ("Cirugía (resecar)", lambda s: s.resect()),
        ("Ablación en disco", lambda s: s.ablate(rim_survival=0.3)),
        ("Ablación guiada por imagen", lambda s: s.ablate(rim_survival=0.3, guided=True, rim=1)),
    )
    plt.style.use("dark_background")
    fig, panels = plt.subplots(1, 4, figsize=(16, 4.9))
    for ax, (title, act) in zip(panels, cases):
        v.sim = copy.deepcopy(base)
        sub = ""
        if act:
            r = act(v.sim)
            advance(v, 1)
            sub = f" · sanas destruidas: {r['healthy']}"
        ax.imshow(v._rgb(), interpolation="nearest")
        ax.set_title(title + sub, fontsize=11)
        ax.set_xticks([]), ax.set_yticks([])
    fig.suptitle(
        "Misma paciente, el día después del tratamiento local (rojo: carcinoma · "
        "gris oscuro: tejido destruido · hueco: tejido extirpado)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT / "03_local_treatment_comparison.png", dpi=90)
    plt.close(fig)
    print("  · 03_local_treatment_comparison.png")
    plt.close(v.fig)


def results_chart():
    """Dos paneles (dos escalas distintas → dos gráficas, nunca doble eje)."""
    rows = list(csv.DictReader(open(RESULTS_DIR / "2026-11-05" / "summary.csv", encoding="utf-8")))
    rows = [f for f in rows if not f["arm"].startswith("3D")]  # solo el banco 2D (n=24)
    rows.reverse()  # la referencia arriba
    tags = [f["arm"] for f in rows]
    control = [int(f["controlled"]) for f in rows]
    eradicates = [int(f["eradicated"]) for f in rows]
    healthy = [float(f["healthy_destroyed"]) for f in rows]
    n = int(rows[0]["n"])

    plt.style.use("default")
    ink, ink2, lattice, series = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
    fig, (a1, a2) = plt.subplots(
        1, 2, figsize=(14, 0.5 * len(rows) + 2.2), sharey=True, gridspec_kw=dict(wspace=0.08)
    )
    fig.patch.set_facecolor("#fcfcfb")
    y = range(len(rows))
    a1.barh(y, control, color="#b7d3f6", height=0.62, label="controladas (<100 peligrosas)")
    a1.barh(y, eradicates, color=series, height=0.62, label="erradicadas (0)")
    for i, (c, e) in enumerate(zip(control, eradicates)):
        a1.text(c + 0.3, i, f"{c}/{n} · {e} a 0", va="center", fontsize=8.5, color=ink2)
    a1.set_xlim(0, n + 7)
    a1.set_title("Pacientes controladas y erradicadas (de 24)", loc="left", fontsize=11, color=ink, pad=10)
    a1.legend(loc="upper left", bbox_to_anchor=(0, -0.06), ncol=2, frameon=False, fontsize=9)
    a2.barh(y, healthy, color=series, height=0.62)
    for i, s in enumerate(healthy):
        a2.text(s + 4, i, f"{s:.0f}", va="center", fontsize=8.5, color=ink2)
    a2.set_xlim(0, max(healthy) * 1.15)
    a2.set_title("Células sanas destruidas (mediana)", loc="left", fontsize=11, color=ink, pad=10)
    a1.set_yticks(list(y), tags, fontsize=9, color=ink)
    for a in (a1, a2):
        a.set_facecolor("#fcfcfb")
        a.grid(axis="x", color=lattice, linewidth=0.8)
        a.set_axisbelow(True)
        for side in ("top", "right", "left"):
            a.spines[side].set_visible(False)
        a.spines["bottom"].set_color(lattice)
        a.tick_params(axis="x", colors=ink2, labelsize=8.5)
        a.tick_params(axis="y", length=0)
    fig.suptitle(
        "¿Se puede tratar sin operar? · luminal, tratamiento día 120, seguimiento "
        "500 días, 24 pacientes simuladas",
        x=0.01,
        ha="left",
        fontsize=12,
        color=ink,
    )
    fig.subplots_adjust(left=0.25, right=0.98, top=0.86, bottom=0.12, wspace=0.08)
    fig.savefig(OUTPUT / "04_no_surgery_results.png", dpi=110)
    plt.close(fig)
    print("  · 04_no_surgery_results.png")


def _state3d(sim, pos, c) -> int:
    """El mismo código de colores que el visor 2D."""
    from breast_cancer_abm.model.tissue3d import MYO, STROMA

    if c.kind == FIBROBLAST:
        return -1
    if c.kind == MYOEPITHELIAL:
        return -2
    if c.necrosis:
        return NECROSIS
    if c.n_drivers == 0:
        return NORMAL
    if c.stem:
        return QUEEN
    if sim.zone[pos] in (MYO, STROMA):
        return INVASIVE
    if c.anoikis_res:
        return IN_SITU
    return HYPERPLASIA if c.prolif_signal else INITIATED


def _cutaway3d(ax, sim, title: str):
    """Medio cilindro (x <= eje) para ver dentro: lo que una célula terapéutica ve."""
    axis_mid = sim.center
    groups: dict[int, list] = {}
    for pos, c in sim.cells.items():
        z, y, x = pos
        if x > axis_mid + 0.5:
            continue
        e = _state3d(sim, pos, c)
        if e == -1:
            continue
        groups.setdefault(e, []).append((x, y, z))
    styles = {  # (color, tamaño, opacidad, etiqueta)
        -2: ("#8e7cc3", 3, 0.10, "mioepitelial"),
        NORMAL: ("#f4b6c8", 3, 0.18, "normal"),
        INITIATED: (COLORS[INITIATED], 5, 0.6, "iniciada"),
        HYPERPLASIA: (COLORS[HYPERPLASIA], 5, 0.6, "hiperplasia"),
        IN_SITU: (COLORS[IN_SITU], 6, 0.55, "carcinoma in situ"),
        INVASIVE: (COLORS[INVASIVE], 7, 0.8, "carcinoma invasivo"),
        QUEEN: ("#ffffff", 10, 1.0, "célula madre tumoral"),
        NECROSIS: (COLORS[NECROSIS], 6, 0.7, "necrosis"),
    }
    for e, (color, size, opacity, name) in styles.items():
        if e in groups:
            xs, ys, zs = zip(*groups[e])
            ax.scatter(
                xs,
                ys,
                zs,
                s=size,
                c=color,
                alpha=opacity,
                linewidths=0,
                label=f"{name} ({len(groups[e])})",
                depthshade=False,
            )
    ai = [a.pos for a in sim.agents if a.kind == "AI" and a.pos[2] <= axis_mid + 0.5]
    if ai:
        zs, ys, xs = zip(*ai)
        ax.scatter(
            xs,
            ys,
            zs,
            s=26,
            c="#f72585",
            marker="*",
            linewidths=0,
            label=f"célula terapéutica ({len(ai)})",
            depthshade=False,
        )
    z_len, side, _ = sim.grid_shape
    ax.set_xlim(0, side), ax.set_ylim(0, side), ax.set_zlim(0, z_len)
    ax.set_box_aspect((side / 2, side, z_len))
    ax.view_init(elev=18, azim=-28)
    ax.set_axis_off()
    r = sim.summary()
    ax.set_title(
        f"{title}\ntumor {r['tumor_cells']:,} · peligrosas {r['dangerous']:,} · "
        f"accesible {sim.accessible_fraction():.0%}".replace(",", "."),
        fontsize=10,
    )
    ax.legend(loc="upper left", fontsize=7, frameon=False, markerscale=1.6)


def figures_3d():
    from breast_cancer_abm.model.tissue3d import Simulation3D

    def advance3d(sim, days):
        for _ in range(int(days * STEPS_PER_DAY)):
            sim.step()

    plt.style.use("dark_background")
    sim = Simulation3D("luminal", 1)
    fig = plt.figure(figsize=(18, 7.2))
    advance3d(sim, 45)
    _cutaway3d(fig.add_subplot(1, 3, 1, projection="3d"), sim, "Día 45 · lesión inicial")
    advance3d(sim, 35)
    _cutaway3d(fig.add_subplot(1, 3, 2, projection="3d"), sim, "Día 80 · la masa llena el conducto")
    sim.release_therapy("trojan+niche_guardian+memory", 300)
    advance3d(sim, 20)
    _cutaway3d(fig.add_subplot(1, 3, 3, projection="3d"), sim, "Día 100 · 20 días después del arsenal")
    fig.suptitle(
        "Modelo 3D (tejido3d.py): medio conducto mamario cortado a lo largo. "
        "«Accesible» = células tumorales con un hueco al lado, por donde puede "
        "entrar una terapia",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT / "05_3d_duct_cutaway.png", dpi=85)
    plt.close(fig)
    print("  · 05_3d_duct_cutaway.png")

    # Secciones transversales del mismo día 100: dónde están las terapéuticas.
    z_len, side, _ = sim.grid_shape
    import numpy as np
    from matplotlib.colors import to_rgb

    fig, panels = plt.subplots(1, 4, figsize=(16, 4.6))
    for ax, z in zip(panels, (z_len // 8, z_len // 3, z_len // 2, 3 * z_len // 4)):
        img = np.full((side, side, 3), to_rgb("#231b17"))
        img[sim.zone[z] == 0] = to_rgb("#120c16")
        img[sim.vessels[z]] = to_rgb("#ff2d55")
        for (zz, y, x), c in sim.cells.items():
            if zz != z:
                continue
            e = _state3d(sim, (zz, y, x), c)
            img[y, x] = to_rgb({-1: COLORS[1], -2: "#8e7cc3"}.get(e, COLORS.get(e, "#ffffff")))
        ax.imshow(img, interpolation="nearest")
        pts = [(a.pos[2], a.pos[1]) for a in sim.agents if a.kind == "AI" and a.pos[0] == z]
        if pts:
            ax.scatter(*zip(*pts), s=40, marker="*", c="#f72585")
        ax.set_title(f"sección z = {z} · terapéuticas en el plano: {len(pts)}", fontsize=10)
        ax.set_xticks([]), ax.set_yticks([])
    fig.suptitle(
        "Día 100 en 3D, secciones transversales del conducto (rojo: carcinoma · "
        "estrellas: células terapéuticas)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT / "06_3d_sections.png", dpi=85)
    plt.close(fig)
    print("  · 06_3d_sections.png")


def ablation_3d_figure():
    """07: lo que hace cada ablación en 3D el día del rescate por ctDNA (día 90)."""
    from breast_cancer_abm.model.biology import CALIBRATIONS
    from breast_cancer_abm.studies.no_surgery_3d import ARSENAL, _new_sim

    plt.style.use("dark_background")
    sim = _new_sim(4)
    for step in range(int(90 * STEPS_PER_DAY)):
        if step == int(60 * STEPS_PER_DAY):
            sim.p.mutation_rate = CALIBRATIONS["realistic"]
            sim.release_therapy(ARSENAL, 300)
        sim.step()
    disc, guided = copy.deepcopy(sim), copy.deepcopy(sim)
    rd = disc.ablate(guided=False, rim=2)
    rg = guided.ablate(guided=True, rim=1)
    fig = plt.figure(figsize=(18, 7.2))
    _cutaway3d(fig.add_subplot(1, 3, 1, projection="3d"), sim, "Día 90 · ctDNA positivo (arsenal el día 60)")
    _cutaway3d(
        fig.add_subplot(1, 3, 2, projection="3d"),
        disc,
        f"Ablación en disco · {rd['healthy']:,} sanas destruidas".replace(",", "."),
    )
    _cutaway3d(
        fig.add_subplot(1, 3, 3, projection="3d"),
        guided,
        f"Ablación guiada · {rg['healthy']:,} sanas destruidas".replace(",", "."),
    )
    fig.suptitle(
        "Modelo 3D, paciente 4: la misma paciente copiada tres veces. El disco "
        "apunta al centro de todo el tumor y deja focos fuera; la guiada va a "
        "cada foco visible",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT / "07_3d_disc_vs_guided_ablation.png", dpi=85)
    plt.close(fig)
    print("  · 07_3d_disc_vs_guided_ablation.png")


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=("2d", "3d", "chart", "ablation3d"), default=None)
    a = ap.parse_args(argv)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    print("Capturas en", OUTPUT)
    if a.only in (None, "2d"):
        viewer_figures()
    if a.only in (None, "chart"):
        results_chart()
    if a.only in (None, "3d"):
        figures_3d()
    if a.only in (None, "ablation3d"):
        ablation_3d_figure()


if __name__ == "__main__":
    main()
