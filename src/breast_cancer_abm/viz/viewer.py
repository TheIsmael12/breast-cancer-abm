"""Viewer interactivo: Play/Pausa, paso a paso, terapias y clic para inspeccionar células."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.widgets import Button, RadioButtons, Slider

from breast_cancer_abm.model.biology import STEPS_PER_DAY
from breast_cancer_abm.model.policies import THERAPIES
from breast_cancer_abm.model.tissue import (
    FIBRO,
    HYPERPLASIA,
    IN_SITU,
    INITIATED,
    INVASIVE,
    LUMEN,
    MYOEP,
    NECROSIS,
    NORMAL,
    QUEEN,
    REPAIRED,
    STATE_NAME,
    STROMA,
    Simulation,
)
from breast_cancer_abm.palette import STATE_COLORS, VESSEL_COLOR, ZONE_BACKGROUND

COLORS = STATE_COLORS
PALETTE = np.zeros((11, 3))
for code, color in COLORS.items():
    PALETTE[code] = to_rgb(color)
CURVES = (INITIATED, HYPERPLASIA, IN_SITU, INVASIVE, QUEEN, NECROSIS, REPAIRED)
AGENT_MARKERS = {"T": ("^", "#4cc9f0"), "NK": ("D", "#f9f871"), "AI": ("*", "#f72585")}


class Viewer:
    def __init__(self, scenario: str, seed: int | None, therapy: str | None, therapy_day: float | None):
        self.scenario, self.seed = scenario, seed
        self.auto_therapy, self.therapy_day = therapy, therapy_day
        self.playing = False
        self.show_o2 = False
        self.selection = None

        plt.style.use("dark_background")
        self.fig = plt.figure(figsize=(15, 8.5))
        self.fig.canvas.manager.set_window_title("Simulador de cáncer de mama")
        self.ax = self.fig.add_axes([0.01, 0.12, 0.50, 0.84])
        self.ax_curves = self.fig.add_axes([0.57, 0.62, 0.41, 0.33])
        self.ax_info = self.fig.add_axes([0.53, 0.19, 0.46, 0.38])
        self.ax_info.axis("off")
        self._widgets()
        self.fig.canvas.mpl_connect("button_press_event", self._click)
        self._reset()
        self.anim = FuncAnimation(self.fig, self._frame, interval=60, cache_frame_data=False)

    # ------------------------------------------------------------------ montaje
    def _widgets(self):
        f = self.fig
        self.b_play = Button(f.add_axes([0.03, 0.03, 0.09, 0.05]), "Play")
        self.b_step = Button(f.add_axes([0.13, 0.03, 0.09, 0.05]), "Paso (6 h)")
        self.b_reset = Button(f.add_axes([0.23, 0.03, 0.09, 0.05]), "Reiniciar")
        self.b_o2 = Button(f.add_axes([0.33, 0.03, 0.09, 0.05]), "Oxígeno")
        self.b_play.on_clicked(lambda _: self._toggle())
        self.b_step.on_clicked(lambda _: self._advance(1))
        self.b_reset.on_clicked(lambda _: self._reset())
        self.b_o2.on_clicked(lambda _: self._toggle_o2())
        self.b_vessels = Button(f.add_axes([0.43, 0.03, 0.08, 0.05]), "Vasos")
        self.b_vessels.on_clicked(lambda _: self._next_vessel_mode())

        names = list(THERAPIES)
        self.radius = RadioButtons(
            f.add_axes([0.55, 0.01, 0.13, 0.17]),
            names,
            active=names.index(self.auto_therapy.split("+")[0]) if self.auto_therapy else 0,
        )
        for tag in self.radius.labels:
            tag.set_fontsize(8)
        self.b_release = Button(
            f.add_axes([0.70, 0.03, 0.12, 0.05]), "Liberar terapia", color="#5a189a", hovercolor="#7b2cbf"
        )
        self.b_release.on_clicked(lambda _: self._release(self.radius.value_selected))
        self.speed = Slider(
            f.add_axes([0.70, 0.11, 0.25, 0.025]), "Pasos/fotograma", 1, 40, valinit=8, valstep=1
        )

    def _reset(self):
        self.sim = Simulation(self.scenario, self.seed)
        if hasattr(self, "b_vessels"):
            self.b_vessels.label.set_text("Vasos")
        self.therapy_released = False
        self.selection = None
        n = self.sim.p.grid_size
        self.background = np.zeros((n, n, 3))
        self.background[:] = to_rgb(ZONE_BACKGROUND[STROMA])
        self.background[self.sim.zone <= 1] = to_rgb(ZONE_BACKGROUND[LUMEN])

        self.ax.clear()
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.img = self.ax.imshow(self._rgb(), interpolation="nearest")
        self.img_o2 = self.ax.imshow(
            self.sim.O2, cmap="inferno", vmin=0, vmax=1, alpha=0.75, visible=self.show_o2
        )
        self.points = {
            t: self.ax.scatter([], [], marker=m, s=22, c=c, edgecolors="black", linewidths=0.4, zorder=3)
            for t, (m, c) in AGENT_MARKERS.items()
        }
        self.sel_mark = self.ax.scatter([], [], s=140, facecolors="none", edgecolors="white", zorder=4)
        legend_items = [
            Patch(color=COLORS[e], label=STATE_NAME[e])
            for e in (
                NORMAL,
                INITIATED,
                HYPERPLASIA,
                IN_SITU,
                INVASIVE,
                QUEEN,
                NECROSIS,
                REPAIRED,
                MYOEP,
                FIBRO,
            )
        ]
        legend_items.append(Patch(color=VESSEL_COLOR, label="vaso sanguíneo"))
        legend_items += [
            plt.Line2D([], [], marker=m, color=c, ls="", label=t)
            for t, (m, c) in (
                ("linfocito T", AGENT_MARKERS["T"]),
                ("célula NK", AGENT_MARKERS["NK"]),
                ("célula IA", AGENT_MARKERS["AI"]),
            )
        ]
        self.ax.legend(handles=legend_items, loc="upper left", fontsize=7, framealpha=0.6, ncol=2)

        ac = self.ax_curves
        ac.clear()
        ac.set_title("Evolución (escala logarítmica)", fontsize=9)
        ac.set_xlabel("días", fontsize=8)
        ac.set_yscale("symlog", linthresh=10)
        ac.tick_params(labelsize=7)
        ac.grid(alpha=0.15)
        self.lines = {e: ac.plot([], [], color=COLORS[e], lw=1.6, label=STATE_NAME[e])[0] for e in CURVES}
        self.lines["immune"] = ac.plot(
            [], [], color=AGENT_MARKERS["T"][1], lw=1, ls="--", label="inmunes naturales"
        )[0]
        self.lines["AI"] = ac.plot([], [], color=AGENT_MARKERS["AI"][1], lw=1.4, ls="--", label="células IA")[
            0
        ]
        ac.legend(fontsize=7, loc="upper left", ncol=2, framealpha=0.5)
        self.message = self.ax_info.text(0, 1, "", va="top", ha="left", family="monospace", fontsize=8)
        self._draw()

    def _rgb(self):
        rgb = self.background.copy()
        code = self.sim.code
        m = code > 0
        rgb[m] = PALETTE[code[m]]
        rgb[self.sim.vessels] = to_rgb(VESSEL_COLOR)
        return rgb

    # ------------------------------------------------------------------ control
    def _toggle(self):
        self.playing = not self.playing
        self.b_play.label.set_text("Pausa" if self.playing else "Play")

    def _toggle_o2(self):
        self.show_o2 = not self.show_o2
        self.img_o2.set_visible(self.show_o2)
        self._draw()

    def _next_vessel_mode(self):
        modes = ("normal", "reroute", "normalize", "antiangiogenic", "starve")
        self.sim.treat_vessels(modes[(modes.index(self.sim.vessel_mode) + 1) % len(modes)])
        self.b_vessels.label.set_text(self.sim.vessel_mode[:10])
        self._draw()

    def _release(self, name: str):
        self.sim.release_therapy(name)
        self.therapy_released = True
        self._draw()

    def _advance(self, steps: int):
        for _ in range(steps):
            if self.auto_therapy and not self.therapy_released and self.sim.day >= self.therapy_day:
                self._release(self.auto_therapy)
            self.sim.step()
        self._draw()

    def _frame(self, _):
        if self.playing:
            self._advance(int(self.speed.val))

    def _click(self, ev):
        if ev.inaxes is self.ax and ev.xdata is not None:
            self.selection = (int(round(ev.ydata)), int(round(ev.xdata)))
            self._draw()

    # ------------------------------------------------------------------ dibujo
    def _draw(self):
        s = self.sim
        self.img.set_data(self._rgb())
        if self.show_o2:
            self.img_o2.set_data(s.O2)
        for kind, sc in self.points.items():
            pos = [(x, y) for a in s.agents if a.kind == kind for y, x in (a.pos,)]
            sc.set_offsets(pos if pos else np.empty((0, 2)))
        self.sel_mark.set_offsets([self.selection[::-1]] if self.selection else np.empty((0, 2)))
        hours = (s.t % STEPS_PER_DAY) * 24 // STEPS_PER_DAY
        self.ax.set_title(
            f"Conducto mamario · día {int(s.day)} · {hours:02d}:00 h · escenario {s.scenario}", fontsize=10
        )

        h = s.history
        for key, line in self.lines.items():
            if key == "immune":
                dataset = [a + b for a, b in zip(h["T"], h["NK"])]
            else:
                dataset = h[key]
            line.set_data(h["day"], dataset)
        self.ax_curves.relim()
        self.ax_curves.autoscale_view()

        self.message.set_text(self._report())
        self.fig.canvas.draw_idle()

    def _report(self) -> str:
        s, r, c = self.sim, self.sim.summary(), self.sim.counts
        n_kind = {t: sum(1 for a in s.agents if a.kind == t) for t in ("T", "NK", "AI")}
        lines = [
            f"Tumor: {r['tumor_cells']:>5}   iniciadas {c[INITIATED]} · hiperplasia {c[HYPERPLASIA]} · "
            f"in situ {c[IN_SITU]} · invasivas {c[INVASIVE]}",
            f"Normales {r['normal_cells']} · reparadas {r['repaired_cells']} · necrosis {c[NECROSIS]} · "
            f"vasos nuevos {r['new_vessels']} · roturas de membrana {r['membrane_breaches']}",
            f"Reinas (madres tumorales) {r['queens']} · obreras que se volvieron reina "
            f"{r['dedifferentiations']} · guardianas {r['guardians']} · se han comido "
            f"{r['phagocytosed']} (sanas {r['phagocytosed_healthy']}) · vasos {s.vessel_mode}",
            f"Inmunes: T {n_kind['T']} · NK {n_kind['NK']} · IA {n_kind['AI']}   "
            f"muertes por inmunidad {r['immune_kills']} · por terapia {r['therapy_kills']}",
            f"Daño colateral (sanas) {r['collateral']} · reparaciones {r['repairs']} "
            f"(inútiles {r['repairs_useless']}) · apoptosis p53 {r['apoptosis_p53']}",
            f"  de las reparaciones: curan {r['repairs_that_cure']} · "
            f"matan vía p53 {r['repairs_p53_lethal']} · "
            f"estériles {r['repairs_sterile']}",
        ]
        freq = s.tumor_frequencies()
        if freq:
            top = sorted(((v, g) for g, v in freq.items() if v > 0), reverse=True)
            lines.append("Genes alterados en el tumor: " + " · ".join(f"{g} {v:.0%}" for v, g in top))
        lines.append("")
        lines.append("Eventos:")
        lines += [e[:95] for e in list(s.events)[-8:]]
        if self.selection:
            lines.append("")
            lines += s.describe(self.selection)
        else:
            lines += ["", "(haz clic en una célula para ver su genoma)"]
        return "\n".join(lines)

    def show(self):
        plt.show()
