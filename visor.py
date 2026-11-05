"""Visor interactivo: Play/Pausa, paso a paso, terapias y clic para inspeccionar células."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.widgets import Button, RadioButtons, Slider

from biologia import PASOS_POR_DIA
from politicas import TERAPIAS
from tejido import (ESTROMA, FIBRO, HIPERPLASIA, IN_SITU, INICIADA, INVASIVA, LUZ, MIOEP,
                    NECROSIS, NOMBRE_ESTADO, NORMAL, REINA, REPARADA, Simulacion)

COLORES = {
    FIBRO: "#6b5a4a", MIOEP: "#8e7cc3", NORMAL: "#f4b6c8", REPARADA: "#3ddc84",
    INICIADA: "#ffe066", HIPERPLASIA: "#ff9f1c", IN_SITU: "#e63946", INVASIVA: "#9d0208",
    NECROSIS: "#4a4a4a", REINA: "#ffffff",
}
FONDO_ZONA = {LUZ: "#120c16", ESTROMA: "#231b17"}
COLOR_VASO = "#ff2d55"
PALETA = np.zeros((11, 3))
for codigo, color in COLORES.items():
    PALETA[codigo] = to_rgb(color)
CURVAS = (INICIADA, HIPERPLASIA, IN_SITU, INVASIVA, REINA, NECROSIS, REPARADA)
MARCAS_AGENTE = {"T": ("^", "#4cc9f0"), "NK": ("D", "#f9f871"), "IA": ("*", "#f72585")}


class Visor:
    def __init__(self, escenario: str, semilla: int | None, terapia: str | None, dia_terapia: float | None):
        self.escenario, self.semilla = escenario, semilla
        self.terapia_auto, self.dia_terapia = terapia, dia_terapia
        self.jugando = False
        self.mostrar_o2 = False
        self.seleccion = None

        plt.style.use("dark_background")
        self.fig = plt.figure(figsize=(15, 8.5))
        self.fig.canvas.manager.set_window_title("Simulador de cáncer de mama")
        self.ax = self.fig.add_axes([0.01, 0.12, 0.50, 0.84])
        self.ax_curvas = self.fig.add_axes([0.57, 0.62, 0.41, 0.33])
        self.ax_info = self.fig.add_axes([0.53, 0.19, 0.46, 0.38])
        self.ax_info.axis("off")
        self._widgets()
        self.fig.canvas.mpl_connect("button_press_event", self._clic)
        self._reiniciar()
        self.anim = FuncAnimation(self.fig, self._fotograma, interval=60, cache_frame_data=False)

    # ------------------------------------------------------------------ montaje
    def _widgets(self):
        f = self.fig
        self.b_play = Button(f.add_axes([0.03, 0.03, 0.09, 0.05]), "Play")
        self.b_paso = Button(f.add_axes([0.13, 0.03, 0.09, 0.05]), "Paso (6 h)")
        self.b_reset = Button(f.add_axes([0.23, 0.03, 0.09, 0.05]), "Reiniciar")
        self.b_o2 = Button(f.add_axes([0.33, 0.03, 0.09, 0.05]), "Oxígeno")
        self.b_play.on_clicked(lambda _: self._alternar())
        self.b_paso.on_clicked(lambda _: self._avanzar(1))
        self.b_reset.on_clicked(lambda _: self._reiniciar())
        self.b_o2.on_clicked(lambda _: self._alternar_o2())
        self.b_vasos = Button(f.add_axes([0.43, 0.03, 0.08, 0.05]), "Vasos")
        self.b_vasos.on_clicked(lambda _: self._siguiente_modo_vasos())

        nombres = list(TERAPIAS)
        self.radio = RadioButtons(f.add_axes([0.55, 0.01, 0.13, 0.17]), nombres,
                                  active=nombres.index(self.terapia_auto.split("+")[0]) if self.terapia_auto else 0)
        for etiqueta in self.radio.labels:
            etiqueta.set_fontsize(8)
        self.b_liberar = Button(f.add_axes([0.70, 0.03, 0.12, 0.05]), "Liberar terapia",
                                color="#5a189a", hovercolor="#7b2cbf")
        self.b_liberar.on_clicked(lambda _: self._liberar(self.radio.value_selected))
        self.velocidad = Slider(f.add_axes([0.70, 0.11, 0.25, 0.025]), "Pasos/fotograma", 1, 40,
                                valinit=8, valstep=1)

    def _reiniciar(self):
        self.sim = Simulacion(self.escenario, self.semilla)
        if hasattr(self, "b_vasos"):
            self.b_vasos.label.set_text("Vasos")
        self.terapia_liberada = False
        self.seleccion = None
        n = self.sim.p.tamano
        self.fondo = np.zeros((n, n, 3))
        self.fondo[:] = to_rgb(FONDO_ZONA[ESTROMA])
        self.fondo[self.sim.zona <= 1] = to_rgb(FONDO_ZONA[LUZ])

        self.ax.clear()
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.img = self.ax.imshow(self._rgb(), interpolation="nearest")
        self.img_o2 = self.ax.imshow(self.sim.O2, cmap="inferno", vmin=0, vmax=1, alpha=0.75,
                                     visible=self.mostrar_o2)
        self.puntos = {t: self.ax.scatter([], [], marker=m, s=22, c=c, edgecolors="black",
                                          linewidths=0.4, zorder=3)
                       for t, (m, c) in MARCAS_AGENTE.items()}
        self.marca_sel = self.ax.scatter([], [], s=140, facecolors="none", edgecolors="white", zorder=4)
        leyenda = [Patch(color=COLORES[e], label=NOMBRE_ESTADO[e]) for e in
                   (NORMAL, INICIADA, HIPERPLASIA, IN_SITU, INVASIVA, REINA, NECROSIS, REPARADA, MIOEP, FIBRO)]
        leyenda.append(Patch(color=COLOR_VASO, label="vaso sanguíneo"))
        leyenda += [plt.Line2D([], [], marker=m, color=c, ls="", label=t)
                    for t, (m, c) in (("linfocito T", MARCAS_AGENTE["T"]), ("célula NK", MARCAS_AGENTE["NK"]),
                                      ("célula IA", MARCAS_AGENTE["IA"]))]
        self.ax.legend(handles=leyenda, loc="upper left", fontsize=7, framealpha=0.6, ncol=2)

        ac = self.ax_curvas
        ac.clear()
        ac.set_title("Evolución (escala logarítmica)", fontsize=9)
        ac.set_xlabel("días", fontsize=8)
        ac.set_yscale("symlog", linthresh=10)
        ac.tick_params(labelsize=7)
        ac.grid(alpha=0.15)
        self.lineas = {e: ac.plot([], [], color=COLORES[e], lw=1.6, label=NOMBRE_ESTADO[e])[0] for e in CURVAS}
        self.lineas["inmunes"] = ac.plot([], [], color=MARCAS_AGENTE["T"][1], lw=1, ls="--",
                                         label="inmunes naturales")[0]
        self.lineas["IA"] = ac.plot([], [], color=MARCAS_AGENTE["IA"][1], lw=1.4, ls="--",
                                    label="células IA")[0]
        ac.legend(fontsize=7, loc="upper left", ncol=2, framealpha=0.5)
        self.texto = self.ax_info.text(0, 1, "", va="top", ha="left", family="monospace", fontsize=8)
        self._dibujar()

    def _rgb(self):
        rgb = self.fondo.copy()
        codigo = self.sim.codigo
        m = codigo > 0
        rgb[m] = PALETA[codigo[m]]
        rgb[self.sim.vasos] = to_rgb(COLOR_VASO)
        return rgb

    # ------------------------------------------------------------------ control
    def _alternar(self):
        self.jugando = not self.jugando
        self.b_play.label.set_text("Pausa" if self.jugando else "Play")

    def _alternar_o2(self):
        self.mostrar_o2 = not self.mostrar_o2
        self.img_o2.set_visible(self.mostrar_o2)
        self._dibujar()

    def _siguiente_modo_vasos(self):
        modos = ("normal", "reconducir", "normalizar", "antiangiogenico", "asfixiar")
        self.sim.tratar_vasos(modos[(modos.index(self.sim.modo_vasos) + 1) % len(modos)])
        self.b_vasos.label.set_text(self.sim.modo_vasos[:10])
        self._dibujar()

    def _liberar(self, nombre: str):
        self.sim.liberar_terapia(nombre)
        self.terapia_liberada = True
        self._dibujar()

    def _avanzar(self, pasos: int):
        for _ in range(pasos):
            if self.terapia_auto and not self.terapia_liberada and self.sim.dia >= self.dia_terapia:
                self._liberar(self.terapia_auto)
            self.sim.paso()
        self._dibujar()

    def _fotograma(self, _):
        if self.jugando:
            self._avanzar(int(self.velocidad.val))

    def _clic(self, ev):
        if ev.inaxes is self.ax and ev.xdata is not None:
            self.seleccion = (int(round(ev.ydata)), int(round(ev.xdata)))
            self._dibujar()

    # ------------------------------------------------------------------ dibujo
    def _dibujar(self):
        s = self.sim
        self.img.set_data(self._rgb())
        if self.mostrar_o2:
            self.img_o2.set_data(s.O2)
        for tipo, sc in self.puntos.items():
            pos = [(x, y) for a in s.agentes if a.tipo == tipo for y, x in (a.pos,)]
            sc.set_offsets(pos if pos else np.empty((0, 2)))
        self.marca_sel.set_offsets([self.seleccion[::-1]] if self.seleccion else np.empty((0, 2)))
        horas = (s.t % PASOS_POR_DIA) * 24 // PASOS_POR_DIA
        self.ax.set_title(f"Conducto mamario · día {int(s.dia)} · {horas:02d}:00 h · escenario {s.escenario}",
                          fontsize=10)

        h = s.historia
        for clave, linea in self.lineas.items():
            if clave == "inmunes":
                datos = [a + b for a, b in zip(h["T"], h["NK"])]
            else:
                datos = h[clave]
            linea.set_data(h["dia"], datos)
        self.ax_curvas.relim()
        self.ax_curvas.autoscale_view()

        self.texto.set_text(self._informe())
        self.fig.canvas.draw_idle()

    def _informe(self) -> str:
        s, r, c = self.sim, self.sim.resumen(), self.sim.conteo
        n_tipo = {t: sum(1 for a in s.agentes if a.tipo == t) for t in ("T", "NK", "IA")}
        lineas = [
            f"Tumor: {r['tumorales']:>5}   iniciadas {c[INICIADA]} · hiperplasia {c[HIPERPLASIA]} · "
            f"in situ {c[IN_SITU]} · invasivas {c[INVASIVA]}",
            f"Normales {r['normales']} · reparadas {r['reparadas']} · necrosis {c[NECROSIS]} · "
            f"vasos nuevos {r['vasos_nuevos']} · roturas de membrana {r['roturas_membrana']}",
            f"Reinas (madres tumorales) {r['reinas']} · obreras que se volvieron reina "
            f"{r['desdiferenciaciones']} · guardianas {r['guardianas']} · se han comido "
            f"{r['fagocitadas']} (sanas {r['fagocitadas_sanas']}) · vasos {s.modo_vasos}",
            f"Inmunes: T {n_tipo['T']} · NK {n_tipo['NK']} · IA {n_tipo['IA']}   "
            f"muertes por inmunidad {r['muertes_inmunes']} · por terapia {r['muertes_terapia']}",
            f"Daño colateral (sanas) {r['colateral']} · reparaciones {r['reparaciones']} "
            f"(inútiles {r['reparaciones_inutiles']}) · apoptosis p53 {r['apoptosis_p53']}",
            f"  de las reparaciones: curan {r['reparaciones_que_curan']} · "
            f"matan vía p53 {r['reparaciones_p53_letales']} · "
            f"estériles {r['reparaciones_esteriles']}",
        ]
        frec = s.frecuencias_tumor()
        if frec:
            top = sorted(((v, g) for g, v in frec.items() if v > 0), reverse=True)
            lineas.append("Genes alterados en el tumor: " + " · ".join(f"{g} {v:.0%}" for v, g in top))
        lineas.append("")
        lineas.append("Eventos:")
        lineas += [e[:95] for e in list(s.eventos)[-8:]]
        if self.seleccion:
            lineas.append("")
            lineas += s.describir(self.seleccion)
        else:
            lineas += ["", "(haz clic en una célula para ver su genoma)"]
        return "\n".join(lineas)

    def mostrar(self):
        plt.show()
