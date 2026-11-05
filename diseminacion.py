"""Diseminación: las células que se escapan antes de operar.

    python diseminacion.py                 # valida la fase 1 (D1-D4)
    python diseminacion.py --semillas 10
    python diseminacion.py --d4            # D4: varía los SUPUESTO (32 semillas)

QUÉ AÑADE AL MODELO. El simulador base es un conducto mamario: un compartimento.
Pero lo que mata al sujeto no es el tumor que se extirpa, es la célula que ya se
había ido — y está en otro órgano, no se divide, y nadie la ve. Este módulo añade
los dos compartimentos que faltan:

    mama (tejido.py)  ──intravasación──>  circulación  ──siembra──>  médula ósea
                                          (casi todas              (dormancia,
                                           mueren)                  años)

NADA DE ESTO SE LE PROGRAMA AL MODELO como objetivo. Se deja que el tumor crezca,
invada y siembre por su cuenta, y después se mide contra los criterios de
`DISEMINACION.md`. Si no los pasa, la fase 1 no está lista.

Mecanismo de la dormancia, que NO es inventado (Ghajar, Nat Cell Biol 2013):
  · el microvaso ESTABLE mantiene dormida a la célula (trombospondina-1)
  · el BROTE vascular la despierta (TGF-β1, periostina)
Desde la v2 (2026-10-28) esto es literal en el código: la médula es un conjunto de
microvasos que brotan y se estabilizan; la DTC se aloja junto a uno y su estado
(dormida o no) es el de su vaso, y cambia cuando el vaso cambia. Ya no hay una
probabilidad de «entrar dormida». Diseño y supuestos: DISEMINACION.md, «Fase 1 · v2».
Y lo que NO hay que modelar como protección: la quiescencia no protege del ataque
de células T (93 % de eliminación, Goddard, Cancer Cell 2024). Lo que protege a la
DTC es la ESCASEZ RELATIVA — son poquísimas y raramente se encuentran con un
linfocito que las reconozca. Es una propiedad espacial, no bioquímica.
"""
from __future__ import annotations

import argparse
import math
import random
import statistics
import sys
from dataclasses import dataclass, field

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from biologia import PASOS_POR_DIA, Parametros
from tejido import Simulacion, LUMINAL


# --------------------------------------------------------------------- parámetros
@dataclass
class ParamsDiseminacion:
    """Todo lo que aquí no tiene fuente va marcado como SUPUESTO y entra en D4."""
    # DOS VÍAS DE ESCAPE, no una. El modelo original solo tenía la segunda, y por eso
    # D1 daba un interruptor (0 % con tumor pequeño, 100 % con tumor grande) en vez del
    # 30-40 % real: suponía, sin decirlo, que solo disemina el tumor invasivo.
    #
    # (1) DISEMINACIÓN TEMPRANA. Células que escapan desde lesiones AÚN NO INVASIVAS.
    #     Demostrado en modelos murinos (Hosseini et al. y Harper et al., Nature 2016):
    #     cambió cómo se entiende la cascada metastásica. Explica que la fracción de
    #     pacientes con DTC sea intermedia y relativamente independiente del tamaño:
    #     unas siembran pronto y otras no, en vez de sembrar todas al crecer.
    prob_diseminacion_temprana: float = 5e-6      # CALIBRADO (ver nota al final)
    competencia_temprana: float = 0.25    # SUPUESTO: fracción de tumores capaces de ello
    # (2) Intravasación clásica: una célula invasiva que toca un vaso se suelta.
    prob_intravasacion: float = 3e-4      # CALIBRADO (ver nota al final)
    # Tránsito: la cascada metastásica es MUY ineficiente (criterio D3).
    prob_sobrevive_transito: float = 0.02  # SUPUESTO dentro de un rango muy aceptado
    prob_muere_en_nicho: float = 0.0015    # SUPUESTO, por paso: el nicho no es gratis
    # NICHO MEDULAR (v2). La dormancia ya NO es un parámetro (`prob_entra_dormida` y
    # `prob_brote_vascular` eliminados el 2026-10-28): sale del vaso junto al que se
    # aloja la DTC. Vaso estable -> dormida (trombospondina-1); vaso en brote ->
    # despierta (TGF-β1, periostina). Ghajar et al., Nat Cell Biol 2013.
    # MEDIDO: volumen de vaso tipo H (angiogénico) / volumen vascular total en fémur
    # humano, 18,9 % a 45-59 años y 13,9 % a 60-75 (Wang et al., Cell Death Dis 2017).
    fraccion_vasos_brote: float = 0.16
    duracion_brote_dias: float = 14.0      # SUPUESTO: vida media de un brote
    # Peso del vaso estable al alojarse (1 = proporcional al volumen vascular). Price
    # 2016 (nicho perisinusoidal E-selectina/SDF-1) sugiere > 1; magnitud SUPUESTA.
    sesgo_homing_estable: float = 1.0
    n_vasos: int = 500                     # SUPUESTO de discretización
    # La DTC vuelve a dormirse si su vaso se estabiliza (Sosa et al., Nat Rev Cancer
    # 2014: dormancia reversible dictada por el entorno). Respaldo cualitativo.
    dormancia_reversible: bool = True
    # Detección: sensibilidad imperfecta, y esto SÍ está medido.
    sensibilidad_ihc: float = 0.24         # biopsia de médula por IHC
    sensibilidad_citometria: float = 0.47  # la misma población, otro método


@dataclass
class DTC:
    """Una célula tumoral diseminada, viviendo en la médula."""
    id: int
    n_drivers: int
    madre: bool
    vaso: int = 0
    dormida: bool = True
    dia_llegada: float = 0.0
    dia_despertar: float | None = None


@dataclass
class Organismo:
    """El sujeto completo: la mama más los compartimentos que faltaban."""
    sim: Simulacion
    p: ParamsDiseminacion = field(default_factory=ParamsDiseminacion)
    medula: list[DTC] = field(default_factory=list)
    en_transito: int = 0
    stats: dict = field(default_factory=lambda: dict(
        intravasaciones=0, muertas_en_transito=0, sembradas=0,
        despertadas=0, redormidas=0, muertas_en_nicho=0, via_temprana=0,
        via_clasica=0))
    # La médula tiene su propio generador: variar el nicho no cambia tumor ni siembra.
    rng_medula: random.Random | None = None
    vasos_brote: list[bool] = field(default_factory=list)
    _siguiente_id: int = 0
    _competente: bool | None = None

    def __post_init__(self):
        if self.rng_medula is None:
            self.rng_medula = random.Random(0)
        f = self.p.fraccion_vasos_brote
        # La médula adulta ya está en equilibrio: estado inicial estacionario.
        self.vasos_brote = [self.rng_medula.random() < f for _ in range(self.p.n_vasos)]
        dt = 1 / PASOS_POR_DIA
        k_menos = 1 / self.p.duracion_brote_dias
        k_mas = f * k_menos / (1 - f) if f < 1 else float("inf")
        self._p_brota = 1 - math.exp(-k_mas * dt)
        self._p_estabiliza = 1 - math.exp(-k_menos * dt)

    # ---------------------------------------------------------------- un paso
    def paso(self):
        self.sim.paso()
        self._intravasar()
        self._nicho()

    def _competente_temprana(self) -> bool:
        """¿Este sujeto tiene un tumor capaz de diseminar antes de ser invasivo?

        Se decide una vez por sujeto, no por célula: es una propiedad del clon, y es
        lo que produce una fracción intermedia de pacientes con DTC en vez de un
        todo-o-nada ligado al tamaño."""
        if self._competente is None:
            self._competente = self.sim.rng.random() < self.p.competencia_temprana
        return self._competente

    def _intravasar(self):
        """Dos vías: la temprana (no invasiva) y la clásica (invasiva junto a vaso)."""
        sim, p = self.sim, self.p
        # (1) VÍA TEMPRANA: no necesita invasión ni contacto con vaso
        if p.prob_diseminacion_temprana and self._competente_temprana():
            for pos, c in list(sim.celulas.items()):
                if c.tipo != LUMINAL or not c.n_drivers or c.necrosis:
                    continue
                if sim.rng.random() >= p.prob_diseminacion_temprana:
                    continue
                self.stats["intravasaciones"] += 1
                self.stats["via_temprana"] += 1
                if sim.rng.random() >= p.prob_sobrevive_transito:
                    self.stats["muertas_en_transito"] += 1
                    continue
                self._siembra(c)
        # (2) VÍA CLÁSICA
        for pos, c in list(sim.celulas.items()):
            if c.tipo != LUMINAL or not c.n_drivers or c.necrosis or not c.invasiva:
                continue
            if not any(sim.vasos[q] for q in sim._vecinos(pos)):
                continue
            if sim.rng.random() >= p.prob_intravasacion:
                continue
            self.stats["intravasaciones"] += 1
            # Tránsito: la inmensa mayoría no sobrevive (criterio D3)
            if sim.rng.random() >= p.prob_sobrevive_transito:
                self.stats["muertas_en_transito"] += 1
                continue
            self.stats["via_clasica"] += 1
            self._siembra(c)

    def _siembra(self, c):
        """La célula sobrevivió al tránsito y se aloja junto a un microvaso.

        Su estado NO se sortea: es el del vaso donde cae (estable -> dormida)."""
        p, rng = self.p, self.rng_medula
        pesos = [1.0 if b else p.sesgo_homing_estable for b in self.vasos_brote]
        v = rng.choices(range(len(pesos)), weights=pesos)[0]
        self._siguiente_id += 1
        brote = self.vasos_brote[v]
        self.medula.append(DTC(
            id=self._siguiente_id, n_drivers=c.n_drivers, madre=c.madre, vaso=v,
            dormida=not brote, dia_llegada=self.sim.dia,
            dia_despertar=self.sim.dia if brote else None))
        self.stats["sembradas"] += 1

    def _nicho(self):
        """La médula: los vasos brotan y se estabilizan, y las DTC siguen a su vaso."""
        p, rng = self.p, self.rng_medula
        vb = self.vasos_brote
        for i, b in enumerate(vb):
            if rng.random() < (self._p_estabiliza if b else self._p_brota):
                vb[i] = not b
        vivas = []
        for d in self.medula:
            if rng.random() < p.prob_muere_en_nicho:
                self.stats["muertas_en_nicho"] += 1
                continue
            brote = vb[d.vaso]
            if d.dormida and brote:
                # El brote vascular rompe la dormancia (TGF-β1, periostina)
                d.dormida = False
                d.dia_despertar = self.sim.dia
                self.stats["despertadas"] += 1
            elif not d.dormida and not brote and p.dormancia_reversible:
                # El vaso se estabiliza: vuelve la trombospondina-1
                d.dormida = True
                self.stats["redormidas"] += 1
            vivas.append(d)
        self.medula = vivas

    # ------------------------------------------------------------- observación
    def biopsia(self, metodo: str = "ihc") -> int:
        """Mirar la médula con sensibilidad IMPERFECTA, que es lo que pasa de verdad.

        Una terapia perfecta sobre una población que no sabes detectar no sirve."""
        s = (self.p.sensibilidad_ihc if metodo == "ihc"
             else self.p.sensibilidad_citometria)
        return sum(1 for _ in self.medula if self.sim.rng.random() < s)

    def resumen(self) -> dict:
        dormidas = sum(1 for d in self.medula if d.dormida)
        return dict(dtc=len(self.medula), dormidas=dormidas,
                    despiertas=len(self.medula) - dormidas,
                    frac_dormidas=dormidas / len(self.medula) if self.medula else 0.0,
                    tumor=sum(1 for c in self.sim.celulas.values()
                              if c.tipo == LUMINAL and c.n_drivers and not c.necrosis),
                    **self.stats)


# ------------------------------------------------------------------ validación
def correr(semilla: int, dias: float, p: ParamsDiseminacion | None = None) -> dict:
    sim = Simulacion("luminal", semilla, Parametros(prob_neoantigeno=0.1))
    org = Organismo(sim, p or ParamsDiseminacion(),
                    rng_medula=random.Random(1_000_003 + semilla))
    for _ in range(int(dias * PASOS_POR_DIA)):
        org.paso()
    r = org.resumen()
    r["biopsia_ihc"] = org.biopsia("ihc")
    r["biopsia_citometria"] = org.biopsia("citometria")
    return r


def _correr_tarea(args):
    return correr(*args)


def correr_lote(semillas, dias, p=None, procesos=None) -> list[dict]:
    """Las mismas corridas que `correr`, repartidas por semillas entre procesos."""
    tareas = [(s, dias, p) for s in semillas]
    if procesos == 1 or len(tareas) == 1:
        return [_correr_tarea(t) for t in tareas]
    import multiprocessing as mp
    with mp.Pool(procesos or max(1, (mp.cpu_count() or 2) - 4)) as pool:
        return pool.map(_correr_tarea, tareas)


def _d2_agregada(res: list[dict]) -> float:
    tot = sum(r["dtc"] for r in res)
    return sum(r["dormidas"] for r in res) / tot if tot else 0.0


def ee_d2_bootstrap(res: list[dict], n: int = 1000) -> float:
    """Error estándar de la D2 agregada remuestreando PACIENTES (generador fijo)."""
    rng = random.Random(0)
    vals = [_d2_agregada(rng.choices(res, k=len(res))) for _ in range(n)]
    return statistics.stdev(vals)


def metricas(res: list[dict]) -> dict:
    """D1, D2 y D3 tal como los define DISEMINACION.md (v2).

    D2 es la fracción AGREGADA (dormidas de todos / DTC de todos), anunciada antes de
    correr la v2; la mediana por paciente (v1) se devuelve como `d2_mediana`, sin
    veredicto: con 1-2 DTC por paciente solo vale 0, ½ o 1 y se satura."""
    con_dtc = sum(r["dtc"] > 0 for r in res)
    fracs = [r["frac_dormidas"] for r in res if r["dtc"] > 0]
    tot_intra = sum(r["intravasaciones"] for r in res)
    tot_muertas = sum(r["muertas_en_transito"] for r in res)
    d1 = con_dtc / len(res)
    d2 = _d2_agregada(res)
    d3 = tot_muertas / tot_intra if tot_intra else 0.0
    return dict(d1=d1, d2=d2, d3=d3, con_dtc=con_dtc, n=len(res), intra=tot_intra,
                dtc=sum(r["dtc"] for r in res),
                d2_mediana=statistics.median(fracs) if fracs else 0.0,
                ok1=0.25 <= d1 <= 0.50, ok2=d2 > 0.80, ok3=d3 > 0.90)


# Supuestos que entran en D4, y cómo se varía cada uno (fijado antes de correr).
SUPUESTOS_D4 = ("competencia_temprana", "prob_sobrevive_transito",
                "prob_muere_en_nicho", "duracion_brote_dias",
                "sesgo_homing_estable", "n_vasos", "fraccion_vasos_brote")
# Los parámetros MEDIDOS no se varían ×0,5/×2 sino a los extremos de su rango medido.
RANGOS_MEDIDOS = {"fraccion_vasos_brote": [("60-75 años", 0.139), ("45-59 años", 0.189)]}


def _variantes(nombre: str, valor: float) -> list[tuple[str, float]]:
    if nombre in RANGOS_MEDIDOS:
        return RANGOS_MEDIDOS[nombre]
    if isinstance(valor, int):
        return [("×0,5", max(1, round(valor * 0.5))), ("×2", valor * 2)]
    tope = 1.0 if nombre.startswith("prob_") or nombre == "competencia_temprana" else None
    doble = valor * 2 if tope is None else min(tope, valor * 2)
    return [("×0,5", valor * 0.5), ("×2", doble)]


def control_d4(semillas, dias, procesos=None):
    """D4: ¿depende el resultado de un parámetro inventado sin respaldo?

    CRITERIO DE ACEPTACIÓN v2 (escrito antes de ejecutar; copia en DISEMINACION.md,
    apartado «D4 v2»):
      · Cada supuesto de SUPUESTOS_D4 se corre a ×0,5 y ×2 (probabilidades recortadas
        a [0, 1]); los MEDIDOS (fraccion_vasos_brote) a los extremos de su rango
        medido (0,139 y 0,189).
      · Mismas semillas en la base y en todas las variantes.
      · Para cada supuesto y cada Dk (D1-D3) el veredicto PASA/FALLA «se mantiene» si
        es igual que en la base en ambas variantes; «cambia» si no, y la diferencia
        supera el ruido (D1: 2·√(2·p(1−p)/n); D2 agregada: 2·√2·EE bootstrap por
        pacientes en la base; D3 tal cual); «no concluyente» si cambia dentro del ruido.
      · D4 PASA si ningún supuesto cambia ningún veredicto; FALLA si alguno lo cambia
        por encima del ruido; NO CONCLUYENTE si solo hay cambios dentro del ruido.
      · Una dependencia por construcción que cambie el veredicto cuenta como FALLO.
      · Control estructural INFORMATIVO (no cuenta para D4): despertar irreversible.
    """
    base_p = ParamsDiseminacion()
    n = len(semillas)
    print(f"D4 v2 · control de dependencia · {n} semillas · {dias:g} días\n")
    res_base = correr_lote(semillas, dias, base_p, procesos)
    base = metricas(res_base)

    def fila(etq, m):
        marca = lambda ok: "ok" if ok else "XX"
        print(f"  {etq:42s} D1 {m['d1']:5.0%} [{marca(m['ok1'])}]  "
              f"D2 {m['d2']:5.1%} [{marca(m['ok2'])}]  "
              f"D3 {m['d3']:6.1%} [{marca(m['ok3'])}]  "
              f"(DTC {m['dtc']}, mediana D2 {m['d2_mediana']:.0%})", flush=True)

    fila("BASE", base)
    ruido_d1 = 2 * (2 * max(base["d1"] * (1 - base["d1"]), 1 / n) / n) ** 0.5
    ruido_d2 = 2 * 2 ** 0.5 * ee_d2_bootstrap(res_base)
    ruido = {"d1": ruido_d1, "d2": ruido_d2}
    peor = "PASA"
    for nombre in SUPUESTOS_D4:
        v0 = getattr(base_p, nombre)
        print(f"\n  {nombre} (base {v0:g})")
        estados = {}
        for etq, v in _variantes(nombre, v0):
            m = metricas(correr_lote(
                semillas, dias, ParamsDiseminacion(**{nombre: v}), procesos))
            fila(f"  {etq:10s} = {v:g}", m)
            for k, dk in (("ok1", "d1"), ("ok2", "d2"), ("ok3", "d3")):
                if m[k] == base[k]:
                    continue
                dentro = dk in ruido and abs(m[dk] - base[dk]) <= ruido[dk]
                est = "no concluyente" if dentro else "CAMBIA"
                if estados.get(dk) != "CAMBIA":
                    estados[dk] = est
        for dk in ("d1", "d2", "d3"):
            est = estados.get(dk, "se mantiene")
            if est != "se mantiene":
                print(f"    -> {dk.upper()}: veredicto {est}")
            if est == "CAMBIA":
                peor = "FALLA"
            elif est == "no concluyente" and peor == "PASA":
                peor = "NO CONCLUYENTE"
    print(f"\n  (ruido tolerado: D1 ±{ruido_d1:.0%}, D2 ±{ruido_d2:.1%})")
    print(f"\n  D4 {peor}")
    print("\n  Control estructural (INFORMATIVO, no cuenta para D4):")
    m = metricas(correr_lote(semillas, dias,
                             ParamsDiseminacion(dormancia_reversible=False), procesos))
    fila("  despertar irreversible", m)
    return peor


def validar(semillas, dias, procesos=None):
    """Criterios D1-D4 de DISEMINACION.md. Se pasan o no se pasan."""
    print(f"FASE 1 · {len(semillas)} semillas · {dias:g} días (momento del diagnóstico)\n")
    res = correr_lote(semillas, dias, None, procesos)
    for s, r in zip(semillas, res):
        print(f"  s{s:2d}  tumor {r['tumor']:6d} · DTC {r['dtc']:4d} "
              f"({r['frac_dormidas']:.0%} dormidas) · intravasan {r['intravasaciones']:5d} "
              f"· mueren en tránsito {r['muertas_en_transito']:5d} · "
              f"biopsia IHC {r['biopsia_ihc']:3d}", flush=True)

    m = metricas(res)
    con_dtc, d1, d2, d3, tot_intra = m["con_dtc"], m["d1"], m["d2"], m["d3"], m["intra"]

    print(f"\n{'='*70}")
    def linea(clave, texto, valor, ok, diana):
        print(f"  [{'PASA ' if ok else 'FALLA'}] {clave} {texto}")
        print(f"          medido {valor} · diana {diana}")
        return ok

    ok1 = linea("D1", "pacientes con DTC al diagnóstico",
                f"{d1:.0%} ({con_dtc}/{len(res)})", 0.25 <= d1 <= 0.50, "30-40 % (Braun 2005)")
    ok2 = linea("D2", "fracción de DTC quiescentes (agregada)",
                f"{d2:.1%} de {m['dtc']} DTC (mediana por paciente {m['d2_mediana']:.0%}, "
                f"sin veredicto)", d2 > 0.80, ">80 % (Pantel 2004)")
    ok3 = linea("D3", "mueren en el trayecto",
                f"{d3:.1%} de {tot_intra} intravasaciones", d3 > 0.90, "la gran mayoría")
    print("\n  [ -- ] D4 control de dependencia: ejecutar aparte "
          "(--d4) antes de dar la fase por cerrada")
    print(f"\n  Detección: IHC ve {statistics.median(r['biopsia_ihc'] for r in res):.0f} "
          f"y citometría {statistics.median(r['biopsia_citometria'] for r in res):.0f} "
          f"de {statistics.median(r['dtc'] for r in res):.0f} reales "
          f"(el 24 % y el 47 % medidos en clínica)")
    print(f"\n{'FASE 1 LISTA' if (ok1 and ok2 and ok3) else 'FASE 1 NO LISTA: no se construye encima'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--semillas", type=int, default=None,
                    help="por defecto 8 (validación) o 32 (--d4)")
    ap.add_argument("--dias", type=float, default=250,
                    help="momento del 'diagnóstico' en el que se mide")
    ap.add_argument("--d4", action="store_true",
                    help="control de dependencia de los parámetros SUPUESTO")
    ap.add_argument("--procesos", type=int, default=None,
                    help="procesos en paralelo (1 = en serie)")
    a = ap.parse_args()
    n = a.semillas or (32 if a.d4 else 8)
    if a.d4:
        control_d4(tuple(range(1, n + 1)), a.dias, a.procesos)
    else:
        validar(tuple(range(1, n + 1)), a.dias, a.procesos)


if __name__ == "__main__":
    main()


# --------------------------------------------------------------- NOTA DE CALIBRACIÓN
# [HISTÓRICO, v1] prob_brote_vascular ya no existe (2026-10-28): el despertar sale
# del recambio vascular del nicho. prob_intravasacion y prob_diseminacion_temprana
# siguen siendo los elegidos aquí; NO se han recalibrado para la v2.
#
# Medido el 2026-10-23, 10 semillas, día 250. Barrido de prob_intravasacion,
# prob_brote_vascular y prob_diseminacion_temprana contra D1 y D2:
#
#   intravasación   brote   temprana | D1 % con DTC | D2 % dormidas
#       4e-03       8e-04    2e-05   |     80 %     |     76 %
#       1e-03       8e-04    2e-05   |     70 %     |    100 %
#    >> 3e-04       8e-04    5e-06   |     40 %     |     83 %   <- ELEGIDA
#       1e-04       3e-04    2e-06   |     30 %     |    100 %
#       3e-05       3e-04    1e-06   |     10 %     |    100 %
#       1e-05       1e-04    5e-07   |      0 %     |      0 %
#
# POR QUÉ LA ELEGIDA Y NO LA DE 1e-04, que también pasa los dos criterios: aquella
# da D2 = 100 %, es decir NINGUNA célula despierta nunca. Un modelo así no puede
# producir recurrencias, y la fase 2 consiste exactamente en reproducir las recaídas
# repartidas entre los años 5 y 20. El ajuste elegido deja un 17 % de células
# despiertas, que es el material con el que la fase 2 tiene que trabajar.
#
# ESTO ES CALIBRACIÓN, NO VALIDACIÓN, Y LA DIFERENCIA IMPORTA. Un parámetro ajustado
# contra un dato no predice ese dato: lo reproduce por construcción. Que el modelo dé
# ahora un 40 % de pacientes con DTC NO es un acierto — es el ajuste funcionando.
# La validación, si llega, vendrá de lo que el modelo haga DESPUÉS con estos números
# fijos: la curva de recurrencias (fase 2) y la asimetría célula T / ADC contra células
# dormidas (fase 3). Compárese con T1 en validacion.py, donde la geometría del oxígeno
# ACERTÓ los 150-180 µm de Thomlinson & Gray sin que nadie la ajustara: eso sí es
# validación.
