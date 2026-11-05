"""Motor de la simulación: corte transversal de un conducto mamario.

Geometría (de dentro a fuera): luz del conducto · epitelio luminal (donde nace
casi todo el cáncer de mama) · capa mioepitelial + membrana basal (la barrera
que separa in situ de invasivo) · estroma con fibroblastos y vasos sanguíneos.
"""
from __future__ import annotations

import random
from collections import deque
from dataclasses import replace

import numpy as np

from agentes import (Agente, Analizar, Atacar, Mover, Observacion, Percepcion, Reparar,
                     politica_linfocito_t, politica_macrofago, politica_nk)
from biologia import (ESCENARIOS, GENES, GENES_INMUNES, PASOS_POR_DIA, Parametros,
                      PREVALENCIA_DIANA)
from politicas import TERAPIAS

LUZ, EPITELIO, MIO, ESTROMA = range(4)
NOMBRE_ZONA = ("luz", "epitelio", "mioepitelio", "estroma")

LUMINAL, MIOEPITELIAL, FIBROBLASTO = range(3)
NOMBRE_TEJIDO = ("luminal", "mioepitelial", "fibroblasto")

VACIO, FIBRO, MIOEP, NORMAL, REPARADA, INICIADA, HIPERPLASIA, IN_SITU, INVASIVA, NECROSIS, REINA = range(11)
NOMBRE_ESTADO = {
    VACIO: "vacío", FIBRO: "fibroblasto", MIOEP: "mioepitelial", NORMAL: "normal",
    REPARADA: "reparada", INICIADA: "iniciada", HIPERPLASIA: "hiperplasia/atipia",
    IN_SITU: "carcinoma in situ", INVASIVA: "carcinoma invasivo", NECROSIS: "necrosis",
    REINA: "célula madre tumoral (reina)",
}
TUMORALES = (INICIADA, HIPERPLASIA, IN_SITU, INVASIVA, REINA)
CONSUMO_O2 = np.array([0.0005, 0.002, 0.004, 0.004, 0.004, 0.005, 0.010, 0.016, 0.016, 0.0, 0.016])
FUENTE_QUIMIOCINA = np.array([0, 0, 0, 0, 0, 0, 0.3, 1.0, 1.0, 1.0, 1.0])

VECINOS = tuple((dy, dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx)
MAX_AGENTES = 500
VIDA_NECROSIS = 40


def _signo(v: float) -> int:
    return (v > 0) - (v < 0)


class Celula:
    __slots__ = ("id", "tipo", "genoma", "er", "mutaciones", "reparada", "necrosis",
                 "ultima_division", "alterados", "n_drivers", "senal_prolif", "p53_ok",
                 "hrd", "cdh1_perdido", "mhc", "pdl1", "her2", "anoikis_res", "invasiva",
                 "factor", "madre", "divisiones", "guardiana", "neo", "inmortal",
                 "telomerasa_ok", "senescente", "cebo", "aneuploide", "cd47", "caf",
                 "corrector", "pdl1_hasta")

    def __init__(self, id_: int, tipo: int, genoma: dict | None = None, er: bool = False,
                 mutaciones: int = 0, madre: bool = False):
        self.id = id_
        self.madre = madre        # célula madre: se divide sin límite (la "reina")
        self.divisiones = None    # divisiones que le quedan a una obrera tumoral (None = sin límite)
        self.guardiana = None     # reina curada y reprogramada: "simple" | "doble"
        self.neo = 0              # neoantígenos que su MHC-I puede llegar a presentar
        self.telomerasa_ok = True # False si un inhibidor de telomerasa la ha alcanzado
        self.senescente = False   # jubilada: no se divide, pero sigue viva y segrega SASP
        self.cebo = False         # cargada con p53 + 6-tio-dG: envenena a quien se la coma
        self.corrector = False    # lleva la maquinaria de corrección y la pasa a sus vecinas
        self.pdl1_hasta = -1      # PD-L1 ADAPTATIVO (inducido por IFN-γ) hasta este paso
        self.caf = False          # fibroblasto activado (myCAF): el que fabrica la malla
        self.cd47 = True          # "no me comas": lo lleva de serie TODA célula propia
        self.aneuploide = False   # cromosomas descuadrados: casi exclusivo del tumor, pero INVISIBLE desde fuera
        self.tipo = tipo
        self.genoma = genoma or {}
        self.er = er
        self.mutaciones = mutaciones   # pasajeras: alimentan los neoantígenos
        self.reparada = False
        self.necrosis = 0
        self.ultima_division = -999
        self.actualizar_fenotipo()

    def actualizar_fenotipo(self):
        """Traduce el genotipo (golpes por gen) a rasgos de comportamiento."""
        alt = set()
        for nombre, gen in GENES.items():
            h = self.genoma.get(nombre, 0)
            if (gen.tipo == "oncogen" and h >= 1) or h >= 2:
                alt.add(nombre)
        self.alterados = frozenset(alt)
        self.n_drivers = len(alt - GENES_INMUNES)
        self.senal_prolif = bool(alt & {"PIK3CA", "ERBB2", "MYC", "PTEN"})
        self.p53_ok = "TP53" not in alt
        self.hrd = bool(alt & {"BRCA1", "BRCA2"})
        self.cdh1_perdido = "CDH1" in alt
        self.mhc = "HLA" not in alt
        self.pdl1 = "CD274" in alt
        self.her2 = "ERBB2" in alt
        # Sobrevivir sin anclaje (llenar la luz) exige señal de crecimiento + otro fallo.
        self.inmortal = "TERT" in alt
        self.anoikis_res = self.senal_prolif and self.n_drivers >= 2
        self.invasiva = self.anoikis_res and (self.cdh1_perdido or self.n_drivers >= 4)
        f = 1.0
        if "PIK3CA" in alt:
            f *= 2.5
        if "ERBB2" in alt:
            f *= 3.0
        if "MYC" in alt:
            f *= 2.0
        if "PTEN" in alt:
            f *= 1.5
        if not self.p53_ok:
            f *= 1.3
        self.factor = min(f, 5.0)   # el ciclo celular tiene un mínimo físico (~1 día)


class Simulacion:
    def __init__(self, escenario: str = "luminal", semilla: int | None = None,
                 params: Parametros | None = None):
        esc = ESCENARIOS[escenario]
        self.escenario = escenario
        self.p = replace(params or Parametros(), estrogeno=esc.estrogeno)
        self.rng = random.Random(semilla)
        self.t = 0
        self._siguiente_id = 0
        self.celulas: dict[tuple, Celula] = {}
        self.agentes: list[Agente] = []
        self.eventos: deque[str] = deque(maxlen=12)
        self._hitos: set[str] = set()
        self.clones: dict[frozenset, int] = {}
        self.stats = dict(apoptosis_p53=0, anoikis=0, necrosis=0, muertes_inmunes=0,
                          muertes_terapia=0, colateral=0, reparaciones=0,
                          reparaciones_inutiles=0, roturas_membrana=0, vasos_nuevos=0,
                          desdiferenciaciones=0, guardianas=0, fagocitadas=0, fagocitadas_sanas=0,
                          dosis=0, senescentes=0, madres_sanas_tocadas=0, propagaciones=0,
                          cebos=0, entosis=0, muertes_cebo=0, pasos_bloqueados=0,
                          cd47_bloqueado=0, cd47_sanas=0, comidas_macrofago=0,
                          comidas_macrofago_sanas=0, cafs_activados=0, invasiones_frenadas=0,
                          entregas_troyano=0, castradas=0, castradas_sanas=0,
                          guardianas_revertidas=0, correcciones_propagadas=0,
                          corregidas_completas=0, errores_edicion=0,
                          extirpadas=0, sanas_extirpadas=0,
                          # Desglose de en qué acaba cada reparación génica (E0).
                          # "reparaciones" cuenta solo las de _reparar; el fármaco p53 y el
                          # cebo llevan contador propio para no mezclar tres mecanismos.
                          reparaciones_que_curan=0, reparaciones_p53_letales=0,
                          reparaciones_esteriles=0, reparaciones_farmaco=0,
                          reparaciones_cebo=0,
                          # ADC: molécula, no agente. Ver liberar_adc().
                          adc_alcanzadas=0, adc_muertas=0, adc_sanas=0,
                          adc_quiescentes_salvadas=0)
        self.historia: dict[str, list] = {k: [] for k in
                                          ("dia", *TUMORALES, REPARADA, NECROSIS, "T", "NK", "IA")}
        self._evento(f"Escenario '{escenario}': {esc.descripcion}")
        self._construir_tejido(esc)
        self._actualizar_mapa()
        for _ in range(60):
            self._difundir()
        self._registrar()

    # ------------------------------------------------------------------ montaje
    def _nuevo_id(self) -> int:
        self._siguiente_id += 1
        return self._siguiente_id

    def _construir_tejido(self, esc):
        n, r = self.p.tamano, self.p.radio_conducto
        self.centro = (n - 1) / 2
        yy, xx = np.mgrid[0:n, 0:n]
        self.dist = np.hypot(yy - self.centro, xx - self.centro)
        self.zona = np.full((n, n), ESTROMA, np.uint8)
        self.zona[self.dist < r] = MIO
        self.zona[self.dist < r - 1.5] = EPITELIO
        self.zona[self.dist < r - 3.2] = LUZ

        self.vasos = np.zeros((n, n), bool)
        periductal = (self.dist >= r + 2) & (self.dist < r + 3)
        lejos = self.dist >= r + 3
        azar = np.array([[self.rng.random() for _ in range(n)] for _ in range(n)])
        self.vasos[(periductal & (azar < 0.15)) | (lejos & (azar < 0.025))] = True

        for y in range(n):
            for x in range(n):
                z = self.zona[y, x]
                if z == EPITELIO:
                    self.celulas[(y, x)] = Celula(self._nuevo_id(), LUMINAL, dict(esc.germinal),
                                                  self.rng.random() < esc.fraccion_er,
                                                  madre=self.rng.random() < self.p.fraccion_madre_normal)
                elif z == MIO:
                    self.celulas[(y, x)] = Celula(self._nuevo_id(), MIOEPITELIAL)
                elif z == ESTROMA and not self.vasos[y, x] and self.rng.random() < 0.08:
                    self.celulas[(y, x)] = Celula(self._nuevo_id(), FIBROBLASTO)

        if esc.semilla:
            pos = self.rng.choice([q for q, c in self.celulas.items() if c.tipo == LUMINAL and c.madre])
            c = self.celulas[pos]
            c.genoma.update(esc.semilla)
            c.neo = sum(self.rng.random() < self.p.prob_neoantigeno for _ in esc.semilla)
            c.actualizar_fenotipo()
            self._registrar_clon(c, pos)
        # La madre SANA solo es madre en su nicho basal (ver _madre_sana_sin_nicho): las que
        # el sorteo inicial dejó en filas suprabasales nacen ya diferenciadas.
        for pos, c in self.celulas.items():
            if self._madre_sana_sin_nicho(pos, c):
                c.madre = False

        self.vasos_tumorales = np.zeros((n, n), bool)
        self.modo_vasos = "normal"
        self.O2 = np.full((n, n), 0.8)
        self.Q = np.zeros((n, n))   # quimiocinas (CXCL9/10...): atraen a las células inmunes
        self.M = np.zeros((n, n))   # densidad de malla: colágeno I + ácido hialurónico
        self.Ms = np.zeros((n, n))  # la misma, extendida: la malla no acaba en la celda del CAF
        self.P = np.zeros((n, n))   # presión intersticial: empuja hacia fuera de la masa
        self.factor_presion = 1.0   # lo baja degradar la matriz o normalizar los vasos
        self.cd47_hasta = -1        # hasta qué paso sigue circulando el anti-CD47
        self.pd1_hasta = -1         # hasta qué paso sigue circulando el anti-PD-1 (ver bloquear_pd1)

    # ------------------------------------------------------------------ bucle
    def paso(self):
        self.t += 1
        orden = list(self.celulas.items())
        self.rng.shuffle(orden)
        for pos, c in orden:
            if self.celulas.get(pos) is c:
                self._vivir(pos, c)
        self._actualizar_mapa()
        self._difundir()
        self._mover_agentes()
        self._reclutar()
        self._registrar()

    @property
    def dia(self) -> float:
        return self.t / PASOS_POR_DIA

    def _vecinos(self, pos):
        n = self.p.tamano
        y, x = pos
        for dy, dx in VECINOS:
            ny, nx = y + dy, x + dx
            if 0 <= ny < n and 0 <= nx < n:
                yield ny, nx

    def _vivir(self, pos, c: Celula):
        rng, p = self.rng, self.p
        if c.tipo == FIBROBLASTO:
            if not c.caf:
                # Dos vías de activación: el TGF-β del tumor vecino, y la propia rigidez
                # de la malla (YAP/TAZ). La segunda es la que la mantiene sola.
                tumor_cerca = any((v := self.celulas.get(q)) is not None and v.n_drivers
                                  for q in self._vecinos(pos))
                riesgo = (p.activacion_caf_tgfb if tumor_cerca else 0.0) +                     p.activacion_caf_rigidez * float(self.M[pos]) * p.lox
                if riesgo and rng.random() < riesgo:
                    c.caf = True
                    self.stats["cafs_activados"] += 1
                    self._hito("caf", "Un fibroblasto se activa (myCAF) y empieza a tejer malla")
            elif self.M[pos] < 1.0:
                self.M[pos] = min(1.0, self.M[pos] + p.produccion_matriz)
            return
        if c.tipo == MIOEPITELIAL:
            if rng.random() < 0.003:   # la capa mioepitelial se autorrepara si puede
                huecos = [q for q in self._vecinos(pos)
                          if self.zona[q] == MIO and q not in self.celulas]
                if huecos:
                    self.celulas[rng.choice(huecos)] = Celula(self._nuevo_id(), MIOEPITELIAL)
            return

        if c.necrosis:
            c.necrosis += 1
            if c.necrosis > VIDA_NECROSIS:   # los macrófagos retiran los restos
                del self.celulas[pos]
            return
        o2 = self.O2[pos]
        # MUERTE POR O2 GRADUAL (variante B, elegida el 2026-10-28). Antes la única muerte por
        # O2 era O2 < o2_necrosis con p = 0,08/paso (vida 3,1 d), y entre o2_necrosis y
        # o2_hipoxia la tasa era ~0,0001/paso (vida ~1500-2500 d): las hipóxicas eran casi
        # inmortales y acumulaban el 28 % de las célula-pasos tumorales. Datos:
        #  - las hipóxicas humanas viven 4-10 d en xenoinjerto y 3-5 d en esferoide (Durand y
        #    Sham, IJROBP 1998;42:711): la hipoxia crónica SÍ mata, en días;
        #  - con hipoxia moderada (~0,5 % O2) no hay efecto apreciable sobre la supervivencia
        #    (Papandreou et al., Cancer Res 2005;65:3171): en el borde de la hipoxia ~0.
        # Por eso la probabilidad por paso crece LINEALMENTE de 0 en o2_hipoxia a 0,08 en
        # o2_necrosis (y 0,08 por debajo). La rampa es la interpolación más simple entre esos
        # dos extremos, no un ajuste a una curva medida. Vida en la banda: ~1/p (O2 = 0,095
        # -> p = 0,04 -> 6,25 d; ver biologia.PASOS_POR_DIA).
        if o2 < p.o2_necrosis:
            p_necrosis = 0.08
        elif o2 < p.o2_hipoxia:
            p_necrosis = 0.08 * (p.o2_hipoxia - o2) / (p.o2_hipoxia - p.o2_necrosis)
        else:
            p_necrosis = 0.0
        if p_necrosis and rng.random() < p_necrosis:
            c.necrosis = 1
            self.stats["necrosis"] += 1
            self._hito("necrosis", "Primera necrosis: el centro del tumor se queda sin oxígeno (comedonecrosis)")
            return
        if c.n_drivers == 0 and rng.random() < p.prob_muerte_normal * (0.1 if c.madre else 1):
            del self.celulas[pos]
            return
        if c.n_drivers and not c.madre:
            if c.divisiones == 0 and rng.random() < p.prob_muerte_diferenciada:
                del self.celulas[pos]   # obrera agotada: envejece y muere sola
                return
            sasp = any((v := self.celulas.get(q)) is not None and v.senescente
                       for q in self._vecinos(pos))
            if rng.random() < (p.plasticidad * (3 if o2 < p.o2_hipoxia else 1)
                               + (p.prob_sasp_plasticidad * (1 - p.farmaco_anti_sasp)
                                  if sasp else 0)):
                # Plasticidad: una obrera se desdiferencia y se convierte en nueva reina.
                c.madre, c.divisiones = True, None
                self.stats["desdiferenciaciones"] += 1
                self._hito("plasticidad", "Plasticidad: una obrera tumoral se convierte en nueva reina")

        if self._madre_sana_sin_nicho(pos, c):
            c.madre = False      # ha perdido el nicho (o dejó de ser guardiana fuera de él)
        z = self.zona[pos]
        fuera_de_sitio = (z == LUZ and not c.anoikis_res) or (z in (MIO, ESTROMA) and not c.invasiva)
        if fuera_de_sitio and rng.random() < p.prob_anoikis:
            del self.celulas[pos]
            self.stats["anoikis"] += 1
            return

        if c.invasiva and z != ESTROMA and rng.random() < p.prob_invasion * 8:
            q = rng.choice(VECINOS)
            q = (pos[0] + q[0], pos[1] + q[1])
            v = self.celulas.get(q)
            if v is not None and v.tipo == MIOEPITELIAL:
                del self.celulas[q]
                self.stats["roturas_membrana"] += 1
                self._hito("rotura", "Una célula invasiva rompe la capa mioepitelial (MMPs): empieza la invasión")

        if c.corrector and p.propagacion_correccion and p.alcance_correccion:
            self._corregir_vecinas(pos, c)
        if c.guardiana and p.reversion_guardiana and rng.random() < p.reversion_guardiana:
            # El estado impuesto se apaga: la célula sigue sana, pero deja de vigilar.
            c.guardiana = None
            self.stats["guardianas_revertidas"] += 1
            self._hito("reversion", "Una guardiana pierde el estado reprogramado y deja de vigilar")
            return
        if c.guardiana and rng.random() < p.prob_fagocitosis * (
                p.coste_lectura_interna if c.guardiana in self.PUERTAS_INTERNAS else 1.0):
            self._fagocitar(pos, c, c.guardiana == "cebo")
            if self.celulas.get(pos) is not c:
                return
        if c.n_drivers and not c.cebo and not c.senescente and rng.random() < p.prob_entosis:
            self._entosis(pos, c)
            if self.celulas.get(pos) is not c:
                return
        if o2 < p.o2_hipoxia and c.n_drivers >= 2 and rng.random() < p.prob_angiogenesis:
            self._angiogenesis(pos)

        pdiv = p.prob_division * c.factor * ((1 + p.estrogeno) if c.er else 1)
        if p.fuerza_contexto and c.n_drivers:
            # El tejido recupera la autoridad: la célula sigue mutada, pero obedece.
            pdiv *= 1 - p.fuerza_contexto * self._contexto(pos)
        if o2 < p.o2_hipoxia:
            pdiv *= 0.15
        if c.madre and c.n_drivers:
            pdiv *= p.factor_quiescencia_madre   # las reinas se dividen poco: la quimio no las ve
        if rng.random() < pdiv:
            self._dividir(pos, c)

    def _contexto(self, pos) -> float:
        """Cuánta autoridad tiene el tejido sobre la célula que ocupa `pos` (0 a 1).

        No es una señal que fabrique nadie: es la ARQUITECTURA. Una célula obedece
        cuando está en su sitio (epitelio), rodeada de vecinas sanas que la sujetan, y
        anclada a una membrana basal en condiciones. Es el contexto que Weaver y Bissell
        restauraron con el anticuerpo anti-β1-integrina para revertir el fenotipo
        maligno sin corregir un solo gen.

        Y por eso el tumor grande se defiende solo: a medida que crece destruye la
        arquitectura, y con ella la autoridad que podría frenarlo."""
        if self.zona[pos] not in (EPITELIO, MIO):
            return 0.0                      # fuera de sitio no hay arquitectura que valga
        vecinas = [self.celulas.get(q) for q in self._vecinos(pos)]
        presentes = [v for v in vecinas if v is not None and not v.necrosis]
        if not presentes:
            return 0.0                      # sin vecinas no hay tejido: no hay quien mande
        sanas = sum(1 for v in presentes if v.tipo != LUMINAL or v.n_drivers == 0)
        orden = sanas / len(presentes)
        anclaje = 1.0 if any(self.zona[q] == MIO for q in self._vecinos(pos)) else 0.4
        if self.p.contexto_necesita_matriz:
            # Una malla rígida y rota no sostiene la arquitectura: la destruye. Por eso
            # la reversión funcionó en matriz de laboratorio y no se ha logrado in vivo.
            anclaje *= max(0.0, 1.0 - float(self.Ms[pos]) * self.p.lox)
        return orden * anclaje

    def _permitido(self, c: Celula, origen, destino) -> bool:
        zd = self.zona[destino]
        if self.vasos[destino]:
            return False
        if zd == EPITELIO:
            return True
        if zd == LUZ:
            return c.anoikis_res
        if zd == MIO:
            return c.invasiva and not self._frena_la_malla(destino)
        if not (c.invasiva and self.zona[origen] in (MIO, ESTROMA)):
            return False
        return not self._frena_la_malla(destino)

    def _frena_la_malla(self, destino) -> bool:
        """La malla no solo estorba a la terapia: CONTIENE al tumor. Cuanto más densa y
        más entrecruzada (LOX), más cuesta atravesarla. Esta es la pieza sin la cual el
        modelo no podría reproducir por qué eliminar los CAFs empeoró los tumores."""
        freno = self.p.contencion_matriz * float(self.Ms[destino]) * self.p.lox
        if freno and self.rng.random() < min(1.0, freno):
            self.stats["invasiones_frenadas"] += 1
            return True
        return False

    def _dividir(self, pos, c: Celula):
        if c.divisiones == 0:
            return   # obrera sin divisiones: senescencia
        huecos = [q for q in self._vecinos(pos)
                  if q not in self.celulas and self._permitido(c, pos, q)]
        if not huecos and c.senal_prolif:
            # Competición celular: el clon hiperproliferativo expulsa a una vecina sana.
            huecos = [q for q in self._vecinos(pos)
                      if (v := self.celulas.get(q)) is not None and v.tipo == LUMINAL
                      and v.n_drivers == 0 and not v.necrosis and self._permitido(c, pos, q)]
            if huecos and self.rng.random() < 0.3:
                del self.celulas[huecos[0]]
                huecos = huecos[:1]
            else:
                huecos = []
        if not huecos:
            return   # inhibición por contacto: sin sitio no hay división
        q = self.rng.choice(huecos)
        hija = Celula(self._nuevo_id(), LUMINAL, dict(c.genoma), c.er, c.mutaciones)
        hija.reparada = c.reparada
        hija.neo = c.neo
        # Mientras el anticuerpo circula, también alcanza a las células recién nacidas.
        hija.cd47 = c.cd47 if self.t < self.cd47_hasta else True
        hija.corrector = c.corrector
        # La aneuploidía se hereda, y cada mitosis sin control (p53 roto) puede generarla.
        pa = self.p.prob_aneuploidia if not c.p53_ok else self.p.prob_aneuploidia_p53ok
        hija.aneuploide = c.aneuploide or self.rng.random() < pa
        if not c.aneuploide and hija.aneuploide:
            self._hito("aneuploidia", "Una mitosis sin control reparte mal los cromosomas: aneuploidía")
        inmortal = (c.inmortal or c.madre) and c.telomerasa_ok
        golpe_farmaco = (getattr(self, "farmaco_activo", False) and self.p.farmaco_persistente
                         and inmortal and self.rng.random() < self.p.farmaco_persistente)
        if c.madre and c.n_drivers:
            # División asimétrica: casi siempre reina + obrera con divisiones contadas.
            hija.madre = self.rng.random() < self.p.prob_division_simetrica
            hija.telomerasa_ok = c.telomerasa_ok
            if hija.madre and inmortal:
                hija.divisiones = None
            else:
                hija.divisiones = self.p.divisiones_obrera if not hija.inmortal or not c.telomerasa_ok else None
            if not inmortal and c.divisiones is not None:
                c.divisiones -= 1
        elif c.madre:
            # División orientada: la hija que sigue siendo madre es la que queda pegada a la
            # membrana basal; la que cae hacia la luz sale diferenciada. La guardiana es un
            # estado IMPUESTO y se hereda aunque la hija caiga fuera del nicho.
            hija.madre = (self.rng.random() < self.p.fraccion_madre_normal
                          and (bool(c.guardiana and self.p.herencia_guardiana) or self._en_nicho(q)))
            hija.guardiana = c.guardiana if (hija.madre and self.p.herencia_guardiana) else None
        else:
            hija.telomerasa_ok = c.telomerasa_ok
            if inmortal:
                hija.divisiones = c.divisiones = None
            else:
                if c.divisiones is None:
                    c.divisiones = self.p.divisiones_obrera
                c.divisiones -= 1
                hija.divisiones = c.divisiones
        c.ultima_division = hija.ultima_division = self.t
        if golpe_farmaco:
            for cel in (c, hija):
                if not cel.p53_ok and cel.genoma.get("TP53") and \
                        self.rng.random() < self.p.farmaco_p53:
                    cel.genoma["TP53"] = 0
                    cel.actualizar_fenotipo()
                if cel.p53_ok:
                    cel.senescente = True   # se marca; la apoptosis la resuelve el paso siguiente
                    cel.divisiones = 0
                else:
                    cel.senescente, cel.divisiones = True, 0
                    if self.p.farmaco_cebo:
                        cel.cebo = True
                        self.stats["cebos"] += 1
        for cel, sitio in ((c, pos), (hija, q)):
            golpes = self._mutar(cel)
            dano = golpes + (1 if cel.hrd and self.rng.random() < self.p.prob_dano_hrd else 0)
            if dano and cel.p53_ok and self.rng.random() < 1 - (1 - self.p.prob_apoptosis_p53) ** dano:
                self.stats["apoptosis_p53"] += 1   # p53 detecta el daño y la célula se suicida
                if cel is c:
                    del self.celulas[pos]
                continue
            if cel is hija:
                self.celulas[q] = hija
            if golpes:
                self._registrar_clon(cel, sitio)

    def _mutar(self, c: Celula) -> int:
        mult = self.p.tasa_mutacion * (self.p.factor_hrd if c.hrd else 1) * (1 if c.p53_ok else 3)
        golpes = 0
        for nombre, gen in GENES.items():
            if self.rng.random() < mult * gen.peso:
                h = c.genoma.get(nombre, 0)
                if h < (1 if gen.tipo == "oncogen" else 2):
                    c.genoma[nombre] = h + 1
                    golpes += 1
                    c.neo += self.rng.random() < self.p.prob_neoantigeno
        lam = 0.3 * (10 if c.hrd else 1)
        nuevas = int(lam) + (self.rng.random() < lam - int(lam))
        c.mutaciones += nuevas
        c.neo += sum(self.rng.random() < self.p.prob_neoantigeno / 15 for _ in range(nuevas))
        if golpes:
            c.actualizar_fenotipo()
            if c.n_drivers and not c.madre and c.divisiones is None and not c.inmortal:
                c.divisiones = self.p.divisiones_obrera   # mutación en célula de vida corta
        return golpes

    # Puertas lógicas de la guardiana. Cada una es un IDENTIFICADOR distinto del tumor;
    # ver INFORME_IDENTIFICADORES.md: ninguna señal tiene a la vez cobertura y exclusividad,
    # así que la especificidad tiene que venir de exigir DOS condiciones a la vez.
    PUERTAS_INTERNAS = ("tert", "tert+nicho")

    def _puerta(self, modo: str, v: Celula, o) -> bool:
        # Señales de SUPERFICIE (lo único que un receptor tipo synNotch puede leer de verdad).
        senal_tumoral = o.neoantigenos > 0 or not o.mhc_i or o.senescente
        senal_actividad = o.proliferando or o.hipoxica or o.ligandos_estres
        if modo == "simple":
            return senal_actividad                      # una sola señal: se equivoca
        if modo in ("doble", "cebo"):
            return senal_tumoral and senal_actividad    # neoantígeno/MHC/senescencia + actividad
        if modo == "nicho":
            # Puerta Y posicional: marcador de célula madre Y fuera de la capa basal.
            # La madre SANA lleva el mismo marcador, pero vive pegada a la membrana basal.
            return o.marcador_madre and o.fuera_de_nicho
        if modo in self.PUERTAS_INTERNAS:
            # Puerta Y interna: telomerasa activa Y cromosomas descuadrados. La madre sana
            # tiene telomerasa pero NO es aneuploide; el resto del tejido no tiene telomerasa.
            # Exige leer dentro del núcleo: más lenta (coste_lectura_interna) y con falsos
            # positivos (error_lectura_interna).
            real = (v.inmortal or v.madre) and v.telomerasa_ok and v.aneuploide
            if modo == "tert+nicho":
                real = real or (o.marcador_madre and o.fuera_de_nicho)
            return real or self.rng.random() < self.p.error_lectura_interna
        return senal_tumoral and senal_actividad

    def _corregir_vecinas(self, pos, c: Celula):
        """Efecto espectador: la célula corregida pasa la corrección a sus vecinas.

        Corrige UN gen por intento, como haría un vector que transfiere una copia sana
        por las uniones comunicantes. Cuando una vecina queda limpia, se vuelve corrector
        a su vez: de ahí el crecimiento exponencial, y de ahí también el peligro."""
        p, rng = self.p, self.rng
        vecinas = [q for q in self._vecinos(pos)][:p.alcance_correccion]
        for q in vecinas:
            v = self.celulas.get(q)
            if v is None or v.tipo != LUMINAL or v.necrosis or v.corrector:
                continue
            if rng.random() >= p.propagacion_correccion:
                continue
            rotos = [g for g in v.genoma if v.genoma[g]]
            if not rotos:
                v.corrector = True        # ya estaba limpia: ahora también propaga
                continue
            gen = min(rotos, key=lambda g: ("TP53" != g, g))   # p53 primero
            if rng.random() < p.eficacia_reparacion:
                v.genoma[gen] -= 1
                self.stats["correcciones_propagadas"] += 1
                if rng.random() < p.error_edicion:
                    # Error fuera de diana: la edición rompe otro gen. El precio de que
                    # algo se copie solo es que los errores también se copian.
                    otro = rng.choice([g for g in GENES if g != gen])
                    v.genoma[otro] = v.genoma.get(otro, 0) + 1
                    self.stats["errores_edicion"] += 1
                v.actualizar_fenotipo()
                if not any(v.genoma.values()):
                    v.corrector = True
                    v.reparada = True
                    self.stats["corregidas_completas"] += 1
                    self._hito("propagacion",
                               "Una célula corregida transmite la corrección a una vecina")

    def _fagocitar(self, pos, c: Celula, ag_ceba: bool = False):
        """La guardiana se come a una vecina que pasa su puerta lógica y ocupa el hueco
        con una hija sana (competición celular / extrusión epitelial). Si ag_ceba, a las
        senescentes las deja cargadas como cebo en lugar de retirarlas."""
        for q in self._vecinos(pos):
            v = self.celulas.get(q)
            if v is None or v.tipo != LUMINAL or v.necrosis or v.guardiana:
                continue
            o = self._observar(v, q, 0, 0)
            diana = self._puerta(c.guardiana, v, o)
            if not diana:
                continue
            if ag_ceba and v.senescente and not v.cebo:
                v.cebo = True       # trampa cargada: el SASP deja de ser abono y pasa a ser anzuelo
                self.stats["cebos"] += 1
                self._hito("cebo", "Una célula senescente queda CARGADA como cebo (p53 + 6-tio-dG)")
                return
            del self.celulas[q]
            self.stats["fagocitadas"] += 1
            if v.n_drivers == 0:
                self.stats["fagocitadas_sanas"] += 1
                self.stats["colateral"] += 1
                self._hito("guardiana_error", "ERROR: una guardiana se come a una célula SANA")
            else:
                self._hito("guardiana", "Una guardiana se come a una obrera corrupta y ocupa su sitio")
            if self.rng.random() < 0.5:
                self._dividir(pos, c)
            return

    def resecar(self, cobertura: float = 0.95, margen: int = 2) -> dict:
        """Cirugía: retira un BLOQUE de tejido que cubre la mayor parte de la masa.

        No es una terapia elegante y no pretende serlo. Es la única forma conocida de
        convertir un problema de VOLUMEN en uno de superficie, que es lo que este
        simulador ha acabado señalando como el muro de los tumores grandes.

        Se modela como es: un corte contiguo con margen, no una retirada selectiva.
        Se lleva por delante tejido sano (eso es el precio) y deja fuera lo que ya había
        escapado del bloque: eso es la ENFERMEDAD RESIDUAL, y es el problema real de
        después, no el tumor que se quitó.

        CUIDADO AL INTERPRETARLA (medido 2026-10-29, ver T18 en validacion.py). Esto NO
        es una tumorectomía: el tumor acaba rodeando la luz (ocupa ~90 % del anillo en
        ~5 focos), así que el centroide cae en el centro del conducto y el disco engloba
        el conducto ENTERO. Con cobertura 0,95 y margen 2 se extirpa el 93-96 % del
        conducto, sea a 400 células o el día 120. En términos clínicos es la resección
        completa del único órgano del campo, no una cirugía conservadora, y no se puede
        comparar con las series de conservadora. Bajar la cobertura NO lo arregla: un
        margen quirúrgico real (1-2 mm = 67-133 celdas) no cabe en un campo de 110.
        """
        tumor = [q for q, c in self.celulas.items()
                 if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]
        if not tumor:
            return dict(extirpadas=0, sanas=0, residuo=0)
        cy = sum(q[0] for q in tumor) / len(tumor)
        cx = sum(q[1] for q in tumor) / len(tumor)
        dist = sorted(((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5 for q in tumor)
        radio = dist[min(len(dist) - 1, int(len(dist) * cobertura))] + margen
        extirpadas = sanas = 0
        for q, c in list(self.celulas.items()):
            if ((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5 > radio:
                continue
            if c.tipo == LUMINAL and c.n_drivers:
                extirpadas += 1
            else:
                sanas += 1
            del self.celulas[q]
        self._actualizar_mapa()
        residuo = sum(1 for c in self.celulas.values()
                      if c.tipo == LUMINAL and c.n_drivers and not c.necrosis)
        self.stats["extirpadas"] = extirpadas
        self.stats["sanas_extirpadas"] = sanas
        self._evento(f"CIRUGÍA: se extirpan {extirpadas} células tumorales y {sanas} sanas; "
                     f"queda un residuo de {residuo}")
        return dict(extirpadas=extirpadas, sanas=sanas, residuo=residuo)

    def ablacionar(self, cobertura: float = 0.95, margen: int = 2, borde: int = 2,
                   supervivencia_borde: float = 0.0, guiada: bool = False,
                   umbral_visible: int = 1) -> dict:
        """Ablación (crio o radioterapia ablativa): destruir el tumor SIN sacarlo.

        Inspirada en el ensayo FIRST (crioablación, 48 pacientes, 2026: ablación completa
        89,6 %, CDIS residual en el BORDE en el 12,5 %) y en la radioterapia ablativa de
        MD Anderson (ESTRO 2026). Misma geometría que resecar() para que sean comparables,
        con dos diferencias que son las reales:
          · lo destruido NO se retira: queda como necrosis hasta que los macrófagos lo
            limpian (VIDA_NECROSIS), y esa necrosis ya cuenta como inflamación en
            _reclutar() — el efecto "vacuna in situ" sale del modelo, no se le añade;
          · el filo de la bola de hielo no es letal del todo: en la corona exterior de
            `borde` celdas cada célula sobrevive con probabilidad `supervivencia_borde`.
            SUPUESTO: no hay dato que traduzca isotermas a celdas de 15 µm; se barre.
        Igual que resecar(), a esta escala el volumen engloba el conducto casi entero:
        destruye tejido sano y lo dice.

        `guiada=True` (2026-11-05): ablación guiada por imagen. Apunta SOLO al carcinoma
        (in situ, invasivo, reina), que es lo que se ve; las células iniciadas o
        hiperplásicas son invisibles y no arrastran el disco. Destruye todo lo que esté a
        ≤ `margen` celdas de algún foco, con varias sondas si es multifocal.
        SUPUESTO OPTIMISTA, declarado: se asume que la imagen localiza focos de muy pocas
        células; la resonancia real no ve por debajo de ~mm."""
        if guiada:
            return self._ablacion_guiada(margen, borde, supervivencia_borde, umbral_visible)
        tumor = [q for q, c in self.celulas.items()
                 if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]
        if not tumor:
            return dict(destruidas=0, sanas=0, residuo=0)
        cy = sum(q[0] for q in tumor) / len(tumor)
        cx = sum(q[1] for q in tumor) / len(tumor)
        self._centro_lesion = (cy, cx)
        dist = sorted(((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5 for q in tumor)
        radio = dist[min(len(dist) - 1, int(len(dist) * cobertura))] + margen
        self._radio_lesion = radio
        destruidas = sanas = 0
        for q, c in self.celulas.items():
            d = ((q[0] - cy) ** 2 + (q[1] - cx) ** 2) ** 0.5
            if d > radio or c.necrosis:
                continue
            if d > radio - borde and self.rng.random() < supervivencia_borde:
                continue
            if c.tipo == LUMINAL and c.n_drivers:
                destruidas += 1
            else:
                sanas += 1
            c.necrosis = 1
        self._actualizar_mapa()
        residuo = sum(1 for c in self.celulas.values()
                      if c.tipo == LUMINAL and c.n_drivers and not c.necrosis)
        self.stats["ablacionadas"] = destruidas
        self.stats["sanas_ablacionadas"] = sanas
        self._evento(f"ABLACIÓN: se destruyen in situ {destruidas} células tumorales y "
                     f"{sanas} sanas; queda un residuo de {residuo}")
        return dict(destruidas=destruidas, sanas=sanas, residuo=residuo)

    def focos_carcinoma(self) -> list[list]:
        """Focos de carcinoma = componentes conexas (8 vecinas) de células de carcinoma."""
        carcinoma = (IN_SITU, INVASIVA, REINA)
        celdas = {q for q, c in self.celulas.items()
                  if c.tipo == LUMINAL and self.estado(c, q) in carcinoma}
        focos, vistos = [], set()
        for q in celdas:
            if q in vistos:
                continue
            foco, pila = [], [q]
            vistos.add(q)
            while pila:
                y, x = pila.pop()
                foco.append((y, x))
                for dy, dx in VECINOS:
                    v = (y + dy, x + dx)
                    if v in celdas and v not in vistos:
                        vistos.add(v)
                        pila.append(v)
            focos.append(foco)
        return focos

    def _ablacion_guiada(self, margen: int, borde: int, supervivencia_borde: float,
                         umbral_visible: int = 1) -> dict:
        """`umbral_visible`: tamaño mínimo (células) de un foco para que la imagen lo vea.
        A la escala del modelo (15 µm/celda) una resonancia real (≥2 mm) no vería NINGÚN
        foco: el valor realista es 'infinito'. Se barre para ver cuánto depende de ello."""
        carcinoma = (IN_SITU, INVASIVA, REINA)
        focos = self.focos_carcinoma()
        visibles = [f for f in focos if len(f) >= umbral_visible]
        dianas = [q for f in visibles for q in f]
        if not dianas:
            self._evento(f"IMAGEN: {len(focos)} focos de carcinoma, ninguno visible "
                         f"(umbral {umbral_visible} células): no hay nada que ablacionar")
            return dict(destruidas=0, sanas=0, residuo=sum(map(len, focos)), focos_vistos=0)
        self._centro_lesion = (sum(q[0] for q in dianas) / len(dianas),
                               sum(q[1] for q in dianas) / len(dianas))
        self._radio_lesion = max(3, max(((q[0] - self._centro_lesion[0]) ** 2 +
                                         (q[1] - self._centro_lesion[1]) ** 2) ** 0.5
                                        for q in dianas))
        dmin: dict = {}
        for (y, x) in dianas:
            for dy in range(-margen, margen + 1):
                for dx in range(-margen, margen + 1):
                    d = (dy * dy + dx * dx) ** 0.5
                    if d <= margen:
                        q = (y + dy, x + dx)
                        if d < dmin.get(q, 1e9):
                            dmin[q] = d
        destruidas = sanas = 0
        for q, d in dmin.items():
            c = self.celulas.get(q)
            if c is None or c.necrosis:
                continue
            if d > margen - borde and self.rng.random() < supervivencia_borde:
                continue
            if c.tipo == LUMINAL and c.n_drivers:
                destruidas += 1
            else:
                sanas += 1
            c.necrosis = 1
        self._actualizar_mapa()
        residuo = sum(1 for q, c in self.celulas.items()
                      if c.tipo == LUMINAL and self.estado(c, q) in carcinoma)
        self.stats["ablacionadas"] = self.stats.get("ablacionadas", 0) + destruidas
        self.stats["sanas_ablacionadas"] = self.stats.get("sanas_ablacionadas", 0) + sanas
        self._evento(f"ABLACIÓN GUIADA: {len(visibles)}/{len(focos)} focos visibles "
                     f"({len(dianas)} células de carcinoma); se destruyen {destruidas} células "
                     f"con driver y {sanas} sanas; carcinoma residual {residuo}")
        return dict(destruidas=destruidas, sanas=sanas, residuo=residuo,
                    focos_vistos=len(visibles))

    def biopsia(self, cilindros: int = 6, celulas_por_cilindro: int = 12) -> dict:
        """Biopsia por vacío del lecho tumoral (la llave para NO operar, Kuerer 2022/2025).

        Toma `cilindros` muestras de `celulas_por_cilindro` celdas contiguas alrededor de
        donde estaba la lesión (el clip). Positiva si sale UNA sola célula tumoral viva.
        No se le programa ninguna tasa de falsos negativos: sale del muestreo, y se compara
        después con la de los ensayos (RESPONDER 17,8 %)."""
        if getattr(self, "_centro_lesion", None) is None:
            tumor = [q for q, c in self.celulas.items()
                     if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]
            if not tumor:
                return dict(positiva=False, muestreadas=0, real=0)
            self._centro_lesion = (sum(q[0] for q in tumor) / len(tumor),
                                   sum(q[1] for q in tumor) / len(tumor))
            dist = sorted(((q[0] - self._centro_lesion[0]) ** 2 +
                           (q[1] - self._centro_lesion[1]) ** 2) ** 0.5 for q in tumor)
            self._radio_lesion = dist[int(len(dist) * 0.95)] if dist else 5
        cy, cx = self._centro_lesion
        r = max(self._radio_lesion, 3)
        vistas = set()
        for _ in range(cilindros):
            ang = self.rng.uniform(0, 2 * np.pi)
            rad = r * self.rng.random() ** 0.5
            oy, ox = cy + rad * np.sin(ang), cx + rad * np.cos(ang)
            for k in range(celulas_por_cilindro):   # el cilindro: una fila de celdas
                vistas.add((int(round(oy + (k - celulas_por_cilindro / 2) * np.sin(ang + 1.57))),
                            int(round(ox + (k - celulas_por_cilindro / 2) * np.cos(ang + 1.57)))))
        # Lo que ve un patólogo: CARCINOMA (in situ, invasivo o célula madre tumoral). Una
        # célula "iniciada" (un driver, morfología normal) es invisible al microscopio.
        # Corregido 2026-11-05: la primera versión contaba cualquier célula con driver
        # como positiva y como "tumor real", e inflaba los falsos negativos.
        carcinoma = (IN_SITU, INVASIVA, REINA)
        positiva = any((c := self.celulas.get(q)) is not None and c.tipo == LUMINAL
                       and self.estado(c, q) in carcinoma for q in vistas)
        real = sum(1 for q, c in self.celulas.items()
                   if c.tipo == LUMINAL and self.estado(c, q) in carcinoma)
        self._evento(f"BIOPSIA: {len(vistas)} celdas muestreadas → "
                     f"{'POSITIVA' if positiva else 'negativa'} (carcinoma real: {real})")
        return dict(positiva=positiva, muestreadas=len(vistas), real=real)

    def ctdna(self, fragmentos_por_celula: float, peso_necrosis: float = 3.0) -> dict:
        """Biopsia líquida: ADN tumoral circulante en una extracción de sangre.

        Cada célula de carcinoma viva vierte ADN al morir y renovarse; las necróticas
        tumorales lo vierten de golpe (`peso_necrosis`). Los fragmentos mutantes que caen
        en el tubo siguen una Poisson; positivo si aparece ≥1 (ensayo tumor-informado).
        `fragmentos_por_celula` es la única incógnita y NO se inventa: se calibra fuera
        contra la sensibilidad clínica antes del tratamiento (ver prueba_residuo.py).
        La célula 'iniciada' también vierte ADN, pero el ensayo busca las mutaciones del
        TUMOR ORIGINAL: cuentan las de carcinoma, el mismo criterio que biopsia()."""
        carcinoma = (IN_SITU, INVASIVA, REINA)
        vivas = sum(1 for q, c in self.celulas.items() if c.tipo == LUMINAL
                    and self.estado(c, q) in carcinoma)
        muertas = sum(1 for c in self.celulas.values() if c.tipo == LUMINAL
                      and c.n_drivers and c.necrosis)
        lam = fragmentos_por_celula * (vivas + peso_necrosis * muertas)
        n = int(np.random.default_rng(self.rng.getrandbits(32)).poisson(lam))
        return dict(positiva=n > 0, fragmentos=n, esperado=lam, vivas=vivas)

    def restaurar_arquitectura(self, fuerza: float = 0.8):
        """Anti-β1-integrina: devuelve al tejido la autoridad sobre sus células.

        No corrige genes, no mata, no entra en ninguna célula. Es una señal que baña el
        tejido entero, así que su coste NO crece con el tamaño del tumor: es la única
        intervención de este simulador que trata igual de barato 10 células que 10.000.

        Es reversible: `fuerza=0` la retira y las células vuelven a mandar sobre el tejido."""
        self.p.fuerza_contexto = fuerza
        self._evento(f"ARQUITECTURA: el tejido recupera autoridad sobre sus células "
                     f"(fuerza {fuerza:.0%}); las mutaciones siguen ahí")

    def bloquear_pd1(self, dias: float | None = None):
        """Anti-PD-1 / anti-PD-L1 (dostarlimab, pembrolizumab...): quita el freno.

        Intervención MÍNIMA y aislada (añadida 2026-10-28 para el test T15 de
        validacion.py). Lo único que hace es anular, mientras el anticuerpo circula, el
        efecto que ya existía en el modelo de que una célula con PD-L1 (CD274 alterado)
        agote al linfocito T que la reconoce por TCR (o CAR). No toca el reclutamiento,
        ni los neoantígenos, ni el MHC-I, ni las NK, ni ningún parámetro: si el linfocito
        no veía a la célula, sigue sin verla. `dias=None` = bloqueo hasta el final.
        Por defecto (sin llamarla) el modelo se comporta exactamente igual que antes."""
        self.pd1_hasta = 10**12 if dias is None else self.t + int(dias * PASOS_POR_DIA)
        self._evento("TERAPIA: anti-PD-1 — el PD-L1 deja de agotar a los linfocitos T")

    def bloquear_cd47(self, eficacia: float | None = None, dirigido: bool = False):
        """Anti-CD47 (como el magrolimab): le quita a las células el cartel de 'no me comas'.

        No es selectivo y no puede serlo: CD47 lo lleva todo el cuerpo. Lo que impide que
        los macrófagos se coman el tejido sano NO es este fármaco, es la segunda señal:
        la célula sana no expone calreticulina. En la clínica el precio de esto es la
        anemia, porque los glóbulos rojos sí llevan mucho CD47 y son frágiles."""
        ef = self.p.eficacia_anti_cd47 if eficacia is None else eficacia
        self.cd47_hasta = self.t + self.p.duracion_anti_cd47 * PASOS_POR_DIA
        tocadas = sanas = 0
        for pos, c in self.celulas.items():
            if c.tipo != LUMINAL or self.rng.random() >= ef:
                continue
            if c.n_drivers == 0 and not self.p.cd47_alcanza_sanas:
                continue
            if dirigido:
                # Segunda diana del biespecífico: algo que la célula sana no muestra.
                o = self._observar(c, pos, 0, 0)
                if not (o.her2_alto or o.ligandos_estres or o.neoantigenos > 0 or not o.mhc_i):
                    continue
            c.cd47 = False
            tocadas += 1
            sanas += c.n_drivers == 0
        self.stats["cd47_bloqueado"] = tocadas
        self.stats["cd47_sanas"] = sanas
        self._evento(f"ANTI-CD47: se retira el 'no me comas' a {tocadas} células "
                     f"({sanas} de ellas SANAS: el fármaco no distingue)")

    def danar_telomeros(self):
        """Rompe el telómero AL INSTANTE en las células con telomerasa activa.

        Si p53 funciona, la célula se autodestruye. Si no, solo se jubila y se queda
        gritando (SASP). Las cualidades de Parametros.farmaco_* modifican el resultado:
        restaurar p53, dejar la carga dentro (cebo), silenciar el SASP o persistir."""
        p = self.p
        self.farmaco_activo = True
        muertas = jubiladas = sanas = 0
        for pos, c in list(self.celulas.items()):
            telomerasa = c.inmortal or c.madre
            if c.tipo != LUMINAL or not telomerasa or not c.telomerasa_ok:
                continue
            if self.rng.random() >= self.p.eficacia_6tiodg:
                continue
            if c.n_drivers == 0 and c.madre:
                sanas += 1
            if not c.p53_ok and c.genoma.get("TP53") and self.rng.random() < p.farmaco_p53:
                c.genoma["TP53"] = 0        # el propio fármaco restaura p53 (vector tipo Gendicine)
                c.actualizar_fenotipo()
                self.stats["reparaciones_farmaco"] += 1
            if c.p53_ok:
                del self.celulas[pos]       # daño irreparable + p53 = apoptosis
                self.stats["apoptosis_p53"] += 1
                muertas += 1
            else:
                c.senescente = True         # sin p53 solo se jubila: hay que retirarla aparte
                c.divisiones = 0
                if p.farmaco_cebo:
                    c.cebo = True           # la carga se queda dentro: la jubilada es una mina
                    self.stats["cebos"] += 1
                jubiladas += 1
        self.stats["madres_sanas_tocadas"] += sanas
        self._evento(f"6-TIO-dG: {muertas} células se autodestruyen (p53 intacto), "
                     f"{jubiladas} solo se jubilan (p53 roto), {sanas} madres sanas afectadas")

    def inhibir_telomerasa(self):
        """La 'castración': bloquea la telomerasa. Toda célula alcanzada recupera un contador
        de divisiones finito y se jubila al agotarlo. Daño colateral: las células madre SANAS
        del tejido también usan telomerasa."""
        tocadas = sanas = 0
        for c in self.celulas.values():
            if c.tipo != LUMINAL or self.rng.random() >= self.p.eficacia_antitelomerasa:
                continue
            c.telomerasa_ok = False
            tocadas += 1
            if c.divisiones is None:
                c.divisiones = self.p.divisiones_obrera
            if c.n_drivers == 0 and c.madre:
                sanas += 1
        self.stats["madres_sanas_tocadas"] += sanas
        self._evento(f"INHIBIDOR DE TELOMERASA: {tocadas} células pierden el rebobinado "
                     f"({sanas} madres sanas afectadas)")

    def _entosis(self, pos, c: Celula):
        """La célula tumoral se traga a una vecina (entosis, descrita en cáncer de mama).
        Si era un cebo, recibe p53 reparado + telómero roto: el banquete la mata."""
        cebos = [q for q in self._vecinos(pos)
                 if (v := self.celulas.get(q)) is not None and v.cebo]
        if not cebos:
            return
        q = self.rng.choice(cebos)
        del self.celulas[q]                     # la digiere en su lisosoma
        self.stats["entosis"] += 1
        if self.rng.random() >= self.p.eficacia_cebo:
            return
        if c.genoma.get("TP53"):
            c.genoma["TP53"] = 0                # carga 1: p53 restaurado
            c.actualizar_fenotipo()
            self.stats["reparaciones_cebo"] += 1
        if c.p53_ok and (c.inmortal or c.madre):
            del self.celulas[pos]               # carga 2: telómero roto + p53 = autodestrucción
            self.stats["muertes_cebo"] += 1
            self.stats["apoptosis_p53"] += 1
            self._hito("cebo_mata", "Una célula tumoral se come un cebo y se autodestruye")
            return
        c.senescente, c.divisiones = True, 0
        if self.rng.random() < self.p.prob_cebo_encadena:
            c.cebo = True                       # la comedora se vuelve cebo: la trampa se propaga
            self.stats["cebos"] += 1

    def liberar_adc(self, diana: str = "trop2", afinidad: str = "alta",
                    prevalencia: float | None = None, dosis: float = 0.35,
                    radio_espectador: int = 8, depende_ciclo: bool = True) -> dict:
        """Anticuerpo conjugado (ADC). NO es un agente que camina: es una molécula.

        Tres diferencias con las terapias celulares del modelo, y las tres importan:

        1. PENETRA POCO Y DESDE EL VASO. Medido en humanos con variantes de afinidad
           de un anti-HER2: alta afinidad se queda a ~20 µm del vaso, baja afinidad
           llega a ~85 µm. Es la "barrera del sitio de unión": agarrarse fuerte a la
           primera célula impide llegar al fondo. A 15 µm/celda → 1,3 y 5,7 celdas.
        2. EL VENENO SE DERRAMA 50-200 µm, así que mata también a la vecina que NO
           expresa la diana. Es lo que resolvió la heterogeneidad antigénica, y la
           razón de que la FDA apruebe estos fármacos SIN exigir expresión de diana.
        3. EL PAYLOAD DEPENDE DE LA DIVISIÓN. Auristatinas, maitansinoides e
           inhibidores de topoisomerasa I necesitan una célula en ciclo. Nunca se han
           probado contra células dormidas. Por eso `depende_ciclo`.

        `prevalencia` es la fracción de células tumorales que expresan la diana EN
        CRIBADO PROSPECTIVO, que no es la de la literatura: con folato-α se pasó de un
        80 % publicado a un 10 % real medido en 80 pacientes consecutivas, y la tasa de
        respuesta fue 0. Por defecto se usa la del cribado, no la publicada.
        """
        p, rng = self.p, self.rng
        alcance = {"alta": 1, "media": 3, "baja": 6}.get(afinidad, 3)   # celdas desde el vaso
        prev = PREVALENCIA_DIANA.get(diana, 0.56) if prevalencia is None else prevalencia

        # Distancia al vaso más cercano (la molécula sale de ahí y difunde)
        dist = np.full((p.tamano, p.tamano), 1e9)
        frente = [tuple(q) for q in np.argwhere(self.vasos)]
        for q in frente:
            dist[q] = 0
        d = 0
        while frente and d <= alcance + radio_espectador:
            d += 1
            nuevo_frente = []
            for y, x in frente:
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    q = (y + dy, x + dx)
                    if 0 <= q[0] < p.tamano and 0 <= q[1] < p.tamano and dist[q] > d:
                        dist[q] = d
                        nuevo_frente.append(q)
            frente = nuevo_frente

        # 1) Células ALCANZADAS: dentro del alcance Y que expresan la diana
        alcanzadas = []
        for pos, c in self.celulas.items():
            if c.tipo != LUMINAL or not c.n_drivers or c.necrosis:
                continue
            if dist[pos] > alcance:
                continue
            if rng.random() < prev:              # ¿expresa la diana?
                alcanzadas.append(pos)

        # 2) El veneno se derrama: zona letal alrededor de cada célula alcanzada
        zona = set()
        for pos in alcanzadas:
            for dy in range(-radio_espectador, radio_espectador + 1):
                for dx in range(-radio_espectador, radio_espectador + 1):
                    if dy * dy + dx * dx <= radio_espectador * radio_espectador:
                        zona.add((pos[0] + dy, pos[1] + dx))

        muertas = sanas = quiescentes_salvadas = 0
        for q in list(zona):
            c = self.celulas.get(q)
            if c is None or c.tipo != LUMINAL or c.necrosis:
                continue
            if depende_ciclo and (c.divisiones == 0 or c.senescente
                                  or self.t - c.ultima_division >= 40):
                quiescentes_salvadas += c.n_drivers > 0
                continue                          # el payload no la toca: no se divide
            if rng.random() >= dosis:
                continue
            if c.n_drivers == 0:
                sanas += 1
                self.stats["colateral"] += 1
            del self.celulas[q]
            muertas += 1
        self._actualizar_mapa()
        self.stats["adc_alcanzadas"] += len(alcanzadas)
        self.stats["adc_muertas"] += muertas
        self.stats["adc_sanas"] += sanas
        self.stats["adc_quiescentes_salvadas"] += quiescentes_salvadas
        self._evento(f"ADC anti-{diana} (afinidad {afinidad}): alcanza {len(alcanzadas)} "
                     f"células con diana, mata {muertas} ({sanas} sanas); "
                     f"{quiescentes_salvadas} tumorales quiescentes SOBREVIVEN")
        return dict(alcanzadas=len(alcanzadas), muertas=muertas, sanas=sanas,
                    quiescentes=quiescentes_salvadas)

    def tratar_matriz(self, modo: str, intensidad: float = 0.6):
        """Las cinco maneras de atacar la malla. NO son equivalentes, y ahí está el asunto:
        en la clínica todo lo que la DESTRUYE ha fracasado y lo que la ABLANDA funciona.

        'hialuronidasa'  degrada el ácido hialurónico (PEGPH20). Baja la presión de golpe,
                         pero los CAFs siguen ahí y vuelven a tejer. Fase III falló (2019).
        'anti_lox'       impide el entrecruzamiento (simtuzumab). No retira material:
                         lo destensa. Baja presión Y contención a la vez.
        'losartan'       reduce la PRODUCCIÓN de colágeno y HA. Barato, positivo en
                         páncreas con quimio (Murphy 2019).
        'anti_tgfb'      corta la orden de activación: los fibroblastos no se vuelven CAF.
        'matar_caf'      elimina al fabricante. Es lo que empeoró los tumores en ratones.
        """
        if modo == "hialuronidasa":
            self.M *= 1 - intensidad
            texto = f"se degrada el {intensidad:.0%} de la malla; los CAFs siguen ahí"
        elif modo == "anti_lox":
            self.p.lox *= 1 - intensidad
            texto = f"entrecruzamiento al {self.p.lox:.0%}: la malla se destensa sin retirarse"
        elif modo == "losartan":
            self.p.produccion_matriz *= 1 - intensidad
            texto = f"producción de malla al {1 - intensidad:.0%}"
        elif modo == "anti_tgfb":
            self.p.activacion_caf_tgfb *= 1 - intensidad
            self.p.activacion_caf_rigidez *= 1 - intensidad
            texto = "se corta la activación de nuevos CAFs"
        elif modo == "matar_caf":
            return self.degradar_matriz(intensidad)
        else:
            raise ValueError(f"modo de matriz desconocido: {modo}")
        self._evento(f"MALLA ({modo}): {texto}")

    def degradar_matriz(self, fraccion: float = 0.25):
        """Enzima tipo PEGPH20/hialuronidasa: deshace la matriz densa, baja la presión y
        abre huecos por los que puede circular una célula terapéutica. No mata nada:
        retira fibroblastos y descomprime el estroma."""
        quitados = 0
        for pos, c in list(self.celulas.items()):
            if c.tipo == FIBROBLASTO and self.rng.random() < fraccion:
                del self.celulas[pos]
                quitados += 1
        self.factor_presion *= 1 - self.p.caida_presion_matriz * fraccion
        self._evento(f"MATRIZ: se degradan {quitados} fibroblastos; la presión baja al "
                     f"{self.factor_presion:.0%} y se abren huecos de paso")

    def tratar_vasos(self, modo: str):
        """normal · antiangiogenico (bloquea VEGF: retira vasos tumorales y no deja crecer más)
        · normalizar (los vasos tumorales riegan bien) · asfixiar (disruptor vascular cerca del tumor)
        · reconducir (desvía oxígeno del tejido sano lejano a la zona tumoral: menos hipoxia)."""
        self.modo_vasos = modo
        if modo == "normalizar":
            # Normalizar la vasculatura cierra las fugas: baja la presión Y mejora el riego.
            self.factor_presion *= 1 - self.p.caida_presion_normalizar
        if modo == "antiangiogenico":
            self.vasos &= ~self.vasos_tumorales
            self.vasos_tumorales[:] = False
        self._evento(f"VASOS: modo {modo}")

    def _angiogenesis(self, pos):
        if self.modo_vasos == "antiangiogenico":
            return
        y, x = pos
        ny, nx = y + self.rng.randint(-5, 5), x + self.rng.randint(-5, 5)
        n = self.p.tamano
        if 0 <= ny < n and 0 <= nx < n and self.zona[ny, nx] == ESTROMA and not self.vasos[ny, nx]:
            v = self.celulas.get((ny, nx))
            if v is None or v.tipo == FIBROBLASTO:
                self.celulas.pop((ny, nx), None)
                self.vasos[ny, nx] = True
                self.vasos_tumorales[ny, nx] = True
                self.stats["vasos_nuevos"] += 1
                self._hito("vegf", "El tumor hipóxico segrega VEGF: crece un vaso nuevo (angiogénesis)")

    # ------------------------------------------------------------------ campos
    def estado(self, c: Celula, pos) -> int:
        if c.tipo == FIBROBLASTO:
            return FIBRO
        if c.tipo == MIOEPITELIAL:
            return MIOEP
        if c.necrosis:
            return NECROSIS
        if c.n_drivers == 0:
            return REPARADA if c.reparada else NORMAL
        if c.madre:
            return REINA
        if self.zona[pos] in (MIO, ESTROMA):
            return INVASIVA
        if c.anoikis_res:
            return IN_SITU
        if c.senal_prolif:
            return HIPERPLASIA
        return INICIADA

    def _actualizar_mapa(self):
        self.codigo = np.zeros(self.zona.shape, np.uint8)
        for pos, c in self.celulas.items():
            self.codigo[pos] = self.estado(c, pos)
        self.conteo = np.bincount(self.codigo.ravel(), minlength=11)

    def _difundir(self, iteraciones: int = 6):
        k = CONSUMO_O2[self.codigo]
        fuente = FUENTE_QUIMIOCINA[self.codigo]
        riego = self.vasos.astype(float)
        if self.modo_vasos not in ("normalizar", "reconducir"):
            riego[self.vasos_tumorales] = self.p.o2_vaso_tumoral
        dirigido = None
        if self.modo_vasos in ("asfixiar", "reconducir"):
            cerca = np.isin(self.codigo, TUMORALES)
            for _ in range(4):
                cerca = (cerca | np.roll(cerca, 1, 0) | np.roll(cerca, -1, 0)
                         | np.roll(cerca, 1, 1) | np.roll(cerca, -1, 1))
            if self.modo_vasos == "asfixiar":
                riego[cerca & self.vasos] = self.p.o2_vaso_asfixiado
            else:
                # Reconducir: el flujo se desvía del tejido sano lejano hacia la zona tumoral.
                riego[~cerca & self.vasos] = 0.85
                dirigido = cerca.astype(float)
        self.M *= 1 - self.p.recambio_matriz
        # La malla que teje un CAF no se queda en su celda: impregna el entorno, incluida
        # la membrana basal que el tumor tiene que atravesar para invadir.
        self.Ms = self.M.copy()
        for eje in (0, 1):
            self.Ms = np.maximum(self.Ms, 0.7 * np.roll(self.M, 1, eje))
            self.Ms = np.maximum(self.Ms, 0.7 * np.roll(self.M, -1, eje))
        self._difundir_presion(iteraciones)
        for campo, sumidero, aporte in ((self.O2, k, None), (self.Q, 0.08, fuente)):
            for _ in range(iteraciones):
                pad = np.pad(campo, 1, mode="edge")
                lap = pad[:-2, 1:-1] + pad[2:, 1:-1] + pad[1:-1, :-2] + pad[1:-1, 2:] - 4 * campo
                campo += 0.2 * lap - sumidero * campo
                if aporte is not None:
                    campo += 0.05 * aporte
                else:
                    campo[self.vasos] = riego[self.vasos]
                    if dirigido is not None:
                        campo += 0.03 * dirigido * (1 - campo)
                np.clip(campo, 0.0, 1.0, out=campo)

    def _difundir_presion(self, iteraciones: int = 6):
        """Campo de presión intersticial. Lo generan las células tumorales (vasos porosos
        que gotean plasma) y los fibroblastos (matriz densa que impide el drenaje); fuera
        de la masa los linfáticos sanos lo drenan deprisa. El resultado es una meseta de
        presión alta dentro y una caída brusca en el borde: justo el perfil que se mide
        en tumores reales, y la razón de que el borde se trate solo y el centro no."""
        # El ácido hialurónico de la malla es un gel que retiene agua: la presión sale
        # de ahí, no de que haya fibroblastos. Por eso se puede ABLANDAR sin destruir.
        fuente = np.isin(self.codigo, TUMORALES).astype(float)
        fuente += 1.2 * self.Ms * self.p.lox
        fuente *= self.p.presion_intersticial * self.factor_presion
        drenaje = np.where(fuente > 0, 0.02, 0.25)   # linfáticos: colapsados dentro, bien fuera
        for _ in range(iteraciones):
            pad = np.pad(self.P, 1, mode="constant")
            lap = pad[:-2, 1:-1] + pad[2:, 1:-1] + pad[1:-1, :-2] + pad[1:-1, 2:] - 4 * self.P
            self.P += 0.2 * lap + 0.08 * fuente - drenaje * self.P
            np.clip(self.P, 0.0, 1.0, out=self.P)

    def _frena_la_presion(self, origen, destino, troyano: float = 0.0) -> bool:
        """¿La corriente de salida impide este paso? Solo frena a quien intenta ENTRAR
        (moverse hacia donde la presión sube); salir es cuesta abajo y es gratis.

        A un vehículo troyano casi no le afecta: no nada contra la corriente, lo lleva
        dentro una célula que el tumor ha reclutado y a la que le han abierto el paso."""
        if not self.p.presion_intersticial or troyano >= 1.0:
            return False
        subida = float(self.P[destino] - self.P[origen])
        if subida <= 0:
            return False
        freno = self.p.efecto_presion * subida * (1 - troyano)
        return self.rng.random() < min(1.0, freno)

    # ------------------------------------------------------------------ agentes
    def _crear_agente(self, tipo: str, pos) -> Agente:
        if tipo == "T":
            return Agente("T", "linfocito T CD8+", politica_linfocito_t, frozenset({"TCR"}),
                          pos=pos, expansion=0.2, vida=160)
        if tipo == "MAC":
            # El tumor polariza a la mayoría hacia M2: presentes, pero no fagocitan.
            m2 = self.rng.random() < self.p.fraccion_macrofagos_m2
            return Agente("MAC", "macrófago M2 (no come)" if m2 else "macrófago",
                          politica_macrofago,
                          frozenset() if m2 else frozenset({"FAGOCITO"}), pos=pos, vida=240)
        return Agente("NK", "célula NK", politica_nk, frozenset({"NK"}), pos=pos, vida=80)

    def _reclutar(self):
        naturales = sum(1 for a in self.agentes if a.tipo in ("T", "NK") and not a.clonal)
        inflamacion = (self.conteo[IN_SITU] + self.conteo[INVASIVA] + self.conteo[NECROSIS]) / 60
        if naturales < self.p.max_inmunes_naturales and \
                self.rng.random() < self.p.reclutamiento_inmune * (1 + inflamacion):
            tipo = "T" if self.rng.random() < 0.6 else "NK"
            self.agentes.append(self._crear_agente(tipo, self._vaso_aleatorio()))
        # Los macrófagos están siempre ahí: no hay que infundirlos, hay que desbloquearlos.
        macros = sum(1 for a in self.agentes if a.tipo == "MAC")
        if macros < self.p.max_macrofagos and                 self.rng.random() < self.p.reclutamiento_macrofagos * (1 + inflamacion):
            self.agentes.append(self._crear_agente("MAC", self._vaso_aleatorio()))

    def _vaso_aleatorio(self, cerca: bool = False):
        ys, xs = np.nonzero(self.vasos)
        if cerca:
            d = self.dist[ys, xs]
            sel = np.argsort(d)[:40]
            ys, xs = ys[sel], xs[sel]
        i = self.rng.randrange(len(ys))
        return int(ys[i]), int(xs[i])

    def _posiciones_entrega(self, n: int, via: str):
        """De dónde salen las células terapéuticas.

        'vena'          por la sangre: aparecen en un vaso y tienen que abrirse paso (lo de siempre)
        'intratumoral'  aguja guiada por ecografía: aparecen DENTRO de la masa, repartidas
        'aguja'         una sola punción: todas en un punto del tumor (lo realista con una aguja)
        """
        tumor = [q for q, c in self.celulas.items()
                 if c.tipo == LUMINAL and c.n_drivers and not c.necrosis]
        if via == "vena" or not tumor:
            return [self._vaso_aleatorio() for _ in range(n)]
        if via == "aguja":
            centro = self.rng.choice(tumor)
            cerca = sorted(tumor, key=lambda q: sum((a - b) ** 2 for a, b in zip(q, centro)))[:n]
            return [self.rng.choice(cerca) for _ in range(n)]
        return [self.rng.choice(tumor) for _ in range(n)]

    def liberar_terapia(self, nombre: str, cantidad: int | None = None, via: str = "vena"):
        """Admite combinaciones: "guardiana+memoria". `via`: vena | intratumoral | aguja."""
        if "+" in nombre:
            for parte in nombre.split("+"):
                self.liberar_terapia(parte, cantidad, via)
            return
        cfg = TERAPIAS[nombre]
        n = cantidad or cfg["n"]
        self.stats["dosis"] += n
        for pos in self._posiciones_entrega(n, via):
            ag = cfg["crear"]()
            ag.pos = pos
            self.agentes.append(ag)
        self._evento(f"TERAPIA: se infunden {cantidad or cfg['n']} × {nombre}")

    def _observar(self, c: Celula, pos, dy, dx) -> Observacion:
        hipoxica = self.O2[pos] < self.p.o2_hipoxia
        neo = c.neo if (c.mhc and c.n_drivers > 0) else 0
        if c.madre:
            neo //= 3   # las madres tumorales presentan menos antígeno (MHC-I bajo)
        return Observacion(dy, dx, c.id, NOMBRE_TEJIDO[c.tipo], c.mhc, neo, self._pdl1(c), c.her2,
                           c.senescente or c.n_drivers >= 2 or (hipoxica and c.n_drivers > 0),
                           self.t - c.ultima_division < 8, hipoxica, c.madre, c.senescente,
                           c.cd47, self._expone_comeme(c, hipoxica),
                           c.madre and not self._en_nicho(pos))

    @staticmethod
    def _expone_comeme(c: Celula, hipoxica: bool) -> bool:
        """Calreticulina en la superficie: la señal de 'cómeme'. La expone la célula
        jubilada (senescente), la dañada y la muy estresada. Una célula sana que sigue
        su ciclo normal NO la expone: esa es la mitad de la puerta que la protege
        aunque le hayan bloqueado el CD47."""
        if c.n_drivers == 0 and not c.senescente:
            return False
        return c.senescente or c.divisiones == 0 or c.cebo or (c.n_drivers >= 2 and hipoxica)

    def _madre_sana_sin_nicho(self, pos, c: Celula) -> bool:
        """La identidad de célula madre SANA depende del nicho: contacto con el mioepitelio y
        la lámina basal (Notch/Wnt). Fuera de él se diferencia y pierde CD44+/CD24-/ALDH1.
        La reina TUMORAL no depende del nicho (independencia de nicho = rasgo del tumor), y la
        guardiana tampoco: su estado es impuesto. Corrige un defecto (2026-10-28): el modelo
        dejaba madres sanas suprabasales indefinidamente con el marcador de madre, que es
        exactamente lo que la puerta "nicho" busca -> falso positivo (T13 con el flag
        adaptativo: semilla 1, día 232; puertas.py sin tratar: 6/301 sanas tocadas el día 120).
        No da a la guardiana ninguna percepción nueva ni toca la puerta.
        MEDIDO: daño a sanas en la terapia de T13, semillas 1-24: antes 3/24 semillas sin flag
        (7, 16, 18) y 2/24 con flag (1, 16); después 0/24 en ambos. Control del tumor igual
        (semillas 9-24: 12/16 → 11/16 sin flag, 12/16 → 12/16 con flag; a igual tamaño el
        día 120 se controla lo mismo). puertas.py día 120: nicho toca 0 sanas (antes 6)."""
        return (c.madre and c.n_drivers == 0 and not c.guardiana and c.tipo == LUMINAL
                and not self._en_nicho(pos))

    def _en_nicho(self, pos) -> bool:
        """El nicho de la célula madre mamaria es la capa basal: epitelio EN CONTACTO con
        el mioepitelio. Una reina tumoral no respeta esa arquitectura y acaba fuera."""
        if self.zona[pos] != EPITELIO:
            return False
        return any(self.zona[q] == MIO for q in self._vecinos(pos))

    def _percibir(self, ag: Agente) -> Percepcion:
        y, x = ag.pos
        n = self.p.tamano
        vecinas = []
        for dy in range(-2, 3):
            for dx in range(-2, 3):
                c = self.celulas.get((y + dy, x + dx))
                if c is not None and not c.necrosis:
                    vecinas.append(self._observar(c, (y + dy, x + dx), dy, dx))
        mejor, quimio = self.Q[y, x] + 1e-4, None
        for dy, dx in VECINOS:
            ny, nx = y + dy, x + dx
            if 0 <= ny < n and 0 <= nx < n and self.Q[ny, nx] > mejor:
                mejor, quimio = self.Q[ny, nx], (dy, dx)
        return Percepcion(vecinas, quimio, (_signo(self.centro - y), _signo(self.centro - x)),
                          NOMBRE_ZONA[self.zona[y, x]], float(self.O2[y, x]), self.dia, self.rng)

    def _mover_agentes(self):
        n = self.p.tamano
        vivos = []
        adaptativa = self.p.inmunidad_adaptativa
        if adaptativa:
            self._n_clonales = sum(1 for a in self.agentes if a.clonal and a.vida > 1)
        for ag in self.agentes:
            ag.vida -= 1
            if ag.vida <= 0 or ag.agotamiento >= 1:
                continue
            vivos.append(ag)
            if adaptativa and ag.tipo == "T" and ag.estado_t:
                self._ciclo_t(ag, vivos)
            if ag.entrega and self._entregar_troyano(ag, vivos):
                continue
            accion = ag.politica(ag, self._percibir(ag))
            if isinstance(accion, Mover):
                destino = (min(n - 1, max(0, ag.pos[0] + _signo(accion.dy))),
                           min(n - 1, max(0, ag.pos[1] + _signo(accion.dx))))
                if not self._frena_la_presion(ag.pos, destino, ag.troyano):
                    ag.pos = destino
                    if ag.oxigeno:
                        # La bombona va soltando O2 por donde pasa.
                        self.O2[ag.pos] = min(1.0, self.O2[ag.pos] + ag.oxigeno)
                else:
                    self.stats["pasos_bloqueados"] += 1
                continue
            if not isinstance(accion, (Atacar, Analizar, Reparar)) or \
                    max(abs(accion.dy), abs(accion.dx)) > 1:
                continue
            q = (ag.pos[0] + accion.dy, ag.pos[1] + accion.dx)
            c = self.celulas.get(q)
            if c is None or c.necrosis:
                continue
            if isinstance(accion, Atacar):
                self._atacar(ag, c, q, vivos)
            elif isinstance(accion, Analizar):
                # Secuenciar: devuelve TODOS los golpes, incluidos alelos heterocigotos.
                ag.memoria.setdefault("genomas", {})[c.id] = {g: h for g, h in c.genoma.items() if h}
            else:
                self._reparar(ag, c, q, accion.gen)
        self.agentes = vivos

    def _entregar_troyano(self, ag: Agente, vivos: list) -> bool:
        """El vehículo ha llegado: suelta la carga y desaparece. La terapia real nace
        AQUÍ DENTRO, que es la gracia — ya no tiene que atravesar nada."""
        rodeado = sum(1 for q in self._vecinos(ag.pos)
                      if (c := self.celulas.get(q)) is not None and c.n_drivers)
        if rodeado < 5:
            return False
        cfg = TERAPIAS[ag.entrega]
        for _ in range(3):
            if len(vivos) >= MAX_AGENTES:
                break
            nuevo = cfg["crear"]()
            nuevo.pos = ag.pos
            vivos.append(nuevo)
        vivos.remove(ag)
        self.stats["entregas_troyano"] += 1
        self._hito("troyano", f"Un vehículo llega al interior del tumor y suelta {ag.entrega}")
        return True

    def _reconoce(self, ag: Agente, c: Celula, q) -> str | None:
        r = ag.receptores
        if "INESPECIFICO" in r:
            return "inespecifico"
        if c.tipo != LUMINAL:
            return None
        obs = self._observar(c, q, 0, 0)
        if "FAGOCITO" in r:
            # Puerta Y natural: hace falta que falte el freno Y que esté el acelerador.
            return "fagocito" if (obs.comeme and not obs.cd47) else None
        if "HER2" in r and c.her2:
            return "car"
        if "TCR" in r and obs.neoantigenos > 0:
            return "tcr"
        if "NK" in r and (not c.mhc or obs.ligandos_estres):
            return "nk"
        return None

    def _atacar(self, ag: Agente, c: Celula, q, vivos: list):
        via = self._reconoce(ag, c, q)
        if via is None:
            return
        adaptativa = self.p.inmunidad_adaptativa
        t_natural = adaptativa and ag.tipo == "T"
        if t_natural and via == "tcr":
            self._estimular_t(ag)          # cada encuentro con su antígeno cuenta (Wherry 2011)
        pdl1 = self._pdl1(c)
        if via in ("tcr", "car") and pdl1 and not ag.resistente_pdl1 and self.t >= self.pd1_hasta:
            self.stats["frenos_pdl1"] = self.stats.get("frenos_pdl1", 0) + 1
            if c.pdl1:
                self.stats["frenos_pdl1_geneticos"] = self.stats.get("frenos_pdl1_geneticos", 0) + 1
            if not t_natural:
                ag.agotamiento += 0.2   # PD-L1 frena al linfocito: se agota
            # Con inmunidad adaptativa el T natural NO desaparece por el freno: queda
            # inhibido este turno. Lo que lo agota es la estimulación crónica (_estimular_t);
            # sin PD-1 se agota igual o más (Odorizzi, J Exp Med 2015;212:1125).
            # DESENGANCHE (corrección 2026-10-28; generalizada el 2026-10-29 a CUALQUIER
            # agente frenado, no sólo a los T endógenos: los de terapia son de tipo "IA" y
            # nunca lo recibían. Censo de politicas.TERAPIAS: el único agente de terapia
            # frenable por PD-L1 es car_her2, el único con receptor y sin resistente_pdl1.
            # Medido antes de corregirlo: 1.556 frenos en 12 semillas, el 45 % en pasos
            # CONSECUTIVOS sobre la misma célula, rachas de hasta 25 pasos = 6,25 días.
            # Antes, la política del T no tenía
            # memoria y elegía siempre la diana visible más cercana: tras un freno volvía
            # a la MISMA célula PD-L1+ en el paso siguiente. Medido (luminal, neo 0,1,
            # semillas 1-3, 250 días): sin el flag, cada T frenado recibía sus 5 frenos
            # en 5 pasos seguidos sobre la misma célula y moría; con el flag adaptativo,
            # el 94 % de los ~114.000 frenos eran repeticiones sobre la misma célula, en
            # rachas de hasta 165 pasos (41 días, la vida entera del T).
            # POR QUÉ ES HONESTO: la señal PD-1 inhibe al CTL pero no lo ancla. Al
            # contrario, PD-1/PD-L1 BLOQUEA la señal de parada del TCR y el T sigue
            # moviéndose (Fife, Nat Immunol 2009;10:1185); los contactos T-diana duran
            # minutos-horas, no semanas. Se conserva TODO lo demás: el freno sigue
            # costando agotamiento (sin flag) o contando como estímulo (con flag), y el
            # T puede volver a esa célula pasado el desenganche si sigue ahí.
            # EXPECTATIVA, escrita ANTES de ejecutar: menos frenos totales (sobre todo
            # con el flag), T frenados que sobreviven más y matan algo más de células
            # PD-L1−; efecto pequeño en el tumor porque las alternativas PD-L1− al
            # alcance eran escasas (adyacentes en 3-8 % de los frenos sin flag; <2 % con
            # flag por el halo de IFN-γ). Veredictos de la batería sin cambios; riesgo
            # en T12 (con flag ya iba justo, recaída 1/3). T15 sigue fallando.
            # MEDIDO DESPUÉS: rachas sobre la misma célula 0 (antes hasta 165 pasos).
            # Sin flag: T1-T13 pasan igual; T15 falla igual. Con flag: T1-T12 pasan y
            # T13 FALLA por 1 célula sana (semilla 1, día 232) comida por una GUARDIANA
            # de puerta "nicho", no por un T: la trayectoria cambia y destapa un falso
            # positivo de esa puerta. No se ha tocado nada aquí para arreglarlo: la
            # causa era una madre SANA suprabasal y se corrigió en el tejido
            # (_madre_sana_sin_nicho, 2026-10-28), no en el T ni en la puerta.
            ag.memoria.setdefault("frenado_por", {})[c.id] = \
                self.dia + self.p.dias_desenganche_pd1
            self._hito("pdl1", "Un clon con PD-L1 alto agota a los linfocitos T (escape inmune)")
            return
        if via in ("tcr", "car") and pdl1 and not ag.resistente_pdl1:
            # Solo se llega aquí con el anti-PD-1 circulando: el freno estaba y se ha quitado.
            self.stats["frenos_pdl1_bloqueados"] = self.stats.get("frenos_pdl1_bloqueados", 0) + 1
        prob = self.p.prob_muerte_ataque
        if t_natural and ag.estado_t == "tex_prog":
            prob *= self.p.eficacia_tex_prog
        elif t_natural and ag.estado_t == "tex_term":
            prob *= self.p.eficacia_tex_term
        mata = self.rng.random() < prob
        if adaptativa and via in ("tcr", "car"):
            # El linfocito que ha reconocido suelta IFN-γ: PD-L1 en la diana y alrededor.
            self._inducir_pdl1(q)
        if not mata:
            return
        if ag.castra:
            # No la mata: la jubila. Telómero roto, sin más divisiones, y marcada con el
            # "cómeme" para que los macrófagos la retiren.
            if c.senescente:
                return
            c.senescente = True
            c.divisiones = 0
            c.telomerasa_ok = False
            ag.muertes += 1
            self.stats["castradas"] += 1
            if c.n_drivers == 0:
                self.stats["castradas_sanas"] += 1
                self.stats["colateral"] += 1
            self._hito("castrada", "Una célula tumoral queda estéril: ya no puede dividirse")
            if ag.expansion and len(vivos) + 1 < MAX_AGENTES and                     self.rng.random() < ag.expansion:
                vivos.append(replace(ag, memoria={}, agotamiento=0.0, muertes=0, reparaciones=0))
            return
        del self.celulas[q]
        ag.muertes += 1
        if ag.tipo == "MAC":
            self.stats["comidas_macrofago"] += 1
            self.stats["comidas_macrofago_sanas"] += c.n_drivers == 0
            self._hito("macrofago", "Un macrófago se come una célula sin 'no me comas' y con 'cómeme'")
        self.stats["muertes_inmunes" if ag.tipo in ("T", "NK", "MAC") else "muertes_terapia"] += 1
        if c.tipo != LUMINAL or c.n_drivers == 0:
            self.stats["colateral"] += 1
            self._hito("colateral", f"Daño colateral: {ag.nombre} mata una célula sana")
        if ag.propaga:
            # Replicación lítica: la progenie salta SOLO a vecinas tumorales en contacto.
            vecinas = [w for w in self._vecinos(q)
                       if (u := self.celulas.get(w)) is not None and u.tipo == LUMINAL
                       and u.n_drivers > 0 and not u.necrosis]
            self.rng.shuffle(vecinas)
            for w in vecinas[:int(ag.propaga)]:
                if len(vivos) + 1 >= MAX_AGENTES:
                    break
                vivos.append(replace(ag, pos=w, memoria={}, agotamiento=0.0, muertes=0,
                                     reparaciones=0, vida=ag.vida))
                self.stats["propagaciones"] += 1
        if ag.expansion and not t_natural and len(vivos) + 1 < MAX_AGENTES and \
                self.rng.random() < ag.expansion:
            vivos.append(replace(ag, memoria={}, agotamiento=0.0, muertes=0, reparaciones=0))

    # ------------------------------------------------- inmunidad adaptativa (T15)
    # Ver biologia.Parametros.inmunidad_adaptativa para fuentes y SUPUESTOS.
    def _pdl1(self, c: Celula) -> bool:
        """PD-L1 en superficie: por mutación (CD274) o inducido por IFN-γ (adaptativo)."""
        return c.pdl1 or self.t < c.pdl1_hasta

    def _inducir_pdl1(self, q):
        """Resistencia adaptativa (Taube 2012; Ribas 2015): el IFN-γ del linfocito que
        ataca enciende PD-L1 en la diana y en las células epiteliales de alrededor."""
        r = self.p.radio_ifng
        hasta = self.t + int(self.p.dias_pdl1_inducido * PASOS_POR_DIA)
        y, x = q
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                c = self.celulas.get((y + dy, x + dx))
                if c is not None and c.tipo == LUMINAL and c.pdl1_hasta < hasta:
                    if c.pdl1_hasta < self.t:
                        self.stats["pdl1_inducidos"] = self.stats.get("pdl1_inducidos", 0) + 1
                    c.pdl1_hasta = hasta
        self._hito("pdl1_adaptativo", "El IFN-γ del ataque enciende PD-L1 alrededor "
                                      "(resistencia adaptativa)")

    def _estimular_t(self, ag: Agente):
        """Un encuentro con el antígeno. El primero arranca la expansión clonal; los
        siguientes, acumulados, agotan (estimulación crónica, Wherry 2011)."""
        p = self.p
        ag.estimulos += 1
        if not ag.estado_t:
            ag.estado_t, ag.programa = "efector", p.divisiones_programa
            self._hito("expansion", "Un linfocito T reconoce un neoantígeno y arranca su "
                                    "expansión clonal")
        elif ag.estado_t == "efector" and ag.estimulos >= p.umbral_agotamiento:
            ag.programa = 0
            if ag.derivado:
                self._a_terminal(ag)
            else:
                ag.estado_t = "tex_prog"
                ag.vida = max(ag.vida, 160)
                self.stats["tex_progenitores"] = self.stats.get("tex_progenitores", 0) + 1
                self._hito("tex", "Estimulación crónica: aparecen linfocitos T agotados "
                                  "progenitores (TCF1+)")
        elif ag.estado_t == "tex_prog" and self.rng.random() < p.prob_tex_terminal:
            self._a_terminal(ag)

    def _a_terminal(self, ag: Agente):
        ag.estado_t, ag.programa = "tex_term", 0
        ag.vida = min(ag.vida, int(self.p.dias_vida_tex_term * PASOS_POR_DIA))
        self.stats["tex_terminales"] = self.stats.get("tex_terminales", 0) + 1

    def _hueco_clonal(self, vivos: list) -> bool:
        return self._n_clonales < self.p.max_t_clonales and len(vivos) + 1 < MAX_AGENTES

    def _nace_clon(self, ag: Agente, vivos: list, **cambios) -> Agente:
        hija = replace(ag, memoria={}, agotamiento=0.0, muertes=0, reparaciones=0,
                       clonal=True, **cambios)
        vivos.append(hija)
        self._n_clonales += 1
        self.stats["t_clonales"] = self.stats.get("t_clonales", 0) + 1
        return hija

    def _contraer(self, ag: Agente):
        """Fin del programa: el 90-95 % de los efectores muere (Kaech 2002)."""
        if self.rng.random() < self.p.contraccion_t:
            ag.vida = min(ag.vida, self.rng.randint(1, int(self.p.dias_contraccion * PASOS_POR_DIA)))

    def _ciclo_t(self, ag: Agente, vivos: list):
        """Un paso de vida de un linfocito T natural ya activado."""
        p, rng = self.p, self.rng
        if ag.estado_t == "efector" and ag.programa > 0:
            if rng.random() < p.prob_division_t:
                ag.programa -= 1
                hija = None
                if self._hueco_clonal(vivos):
                    hija = self._nace_clon(ag, vivos, vida=160)
                if ag.programa == 0:
                    self._contraer(ag)
                    if hija is not None:
                        self._contraer(hija)
        elif ag.estado_t == "tex_prog":
            if self.t < self.pd1_hasta and rng.random() < p.prob_revigorizacion:
                # Anti-PD-1: el progenitor se autorrenueva y suelta progenie efectora
                # que vuelve a matar y acaba terminal (Siddiqui 2019; Miller 2019).
                if self._hueco_clonal(vivos):
                    self._nace_clon(ag, vivos, estado_t="efector",
                                    programa=p.divisiones_revigorizacion,
                                    estimulos=p.umbral_agotamiento // 2, derivado=True, vida=160)
                    self.stats["revigorizaciones"] = self.stats.get("revigorizaciones", 0) + 1
                    self._hito("revigor", "Anti-PD-1: los T agotados progenitores vuelven a "
                                          "dividirse (revigorización)")
            elif rng.random() < p.autorrenovacion_tex and self._hueco_clonal(vivos):
                self._nace_clon(ag, vivos)

    def _reparar(self, ag: Agente, c: Celula, q, gen: str):
        if gen not in ag.vectores or ag.carga <= 0 or c.tipo != LUMINAL:
            return
        ag.carga -= 1
        if not c.genoma.get(gen):
            self.stats["reparaciones_inutiles"] += 1
            return
        if self.rng.random() >= self.p.eficacia_reparacion:
            return
        tenia_p53_roto = not c.p53_ok
        c.genoma[gen] = 0
        c.actualizar_fenotipo()
        ag.reparaciones += 1
        self.stats["reparaciones"] += 1
        if c.n_drivers == 0:
            c.reparada = True
            self.stats["reparaciones_que_curan"] += 1
            # Chispa del efecto espectador: la célula corregida pasa a corregir vecinas.
            if self.p.alcance_correccion:
                c.corrector = True
            self._hito("reparada", f"{ag.nombre} devuelve una célula a la normalidad ({gen} corregido)")
            if c.madre and ag.reprograma:
                c.guardiana = ag.reprograma
                self.stats["guardianas"] += 1
                self._hito("reconducida", f"Reina curada y reprogramada como guardiana ({ag.reprograma})")
        elif tenia_p53_roto and c.p53_ok and self.rng.random() < 0.5:
            # p53 restaurado detecta el resto de daños y ordena la apoptosis.
            del self.celulas[q]
            self.stats["apoptosis_p53"] += 1
            self.stats["reparaciones_p53_letales"] += 1
            self._hito("p53", "p53 restaurado: la célula tumoral reconoce su propio daño y se autodestruye")
        else:
            # Gen corregido pero la célula sigue siendo tumoral: la reparación no cambió nada.
            self.stats["reparaciones_esteriles"] += 1

    # ------------------------------------------------------------------ registro
    def _evento(self, texto: str):
        self.eventos.append(f"d{self.dia:6.1f}  {texto}")

    def _hito(self, clave: str, texto: str):
        if clave not in self._hitos:
            self._hitos.add(clave)
            self._evento(texto)

    def _registrar_clon(self, c: Celula, pos):
        firma = c.alterados
        if firma and firma not in self.clones:
            self.clones[firma] = len(self.clones) + 1
            self._evento(f"clon #{self.clones[firma]}: {' + '.join(sorted(firma))} "
                         f"→ {NOMBRE_ESTADO[self.estado(c, pos)]}")

    def _registrar(self):
        h = self.historia
        h["dia"].append(self.dia)
        for e in (*TUMORALES, REPARADA, NECROSIS):
            h[e].append(int(self.conteo[e]))
        for tipo in ("T", "NK", "IA"):
            h[tipo].append(sum(1 for a in self.agentes if a.tipo == tipo))

    def frecuencias_tumor(self) -> dict[str, float]:
        tumor = [c for c in self.celulas.values()
                 if c.tipo == LUMINAL and c.n_drivers > 0 and not c.necrosis]
        if not tumor:
            return {}
        return {g: sum(g in c.alterados for c in tumor) / len(tumor) for g in GENES}

    def describir(self, pos) -> list[str]:
        c = self.celulas.get(pos)
        if c is None:
            return [f"({pos[0]},{pos[1]}) {NOMBRE_ZONA[self.zona[pos]]}: vacío"
                    + (" · vaso sanguíneo" if self.vasos[pos] else ""), f"O2 {self.O2[pos]:.2f}"]
        golpes = ", ".join(f"{g}×{h}" for g, h in c.genoma.items() if h) or "ninguno"
        return [
            f"Célula #{c.id} · {NOMBRE_TEJIDO[c.tipo]} · {NOMBRE_ESTADO[self.estado(c, pos)]}",
            f"Zona {NOMBRE_ZONA[self.zona[pos]]} · O2 {self.O2[pos]:.2f} · ER{'+' if c.er else '-'}",
            f"Golpes: {golpes}",
            f"Alterados: {', '.join(sorted(c.alterados)) or '-'} · pasajeras {c.mutaciones}",
            f"p53 {'OK' if c.p53_ok else 'PERDIDO'} · HRD {'sí' if c.hrd else 'no'} · "
            f"MHC-I {'sí' if c.mhc else 'NO'} · PD-L1 {'alto' if c.pdl1 else 'bajo'} · "
            f"prolif ×{c.factor:.1f}",
            f"GUARDIANA (puerta {c.guardiana})" if c.guardiana else
            "Célula madre (reina)" if c.madre else
            f"Obrera · divisiones restantes {'∞' if c.divisiones is None else c.divisiones}"
            + (" · TELOMERASA BLOQUEADA" if not c.telomerasa_ok else "")
            + (" · SENESCENTE (jubilada)" if c.senescente else "")
            + (" · CEBO CARGADO" if c.cebo else ""),
        ]

    def resumen(self) -> dict:
        tumor = int(sum(self.conteo[e] for e in TUMORALES))
        self.stats["senescentes"] = sum(
            1 for c in self.celulas.values()
            if c.tipo == LUMINAL and c.n_drivers and (c.senescente or c.divisiones == 0))
        self.stats["aneuploides"] = sum(1 for c in self.celulas.values() if c.aneuploide)
        self.stats["aneuploides_sanas"] = sum(
            1 for c in self.celulas.values() if c.aneuploide and c.n_drivers == 0)
        return dict(dia=round(self.dia, 1), tumorales=tumor, in_situ=int(self.conteo[IN_SITU]),
                    peligrosas=int(self.conteo[IN_SITU] + self.conteo[INVASIVA] + self.conteo[REINA]),
                    reinas=int(self.conteo[REINA]),
                    invasivas=int(self.conteo[INVASIVA]), reparadas=int(self.conteo[REPARADA]),
                    normales=int(self.conteo[NORMAL]), **self.stats)
