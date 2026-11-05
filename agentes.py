"""Células inmunes y células 'inteligentes': qué perciben y qué pueden hacer.

Regla de oro del modelo: una célula NO puede leer el ADN de otra. Solo percibe
lo que la biología expone en la superficie (MHC-I con neoantígenos, PD-L1, HER2,
ligandos de estrés tipo MICA/B) y señales químicas (quimiocinas). Para conocer
el genoma tiene que Analizar (equivale a secuenciar), y eso le cuesta un turno.

Una política es una función  politica(agente, percepcion) -> acción.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Callable


@dataclass(frozen=True)
class Observacion:
    dy: int
    dx: int
    id: int
    tejido: str            # "luminal" | "mioepitelial" | "fibroblasto"
    mhc_i: bool            # presenta antígenos (False si ha perdido HLA)
    neoantigenos: int      # lo que un linfocito T puede reconocer (0 si no hay MHC-I)
    pdl1: bool
    her2_alto: bool
    ligandos_estres: bool  # MICA/B: aparecen con transformación o hipoxia
    proliferando: bool     # Ki-67
    hipoxica: bool
    marcador_madre: bool   # CD44+/CD24-/ALDH1+: también lo llevan las madres SANAS
    senescente: bool       # marcadores de senescencia (uPAR, DPP4): diana de los senolíticos
    cd47: bool             # "NO me comas": freno de la fagocitosis. Lo llevan TODAS las células
                           # sanas (los glóbulos rojos, los primeros). El tumor lo sube.
    comeme: bool           # "CÓMEME": calreticulina expuesta. Aparece con el estrés, el daño
                           # y sobre todo con la SENESCENCIA. Una célula sana no la expone.
    fuera_de_nicho: bool   # lleva marcador de madre pero NO está pegada a la membrana basal.
                           # Es posicional, no molecular: una vecina en contacto sí lo percibe
                           # (adhesión, polaridad, contacto con la lámina basal).


@dataclass
class Percepcion:
    vecinas: list[Observacion]           # células a distancia <= 2
    quimiotaxis: tuple[int, int] | None  # hacia dónde sube la señal inflamatoria
    hacia_conducto: tuple[int, int]
    zona: str                            # "luz" | "epitelio" | "mioepitelio" | "estroma"
    oxigeno: float
    dia: float
    azar: random.Random


@dataclass(frozen=True)
class Mover:
    dy: int
    dx: int


@dataclass(frozen=True)
class Atacar:
    dy: int
    dx: int


@dataclass(frozen=True)
class Analizar:
    dy: int
    dx: int


@dataclass(frozen=True)
class Reparar:
    dy: int
    dx: int
    gen: str


@dataclass(frozen=True)
class Esperar:
    pass


@dataclass
class Agente:
    tipo: str                      # "T" | "NK" | "IA"
    nombre: str
    politica: Callable
    receptores: frozenset          # "TCR", "NK", "HER2" (CAR), "INESPECIFICO"
    pos: tuple = (0, 0)
    vectores: frozenset = frozenset()  # genes que sabe corregir
    carga: int = 0                 # ediciones génicas que le quedan
    resistente_pdl1: bool = False  # PD-1 desactivado (como un CAR-T "blindado")
    expansion: float = 0.0         # probabilidad de dividirse tras matar (expansión clonal)
    vida: int = 200                # pasos
    agotamiento: float = 0.0
    memoria: dict = field(default_factory=dict)
    muertes: int = 0
    reparaciones: int = 0
    reprograma: str | None = None  # al curar una reina la convierte en guardiana ("simple" | "doble")
    propaga: float = 0.0           # al matar, copias que saltan a vecinas tumorales (virus oncolítico)
    ceba: bool = False             # deja cebos cargados en vez de limpiar células senescentes
    # CABALLO DE TROYA: no pelea contra la presión intersticial, entra invitado.
    # El tumor recluta monocitos activamente (CCL2) y les abre el endotelio: el vehículo
    # no nada contra corriente, lo llevan. Es la base de los CAR-M y de las células madre
    # mesenquimales usadas como vehículo. `troyano` = fracción de la corriente que ignora.
    troyano: float = 0.0
    entrega: str | None = None     # terapia en la que se convierte al llegar al tumor
    # BOMBONA DE OXÍGENO: el vehículo suelta O2 por donde pasa (perfluorocarbonos,
    # hemoglobina modificada, cianobacterias oxigenantes). Sirve para que la carga
    # sobreviva en el núcleo hipóxico... pero el O2 también alimenta al tumor.
    oxigeno: float = 0.0
    # CASTRAR en vez de matar: no elimina la célula, la jubila (telómero roto + p53).
    # Erradicar selecciona resistentes; esterilizar no. Y una célula jubilada expone
    # el "cómeme", así que los macrófagos la retiran después.
    castra: bool = False
    # INMUNIDAD ADAPTATIVA (biologia.Parametros.inmunidad_adaptativa). Solo la usan los
    # linfocitos T naturales; para el resto estos campos se quedan en su valor neutro.
    estado_t: str = ""             # "" sin activar · "efector" · "tex_prog" · "tex_term"
    programa: int = 0              # divisiones que le quedan del programa de expansión
    estimulos: int = 0             # encuentros con su antígeno (estimulación crónica)
    clonal: bool = False           # nacido por expansión (no cuenta para el reclutamiento)
    derivado: bool = False         # progenie de un Tex progenitor: acaba terminal, no progenitor


def distancia(o: Observacion) -> int:
    return max(abs(o.dy), abs(o.dx))


def adyacente(o: Observacion) -> bool:
    return distancia(o) <= 1


def hacia(dy: int, dx: int) -> Mover:
    return Mover((dy > 0) - (dy < 0), (dx > 0) - (dx < 0))


def deambular(p: Percepcion) -> Mover:
    return Mover(p.azar.randint(-1, 1), p.azar.randint(-1, 1))


def cazar(dianas: list[Observacion], p: Percepcion):
    """Ataca a la diana más cercana; si no hay, sigue las quimiocinas o deambula."""
    if dianas:
        o = min(dianas, key=distancia)
        return Atacar(o.dy, o.dx) if adyacente(o) else hacia(o.dy, o.dx)
    if p.quimiotaxis and p.azar.random() < 0.7:
        return hacia(*p.quimiotaxis)
    return deambular(p)


def politica_linfocito_t(ag: Agente, p: Percepcion):
    """Linfocito T CD8+: reconoce neoantígenos presentados por MHC-I.

    No reintenta una célula que acaba de frenarlo por PD-1 (ver Simulacion._atacar)."""
    frenado = ag.memoria.get("frenado_por", {})
    return cazar([o for o in p.vecinas
                  if o.tejido == "luminal" and o.mhc_i and o.neoantigenos > 0
                  and frenado.get(o.id, -1.0) <= p.dia], p)


def politica_macrofago(ag: Agente, p: Percepcion):
    """Macrófago: la puerta Y natural del cuerpo. No come por una señal, sino por dos.

    Necesita que falte el freno (CD47, "no me comas") Y que esté el acelerador
    (calreticulina, "cómeme"). Por eso bloquear CD47 no hace que se coma el cuerpo
    entero: a la célula sana le falta la segunda señal. Y por eso la senescencia
    importa tanto aquí: jubilar una célula es ponerle encima el "cómeme".
    """
    return cazar([o for o in p.vecinas if o.comeme and not o.cd47], p)


def politica_vehiculo(ag: Agente, p: Percepcion):
    """Monocito cargado: su único trabajo es LLEGAR. No ataca ni repara por el camino.

    Sigue las quimiocinas del tumor, que es justo lo que el tumor usa para reclutarlo.
    Al llegar al núcleo suelta la carga (ver Simulacion._entregar_troyano)."""
    if p.quimiotaxis:
        return hacia(*p.quimiotaxis)
    return deambular(p)


def politica_castrador(ag: Agente, p: Percepcion):
    """Soldado que esteriliza. No mata: deja a la célula sin poder dividirse.

    Objetivo declarado: que el tumor no pueda CRECER. Lo que queda, jubilado y marcado
    con el 'cómeme', lo retiran los macrófagos o se extirpa después."""
    return cazar([o for o in p.vecinas if o.tejido == "luminal"
                  and (o.neoantigenos > 0 or not o.mhc_i or o.ligandos_estres
                       or o.marcador_madre)], p)


def politica_nk(ag: Agente, p: Percepcion):
    """Célula NK: mata lo que ha perdido MHC-I ('falta de lo propio') o muestra estrés."""
    return cazar([o for o in p.vecinas
                  if o.tejido == "luminal" and (not o.mhc_i or o.ligandos_estres)], p)
