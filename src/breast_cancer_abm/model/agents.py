"""Células inmunes y células 'inteligentes': qué perciben y qué pueden hacer.

Regla de oro del modelo: una célula NO puede leer el ADN de otra. Solo percibe
lo que la biología expone en la superficie (MHC-I con neoantígenos, PD-L1, HER2,
ligandos de estrés tipo MICA/B) y señales químicas (quimiocinas). Para conocer
el genoma tiene que Analyze (equivale a secuenciar), y eso le cuesta un turno.

Una política es una función  policy(agente, percepcion) -> acción.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from typing import Callable


@dataclass(frozen=True)
class Observation:
    dy: int
    dx: int
    id: int
    tissue: str  # "luminal" | "mioepitelial" | "fibroblasto"
    mhc_i: bool  # presenta antígenos (False si ha perdido HLA)
    neoantigens: int  # lo que un linfocito T puede reconocer (0 si no hay MHC-I)
    pdl1: bool
    her2_high: bool
    stress_ligands: bool  # MICA/B: aparecen con transformación o hipoxia
    proliferating: bool  # Ki-67
    hypoxic: bool
    stem_marker: bool  # CD44+/CD24-/ALDH1+: también lo llevan las madres SANAS
    senescent: bool  # marcadores de senescencia (uPAR, DPP4): diana de los senolíticos
    cd47: bool  # "NO me comas": freno de la fagocitosis. Lo llevan TODAS las células
    # sanas (los glóbulos rojos, los primeros). El tumor lo sube.
    eat_me: bool  # "CÓMEME": calreticulina expuesta. Aparece con el estrés, el daño
    # y sobre todo con la SENESCENCIA. Una célula sana no la expone.
    off_niche: bool  # lleva marcador de madre pero NO está pegada a la membrana basal.
    # Es posicional, no molecular: una vecina en contacto sí lo percibe
    # (adhesión, polaridad, contacto con la lámina basal).


@dataclass
class Perception:
    neighbors: list[Observation]  # células a distancia <= 2
    chemotaxis: tuple[int, int] | None  # hacia dónde sube la señal inflamatoria
    toward_duct: tuple[int, int]
    zone: str  # "luz" | "epitelio" | "mioepitelio" | "estroma"
    oxygen: float
    day: float
    rand: random.Random


@dataclass(frozen=True)
class Move:
    dy: int
    dx: int


@dataclass(frozen=True)
class Attack:
    dy: int
    dx: int


@dataclass(frozen=True)
class Analyze:
    dy: int
    dx: int


@dataclass(frozen=True)
class Repair:
    dy: int
    dx: int
    gene: str


@dataclass(frozen=True)
class Wait:
    pass


@dataclass
class Agent:
    kind: str  # "T" | "NK" | "IA"
    name: str
    policy: Callable
    receptors: frozenset  # "TCR", "NK", "HER2" (CAR), "INESPECIFICO"
    pos: tuple = (0, 0)
    vectors: frozenset = frozenset()  # genes que sabe corregir
    payload: int = 0  # ediciones génicas que le quedan
    pdl1_resistant: bool = False  # PD-1 desactivado (como un CAR-T "blindado")
    expansion: float = 0.0  # probabilidad de dividirse tras matar (expansión clonal)
    life: int = 200  # pasos
    exhaustion: float = 0.0
    memory: dict = field(default_factory=dict)
    kills: int = 0
    repairs: int = 0
    reprograms: str | None = None  # al curar una reina la convierte en guardiana ("simple" | "doble")
    spreads: float = 0.0  # al matar, copias que saltan a vecinas tumorales (virus oncolítico)
    leaves_bait: bool = False  # deja cebos cargados en vez de limpiar células senescentes
    # CABALLO DE TROYA: no pelea contra la presión intersticial, entra invitado.
    # El tumor recluta monocitos activamente (CCL2) y les abre el endotelio: el vehículo
    # no nada contra corriente, lo llevan. Es la base de los CAR-M y de las células madre
    # mesenquimales usadas como vehículo. `troyano` = fracción de la corriente que ignora.
    trojan: float = 0.0
    delivers: str | None = None  # terapia en la que se convierte al llegar al tumor
    # BOMBONA DE OXÍGENO: el vehículo suelta O2 por donde pasa (perfluorocarbonos,
    # hemoglobina modificada, cianobacterias oxigenantes). Sirve para que la carga
    # sobreviva en el núcleo hipóxico... pero el O2 también alimenta al tumor.
    oxygen: float = 0.0
    # STERILIZE en vez de matar: no elimina la célula, la jubila (telómero roto + p53).
    # Erradicar selecciona resistentes; esterilizar no. Y una célula jubilada expone
    # el "cómeme", así que los macrófagos la retiran después.
    sterilizes: bool = False
    # INMUNIDAD ADAPTATIVA (biologia.Parameters.adaptive_immunity). Solo la usan los
    # linfocitos T naturales; para el resto estos campos se quedan en su valor neutro.
    t_state: str = ""  # "" sin activar · "efector" · "tex_prog" · "tex_term"
    program: int = 0  # divisiones que le quedan del programa de expansión
    stimuli: int = 0  # encuentros con su antígeno (estimulación crónica)
    clonal: bool = False  # nacido por expansión (no cuenta para el reclutamiento)
    derived: bool = False  # progenie de un Tex progenitor: acaba terminal, no progenitor


def distance(o: Observation) -> int:
    return max(abs(o.dy), abs(o.dx))


def adjacent(o: Observation) -> bool:
    return distance(o) <= 1


def toward(dy: int, dx: int) -> Move:
    return Move((dy > 0) - (dy < 0), (dx > 0) - (dx < 0))


def wander(p: Perception) -> Move:
    return Move(p.rand.randint(-1, 1), p.rand.randint(-1, 1))


def hunt(targets: list[Observation], p: Perception):
    """Ataca a la diana más cercana; si no hay, sigue las quimiocinas o deambula."""
    if targets:
        o = min(targets, key=distance)
        return Attack(o.dy, o.dx) if adjacent(o) else toward(o.dy, o.dx)
    if p.chemotaxis and p.rand.random() < 0.7:
        return toward(*p.chemotaxis)
    return wander(p)


def t_cell_policy(ag: Agent, p: Perception):
    """Linfocito T CD8+: reconoce neoantígenos presentados por MHC-I.

    No reintenta una célula que acaba de frenarlo por PD-1 (ver Simulation._attack)."""
    blocked = ag.memory.get("blocked_by", {})
    return hunt(
        [
            o
            for o in p.neighbors
            if o.tissue == "luminal" and o.mhc_i and o.neoantigens > 0 and blocked.get(o.id, -1.0) <= p.day
        ],
        p,
    )


def macrophage_policy(ag: Agent, p: Perception):
    """Macrófago: la puerta Y natural del cuerpo. No come por una señal, sino por dos.

    Necesita que falte el freno (CD47, "no me comas") Y que esté el acelerador
    (calreticulina, "cómeme"). Por eso bloquear CD47 no hace que se coma el cuerpo
    entero: a la célula sana le falta la segunda señal. Y por eso la senescencia
    importa tanto aquí: jubilar una célula es ponerle encima el "cómeme".
    """
    return hunt([o for o in p.neighbors if o.eat_me and not o.cd47], p)


def vehicle_policy(ag: Agent, p: Perception):
    """Monocito cargado: su único trabajo es LLEGAR. No ataca ni repara por el camino.

    Sigue las quimiocinas del tumor, que es justo lo que el tumor usa para reclutarlo.
    Al llegar al núcleo suelta la carga (ver Simulation._deliver_trojan)."""
    if p.chemotaxis:
        return toward(*p.chemotaxis)
    return wander(p)


def sterilizer_policy(ag: Agent, p: Perception):
    """Soldado que esteriliza. No mata: deja a la célula sin poder dividirse.

    Objetivo declarado: que el tumor no pueda CRECER. Lo que queda, jubilado y marcado
    con el 'cómeme', lo retiran los macrófagos o se extirpa después."""
    return hunt(
        [
            o
            for o in p.neighbors
            if o.tissue == "luminal"
            and (o.neoantigens > 0 or not o.mhc_i or o.stress_ligands or o.stem_marker)
        ],
        p,
    )


def nk_policy(ag: Agent, p: Perception):
    """Célula NK: mata lo que ha perdido MHC-I ('falta de lo propio') o muestra estrés."""
    return hunt([o for o in p.neighbors if o.tissue == "luminal" and (not o.mhc_i or o.stress_ligands)], p)
