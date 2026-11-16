"""Células terapéuticas 'inteligentes'. Cada una = receptores + vectores + una política.

Todas usan la misma API que los linfocitos naturales (ver agents.py): solo ven
la superficie de las células y, si quieren conocer el genoma, deben Analyze.
"""

from __future__ import annotations

from breast_cancer_abm.model import my_cell
from breast_cancer_abm.model.agents import (
    Agent,
    Analyze,
    Attack,
    Repair,
    adjacent,
    distance,
    hunt,
    sterilizer_policy,
    toward,
    vehicle_policy,
    wander,
)
from breast_cancer_abm.model.biology import GENES

# Orden en que conviene corregir: primero el guardián (p53 restaurado puede
# eliminar él solo a la célula), luego la reparación del ADN, luego los aceleradores.
PRIORITY = ("TP53", "BRCA1", "BRCA2", "PTEN", "CDH1", "PIK3CA", "ERBB2", "MYC")


def _suspicious(o) -> bool:
    return o.tissue == "luminal" and (o.stress_ligands or o.neoantigens > 0 or o.proliferating)


def hunter(ag: Agent, p):
    """CAR-T 'blindado': reconoce como un T y como una NK, y es inmune a PD-L1."""
    return hunt(
        [
            o
            for o in p.neighbors
            if o.tissue == "luminal" and ((o.mhc_i and o.neoantigens > 0) or not o.mhc_i or o.stress_ligands)
        ],
        p,
    )


def car_her2(ag: Agent, p):
    """CAR-T anti-HER2: solo ve HER2. Inútil contra clones sin ERBB2.

    Es el único agente de terapia al que PD-L1 puede frenar (no lleva pdl1_resistant),
    así que es el único que necesita el desenganche: sin él volvía a la misma célula
    PD-L1+ que acababa de frenarlo (ver Simulation._attack)."""
    blocked = ag.memory.get("blocked_by", {})
    return hunt([o for o in p.neighbors if o.her2_high and blocked.get(o.id, -1.0) <= p.day], p)


def _repair_or_attack(ag: Agent, p, attack_threshold: int | None, include_heterozygous: bool):
    genomes = ag.memory.setdefault("genomes", {})
    candidate_cells = [
        o for o in p.neighbors if o.tissue == "luminal" and (include_heterozygous or _suspicious(o))
    ]
    for o in sorted(candidate_cells, key=distance):
        known = genomes.get(o.id)
        if known == {}:
            continue  # ya analizada y sana
        if not adjacent(o):
            return toward(o.dy, o.dx)
        if known is None:
            return Analyze(o.dy, o.dx)  # secuenciar antes de tocar nada
        broken = [g for g in PRIORITY if g in known and g in ag.vectors]
        if attack_threshold is not None and len(known) >= attack_threshold:
            del genomes[o.id]
            return Attack(o.dy, o.dx)  # demasiados fallos: ya no compensa reparar
        if broken and ag.payload > 0:
            del genomes[o.id]  # tras editar hay que volver a comprobar
            return Repair(o.dy, o.dx, broken[0])
        genomes[o.id] = {}
    if p.chemotaxis and p.rand.random() < 0.6:
        return toward(*p.chemotaxis)
    if p.zone == "stroma":
        return toward(*p.toward_duct)
    return wander(p)


def repairer(ag: Agent, p):
    """Terapia génica pura: analiza y corrige, nunca mata."""
    return _repair_or_attack(ag, p, attack_threshold=None, include_heterozygous=False)


def hybrid(ag: Agent, p):
    """Repara lo que aún tiene arreglo (< 3 genes) y elimina lo que ya no."""
    return _repair_or_attack(ag, p, attack_threshold=3, include_heterozygous=False)


def sentinel(ag: Agent, p):
    """Prevención: patrulla el epitelio y corrige cualquier golpe, incluso el alelo heredado."""
    if p.zone in ("stroma", "lumen") and p.rand.random() < 0.8:
        act = toward(*p.toward_duct)
        if p.zone == "lumen":  # en la luz hay que ir hacia fuera, al epitelio
            act = toward(-act.dy, -act.dx)
        return act
    return _repair_or_attack(ag, p, attack_threshold=None, include_heterozygous=True)


def chemo(ag: Agent, p):
    """No es una célula: modela un citotóxico que mata lo que prolifera, sea tumor o no."""
    targets = [o for o in p.neighbors if o.proliferating]
    if targets:
        o = min(targets, key=distance)
        return Attack(o.dy, o.dx) if adjacent(o) else toward(o.dy, o.dx)
    return toward(*p.toward_duct) if p.zone == "stroma" else wander(p)


def _is_tumor(genome: dict) -> bool:
    """Un oncogén activado o un supresor con las dos copias rotas."""
    return any(GENES[g].kind == "oncogene" or h >= 2 for g, h in genome.items())


def _queen_strategy(ag: Agent, p, repair: bool, clean: bool):
    """Busca reinas (marcador de célula madre), las secuencia y las arregla o las mata.
    Si no ve ninguna, opcionalmente limpia las obreras corruptas."""
    genomes = ag.memory.setdefault("genomes", {})
    queens = [o for o in p.neighbors if o.tissue == "luminal" and o.stem_marker and genomes.get(o.id) != {}]
    for o in sorted(queens, key=distance):
        if not adjacent(o):
            return toward(o.dy, o.dx)
        known = genomes.get(o.id)
        if known is None:
            return Analyze(o.dy, o.dx)
        broken = [g for g in PRIORITY if g in known and g in ag.vectors]
        if not _is_tumor(known) or (repair and not broken):
            genomes[o.id] = {}  # madre sana: no se toca
            continue
        del genomes[o.id]  # tras actuar, volver a comprobar
        if repair and ag.payload > 0:
            return Repair(o.dy, o.dx, broken[0])
        return Attack(o.dy, o.dx)
    if clean:
        corrupted = [
            o
            for o in p.neighbors
            if o.tissue == "luminal" and not o.stem_marker and (o.stress_ligands or o.neoantigens > 0)
        ]
        if corrupted:
            return hunt(corrupted, p)
    if p.chemotaxis and p.rand.random() < 0.6:
        return toward(*p.chemotaxis)
    if p.zone == "stroma":
        return toward(*p.toward_duct)
    return wander(p)


def propagator(ag: Agent, p):
    """Virus oncolítico: solo se replica dentro de células tumorales y salta a las vecinas.
    Se extingue solo cuando no queda tumor alrededor (no hace falta tratar célula a célula)."""
    targets = [
        o
        for o in p.neighbors
        if o.tissue == "luminal"
        and (o.neoantigens > 0 or not o.mhc_i)
        and (o.proliferating or o.stress_ligands)
    ]
    if targets:
        o = min(targets, key=distance)
        return Attack(o.dy, o.dx) if adjacent(o) else toward(o.dy, o.dx)
    if p.chemotaxis and p.rand.random() < 0.6:
        return toward(*p.chemotaxis)
    return wander(p)


def queen(ag: Agent, p):
    """Arregla el bug en la reina (célula madre tumoral) y extirpa las obreras corruptas."""
    return _queen_strategy(ag, p, repair=True, clean=True)


def queen_only(ag: Agent, p):
    """Solo arregla las reinas; espera a que las obreras corruptas se agoten solas."""
    return _queen_strategy(ag, p, repair=True, clean=False)


def queen_killer(ag: Agent, p):
    """No repara: mata a las reinas tumorales (tras secuenciarlas) y limpia obreras."""
    return _queen_strategy(ag, p, repair=False, clean=True)


ALL_GENES = frozenset(PRIORITY)

THERAPIES = {
    "hunter": dict(
        n=40,
        description=hunter.__doc__,
        create=lambda: Agent(
            "AI", "hunter", hunter, frozenset({"TCR", "NK"}), pdl1_resistant=True, expansion=0.25, life=300
        ),
    ),
    "car_her2": dict(
        n=40,
        description=car_her2.__doc__,
        create=lambda: Agent("AI", "CAR-HER2", car_her2, frozenset({"HER2"}), expansion=0.25, life=300),
    ),
    "repairer": dict(
        n=40,
        description=repairer.__doc__,
        create=lambda: Agent(
            "AI", "repairer", repairer, frozenset(), vectors=ALL_GENES, payload=30, life=400
        ),
    ),
    "hybrid": dict(
        n=40,
        description=hybrid.__doc__,
        create=lambda: Agent(
            "AI",
            "híbrido",
            hybrid,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=30,
            pdl1_resistant=True,
            expansion=0.15,
            life=400,
        ),
    ),
    "sentinel": dict(
        n=30,
        description=sentinel.__doc__,
        create=lambda: Agent(
            "AI", "sentinel", sentinel, frozenset(), vectors=ALL_GENES, payload=80, life=800
        ),
    ),
    "queen": dict(
        n=40,
        description=queen.__doc__,
        create=lambda: Agent(
            "AI",
            "queen",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=60,
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
        ),
    ),
    "queen_only": dict(
        n=40,
        description=queen_only.__doc__,
        create=lambda: Agent(
            "AI", "reina (solo reparar)", queen_only, frozenset(), vectors=ALL_GENES, payload=60, life=500
        ),
    ),
    "queen_killer": dict(
        n=40,
        description=queen_killer.__doc__,
        create=lambda: Agent(
            "AI",
            "queen-killer",
            queen_killer,
            frozenset({"NONSPECIFIC"}),
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
        ),
    ),
    "guardian": dict(
        n=40,
        description="Cura reinas y las convierte en guardianas con puerta Y",
        create=lambda: Agent(
            "AI",
            "reconductor (Y)",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=60,
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
            reprograms="double",
        ),
    ),
    "simple_guardian": dict(
        n=40,
        description="Igual, pero la guardiana ataca con una sola señal",
        create=lambda: Agent(
            "AI",
            "reconductor (1 señal)",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=60,
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
            reprograms="simple",
        ),
    ),
    "niche_guardian": dict(
        n=40,
        description="Guardiana con puerta Y posicional: marcador de madre Y fuera de la capa basal",
        create=lambda: Agent(
            "AI",
            "reconductor (Y nicho)",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=60,
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
            reprograms="niche",
        ),
    ),
    "tert_guardian": dict(
        n=40,
        description="Guardiana con puerta Y interna: telomerasa activa Y aneuploidía",
        create=lambda: Agent(
            "AI",
            "reconductor (Y TERT+aneu)",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=60,
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
            reprograms="tert",
        ),
    ),
    "tert_niche_guardian": dict(
        n=40,
        description="Guardiana con las dos puertas Y en paralelo (O de dos Y)",
        create=lambda: Agent(
            "AI",
            "reconductor (Y TERT ó Y nicho)",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=60,
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
            reprograms="tert+niche",
        ),
    ),
    "propagator": dict(
        n=4,
        description=propagator.__doc__,
        create=lambda: Agent(
            "AI",
            "propagator",
            propagator,
            frozenset({"TCR", "NK"}),
            pdl1_resistant=True,
            life=120,
            spreads=2.0,
        ),
    ),
    "bait": dict(
        n=40,
        description="Guardiana que deja cebos cargados en las senescentes",
        create=lambda: Agent(
            "AI",
            "bait",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=60,
            pdl1_resistant=True,
            expansion=0.15,
            life=500,
            reprograms="bait",
        ),
    ),
    "tert_memory": dict(
        n=15,
        description="Vigilancia residente con la puerta Y interna (TERT + aneuploidía)",
        create=lambda: Agent(
            "AI",
            "memoria (Y TERT+aneu)",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=200,
            pdl1_resistant=True,
            life=10_000,
            reprograms="tert",
        ),
    ),
    "tert_niche_memory": dict(
        n=15,
        description="Vigilancia residente con las dos puertas Y en paralelo",
        create=lambda: Agent(
            "AI",
            "memoria (Y TERT ó Y nicho)",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=200,
            pdl1_resistant=True,
            life=10_000,
            reprograms="tert+niche",
        ),
    ),
    "memory": dict(
        n=15,
        description="Vigilancia residente de larga vida (como los linfocitos T de memoria)",
        create=lambda: Agent(
            "AI",
            "memory",
            queen,
            frozenset({"TCR", "NK"}),
            vectors=ALL_GENES,
            payload=200,
            pdl1_resistant=True,
            life=10_000,
            reprograms="double",
        ),
    ),
    # ---- CABALLO DE TROYA -------------------------------------------------------
    # No pelean contra la presión: entran invitados. El tumor recluta monocitos con CCL2
    # y les abre el endotelio; el vehículo no nada contra corriente, lo llevan. Al llegar
    # al interior suelta la terapia real, que ya nace dentro y no tiene que atravesar nada.
    # Real: CAR-M (macrófagos CAR, en ensayo) y células madre mesenquimales como vehículo.
    "trojan": dict(
        n=40,
        description="Monocito cargado: entra invitado y suelta la guardiana dentro",
        create=lambda: Agent(
            "AI",
            "vehículo troyano",
            vehicle_policy,
            frozenset(),
            life=400,
            trojan=0.9,
            delivers="niche_guardian",
        ),
    ),
    "trojan_bait": dict(
        n=40,
        description="Vehículo que suelta cebos: deja que el tumor se lo coma",
        create=lambda: Agent(
            "AI", "troyano (cebo)", vehicle_policy, frozenset(), life=400, trojan=0.9, delivers="bait"
        ),
    ),
    "trojan_memory": dict(
        n=20,
        description="Vehículo que deja vigilancia residente dentro del tumor",
        create=lambda: Agent(
            "AI", "troyano (memoria)", vehicle_policy, frozenset(), life=400, trojan=0.9, delivers="memory"
        ),
    ),
    # Soldados que ESTERILIZAN en vez de matar. El objetivo no es vaciar el tumor, es
    # que no pueda crecer: lo que queda, jubilado, lo retiran los macrófagos o se extirpa.
    "sterilizer": dict(
        n=60,
        description=sterilizer_policy.__doc__,
        create=lambda: Agent(
            "AI",
            "sterilizer",
            sterilizer_policy,
            frozenset({"TCR", "NK"}),
            pdl1_resistant=True,
            expansion=0.2,
            life=400,
            sterilizes=True,
        ),
    ),
    # El caballo de Troya completo: entra invitado, lleva su propio oxígeno para que la
    # carga no se asfixie en el núcleo, y suelta castradores dentro.
    "trojan_horse": dict(
        n=120,
        description="Vehículo con bombona de O2 que suelta castradores dentro",
        create=lambda: Agent(
            "AI",
            "troya (O2 + castradores)",
            vehicle_policy,
            frozenset(),
            life=400,
            trojan=0.9,
            delivers="sterilizer",
            oxygen=0.25,
        ),
    ),
    "trojan_horse_no_o2": dict(
        n=120,
        description="Igual pero sin bombona: mide cuánto aporta el O2",
        create=lambda: Agent(
            "AI", "troya (sin O2)", vehicle_policy, frozenset(), life=400, trojan=0.9, delivers="sterilizer"
        ),
    ),
    "chemo": dict(
        n=80,
        description=chemo.__doc__,
        create=lambda: Agent("AI", "chemo", chemo, frozenset({"NONSPECIFIC"}), life=40),
    ),
    "my_cell": dict(n=my_cell.AMOUNT, description=my_cell.my_policy.__doc__, create=my_cell.create),
}
