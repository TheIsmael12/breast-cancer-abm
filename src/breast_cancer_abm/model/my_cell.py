"""PLANTILLA: programa aquí tu propia célula y lánzala con la opción 'my_cell'.

Cada paso (6 h) el simulador llama a my_policy(agente, percepcion) y ejecuta la
acción que devuelvas:

    Move(dy, dx)          desplazarse una casilla (dy, dx en -1..1)
    Attack(dy, dx)         intentar matar a una célula adyacente (según tus receptores)
    Analyze(dy, dx)       secuenciarla: guarda {gen: golpes} en agente.memoria["genomas"][id]
    Repair(dy, dx, gen)   editar un gen (debe estar en tus vectores y gasta carga)
    Wait()

percepcion.vecinas trae Observation(dy, dx, id, tejido, mhc_i, neoantigenos, pdl1,
her2_high, stress_ligands, proliferando, hipoxica, stem_marker, senescente). Además:
quimiotaxis, toward_duct, zona, oxigeno, dia y azar (generador aleatorio).
"""

from breast_cancer_abm.model.agents import Agent, Analyze, Attack, Repair, adjacent, distance, hunt, toward

AMOUNT = 30


def my_policy(ag: Agent, p):
    """Ejemplo editable: ataca lo que tiene PD-L1 y repara TP53 en el resto."""
    genomes = ag.memory.setdefault("genomes", {})
    for o in sorted((o for o in p.neighbors if o.tissue == "luminal" and o.stress_ligands), key=distance):
        if not adjacent(o):
            return toward(o.dy, o.dx)
        if o.pdl1:
            return Attack(o.dy, o.dx)
        if o.id not in genomes:
            return Analyze(o.dy, o.dx)
        if genomes.pop(o.id).get("TP53") and ag.payload > 0:
            return Repair(o.dy, o.dx, "TP53")
    return hunt([], p)


def create() -> Agent:
    return Agent(
        "AI",
        "mi célula",
        my_policy,
        receptors=frozenset({"TCR"}),  # "TCR", "NK", "HER2", "INESPECIFICO"
        vectors=frozenset({"TP53"}),
        payload=20,
        pdl1_resistant=True,
        expansion=0.1,
        life=300,
    )
