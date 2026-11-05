"""PLANTILLA: programa aquí tu propia célula y lánzala con la opción 'mi_celula'.

Cada paso (6 h) el simulador llama a mi_politica(agente, percepcion) y ejecuta la
acción que devuelvas:

    Mover(dy, dx)          desplazarse una casilla (dy, dx en -1..1)
    Atacar(dy, dx)         intentar matar a una célula adyacente (según tus receptores)
    Analizar(dy, dx)       secuenciarla: guarda {gen: golpes} en agente.memoria["genomas"][id]
    Reparar(dy, dx, gen)   editar un gen (debe estar en tus vectores y gasta carga)
    Esperar()

percepcion.vecinas trae Observacion(dy, dx, id, tejido, mhc_i, neoantigenos, pdl1,
her2_alto, ligandos_estres, proliferando, hipoxica, marcador_madre, senescente). Además:
quimiotaxis, hacia_conducto, zona, oxigeno, dia y azar (generador aleatorio).
"""
from agentes import Agente, Analizar, Atacar, Reparar, adyacente, cazar, distancia, hacia

CANTIDAD = 30


def mi_politica(ag: Agente, p):
    """Ejemplo editable: ataca lo que tiene PD-L1 y repara TP53 en el resto."""
    genomas = ag.memoria.setdefault("genomas", {})
    for o in sorted((o for o in p.vecinas if o.tejido == "luminal" and o.ligandos_estres),
                    key=distancia):
        if not adyacente(o):
            return hacia(o.dy, o.dx)
        if o.pdl1:
            return Atacar(o.dy, o.dx)
        if o.id not in genomas:
            return Analizar(o.dy, o.dx)
        if genomas.pop(o.id).get("TP53") and ag.carga > 0:
            return Reparar(o.dy, o.dx, "TP53")
    return cazar([], p)


def crear() -> Agente:
    return Agente("IA", "mi célula", mi_politica,
                  receptores=frozenset({"TCR"}),   # "TCR", "NK", "HER2", "INESPECIFICO"
                  vectores=frozenset({"TP53"}),
                  carga=20, resistente_pdl1=True, expansion=0.1, vida=300)
