"""Células terapéuticas 'inteligentes'. Cada una = receptores + vectores + una política.

Todas usan la misma API que los linfocitos naturales (ver agentes.py): solo ven
la superficie de las células y, si quieren conocer el genoma, deben Analizar.
"""
from __future__ import annotations

from agentes import (Agente, Analizar, Atacar, Reparar, adyacente, cazar, deambular, distancia,
                     hacia, politica_castrador, politica_vehiculo)
from biologia import GENES
import mi_celula

# Orden en que conviene corregir: primero el guardián (p53 restaurado puede
# eliminar él solo a la célula), luego la reparación del ADN, luego los aceleradores.
PRIORIDAD = ("TP53", "BRCA1", "BRCA2", "PTEN", "CDH1", "PIK3CA", "ERBB2", "MYC")


def _sospechosa(o) -> bool:
    return o.tejido == "luminal" and (o.ligandos_estres or o.neoantigenos > 0 or o.proliferando)


def cazador(ag: Agente, p):
    """CAR-T 'blindado': reconoce como un T y como una NK, y es inmune a PD-L1."""
    return cazar([o for o in p.vecinas if o.tejido == "luminal" and
                  ((o.mhc_i and o.neoantigenos > 0) or not o.mhc_i or o.ligandos_estres)], p)


def car_her2(ag: Agente, p):
    """CAR-T anti-HER2: solo ve HER2. Inútil contra clones sin ERBB2.

    Es el único agente de terapia al que PD-L1 puede frenar (no lleva resistente_pdl1),
    así que es el único que necesita el desenganche: sin él volvía a la misma célula
    PD-L1+ que acababa de frenarlo (ver Simulacion._atacar)."""
    frenado = ag.memoria.get("frenado_por", {})
    return cazar([o for o in p.vecinas
                  if o.her2_alto and frenado.get(o.id, -1.0) <= p.dia], p)


def _reparar_o_atacar(ag: Agente, p, umbral_atacar: int | None, incluir_heterocigotos: bool):
    genomas = ag.memoria.setdefault("genomas", {})
    candidatas = [o for o in p.vecinas if o.tejido == "luminal" and
                  (incluir_heterocigotos or _sospechosa(o))]
    for o in sorted(candidatas, key=distancia):
        conocido = genomas.get(o.id)
        if conocido == {}:
            continue                       # ya analizada y sana
        if not adyacente(o):
            return hacia(o.dy, o.dx)
        if conocido is None:
            return Analizar(o.dy, o.dx)    # secuenciar antes de tocar nada
        rotos = [g for g in PRIORIDAD if g in conocido and g in ag.vectores]
        if umbral_atacar is not None and len(conocido) >= umbral_atacar:
            del genomas[o.id]
            return Atacar(o.dy, o.dx)      # demasiados fallos: ya no compensa reparar
        if rotos and ag.carga > 0:
            del genomas[o.id]              # tras editar hay que volver a comprobar
            return Reparar(o.dy, o.dx, rotos[0])
        genomas[o.id] = {}
    if p.quimiotaxis and p.azar.random() < 0.6:
        return hacia(*p.quimiotaxis)
    if p.zona == "estroma":
        return hacia(*p.hacia_conducto)
    return deambular(p)


def reparador(ag: Agente, p):
    """Terapia génica pura: analiza y corrige, nunca mata."""
    return _reparar_o_atacar(ag, p, umbral_atacar=None, incluir_heterocigotos=False)


def hibrido(ag: Agente, p):
    """Repara lo que aún tiene arreglo (< 3 genes) y elimina lo que ya no."""
    return _reparar_o_atacar(ag, p, umbral_atacar=3, incluir_heterocigotos=False)


def centinela(ag: Agente, p):
    """Prevención: patrulla el epitelio y corrige cualquier golpe, incluso el alelo heredado."""
    if p.zona in ("estroma", "luz") and p.azar.random() < 0.8:
        accion = hacia(*p.hacia_conducto)
        if p.zona == "luz":   # en la luz hay que ir hacia fuera, al epitelio
            accion = hacia(-accion.dy, -accion.dx)
        return accion
    return _reparar_o_atacar(ag, p, umbral_atacar=None, incluir_heterocigotos=True)


def quimio(ag: Agente, p):
    """No es una célula: modela un citotóxico que mata lo que prolifera, sea tumor o no."""
    dianas = [o for o in p.vecinas if o.proliferando]
    if dianas:
        o = min(dianas, key=distancia)
        return Atacar(o.dy, o.dx) if adyacente(o) else hacia(o.dy, o.dx)
    return hacia(*p.hacia_conducto) if p.zona == "estroma" else deambular(p)


def _es_tumoral(genoma: dict) -> bool:
    """Un oncogén activado o un supresor con las dos copias rotas."""
    return any(GENES[g].tipo == "oncogen" or h >= 2 for g, h in genoma.items())


def _estrategia_reina(ag: Agente, p, reparar: bool, limpiar: bool):
    """Busca reinas (marcador de célula madre), las secuencia y las arregla o las mata.
    Si no ve ninguna, opcionalmente limpia las obreras corruptas."""
    genomas = ag.memoria.setdefault("genomas", {})
    reinas = [o for o in p.vecinas
              if o.tejido == "luminal" and o.marcador_madre and genomas.get(o.id) != {}]
    for o in sorted(reinas, key=distancia):
        if not adyacente(o):
            return hacia(o.dy, o.dx)
        conocido = genomas.get(o.id)
        if conocido is None:
            return Analizar(o.dy, o.dx)
        rotos = [g for g in PRIORIDAD if g in conocido and g in ag.vectores]
        if not _es_tumoral(conocido) or (reparar and not rotos):
            genomas[o.id] = {}             # madre sana: no se toca
            continue
        del genomas[o.id]                  # tras actuar, volver a comprobar
        if reparar and ag.carga > 0:
            return Reparar(o.dy, o.dx, rotos[0])
        return Atacar(o.dy, o.dx)
    if limpiar:
        corruptas = [o for o in p.vecinas if o.tejido == "luminal" and not o.marcador_madre
                     and (o.ligandos_estres or o.neoantigenos > 0)]
        if corruptas:
            return cazar(corruptas, p)
    if p.quimiotaxis and p.azar.random() < 0.6:
        return hacia(*p.quimiotaxis)
    if p.zona == "estroma":
        return hacia(*p.hacia_conducto)
    return deambular(p)


def propagador(ag: Agente, p):
    """Virus oncolítico: solo se replica dentro de células tumorales y salta a las vecinas.
    Se extingue solo cuando no queda tumor alrededor (no hace falta tratar célula a célula)."""
    dianas = [o for o in p.vecinas if o.tejido == "luminal"
              and (o.neoantigenos > 0 or not o.mhc_i) and (o.proliferando or o.ligandos_estres)]
    if dianas:
        o = min(dianas, key=distancia)
        return Atacar(o.dy, o.dx) if adyacente(o) else hacia(o.dy, o.dx)
    if p.quimiotaxis and p.azar.random() < 0.6:
        return hacia(*p.quimiotaxis)
    return deambular(p)


def reina(ag: Agente, p):
    """Arregla el bug en la reina (célula madre tumoral) y extirpa las obreras corruptas."""
    return _estrategia_reina(ag, p, reparar=True, limpiar=True)


def reina_solo(ag: Agente, p):
    """Solo arregla las reinas; espera a que las obreras corruptas se agoten solas."""
    return _estrategia_reina(ag, p, reparar=True, limpiar=False)


def mata_reinas(ag: Agente, p):
    """No repara: mata a las reinas tumorales (tras secuenciarlas) y limpia obreras."""
    return _estrategia_reina(ag, p, reparar=False, limpiar=True)


TODOS_LOS_GENES = frozenset(PRIORIDAD)

TERAPIAS = {
    "cazador": dict(n=40, descripcion=cazador.__doc__, crear=lambda: Agente(
        "IA", "cazador", cazador, frozenset({"TCR", "NK"}), resistente_pdl1=True,
        expansion=0.25, vida=300)),
    "car_her2": dict(n=40, descripcion=car_her2.__doc__, crear=lambda: Agente(
        "IA", "CAR-HER2", car_her2, frozenset({"HER2"}), expansion=0.25, vida=300)),
    "reparador": dict(n=40, descripcion=reparador.__doc__, crear=lambda: Agente(
        "IA", "reparador", reparador, frozenset(), vectores=TODOS_LOS_GENES, carga=30, vida=400)),
    "hibrido": dict(n=40, descripcion=hibrido.__doc__, crear=lambda: Agente(
        "IA", "híbrido", hibrido, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES, carga=30,
        resistente_pdl1=True, expansion=0.15, vida=400)),
    "centinela": dict(n=30, descripcion=centinela.__doc__, crear=lambda: Agente(
        "IA", "centinela", centinela, frozenset(), vectores=TODOS_LOS_GENES, carga=80, vida=800)),
    "reina": dict(n=40, descripcion=reina.__doc__, crear=lambda: Agente(
        "IA", "reina", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES, carga=60,
        resistente_pdl1=True, expansion=0.15, vida=500)),
    "reina_solo": dict(n=40, descripcion=reina_solo.__doc__, crear=lambda: Agente(
        "IA", "reina (solo reparar)", reina_solo, frozenset(), vectores=TODOS_LOS_GENES, carga=60,
        vida=500)),
    "mata_reinas": dict(n=40, descripcion=mata_reinas.__doc__, crear=lambda: Agente(
        "IA", "mata-reinas", mata_reinas, frozenset({"INESPECIFICO"}), resistente_pdl1=True,
        expansion=0.15, vida=500)),
    "guardiana": dict(n=40, descripcion="Cura reinas y las convierte en guardianas con puerta Y",
                      crear=lambda: Agente(
        "IA", "reconductor (Y)", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES,
        carga=60, resistente_pdl1=True, expansion=0.15, vida=500, reprograma="doble")),
    "guardiana_simple": dict(n=40, descripcion="Igual, pero la guardiana ataca con una sola señal",
                             crear=lambda: Agente(
        "IA", "reconductor (1 señal)", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES,
        carga=60, resistente_pdl1=True, expansion=0.15, vida=500, reprograma="simple")),
    "guardiana_nicho": dict(n=40, descripcion="Guardiana con puerta Y posicional: marcador de madre Y fuera de la capa basal",
                            crear=lambda: Agente(
        "IA", "reconductor (Y nicho)", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES,
        carga=60, resistente_pdl1=True, expansion=0.15, vida=500, reprograma="nicho")),
    "guardiana_tert": dict(n=40, descripcion="Guardiana con puerta Y interna: telomerasa activa Y aneuploidía",
                           crear=lambda: Agente(
        "IA", "reconductor (Y TERT+aneu)", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES,
        carga=60, resistente_pdl1=True, expansion=0.15, vida=500, reprograma="tert")),
    "guardiana_tert_nicho": dict(n=40, descripcion="Guardiana con las dos puertas Y en paralelo (O de dos Y)",
                                 crear=lambda: Agente(
        "IA", "reconductor (Y TERT ó Y nicho)", reina, frozenset({"TCR", "NK"}),
        vectores=TODOS_LOS_GENES, carga=60, resistente_pdl1=True, expansion=0.15, vida=500,
        reprograma="tert+nicho")),
    "propagador": dict(n=4, descripcion=propagador.__doc__, crear=lambda: Agente(
        "IA", "propagador", propagador, frozenset({"TCR", "NK"}), resistente_pdl1=True,
        vida=120, propaga=2.0)),
    "cebo": dict(n=40, descripcion="Guardiana que deja cebos cargados en las senescentes",
                 crear=lambda: Agente(
        "IA", "cebo", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES, carga=60,
        resistente_pdl1=True, expansion=0.15, vida=500, reprograma="cebo")),
    "memoria_tert": dict(n=15, descripcion="Vigilancia residente con la puerta Y interna (TERT + aneuploidía)",
                         crear=lambda: Agente(
        "IA", "memoria (Y TERT+aneu)", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES,
        carga=200, resistente_pdl1=True, vida=10_000, reprograma="tert")),
    "memoria_tert_nicho": dict(n=15, descripcion="Vigilancia residente con las dos puertas Y en paralelo",
                               crear=lambda: Agente(
        "IA", "memoria (Y TERT ó Y nicho)", reina, frozenset({"TCR", "NK"}),
        vectores=TODOS_LOS_GENES, carga=200, resistente_pdl1=True, vida=10_000,
        reprograma="tert+nicho")),
    "memoria": dict(n=15, descripcion="Vigilancia residente de larga vida (como los linfocitos T de memoria)",
                    crear=lambda: Agente(
        "IA", "memoria", reina, frozenset({"TCR", "NK"}), vectores=TODOS_LOS_GENES, carga=200,
        resistente_pdl1=True, vida=10_000, reprograma="doble")),
    # ---- CABALLO DE TROYA -------------------------------------------------------
    # No pelean contra la presión: entran invitados. El tumor recluta monocitos con CCL2
    # y les abre el endotelio; el vehículo no nada contra corriente, lo llevan. Al llegar
    # al interior suelta la terapia real, que ya nace dentro y no tiene que atravesar nada.
    # Real: CAR-M (macrófagos CAR, en ensayo) y células madre mesenquimales como vehículo.
    "troyano": dict(n=40, descripcion="Monocito cargado: entra invitado y suelta la guardiana dentro",
                    crear=lambda: Agente(
        "IA", "vehículo troyano", politica_vehiculo, frozenset(), vida=400,
        troyano=0.9, entrega="guardiana_nicho")),
    "troyano_cebo": dict(n=40, descripcion="Vehículo que suelta cebos: deja que el tumor se lo coma",
                         crear=lambda: Agente(
        "IA", "troyano (cebo)", politica_vehiculo, frozenset(), vida=400,
        troyano=0.9, entrega="cebo")),
    "troyano_memoria": dict(n=20, descripcion="Vehículo que deja vigilancia residente dentro del tumor",
                            crear=lambda: Agente(
        "IA", "troyano (memoria)", politica_vehiculo, frozenset(), vida=400,
        troyano=0.9, entrega="memoria")),
    # Soldados que ESTERILIZAN en vez de matar. El objetivo no es vaciar el tumor, es
    # que no pueda crecer: lo que queda, jubilado, lo retiran los macrófagos o se extirpa.
    "castrador": dict(n=60, descripcion=politica_castrador.__doc__, crear=lambda: Agente(
        "IA", "castrador", politica_castrador, frozenset({"TCR", "NK"}),
        resistente_pdl1=True, expansion=0.2, vida=400, castra=True)),
    # El caballo de Troya completo: entra invitado, lleva su propio oxígeno para que la
    # carga no se asfixie en el núcleo, y suelta castradores dentro.
    "troya": dict(n=120, descripcion="Vehículo con bombona de O2 que suelta castradores dentro",
                  crear=lambda: Agente(
        "IA", "troya (O2 + castradores)", politica_vehiculo, frozenset(), vida=400,
        troyano=0.9, entrega="castrador", oxigeno=0.25)),
    "troya_sin_o2": dict(n=120, descripcion="Igual pero sin bombona: mide cuánto aporta el O2",
                         crear=lambda: Agente(
        "IA", "troya (sin O2)", politica_vehiculo, frozenset(), vida=400,
        troyano=0.9, entrega="castrador")),
    "quimio": dict(n=80, descripcion=quimio.__doc__, crear=lambda: Agente(
        "IA", "quimio", quimio, frozenset({"INESPECIFICO"}), vida=40)),
    "mi_celula": dict(n=mi_celula.CANTIDAD, descripcion=mi_celula.mi_politica.__doc__,
                      crear=mi_celula.crear),
}
