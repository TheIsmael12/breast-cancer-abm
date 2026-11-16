"""Catálogo legible para personas: nombres en español de terapias, escenarios y parámetros.

El código usa claves internas en inglés (``trojan``, ``prob_neoantigen``); la aplicación
enseña estos nombres y explicaciones. Los textos salen de los comentarios de
``model/biology.py`` y ``model/policies.py``: si cambias un parámetro allí, revisa aquí.
"""

from __future__ import annotations

from dataclasses import dataclass

SCENARIOS = {
    "luminal": (
        "Luminal ER+",
        "El más frecuente (~70 %): una célula adquiere PIK3CA y los estrógenos la empujan.",
    ),
    "her2": ("HER2+", "Una célula amplifica ERBB2 y prolifera muy rápido."),
    "basal": (
        "Triple negativo / basal",
        "TP53 se pierde pronto; el subtipo más agresivo y sin diana hormonal.",
    ),
    "brca1": ("Portadora de BRCA1", "Todas las células nacen con un alelo roto; basta un segundo golpe."),
}

# clave -> (nombre, grupo)
THERAPIES = {
    "hunter": ("Linfocito cazador (CAR-T blindado)", "Matar"),
    "car_her2": ("CAR-T anti-HER2", "Matar"),
    "chemo": ("Quimioterapia citotóxica", "Matar"),
    "queen_killer": ("Mata-reinas", "Matar"),
    "sterilizer": ("Esterilizador (no mata, castra)", "Matar"),
    "propagator": ("Virus oncolítico", "Matar"),
    "repairer": ("Reparador génico", "Reparar"),
    "hybrid": ("Híbrido: repara o elimina", "Reparar"),
    "sentinel": ("Centinela preventivo", "Reparar"),
    "queen": ("Arreglar la reina", "Reparar"),
    "queen_only": ("Solo reinas", "Reparar"),
    "guardian": ("Guardiana (puerta Y)", "Guardianas"),
    "simple_guardian": ("Guardiana de una señal", "Guardianas"),
    "niche_guardian": ("Guardiana de nicho", "Guardianas"),
    "tert_guardian": ("Guardiana TERT + aneuploidía", "Guardianas"),
    "tert_niche_guardian": ("Guardiana de doble puerta", "Guardianas"),
    "bait": ("Guardiana con cebos", "Guardianas"),
    "memory": ("Vigilancia de memoria", "Vigilancia"),
    "tert_memory": ("Memoria con puerta TERT", "Vigilancia"),
    "tert_niche_memory": ("Memoria de doble puerta", "Vigilancia"),
    "trojan": ("Caballo de Troya (monocito)", "Vehículos"),
    "trojan_bait": ("Troyano con cebos", "Vehículos"),
    "trojan_memory": ("Troyano con memoria", "Vehículos"),
    "trojan_horse": ("Troyano con O₂ y esterilizadores", "Vehículos"),
    "trojan_horse_no_o2": ("Troyano sin O₂", "Vehículos"),
    "my_cell": ("Mi célula (programable)", "Propias"),
}
ARSENAL = "trojan+niche_guardian+memory"


def therapy_label(key: str) -> str:
    """Nombre legible de una terapia o de una combinación ``a+b+c``."""
    return " + ".join(THERAPIES.get(k, (k,))[0] for k in key.split("+"))


VESSEL_MODES = {
    "normal": "Normales",
    "reroute": "Reconducir el O₂ hacia el tumor",
    "normalize": "Normalizar (Jain)",
    "antiangiogenic": "Antiangiogénico",
    "starve": "Asfixiar (disruptor vascular)",
}


@dataclass(frozen=True)
class Param:
    key: str
    label: str
    help: str
    lo: float = 0.0
    hi: float = 1.0
    step: float = 0.01


P = Param
# (título, explicación, parámetros). El orden es el de la pantalla.
PARAM_GROUPS: list[tuple[str, str, list[Param]]] = [
    (
        "Genética y crecimiento",
        "Cómo se acumulan los errores y cuánto empujan a la célula a dividirse.",
        [
            P(
                "mutation_rate",
                "Tasa de mutación",
                "Por gen y división. La real es ~1e-6; el modelo la acelera para "
                "ver en días lo que tarda años. 2e-4 es la calibración «realista».",
                1e-5,
                2e-2,
                1e-5,
            ),
            P(
                "prob_division",
                "Probabilidad de división",
                "Por paso y célula sana con hueco; los oncogenes la multiplican.",
                0.001,
                0.05,
                0.001,
            ),
            P(
                "estrogen",
                "Estrógenos",
                "Multiplica la proliferación de las células ER+. El escenario fija su valor de partida.",
                0.0,
                1.5,
                0.05,
            ),
            P(
                "hrd_factor",
                "Inestabilidad sin BRCA (×)",
                "Cuántas más mutaciones acumula una célula sin recombinación homóloga.",
                1,
                30,
                1,
            ),
            P(
                "prob_hrd_damage",
                "Daño extra por división (HRD)",
                "Daño adicional en células sin BRCA1/2.",
                0,
                1,
            ),
            P(
                "prob_apoptosis_p53",
                "Apoptosis por p53",
                "Probabilidad de que p53 mande suicidarse a una célula dañada.",
                0,
                1,
            ),
            P(
                "prob_normal_death",
                "Recambio fisiológico",
                "Muerte natural por paso (las madres viven 10 veces más).",
                0,
                0.02,
                0.0005,
            ),
            P("prob_aneuploidy", "Aneuploidía sin p53", "Por división de una célula con p53 roto.", 0, 0.5),
            P(
                "prob_aneuploidy_p53ok",
                "Aneuploidía con p53",
                "Con p53 intacto es mucho más rara.",
                0,
                0.05,
                0.001,
            ),
        ],
    ),
    (
        "Células madre tumorales (reinas)",
        "Pocas células madre sostienen el tumor; las obreras se agotan.",
        [
            P(
                "normal_stem_fraction",
                "Fracción de madres sanas",
                "Proporción de células madre en el epitelio sano.",
                0,
                0.5,
            ),
            P(
                "prob_symmetric_division",
                "División simétrica",
                "Una reina da dos reinas (si no, reina + obrera).",
                0,
                1,
            ),
            P(
                "worker_divisions",
                "Divisiones de una obrera",
                "Límite de Hayflick: divisiones antes de la senescencia.",
                1,
                30,
                1,
            ),
            P(
                "plasticity",
                "Plasticidad (obrera → reina)",
                "Por paso; ×3 con hipoxia. La pregunta biológica decisiva del proyecto.",
                0,
                0.005,
                0.0001,
            ),
            P(
                "prob_sasp_plasticity",
                "Plasticidad por SASP",
                "Una senescente vecina empuja a otra a desdiferenciarse.",
                0,
                0.05,
                0.001,
            ),
            P(
                "stem_quiescence_factor",
                "Quiescencia de las madres",
                "Cuánto más despacio se dividen las madres.",
                0,
                1,
            ),
            P("prob_differentiated_death", "Muerte de diferenciadas", "Por paso.", 0, 0.1, 0.001),
        ],
    ),
    (
        "Sistema inmune",
        "Linfocitos T, NK y macrófagos naturales.",
        [
            P(
                "prob_neoantigen",
                "Inmunogenicidad",
                "Probabilidad de que una mutación sea un neoantígeno visible. "
                "1 = optimista; ~0,1 calibrado con TCGA.",
                0,
                1,
            ),
            P(
                "immune_recruitment",
                "Reclutamiento inmune",
                "Llegada de linfocitos T y NK desde los vasos.",
                0,
                0.3,
            ),
            P(
                "max_natural_immune",
                "Máximo de inmunes naturales",
                "Tope de T + NK en la rejilla.",
                0,
                200,
                1,
            ),
            P("prob_attack_kill", "Eficacia de un ataque", "Probabilidad de que un ataque mate.", 0, 1),
            P("macrophage_recruitment", "Reclutamiento de macrófagos", "Por paso.", 0, 0.5),
            P("max_macrophages", "Máximo de macrófagos", "Tope en la rejilla.", 0, 200, 1),
            P(
                "m2_macrophage_fraction",
                "Macrófagos M2",
                "Fracción polarizada a favor del tumor (no fagocitan).",
                0,
                1,
            ),
            P(
                "pd1_disengage_days",
                "Desenganche tras PD-1 (días)",
                "SUPUESTO: días que un T frenado no reintenta la misma célula.",
                0,
                10,
                0.5,
            ),
        ],
    ),
    (
        "Inmunidad adaptativa",
        "Expansión clonal, PD-L1 inducido y agotamiento (Kaech, Wherry, Siddiqui). "
        "Desactivada por defecto: varios supuestos sin calibrar.",
        [
            P(
                "adaptive_immunity",
                "Activar inmunidad adaptativa",
                "Expansión clonal + PD-L1 por IFN-γ + agotamiento y revigorización con anti-PD-1.",
            ),
            P(
                "program_divisions",
                "Divisiones del programa",
                "≥7 divisiones tras un estímulo breve.",
                0,
                15,
                1,
            ),
            P("prob_division_t", "División de un T por paso", "~1 división cada 6-8 h.", 0, 1),
            P("t_contraction", "Contracción tras el pico", "90-95 % de los efectores muere.", 0, 1),
            P("contraction_days", "Duración de la contracción (días)", "SUPUESTO.", 1, 30, 0.5),
            P("max_clonal_t", "Tope de T clonales", "SUPUESTO: 300 ≈ 110 CD8/mm².", 0, 2000, 10),
            P("ifng_radius", "Alcance del IFN-γ (celdas)", "SUPUESTO conservador (~30 µm).", 0, 10, 1),
            P("induced_pdl1_days", "PD-L1 inducido (días)", "SUPUESTO: dura ~48 h sin IFN-γ.", 0, 10, 0.5),
            P(
                "exhaustion_threshold",
                "Umbral de agotamiento",
                "Encuentros con antígeno para agotarse.",
                1,
                100,
                1,
            ),
            P("prob_tex_terminal", "Paso a terminal", "Por encuentro.", 0, 0.5),
            P("tex_prog_efficacy", "Eficacia del progenitor", "El TCF1+ mata poco.", 0, 1),
            P("tex_term_efficacy", "Eficacia del terminal", "Mata más pero no se divide.", 0, 1),
            P("tex_self_renewal", "Autorrenovación", "El progenitor se renueva despacio.", 0, 0.2),
            P("tex_term_lifespan_days", "Vida del terminal (días)", "SUPUESTO.", 1, 60, 1),
            P(
                "prob_reinvigoration",
                "Revigorización con anti-PD-1",
                "Por paso, con el anticuerpo circulando.",
                0,
                1,
            ),
            P("reinvigoration_divisions", "Divisiones de la revigorización", "Ráfaga corta.", 0, 10, 1),
        ],
    ),
    (
        "Microambiente: oxígeno, vasos y presión",
        "Lo que impide llegar al interior del tumor.",
        [
            P("o2_hypoxia", "Umbral de hipoxia", "Por debajo, la célula está hipóxica.", 0, 0.5),
            P("o2_necrosis", "Umbral de necrosis", "Por debajo, la célula muere.", 0, 0.2, 0.005),
            P(
                "prob_angiogenesis",
                "Angiogénesis",
                "VEGF de una célula tumoral hipóxica, por paso.",
                0,
                0.02,
                0.0005,
            ),
            P("o2_tumor_vessel", "O₂ de un vaso tumoral", "Los vasos que induce el tumor riegan mal.", 0, 1),
            P("o2_starved_vessel", "O₂ tras disruptor vascular", "Vasos cercanos al tumor asfixiados.", 0, 1),
            P(
                "interstitial_pressure",
                "Presión intersticial",
                "0 = sin presión (irreal); 1 = la del modelo.",
                0,
                2,
            ),
            P(
                "pressure_effect",
                "Freno de la presión",
                "Cuánto frena el gradiente al que intenta entrar.",
                0,
                20,
                0.5,
            ),
            P("pressure_drop_matrix", "Caída de presión al degradar la matriz", "", 0, 1),
            P("pressure_drop_normalize", "Caída de presión al normalizar vasos", "Jain 2005.", 0, 1),
        ],
    ),
    (
        "Invasión y malla",
        "La capa mioepitelial, la anoikis y la matriz que fabrican los CAF.",
        [
            P("prob_invasion", "Rotura de la membrana", "Por paso (metaloproteasas).", 0, 0.05, 0.0005),
            P("prob_anoikis", "Anoikis", "Muerte de una célula fuera de su sitio sin resistencia.", 0, 1),
            P("matrix_production", "Producción de matriz", "Por CAF activado y paso.", 0, 0.3),
            P("matrix_turnover", "Recambio de matriz", "Degradación fisiológica.", 0, 0.1, 0.001),
            P("lox", "Entrecruzamiento (LOX)", "Multiplica rigidez y contención.", 0, 3, 0.1),
            P(
                "caf_activation_tgfb",
                "Activación de CAF por TGF-β",
                "El tumor activa al fibroblasto vecino.",
                0,
                0.2,
            ),
            P(
                "caf_activation_stiffness",
                "Activación por rigidez (YAP/TAZ)",
                "La rigidez activa más CAF.",
                0,
                0.2,
            ),
            P(
                "matrix_containment",
                "Contención de la malla",
                "Cuánto frena la invasión una malla densa.",
                0,
                1,
            ),
            P(
                "context_strength",
                "Fuerza del contexto tisular",
                "Cuánto frena el tejido a una célula alterada (Weaver y Bissell, 1997).",
                0,
                1,
            ),
            P(
                "context_needs_matrix",
                "El contexto necesita la malla intacta",
                "Si la malla está rota, el contexto no manda.",
            ),
        ],
    ),
    (
        "Terapias celulares y edición génica",
        "Cómo de bien funcionan las células programadas.",
        [
            P(
                "repair_efficacy",
                "Éxito de una edición génica",
                "Probabilidad de que la corrección funcione.",
                0,
                1,
            ),
            P(
                "editing_error",
                "Error de edición",
                "Probabilidad de introducir un driver NUEVO al editar.",
                0,
                0.2,
                0.005,
            ),
            P(
                "correction_spread",
                "Corrección que se propaga",
                "Prob. por paso de corregir un gen de una vecina.",
                0,
                0.2,
            ),
            P(
                "correction_reach",
                "Alcance de la propagación",
                "Vecinas a las que llega (0 = apagado).",
                0,
                8,
                1,
            ),
            P(
                "prob_phagocytosis",
                "Fagocitosis de la guardiana",
                "Intentos por paso de comerse a una vecina.",
                0,
                1,
            ),
            P(
                "guardian_reversion",
                "Reversión de la guardiana",
                "Probabilidad por paso de perder el estado.",
                0,
                0.05,
                0.001,
            ),
            P(
                "guardian_inheritance",
                "Las hijas heredan el estado",
                "¿Lo heredan las hijas de la guardiana?",
            ),
            P(
                "internal_read_cost",
                "Coste de leer el interior",
                "Puertas TERT/aneuploidía: más lentas.",
                0,
                1,
            ),
            P(
                "internal_read_error",
                "Error de lectura interna",
                "Falsos positivos sobre célula sana.",
                0,
                0.2,
            ),
        ],
    ),
    (
        "Fármacos de diseño y CD47",
        "Telómeros, cebos y el bloqueo del «no me comas».",
        [
            P(
                "antitelomerase_efficacy",
                "Inhibidor de telomerasa",
                "Fracción de células en las que funciona.",
                0,
                1,
            ),
            P(
                "efficacy_6tiodg",
                "6-tio-dG",
                "Estropea el telómero de las células con telomerasa activa.",
                0,
                1,
            ),
            P("drug_p53", "El fármaco restaura p53", "La célula alcanzada muere (como Gendicine).", 0, 1),
            P("drug_bait", "Carga que envenena a quien se la come", ""),
            P("drug_anti_sasp", "Silenciar el SASP", "Inhibidores JAK (ruxolitinib).", 0, 1),
            P(
                "drug_persistent",
                "Fármaco persistente",
                "Fracción de células nuevas que lo siguen recibiendo.",
                0,
                1,
            ),
            P("prob_entosis", "Entosis", "Una tumoral se traga un cebo vecino, por paso.", 0, 0.3),
            P("bait_efficacy", "Eficacia del cebo", "Fracción de la carga que llega activa.", 0, 1),
            P("prob_bait_chain", "Cebo en cadena", "Si sobrevive, la comedora se vuelve cebo.", 0, 1),
            P(
                "anti_cd47_efficacy",
                "Anti-CD47: alcance",
                "Fracción de células a las que llega el bloqueo.",
                0,
                1,
            ),
            P(
                "cd47_reaches_healthy",
                "Anti-CD47 también llega a las sanas",
                "CD47 no es exclusivo del tumor.",
            ),
            P(
                "anti_cd47_duration",
                "Anti-CD47: duración (días)",
                "Días que el anticuerpo sigue circulando.",
                0,
                200,
                1,
            ),
        ],
    ),
]
PARAMS: dict[str, Param] = {p.key: p for _, _, ps in PARAM_GROUPS for p in ps}
