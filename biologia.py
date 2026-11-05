"""Conocimiento biológico del modelo: genes, parámetros y escenarios.

Cada número de este fichero es una SIMPLIFICACIÓN documentada. Las frecuencias
reales vienen de la literatura (ver README.md, sección Fuentes); las tasas de la
simulación están ACELERADAS para que un proceso que en una persona tarda 10-30
años se vea en minutos.
"""
from __future__ import annotations

from dataclasses import dataclass, field

HORAS_POR_PASO = 6
PASOS_POR_DIA = 24 // HORAS_POR_PASO


@dataclass(frozen=True)
class Gen:
    nombre: str
    tipo: str        # "oncogen": 1 golpe lo activa · "supresor": hacen falta 2 (hipótesis de Knudson)
    papel: str
    frecuencia: str  # en tumores de mama reales
    peso: float      # probabilidad relativa de recibir un golpe en cada división
    efecto: str      # qué hace en la simulación


GENES: dict[str, Gen] = {g.nombre: g for g in [
    Gen("PIK3CA", "oncogen", "Vía PI3K/AKT: señal de crecimiento y supervivencia siempre encendida",
        "~35 % (sobre todo luminales ER+)", 3.0, "proliferación ×2,5"),
    Gen("ERBB2", "oncogen", "Receptor HER2 amplificado: el acelerador atascado",
        "~15-20 % (amplificación)", 1.2, "proliferación ×3 · HER2 visible en superficie"),
    Gen("MYC", "oncogen", "Factor de transcripción del crecimiento",
        "amplificado en ~15-20 %, más en basales/BRCA1", 1.0, "proliferación ×2"),
    Gen("TP53", "supresor", "Guardián del genoma: frena el ciclo o manda suicidarse a la célula dañada",
        "~30-35 % (~80 % en triple negativo)", 2.5, "sin apoptosis por daño · inestabilidad ×3"),
    Gen("PTEN", "supresor", "Freno de la vía PI3K",
        "~5-10 % (más por pérdida de expresión)", 1.0, "proliferación ×1,5"),
    Gen("CDH1", "supresor", "E-cadherina: mantiene a las células pegadas y en su sitio",
        "~10-15 % (casi todos los lobulillares)", 1.0, "facilita la invasión"),
    Gen("BRCA1", "supresor", "Reparación del ADN por recombinación homóloga",
        "germinal en ~5 % de los casos", 0.6, "mutaciones ×10 (HRD)"),
    Gen("BRCA2", "supresor", "Reparación del ADN por recombinación homóloga",
        "germinal en ~2-3 % de los casos", 0.6, "mutaciones ×10 (HRD)"),
    Gen("HLA", "supresor", "MHC-I/B2M: el 'escaparate' que enseña al sistema inmune lo que hay dentro",
        "pérdida frecuente en tumores avanzados", 0.8, "invisible a los linfocitos T"),
    Gen("TERT", "oncogen", "Telomerasa: rebobina el contador de divisiones (telómeros). "
        "Sin ella ninguna célula puede dividirse sin límite",
        "activa en ~73-90 % de los cánceres; en mama la vía alternativa (ALT) es solo ~3 %",
        2.0, "divisiones ilimitadas: la célula se hace inmortal"),
    Gen("CD274", "oncogen", "PD-L1: la señal de 'no me ataques' para los linfocitos T",
        "alto en ~20 % de triple negativos", 0.6, "agota a los linfocitos T"),
]}

# Genes que, alterados por sí solos, no hacen a una célula neoplásica.
GENES_INMUNES = frozenset({"HLA", "CD274"})


@dataclass
class Parametros:
    tamano: int = 110                 # lado de la rejilla (1 píxel ≈ 1 célula ≈ 15 µm)
    radio_conducto: float = 22.0
    tasa_mutacion: float = 6e-3       # por gen y división. Real ~1e-6. Ver CALIBRACIONES
    factor_hrd: float = 10.0          # BRCA1/2 perdidos → reparación por recombinación rota
    prob_dano_hrd: float = 0.3        # daño extra por división en células HRD
    prob_apoptosis_p53: float = 0.35  # por cada golpe de daño si p53 funciona
    prob_division: float = 0.012      # por paso y célula sana con hueco (se multiplica por oncogenes)
    prob_muerte_normal: float = 0.002  # recambio fisiológico (las madres viven 10 veces más)
    prob_anoikis: float = 0.25        # muerte de una célula fuera de su sitio sin resistencia
    prob_invasion: float = 0.004      # rotura de la capa mioepitelial por paso (MMP)
    prob_angiogenesis: float = 0.002  # VEGF de una célula tumoral hipóxica
    estrogeno: float = 0.5            # multiplica la proliferación de las células ER+
    o2_hipoxia: float = 0.15
    o2_necrosis: float = 0.04
    reclutamiento_inmune: float = 0.05
    # MACRÓFAGOS y la puerta Y natural de la fagocitosis (CD47 / calreticulina).
    reclutamiento_macrofagos: float = 0.10
    max_macrofagos: int = 60
    eficacia_anti_cd47: float = 0.85   # fracción de células a las que llega el bloqueo...
    cd47_alcanza_sanas: bool = True    # ...y sí, también a las sanas: CD47 no es exclusivo
    duracion_anti_cd47: float = 60.0   # días que el anticuerpo sigue circulando (también
                                       # alcanza a las células que nacen durante ese tiempo)
    # Macrófagos asociados al tumor: la mayoría están polarizados a M2 por el propio tumor
    # y NO fagocitan; colaboran con él. Es otra razón de que el cuerpo no se lo coma solo.
    fraccion_macrofagos_m2: float = 0.6
    max_inmunes_naturales: int = 40
    eficacia_reparacion: float = 0.6  # probabilidad de que una edición génica funcione
    prob_muerte_ataque: float = 0.5
    # DESENGANCHE TRAS EL FRENO PD-1 (2026-10-28; ver Simulacion._atacar). El linfocito T
    # natural que recibe la señal PD-1 de una célula no vuelve a intentar matar a ESA célula
    # durante este tiempo; queda libre para buscar otra diana. SUPUESTO: 2 días, la misma
    # escala en la que cambia el PD-L1 inducido (dias_pdl1_inducido); sin calibrar y sin
    # análisis de sensibilidad.
    dias_desenganche_pd1: float = 2.0
    # Jerarquía de células madre (Al-Hajj 2003): pocas "reinas" sostienen el tumor.
    fraccion_madre_normal: float = 0.1
    prob_division_simetrica: float = 0.35  # reina → 2 reinas (si no, reina + obrera)
    divisiones_obrera: int = 8             # límite de Hayflick: divisiones antes de la senescencia
    eficacia_antitelomerasa: float = 0.9   # fracción de células en las que el inhibidor funciona
    # 6-tio-dG: la telomerasa lo incorpora al telómero y lo estropea YA, sin esperar al desgaste.
    # Solo afecta a células con telomerasa activa (tumorales y madres); las demás lo ignoran.
    eficacia_6tiodg: float = 0.85
    prob_sasp_plasticidad: float = 0.004   # una senescente vecina empuja a otra a desdiferenciarse
    # CEBO: en vez de retirar la célula senescente, la guardiana la deja CARGADA (p53 + 6-tio-dG).
    # El tumor se la come por entosis (real en mama: ~30 % de las MCF7 alojan otra célula dentro)
    # y recibe la carga. El veneno de los ADC como el T-DXd funciona así: salta a las vecinas.
    prob_entosis: float = 0.05             # por paso: una tumoral se traga un cebo vecino
    eficacia_cebo: float = 0.7             # fracción de la carga que llega activa
    prob_cebo_encadena: float = 0.5        # si sobrevive, la comedora se vuelve cebo a su vez

    # ---- FÁRMACO DE DISEÑO -------------------------------------------------
    # Cualidades que puede tener el agente que rompe telómeros. Sirven para responder:
    # ¿qué propiedades MÍNIMAS necesita para que ayude en vez de empeorar?
    farmaco_p53: float = 0.0        # también restaura p53 (como Gendicine): la célula muere
    farmaco_cebo: bool = False      # la carga se queda dentro: quien se coma la célula se envenena
    farmaco_anti_sasp: float = 0.0  # silencia el grito inflamatorio (inhibidores JAK: ruxolitinib)
    farmaco_persistente: float = 0.0  # fracción de células NUEVAS que siguen recibiendo el fármaco
    prob_muerte_diferenciada: float = 0.01
    factor_quiescencia_madre: float = 0.6
    plasticidad: float = 0.0003            # obrera → reina por paso (×3 con hipoxia)
    # Inmunogenicidad: probabilidad de que una mutación se presente como neoantígeno reconocible.
    # 1.0 = optimista (todo driver es visible). Calibrado con TCGA (TMB mediana 1,3/Mb): ~0.1.
    prob_neoantigeno: float = 1.0
    # ¿CUÁNTO DURA UN ESTADO CELULAR IMPUESTO? La guardiana es la pieza INVENTADA de
    # este proyecto: una reina curada y reprogramada que vigila para siempre y lo hereda.
    # Esa permanencia es una suposición, no un dato. Y la literatura de terapia de
    # diferenciación dice lo contrario: en tumores sólidos el estado inducido REVIERTE
    # cuando se retira el inductor (a diferencia de la APL con ATRA, donde se mantiene).
    # EL CONTEXTO MANDA SOBRE EL GENOTIPO (Weaver & Bissell, J Cell Biol 1997).
    # Células de cáncer de mama humano, en 3D y con un anticuerpo que bloquea la
    # β1-integrina, VOLVIERON a comportarse como normales: reconstruyeron su membrana
    # basal, recuperaron la E-cadherina, bajaron la ciclina D1, subieron p21 y dejaron
    # de dividirse. Con TODAS sus mutaciones intactas. No se corrigió ni un gen.
    #
    # Importa porque es la única intervención conocida cuyo coste NO crece con el número
    # de células: una señal que baña el tejido entero cuesta lo mismo para 10 que para
    # 10 millones. Todo lo demás en este simulador paga por célula, y por eso se agota.
    #
    # Pegas reales: es REVERSIBLE (al retirarla vuelven), la β1-integrina es esencial en
    # tejido sano, y funcionó en una matriz de laboratorio blanda y bien formada, no en
    # la malla rígida y rota de un tumor de verdad.
    fuerza_contexto: float = 0.0        # cuánto puede el tejido frenar a una célula alterada
    contexto_necesita_matriz: bool = True   # si la malla está rota, el contexto no manda
    # CORRECCIÓN QUE SE PROPAGA (efecto espectador / gene drive terapéutico).
    # Una célula ya corregida corrige a sus vecinas, y las que quedan corregidas corrigen
    # a su vez. Pasa de coste lineal (un agente por célula) a crecimiento exponencial.
    # Real a medias: el efecto espectador existe en terapia génica (HSV-TK por uniones
    # comunicantes) y los gene drives de CRISPR se copian solos en la población.
    # El peligro también es real: si se propaga sola no se puede parar, y la edición
    # génica tiene errores fuera de diana.
    propagacion_correccion: float = 0.0   # prob. por paso de corregir un gen de una vecina
    error_edicion: float = 0.0            # prob. de introducir un driver NUEVO al editar
    alcance_correccion: int = 0           # nº de vecinas a las que puede llegar (0 = apagado)
    reversion_guardiana: float = 0.0    # probabilidad por paso de perder el estado
    herencia_guardiana: bool = True     # ¿lo heredan las hijas?
    prob_fagocitosis: float = 0.15         # guardiana: intentos por paso de comerse a una vecina
    # Aneuploidía: número de cromosomas descuadrado. Aparece cuando se pierde el control
    # mitótico (sobre todo sin p53) y está en ~95 % de los tumores, pero en casi ninguna
    # célula sana. Problema: está DENTRO del núcleo, no se ve desde fuera.
    prob_aneuploidia: float = 0.12         # por división de una célula con p53 roto
    prob_aneuploidia_p53ok: float = 0.004  # con p53 intacto es mucho más rara
    # Puertas que leen el interior (telomerasa, cromosomas) necesitan un paso de lectura
    # (promotor hTERT, sonda de ploidía): más lentas y con algún error.
    coste_lectura_interna: float = 0.5     # multiplica los intentos de fagocitosis
    error_lectura_interna: float = 0.02    # falsos positivos sobre célula sana
    # LA MALLA (matriz desmoplásica) y quién la fabrica.
    # Los CAFs son al menos dos poblaciones (Öhlund 2017; Sahai 2020). Los que hacen la
    # malla son los myCAF, activados por el TGF-β del tumor. Fabrican colágeno I y ácido
    # hialurónico; el HA es un gel que retiene agua, y de ahí sale la presión. La LOX
    # ENTRECRUZA el colágeno: no añade material, lo tensa.
    produccion_matriz: float = 0.06     # por CAF activado y paso
    recambio_matriz: float = 0.012      # degradación fisiológica
    lox: float = 1.0                    # entrecruzamiento: multiplica rigidez y contención
    activacion_caf_tgfb: float = 0.03   # el tumor activa al fibroblasto vecino
    # Bucle YAP/TAZ (Calvo 2013): la rigidez que fabrica el CAF activa MÁS CAFs.
    # Es la razón de que la malla se mantenga sola y de que retirar material no baste.
    activacion_caf_rigidez: float = 0.02
    # La malla no solo estorba: CONTIENE. Eliminar los CAFs empeoró los tumores en ratones
    # (Özdemir 2014; Rhim 2014). Sin esto, el modelo nunca podría reproducir ese fracaso.
    contencion_matriz: float = 0.8      # cuánto frena la invasión una malla densa
    # PRESIÓN INTERSTICIAL: la barrera real de un tumor sólido. Los vasos tumorales son
    # porosos y los linfáticos están colapsados, así que el líquido se acumula dentro y
    # empuja HACIA FUERA. No es una pared: es corriente en contra. Es la razón principal
    # de que lo que se inyecta en vena no llegue al interior de la masa (Jain, 1987-2014).
    presion_intersticial: float = 1.0   # 0 = sin presión (irreal); 1 = la del modelo
    efecto_presion: float = 8.0         # cuánto frena el gradiente al que intenta entrar
    caida_presion_matriz: float = 0.5   # lo que baja la presión al degradar la matriz
    caida_presion_normalizar: float = 0.6  # normalizar los vasos también la baja (Jain 2005)
    o2_vaso_tumoral: float = 0.4           # los vasos que induce el tumor son caóticos y riegan mal
    o2_vaso_asfixiado: float = 0.3         # vasos cercanos al tumor tras un disruptor vascular

    # ---- INMUNIDAD ADAPTATIVA (añadido 2026-10-28, a raíz del fallo de T15) -------------
    # Diagnóstico de T15: el anti-PD-1 quitaba ~35.000 frenos y no cambiaba nada, porque
    # el límite era la CAPACIDAD (≤40 T+NK, sin expansión), no el freno. Faltaban tres
    # piezas de la biología, que se añaden aquí SIN mirar qué hace falta para que T15 pase:
    #
    #  (a) EXPANSIÓN CLONAL. Un CD8 que reconoce su antígeno arranca un programa autónomo
    #      de ≥7 divisiones aunque el estímulo dure solo horas (Kaech & Ahmed, Nat Immunol
    #      2001;2:415 · van Stipdonk, Nat Immunol 2001;2:423), a ~1 división cada 6-8 h
    #      una vez en marcha (Yoon, J Immunol 2010; De Boer, J Immunol 2003). Tras el pico,
    #      CONTRACCIÓN: muere el 90-95 % de los efectores (Kaech, Wherry & Ahmed, Nat Rev
    #      Immunol 2002;2:251 · Badovinac, Nat Immunol 2002).
    #  (b) PD-L1 ADAPTATIVO. El linfocito que reconoce su diana suelta IFN-γ y el IFN-γ
    #      induce PD-L1 en la diana y en sus vecinas por JAK1/2-IRF1: el tumor se pone el
    #      freno justo DONDE le atacan (Taube, Sci Transl Med 2012;4:127ra37 · Ribas,
    #      Cancer Discov 2015;5:915 · Garcia-Diaz, Cell Rep 2017;19:1189).
    #  (c) AGOTAMIENTO Y REVIGORIZACIÓN. La estimulación crónica por antígeno lleva al T a
    #      un estado agotado (Wherry, Nat Immunol 2011;12:492) con dos poblaciones: los
    #      progenitores TCF1+ (se autorrenuevan, matan poco) y los terminales TIM-3+ (no
    #      se dividen, viven poco). El anti-PD-1 actúa sobre los PROGENITORES: provoca una
    #      ráfaga proliferativa cuya progenie es efectora y acaba terminal; los terminales
    #      no responden (Im, Nature 2016;537:417 · Siddiqui, Immunity 2019;50:195 ·
    #      Miller, Nat Immunol 2019;20:326 · Hudson, Immunity 2019;51:1043).
    #
    # FLAG. Las tres piezas cambian el comportamiento POR DEFECTO de toda la inmunidad
    # natural, así que van detrás de `inmunidad_adaptativa`. Con False el modelo es
    # exactamente el de antes (mismo código, misma secuencia de números aleatorios).
    # El valor por defecto se decide tras correr la batería (ver nota más abajo).
    #
    # EXPECTATIVA, escrita ANTES de ejecutar nada (2026-10-28):
    #   · Sin fármaco, el PD-L1 adaptativo debería frenar a la propia expansión: cada
    #     ataque enciende el freno alrededor, y los clones se agotan. Espero que el tumor
    #     siga creciendo en ambos brazos sin tratar (quizá algo menos en el dMMR).
    #   · Con anti-PD-1 espero una mejora REAL pero PARCIAL y mayor en dMMR que en
    #     proficiente. NO espero respuestas completas: las reinas presentan neo//3, y con
    #     una semilla de 1 driver una reina dMMR tiene neo=1 → 0 visible. Ese escondite
    #     no lo toca ningún mecanismo de los añadidos. Probable veredicto de T15: FALLA o
    #     NO CONCLUYENTE por falta de RC, con diferencial dMMR > proficiente.
    #   · Resto de la batería: T13 no debería moverse (los T solo matan neo>0 y drivers>0,
    #     nunca sanas). Riesgos: T12 y T10 (si la inmunidad natural pasara a controlar el
    #     tumor, desaparecería la recaída y la citorreducción "curaría") y T1 (si no hay
    #     tumor grande, no hay núcleo necrótico). Espero que sigan pasando porque el freno
    #     adaptativo limita al T natural sin fármaco, pero es el riesgo a vigilar.
    #
    # RESULTADO MEDIDO (2026-10-28, después de escribir lo anterior):
    #   · Con False el modelo es idéntico al de antes (misma serie, semilla 1, 200 días).
    #   · T15 con True (8 semillas): sigue FALLANDO. Mejora 0/8 y RC 0/8 en los dos brazos;
    #     el día 500 casi todo está en el techo de la rejilla (~10.400) con o sin fármaco.
    #     La expansión sí retrasa el tumor dMMR antes del día 120 (mediana 174 frente a 288).
    #   · T1, T2, T7, T10, T12 y T13 (3 semillas) siguen PASANDO con True, pero T12 pierde
    #     margen: recaída en 1/3 semillas (antes 2/3). Con una semilla menos, fallaría.
    # DECISIÓN: por defecto False. No mejora ningún test, estrecha T12, lleva 9 SUPUESTOS
    # sin calibrar, y diseminacion.py (en desarrollo en paralelo) usa Parametros() por
    # defecto: cambiarlo le movería los resultados por debajo. Para estudiarlo:
    # Parametros(inmunidad_adaptativa=True).
    inmunidad_adaptativa: bool = False
    # (a) expansión y contracción
    divisiones_programa: int = 7          # Kaech & Ahmed 2001: ≥7 divisiones tras estímulo breve
    prob_division_t: float = 0.8          # por paso de 6 h: ~1 división cada 6-8 h (6/7,5 ≈ 0,8)
    contraccion_t: float = 0.925          # 90-95 % de los efectores muere tras el pico (Kaech 2002)
    dias_contraccion: float = 7.0         # SUPUESTO: la contracción se reparte en ~1 semana
                                          # (en ratón dura ~1-3 semanas; se toma el borde corto)
    # SUPUESTO: tope de linfocitos T clonales en la rejilla (2,7 mm²). 300 ≈ 110 CD8/mm²,
    # una densidad de TIL intermedia en mama; el techo real (espacio, IL-2, Treg) no se modela.
    max_t_clonales: int = 300
    # (b) PD-L1 inducido por IFN-γ
    radio_ifng: int = 2                   # SUPUESTO conservador: el IFN-γ llega a ~30 µm.
                                          # Hoekstra (Sci Adv 2020) y Thibaut (Nat Commun 2020)
                                          # lo ven actuar mucho más lejos (cientos de µm).
    dias_pdl1_inducido: float = 2.0       # SUPUESTO: PD-L1 inducido dura ~48 h sin IFN-γ
                                          # (in vitro sube en 24-48 h y baja al retirarlo)
    # (c) agotamiento
    umbral_agotamiento: int = 28          # SUPUESTO: encuentros con antígeno para entrar en
                                          # Tex progenitor (~1 semana a 1 encuentro por paso;
                                          # el programa de agotamiento se ve a los 5-8 días)
    prob_tex_terminal: float = 0.05       # SUPUESTO: por encuentro, Tex progenitor → terminal
    eficacia_tex_prog: float = 0.3        # SUPUESTO: el TCF1+ mata poco (poca granzima B; Miller 2019)
    eficacia_tex_term: float = 0.6        # SUPUESTO: el terminal mata más pero no se divide (Miller 2019)
    autorrenovacion_tex: float = 0.02     # SUPUESTO: el progenitor se autorrenueva despacio (~1/12 días)
    dias_vida_tex_term: float = 10.0      # SUPUESTO: el terminal vive poco
    # Revigorización por anti-PD-1 (solo progenitores; Siddiqui/Miller 2019)
    prob_revigorizacion: float = 0.1      # SUPUESTO: por paso, con el anticuerpo circulando
    divisiones_revigorizacion: int = 3    # SUPUESTO: ráfaga corta (Huang, Nature 2017;545:60:
                                          # pico de Ki-67 en CD8 PD-1+ a las ~3 semanas)


@dataclass(frozen=True)
class Escenario:
    descripcion: str
    germinal: dict = field(default_factory=dict)   # golpes que traen TODAS las células al nacer
    semilla: dict = field(default_factory=dict)    # golpes de una única célula iniciadora
    estrogeno: float = 0.5
    fraccion_er: float = 0.6


ESCENARIOS: dict[str, Escenario] = {
    "luminal": Escenario(
        "Esporádico ER+ (el más frecuente, ~70 %): una célula adquiere PIK3CA y los estrógenos la empujan",
        semilla={"PIK3CA": 1}, estrogeno=0.6, fraccion_er=0.7),
    "brca1": Escenario(
        "Portadora de BRCA1: todas las células nacen con un alelo roto; basta un segundo golpe",
        germinal={"BRCA1": 1}, estrogeno=0.3, fraccion_er=0.2),
    "basal": Escenario(
        "Triple negativo / basal: TP53 se pierde pronto (~83-90 % de estos tumores) y ya no hay "
        "guardián del genoma; es el subtipo más agresivo y sin diana hormonal",
        semilla={"TP53": 2, "MYC": 1}, estrogeno=0.15, fraccion_er=0.05),
    "her2": Escenario(
        "HER2+: una célula amplifica ERBB2 y prolifera muy rápido",
        semilla={"ERBB2": 1}, estrogeno=0.3, fraccion_er=0.5),
}


# ---------------------------------------------------------------- CALIBRACIONES
# Medido el 2026-10-22 sobre TEJIDO VIRGEN (sin célula sembrada), 6 semillas, 700 días,
# contando como "hay tumor" a partir de 100 células tumorales:
#
#     tasa     ×real   incidencia  día medio  premutadas máx  carga final
#     6e-03     6000      6/6         110           8            10.483
#     2e-03     2000      4/6         370          14               584
#     6e-04      600      3/6         470           2               568
#     2e-04      200      1/6         500           0                 2
#     6e-05       60      0/6          --           0                 0
#
# NO hay una tasa correcta única: depende de la pregunta, y el compromiso es real.
#
#   · Con la tasa acelerada el tumor aparece SIEMPRE (6/6) y rápido. Eso es cómodo para
#     estudiar TRATAMIENTO —que es para lo que se escribió el modelo— e INVALIDA
#     cualquier conclusión sobre prevención: si el tumor sale igual hagas lo que hagas,
#     no hay nada que prevenir. Es lo que hace fallar a T14 en validacion.py.
#   · Al bajarla, la incidencia pasa a depender del azar y la prevención se vuelve
#     estudiable, pero hacen falta muchas más semillas y muchos más días.
#   · A 2e-4 la incidencia (17 %) se acerca al riesgo real de por vida (~12 %), y
#     justo ahí deja de poder medirse nada: con 1 evento de 6 no se detecta una
#     reducción del 50 %. Es exactamente por lo que los ensayos de prevención reales
#     necesitan decenas de miles de participantes y décadas.
#
# Un aviso que vale para las tres: el reservorio de células PREMUTADAS (las que parecen
# sanas y llevan un solo golpe en un supresor) nunca pasa de 14 células a NINGUNA tasa.
# El estado intermedio no se sostiene como población, así que "reparar el tejido" no
# tiene a quién reparar, y eso NO se arregla calibrando.
CALIBRACIONES = {
    "tratamiento": 6e-3,   # por defecto: tumores fiables y rápidos (incidencia 6/6)
    "intermedia":  2e-3,   # incidencia 4/6: el mejor compromiso para prevención/campo
    "prevencion":  6e-4,   # incidencia 3/6, tumores al día ~470
    "realista":    2e-4,   # incidencia 1/6 (~17 %, cerca del 12 % real); casi inmedible
}


# ------------------------------------------------- PREVALENCIA EFECTIVA DE DIANAS
# La fracción de tumores que expresa una diana EN CRIBADO PROSPECTIVO, con el ensayo
# y el corte que se usan de verdad — que NO es la cifra que repite la literatura.
#
# El caso que obliga a distinguirlo: de folato-α se decía "~80 % de los TNBC".
# Al cribar 80 pacientes consecutivas la positividad real fue del 10 %, solo 2 de 96
# llegaron a recibir el fármaco y la tasa de respuesta fue 0 (Yam, Invest New Drugs
# 2021). El ensayo se cerró antes de tiempo. Se repite con Trop-2 (82-88 % publicado,
# 56 % con los cortes dicotómicos reales) y con CEACAM5 (24 % en TNBC).
#
# Por eso aquí van las cifras MEDIDAS, no las citadas. Valores para triple negativo.
PREVALENCIA_DIANA = {
    "trop2":    0.56,   # 56 % con los cortes de los ensayos (82-88 % "publicado")
    "b7h3":     0.70,   # 70 % con CPS>=50; hasta 90 % en metástasis cerebrales
    "nectina4": 0.58,
    "her3":     0.50,
    "egfr":     0.50,
    "muc1":     0.94,   # la mejor expresión... y clasificada de ALTO RIESGO por
                        # solapamiento completo con epitelio normal (UPTIDER)
    "gpnmb":    0.40,
    "folato":   0.10,   # el escarmiento: 10 % real frente a 80 % publicado
    "ceacam5":  0.24,
}
