"""Batería de validación: ¿reproduce el simulador los resultados reales?

    python validacion.py              # todos los tests
    python validacion.py --test T2    # uno solo
    python validacion.py --semillas 5 # más semillas (por defecto 3; T1, T16 y T17 usan
                                      # siempre ≥ 12 y T13 ≥ 24)

LA REGLA, y es la razón de que este fichero exista:

    Un fracaso clínico conocido NO se programa como parámetro: se usa como TEST.

Si el modelo, calibrado sólo con biología básica, reproduce por su cuenta lo que
pasó en la clínica, entonces sus predicciones sobre lo que nadie ha probado
empiezan a merecer atención. Si se le programa la conclusión, sólo devuelve el eco
de quien lo escribió — eso ya pasó en este proyecto con `contencion_matriz = 0.8`,
que "reprodujo" el fracaso de eliminar CAFs porque se le había dicho dos horas antes.

Un test que FALLA es tan útil como uno que pasa: dice qué le falta al modelo.
No se toca el modelo para que un test pase. Se anota el fallo y se decide aparte.

Escala: 1 celda ≈ 1 célula ≈ 15 µm (biologia.py:59).
"""
from __future__ import annotations

import argparse
import statistics
import sys

import numpy as np

# En Windows la consola usa cp1252 y rompe al redirigir a fichero (flechas, ·, µ).
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from biologia import PASOS_POR_DIA, Parametros
from tejido import Simulacion, LUMINAL

MICRAS_POR_CELDA = 15.0


class Resultado:
    def __init__(self, clave, titulo, fuente):
        self.clave, self.titulo, self.fuente = clave, titulo, fuente
        self.esperado = self.medido = ""
        self.pasa = None          # None = no concluyente
        self.nota = ""

    def veredicto(self):
        # numpy devuelve np.bool_ en las comparaciones: normalizar a bool de Python.
        return None if self.pasa is None else bool(self.pasa)

    def imprimir(self):
        marca = {True: "PASA", False: "FALLA", None: "¿?"}[self.veredicto()]
        print(f"\n[{marca}] {self.clave} — {self.titulo}")
        print(f"       esperado: {self.esperado}")
        print(f"       medido:   {self.medido}")
        if self.nota:
            print(f"       nota:     {self.nota}")
        print(f"       fuente:   {self.fuente}")


def _correr(semilla, dias, terapia=None, dia_terapia=120, dosis=300, params=None):
    p = params or Parametros(prob_neoantigeno=0.1)
    sim = Simulacion("luminal", semilla, p)
    for paso in range(int(dias * PASOS_POR_DIA)):
        if terapia and paso == int(dia_terapia * PASOS_POR_DIA):
            for t in terapia:
                sim.liberar_terapia(t, dosis)
        sim.paso()
    return sim


def _distancia_a_vaso(sim):
    """Distancia (en celdas) de cada posición al vaso más cercano, por BFS."""
    n = sim.p.tamano
    dist = np.full((n, n), 1e9)
    frente = [tuple(q) for q in np.argwhere(sim.vasos)]
    for q in frente:
        dist[q] = 0
    d = 0
    while frente:
        d += 1
        nuevo = []
        for y, x in frente:
            for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                q = (y + dy, x + dx)
                if 0 <= q[0] < n and 0 <= q[1] < n and dist[q] > d:
                    dist[q] = d
                    nuevo.append(q)
        frente = nuevo
    return dist


def _distancia_euclidea_a_vaso(sim):
    """Distancia EUCLÍDEA exacta (en celdas) de cada posición al centro de la celda de vaso
    más cercana. Fuerza bruta por bloques (rejilla 110×110, pocos cientos de vasos: barato);
    no usa scipy porque no está en requirements.txt."""
    n = sim.p.tamano
    vasos = np.argwhere(sim.vasos).astype(float)
    yy, xx = np.mgrid[0:n, 0:n]
    pts = np.stack([yy.ravel(), xx.ravel()], axis=1).astype(float)
    if len(vasos) == 0:
        return np.full((n, n), 1e9)
    mejor = np.full(len(pts), np.inf)
    for i in range(0, len(vasos), 64):
        v = vasos[i:i + 64]
        d2 = ((pts[:, None, :] - v[None, :, :]) ** 2).sum(axis=2).min(axis=1)
        mejor = np.minimum(mejor, d2)
    return np.sqrt(mejor).reshape(n, n)


# --------------------------------------------------------------------------- T1
T1_SEMILLAS_MIN = 12
T1_MIN_ANILLO = 15      # celdas tumorales mínimas para que un anillo cuente


def _t1_semilla(s):
    """Por anillo de distancia al vaso: (total tumoral, hipóxicas vivas, necróticas),
    y además la lista de distancias de las necróticas para la mediana antigua."""
    sim = _correr(s, 250)
    dist = _distancia_euclidea_a_vaso(sim)
    anillos, nec = {}, []
    for pos, c in sim.celulas.items():
        if not (c.necrosis or c.n_drivers > 0):
            continue
        d = int(round(float(dist[pos])))   # anillo d = distancia euclídea en [d−0,5, d+0,5)
        a = anillos.setdefault(d, [0, 0, 0])
        a[0] += 1
        if c.necrosis:
            a[2] += 1
            nec.append(d)
        elif sim.O2[pos] < sim.p.o2_hipoxia:
            a[1] += 1
    return anillos, nec


def T1(semillas):
    """Geometría del O2, medida CONDICIONADA A LA DISTANCIA y no al tamaño del tumor.

    Versión anterior (hasta 2026-10-28): mediana de la distancia al vaso de las celdas
    necróticas. Mide DÓNDE HAY necrosis, no DÓNDE EMPIEZA: un tumor más grande o más
    lejos de los vasos la sube aunque el perfil de O2 sea idéntico (correlación por
    semilla −0,7 con el tamaño y +0,85 con la distancia media al vaso). Con 3 semillas
    saltó de 240 a 330 µm por una sola semilla, y con 48 semillas daba 255 / 270 µm.

    Métrica actual (criterio escrito ANTES de correrla: scratchpad t1/CRITERIO_T1_ROBUSTO):
    se agrupan las celdas tumorales por anillos de distancia d al vaso y se mira qué
    fracción de cada anillo está hipóxica o necrótica. Thomlinson & Gray miden eso: el
    ESPESOR del anillo viable, la distancia a la que aparece la necrosis.
        frente de hipoxia = menor d con (hipóxicas + necróticas) / total ≥ 0,5
        inicio de necrosis = menor d con necróticas / total ≥ 0,5; si ningún anillo llega
                             (las necróticas se retiran a los 40 pasos y casi nunca son
                             mayoría), se usa 0,25 y se dice.
    Solo anillos con ≥ T1_MIN_ANILLO celdas. Pasa si el inicio de necrosis está en
    75-360 µm (150-180 µm con tolerancia 2×).

    Distancia EUCLÍDEA al vaso (desde 2026-10-28; antes BFS Manhattan). Por qué: la BFS con
    4 vecinos mide |dy| + |dx|, que sobreestima la distancia real hasta √2 en diagonal (~1,2×
    de media); a 15-25 celdas del vaso eso son ~3-5 celdas, ~75 µm de más, justo en la
    magnitud que T1 pretende medir. El O2 difunde en el plano, así que la distancia que
    compara con Thomlinson & Gray es la euclídea. T2 sigue con la BFS (_distancia_a_vaso):
    allí solo importa el orden de las bandas y cambiarla movería otro test."""
    from multiprocessing import Pool
    r = Resultado("T1", "Geometría del O2: ¿a qué distancia del vaso empieza la necrosis?",
                  "Thomlinson & Gray, Br J Cancer 1955;9:539 · Macklin, J Theor Biol 2012")
    r.esperado = ("viable hasta ~150-180 µm del vaso (10-12 celdas): frente de hipoxia "
                  "ahí y necrosis a partir de ahí; falla si el inicio de necrosis se desvía >2×")
    semillas = (tuple(semillas) if len(semillas) >= T1_SEMILLAS_MIN
                else tuple(range(1, T1_SEMILLAS_MIN + 1)))
    with Pool(maxtasksperchild=1) as pool:
        filas = pool.map(_t1_semilla, semillas, chunksize=1)
    anillos, nec = {}, []
    for an, ne in filas:
        nec += ne
        for d, v in an.items():
            a = anillos.setdefault(d, [0, 0, 0])
            for i in range(3):
                a[i] += v[i]
    validos = sorted(d for d, a in anillos.items() if a[0] >= T1_MIN_ANILLO)

    def primero(frac, umbral):
        return next((d for d in validos if frac(anillos[d]) >= umbral), None)

    d_hip = primero(lambda a: (a[1] + a[2]) / a[0], 0.5)
    umbral_nec = 0.5
    d_nec = primero(lambda a: a[2] / a[0], umbral_nec)
    if d_nec is None:
        umbral_nec = 0.25
        d_nec = primero(lambda a: a[2] / a[0], umbral_nec)

    def um(d):
        return "—" if d is None else f"{d} celdas = {d * MICRAS_POR_CELDA:.0f} µm"

    antigua = (f"{statistics.median(nec) * MICRAS_POR_CELDA:.0f} µm" if nec else "—")
    r.medido = (f"{len(semillas)} semillas · frente de hipoxia a {um(d_hip)} · inicio de "
                f"necrosis (≥{umbral_nec:.0%} del anillo) a {um(d_nec)}")
    r.nota = (f"métrica antigua (mediana de las necróticas, depende del tamaño): {antigua}"
              + ("; ningún anillo llega al 50 % necrótico, se usó el umbral 25 %"
                 if umbral_nec < 0.5 else ""))
    if d_nec is None:
        r.nota += "; sin anillo necrótico con muestra suficiente: no evaluable"
        return r
    r.pasa = 75 <= d_nec * MICRAS_POR_CELDA <= 360
    return r


def mejor_que(antes: float, despues: float, factor: float = 0.8, suelo: int = 100) -> bool:
    """¿`despues` es MEJOR que `antes` de forma que importe?

    Un umbral puramente relativo miente en los extremos: pasar de 1 a 8 células
    peligrosas es un +700 % y no significa nada — las dos son control total. Por eso
    hace falta un suelo absoluto por debajo del cual las diferencias no cuentan.
    Este defecto apareció comparando las tres capas y contaba 7 derrotas de 13 que
    en realidad eran empates.
    """
    if antes < suelo and despues < suelo:
        return False            # ambos son control: ni mejor ni peor
    return despues < antes * factor


def peor_que(antes: float, despues: float, factor: float = 1.25, suelo: int = 100) -> bool:
    if antes < suelo and despues < suelo:
        return False
    return despues > antes * factor


def _med(xs):
    if not xs:
        return "n/d"
    d = statistics.median(xs)
    return f"{d:.1f} celdas = {d * MICRAS_POR_CELDA:.0f} µm"


# --------------------------------------------------------------------------- T2
def T2(semillas):
    r = Resultado("T2", "Gradiente de penetración: ¿la terapia llega igual a todas partes?",
                  "Primeau, Tannock, Clin Cancer Res 2005;11:8782 · Minchinton, Nat Rev Cancer 2006")
    r.esperado = ("la letalidad debe caer con la distancia al vaso (a la mitad a ~40-50 µm "
                  "para una molécula pequeña; una CÉLULA terapéutica lo tiene aún peor)")
    # Hace falta un brazo CONTROL sin fármaco: a >150 µm está el núcleo necrótico,
    # donde las células mueren solas. Sin restar esa mortalidad basal, el test
    # atribuye al fármaco muertes por hipoxia y mide justo lo contrario de lo que dice.
    perfiles = []
    for s in semillas:
        anillos_brazo = {}
        for con_farmaco in (False, True):
            sim = _correr(s, 120)
            dist = _distancia_a_vaso(sim)
            antes = {pos: c.id for pos, c in sim.celulas.items()
                     if c.tipo == LUMINAL and c.n_drivers > 0}
            if con_farmaco:
                # Se mide con NUESTRA terapia, no con quimio: el proyecto no la usa,
                # y además una célula terapéutica penetra peor que una molécula pequeña,
                # así que el test queda más duro, no más blando.
                sim.liberar_terapia("guardiana_nicho", 300)
            for _ in range(40 * PASOS_POR_DIA):
                sim.paso()
            vivos = {c.id for c in sim.celulas.values()}
            anillos = {}
            for pos, cid in antes.items():
                k = int(min(dist[pos], 12)) // 2    # anillos de 2 celdas = 30 µm
                muertas, total = anillos.get(k, (0, 0))
                anillos[k] = (muertas + (cid not in vivos), total + 1)
            anillos_brazo[con_farmaco] = {k: m / t for k, (m, t) in anillos.items() if t >= 10}
        # letalidad ATRIBUIBLE al fármaco = con − sin
        comunes = set(anillos_brazo[True]) & set(anillos_brazo[False])
        perfiles.append({k: anillos_brazo[True][k] - anillos_brazo[False][k] for k in comunes})
    claves = sorted(set().union(*[set(p) for p in perfiles])) if perfiles else []
    linea = []
    for k in claves:
        vals = [p[k] for p in perfiles if k in p]
        if vals:
            linea.append(f"{k*2*MICRAS_POR_CELDA:.0f}-{(k*2+2)*MICRAS_POR_CELDA:.0f}µm:"
                         f"{statistics.median(vals):.0%}")
    r.medido = "letalidad ATRIBUIBLE al fármaco (con − sin) → " + (" · ".join(linea) or "sin datos")
    if len(claves) >= 2:
        cerca = statistics.median([p[claves[0]] for p in perfiles if claves[0] in p])
        lejos = statistics.median([p[claves[-1]] for p in perfiles if claves[-1] in p])
        if cerca <= 0.02:
            # No se puede medir el gradiente de algo que no mata: el agente "quimio"
            # resulta ser casi inerte (letalidad atribuible ~0 a toda distancia).
            r.pasa = None
            r.nota = (f"NO CONCLUYENTE: la letalidad atribuible cerca del vaso es {cerca:+.0%}. "
                      f"El fármaco no mata lo bastante para que exista un gradiente que medir. "
                      f"Antes de reevaluar T2 hay que darle efecto a 'quimio'.")
        else:
            r.nota = (f"cerca {cerca:+.0%} vs lejos {lejos:+.0%}; "
                      f"se espera que lejos sea ≤ la mitad de cerca")
            r.pasa = lejos <= cerca / 2
    return r


# --------------------------------------------------------------------------- T7
def T7(semillas):
    r = Resultado("T7", "Matar la 'reina' no termina la partida",
                  "Gupta, Cell 2011;146:633 · fRida NCT02370238 (HR 1,13, p=0,60)")
    r.esperado = "tras eliminar el compartimento madre, repoblación por transición de estado"
    rebrotes = []
    for s in semillas:
        p = Parametros(prob_neoantigeno=0.1)
        sim = Simulacion("luminal", s, p)
        minimo = None
        for paso in range(int(400 * PASOS_POR_DIA)):
            if paso == int(120 * PASOS_POR_DIA):
                sim.liberar_terapia("mata_reinas", 300)
            sim.paso()
            if paso > int(120 * PASOS_POR_DIA) and paso % (10 * PASOS_POR_DIA) == 0:
                reinas = sim.resumen()["reinas"]
                minimo = reinas if minimo is None else min(minimo, reinas)
        rebrotes.append((minimo, sim.resumen()["reinas"]))
    mins = [a for a, _ in rebrotes]
    fins = [b for _, b in rebrotes]
    r.medido = (f"reinas: mínimo tras tratar {statistics.median(mins):.0f} → "
                f"final {statistics.median(fins):.0f}  (por semilla {rebrotes})")
    r.pasa = statistics.median(fins) > statistics.median(mins)
    r.nota = "pasa si el compartimento madre se repuebla tras el mínimo"
    return r


# -------------------------------------------------------------------------- T10
def T10(semillas):
    r = Resultado("T10", "La citorreducción incompleta no debe ayudar",
                  "LION NEJM 2019 (HR 1,06) · CARMENA NEJM 2018 (18,4 vs 13,9 m) · GOG-0213")
    r.esperado = "matar mucho sin cerrar la fuente NO debe mejorar el desenlace"
    sin_t, con_t = [], []
    for s in semillas:
        sin_t.append(_correr(s, 400).resumen()["peligrosas"])
        con_t.append(_correr(s, 400, terapia=("cazador",), dosis=300).resumen()["peligrosas"])
    mejoras = sum(mejor_que(a, c, factor=0.5) for a, c in zip(sin_t, con_t))
    r.medido = (f"sin tratar {sin_t} → citorreducción {con_t}; "
                f"mejora >2× en {mejoras}/{len(semillas)} semillas")
    r.pasa = mejoras <= len(semillas) // 2
    r.nota = "falla si la citorreducción sola cura: eso contradice tres ensayos aleatorizados"
    return r


# -------------------------------------------------------------------------- T12
def T12(semillas):
    r = Resultado("T12", "El encogimiento del tumor NO es supervivencia",
                  "endpoint: ORR vs OS. Dato-DXd en mama: PFS sí, OS HR 1,01")
    r.esperado = "debe existir respuesta inicial seguida de recaída: el mínimo precoz no predice el final"
    filas = []
    for s in semillas:
        p = Parametros(prob_neoantigeno=0.1)
        sim = Simulacion("luminal", s, p)
        serie = []
        for paso in range(int(500 * PASOS_POR_DIA)):
            if paso == int(120 * PASOS_POR_DIA):
                for t in ("troyano", "guardiana_nicho", "memoria"):
                    sim.liberar_terapia(t, 300)
            sim.paso()
            if paso % (20 * PASOS_POR_DIA) == 0:
                serie.append(sim.resumen()["peligrosas"])
        tras = serie[7:]            # a partir del día ~140
        filas.append((min(tras) if tras else 0, serie[-1]))
    mins = [a for a, _ in filas]
    fins = [b for _, b in filas]
    recaen = sum(b > 10 * max(a, 1) for a, b in filas)
    r.medido = (f"mínimo tras tratar {statistics.median(mins):.0f} → final "
                f"{statistics.median(fins):.0f}; recaída >10× en {recaen}/{len(filas)} "
                f"(por semilla {filas})")
    r.pasa = recaen >= 1
    r.nota = "pasa si el modelo produce respuesta y recaída, no sólo un desenlace monótono"
    return r


# --------------------------------------------------------------------------- T13
T13_SEMILLAS_MIN = 24


def _t13_semilla(s):
    sim = _correr(s, 400, terapia=("troyano", "guardiana_nicho", "memoria"), dosis=300)
    st = sim.resumen()
    return dict(colateral=st["colateral"],
                comidas_sanas=st["comidas_macrofago_sanas"],
                castradas_sanas=st["castradas_sanas"],
                madres_sanas=st["madres_sanas_tocadas"],
                extirpadas_sanas=st["sanas_extirpadas"],
                muertes=st["muertes_terapia"] + st["muertes_inmunes"])


def T13(semillas):
    """El requisito del proyecto: eliminar lo tumoral SIN tocar al sujeto.

    SIEMPRE con al menos T13_SEMILLAS_MIN = 24 semillas, aunque --semillas pida menos.
    Por qué: el 2026-10-28 se descubrió que con las 3 semillas por defecto este test PASABA
    POR SUERTE. Con 24 semillas fallaban 3/24 sin el flag adaptativo (semillas 7, 16, 18) y
    2/24 con él (1, 16): una madre SANA suprabasal conservaba el marcador de madre y la
    puerta "nicho" se la comía (corregido en tejido.py, _madre_sana_sin_nicho). Un fallo
    con frecuencia ~1/8 por semilla pasa desapercibido con 3 semillas el 2/3 de las veces;
    con 24 la probabilidad de no verlo baja a ~4 %. Un requisito "cero daño" no se puede
    dar por cumplido con una muestra que no tiene potencia para ver el fallo.
    Se corre en paralelo (una semilla por proceso)."""
    from multiprocessing import Pool
    r = Resultado("T13", "No dañar al sujeto: cero daño colateral en tejido sano",
                  "requisito del proyecto · contraste real: Morgan 2010 (muerte por HER2 "
                  "en epitelio pulmonar) y vantictumab/ipafricept (fracturas óseas)")
    r.esperado = "daño a células sanas = 0 en todas las terapias de la estrategia"
    semillas = (tuple(semillas) if len(semillas) >= T13_SEMILLAS_MIN
                else tuple(range(1, T13_SEMILLAS_MIN + 1)))
    with Pool(maxtasksperchild=1) as pool:
        filas = pool.map(_t13_semilla, semillas, chunksize=1)
    total_dano = sum(f["colateral"] + f["comidas_sanas"] + f["castradas_sanas"]
                     + f["madres_sanas"] + f["extirpadas_sanas"] for f in filas)
    muertes = sum(f["muertes"] for f in filas)
    sucias = [s for s, f in zip(semillas, filas)
              if f["colateral"] + f["comidas_sanas"] + f["castradas_sanas"]
              + f["madres_sanas"] + f["extirpadas_sanas"]]
    r.medido = (f"{len(semillas)} semillas, con daño {len(sucias)}/{len(semillas)} {sucias}; "
                f"daño a sanas = {total_dano} sobre {muertes} muertes totales "
                f"({total_dano/max(muertes,1):.4%}); por semilla "
                f"{[(f['colateral'], f['comidas_sanas'], f['madres_sanas']) for f in filas]}")
    r.pasa = total_dano == 0
    r.nota = ("tuplas = (colateral, sanas comidas, madres sanas tocadas). "
              "Este test NO es negociable: si falla, la estrategia deja de cumplir "
              "su propio requisito, por muy bien que vayan los demás.")
    return r


# --------------------------------------------------------------------------- T14
def T14(semillas):
    """Prevención primaria. ESTE TEST FALLA Y DEBE SEGUIR FALLANDO hasta que se
    calibren las tasas de mutación: es el límite conocido del modelo."""
    r = Resultado("T14", "Prevención: ¿apagar el campo evita que aparezca el tumor?",
                  "tamoxifeno en prevención primaria: ~50 % menos incidencia en alto riesgo · "
                  "riesgo de cáncer de mama a lo largo de la vida ~12 %")
    r.esperado = ("incidencia < 100 % en tejido virgen, y ~50 % menos con el campo apagado")
    r.medido = ("MEDIDO 2026-10-22 (8 semillas, tejido virgen sin célula sembrada, 700 días): "
                "incidencia 8/8 con campo normal y 8/8 con campo apagado = 0 % de reducción; "
                "carga final 10.474 vs 10.478; el día de aparición pasa de 107 a 151 pero "
                "2 de 8 semillas ADELANTAN")
    r.pasa = False
    r.nota = ("FALLO CONOCIDO, EXPLICADO Y CON SALIDA. Causa: la tasa de mutación por "
              "defecto (6e-3, ~6000× la real) hace que CUALQUIER tejido produzca tumor, "
              "así que no hay nada que prevenir. Salida: biologia.CALIBRACIONES — con "
              "'intermedia' (2e-3) la incidencia baja a 4/6 y la prevención pasa a ser "
              "estudiable; con 'realista' (2e-4) es 1/6 (~17 %, cerca del 12 % real) pero "
              "ya no se puede medir una reducción del 50 % sin cientos de semillas, que es "
              "exactamente por lo que los ensayos de prevención reales necesitan decenas de "
              "miles de participantes. El modelo POR DEFECTO sirve para TRATAMIENTO; para "
              "prevención hay que cambiar de calibración y asumir el coste.")
    return r


# --------------------------------------------------------------------------- T15
# Brazos de T15. "dMMR" NO es un mecanismo nuevo: el modelo no tiene genes de
# reparación de desapareamientos (solo BRCA1/2 → HRD). Lo único que distingue a un
# tumor dMMR/MSI-H a efectos inmunes es su carga de neoantígenos, y eso ya existe:
# Parametros.prob_neoantigeno. Calibración real (TMB mediana 1,3/Mb) = 0.1; un MSI-H
# tiene ~20-50× más mutaciones, lo que satura la escala → 1.0 (todo golpe visible).
# NO se toca la tasa de mutación: subirla cambiaría también los drivers y el tamaño,
# y se compararía otra cosa (lección: una sola variable por brazo).
T15_BRAZOS = {"proficiente (neo 0.1)": 0.1, "tipo dMMR (neo 1.0)": 1.0}
T15_DIA_INICIO, T15_DIAS_BLOQUEO, T15_DIA_FINAL = 120, 180, 500   # 6 meses, como AZUR-1


def _t15_semilla(arg):
    """Un par emparejado: la MISMA simulación hasta el día 120 y, desde ahí, dos
    copias — una sigue sin tratar (control) y la otra recibe el anti-PD-1."""
    import copy
    neo, s = arg
    sim = Simulacion("luminal", s, Parametros(prob_neoantigeno=neo))
    for _ in range(int(T15_DIA_INICIO * PASOS_POR_DIA)):
        sim.paso()
    tum = [c for c in sim.celulas.values() if c.tipo == LUMINAL and c.n_drivers]
    base = dict(peligrosas=sim.resumen()["peligrosas"],
                pdl1=sum(c.pdl1 for c in tum) / max(len(tum), 1),
                visibles=sum(c.mhc and c.neo > 0 for c in tum) / max(len(tum), 1))
    tratada = copy.deepcopy(sim)
    tratada.bloquear_pd1(dias=T15_DIAS_BLOQUEO)
    out = {}
    for nombre, x in (("sin", sim), ("con", tratada)):
        for _ in range(int((T15_DIA_FINAL - T15_DIA_INICIO) * PASOS_POR_DIA)):
            x.paso()
        st = x.resumen()
        out[nombre] = dict(peligrosas=st["peligrosas"], frenos=st.get("frenos_pdl1", 0),
                           quitados=st.get("frenos_pdl1_bloqueados", 0))
    return neo, s, base, out


def T15(semillas):
    """Contraste dMMR vs pMMR ante el anti-PD-1. Nadie le ha dicho al modelo que el
    bloqueo funcione en tumores con muchas mutaciones: solo sabe que un linfocito T
    mata lo que presenta neoantígeno y que el PD-L1 lo agota."""
    r = Resultado("T15", "Anti-PD-1: ¿funciona en el tumor con muchos neoantígenos y no en el proficiente?",
                  "AZUR-1, dostarlimab en recto dMMR/MSI-H (FDA 19-10-2026; n=154): RC clínica "
                  "mantenida 12 m 47 %, ORR 79 %. pMMR/MSS: ORR ~0-5 % (Le, NEJM 2015;372:2509)")
    # DIANA Y CRITERIO, fijados ANTES de ejecutar (2026-10-28):
    #   mejora(semilla) = mejor_que(sin, con, factor=0.5, suelo=100) sobre las peligrosas
    #                     el día 500 (≥ 2× menos, con suelo absoluto de 100 células).
    #   RC(semilla)     = con anti-PD-1 quedan 0 peligrosas el día 500 Y sin él había ≥100
    #                     (respuesta completa ATRIBUIBLE al fármaco, mantenida ≥ 6 meses
    #                     después de retirarlo).
    #   PASA   si en el brazo dMMR mejora ≥ 50 % de semillas y hay RC en ≥ 1/3, y la
    #          fracción de mejora dMMR supera a la proficiente en ≥ 0,4 (proficiente ≤ 20 %).
    #   FALLA  si el freno PD-L1 nunca actúa en los controles (el fármaco no tiene diana:
    #          al modelo le falta el mecanismo), o si hay margen (dMMR sin tratar ≥ 100 en
    #          ≥ la mitad) y el diferencial es < 0,2 o no hay ninguna RC en dMMR.
    #   NO CONCLUYENTE si el dMMR sin tratar ya está controlado en más de la mitad de
    #          semillas (no hay nada que rescatar), o en cualquier caso intermedio.
    r.esperado = ("anti-PD-1 (6 meses desde el día 120): mejora ≥2× en ≥50 % de semillas dMMR "
                  "con RC en ≥1/3, frente a ≤20 % en el proficiente (diferencial ≥0,4)")
    from multiprocessing import Pool
    semillas = tuple(semillas) if len(semillas) >= 6 else tuple(range(1, 7))
    tareas = [(neo, s) for neo in T15_BRAZOS.values() for s in semillas]
    with Pool() as pool:
        filas = pool.map(_t15_semilla, tareas)
    lineas, frac, rcs, margen, frenos_ctrl = [], {}, {}, {}, 0
    for nombre, neo in T15_BRAZOS.items():
        fs = [(s, b, o) for n_, s, b, o in filas if n_ == neo]
        k = len(fs)
        mej = sum(mejor_que(o["sin"]["peligrosas"], o["con"]["peligrosas"], factor=0.5, suelo=100)
                  for _, _, o in fs)
        rc = sum(o["con"]["peligrosas"] == 0 and o["sin"]["peligrosas"] >= 100 for _, _, o in fs)
        frac[nombre], rcs[nombre] = mej / k, rc / k
        margen[nombre] = sum(o["sin"]["peligrosas"] >= 100 for _, _, o in fs) / k
        frenos_ctrl += sum(o["sin"]["frenos"] for _, _, o in fs)
        sin_ = [o["sin"]["peligrosas"] for _, _, o in fs]
        con_ = [o["con"]["peligrosas"] for _, _, o in fs]
        lineas.append(
            f"\n         {nombre}: día120 peligrosas med {statistics.median(b['peligrosas'] for _, b, _ in fs):.0f}, "
            f"PD-L1 {statistics.mean(b['pdl1'] for _, b, _ in fs):.0%}, visibles al TCR "
            f"{statistics.mean(b['visibles'] for _, b, _ in fs):.0%} · día500 sin {sin_} → con {con_} · "
            f"mejora {mej}/{k} · RC {rc}/{k} · frenos PD-L1 quitados "
            f"{sum(o['con']['quitados'] for _, _, o in fs)}")
    d, p = "tipo dMMR (neo 1.0)", "proficiente (neo 0.1)"
    dif = frac[d] - frac[p]
    r.medido = f"{len(semillas)} semillas por brazo, pares emparejados hasta el día 120:" + "".join(lineas)
    if frenos_ctrl == 0:
        r.pasa = False
        r.nota = ("el PD-L1 no frenó a ningún linfocito en los controles: el anti-PD-1 no tiene "
                  "diana. En el modelo el PD-L1 solo aparece por mutación de CD274; falta el PD-L1 "
                  "ADAPTATIVO (inducido por IFN-γ al ser atacada), que es el que existe en la clínica.")
    elif margen[d] <= 0.5:
        r.pasa = None
        r.nota = (f"NO CONCLUYENTE: sin tratar, el brazo dMMR ya está controlado en "
                  f"{1 - margen[d]:.0%} de semillas; no queda tumor que el bloqueo pueda rescatar.")
    elif frac[d] >= 0.5 and rcs[d] >= 1 / 3 and frac[p] <= 0.2 and dif >= 0.4:
        r.pasa = True
        r.nota = f"diferencial de mejora {dif:+.0%} (dMMR {frac[d]:.0%} vs proficiente {frac[p]:.0%})"
    elif dif < 0.2 or rcs[d] == 0:
        r.pasa = False
        r.nota = (f"sin contraste: mejora dMMR {frac[d]:.0%} vs proficiente {frac[p]:.0%} "
                  f"(diferencial {dif:+.0%}), RC dMMR {rcs[d]:.0%}")
    else:
        r.pasa = None
        r.nota = (f"NO CONCLUYENTE (intermedio): mejora dMMR {frac[d]:.0%} vs proficiente "
                  f"{frac[p]:.0%} (diferencial {dif:+.0%}), RC dMMR {rcs[d]:.0%}")
    return r


# --------------------------------------------------------------------------- T16
T16_SEMILLAS_MIN = 12
T16_N1, T16_N2, T16_DIAS = 100, 1000, 500


def _t16_semilla(s):
    sim = Simulacion("luminal", s, Parametros(prob_neoantigeno=0.1))
    t1 = t2 = None
    for d in range(1, T16_DIAS + 1):
        for _ in range(PASOS_POR_DIA):
            sim.paso()
        n = sim.resumen()["tumorales"]
        if t1 is None and n >= T16_N1:
            t1 = d
        if n >= T16_N2:
            t2 = d
            break
    return t1, t2


def T16(semillas):
    """Historia natural: ¿a qué ritmo crece el tumor sin tratar?

    Es el único contraste TEMPORAL que este modelo puede intentar honestamente. El resto
    no: la tasa de mutación por defecto está acelerada ~6000× (biologia.CALIBRACIONES) y
    el tejido virgen produce carcinoma de novo en ~110 días. Pero el RITMO de crecimiento
    de una masa ya formada no depende de esa tasa, sino de prob_division, del estrógeno y
    del O2, así que sí se puede comparar con el tiempo de duplicación medido en mama.

    El modelo es un corte 2D: N ∝ r², un tumor real V ∝ r³. Si N se dobla en T_celulas,
    el radio se multiplica por √2 y el volumen por 2^1,5, así que el equivalente de
    volumen es T_volumen = T_celulas / 1,5. El criterio va sobre T_volumen.

    CRITERIO, escrito ANTES de ejecutarlo (2026-10-28, scratchpad optimismo/CRITERIOS.md):
      PASA   si la mediana de T_volumen cae en 50-520 días (banda clínica 100-260 con
             tolerancia 2×).
      FALLA  'demasiado rápido' si < 50 d; 'demasiado lento' (optimista) si > 520 d.
      NO CONCLUYENTE si menos de 2/3 de las semillas llegan a 1000 células en 500 días."""
    from multiprocessing import Pool
    r = Resultado("T16", "Tiempo de duplicación del tumor sin tratar",
                  "Spratt, Cancer 1993;71:2020 · Kuroishi, Jpn J Cancer Res 1990;81:454 · "
                  "Peer, Cancer 1993;71:3547 (mediana ~100-260 días; rango 20->1000)")
    r.esperado = ("tiempo de duplicación de VOLUMEN con mediana en 50-520 días "
                  "(100-260 d con tolerancia 2×)")
    semillas = (tuple(semillas) if len(semillas) >= T16_SEMILLAS_MIN
                else tuple(range(1, T16_SEMILLAS_MIN + 1)))
    with Pool(maxtasksperchild=1) as pool:
        filas = pool.map(_t16_semilla, semillas, chunksize=1)
    import math
    tv = [((t2 - t1) * math.log(2) / math.log(T16_N2 / T16_N1)) / 1.5
          for t1, t2 in filas if t1 and t2 and t2 > t1]
    if len(tv) < 2 * len(semillas) / 3:
        r.medido = f"solo {len(tv)}/{len(semillas)} semillas llegan a {T16_N2} células"
        r.nota = "NO CONCLUYENTE: sin fase de crecimiento medible en la mayoría de semillas"
        return r
    med = statistics.median(tv)
    r.medido = (f"{len(semillas)} semillas · T_volumen mediana {med:.0f} d "
                f"(min {min(tv):.0f}, max {max(tv):.0f}) · T_celulas mediana {med*1.5:.0f} d")
    r.pasa = 50 <= med <= 520
    if not r.pasa:
        r.nota = ("el modelo crece DEMASIADO RÁPIDO: el tiempo del simulador está "
                  f"comprimido ~{200/med:.0f}× frente al de una paciente. No es optimismo, "
                  "es lo contrario; pero obliga a leer cualquier 'día' del modelo como un "
                  "tiempo comprimido, también en los seguimientos tras tratar."
                  if med < 50 else
                  "el modelo crece DEMASIADO LENTO: eso sí es optimismo y afecta a todo "
                  "desenlace medido a día fijo.")
    return r


# --------------------------------------------------------------------------- T17
T17_SEMILLAS_MIN = 12
T17_UMBRAL, T17_DIA_MAX = 200, 400
T17_CICLOS, T17_INTERVALO, T17_DOSIS, T17_ESPERA = 4, 21, 300, 21


def _t17_semilla(s):
    import copy
    sim = Simulacion("luminal", s, Parametros(prob_neoantigeno=0.1))
    d = 0
    while d < T17_DIA_MAX and sim.resumen()["tumorales"] < T17_UMBRAL:
        for _ in range(PASOS_POR_DIA):
            sim.paso()
        d += 1
    base = sim.resumen()
    if base["tumorales"] < T17_UMBRAL:
        return None
    tratada = copy.deepcopy(sim)
    for dia in range((T17_CICLOS - 1) * T17_INTERVALO + T17_ESPERA):
        if dia % T17_INTERVALO == 0 and dia // T17_INTERVALO < T17_CICLOS:
            tratada.liberar_terapia("quimio", T17_DOSIS)
        for _ in range(PASOS_POR_DIA):
            tratada.paso()
            sim.paso()
    return base, tratada.resumen(), sim.resumen()


def T17(semillas):
    """Quimio neoadyuvante en luminal (ER+): respuesta baja, pCR muy baja.

    El proyecto RECHAZA la quimioterapia como tratamiento, pero el agente ya existe
    (politicas.quimio) y es el único brazo del modelo que imita algo con dato clínico
    directo. Sirve como test de optimismo: si el modelo "cura" un ER+ con quimio, miente.
    Pares emparejados: la misma simulación hasta el disparo, y desde ahí copia tratada y
    copia control.

    CRITERIO, escrito ANTES de ejecutarlo (2026-10-28, scratchpad optimismo/CRITERIOS.md):
      pCR (modelo) = 0 invasivas y 0 reinas en el brazo tratado (equivale a ypT0/is: el
                     carcinoma in situ residual no cuenta como fallo de pCR).
      Respuesta objetiva = caída ≥30 % del DIÁMETRO equivalente, es decir N ≤ 0,49·N_basal
                     (RECIST mide diámetro, no área; en 2D 0,7² = 0,49).
      PASA   si pCR en 0-30 % y respuesta objetiva en 40-90 %.
      FALLA  'demasiado optimista' si pCR > 30 % o respuesta > 90 %.
      FALLA  'demasiado pesimista / inerte' si respuesta < 20 %, o si la letalidad
             atribuible (muertes_terapia tratado − control) es < 1 % del tumor basal."""
    from multiprocessing import Pool
    r = Resultado("T17", "Quimio neoadyuvante en ER+: pCR baja y respuesta parcial",
                  "Cortazar, Lancet 2014;384:164 · von Minckwitz, JCO 2012;30:1796 "
                  "(HR+/HER2−: pCR ~7-15 %; respuesta clínica ~60-70 %)")
    r.esperado = ("4 ciclos de quimio: pCR 0-30 % (clínica 7-15 %) y respuesta objetiva "
                  "40-90 % (clínica 60-70 %)")
    semillas = (tuple(semillas) if len(semillas) >= T17_SEMILLAS_MIN
                else tuple(range(1, T17_SEMILLAS_MIN + 1)))
    with Pool(maxtasksperchild=1) as pool:
        filas = [f for f in pool.map(_t17_semilla, semillas, chunksize=1) if f]
    if not filas:
        r.nota = "NO CONCLUYENTE: ninguna semilla alcanza el umbral de tumor tratable"
        r.medido = "sin datos"
        return r
    k = len(filas)
    pcr = sum(t["invasivas"] == 0 and t["reinas"] == 0 for _, t, _ in filas)
    resp = sum(t["tumorales"] <= 0.49 * b["tumorales"] for b, t, _ in filas)
    letal = statistics.median((t["muertes_terapia"] - c["muertes_terapia"]) / max(b["tumorales"], 1)
                              for b, t, c in filas)
    ratio = statistics.median(t["tumorales"] / max(c["tumorales"], 1) for _, t, c in filas)
    r.medido = (f"{k} semillas · pCR {pcr}/{k} = {pcr/k:.0%} · respuesta objetiva "
                f"{resp}/{k} = {resp/k:.0%} · tratado/control (mediana) {ratio:.2f}× · "
                f"letalidad atribuible {letal:.2f}× el tumor basal")
    if pcr / k > 0.30 or resp / k > 0.90:
        r.pasa = False
        r.nota = "DEMASIADO OPTIMISTA: la quimio no hace esto en un ER+."
    elif resp / k < 0.20 or letal < 0.01:
        r.pasa = False
        r.nota = ("DEMASIADO PESIMISTA, y no por falta de muertes: la quimio mata del orden "
                  "del tumor basal entero y aun así la masa no encoge, porque el tumor "
                  "repuebla a la vez (ver T16: duplicación ~30 d en células). Lo que falta "
                  "no es letalidad, es que la letalidad gane a la repoblación.")
    elif 0 <= pcr / k <= 0.30 and 0.40 <= resp / k <= 0.90:
        r.pasa = True
        r.nota = "respuesta parcial frecuente y pCR rara, como en HR+ real"
    else:
        r.nota = "NO CONCLUYENTE: caso intermedio"
    return r


# --------------------------------------------------------------------------- T18
def T18(semillas):
    """Efecto del MARGEN quirúrgico sobre la recaída: una RAZÓN, no una tasa.

    NO SE EJECUTA en la batería: necesita cambiar de calibración (2e-3, 'intermedia') y
    ~1.150 días por semilla y brazo. Se informa el resultado MEDIDO, como hace T14.

    POR QUÉ ESTE TEST NO ES EL DEL 2026-10-28. La primera versión comparaba la tasa
    ABSOLUTA de recaída del modelo con EBCTCG de cirugía conservadora (~35 %) y falló
    "por optimismo"; la causa que se escribió entonces fue que resecar() cortaba un disco
    sobre un foco único y compacto. Tres sondas de geometría del 2026-10-29 la refutaron,
    y ninguna de ellas mira el desenlace (scratchpad cirugia/):
      · geometria.py — el tumor ocupa el 90 % del anillo ductal en 5 componentes conexas:
        el modelo SÍ tiene multifocalidad, no había que añadirla.
      · que_corta.py — resecar(0.95, 2) extirpa el 93 % del conducto (49 % del campo):
        como el tumor rodea la luz, el centroide cae en el centro del conducto y el disco
        (radio ~28 celdas) engloba el conducto entero (radio 22).
      · sonda3_protocolo22.py — lo mismo el día 120 con la tasa por defecto: 96 %.
    O sea: resecar() NO es una tumorectomía en ninguna configuración, es la resección
    completa del único órgano del campo. Compararla con conservadora era compararla con
    otra operación, y el fallo "por optimismo" estaba mal atribuido.

    LÍMITE DE ESCALA, declarado: el campo mide 110 celdas × 15 µm = 1,65 mm de lado y el
    conducto 330 µm de radio; un margen quirúrgico clínico (1-2 mm) son 67-133 celdas,
    más que el semilado del campo. Ninguna cirugía con margen realista cabe aquí, así que
    la TASA absoluta de recaída no es contrastable. Una RAZÓN entre dos brazos del propio
    modelo sí lo es: el sesgo común de escala y de tiempo comprimido se cancela.

    CRITERIO, escrito ANTES de ejecutarlo (2026-10-29, scratchpad cirugia/
    CRITERIOS_T18_v2.md) y medido con las semillas 101-148, NUNCA usadas antes para que
    el contraste sea ciego:
      Brazo principal resecar(0.95, 2) sobre tumor de 400 células, seguimiento 250 días.
      La MISMA operación deja margen negativo (residuo 0) o positivo (residuo > 0) según
      el azar de la geometría: ésa es la situación de la que habla Houssami.
      RR = frac_recaída(margen positivo) / frac_recaída(margen negativo).
      PASA si RR ∈ [1,5 , 6,0] · FALLA 'el margen lo decide todo' si RR > 6 o si el margen
      negativo no recae NUNCA con ≥20 casos · FALLA 'el margen no importa' si RR < 1,5 ·
      NO CONCLUYENTE si el vaciado total recae ≥50 % o quedan <8 casos en un grupo."""
    r = Resultado("T18", "¿Pesa el margen quirúrgico lo que pesa en la clínica?",
                  "Houssami, Eur J Cancer 2010;46:3219 — meta-análisis, 28.162 pacientes: "
                  "margen positivo vs negativo, OR de recurrencia local 2,44 (1,97-3,03). "
                  "Secundaria: EBCTCG Lancet 2011;378:1707 (conservadora SIN RT ~35 %)")
    r.esperado = "RR de recaída entre margen positivo y negativo de 1,5 a 6,0 (ref. 2,44)"
    r.medido = (
        "MEDIDO 2026-10-29 · tasa 2e-3, semillas 101-148 (42 válidas), umbral 400 células, "
        "seguimiento 250 d · brazo resecar(0.95, 2): margen POSITIVO n=19 → recaída 32 %; "
        "margen NEGATIVO n=23 → recaída 0 % → RR = ∞ (15,6 con corrección de continuidad) "
        "· control de de novo (vaciado total) 7 %, así que la señal NO está tapada · "
        "informativo: resecar(0.90, 1) deja residuo 31 y recae 57 %")
    r.pasa = False
    r.nota = (
        "FALLA 'EL MARGEN LO DECIDE TODO', y la causa es estructural y se puede nombrar: "
        "en el modelo un margen negativo CURA SIEMPRE (0/23), mientras que en la clínica "
        "una conservadora sin radioterapia y con margen negativo recae igualmente en el "
        "20-25 % a 10 años. Esa diferencia es exactamente lo que no está en el campo: "
        "FUERA DEL CONDUCTO NO HAY MAMA donde esconder enfermedad oculta. El campo tiene "
        "UN solo conducto de 412 células (188 luminales + 224 mioepiteliales) sobre 1.251 "
        "del campo; la mama real tiene 15-20 sistemas ductales y la recaída tras "
        "conservadora viene en buena parte del resto de la mama, no del lecho. El modelo "
        "no es optimista ni pesimista aquí: le falta el órgano. Lo que SÍ reproduce, y no "
        "se le enseñó, es el SIGNO y la dosis-respuesta: más residuo, más recaída "
        "(0 → 0 % · 31 → 57 %), y una tasa de márgenes positivos del 45 % frente al "
        "20-25 % clínico, el mismo orden de magnitud. NO SE CORRIGE calibrando la "
        "cobertura ni sembrando focos fuera del conducto: lo primero es programar el "
        "resultado y lo segundo exige un campo mayor (varios conductos), que es un cambio "
        "de modelo, no un parámetro. Queda declarado como límite, igual que T14.")
    return r


TESTS = {"T1": T1, "T2": T2, "T7": T7, "T10": T10, "T12": T12, "T13": T13, "T14": T14,
         "T15": T15, "T16": T16, "T17": T17, "T18": T18}
PENDIENTES = {
    "T3": "ADC con payload permeable — requiere implementar el efecto espectador del citotóxico",
    "T4": "degradar el estroma no es gratis — ya medido en malla.py, falta formalizarlo aquí",
    "T5": "equilibrio inmune (Koebel 2007): requiere depleción selectiva de CD4/CD8/NK",
    "T6": "checkpoint en tumor frío — el anti-PD-1 ya existe (Simulacion.bloquear_pd1, ver T15); falta formalizar el brazo frío",
    "T8": "resistencia evolutiva y calendario — requiere dosificación adaptativa on/off",
    "T9": "apilar debe COSTAR — requiere atenuación de dosis (Kurzrock: 51 %/38 %/28 %)",
    "T11": "dormancia y recurrencia tardía — requiere células quiescentes persistentes",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", nargs="*", default=list(TESTS))
    ap.add_argument("--semillas", type=int, default=3)
    a = ap.parse_args()
    semillas = tuple(range(1, a.semillas + 1))

    print(f"BATERÍA DE VALIDACIÓN · {a.semillas} semillas · escala {MICRAS_POR_CELDA:.0f} µm/celda")
    print("Un test que falla dice qué le falta al modelo. No se toca el modelo para que pase.")
    resultados = []
    for clave in a.test:
        if clave not in TESTS:
            print(f"\n[--] {clave}: {PENDIENTES.get(clave, 'desconocido')}")
            continue
        r = TESTS[clave](semillas)
        r.imprimir()
        resultados.append(r)

    pasa = sum(r.veredicto() is True for r in resultados)
    falla = sum(r.veredicto() is False for r in resultados)
    duda = sum(r.veredicto() is None for r in resultados)
    print(f"\n{'='*70}\nRESUMEN: {pasa} pasan · {falla} fallan · {duda} no concluyentes")
    if PENDIENTES:
        print(f"\nSin implementar ({len(PENDIENTES)}):")
        for k, v in PENDIENTES.items():
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
