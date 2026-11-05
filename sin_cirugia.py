"""¿Se puede llegar al mismo sitio que la cirugía SIN operar?

Lo que ha aprendido el mundo entre 2022 y 2026, y que este banco traslada al modelo:

  · Kuerer (MD Anderson; Lancet Oncol 2022, JAMA Oncol 2025): tratamiento sistémico
    PRIMERO, biopsia por vacío del lecho, y si sale limpia NO se opera. 31 pacientes
    HER2+ o triple negativo sin cirugía: 0 % de recaída en la mama a ~5 años.
  · FIRST (crioablación, 2026): destruir el tumor in situ con frío. Ablación completa
    89,6 % (96,9 % si ≤2 cm), pero CDIS residual en el BORDE en el 12,5 %.
  · MD Anderson, ESTRO 2026: hormonoterapia + radioterapia ablativa en luminales
    pequeños → 53 % de respuestas completas en biopsia, sin operar.
  · El punto débil común: la BIOPSIA de comprobación (RESPONDER: falsos negativos 17,8 %).

Brazos (todo el día 120, seguimiento hasta el día 500, misma semilla = misma paciente):
  nada · solo arsenal · cirugía + arsenal (la referencia del 22-10) · ablación + arsenal
  (borde letal y borde con supervivientes) · vía Kuerer (arsenal → biopsia → ablación
  solo si la biopsia sale positiva).

    python sin_cirugia.py --semillas 12 --procesos 6
"""
from __future__ import annotations

import argparse
import statistics
import time
from multiprocessing import Pool

from biologia import CALIBRACIONES, PASOS_POR_DIA, Parametros
from tejido import Simulacion

ARSENAL = "troyano+guardiana_nicho+memoria"
UMBRAL_CONTROL = 100   # suelo absoluto de la métrica (corrección del 2026-10-22 (7))

# (etiqueta, local: None | "cirugia" | "ablacion" | "kuerer", supervivencia_borde, arsenal)
BRAZOS = (
    ("nada", None, 0.0, False),
    ("solo arsenal (sin tocar)", None, 0.0, True),
    ("cirugía + arsenal", "cirugia", 0.0, True),
    ("ablación + arsenal (borde letal)", "ablacion", 0.0, True),
    ("ablación + arsenal (borde 30 % vivo)", "ablacion", 0.3, True),
    ("Kuerer: arsenal → biopsia → ablación", "kuerer", 0.3, True),
    # Puertas de prueba_residuo.py (2026-11-05): el ctDNA seriado, no de una vez.
    ("Kuerer + ctDNA mensual", "kuerer+ctdna", 0.3, True),
    ("solo ctDNA mensual (sin pinchar nunca)", "ctdna", 0.3, True),
    # 2026-11-05 (3): ablación guiada por imagen y tasa de mutación REALISTA tras tratar
    # (acelerada hasta el día del tratamiento para que haya tumor; 2e-4 después, para que
    # el campo no fabrique lesiones nuevas a ritmo de laboratorio).
    ("ctDNA + ablación GUIADA", "ctdna", 0.3, True, dict(guiada=True)),
    ("REAL: solo arsenal", None, 0.0, True, dict(tasa_tras="realista")),
    ("REAL: ctDNA + ablación disco", "ctdna", 0.3, True, dict(tasa_tras="realista")),
    ("REAL: ctDNA + ablación GUIADA", "ctdna", 0.3, True,
     dict(guiada=True, tasa_tras="realista")),
    # 2026-11-05 (4): la imagen solo ve focos de >= N células. Si el ctDNA da positivo y
    # no se ve nada, solo se repite el arsenal (tratamiento sistémico, sin tocar la mama).
    ("IMAGEN≥5: ctDNA + guiada", "ctdna", 0.3, True,
     dict(guiada=True, tasa_tras="realista", umbral_visible=5)),
    ("IMAGEN≥25: ctDNA + guiada", "ctdna", 0.3, True,
     dict(guiada=True, tasa_tras="realista", umbral_visible=25)),
    ("IMAGEN≥100: ctDNA + guiada", "ctdna", 0.3, True,
     dict(guiada=True, tasa_tras="realista", umbral_visible=100)),
    ("IMAGEN ciega: ctDNA → solo arsenal", "ctdna", 0.3, True,
     dict(guiada=True, tasa_tras="realista", umbral_visible=10 ** 9)),
)
# Fragmentos/célula calibrados en prueba_residuo.py para sensibilidad pre-tratamiento 75 %.
CTDNA_FRAGMENTOS = 9.91e-2
CTDNA_CADA = 30          # días entre extracciones
CTDNA_HASTA = 480


def correr(args) -> dict:
    brazo, semilla, dias, dia, dosis, espera_biopsia = args
    _, local, sup_borde, arsenal = brazo[:4]
    tasa_tras = (brazo[4] if len(brazo) > 4 else {}).get("tasa_tras")
    if local in ("kuerer+ctdna", "ctdna"):
        return correr_ctdna(brazo, semilla, dias, dia, dosis, espera_biopsia)
    sim = Simulacion("luminal", semilla, Parametros(prob_neoantigeno=0.1, error_edicion=0.01))
    p0 = int(dia * PASOS_POR_DIA)
    p_terapia = p0 + (4 * PASOS_POR_DIA if local in ("cirugia", "ablacion") else 0)
    p_biopsia = p0 + int(espera_biopsia * PASOS_POR_DIA)
    out = dict(sanas_destruidas=0, biopsia=None, real_en_biopsia=None, operado=0)
    for paso in range(int(dias * PASOS_POR_DIA)):
        if paso == p0 and tasa_tras:
            sim.p.tasa_mutacion = CALIBRACIONES[tasa_tras]
        if paso == p0 and local == "cirugia":
            out["sanas_destruidas"] = sim.resecar()["sanas"]
        elif paso == p0 and local == "ablacion":
            out["sanas_destruidas"] = sim.ablacionar(supervivencia_borde=sup_borde)["sanas"]
        if paso == p_terapia and arsenal:
            sim.liberar_terapia(ARSENAL, dosis)
        if local == "kuerer" and paso == p_biopsia:
            b = sim.biopsia()
            out["biopsia"], out["real_en_biopsia"] = b["positiva"], b["real"]
            if b["positiva"]:   # rescate local SIN bisturí, y segunda dosis del arsenal
                out["sanas_destruidas"] = sim.ablacionar(supervivencia_borde=sup_borde)["sanas"]
                out["operado"] = 1
        if local == "kuerer" and out["operado"] and paso == p_biopsia + 4 * PASOS_POR_DIA:
            sim.liberar_terapia(ARSENAL, dosis)
        sim.paso()
    r = sim.resumen()
    out.update(peligrosas=r["peligrosas"], reinas=r["reinas"], invasivas=r["invasivas"])
    return out


def correr_ctdna(brazo, semilla, dias, dia, dosis, espera_biopsia) -> dict:
    """Arsenal el día `dia`; a los `espera_biopsia` días primera comprobación (biopsia o
    ctDNA según el brazo) y después ctDNA cada CTDNA_CADA días. El primer positivo
    dispara ablación + arsenal UNA vez; después se sigue vigilando sin más rescates."""
    _, local, sup_borde, _ = brazo[:4]
    op = brazo[4] if len(brazo) > 4 else {}
    guiada, tasa_tras = op.get("guiada", False), op.get("tasa_tras")
    umbral = op.get("umbral_visible", 1)
    sim = Simulacion("luminal", semilla, Parametros(prob_neoantigeno=0.1, error_edicion=0.01))
    p0 = int(dia * PASOS_POR_DIA)
    primera = dia + espera_biopsia
    controles = {int(d * PASOS_POR_DIA) for d in range(int(primera), CTDNA_HASTA + 1, CTDNA_CADA)}
    out = dict(sanas_destruidas=0, biopsia=None, real_en_biopsia=None, operado=0,
               dia_rescate=None, extracciones=0, ablacion=0)
    p_rearsenal = None
    for paso in range(int(dias * PASOS_POR_DIA)):
        if paso == p0:
            sim.liberar_terapia(ARSENAL, dosis)
            if tasa_tras:
                sim.p.tasa_mutacion = CALIBRACIONES[tasa_tras]
        if paso in controles and not out["operado"]:
            if local == "kuerer+ctdna" and paso == int(primera * PASOS_POR_DIA):
                b = sim.biopsia()
                out["biopsia"], out["real_en_biopsia"] = b["positiva"], b["real"]
                positivo = b["positiva"] or sim.ctdna(CTDNA_FRAGMENTOS)["positiva"]
            else:
                positivo = sim.ctdna(CTDNA_FRAGMENTOS)["positiva"]
            out["extracciones"] += 1
            if positivo:
                ab = sim.ablacionar(supervivencia_borde=sup_borde, guiada=guiada,
                                    borde=1 if guiada else 2, umbral_visible=umbral)
                out["sanas_destruidas"] = ab["sanas"]
                out["ablacion"] = int(ab.get("focos_vistos", 1) > 0)
                out["operado"], out["dia_rescate"] = 1, round(sim.dia)
                p_rearsenal = paso + 4 * PASOS_POR_DIA
        if paso == p_rearsenal:
            sim.liberar_terapia(ARSENAL, dosis)
        sim.paso()
    r = sim.resumen()
    out.update(peligrosas=r["peligrosas"], reinas=r["reinas"], invasivas=r["invasivas"])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--semillas", type=int, default=12)
    ap.add_argument("--primera", type=int, default=1)
    ap.add_argument("--dias", type=float, default=500)
    ap.add_argument("--dia", type=float, default=120)
    ap.add_argument("--dosis", type=int, default=300)
    ap.add_argument("--espera-biopsia", type=float, default=60,
                    help="días entre el arsenal y la biopsia en la vía Kuerer")
    ap.add_argument("--procesos", type=int, default=6,
                    help="procesos en paralelo (deja núcleos libres)")
    ap.add_argument("--solo", default="",
                    help="ejecutar solo los brazos cuyo nombre contenga este texto")
    a = ap.parse_args()
    semillas = range(a.primera, a.primera + a.semillas)
    global BRAZOS
    BRAZOS = tuple(b for b in BRAZOS if a.solo.lower() in b[0].lower())

    print(f"SIN CIRUGÍA · luminal · tratamiento día {a.dia:g} · seguimiento a {a.dias:g} d · "
          f"{a.semillas} semillas · control = <{UMBRAL_CONTROL} peligrosas\n")
    print(f"{'brazo':<40}{'mediana':>9}{'control':>9}{'sanas dest.':>13}   por semilla")
    tareas = [(b, s, a.dias, a.dia, a.dosis, a.espera_biopsia) for b in BRAZOS for s in semillas]
    t0 = time.time()
    with Pool(a.procesos, maxtasksperchild=1) as pool:
        res = pool.map(correr, tareas, chunksize=1)
    n = len(semillas)
    for i, b in enumerate(BRAZOS):
        rs = res[i * n:(i + 1) * n]
        pel = [r["peligrosas"] for r in rs]
        ctrl = sum(p < UMBRAL_CONTROL for p in pel)
        sanas = statistics.median(r["sanas_destruidas"] for r in rs)
        print(f"{b[0]:<40}{statistics.median(pel):>9.0f}{ctrl:>6}/{n:<2}{sanas:>13.0f}   {pel}")
        if b[1] == "kuerer":
            neg = [r for r in rs if r["biopsia"] is False]
            fn = sum(r["real_en_biopsia"] > 0 for r in neg)
            vivos = sum(r["real_en_biopsia"] > 0 for r in rs)
            print(f"{'':<40}biopsia negativa (no se toca) en {len(neg)}/{n}; "
                  f"de ellas con tumor vivo en el tejido: {fn} "
                  f"(falsos negativos {fn}/{max(vivos, 1)} de las que tenían tumor)")
            print(f"{'':<40}desenlace de las NO tocadas: "
                  f"{[r['peligrosas'] for r in neg]}")
        if b[1] in ("kuerer+ctdna", "ctdna"):
            tocadas = [r for r in rs if r["operado"]]
            libres = [r for r in rs if not r["ablacion"]]
            print(f"{'':<40}sin tocar nunca la mama: {len(libres)}/{n} → desenlace "
                  f"{[r['peligrosas'] for r in libres]}")
            print(f"{'':<40}ctDNA positivo pero nada visible (solo arsenal): "
                  f"{sum(r['operado'] and not r['ablacion'] for r in rs)}/{n}")
            print(f"{'':<40}rescatadas: días {sorted(r['dia_rescate'] for r in tocadas)}")
    print(f"\n({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
