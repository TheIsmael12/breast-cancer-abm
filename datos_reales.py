"""Contraste con pacientes reales: TCGA PanCancer Atlas (~1.000 cánceres de mama) vía cBioPortal.

Responde a la pregunta: en tumores REALES, ¿se cumplen las suposiciones de nuestra estrategia
(guardiana + memoria)? Descarga los datos una vez a datos/ y luego trabaja sin conexión.

    python datos_reales.py            # usa la caché si existe
    python datos_reales.py --refrescar
"""
from __future__ import annotations

import argparse
import json
import statistics
import urllib.request
from collections import Counter
from pathlib import Path

API = "https://www.cbioportal.org/api"
ESTUDIO = "brca_tcga_pan_can_atlas_2018"
LISTA = f"{ESTUDIO}_all"
CACHE = Path(__file__).parent / "datos"

# Genes del modelo y otros drivers conocidos de mama (TCGA 2012; Nik-Zainal 2016).
ONCOGENES = {"PIK3CA": 5290, "ERBB2": 2064, "MYC": 4609, "CD274": 29126, "AKT1": 207,
             "CCND1": 595, "FGFR1": 2260, "ESR1": 2099, "FOXA1": 3169, "SF3B1": 23451}
SUPRESORES = {"TP53": 7157, "PTEN": 5728, "CDH1": 999, "BRCA1": 672, "BRCA2": 675,
              "GATA3": 2625, "MAP3K1": 4214, "MAP2K4": 6416, "KMT2C": 58508, "CBFB": 865,
              "RUNX1": 861, "NF1": 4763, "ARID1A": 8289, "TBX3": 6926, "NCOR1": 9611,
              "RB1": 5925, "PIK3R1": 5295, "CDKN2A": 1029}
PRESENTACION = {"HLA-A": 3105, "HLA-B": 3106, "HLA-C": 3107, "B2M": 567}
GENES = {**ONCOGENES, **SUPRESORES, **PRESENTACION}
NOMBRE = {v: k for k, v in GENES.items()}

GENES_MODELO = {"PIK3CA", "ERBB2", "MYC", "TP53", "PTEN", "CDH1", "BRCA1", "BRCA2"}
SIN_EFECTO = {"Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA", "IGR", "Splice_Region"}


def _pedir(ruta: str, cuerpo: dict | None = None):
    req = urllib.request.Request(
        API + ruta, data=json.dumps(cuerpo).encode() if cuerpo else None,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST" if cuerpo else "GET")
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def descargar(refrescar: bool = False) -> dict:
    CACHE.mkdir(exist_ok=True)
    fichero = CACHE / "tcga_brca.json"
    if fichero.exists() and not refrescar:
        return json.loads(fichero.read_text(encoding="utf-8"))
    cuerpo = {"entrezGeneIds": list(GENES.values()), "sampleListId": LISTA}
    datos = {
        "muestras": _pedir(f"/sample-lists/{LISTA}/sample-ids"),
        "mutaciones": [
            {"s": m["sampleId"], "g": m["entrezGeneId"], "t": m["mutationType"]}
            for m in _pedir(f"/molecular-profiles/{ESTUDIO}_mutations/mutations/fetch", cuerpo)],
        "cna": [
            {"s": c["sampleId"], "g": c["entrezGeneId"], "a": c["alteration"]}
            for c in _pedir(f"/molecular-profiles/{ESTUDIO}_gistic/discrete-copy-number/fetch"
                            "?discreteCopyNumberEventType=HOMDEL_AND_AMP", cuerpo)],
        "subtipo": {d["patientId"]: d["value"] for d in _pedir(
            f"/studies/{ESTUDIO}/clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE")},
        "tmb": {d["sampleId"]: float(d["value"]) for d in _pedir(
            f"/studies/{ESTUDIO}/clinical-data?clinicalDataType=SAMPLE&attributeId=TMB_NONSYNONYMOUS")},
    }
    fichero.write_text(json.dumps(datos), encoding="utf-8")
    return datos


def alteraciones(datos: dict) -> dict[str, set[str]]:
    """Por muestra: genes alterados con criterio funcional (oncogén: mutación o amplificación;
    supresor: mutación no silenciosa o deleción homocigota)."""
    alt: dict[str, set[str]] = {s: set() for s in datos["muestras"]}
    for m in datos["mutaciones"]:
        if m["t"] not in SIN_EFECTO and m["s"] in alt:
            alt[m["s"]].add(NOMBRE[m["g"]])
    for c in datos["cna"]:
        gen = NOMBRE[c["g"]]
        if c["s"] in alt and ((c["a"] == 2 and gen in ONCOGENES) or (c["a"] == -2 and gen not in ONCOGENES)):
            alt[c["s"]].add(gen)
    return alt


def pct(n: int, total: int) -> str:
    return f"{100 * n / total:5.1f} %"


def informe(datos: dict) -> str:
    alt = alteraciones(datos)
    muestras = list(alt)
    n = len(muestras)
    paciente = {s: s[:12] for s in muestras}
    subtipo = {s: datos["subtipo"].get(paciente[s], "?").replace("BRCA_", "") for s in muestras}
    tmb = datos["tmb"]
    drivers = {s: alt[s] - set(PRESENTACION) for s in muestras}
    L = [f"# Contraste con pacientes reales (TCGA PanCancer Atlas, {n} tumores de mama)", ""]

    L += ["## 1. ¿Están bien elegidos los genes del modelo?", "",
          "| Gen | Pacientes reales | Tipo |", "|---|---|---|"]
    frec = Counter(g for s in muestras for g in alt[s])
    for g, k in frec.most_common():
        if k / n >= 0.03 or g in GENES_MODELO:
            L.append(f"| {g}{' ★' if g in GENES_MODELO else ''} | {pct(k, n)} | "
                     f"{'oncogén' if g in ONCOGENES else 'presentación de antígeno' if g in PRESENTACION else 'supresor'} |")
    L += ["", "★ = gen incluido en el simulador.", ""]

    con_modelo = sum(bool(drivers[s] & GENES_MODELO) for s in muestras)
    todos_en_modelo = sum(bool(drivers[s]) and drivers[s] <= GENES_MODELO for s in muestras)
    sin_driver = sum(not drivers[s] for s in muestras)
    nd = [len(drivers[s]) for s in muestras]
    L += ["## 2. ¿Podrían las guardianas reparar el bug de cada paciente?", "",
          f"- Con al menos un error en los genes del modelo: **{pct(con_modelo, n)}**",
          f"- Con TODOS sus errores dentro de los genes del modelo (reparables con nuestros vectores): "
          f"**{pct(todos_en_modelo, n)}**",
          f"- Sin ningún driver de esta lista (su bug está en otra parte): {pct(sin_driver, n)}",
          f"- Errores driver por paciente: mediana {statistics.median(nd):g}, "
          f"máximo {max(nd)}; con 3 o más: {pct(sum(x >= 3 for x in nd), n)}",
          "", "Ojo: estos son solo ~30 genes; el genoma de un tumor de mama real tiene además "
          "miles de mutaciones pasajeras y reordenamientos.", ""]

    valores = [tmb[s] for s in muestras if s in tmb]
    altos = sum(v >= 10 for v in valores)
    sin_mhc = sum(bool(alt[s] & set(PRESENTACION)) for s in muestras)
    L += ["## 3. ¿Podría la puerta Y reconocer al tumor sin equivocarse?", "",
          "La guardiana necesita una señal EXCLUSIVA del tumor: neoantígenos (dependen de la carga",
          "mutacional, TMB) o pérdida de MHC-I.", "",
          f"- TMB mediana: **{statistics.median(valores):.2f} mutaciones/Mb** (melanoma o pulmón: ~10-20)",
          f"- Tumores con TMB alta (>= 10/Mb): **{pct(altos, len(valores))}**",
          f"- Con mutación o deleción en HLA-A/B/C o B2M: {pct(sin_mhc, n)}",
          "", "| Subtipo | Pacientes | TMB mediana | TP53 | PIK3CA | ERBB2 |", "|---|---|---|---|---|---|"]
    for st in ("LumA", "LumB", "Her2", "Basal", "Normal"):
        grupo = [s for s in muestras if subtipo[s] == st]
        if grupo:
            t = [tmb[s] for s in grupo if s in tmb]
            L.append(f"| {st} | {len(grupo)} | {statistics.median(t):.2f} | "
                     + " | ".join(pct(sum(g in alt[s] for s in grupo), len(grupo)) for g in ("TP53", "PIK3CA", "ERBB2"))
                     + " |")
    L += ["", "## Conclusión automática", ""]
    L.append(f"- La puerta Y del simulador da por hecho que todo tumor muestra neoantígenos. En la realidad, "
             f"la mayoría de los cánceres de mama tienen pocos (TMB mediana {statistics.median(valores):.1f}/Mb): "
             "**es el cuello de botella nº 1**. Hace falta otra señal exclusiva (p. ej. una mutación "
             "concreta como PIK3CA H1047R presentada en HLA, o combinaciones de proteínas de superficie).")
    L.append(f"- Reparar gen a gen solo cubriría por completo a ~{100 * todos_en_modelo / n:.0f} % de las "
             "pacientes con nuestro catálogo: **cada paciente necesitaría su propio juego de vectores** "
             "(medicina personalizada), igual que las vacunas de ARNm personalizadas.")
    L.append("- ERBB2/HER2 amplificado es la señal más limpia (proteína de superficie en exceso) y ya hay "
             "fármacos reales contra ella (trastuzumab, T-DXd).")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refrescar", action="store_true")
    a = ap.parse_args()
    texto = informe(descargar(a.refrescar))
    (Path(__file__).parent / "docs" / "informes" / "INFORME_DATOS_REALES.md").write_text(texto, encoding="utf-8")
    print(texto)


if __name__ == "__main__":
    main()
