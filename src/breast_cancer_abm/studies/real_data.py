"""Contraste con pacientes reales: TCGA PanCancer Atlas (~1.000 cánceres de mama) vía cBioPortal.

Responde a la pregunta: en tumores REALES, ¿se cumplen las suposiciones de nuestra estrategia
(guardiana + memoria)? Descarga los datos una vez a datos/ y luego trabaja sin conexión.

    bcabm study real-data            # usa la caché si existe
    bcabm study real-data --refresh
"""

from __future__ import annotations

import argparse
import json
import statistics
import urllib.request
from collections import Counter

from breast_cancer_abm.paths import DATA_DIR, REPORTS_DIR

API = "https://www.cbioportal.org/api"
STUDY = "brca_tcga_pan_can_atlas_2018"
LISTING = f"{STUDY}_all"
CACHE = DATA_DIR

# Genes del modelo y otros drivers conocidos de mama (TCGA 2012; Nik-Zainal 2016).
ONCOGENES = {
    "PIK3CA": 5290,
    "ERBB2": 2064,
    "MYC": 4609,
    "CD274": 29126,
    "AKT1": 207,
    "CCND1": 595,
    "FGFR1": 2260,
    "ESR1": 2099,
    "FOXA1": 3169,
    "SF3B1": 23451,
}
SUPPRESSORS = {
    "TP53": 7157,
    "PTEN": 5728,
    "CDH1": 999,
    "BRCA1": 672,
    "BRCA2": 675,
    "GATA3": 2625,
    "MAP3K1": 4214,
    "MAP2K4": 6416,
    "KMT2C": 58508,
    "CBFB": 865,
    "RUNX1": 861,
    "NF1": 4763,
    "ARID1A": 8289,
    "TBX3": 6926,
    "NCOR1": 9611,
    "RB1": 5925,
    "PIK3R1": 5295,
    "CDKN2A": 1029,
}
PRESENTATION = {"HLA-A": 3105, "HLA-B": 3106, "HLA-C": 3107, "B2M": 567}
GENES = {**ONCOGENES, **SUPPRESSORS, **PRESENTATION}
NAME = {v: k for k, v in GENES.items()}

MODEL_GENES = {"PIK3CA", "ERBB2", "MYC", "TP53", "PTEN", "CDH1", "BRCA1", "BRCA2"}
NO_EFFECT = {"Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "RNA", "IGR", "Splice_Region"}


def _request(file_path: str, body: dict | None = None):
    req = urllib.request.Request(
        API + file_path,
        data=json.dumps(body).encode() if body else None,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST" if body else "GET",
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def download(refresh: bool = False) -> dict:
    CACHE.mkdir(exist_ok=True)
    filename = CACHE / "tcga_brca.json"
    if filename.exists() and not refresh:
        return json.loads(filename.read_text(encoding="utf-8"))
    body = {"entrezGeneIds": list(GENES.values()), "sampleListId": LISTING}
    dataset = {
        "samples": _request(f"/sample-lists/{LISTING}/sample-ids"),
        "mutations": [
            {"s": m["sampleId"], "g": m["entrezGeneId"], "t": m["mutationType"]}
            for m in _request(f"/molecular-profiles/{STUDY}_mutations/mutations/fetch", body)
        ],
        "cna": [
            {"s": c["sampleId"], "g": c["entrezGeneId"], "a": c["alteration"]}
            for c in _request(
                f"/molecular-profiles/{STUDY}_gistic/discrete-copy-number/fetch"
                "?discreteCopyNumberEventType=HOMDEL_AND_AMP",
                body,
            )
        ],
        "subtype": {
            d["patientId"]: d["value"]
            for d in _request(f"/studies/{STUDY}/clinical-data?clinicalDataType=PATIENT&attributeId=SUBTYPE")
        },
        "tmb": {
            d["sampleId"]: float(d["value"])
            for d in _request(
                f"/studies/{STUDY}/clinical-data?clinicalDataType=SAMPLE&attributeId=TMB_NONSYNONYMOUS"
            )
        },
    }
    filename.write_text(json.dumps(dataset), encoding="utf-8")
    return dataset


def alterations(dataset: dict) -> dict[str, set[str]]:
    """Por muestra: genes alterados con criterio funcional (oncogén: mutación o amplificación;
    supresor: mutación no silenciosa o deleción homocigota)."""
    alt: dict[str, set[str]] = {s: set() for s in dataset["samples"]}
    for m in dataset["mutations"]:
        if m["t"] not in NO_EFFECT and m["s"] in alt:
            alt[m["s"]].add(NAME[m["g"]])
    for c in dataset["cna"]:
        gene = NAME[c["g"]]
        if c["s"] in alt and (
            (c["a"] == 2 and gene in ONCOGENES) or (c["a"] == -2 and gene not in ONCOGENES)
        ):
            alt[c["s"]].add(gene)
    return alt


def pct(n: int, total: int) -> str:
    return f"{100 * n / total:5.1f} %"


def report(dataset: dict) -> str:
    alt = alterations(dataset)
    samples = list(alt)
    n = len(samples)
    patient = {s: s[:12] for s in samples}
    subtype = {s: dataset["subtype"].get(patient[s], "?").replace("BRCA_", "") for s in samples}
    tmb = dataset["tmb"]
    drivers = {s: alt[s] - set(PRESENTATION) for s in samples}
    L = [f"# Contraste con pacientes reales (TCGA PanCancer Atlas, {n} tumores de mama)", ""]

    L += [
        "## 1. ¿Están bien elegidos los genes del modelo?",
        "",
        "| Gen | Pacientes reales | Tipo |",
        "|---|---|---|",
    ]
    freq = Counter(g for s in samples for g in alt[s])
    for g, k in freq.most_common():
        if k / n >= 0.03 or g in MODEL_GENES:
            L.append(
                f"| {g}{' ★' if g in MODEL_GENES else ''} | {pct(k, n)} | "
                f"{'oncogén' if g in ONCOGENES else 'presentación de antígeno' if g in PRESENTATION else 'suppressor'} |"
            )
    L += ["", "★ = gen incluido en el simulador.", ""]

    in_model = sum(bool(drivers[s] & MODEL_GENES) for s in samples)
    all_in_model = sum(bool(drivers[s]) and drivers[s] <= MODEL_GENES for s in samples)
    no_driver = sum(not drivers[s] for s in samples)
    nd = [len(drivers[s]) for s in samples]
    L += [
        "## 2. ¿Podrían las guardianas reparar el bug de cada paciente?",
        "",
        f"- Con al menos un error en los genes del modelo: **{pct(in_model, n)}**",
        f"- Con TODOS sus errores dentro de los genes del modelo (reparables con nuestros vectores): "
        f"**{pct(all_in_model, n)}**",
        f"- Sin ningún driver de esta lista (su bug está en otra parte): {pct(no_driver, n)}",
        f"- Errores driver por paciente: mediana {statistics.median(nd):g}, "
        f"máximo {max(nd)}; con 3 o más: {pct(sum(x >= 3 for x in nd), n)}",
        "",
        "Ojo: estos son solo ~30 genes; el genoma de un tumor de mama real tiene además "
        "miles de mutaciones pasajeras y reordenamientos.",
        "",
    ]

    values_list = [tmb[s] for s in samples if s in tmb]
    high = sum(v >= 10 for v in values_list)
    no_mhc = sum(bool(alt[s] & set(PRESENTATION)) for s in samples)
    L += [
        "## 3. ¿Podría la puerta Y reconocer al tumor sin equivocarse?",
        "",
        "La guardiana necesita una señal EXCLUSIVA del tumor: neoantígenos (dependen de la carga",
        "mutacional, TMB) o pérdida de MHC-I.",
        "",
        f"- TMB mediana: **{statistics.median(values_list):.2f} mutaciones/Mb** (melanoma o pulmón: ~10-20)",
        f"- Tumores con TMB alta (>= 10/Mb): **{pct(high, len(values_list))}**",
        f"- Con mutación o deleción en HLA-A/B/C o B2M: {pct(no_mhc, n)}",
        "",
        "| Subtipo | Pacientes | TMB mediana | TP53 | PIK3CA | ERBB2 |",
        "|---|---|---|---|---|---|",
    ]
    for st in ("LumA", "LumB", "Her2", "Basal", "Normal"):
        group = [s for s in samples if subtype[s] == st]
        if group:
            t = [tmb[s] for s in group if s in tmb]
            L.append(
                f"| {st} | {len(group)} | {statistics.median(t):.2f} | "
                + " | ".join(
                    pct(sum(g in alt[s] for s in group), len(group)) for g in ("TP53", "PIK3CA", "ERBB2")
                )
                + " |"
            )
    L += ["", "## Conclusión automática", ""]
    L.append(
        f"- La puerta Y del simulador da por hecho que todo tumor muestra neoantígenos. En la realidad, "
        f"la mayoría de los cánceres de mama tienen pocos (TMB mediana {statistics.median(values_list):.1f}/Mb): "
        "**es el cuello de botella nº 1**. Hace falta otra señal exclusiva (p. ej. una mutación "
        "concreta como PIK3CA H1047R presentada en HLA, o combinaciones de proteínas de superficie)."
    )
    L.append(
        f"- Reparar gen a gen solo cubriría por completo a ~{100 * all_in_model / n:.0f} % de las "
        "pacientes con nuestro catálogo: **cada paciente necesitaría su propio juego de vectores** "
        "(medicina personalizada), igual que las vacunas de ARNm personalizadas."
    )
    L.append(
        "- ERBB2/HER2 amplificado es la señal más limpia (proteína de superficie en exceso) y ya hay "
        "fármacos reales contra ella (trastuzumab, T-DXd)."
    )
    return "\n".join(L)


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args(argv)
    message = report(download(a.refresh))
    (REPORTS_DIR / "REAL_DATA.md").write_text(message, encoding="utf-8")
    print(message)


if __name__ == "__main__":
    main()
