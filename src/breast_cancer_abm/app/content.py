"""Lectura de documentos, resultados y figuras del repositorio para la aplicación."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from breast_cancer_abm.paths import DOCS_DIR, FIGURES_DIR, LOGBOOK, PROGRESS_DIR, REPORTS_DIR, RESULTS_DIR

REPO_BLOB = "https://github.com/TheIsmael12/breast-cancer-abm/blob/main/"
DOCS_URL = "/files/docs"  # ruta estática servida por la aplicación

FIGURE_CAPTIONS = {
    "01_day_060_initial_lesion": "Día 60: la lesión inicial en el epitelio del conducto",
    "02_day_120_tumor": "Día 120: el tumor que se va a tratar",
    "03_local_treatment_comparison": "Tratamiento local: tumor, cirugía, ablación en disco y ablación guiada",
    "04_no_surgery_results": "Resultados de tratar sin operar (24 pacientes)",
    "05_3d_duct_cutaway": "Modelo 3D: el conducto cortado a lo largo",
    "06_3d_sections": "Modelo 3D: secciones del día 100, de la superficie al núcleo",
    "07_3d_disc_vs_guided_ablation": "3D: ablación en disco frente a ablación guiada",
}

SUMMARY_HEADERS = {
    "arm": "Brazo",
    "controlled": "Controladas",
    "eradicated": "Erradicadas",
    "healthy_destroyed": "Sanas destruidas (mediana)",
    "breast_untouched": "Mama sin tocar",
    "n": "Pacientes",
    "rate_after_treatment": "Tasa tras tratar",
    "file": "Salida",
}


@dataclass(frozen=True)
class Doc:
    key: str  # ruta relativa a docs/, p. ej. "progress/2026-11-05-...md"
    title: str
    group: str
    date: str = ""
    summary: str = ""


def _title(path: Path) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return path.stem


def progress_docs() -> list[Doc]:
    """Avances, del más reciente al más antiguo, con el resumen de su índice."""
    index = {}
    readme = PROGRESS_DIR / "README.md"
    if readme.exists():
        for line in readme.read_text(encoding="utf-8").splitlines():
            m = re.match(r"\|\s*([\d-]{10})\s*\|\s*\[([^\]]+)\]\(([^)]+\.md)\)\s*\|\s*(.+?)\s*\|$", line)
            if m:
                index[m.group(3)] = (m.group(1), m.group(2), m.group(4))
    docs = []
    for p in sorted(PROGRESS_DIR.glob("20*.md"), reverse=True):
        date, title, summary = index.get(p.name, (p.name[:10], _title(p), ""))
        docs.append(Doc(f"progress/{p.name}", title, "Avances", date, summary))
    return docs


def report_docs() -> list[Doc]:
    return [Doc(f"reports/{p.name}", _title(p), "Informes") for p in sorted(REPORTS_DIR.glob("*.md"))]


def other_docs() -> list[Doc]:
    out = [Doc("LOGBOOK.md", "Bitácora de investigación", "Proyecto")]
    if (DOCS_DIR / "NAMING.md").exists():
        out.append(Doc("NAMING.md", "Equivalencia de nombres", "Proyecto"))
    return out


def all_docs() -> list[Doc]:
    return progress_docs() + report_docs() + other_docs()


def read_doc(key: str) -> str:
    """Markdown de un documento con las rutas reescritas para la aplicación."""
    path = (DOCS_DIR / key).resolve()
    if DOCS_DIR.resolve() not in path.parents or not path.exists():
        return "# Documento no encontrado"
    text = path.read_text(encoding="utf-8")
    base = path.parent

    def link(m: re.Match) -> str:
        bang, label, target = m.group(1), m.group(2), m.group(3)
        if re.match(r"^(https?:|mailto:|#)", target):
            return m.group(0)
        anchor = ""
        if "#" in target:
            target, anchor = target.split("#", 1)
            anchor = "#" + anchor
        resolved = (base / target).resolve()
        try:
            rel_docs = resolved.relative_to(DOCS_DIR.resolve()).as_posix()
        except ValueError:
            rel_docs = None
        if rel_docs is not None:
            if bang or not rel_docs.endswith(".md"):
                return f"{bang}[{label}]({DOCS_URL}/{rel_docs})"
            return f"[{label}](/docs?doc={rel_docs})"
        try:
            rel_repo = resolved.relative_to(DOCS_DIR.parent.resolve()).as_posix()
        except ValueError:
            return m.group(0)
        if rel_repo.startswith("results/"):
            return f"[{label}](/results?set={rel_repo.split('/')[1] if '/' in rel_repo else ''})"
        return f"{bang}[{label}]({REPO_BLOB}{rel_repo}{anchor})"

    return re.sub(r"(!?)\[([^\]]*)\]\(([^)\s]+)\)", link, text)


def logbook_entries() -> list[tuple[str, str]]:
    """Bitácora partida en entradas (encabezado, cuerpo), de la más reciente a la más antigua."""
    text = LOGBOOK.read_text(encoding="utf-8")
    parts = re.split(r"^### ", text, flags=re.M)
    entries = []
    for part in parts[1:]:
        head, _, body = part.partition("\n")
        entries.append((head.strip(), body.strip()))
    return entries[::-1]


def logbook_intro() -> str:
    text = LOGBOOK.read_text(encoding="utf-8")
    return re.split(r"^### ", text, maxsplit=1, flags=re.M)[0]


@dataclass(frozen=True)
class ResultSet:
    date: str
    summary: list[dict]
    outputs: list[Path]


def result_sets() -> list[ResultSet]:
    sets = []
    for d in sorted((p for p in RESULTS_DIR.iterdir() if p.is_dir() and p.name[:2] == "20"), reverse=True):
        rows = []
        summary = d / "summary.csv"
        if summary.exists():
            with summary.open(encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh))
        sets.append(ResultSet(d.name, rows, sorted(d.glob("*.txt"))))
    return sets


def figures() -> list[tuple[str, str]]:
    """(url, pie de foto) de cada figura de docs/figures."""
    return [
        (f"{DOCS_URL}/figures/{p.name}", FIGURE_CAPTIONS.get(p.stem, p.stem))
        for p in sorted(FIGURES_DIR.glob("*.png"))
    ]
