"""Documentación: avances, informes técnicos, bitácora y equivalencia de nombres."""

from __future__ import annotations

from itertools import groupby

from nicegui import ui

from breast_cancer_abm.app import content
from breast_cancer_abm.app.theme import NavLink, card, shell, sidepanel

EXTRAS = ["fenced-code-blocks", "tables", "header-ids"]


def page(doc: str = "") -> None:
    docs = content.all_docs()
    current = next((d for d in docs if d.key == doc), docs[0] if docs else None)
    with shell(
        "/docs",
        "Documentación",
        "Un documento por avance, con fuentes, resultados, errores corregidos y limitaciones; los informes "
        "técnicos; y la bitácora completa, con los tropiezos incluidos.",
    ):
        with ui.row().classes("w-full no-wrap items-start gap-6"):
            sidepanel(
                [
                    (group, [NavLink(f"/docs?doc={d.key}", d.title, d is current, d.date) for d in items])
                    for group, items in groupby(docs, key=lambda d: d.group)
                ]
            )
            with ui.column().classes("col gap-4").style("min-width: 0"):
                if current is None:
                    ui.label("No hay documentos.").classes("muted")
                elif current.key == "LOGBOOK.md":
                    _logbook()
                else:
                    with card():
                        ui.markdown(content.read_doc(current.key), extras=EXTRAS).classes("doc q-pa-sm")


def _logbook() -> None:
    entries = content.logbook_entries()
    with card():
        ui.markdown(content.logbook_intro(), extras=EXTRAS).classes("doc q-pa-sm")
        search = (
            ui.input(placeholder="Buscar en la bitácora (p. ej. ctDNA, guardiana, T15)")
            .props("outlined dense clearable")
            .classes("w-full")
        )
        count = ui.label().classes("hint")

    @ui.refreshable
    def listing() -> None:
        q = (search.value or "").strip().lower()
        shown = [(h, b) for h, b in entries if not q or q in h.lower() or q in b.lower()]
        count.set_text(f"{len(shown)} de {len(entries)} entradas, de la más reciente a la más antigua")
        for head, body in shown[:60]:
            with card(head, icon="event"):
                ui.markdown(body, extras=EXTRAS).classes("doc text-sm")
        if len(shown) > 60:
            ui.label("Afina la búsqueda para ver el resto.").classes("hint")

    search.on_value_change(lambda _: listing.refresh())
    listing()
