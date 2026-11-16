"""Estudios: el catálogo de bancos de pruebas, con su pregunta y cómo ejecutarlos."""

from __future__ import annotations

from itertools import groupby

from nicegui import ui

from breast_cancer_abm.app.theme import REPO, card, ghost, primary, shell
from breast_cancer_abm.studies import STUDIES


def page() -> None:
    with shell(
        "/studies",
        "Estudios",
        "Cada estudio es un banco de simulaciones con su criterio de éxito escrito antes de ejecutar. Se lanzan "
        "desde la terminal con el comando que aparece en cada tarjeta. Para diseñar uno propio sin programar, "
        "usa el Laboratorio.",
    ) as pg:
        with pg.actions:
            primary("Ir al laboratorio", lambda: ui.navigate.to("/lab"), icon="biotech")
        studies = sorted(STUDIES.values(), key=lambda s: (s.topic, s.name))
        for topic, group in groupby(studies, key=lambda s: s.topic):
            ui.label(topic).classes("text-base font-semibold q-mt-sm")
            with (
                ui.grid()
                .classes("w-full gap-4")
                .style("grid-template-columns: repeat(auto-fill, minmax(440px, 1fr))")
            ):
                for s in group:
                    _card(s)


def _card(s) -> None:
    doc = s.doc.split("\n\n")
    paragraphs = [" ".join(line.strip() for line in p.splitlines()) for p in doc[1:] if p.strip()]
    body = " ".join(p for p in paragraphs if "bcabm study" not in p and not p.startswith("Diseño"))
    if len(body) > 360:
        body = body[:360].rsplit(" ", 1)[0] + "…"
    command = f"bcabm study {s.name}"
    source = s.module.replace("breast_cancer_abm.", "src/breast_cancer_abm/").replace(".", "/") + ".py"

    def actions() -> None:
        ghost("", lambda: ui.navigate.to(f"{REPO}/blob/main/{source}", new_tab=True), icon="code").tooltip(
            "Ver el código"
        )

    with card(s.title, actions=actions, classes="h-full"):
        if body:
            ui.label(body).classes("text-sm muted").style("line-height: 1.55")
        ui.space()
        with ui.row().classes("w-full items-center no-wrap gap-2"):
            ui.label(command).classes("mono command col")
            ghost(
                "", lambda: (ui.clipboard.write(command), ui.notify("Comando copiado")), icon="content_copy"
            ).tooltip("Copiar comando")
