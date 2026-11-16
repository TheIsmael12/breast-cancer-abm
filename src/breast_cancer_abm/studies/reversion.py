"""¿Cuánto tiene que durar la reprogramación para que la estrategia siga valiendo?

La guardiana —reina curada y reprogramada que vigila el tejido— es la pieza INVENTADA
de este proyecto: no existe hoy en ningún laboratorio (ver HOJA_DE_RUTA.md, fase B). Y
en el modelo tenía dos propiedades que nadie ha justificado: dura para siempre y se
hereda al dividirse.

La literatura de terapia de diferenciación dice lo contrario. En la leucemia
promielocítica (APL) el ATRA impone un estado que se mantiene y cura a >90 %. Pero en
tumores SÓLIDOS la diferenciación inducida REVIERTE: al retirar el inductor las células
recuperan su fenotipo proliferativo, y esa es una de las razones por las que ningún
agente de diferenciación ha curado un tumor sólido.

Así que la pregunta no es "¿funciona la guardiana?", sino "¿cuánto tiene que durar?".
Eso sí es contestable en un laboratorio, y es barato comparado con todo lo demás.

    bcabm study reversion --seeds 3
"""

from __future__ import annotations

import argparse
import time

from breast_cancer_abm.model.biology import SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation

# (etiqueta, probabilidad de revertir por paso) — la vida media sale de 1/p pasos.
LEVELS = (
    ("permanente (lo que asumíamos)", 0.0),
    ("vida media ~500 días", 0.0005),
    ("vida media ~125 días", 0.002),
    ("vida media ~25 días", 0.01),
    ("vida media ~5 días", 0.05),
)


def run(
    therapy: str,
    dose: int | None,
    scenario: str,
    seed: int,
    days: float,
    therapy_day: float,
    reversion: float,
    inheritance: bool,
    neoantigens: float,
) -> dict:
    p = Parameters(
        prob_neoantigen=neoantigens, guardian_reversion=reversion, guardian_inheritance=inheritance
    )
    sim = Simulation(scenario, seed, p)
    for step in range(int(days * STEPS_PER_DAY)):
        if step == int(therapy_day * STEPS_PER_DAY):
            sim.release_therapy(therapy, dose)
        sim.step()
    return sim.summary()


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--days", type=float, default=450)
    ap.add_argument("--therapy-day", type=float, default=250)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--neoantigens", type=float, default=0.1)
    ap.add_argument("--therapy", default="trojan+niche_guardian+memory")
    ap.add_argument("--dose", type=int, default=300)
    a = ap.parse_args(argv)

    print(f"REVERSIÓN · {a.therapy} (dosis {a.dose}) · {a.scenario} · {a.days:g} días · {a.seeds} semillas")
    print("La guardiana es la pieza inventada del proyecto. Esto mide cuánto tiene que durar.\n")
    print(
        f"{'duración del estado':<32}{'herencia':>10}{'peligrosas':>12}{'reinas':>8}"
        f"{'revertidas':>12}   por semilla"
    )
    for inheritance in (True, False):
        for tag, rev in LEVELS:
            t0 = time.time()
            res = [
                run(a.therapy, a.dose, a.scenario, s, a.days, a.therapy_day, rev, inheritance, a.neoantigens)
                for s in range(a.seeds)
            ]
            m = {k: sum(r[k] for r in res) / len(res) for k in ("dangerous", "queens", "guardians_reverted")}
            dang = sorted(r["dangerous"] for r in res)
            print(
                f"{tag:<32}{'sí' if inheritance else 'NO':>10}{m['dangerous']:>12.0f}"
                f"{m['queens']:>8.0f}{m['guardians_reverted']:>12.0f}"
                f"   {dang}   ({time.time() - t0:.0f}s)",
                flush=True,
            )
        print()


if __name__ == "__main__":
    main()
