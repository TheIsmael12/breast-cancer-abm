"""Corrección que se propaga: ¿cuánta fidelidad necesita para merecer la pena?

El cuello de botella medido en este proyecto es que cada célula corregida cuesta un
agente que la encuentre, la secuencie y acierte con todos sus genes rotos. En 450 días
solo se consiguen ~23. Es coste lineal.

La idea: que una célula ya corregida corrija a sus vecinas, y que las que quedan limpias
corrijan a su vez. Pasa a crecimiento exponencial, y basta conseguirlo una vez.

Equivalentes reales a medias:
  - EFECTO ESPECTADOR: existe en terapia génica. El gen suicida HSV-TK pasa el tóxico a
    las vecinas por las uniones comunicantes, sin necesidad de transducirlas todas.
  - GENE DRIVE: un elemento de CRISPR que se copia a sí mismo en la población. Funciona
    en insectos; por eso mismo está sujeto a moratorias.

Y el peligro, que es la razón de ser de este banco: si algo se copia solo, **los errores
también se copian**. La edición génica tiene errores fuera de diana. Una corrección que
se auto-propaga con un 1 % de error va sembrando mutaciones nuevas justo donde más
células hay. Por eso aquí no se pregunta "¿funciona?", sino:

    ¿qué fidelidad MÍNIMA necesita la edición para que propagarla sea mejor que no hacerlo?

Eso sí es una especificación que un laboratorio puede medir.

    bcabm study propagation --seeds 3
"""

from __future__ import annotations

import argparse
import time

from breast_cancer_abm.model.biology import SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation

# (propagación por paso, error fuera de diana, alcance en vecinas)
CASES = (
    (0.00, 0.00, 0),  # control: sin propagación, como hasta ahora
    (0.02, 0.00, 8),
    (0.05, 0.00, 8),
    (0.15, 0.00, 8),
    (0.05, 0.01, 8),  # 1 % de error fuera de diana
    (0.05, 0.05, 8),  # 5 %
    (0.15, 0.05, 8),  # propagación alta Y error alto: el peor caso
)


def run(
    prop: float,
    err: float,
    reach: int,
    scenario: str,
    seed: int,
    days: float,
    therapy_day: float,
    therapy: str,
    dose: int,
    neoantigens: float,
) -> dict:
    p = Parameters(
        prob_neoantigen=neoantigens, correction_reach=reach, correction_spread=prop, editing_error=err
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

    print(f"PROPAGACIÓN · {a.therapy} (dosis {a.dose}) · {a.scenario} · {a.days:g} días · {a.seeds} semillas")
    print("La pregunta no es si funciona, sino con cuánto error fuera de diana deja de compensar.\n")
    print(
        f"{'propag.':>9}{'error':>7}{'alcance':>9}{'peligrosas':>12}{'reinas':>8}"
        f"{'propagadas':>12}{'limpias':>9}{'errores':>9}   por semilla"
    )
    for prop, err, reach in CASES:
        t0 = time.time()
        res = [
            run(prop, err, reach, a.scenario, s, a.days, a.therapy_day, a.therapy, a.dose, a.neoantigens)
            for s in range(a.seeds)
        ]
        m = {
            k: sum(r[k] for r in res) / len(res)
            for k in ("dangerous", "queens", "propagated_corrections", "fully_corrected", "editing_errors")
        }
        dang = sorted(r["dangerous"] for r in res)
        print(
            f"{prop:>9.2f}{err:>7.2f}{reach:>9}{m['dangerous']:>12.0f}{m['queens']:>8.0f}"
            f"{m['propagated_corrections']:>12.0f}{m['fully_corrected']:>9.0f}"
            f"{m['editing_errors']:>9.0f}   {dang}   ({time.time() - t0:.0f}s)",
            flush=True,
        )


if __name__ == "__main__":
    main()
