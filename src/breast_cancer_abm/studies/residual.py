"""El escenario real: quitar la masa primero y pelear por el residuo.

Este banco existe por un resultado, no por una idea. Todo lo construido en este proyecto
—el caballo de Troya, la corrección que se propaga, la vigilancia residente, la autoridad
del tejido— actúa sobre la SUPERFICIE del tumor:

  acceso en 3D ................. 26 % accesible, 1 % del núcleo alcanzado
  macrófagos ................... no atraviesan la corriente de salida
  autoridad del tejido ......... 148 células con contexto, tenga el tumor 1.089 o 8.832

La superficie crece con r², el volumen con r³. Por eso existe un umbral (~4.000 células
en este modelo) a partir del cual el tumor gana por geometría, sin necesidad de evolucionar.

La cirugía no es elegante, pero es la única forma conocida de convertir un problema de
volumen en uno de superficie. La pregunta que contesta este banco no es si la cirugía
funciona —eso ya se sabe—, sino si el arsenal de este proyecto sirve para lo que de verdad
importa después: la ENFERMEDAD RESIDUAL y la recaída.

    bcabm study residual --seeds 3
"""

from __future__ import annotations

import argparse
import time

from breast_cancer_abm.model.biology import SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation

# (etiqueta, ¿opera?, terapia tras la cirugía, propagación de la corrección)
ARMS = (
    ("nothing", False, None, 0.0),
    ("solo el arsenal (sin operar)", False, "trojan+niche_guardian+memory", 0.0),
    ("solo cirugía", True, None, 0.0),
    ("cirugía + vigilancia", True, "memory", 0.0),
    ("cirugía + arsenal", True, "trojan+niche_guardian+memory", 0.0),
    ("cirugía + arsenal + propagación", True, "trojan+niche_guardian+memory", 0.15),
)


def run(
    arm, scenario: str, seed: int, days: float, day: float, dose: int, coverage: float, neoantigens: float
) -> dict:
    _, operates, therapy, prop = arm
    p = Parameters(
        prob_neoantigen=neoantigens,
        correction_reach=8 if prop else 0,
        correction_spread=prop,
        editing_error=0.01,
    )
    sim = Simulation(scenario, seed, p)
    step_0 = int(day * STEPS_PER_DAY)
    for step in range(int(days * STEPS_PER_DAY)):
        if step == step_0 and operates:
            sim.resect(coverage)
        # La terapia va DESPUÉS de operar, como la adyuvancia real.
        if step == step_0 + (4 * STEPS_PER_DAY if operates else 0) and therapy:
            sim.release_therapy(therapy, dose)
        sim.step()
    return sim.summary()


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--days", type=float, default=500)
    ap.add_argument("--day", type=float, default=250, help="día de la cirugía")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--neoantigens", type=float, default=0.1)
    ap.add_argument("--dose", type=int, default=300)
    ap.add_argument(
        "--coverage", type=float, default=0.95, help="fracción de la masa que se lleva el bloque extirpado"
    )
    a = ap.parse_args(argv)

    print(
        f"RESIDUAL · {a.scenario} · {a.days:g} días · cirugía el día {a.day:g} "
        f"(cobertura {a.coverage:.0%}) · terapia 4 días después · {a.seeds} semillas"
    )
    print("La pregunta: ¿sirve nuestro arsenal para la recaída, que es el problema real?\n")
    print(f"{'brazo':<34}{'peligrosas':>12}{'reinas':>8}{'extirp.':>9}{'sanas ext.':>12}   por semilla")
    for arm in ARMS:
        t0 = time.time()
        res = [
            run(arm, a.scenario, s, a.days, a.day, a.dose, a.coverage, a.neoantigens) for s in range(a.seeds)
        ]
        m = {
            k: sum(r[k] for r in res) / len(res)
            for k in ("dangerous", "queens", "excised", "healthy_excised")
        }
        dang = sorted(r["dangerous"] for r in res)
        print(
            f"{arm[0]:<34}{m['dangerous']:>12.0f}{m['queens']:>8.0f}"
            f"{m['excised']:>9.0f}{m['healthy_excised']:>12.0f}"
            f"   {dang}   ({time.time() - t0:.0f}s)",
            flush=True,
        )


if __name__ == "__main__":
    main()
