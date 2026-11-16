"""Las cinco maneras de atacar la malla: ¿cuál abre el paso sin soltar al tumor?

La malla desmoplásica hace DOS cosas a la vez, y son opuestas para nosotros:

  ESTORBA   genera la presión intersticial que impide entrar a cualquier terapia.
  CONTIENE  frena la invasión. Es, en parte, lo que mantiene al tumor donde está.

Por eso en la clínica todo lo que la DESTRUYE ha fracasado (PEGPH20 en fase III, 2019;
eliminar CAFs empeoró los tumores en ratones) y lo que la ABLANDA o normaliza funciona
(losartán con quimio, Murphy 2019; normalización vascular, Jain).

Este banco mide las dos caras por separado, para cada intervención:

  acceso     cuántas fagocitosis consiguen los macrófagos (que ya están ahí de serie)
             y cuántos pasos se bloquean contra la corriente de salida
  contención cuántas células invasivas hay al final y cuántas invasiones frenó la malla

Una intervención solo sirve si sube el acceso SIN disparar las invasivas.

    bcabm study mesh --seeds 3 --days 450
"""

from __future__ import annotations

import argparse
import time

from breast_cancer_abm.model.biology import SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation

MODES = ("none", "hyaluronidase", "anti_lox", "losartan", "anti_tgfb", "kill_caf")

EXPLANATION = {
    "none": "sin tocar la malla",
    "hyaluronidase": "degrada el HA (PEGPH20): quita material, los CAFs siguen",
    "anti_lox": "destensa sin retirar (simtuzumab)",
    "losartan": "reduce la producción (Murphy 2019)",
    "anti_tgfb": "corta la activación de nuevos CAFs",
    "kill_caf": "elimina al fabricante (lo que empeoró en ratones)",
}


def run(
    mode: str,
    scenario: str,
    seed: int,
    days: float,
    day: float,
    params: Parameters,
    intensity: float,
    chain: bool,
) -> dict:
    sim = Simulation(scenario, seed, params)
    step_0 = int(day * STEPS_PER_DAY)
    for step in range(int(days * STEPS_PER_DAY)):
        if step == step_0 and mode != "none":
            sim.treat_matrix(mode, intensity)
        if chain:
            if step == step_0 + int(5 * STEPS_PER_DAY):
                sim.block_cd47(targeted=True)
            elif step == step_0 + int(15 * STEPS_PER_DAY):
                sim.damage_telomeres()
            elif step == step_0 + int(20 * STEPS_PER_DAY):
                sim.release_therapy("niche_guardian+memory")
        sim.step()
    r = sim.summary()
    r["mesh"] = float(sim.Ms.max())
    return r


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--days", type=float, default=450)
    ap.add_argument("--day", type=float, default=250, help="día en que se trata la malla")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--neoantigens", type=float, default=0.1)
    ap.add_argument("--intensity", type=float, default=0.6)
    ap.add_argument("--chain", action="store_true", help="además de la malla, aplica la cadena completa")
    ap.add_argument("--modes", nargs="*", default=list(MODES))
    a = ap.parse_args(argv)
    params = Parameters(prob_neoantigen=a.neoantigens)

    print(
        f"MALLA · {a.scenario} · {a.days:g} días · intervención el día {a.day:g} al "
        f"{a.intensity:.0%} · {a.seeds} semillas{' · + cadena completa' if a.chain else ''}"
    )
    print("sirve solo si sube ACCESO sin disparar INVASIVAS\n")
    print(
        f"{'intervención':<16}{'peligros':>10}{'invasivas':>11}{'comidas':>9}"
        f"{'bloqueados':>12}{'frenadas':>10}{'malla':>8}   qué es"
    )
    for mode in a.modes:
        t0 = time.time()
        res = [run(mode, a.scenario, s, a.days, a.day, params, a.intensity, a.chain) for s in range(a.seeds)]
        m = {
            k: sum(r.get(k, 0) for r in res) / len(res)
            for k in (
                "dangerous",
                "invasive_cells",
                "macrophage_eaten",
                "steps_blocked",
                "invasions_blocked",
                "mesh",
            )
        }
        print(
            f"{mode:<16}{m['dangerous']:>10.0f}{m['invasive_cells']:>11.0f}"
            f"{m['macrophage_eaten']:>9.0f}{m['steps_blocked']:>12.0f}"
            f"{m['invasions_blocked']:>10.0f}{m['mesh']:>8.2f}   {EXPLANATION[mode]}"
        )
        danger = sorted(r["dangerous"] for r in res)
        print(f"{'':<16}mediana {danger[len(danger) // 2]:>7} · {danger}   ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
