"""La cadena completa: abrir · desenmascarar · castrar · que se lo coman · cerrar la fuente.

La idea, en los términos en que la planteó el usuario: si al tumor le quitas la careta
(CD47) y además no puede duplicarse (telómeros), se convierte en objetivo único y el
cuerpo se lo come solo. Lo natural, sin cirugía y sin quimio.

El modelo añade dos pasos que la versión original no tenía, y los dos salieron de medir:

  PASO 0 (abrir)   sin esto no llega nadie. Con presión intersticial los macrófagos
                   fagocitaron 3 células; sin ella, 920. Es el cuello de botella real.
  PASO 4 (cerrar)  comerse 920 células bajó las peligrosas de 957 a 932. Limpiar la masa
                   no cura si las células madre tumorales siguen reponiéndola.

Los brazos de control no están de adorno: sirven para saber si la cadena funciona por lo
que creemos que funciona, o si en realidad la hace sola una de sus piezas.

    bcabm study chain --seeds 3 --days 500
"""

from __future__ import annotations

import argparse
import time

from breast_cancer_abm.model.biology import SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import Simulation

# Cada paso es (días después del inicio del tratamiento, qué se hace).
OPEN = "open"  # hialuronidasa + normalizar vasos: baja la presión intersticial
MASK = "cd47"  # biespecífico anti-CD47: quita el "no me comas" solo donde hay sospecha
STERILIZE = "6tiodg"  # 6-tio-dG: rompe el telómero YA en las células con telomerasa
SOURCE = "guardian"  # cierra la fuente: cura las reinas y deja memoria residente

ARMS = {
    "nothing": (),
    "solo abrir": ((0, OPEN),),
    "solo guardiana+memoria": ((20, SOURCE),),
    "cadena SIN abrir": ((5, MASK), (15, STERILIZE), (20, SOURCE)),
    "cadena SIN cerrar": ((0, OPEN), (5, MASK), (15, STERILIZE)),
    "CADENA COMPLETA": ((0, OPEN), (5, MASK), (15, STERILIZE), (20, SOURCE)),
}

COLUMNS = (
    ("dangerous", "peligrosas"),
    ("queens", "reinas"),
    ("macrophage_eaten", "comidas"),
    ("macrophage_eaten_healthy", "comidas SANAS"),
    ("collateral", "colateral"),
    ("steps_blocked", "bloqueados"),
)


def run(
    arm, scenario: str, seed: int, days: float, therapy_day: float, params: Parameters, gate: str
) -> dict:
    sim = Simulation(scenario, seed, params)
    schedule = {int((therapy_day + d) * STEPS_PER_DAY): step for d, step in arm}
    for step in range(int(days * STEPS_PER_DAY)):
        act = schedule.get(step)
        if act == OPEN:
            sim.degrade_matrix(0.5)
            sim.treat_vessels("normalize")
        elif act == MASK:
            sim.block_cd47(targeted=True)
        elif act == STERILIZE:
            sim.damage_telomeres()
        elif act == SOURCE:
            sim.release_therapy(gate)
        sim.step()
    return sim.summary()


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--days", type=float, default=500)
    ap.add_argument("--therapy-day", type=float, default=250)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--neoantigens", type=float, default=0.1)
    ap.add_argument("--gate", default="niche_guardian+memory", help="terapia del paso 'cerrar la fuente'")
    ap.add_argument("--arms", nargs="*", default=list(ARMS))
    a = ap.parse_args(argv)
    params = Parameters(prob_neoantigen=a.neoantigens)

    print(
        f"CADENA · {a.scenario} · {a.days:g} días · inicio el día {a.therapy_day:g} · "
        f"{a.seeds} semillas · neoantígenos {a.neoantigens:g} · fuente: {a.gate}"
    )
    print("abrir(0) - careta(+5) - castrar(+15) - cerrar la fuente(+20)")
    print("peligrosas = in situ + invasivas + reinas · entre paréntesis, el rango de semillas\n")
    print(f"{'brazo':<24}" + "".join(f"{e:>16}" for _, e in COLUMNS))
    for name in a.arms:
        t0 = time.time()
        res = [run(ARMS[name], a.scenario, s, a.days, a.therapy_day, params, a.gate) for s in range(a.seeds)]
        row = f"{name:<24}"
        for key, _ in COLUMNS:
            v = [r.get(key, 0) for r in res]
            row += f"{sum(v) / len(v):>16.0f}"
        print(row + f"   ({time.time() - t0:.0f}s)")
        danger = sorted(r["dangerous"] for r in res)
        median = danger[len(danger) // 2]
        print(f"{'':<24}mediana {median:>8} · por semilla {danger}")


if __name__ == "__main__":
    main()
