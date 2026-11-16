"""Banco de pruebas sin ventana: compara terapias con varias semillas.

bcabm study experiment --scenario luminal --days 200 --therapy-day 120 --seeds 3
bcabm study experiment --scenario her2 --plasticity 0   # colmena sin plasticidad
bcabm study experiment --therapies guardian+memory --surveillance --hormonal
"""

from __future__ import annotations

import argparse
import time

from breast_cancer_abm.model.biology import CALIBRATIONS, SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.policies import THERAPIES
from breast_cancer_abm.model.tissue import Simulation

COLUMNS = (
    "tumor_cells",
    "dangerous",
    "queens",
    "senescent_cells",
    "macrophage_eaten",
    "macrophage_eaten_healthy",
    "cd47_healthy",
    "collateral",
)


def run(
    scenario: str,
    therapy: str | None,
    seed: int,
    days: float,
    therapy_day: float,
    params: Parameters,
    vessels: str = "normal",
    surveillance: bool = False,
    hormonal: bool = False,
    sterilize: bool = False,
    tiodg: bool = False,
    tiodg_day: float | None = None,
    anti_cd47: bool = False,
    cd47_day: float | None = None,
    sterilize_day: float | None = None,
) -> dict:
    sim = Simulation(scenario, seed, params)
    therapy_step = int(therapy_day * STEPS_PER_DAY)
    for step in range(int(days * STEPS_PER_DAY)):
        if therapy and step == therapy_step:
            sim.release_therapy(therapy)
        if vessels != "normal" and step == therapy_step:
            sim.treat_vessels(vessels)
        if anti_cd47 and step == int((cd47_day if cd47_day is not None else therapy_day) * STEPS_PER_DAY):
            sim.block_cd47()
        if sterilize and step == int(
            (sterilize_day if sterilize_day is not None else therapy_day) * STEPS_PER_DAY
        ):
            sim.inhibit_telomerase()
        if tiodg and step == int((tiodg_day if tiodg_day is not None else therapy_day) * STEPS_PER_DAY):
            sim.damage_telomeres()
        if hormonal and step == therapy_step:
            sim.p.estrogen = 0.1  # tamoxifeno / inhibidor de aromatasa
            sim._event("PREVENCIÓN HORMONAL: estrógenos bloqueados")
        if surveillance and step > therapy_step and step % (30 * STEPS_PER_DAY) == 0:
            # Biopsia líquida mensual: si detecta ADN tumoral (>= 5 células peligrosas), dosis pequeña.
            if sim.summary()["dangerous"] >= 5:
                sim.release_therapy("guardian", 10)
        sim.step()
    return sim.summary()


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument("--days", type=float, default=200)
    ap.add_argument("--therapy-day", type=float, default=120)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument(
        "--therapies",
        nargs="*",
        default=["none", *THERAPIES],
        help="admite combinaciones con '+', p. ej. guardiana+memoria",
    )
    ap.add_argument(
        "--plasticity",
        type=float,
        default=Parameters.plasticity,
        help="probabilidad por paso de que una obrera tumoral se vuelva reina",
    )
    ap.add_argument(
        "--vessels", default="normal", choices=("normal", "antiangiogenic", "normalize", "starve", "reroute")
    )
    ap.add_argument("--surveillance", action="store_true", help="biopsia líquida mensual + dosis si positiva")
    ap.add_argument("--hormonal", action="store_true", help="bloqueo de estrógenos desde el día de terapia")
    ap.add_argument(
        "--sterilize",
        action="store_true",
        help="inhibidor de telomerasa el día de la terapia (bloquea el contador de divisiones)",
    )
    ap.add_argument(
        "--anti-cd47",
        action="store_true",
        help="bloquea CD47 ('no me comas'); afecta también a las células sanas",
    )
    ap.add_argument(
        "--cd47-day", type=float, default=None, help="día del anti-CD47; por defecto, el mismo que la terapia"
    )
    ap.add_argument(
        "--sterilize-day",
        type=float,
        default=None,
        help="día del inhibidor de telomerasa; por defecto, el mismo que la terapia",
    )
    ap.add_argument(
        "--tiodg",
        action="store_true",
        help="6-tio-dG: rompe el telómero al instante en células con telomerasa",
    )
    ap.add_argument(
        "--drug-p53", type=float, default=0.0, help="cualidad: el fármaco también restaura p53 (0-1)"
    )
    ap.add_argument(
        "--drug-bait",
        action="store_true",
        help="cualidad: la carga se queda dentro; quien se coma la célula se envenena",
    )
    ap.add_argument(
        "--drug-anti-sasp",
        type=float,
        default=0.0,
        help="cualidad: silencia el SASP (0-1), como un inhibidor de JAK",
    )
    ap.add_argument(
        "--drug-persistent",
        type=float,
        default=0.0,
        help="cualidad: sigue alcanzando a las células que nacen después (0-1)",
    )
    ap.add_argument(
        "--tiodg-day",
        type=float,
        default=None,
        help="día del 6-tio-dG; por defecto, el mismo que la terapia celular",
    )
    ap.add_argument(
        "--neoantigens",
        type=float,
        default=Parameters.prob_neoantigen,
        help="1 = optimista; 0.1 = calibrado con pacientes reales (TCGA)",
    )
    ap.add_argument(
        "--calibration",
        choices=list(CALIBRATIONS),
        default="treatment",
        help="tasa de mutación. 'tratamiento' (por defecto) da tumores siempre "
        "y rápido; las otras bajan la incidencia y hacen estudiable la "
        "prevención, a costa de necesitar más semillas y más días. "
        "Ver la tabla en biologia.py",
    )
    a = ap.parse_args(argv)
    params = Parameters(
        plasticity=a.plasticity,
        prob_neoantigen=a.neoantigens,
        mutation_rate=CALIBRATIONS[a.calibration],
        drug_p53=a.drug_p53,
        drug_bait=a.drug_bait,
        drug_anti_sasp=a.drug_anti_sasp,
        drug_persistent=a.drug_persistent,
    )

    print(
        f"Escenario {a.scenario} · {a.days:g} días · terapia el día {a.therapy_day:g} · "
        f"{a.seeds} semillas (media) · plasticidad {a.plasticity:g} · neoantígenos {a.neoantigens:g} · "
        f"vasos {a.vessels}"
        f"{' · vigilancia' if a.surveillance else ''}{' · hormonal' if a.hormonal else ''}"
        f"{' · ANTI-CD47' if a.anti_cd47 else ''}"
        f"{' · ANTITELOMERASA' if a.sterilize else ''}{f' · 6-TIO-dG dia {a.tiodg_day or a.therapy_day:g}' if a.tiodg else ''}"
    )
    print("peligrosas = in situ + invasivas + reinas (lo que de verdad es cáncer)\n")
    print(f"{'terapia':<20}" + "".join(f"{c[:12]:>13}" for c in COLUMNS))
    for therapy in a.therapies:
        t0 = time.time()
        res = [
            run(
                a.scenario,
                None if therapy == "none" else therapy,
                s,
                a.days,
                a.therapy_day,
                params,
                a.vessels,
                a.surveillance,
                a.hormonal,
                a.sterilize,
                a.tiodg,
                a.tiodg_day,
                a.anti_cd47,
                a.cd47_day,
                a.sterilize_day,
            )
            for s in range(a.seeds)
        ]
        means = {c: sum(r.get(c, 0) for r in res) / len(res) for c in COLUMNS}
        print(
            f"{therapy:<20}"
            + "".join(f"{means[c]:>13.1f}" for c in COLUMNS)
            + f"   ({time.time() - t0:.0f}s)"
        )


if __name__ == "__main__":
    main()
