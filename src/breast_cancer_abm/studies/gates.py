"""¿Qué puerta lógica identifica mejor a la célula tumoral y SOLO a ella?

INFORME_IDENTIFICADORES.md concluyó, sobre 1.084 pacientes reales, que ningún marcador
tiene a la vez cobertura y exclusividad, y que la salida son las puertas Y. Esto lo
comprueba DENTRO del simulador: deja crecer un tumor sin tratar y, en ese tejido, mide
para cada puerta dos números que son los de un test diagnóstico:

  SENSIBILIDAD  de las células tumorales, ¿a cuántas ve?        (si falla, el tumor rebrota)
  ESPECIFICIDAD de las células SANAS, ¿a cuántas respeta?       (si falla, hay daño colateral)

No simula tratamiento: solo pregunta a cada puerta qué haría con cada célula del tejido.

    bcabm study gates --days 400 --seeds 3
"""

from __future__ import annotations

import argparse

from breast_cancer_abm.model.biology import SCENARIOS, STEPS_PER_DAY, Parameters
from breast_cancer_abm.model.tissue import LUMINAL, Simulation

GATES = ("simple", "double", "niche", "tert", "tert+niche")
LABEL = {
    "simple": "1 señal (actividad)",
    "double": "Y: tumoral + actividad",
    "niche": "Y: madre + fuera del nicho",
    "tert": "Y: telomerasa + aneuploidía",
    "tert+niche": "las dos Y en paralelo",
}


def evaluate(sim: Simulation) -> dict[str, tuple[int, int, int, int]]:
    """Por puerta: (tumorales vistas, tumorales totales, sanas tocadas, sanas totales)."""
    r = {p: [0, 0, 0, 0] for p in GATES}
    sim.p = sim.p.__class__(**{**sim.p.__dict__, "internal_read_error": 0.0})
    for pos, c in sim.cells.items():
        if c.kind != LUMINAL or c.necrosis or c.guardian:
            continue
        o = sim._observe(c, pos, 0, 0)
        is_tumoral = c.n_drivers > 0
        for p in GATES:
            d = r[p]
            d[1 if is_tumoral else 3] += 1
            if sim._gate(p, c, o):
                d[0 if is_tumoral else 2] += 1
    return {p: tuple(v) for p, v in r.items()}


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenario", choices=SCENARIOS, default="luminal")
    ap.add_argument(
        "--days",
        nargs="*",
        type=float,
        default=[120, 200, 300],
        help="momentos en los que se mide (el tejido sigue creciendo entre ellos)",
    )
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--neoantigens", type=float, default=0.1)
    a = ap.parse_args(argv)

    days = sorted(a.days)
    total = {d: {p: [0, 0, 0, 0] for p in GATES} for d in days}
    for seed in range(a.seeds):
        sim = Simulation(a.scenario, seed, Parameters(prob_neoantigen=a.neoantigens))
        step = 0
        for d in days:
            while step < int(d * STEPS_PER_DAY):
                sim.step()
                step += 1
            for p, v in evaluate(sim).items():
                for i in range(4):
                    total[d][p][i] += v[i]

    print(f"Escenario {a.scenario} · sin tratar · {a.seeds} semillas · neoantígenos {a.neoantigens:g}")
    print("sensibilidad = tumorales que ve · especificidad = sanas que respeta")
    for d in days:
        t = total[d][GATES[0]]
        print(f"\n— día {d:g}: {t[1]} tumorales · {t[3]} luminales sanas")
        print(f"  {'puerta lógica':<30}{'sensib.':>10}{'especif.':>11}{'sanas tocadas':>15}")
        for p in GATES:
            seen_set, tum_count, touched, healthy = total[d][p]
            sens = seen_set / tum_count if tum_count else float("nan")
            spec = 1 - touched / healthy if healthy else float("nan")
            print(f"  {LABEL[p]:<30}{sens:>9.1%}{spec:>11.2%}{touched:>15}")


if __name__ == "__main__":
    main()
