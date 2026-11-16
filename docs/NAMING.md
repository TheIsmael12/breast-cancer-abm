# Equivalencia de nombres (español → inglés)

El 2026-11-16 el código pasó a inglés y a un paquete de Python (`src/breast_cancer_abm`).
La bitácora, los resultados brutos y los issues anteriores usan los nombres antiguos; esta
tabla los traduce. La traducción se hizo por *tokens* con un diccionario único y se comprobó
con 32 simulaciones deterministas: el modelo da **exactamente** los mismos números.

## Ficheros y comandos

| Antes | Ahora | Comando |
|---|---|---|
| `agentes.py` | `src/breast_cancer_abm/model/agents.py` |  |
| `barrera.py` | `src/breast_cancer_abm/studies/barrier.py` | `bcabm study barrier` |
| `biologia.py` | `src/breast_cancer_abm/model/biology.py` |  |
| `cadena.py` | `src/breast_cancer_abm/studies/chain.py` | `bcabm study chain` |
| `capturas.py` | `src/breast_cancer_abm/viz/figures.py` | `bcabm figures` |
| `datos_reales.py` | `src/breast_cancer_abm/studies/real_data.py` | `bcabm study real-data` |
| `diseminacion.py` | `src/breast_cancer_abm/model/dissemination.py` | `bcabm study dissemination` |
| `experimento.py` | `src/breast_cancer_abm/studies/experiment.py` | `bcabm study experiment` |
| `identificadores.py` | `src/breast_cancer_abm/studies/identifiers.py` | `bcabm study identifiers` |
| `main.py` | `src/breast_cancer_abm/cli.py` | `bcabm viewer` |
| `malla.py` | `src/breast_cancer_abm/studies/mesh.py` | `bcabm study mesh` |
| `mi_celula.py` | `src/breast_cancer_abm/model/my_cell.py` |  |
| `politicas.py` | `src/breast_cancer_abm/model/policies.py` |  |
| `propagacion.py` | `src/breast_cancer_abm/studies/propagation.py` | `bcabm study propagation` |
| `prueba_residuo.py` | `src/breast_cancer_abm/studies/residual_test.py` | `bcabm study residual-test` |
| `puertas.py` | `src/breast_cancer_abm/studies/gates.py` | `bcabm study gates` |
| `residual.py` | `src/breast_cancer_abm/studies/residual.py` | `bcabm study residual` |
| `reversion.py` | `src/breast_cancer_abm/studies/reversion.py` | `bcabm study reversion` |
| `sin_cirugia.py` | `src/breast_cancer_abm/studies/no_surgery.py` | `bcabm study no-surgery` |
| `sin_cirugia_3d.py` | `src/breast_cancer_abm/studies/no_surgery_3d.py` | `bcabm study no-surgery-3d` |
| `tejido.py` | `src/breast_cancer_abm/model/tissue.py` |  |
| `tejido3d.py` | `src/breast_cancer_abm/model/tissue3d.py` | `bcabm study tissue3d` |
| `validacion.py` | `src/breast_cancer_abm/studies/validation.py` | `bcabm study validation` |
| `visor.py` | `src/breast_cancer_abm/viz/viewer.py` |  |
| `tests/test_humo.py` | `tests/test_smoke.py` | `pytest` |
| `docs/avances/`, `docs/informes/`, `docs/capturas/` | `docs/progress/`, `docs/reports/`, `docs/figures/` | |
| `docs/BITACORA.md` | `docs/LOGBOOK.md` | |
| `datos/`, `resultados/` | `data/`, `results/` | |
| `requirements.txt` | `pyproject.toml` | `pip install -e ".[app]"` |

## Clases

| Antes | Ahora |
|---|---|
| `Agente` | `Agent` |
| `Analizar` | `Analyze` |
| `Atacar` | `Attack` |
| `Celula` | `Cell` |
| `Escenario` | `Scenario` |
| `Esperar` | `Wait` |
| `Gen` | `Gene` |
| `Mover` | `Move` |
| `Observacion` | `Observation` |
| `Organismo` | `Organism` |
| `Parametros` | `Parameters` |
| `ParamsDiseminacion` | `DisseminationParams` |
| `Percepcion` | `Perception` |
| `Reparar` | `Repair` |
| `Resultado` | `Result` |
| `Simulacion` | `Simulation` |
| `Simulacion3D` | `Simulation3D` |
| `Visor` | `Viewer` |

## Funciones y métodos públicos

| Antes | Ahora |
|---|---|
| `Celula.actualizar_fenotipo` | `Cell.update_phenotype` |
| `Organismo.biopsia` | `Organism.biopsy` |
| `Organismo.paso` | `Organism.step` |
| `Organismo.resumen` | `Organism.summary` |
| `Resultado.imprimir` | `Result.print_table` |
| `Resultado.veredicto` | `Result.verdict` |
| `Simulacion.ablacionar` | `Simulation.ablate` |
| `Simulacion.biopsia` | `Simulation.biopsy` |
| `Simulacion.bloquear_cd47` | `Simulation.block_cd47` |
| `Simulacion.bloquear_pd1` | `Simulation.block_pd1` |
| `Simulacion.ctdna` | `Simulation.ctdna` |
| `Simulacion.danar_telomeros` | `Simulation.damage_telomeres` |
| `Simulacion.degradar_matriz` | `Simulation.degrade_matrix` |
| `Simulacion.describir` | `Simulation.describe` |
| `Simulacion.dia` | `Simulation.day` |
| `Simulacion.estado` | `Simulation.state` |
| `Simulacion.focos_carcinoma` | `Simulation.carcinoma_foci` |
| `Simulacion.frecuencias_tumor` | `Simulation.tumor_frequencies` |
| `Simulacion.inhibir_telomerasa` | `Simulation.inhibit_telomerase` |
| `Simulacion.liberar_adc` | `Simulation.release_adc` |
| `Simulacion.liberar_terapia` | `Simulation.release_therapy` |
| `Simulacion.paso` | `Simulation.step` |
| `Simulacion.resecar` | `Simulation.resect` |
| `Simulacion.restaurar_arquitectura` | `Simulation.restore_architecture` |
| `Simulacion.resumen` | `Simulation.summary` |
| `Simulacion.tratar_matriz` | `Simulation.treat_matrix` |
| `Simulacion.tratar_vasos` | `Simulation.treat_vessels` |
| `Simulacion3D.ablacionar` | `Simulation3D.ablate` |
| `Simulacion3D.ctdna` | `Simulation3D.ctdna` |
| `Simulacion3D.danar_telomeros` | `Simulation3D.damage_telomeres` |
| `Simulacion3D.degradar_matriz` | `Simulation3D.degrade_matrix` |
| `Simulacion3D.dia` | `Simulation3D.day` |
| `Simulacion3D.es_carcinoma` | `Simulation3D.is_carcinoma` |
| `Simulacion3D.focos_carcinoma` | `Simulation3D.carcinoma_foci` |
| `Simulacion3D.fraccion_accesible` | `Simulation3D.accessible_fraction` |
| `Simulacion3D.liberar_terapia` | `Simulation3D.release_therapy` |
| `Simulacion3D.paso` | `Simulation3D.step` |
| `Simulacion3D.resumen` | `Simulation3D.summary` |
| `Simulacion3D.tratar_vasos` | `Simulation3D.treat_vessels` |
| `Visor.mostrar` | `Viewer.show` |
| `agentes.adyacente` | `agents.adjacent` |
| `agentes.cazar` | `agents.hunt` |
| `agentes.deambular` | `agents.wander` |
| `agentes.distancia` | `agents.distance` |
| `agentes.hacia` | `agents.toward` |
| `agentes.politica_castrador` | `agents.sterilizer_policy` |
| `agentes.politica_linfocito_t` | `agents.t_cell_policy` |
| `agentes.politica_macrofago` | `agents.macrophage_policy` |
| `agentes.politica_nk` | `agents.nk_policy` |
| `agentes.politica_vehiculo` | `agents.vehicle_policy` |
| `politicas.car_her2` | `policies.car_her2` |
| `politicas.cazador` | `policies.hunter` |
| `politicas.centinela` | `policies.sentinel` |
| `politicas.hibrido` | `policies.hybrid` |
| `politicas.mata_reinas` | `policies.queen_killer` |
| `politicas.propagador` | `policies.propagator` |
| `politicas.quimio` | `policies.chemo` |
| `politicas.reina` | `policies.queen` |
| `politicas.reina_solo` | `policies.queen_only` |
| `politicas.reparador` | `policies.repairer` |

## Parámetros del modelo (`Parametros` → `Parameters`)

| Antes | Ahora |
|---|---|
| `tamano` | `grid_size` |
| `radio_conducto` | `duct_radius` |
| `tasa_mutacion` | `mutation_rate` |
| `factor_hrd` | `hrd_factor` |
| `prob_dano_hrd` | `prob_hrd_damage` |
| `prob_muerte_normal` | `prob_normal_death` |
| `estrogeno` | `estrogen` |
| `o2_hipoxia` | `o2_hypoxia` |
| `reclutamiento_inmune` | `immune_recruitment` |
| `reclutamiento_macrofagos` | `macrophage_recruitment` |
| `max_macrofagos` | `max_macrophages` |
| `eficacia_anti_cd47` | `anti_cd47_efficacy` |
| `cd47_alcanza_sanas` | `cd47_reaches_healthy` |
| `duracion_anti_cd47` | `anti_cd47_duration` |
| `fraccion_macrofagos_m2` | `m2_macrophage_fraction` |
| `max_inmunes_naturales` | `max_natural_immune` |
| `eficacia_reparacion` | `repair_efficacy` |
| `prob_muerte_ataque` | `prob_attack_kill` |
| `dias_desenganche_pd1` | `pd1_disengage_days` |
| `fraccion_madre_normal` | `normal_stem_fraction` |
| `prob_division_simetrica` | `prob_symmetric_division` |
| `divisiones_obrera` | `worker_divisions` |
| `eficacia_antitelomerasa` | `antitelomerase_efficacy` |
| `eficacia_6tiodg` | `efficacy_6tiodg` |
| `prob_sasp_plasticidad` | `prob_sasp_plasticity` |
| `eficacia_cebo` | `bait_efficacy` |
| `prob_cebo_encadena` | `prob_bait_chain` |
| `farmaco_p53` | `drug_p53` |
| `farmaco_cebo` | `drug_bait` |
| `farmaco_anti_sasp` | `drug_anti_sasp` |
| `farmaco_persistente` | `drug_persistent` |
| `prob_muerte_diferenciada` | `prob_differentiated_death` |
| `factor_quiescencia_madre` | `stem_quiescence_factor` |
| `plasticidad` | `plasticity` |
| `prob_neoantigeno` | `prob_neoantigen` |
| `fuerza_contexto` | `context_strength` |
| `contexto_necesita_matriz` | `context_needs_matrix` |
| `propagacion_correccion` | `correction_spread` |
| `error_edicion` | `editing_error` |
| `alcance_correccion` | `correction_reach` |
| `reversion_guardiana` | `guardian_reversion` |
| `herencia_guardiana` | `guardian_inheritance` |
| `prob_fagocitosis` | `prob_phagocytosis` |
| `prob_aneuploidia` | `prob_aneuploidy` |
| `prob_aneuploidia_p53ok` | `prob_aneuploidy_p53ok` |
| `coste_lectura_interna` | `internal_read_cost` |
| `error_lectura_interna` | `internal_read_error` |
| `produccion_matriz` | `matrix_production` |
| `recambio_matriz` | `matrix_turnover` |
| `activacion_caf_tgfb` | `caf_activation_tgfb` |
| `activacion_caf_rigidez` | `caf_activation_stiffness` |
| `contencion_matriz` | `matrix_containment` |
| `presion_intersticial` | `interstitial_pressure` |
| `efecto_presion` | `pressure_effect` |
| `caida_presion_matriz` | `pressure_drop_matrix` |
| `caida_presion_normalizar` | `pressure_drop_normalize` |
| `o2_vaso_tumoral` | `o2_tumor_vessel` |
| `o2_vaso_asfixiado` | `o2_starved_vessel` |
| `inmunidad_adaptativa` | `adaptive_immunity` |
| `divisiones_programa` | `program_divisions` |
| `contraccion_t` | `t_contraction` |
| `dias_contraccion` | `contraction_days` |
| `max_t_clonales` | `max_clonal_t` |
| `radio_ifng` | `ifng_radius` |
| `dias_pdl1_inducido` | `induced_pdl1_days` |
| `umbral_agotamiento` | `exhaustion_threshold` |
| `eficacia_tex_prog` | `tex_prog_efficacy` |
| `eficacia_tex_term` | `tex_term_efficacy` |
| `autorrenovacion_tex` | `tex_self_renewal` |
| `dias_vida_tex_term` | `tex_term_lifespan_days` |
| `prob_revigorizacion` | `prob_reinvigoration` |
| `divisiones_revigorizacion` | `reinvigoration_divisions` |

## Terapias (`TERAPIAS` → `THERAPIES`)

| Antes | Ahora |
|---|---|
| `cazador` | `hunter` |
| `car_her2` | `car_her2` |
| `reparador` | `repairer` |
| `hibrido` | `hybrid` |
| `centinela` | `sentinel` |
| `reina` | `queen` |
| `reina_solo` | `queen_only` |
| `mata_reinas` | `queen_killer` |
| `guardiana` | `guardian` |
| `guardiana_simple` | `simple_guardian` |
| `guardiana_nicho` | `niche_guardian` |
| `guardiana_tert` | `tert_guardian` |
| `guardiana_tert_nicho` | `tert_niche_guardian` |
| `propagador` | `propagator` |
| `cebo` | `bait` |
| `memoria_tert` | `tert_memory` |
| `memoria_tert_nicho` | `tert_niche_memory` |
| `memoria` | `memory` |
| `troyano` | `trojan` |
| `troyano_cebo` | `trojan_bait` |
| `troyano_memoria` | `trojan_memory` |
| `castrador` | `sterilizer` |
| `troya` | `trojan_horse` |
| `troya_sin_o2` | `trojan_horse_no_o2` |
| `quimio` | `chemo` |
| `mi_celula` | `my_cell` |

## Opciones de línea de comandos

| Antes | Ahora |
|---|---|
| `--brazos` | `--arms` |
| `--cadena` | `--chain` |
| `--calibracion` | `--calibration` |
| `--castrar` | `--sterilize` |
| `--cobertura` | `--coverage` |
| `--dia` | `--day` |
| `--dia-castrar` | `--sterilize-day` |
| `--dia-cd47` | `--cd47-day` |
| `--dia-terapia` | `--therapy-day` |
| `--dia-tiodg` | `--tiodg-day` |
| `--dias` | `--days` |
| `--dosis` | `--dose` |
| `--escenario` | `--scenario` |
| `--espera-biopsia` | `--biopsy-wait` |
| `--farmaco-anti-sasp` | `--drug-anti-sasp` |
| `--farmaco-cebo` | `--drug-bait` |
| `--farmaco-p53` | `--drug-p53` |
| `--farmaco-persistente` | `--drug-persistent` |
| `--intensidad` | `--intensity` |
| `--matriz` | `--matrix` |
| `--modos` | `--modes` |
| `--neoantigenos` | `--neoantigens` |
| `--plasticidad` | `--plasticity` |
| `--primera` | `--first` |
| `--procesos` | `--processes` |
| `--puerta` | `--gate` |
| `--refrescar` | `--refresh` |
| `--semilla` | `--seed` |
| `--semillas` | `--seeds` |
| `--solo` | `--only` |
| `--terapia` | `--therapy` |
| `--terapias` | `--therapies` |
| `--umbrales` | `--thresholds` |
| `--vasos` | `--vessels` |
| `--via` | `--route` |
| `--vigilancia` | `--surveillance` |

## Claves del resumen y otras

Las claves de `resumen()` (ahora `summary()`) y de los contadores siguen el mismo diccionario:
`peligrosas` → `dangerous`, `tumorales` → `tumor_cells`, `reinas` → `queens`, `sanas` → `healthy`,
`colateral` → `collateral`, `reparaciones` → `repairs`, `guardianas` → `guardians`, etc.
