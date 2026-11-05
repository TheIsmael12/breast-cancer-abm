# Contraste con pacientes reales (TCGA PanCancer Atlas, 1084 tumores de mama)

## 1. ¿Están bien elegidos los genes del modelo?

| Gen | Pacientes reales | Tipo |
|---|---|---|
| PIK3CA ★ |  33.9 % | oncogén |
| TP53 ★ |  32.6 % | supresor |
| MYC ★ |  15.0 % | oncogén |
| CCND1 |  14.9 % | oncogén |
| ERBB2 ★ |  13.6 % | oncogén |
| CDH1 ★ |  13.1 % | supresor |
| GATA3 |  11.7 % | supresor |
| FGFR1 |  11.3 % | oncogén |
| PTEN ★ |  10.3 % | supresor |
| KMT2C |  10.1 % | supresor |
| MAP3K1 |   9.2 % | supresor |
| MAP2K4 |   6.3 % | supresor |
| RB1 |   5.8 % | supresor |
| NCOR1 |   5.1 % | supresor |
| FOXA1 |   4.3 % | oncogén |
| RUNX1 |   4.3 % | supresor |
| CDKN2A |   4.2 % | supresor |
| ARID1A |   4.1 % | supresor |
| BRCA2 ★ |   4.0 % | supresor |
| NF1 |   3.9 % | supresor |
| AKT1 |   3.6 % | oncogén |
| CBFB |   3.5 % | supresor |
| PIK3R1 |   3.5 % | supresor |
| TBX3 |   3.2 % | supresor |
| BRCA1 ★ |   3.0 % | supresor |

★ = gen incluido en el simulador.

## 2. ¿Podrían las guardianas reparar el bug de cada paciente?

- Con al menos un error en los genes del modelo: ** 74.4 %**
- Con TODOS sus errores dentro de los genes del modelo (reparables con nuestros vectores): ** 22.3 %**
- Sin ningún driver de esta lista (su bug está en otra parte):  10.4 %
- Errores driver por paciente: mediana 2, máximo 15; con 3 o más:  43.3 %

Ojo: estos son solo ~30 genes; el genoma de un tumor de mama real tiene además miles de mutaciones pasajeras y reordenamientos.

## 3. ¿Podría la puerta Y reconocer al tumor sin equivocarse?

La guardiana necesita una señal EXCLUSIVA del tumor: neoantígenos (dependen de la carga
mutacional, TMB) o pérdida de MHC-I.

- TMB mediana: **1.30 mutaciones/Mb** (melanoma o pulmón: ~10-20)
- Tumores con TMB alta (>= 10/Mb): **  2.7 %**
- Con mutación o deleción en HLA-A/B/C o B2M:   1.1 %

| Subtipo | Pacientes | TMB mediana | TP53 | PIK3CA | ERBB2 |
|---|---|---|---|---|---|
| LumA | 499 | 1.03 |  11.0 % |  48.1 % |   9.2 % |
| LumB | 197 | 1.67 |  37.6 % |  31.0 % |  13.2 % |
| Her2 | 78 | 2.42 |  71.8 % |  37.2 % |  73.1 % |
| Basal | 171 | 2.37 |  89.5 % |  13.5 % |   1.8 % |
| Normal | 36 | 1.28 |  30.6 % |  22.2 % |  13.9 % |

## Conclusión automática

- La puerta Y del simulador da por hecho que todo tumor muestra neoantígenos. En la realidad, la mayoría de los cánceres de mama tienen pocos (TMB mediana 1.3/Mb): **es el cuello de botella nº 1**. Hace falta otra señal exclusiva (p. ej. una mutación concreta como PIK3CA H1047R presentada en HLA, o combinaciones de proteínas de superficie).
- Reparar gen a gen solo cubriría por completo a ~22 % de las pacientes con nuestro catálogo: **cada paciente necesitaría su propio juego de vectores** (medicina personalizada), igual que las vacunas de ARNm personalizadas.
- ERBB2/HER2 amplificado es la señal más limpia (proteína de superficie en exceso) y ya hay fármacos reales contra ella (trastuzumab, T-DXd).