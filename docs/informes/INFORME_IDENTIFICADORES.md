# ¿Qué identifica a una célula de cáncer de mama?

Datos: 1084 tumores de TCGA PanCancer Atlas. TMB mediana 1.30 mutaciones/Mb.

| Marcador | Cobertura | ¿Lo tiene alguna célula SANA? | Estado |
|---|---|---|---|
| Telomerasa (TERT) activa | 90% | Células madre: médula ósea, intestino, piel, folículo piloso y mama | Inhibidores probados; fracasaron solos y dan toxicidad medular |
| Neoantígenos (proteínas nuevas por mutación) | 3% | NINGUNA: no existen en una célula sana. Exclusividad perfecta | Base de las vacunas de ARNm personalizadas (en ensayo) |
| HER2 muy amplificado | 14% | Corazón y epitelios tienen algo de HER2: es diferencia de CANTIDAD, no de tipo | Trastuzumab y T-DXd aprobados; cardiotoxicidad por esa misma razón |
| Pérdida de MHC-I (HLA/B2M) | 1% | Algunas células sanas lo bajan; las NK ya lo vigilan de serie | Lo explotan las terapias con células NK |
| Marcadores de célula madre (CD44+/CD24-/ALDH1+) | 100% | Las células madre SANAS de la mama los llevan igual | Sin terapia aprobada que los use como diana única |
| Aneuploidía (número de cromosomas alterado) | 95% | Prácticamente ninguna célula sana. Casi exclusivo | INVISIBLE desde fuera: está dentro del núcleo, no en la superficie |

## El problema, en una frase

Ningún marcador tiene a la vez COBERTURA alta y EXCLUSIVIDAD total:

- Los **exclusivos** (neoantígenos) solo están en una minoría: 3% de los tumores de mama tienen carga mutacional alta.
- Los **universales** (telomerasa, marcadores de célula madre, aneuploidía) los comparte el tejido sano regenerativo, o no se ven desde fuera de la célula.

## La salida: combinaciones

Ninguna célula sana cumple VARIAS condiciones a la vez. Ejemplos de combinación con cobertura alta y, en teoría, sin equivalente sano:

- **telomerasa activa  Y  marcador de célula madre  Y  está en el conducto, no en el nicho basal**  
  la madre sana vive en un sitio concreto; la tumoral no respeta la arquitectura
- **telomerasa activa  Y  aneuploidía**  
  una célula madre sana tiene telomerasa pero NO tiene los cromosomas descuadrados
- **HER2 muy alto  Y  proliferando**  
  el corazón tiene algo de HER2 pero no se divide: base de por qué el T-DXd funciona
- **pérdida de MHC-I  Y  marcadores de estrés (MICA/B)**  
  es exactamente lo que ya hacen las células NK del cuerpo

## Lo que esto implica para la estrategia

1. No busques el marcador perfecto: **no existe**. Es un problema abierto de investigación, no un hueco de este proyecto.
2. La especificidad tiene que venir de **combinar condiciones** (puerta AND), que es lo que ya hacen los receptores synNotch en los ensayos reales.
3. La combinación más prometedora para mama no usa neoantígenos (son escasos) sino **telomerasa + algo que la célula madre sana no tenga**: aneuploidía o estar fuera de su sitio.
4. Y sigue pendiente el problema que destapó el 3D: aunque la identifiques bien, hay que **llegar** hasta ella.