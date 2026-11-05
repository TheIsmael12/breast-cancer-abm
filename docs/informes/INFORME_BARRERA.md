# El cuello de botella no es identificar al tumor: es llegar hasta él

Resultado principal de la sesión del 2026-10-20. Cambia la lectura de todo lo anterior.

## Qué se añadió al modelo

Hasta ahora el simulador dejaba entrar a las células terapéuticas sin coste. Se han
añadido las tres barreras reales de un tumor sólido, que no son la misma cosa:

| Barrera | Qué es | Cómo se levanta |
|---|---|---|
| Geometría | una esfera esconde su interior (accesible ∝ 1/r) | no se arregla con fármacos |
| **Presión intersticial** | vasos porosos + linfáticos colapsados: el líquido empuja hacia fuera | degradar matriz, normalizar vasos |
| Oxígeno | difunde ~100-150 µm desde un capilar: el centro está lejos | acercar vasos ("reconducir") |

No hay ninguna membrana que envuelva el tumor. Lo que impide entrar es **corriente en
contra**, no una pared: por eso en el modelo la presión solo frena a quien intenta
ENTRAR, y salir es gratis (`Simulacion._frena_la_presion`).

También se añadieron **macrófagos, CD47 y calreticulina**, porque la fagocitosis es una
**puerta Y natural**: el macrófago necesita que falte el freno ("no me comas", CD47) Y
que esté el acelerador ("cómeme", calreticulina). Es la misma estructura lógica que
`INFORME_IDENTIFICADORES.md` proponía construir a mano.

## El número que lo decide

Mismo tumor, mismo anti-CD47, misma biología. Lo único que cambia es la presión:

| presión | células con "cómeme" | fagocitadas | pasos bloqueados |
|---|---|---|---|
| 1,0 (real) | 621 | **3** | 17.139 |
| 0,0 (irreal) | 559 | **920** | 0 |

**×300.** La señal estaba puesta en 620 células todo el tiempo. Lo que fallaba era llegar.

## Lo que esto le hace al trabajo anterior

`guardiana+memoria`, que era la estrategia del proyecto, **deja de funcionar**:

| | antes (sin barrera) | ahora (con barrera) |
|---|---|---|
| guardiana+memoria | 487 peligrosas | **7.230** (sin tratar: 7.286) |

Y acumula 189.188 pasos bloqueados frente a 63.318 del tejido sin tratar: son sus propias
células terapéuticas chocando contra la corriente de salida sin entrar nunca.

**Todo lo medido antes del 2026-10-20 asumía un tumor al que se entraba gratis.**

## La cadena de cinco tiempos

Secuencia propuesta: abrir (matriz + normalizar vasos) → anti-CD47 dirigido (biespecífico)
→ 6-tio-dG → macrófagos → cerrar la fuente (`guardiana_nicho+memoria`). Ver `cadena.py`.

| brazo | media | por semilla |
|---|---|---|
| nada | 7.286 | 976 · 10.487 · 10.394 |
| solo abrir | 4.302 | 981 · 10.924 · 1.002 |
| solo guardiana+memoria | 7.230 | 955 · 10.457 · 10.279 |
| cadena SIN abrir | 9.258 | 7.178 · 10.473 · 10.122 |
| cadena SIN cerrar | **10.828** | 10.760 · 10.914 · 10.811 |
| CADENA COMPLETA | 5.658 | **2** · 10.922 · 6.050 |

Los dos controles son lo valioso:

- **SIN abrir (9.258) es peor que no tratar.** Desenmascarar y castrar sin poder entrar
  es gastar munición fuera.
- **SIN cerrar (10.828) es el peor brazo de todos.** Se comió 423 células y acabó peor
  que el tumor sin tratar: limpiar la masa sin cerrar la fuente no es neutro, es
  contraproducente, porque las células madre tumorales reponen mientras tanto.

La puerta Y sí protegió al tejido sano: **27 células sanas con el CD47 bloqueado, 0
comidas**. En el modelo la secuencia no produce anemia, porque la célula sana no expone
calreticulina. (En una persona sí la habría: los glóbulos rojos viejos sí la exponen.
De ahí la dosis de cebado, el anti-SIRPα y los biespecíficos.)

## El techo: hasta qué tamaño funciona

Misma cadena, 7 semillas, ordenadas por carga tumoral el día del tratamiento:

| tumor el día 250 | peligrosas al final | |
|---|---|---|
| 213 | 0 | ✔ |
| 471 | 5 | ✔ |
| 861 | 6.050 | ✘ |
| 917 | 2 | ✔ |
| 4.451 | 10.866 | ✘ |
| 7.256 | 10.935 | ✘ |
| 10.079 | 10.922 | ✘ |

**Por encima de ~4.000 células no funciona nunca** (3 de 3). Por debajo de ~1.000 funciona
en 3 de cada 4. La diferencia entre curar y no curar no fue el tratamiento: fue el tamaño
del tumor el día que se empezó.

## Qué significa y qué no

Esto es un modelo educativo: rejilla de 110×110, tasas de mutación aceleradas más de mil
veces, 8 genes de ~20.000. **"Funciona en el simulador" no es un resultado médico.** La
prueba está en esta misma sesión: la estrategia que llevaba dos días dándose por buena se
cayó al añadir un solo mecanismo que faltaba.

Lo que sí queda es una pregunta de investigación bien planteada, y es la que habría que
llevar a un laboratorio:

> El cuello de botella de las terapias celulares contra tumores sólidos, ¿es el
> reconocimiento de la diana o el transporte hasta ella? Y si es el transporte, ¿hasta
> qué carga tumoral sigue siendo reversible?

La segunda pregunta pendiente sigue siendo la de `HOJA_DE_RUTA.md` (A7): con qué facilidad
se desdiferencia una célula tumoral de mama *in vivo*. Las dos son baratas comparadas con
cualquier otra cosa de la fase B.

## Aviso de método

La varianza entre semillas domina al tratamiento: la misma configuración da 976 o 10.487.
**Ninguna media de 3 semillas significa nada aquí.** Usar medianas y 5-7 semillas como
mínimo; `cadena.py` ya informa ambas.
