# Modelo de diseminación: la célula que se escapó

> **Por qué existe este módulo.** El simulador original modela **un conducto mamario**.
> Cuatro días de experimentos terminaron siempre en el mismo sitio: lo que mata no es el
> tumor que se extirpa, es **la célula que ya se había escapado antes de operar** — lejos
> de cualquier vaso, sin dividirse, y en otro órgano. Ese problema **no se puede estudiar
> en un modelo de un solo compartimento**, por construcción. De ahí este módulo.

## Regla del proyecto: cada fase se valida antes de construir encima

No se usa la palabra «perfecto», porque no se puede comprobar. Cada fase tiene
**criterios de aceptación explícitos, escritos antes de programar**, que se pasan o no se
pasan. Y cada fase tiene un **dominio de validez escrito**: para qué sirve y para qué no.

El motivo es concreto: el 2026-10-22 se descubrió que las tasas de mutación aceleradas
invalidaban tres líneas de trabajo a la vez. El modelo llevaba días pareciendo correcto.
Se detectó porque un test falló, no porque el modelo fuera bueno. **Un modelo sin
limitaciones conocidas no es mejor: es uno cuyas limitaciones no se han encontrado aún.**

---

## FASE 1 — Siembra y dormancia. Criterios de aceptación

Nada de esto se le programa al modelo. Se deja que el tumor crezca, invada y siembre por
su cuenta, y **después** se mide.

| # | criterio | diana | fuente |
|---|---|---|---|
| **D1** | fracción de pacientes con DTC en médula **al diagnóstico** | **30-40 %** | Braun et al., *NEJM* 2005, n=4.703 |
| **D2** | fracción de DTC **quiescentes** (Ki67−) | **>80 %** | Pantel & Brakenhoff, *Nat Rev Cancer* 2004 |
| **D3** | la **mayoría** de las células que se sueltan deben morir en el trayecto | ineficiencia metastásica | literatura clásica de cascada metastásica |
| **D4** | el resultado **no debe depender** de un parámetro inventado sin respaldo | control explícito | lección del 2026-10-23 (ciclo celular del ADC) |

**Si no pasa los cuatro, la fase 1 no está lista y no se construye la fase 2.**

D4 merece explicación: se comprueba corriendo el modelo con la suposición dudosa activada
y desactivada. Si el desenlace cambia mucho, el modelo no está prediciendo — está
devolviendo lo que se asumió. Se hizo así con el ADC y la dependencia del ciclo celular, y
resultó que no cambiaba nada (978 muertes frente a 13.173 dieron el mismo tumor final).

### D4 · criterio de aceptación (escrito el 2026-10-28, ANTES de ver los resultados)

Se ejecuta con `python diseminacion.py --d4` (32 semillas por defecto, las mismas en
todas las variantes; se paraleliza por semillas).

**Qué se varía.** Los parámetros sin fuente de `ParamsDiseminacion`:
`competencia_temprana`, `prob_sobrevive_transito`, `prob_entra_dormida` y
`prob_muere_en_nicho`. Cada uno a ×0,5 y ×2 de su valor por defecto, recortado a [0, 1].
Excepción fijada de antemano: en `prob_entra_dormida` (0,85) doblar no tiene sentido, así
que se varía la fracción complementaria (1 − p = 0,15) a ×0,5 y ×2 → p = 0,925 y 0,70.

**Qué se mide.** D1, D2 y D3 en cada variante, con el mismo veredicto PASA/FALLA que la
validación normal (D1 en 25-50 %, D2 > 80 %, D3 > 90 %).

**Regla.** Para cada supuesto y cada criterio Dk se compara el veredicto de las variantes
con el de la base (mismas 32 semillas):
- **se mantiene** si el veredicto es el mismo en las dos variantes;
- **cambia** si en alguna variante cambia y la diferencia con la base supera el ruido
  (para D1, más de 2 errores estándar binomiales de la diferencia,
  2·√(2·p·(1−p)/n), unos ±24 puntos con n = 32 y p = 0,4; D2 y D3 se comparan tal cual);
- **no concluyente** si cambia pero dentro de ese ruido.

**D4 PASA** si ningún supuesto cambia el veredicto de D1-D3. **D4 FALLA** si alguno lo
cambia por encima del ruido. **D4 NO CONCLUYENTE** si solo hay cambios dentro del ruido.

**Dependencias por construcción, declaradas antes.** D3 es casi exactamente
1 − `prob_sobrevive_transito`, y D2 es `prob_entra_dormida` menos los despertares. Que el
desenlace se mueva con ellas es esperable; la pregunta es si el **veredicto** cambia. Si
cambia, **cuenta como fallo de D4** igualmente: significa que ese criterio se cumple
porque se supuso, no porque el modelo lo produzca. La explicación no lo excusa.

**Predicción a priori.** D3 no debería cambiar (0,96 y 0,99 siguen por encima de 0,90).
D2 con `prob_entra_dormida` = 0,70 debería caer por debajo del 80 %, y D1 debería subir
mucho con `competencia_temprana` = 0,5. Si se cumple, D4 falla.

---

## Lo que la fase 1 debe contener, y con qué respaldo

### Los tres compartimentos

| compartimento | contenido | estado |
|---|---|---|
| **mama** | el tumor primario | ya existe (`tejido.py`) |
| **circulación** | células en tránsito, la mayoría muere | nuevo |
| **médula ósea** | el nicho donde esperan años | nuevo |

### Los procesos, y lo que se sabe de cada uno

**Intravasación.** La célula invasiva entra en un vaso. El modelo ya tiene invasión y
rotura de la membrana mioepitelial, así que hay de dónde partir: una célula `invasiva` que
alcanza un vaso puede soltarse.

**Dormancia, y esto tiene mecanismo medido — no hay que inventarlo:**
- Los **microvasos estables mantienen dormidas** a las células, vía **trombospondina-1**
  endotelial.
- Los **brotes vasculares las despiertan**, vía **TGF-β1 y periostina**.
  (Ghajar et al., *Nat Cell Biol* 2013;15:807-817)
- Nicho perisinusoidal rico en **E-selectina y SDF-1**; movilizarlas las sensibiliza
  (Price et al., *Sci Transl Med* 2016)

**Lo que NO hay que modelar como protección:** la quiescencia **no** protege del ataque de
células T. Medido: CAR-T eliminaron el **93 %** de las dormidas, igual que de las
proliferantes (Goddard et al., *Cancer Cell* 2024). Lo que las protege es la **escasez
relativa**: son poquísimas y hay poquísimos linfocitos que las reconozcan, así que
raramente se encuentran. **Esa es la propiedad que el modelo debe reproducir**, y es
espacial/estadística, no bioquímica.

**Detección — el drama clínico no es que estén, es que no se ven:**
- biopsia de médula por inmunohistoquímica: detecta el **24 %**
- por citometría: detecta el **47 %** en la misma población
- ADN tumoral en sangre: cuando se positiviza, **el 55-72 % ya tiene metástasis visibles**

El modelo debe permitir «mirar» con sensibilidad imperfecta, porque una terapia perfecta
sobre una población indetectable no sirve de nada.

---

## Dominio de validez (se rellena al cerrar cada fase)

**Fase 1 · v1** — medido el 2026-10-28 (`python diseminacion.py --d4`, semillas 1-32, día 250).
**Estado: NO LISTA. D4 FALLA.** (Sustituida por la v2, más abajo; se conserva como registro.)

| variante | D1 (25-50 %) | D2 (>80 %) | D3 (>90 %) |
|---|---|---|---|
| **base** | 38 % | 100 % | 97,7 % |
| competencia_temprana 0,125 / 0,5 | 41 % / 41 % | 100 % / 100 % | 97,8 % / 97,3 % |
| prob_sobrevive_transito 0,01 / 0,04 | 25 % / 50 % | **65 %** / 88 % | 99,0 % / 96,2 % |
| prob_entra_dormida 0,925 / 0,70 | 38 % / 38 % | 100 % / **58 %** | 97,8 % / 97,9 % |
| prob_muere_en_nicho 0,00075 / 0,003 | 38 % / 34 % | 100 % / 100 % | 97,7 % / 97,8 % |

Lo que dice, sin adornos:

- **D2 no es un resultado del modelo, es el supuesto.** Con `prob_entra_dormida` = 0,70
  la fracción de quiescentes cae al 58 % y el criterio falla. Estaba predicho antes de
  correr, y precisamente por eso es un fallo: D2 se cumple porque se asumió que la
  mayoría entra dormida. Hasta que la fracción de dormidas salga de un mecanismo (el
  nicho de Ghajar: microvaso estable frente a brote), D2 no valida nada.
- **D2 también se cae con `prob_sobrevive_transito` a la mitad (65 %).** Aquí el
  mecanismo es otro y más incómodo: con menos siembra hay 8 pacientes con DTC, casi todos
  con 1-2 células, y la mediana de una fracción tan discreta (0, ½, 1) salta mucho. Es
  probablemente ruido de la métrica, pero la regla escrita antes decía que D2 se compara
  tal cual, así que cuenta como cambio. Además la propia base da D2 = 100 %: la mediana
  por paciente está saturada y es un mal estimador; habría que medir D2 sobre el total
  de DTC o con más células por paciente. Eso queda como limitación, no se cambia ahora.
- **D3 es robusto, pero por construcción**: es casi exactamente 1 − `prob_sobrevive_transito`.
  Pasa en todas las variantes y no demuestra nada sobre el modelo.
- **D1 mantiene el veredicto en todas las variantes**, pero en `prob_sobrevive_transito`
  recorre todo el rango de aceptación (25 % → 50 %): está en el borde, no holgado.
- **Sorpresa: `competencia_temprana` casi no mueve D1** (38 % → 41 % con 0,125 y con
  0,5). La vía temprana, que se introdujo precisamente para explicar la fracción
  intermedia de pacientes con DTC, apenas contribuye a ella a día 250: D1 lo decide la
  vía clásica. La justificación del docstring de `ParamsDiseminacion` no queda respaldada
  por esta medida.
- `prob_muere_en_nicho` (sin fuente) no cambia ningún veredicto en ×0,5-×2.

**Sobre D1 y el ruido.** Con 8 semillas la validación daba D1 = 50 % (4/8): dentro del
rango de aceptación 25-50 % pero **fuera de la diana real 30-40 %** de Braun 2005, y justo
en el borde. Con 32 semillas la base da **38 %** (12/32), dentro de la diana. El 50 % era
ruido de 8 semillas (un caso = 12,5 puntos); el error estándar con 32 sigue siendo de
unos ±9 puntos. Recuérdese además que D1 es **calibración**, no validación (ver la nota
al final de `diseminacion.py`).

**Para qué sirve hoy la fase 1:** genera una población de DTC con una fracción de
pacientes sembrados plausible y con la ineficiencia del tránsito impuesta. **Para qué no
sirve:** para afirmar nada sobre cuántas DTC están dormidas (es el supuesto devuelto), ni
para atribuir la siembra a la vía temprana. La fase 2 no debe empezar hasta que D2 salga
de un mecanismo y no de `prob_entra_dormida`.

---

## FASE 1 · v2 — la dormancia sale del nicho vascular, no de un parámetro

> Diseño, supuestos y criterio **escritos el 2026-10-28, ANTES de ejecutar la v2.**
> Motivo: la v1 falló D4 porque D2 era `prob_entra_dormida` devuelto (ver arriba).
> `prob_entra_dormida` y `prob_brote_vascular` se **eliminan**.

### Mecanismo

1. **Geometría de la médula.** `n_vasos` segmentos de microvaso. Cada uno está en uno de
   dos estados: **estable** (sinusoide/capilar maduro, tipo L) o **en brote**
   (endotelio angiogénico, tipo H). Recambio de dos estados por segmento:
   estable → brote con tasa *k*₊, brote → estable con tasa *k*₋ = 1/`duracion_brote_dias`,
   y *k*₊ se fija para que la fracción estacionaria en brote sea `fraccion_vasos_brote`:
   *k*₊ = *f*·*k*₋/(1 − *f*). El estado inicial de cada segmento se sortea de esa
   distribución estacionaria (la médula adulta ya está en equilibrio al empezar).
2. **Alojamiento.** La DTC que sobrevive al tránsito se aloja en un segmento elegido con
   probabilidad proporcional a su peso: 1 si está en brote, `sesgo_homing_estable` si es
   estable (base 1 = sin preferencia, es decir, proporcional al volumen vascular).
3. **Estado emergente.** La DTC está **dormida si y solo si su vaso es estable**
   (trombospondina-1 endotelial) y **despierta (Ki67+) si su vaso está en brote**
   (TGF-β1, periostina) — Ghajar et al., *Nat Cell Biol* 2013. Cuando el vaso cambia de
   estado, la DTC cambia con él, en los dos sentidos: la dormancia es un estado
   **reversible** dictado por el entorno (Sosa, Bragado y Aguirre-Ghiso,
   *Nat Rev Cancer* 2014). En la fase 1 la DTC despierta no se divide.
4. **Muerte en el nicho**: `prob_muere_en_nicho` por paso, igual que en la v1.
5. **Aleatoriedad separada.** La médula usa su propio generador (derivado de la semilla),
   así que variar un parámetro del nicho **no cambia el tumor ni la siembra**: las
   variantes de D4 de la médula comparan las mismas DTC en distinto nicho.

### Parámetros nuevos y su respaldo

| parámetro | base | respaldo | ¿D4? |
|---|---|---|---|
| `fraccion_vasos_brote` | **0,16** | **MEDIDO**: vasos tipo H (CD31^hi Emcn^hi, el endotelio angiogénico del hueso) sobre el volumen vascular total en fémur proximal humano: **18,9 %** a los 45-59 años, **13,9 %** a los 60-75 (Wang et al., *Cell Death Dis* 2017, n = 39; Kusumbe et al., *Nature* 2014 para la definición). Base = punto medio, porque la mediana de edad al diagnóstico de cáncer de mama ronda los 60. | sí, **dentro del rango medido**: 0,139 y 0,189 |
| `duracion_brote_dias` | 14 | **SUPUESTO.** No he encontrado una medida de la vida de un segmento tipo H o de un brote en médula adulta. Orden de magnitud: maduración de un brote (pericitos) de días a semanas; el endotelio normal se recambia en 47-23.000 días (Hobson & Denekamp, *Br J Cancer* 1984), lo que acota por arriba. | sí, ×0,5 y ×2 |
| `sesgo_homing_estable` | 1,0 | **SUPUESTO en la magnitud.** Price et al. (*Sci Transl Med* 2016) sitúan a las DTC en el nicho perisinusoidal E-selectina/SDF-1 (vaso estable), lo que apunta a un sesgo > 1, pero no hay un cociente medido. La base 1 es la opción conservadora: **no ayuda a D2**. | sí, ×0,5 y ×2 (0,5 y 2) |
| `n_vasos` | 500 | **SUPUESTO** de discretización: controla cuántas DTC comparten vaso (y por tanto se despiertan a la vez). | sí, ×0,5 y ×2 |
| reversibilidad | sí | Sosa et al. 2014 (cualitativo). | **no** como parámetro; se corre aparte como control estructural (abajo) |

### Lo que este diseño implica, dicho antes de correr

**D2 tampoco aquí es libre.** Con estado reversible y alojamiento sin sesgo, la
fracción esperada de DTC dormidas en equilibrio es **1 − `fraccion_vasos_brote`**
= 0,84 (0,81-0,86 en el rango medido). La diferencia con la v1 no es que D2 deje de
depender de un número, sino **de qué número depende**: antes era una probabilidad
inventada para que saliera >80 %; ahora es una **medida independiente** (fracción de
endotelio angiogénico en hueso humano, histología vascular) que no se tomó de ningún
estudio de DTC. Si D2 pasa, lo que se valida es la **consistencia de dos medidas
independientes a través del mecanismo de Ghajar**, no el modelo de tumor. Y pasa por
poco: 84 % frente a una diana >80 %. Con 0,189 la predicción (81 %) está en el borde.

**Control estructural (informativo, no cuenta para D4 porque no es un parámetro, pero se
publica):** con despertar **irreversible** (la DTC que ve un brote queda despierta para
siempre) la fracción de dormidas cae con el tiempo de residencia,
≈ 0,84 · e^(−*k*₊·*t*), y con *k*₊ ≈ 0,014/día y residencias de decenas a cientos de
días **D2 debería hundirse muy por debajo del 80 %**. Es decir: el resultado depende
de que la dormancia sea reversible, y eso solo tiene respaldo cualitativo.

### Métrica de D2: cambio anunciado antes de correr

La v1 medía D2 como **mediana por paciente** de la fracción de dormidas. Con 1-2 DTC por
paciente esa fracción solo vale 0, ½ o 1, la mediana se satura (la base daba 100 %) y
salta con una sola célula (el 65 % de `prob_sobrevive_transito` ×0,5 en la v1).
**Nueva métrica primaria: fracción agregada** = DTC dormidas de todos los pacientes /
DTC totales. Justificación: es como se reporta la Ki67 en DTC (porcentaje de células
positivas sobre las DTC detectadas en la serie), y no está acotada a {0, ½, 1}.
Limitación: pesa más a los pacientes con muchas DTC. La mediana por paciente se sigue
imprimiendo como dato secundario, sin veredicto.

Este cambio **no se hace para aprobar**: con el mecanismo nuevo la mediana por paciente
también rondaría el 84-100 %. Se hace porque la métrica vieja no distingue 84 % de 100 %.

### D4 v2 · criterio de aceptación (escrito ANTES de ejecutar)

- **Qué se varía** (`SUPUESTOS_D4`): `competencia_temprana`, `prob_sobrevive_transito`,
  `prob_muere_en_nicho`, `duracion_brote_dias`, `sesgo_homing_estable` y `n_vasos`, cada
  uno a ×0,5 y ×2 (recortado a [0, 1] si es probabilidad); y `fraccion_vasos_brote` a
  los dos extremos de su rango medido (0,139 y 0,189). Mismas 32 semillas en todo.
- **Qué se mide.** D1 (25-50 %), D2 agregada (> 80 %), D3 (> 90 %).
- **Ruido.** D1: 2·√(2·p(1−p)/n) como en la v1. **D2: ahora también tiene banda de
  ruido**, porque ya no es una mediana discreta: 2·√2·EE, con EE el error estándar
  bootstrap de la fracción agregada **remuestreando pacientes** (1.000 remuestreos,
  generador fijo), calculado en la base. D3 se compara tal cual (sigue siendo ≈ 1 −
  `prob_sobrevive_transito`, dependencia por construcción ya declarada).
- **Regla** idéntica a la v1: por supuesto y criterio, «se mantiene» / «CAMBIA» (veredicto
  distinto y diferencia fuera del ruido) / «no concluyente» (distinto pero dentro).
  **D4 PASA** si nada cambia, **FALLA** si algo CAMBIA, **NO CONCLUYENTE** si solo hay
  cambios dentro del ruido.
- Un cambio por construcción cuenta como fallo, igual que antes.

**Predicción a priori.** `duracion_brote_dias`, `n_vasos` y `prob_muere_en_nicho` no
deberían mover D2 (la fracción estacionaria no depende de ellos). `sesgo_homing_estable`
a 0,5 baja la fracción dormida **al llegar** a ≈ 1 − 0,16/(0,16 + 0,5·0,84) ≈ 72 %, pero
como el estado sigue al vaso, ese sesgo se olvida con el tiempo de relajación del
recambio, τ = 1/(*k*₊ + *k*₋) = `duracion_brote_dias`·(1 − *f*) ≈ 12 días. Así que D2
quedará entre 72 % y 84 % según cuántas DTC hayan llegado en las últimas ~2 semanas;
como el tumor crece y siembra cada vez más, puede ser una fracción grande. **Posible
FALLA de D4 ahí**, y no lo sé de antemano. Lo mismo, más suave, con
`duracion_brote_dias` ×2 combinado con la siembra tardía. `fraccion_vasos_brote` = 0,189 deja D2 en ≈ 81 %, en el borde: posible «no
concluyente». D3 igual que en la v1.

Nota: `tejido.py` lo está modificando otro trabajo en paralelo (última modificación
2026-10-28 08:40); la base de la v2 puede no coincidir con la de la v1 en D1 por eso.

### Resultados v2 (medido el 2026-10-28, después de escribir lo anterior)

`python diseminacion.py --d4`, semillas 1-32, día 250. **Estado: NO LISTA.
D2 FALLA en la base y D4 es NO CONCLUYENTE.**

| variante | D1 (25-50 %) | D2 agregada (>80 %) | D3 (>90 %) | DTC totales |
|---|---|---|---|---|
| **base** | 34 % | **75,0 %** ✗ | 98,3 % | 24 |
| competencia_temprana 0,125 / 0,5 | 38 % / 38 % | 73,1 % / 85,7 % | 98,3 % / 98,0 % | 26 / 28 |
| prob_sobrevive_transito 0,01 / 0,04 | 28 % / 47 % | 84,6 % / 86,8 % | 99,1 % / 95,7 % | 13 / 53 |
| prob_muere_en_nicho 0,00075 / 0,003 | 38 % / 34 % | 76,9 % / 76,5 % | 98,3 % / 98,3 % | 26 / 17 |
| duracion_brote_dias 7 / 28 | 34 % / 34 % | 83,3 % / 75,0 % | 98,3 % | 24 |
| sesgo_homing_estable 0,5 / 2 | 34 % / 34 % | 91,7 % / 91,7 % | 98,3 % | 24 |
| n_vasos 250 / 1000 | 34 % / 34 % | 78,3 % / 95,0 % | 98,3 % | 23 / 20 |
| fraccion_vasos_brote 0,139 / 0,189 (medido) | 34 % / 34 % | 83,3 % / 79,2 % | 98,3 % | 24 |
| *informativo:* despertar irreversible | 34 % | **37,5 %** | 98,3 % | 24 |

Ruido tolerado (fijado antes): D1 ±24 puntos, **D2 ±21,7 puntos** (bootstrap por pacientes).

Veredictos por criterio:
- **D1 PASA** (34 %) y se mantiene en todas las variantes. Sigue siendo calibración.
- **D2 FALLA en la base: 75,0 % (18 de 24 DTC).** La predicción a priori era 84 %; la
  diferencia son unas 2 células de 24 (18 dormidas frente a ~20 esperadas) y cae de lleno dentro del ruido,
  pero la regla dice que se compara con la diana y no se pasa. No se ajusta nada.
- **D3 PASA** (98,3 %), por construcción, como en la v1.
- **D4 NO CONCLUYENTE**: D2 cambia de veredicto con 5 de los 7 supuestos, pero **todos**
  los cambios caen dentro de la banda de ±21,7 puntos. Ningún cambio por encima del ruido.

Lo que dice, sin adornos:
- **La métrica no tiene resolución para responder.** Con 32 pacientes hay 24 DTC en
  total; cada célula vale 4 puntos de D2. El síntoma está a la vista:
  `sesgo_homing_estable` a 0,5 y a 2 dan **el mismo** 91,7 %, cuando el mecanismo
  predice efectos opuestos — es ruido de cuál vaso toca a cada una de 24 células, no
  efecto del parámetro. «No concluyente» aquí no es «casi pasa»: es «con este tamaño
  de muestra D4 no puede detectar nada en D2 que no sea enorme». Es honesto decir que
  el ×0,5/×2 de la v1 (58 % con `prob_entra_dormida` 0,70) sí se habría detectado.
- **El mecanismo sí funciona como se diseñó**: todas las variantes rondan el 1 − *f* ≈
  84 % predicho (media de las 14 variantes ≈ 83 %), y la mediana por paciente — la
  métrica vieja — sigue saturada en 100 % en casi todas, confirmando que no servía.
- **La conclusión fuerte es la del control estructural**: con despertar irreversible D2
  cae al **37,5 %**, muy por debajo de la diana y del ruido. **El resultado depende
  enteramente de que la dormancia sea reversible**, un supuesto con respaldo solo
  cualitativo (Sosa 2014). Eso se había predicho antes de correr. No cuenta para D4
  porque se excluyó de antemano, pero es la dependencia más importante del modelo y
  queda escrita aquí para que nadie la olvide.
- `prob_brote_vascular` (calibrado en la v1) ya no existe; `prob_entra_dormida` tampoco.
  El D2 de la v2 depende de un número **medido** (Wang 2017), no inventado.

**Para avanzar sin ajustar el modelo** (propuesta, no hecha): D2 necesita cientos de
DTC para tener una banda de ruido útil. O más semillas (≈ 128-256), o medir D2 en un
momento posterior o sobre todos los ingresos al nicho; cualquiera de las dos cosas
habría que escribirla como criterio nuevo **antes** de correrla. Mientras tanto la fase
1 no está lista y la fase 2 no se empieza.

---

Lo que ya se sabe que **no** cubrirá: prevención primaria (las tasas de mutación están
aceleradas ~6000×, ver `biologia.CALIBRACIONES`), toxicidad sistémica, y el coste de
apilar intervenciones (sigue siendo gratis en el modelo — es el test T9 pendiente).

---

## Fases siguientes (no se empiezan hasta cerrar la anterior)

**Fase 2 — reactivación.** Reproducir la curva de recurrencias repartidas entre los años
5 y 20 (EBCTCG, *NEJM* 2017, n=62.923). Criterio: que las recaídas no se concentren
todas al principio.

**Fase 3 — terapias sobre las dormidas.** Aquí está la pregunta que importa, y hay una
asimetría que el modelo podría poner a prueba: **la célula T llega a la dormida y el ADC
probablemente no** (los payloads dependen de la división y nunca se han probado contra
células dormidas; la literatura de revisión se contradice y no hay dato primario). Si el
modelo reprodujera esa asimetría sería la primera cosa accionable que dice sobre el
problema que de verdad mata.
