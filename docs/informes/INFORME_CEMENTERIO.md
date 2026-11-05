# El cementerio: qué se ha intentado en cáncer de mama y contra qué murió

> Investigación de literatura y registros de ensayos, 2026-10-21.
> **Idea rectora**: casi ningún fracaso lo fue por un mecanismo falso. Casi todos
> resolvieron *una* barrera y murieron contra otra que no habían medido. El patrón
> repetido no es «el fármaco no funciona», sino **«el fármaco hizo exactamente lo que
> prometía y eso no bastaba»**.

Para qué sirve este documento: **no para copiar soluciones, sino para no repetir errores
y para saber qué barrera sigue abierta.** Y, sobre todo, como **banco de pruebas** — ver
la sección final.

---

## 1. Las once barreras, y quién murió contra cada una

| # | Barrera | ¿Resuelta? | Quién y cómo | Quién murió por ella |
|---|---|---|---|---|
| B1 | Diana con exclusividad suficiente | Parcial: **redefinida**, no resuelta | ADC, bajando el requisito de *exclusividad* a *enriquecimiento* (HER2-ultralow) | Morgan 2010 (muerte por HER2 en epitelio pulmonar); FAP CAR-T (toxicidad ósea letal); **Notch/Wnt/Hh: diana compartida con la célula madre normal** |
| B2 | **Llegar físicamente al tumor** | **NO** | Nadie. Lo mejor del campo: ~1-2 % de dosis inyectada para ADC, ~0,7 % para nanopartículas | Nanopartículas (Wilhelm 2016); y **todo lo validado en ratón a 20-40 %ID/g — un factor ~3.000×** |
| B3 | Penetrar el interior del tumor | Parcial, **sólo los ADC** | Efecto espectador: el payload permeable alcanza **50-200 µm** más allá del anticuerpo (Khera 2021) | Quimio: la **doxorrubicina cae a la mitad a 40-50 µm del vaso** y la hipoxia empieza a 90-140 µm (Primeau 2005) |
| B4 | **Heterogeneidad antigénica** | **SÍ** | **ADC con payload permeable.** Brentuximab: respuestas en **9/14 (64 %) con CD30 escaso o ausente** | CAR-T (escape antigénico); trastuzumab desnudo (inerte en HER2-low) |
| B5 | **Matar sin tocar tejido sano** | **SÍ** | **T-DM1 en EMILIA/TH3RESA: MENOS EA grado ≥3 que el comparador (41 % vs 57 %) CON mejor eficacia.** Nada más en mama mueve ambas direcciones a la vez | SYD985 (queratitis 38 %); vantictumab/ipafricept (**fracturas**); Morgan 2010 |
| B6 | **Matar células que NO se dividen** | **NO** | Nadie. El intento más directo (profármacos de hipoxia) dio **HR de OS 1,06** en sarcoma | Toda la quimio; todos los antimitóticos; **los inhibidores de telomerasa, que REQUIEREN divisiones** |
| B7 | Evitar la resistencia evolutiva | **NO**, sólo retrasada | Terapia adaptativa (Gatenby 2017) y **dosis densa** (CALGB 9741: −26 % recaída con el mismo fármaco y distinto calendario) | **ESR1: 0 % tras IA adyuvante → 36,4 % tras exposición en metastásico.** Evolución pura, observable |
| B8 | **Alcanzar la enfermedad diseminada y oculta** | **NO** | Nadie elimina DTCs. Lo más cercano, bifosfonatos adyuvantes (EBCTCG 2015, n=18.766), que actúan sobre el **nicho**, no sobre la célula | Toda la inmunoterapia (desierto inmune en médula); CAR-T; terapia génica intratumoral |
| B9 | Durar en el tiempo (recurrencia 5-20 años) | Parcial, **sólo hormonoterapia** | **Tamoxifeno (EBCTCG 2011, ~21.000 mujeres): RR de recurrencia 0,53 / 0,68 / 0,97 en los quinquenios.** La ganancia de la primera década no se pierde | Todo lo demás. **Y sólo funciona en ER+** |
| B10 | **Tolerar la COMBINACIÓN** | **NO — e ignorada sistemáticamente** | Sólo una regla de diseño (DeVita/Frei) | magrolimab+aza+ven (EA fatales 19,0 % vs 11,4 %); nivolumab+crizotinib (2 muertes de los primeros 13); osimertinib+durvalumab (ILD 35 %) |
| B11 | Movilizar inmunidad en tumor frío | **NO** en metastásico; parcial en precoz | **KEYNOTE-522: EFS HR 0,65 — y el delta de pCR es IDÉNTICO en PD-L1+ y PD-L1-negativo** | KEYNOTE-119 (ORR 9,6 % vs 10,6 %); todo el checkpoint en ER+; STING; vacunas |

**Lo que la tabla hace visible**: sólo **B4 y B5** tienen un «sí» limpio, las dos las
resolvió **la misma tecnología** (ADC con payload permeable) y por **el mismo mecanismo**
(derramar el citotóxico ~100 µm desde donde se unió). Seis barreras siguen abiertas.

### El matiz que más importa para este proyecto

**Los ADC NO resolvieron «llegar».** Entregan el 1-2 % de la dosis inyectada, el mismo
orden que el 0,7 % de las nanopartículas que llamamos fracaso. Lo que hicieron fue **hacer
irrelevante el mal transporte subiendo la potencia del payload tres órdenes de magnitud**.
Es una solución de **potencia**, no de **transporte**.

Para un simulador que modela el transporte explícitamente, la distinción es crítica: un
ADC no aumenta el flujo al tumor, aumenta la letalidad por molécula entregada.

---

## 2. La barrera que no ha resuelto nadie

> **La conjunción de tres propiedades en la misma célula: está LEJOS del vaso, está
> QUIESCENTE, y está en OTRO ÓRGANO.**

Las tres **covarían**, y esa covarianza es lo que la hace irreducible. La célula a >100 µm
del capilar es la que recibe menos fármaco, es la hipóxica, es la que no se divide, y es
por tanto la que sobrevive a todo lo que selecciona por proliferación. Y la versión
diseminada de esa misma célula vive en un nicho donde ni la quimio, ni el ADC, ni la
célula T llegan en concentración útil.

Lo que lo sostiene:

- **>80 % de las pacientes de mama albergan sólo DTCs Ki67-negativas**, no proliferantes
  (Pantel & Brakenhoff 2004).
- DTCs en **~30-40 % de las pacientes con mama precoz ya al diagnóstico** (Braun, *NEJM*
  2005, n=4.703).
- Todo el arsenal citotóxico presupone proliferación: **la hipótesis log-kill de Skipper se
  derivó de L1210, una leucemia murina de crecimiento exponencial.**
- El intento más directo —profármacos activados por hipoxia— murió en tres fase III, y la
  explicación publicada es «pacientes no seleccionados por hipoxia»: es decir, **nadie sabe
  medir la barrera en el paciente**, que es una forma más profunda de no haberla resuelto.
- Resultado agregado tras tres décadas: **mediana de supervivencia de 25 meses** en 47.034
  pacientes con enfermedad metastásica de novo (SEER).

**La excepción lo confirma por contraste**: el tamoxifeno es lo único con efecto
demostrado a 15 años, y **no mata la célula dormida — suprime la señal que la despierta**.
Es citostático, exige 5 años de exposición continua, y sólo sirve en ER+.

**Lo más nuevo, y por qué no cierra la pregunta**: el ensayo **CLEVER** (*Nature Medicine*
2025, NCT03032406) redujo las células dormidas un 78-87 % con hidroxicloroquina y/o
everolimus, con RFS a 3 años del 91,7-100 %. Pero es **n=51 y sin brazo control**, y los
propios autores lo llaman prueba de concepto. Es la primera vez que alguien demuestra que
se puede reducir farmacológicamente la carga de DTCs en humanos. No es una barrera
resuelta.

---

## 3. Lo que este cementerio dice sobre nuestra estrategia

**Los inhibidores de vías de célula madre tumoral están muertos.** Notch, Wnt y Hedgehog
—es decir, atacar «la reina»— murieron por B1: **la diana es compartida con la célula madre
normal**. Vantictumab e ipafricept produjeron fracturas óseas. Es exactamente el problema
que `identificadores.py` encontró por su cuenta: ningún marcador tiene cobertura y
exclusividad a la vez.

Y la convergencia que importa: nosotros medimos que **matar reinas no cierra nada mientras
queden obreras** (90 días con cero reinas y el tumor se rehízo por plasticidad), y que la
enfermedad residual controlada es **compartimento madre puro**. El mapa dice que la barrera
abierta es la célula quiescente, lejos del vaso y en otro órgano. **Son la misma cosa vista
desde dos lados.**

---

## 4. Apilar intervenciones tiene precio (y el simulador hoy lo regala)

El modelo permite dar troyano + reconducir + anti-SASP + memoria sin pagar nada. Eso es
**una mentira**, y ahora tiene números.

**Serie de Kurzrock** — porcentaje de ensayos en que **todos** los fármacos pudieron darse
a dosis plena:

| combinación | pacientes | a dosis plena |
|---|---|---|
| 2 dirigidos | 8.568 | **51 %** |
| dirigido + citotóxico | 24.326 | **38 %** |
| 3 fármacos | 37.763 | **28 %** |

Heurística de los propios autores: dos dirigidos de clases **no** solapantes admiten ~70 %
de dosis cada uno; si **solapan**, ~30 % cada uno. Si uno se mantiene al 100 %, el segundo
empieza al ~25 %.

**La toxicidad de combinación no sólo es más frecuente: es más rápida.** Metaanálisis de
112 ensayos y 19.217 pacientes (*JAMA Oncol* 2018): mortalidad por toxicidad 0,36 %
(anti-PD-1) → **1,23 %** (PD-1+CTLA-4), y **la mediana de tiempo hasta la muerte tóxica
baja de 40 días a 14,5**.

**Cuidado con atribuir mal la asimetría** «eficacia independiente / toxicidad aditiva»:
**ningún trabajo único la afirma**. Palmer & Sorger 2017 aporta la mitad de eficacia
(y dice literalmente que la toxicidad «no se estudia aquí por datos insuficientes»);
Küçükosmanoglu 2024 (~15 millones de reportes FAERS) aporta la mitad de toxicidad, con
**independencia de Bliss** como nulo empírico.

**Lo que no existe** y conviene no inventarse: no hay fórmula publicada de «cuántos agentes
antes de que la toxicidad sea limitante», ni metaanálisis de EA grado ≥3 por número de
agentes concurrentes. El patrón aproximado por agente añadido es RR 1,1-1,4 con varianza
enorme — heurística de orden de magnitud, no ley.

---

## 5. Batería de validación: el simulador debe reproducir los fracasos

**La regla contra el eco**: un fracaso conocido **no se programa como parámetro**, se usa
como **test**. Si el modelo, calibrado sólo con biología básica, reproduce por su cuenta lo
que pasó en la clínica, entonces sus predicciones sobre lo que nadie ha probado empiezan a
merecer atención. Si se le programa la conclusión, sólo devuelve el eco de quien lo
escribió — eso ya pasó aquí con `contencion_matriz = 0.8`.

| test | debería dar |
|---|---|
| T1 | geometría del O₂ y núcleo necrótico (calibración base) |
| **T2** | **gradiente de penetración del fármaco** — *el que más probablemente falle hoy* |
| T3 | el ADC derrota la heterogeneidad; el anticuerpo desnudo no |
| T4 | degradar el estroma **no** debe ser gratis |
| T5 | equilibrio inmune, y quién lo mantiene |
| T6 | el checkpoint fracasa en tumor frío y funciona sólo donde hay sustrato |
| **T7** | **matar la «reina» no termina la partida** |
| T8 | resistencia evolutiva y el valor del **calendario** |
| **T9** | **apilar intervenciones debe COSTAR** — *corrige la mentira actual* |
| T10 | citorreducción incompleta no debe ayudar |
| T11 | dormancia, recurrencia tardía, y qué la toca |
| T12 | el encogimiento del tumor **no** es supervivencia |

### Lo que el simulador ya acierta sin que se lo enseñáramos

Seis, y es el mejor argumento a favor de seguir:

1. **«Limpiar sin cerrar la fuente no sirve»** — coincide con LION, CARMENA y GOG-0213,
   tres ensayos aleatorizados.
2. **La curva de respuesta y recaída** — el tumor baja a menos de la mitad hacia el día 200
   y rebrota hasta el día 500, porque la oleada terapéutica se agota justo cuando el tumor
   se construye vasos y se arregla la hipoxia.
3. **«Matar la reina no termina la partida»** (T7) — medido aquí antes de leerlo en la
   literatura.
4. **Dar visión al sistema inmune no rescata** — hacer visibles las reinas produjo un 25 %
   más de muertes y cero beneficio, igual que KEYNOTE-119 (ORR 9,6 % vs 10,6 % del control:
   el checkpoint rindió *peor* que la quimio en mama metastásica).
5. **La hormonoterapia sólo funciona como adyuvante** — apagar el campo en solitario da
   10.394 frente a 10.400 sin tratar, es decir nada; pero tras vaciar el tumor sí aporta.
   Es exactamente cómo se usa el tamoxifeno: después de la cirugía, nunca en lugar de ella.
6. **Retirar la hormonoterapia devuelve el riesgo** — de 6.087 a 9.674 peligrosas (5 de 7
   semillas), que es el patrón del RR 0,97 en los años 10-14 del EBCTCG.

Los aciertos 5 y 6 son los más exigentes de los seis, porque ocurren en un régimen
(prevención de la recaída a largo plazo) distinto de aquel en el que se calibró el modelo.

---

## 6. Honestidad sobre las fuentes

El informe original marca lo verificado contra fuente primaria, lo secundario y lo
inferido, y lista once cosas que se buscaron y **no** se encontraron. Las más relevantes:

- **No existe** una cifra única y autoritativa de % de dosis inyectada de ADC que llega a
  tumor humano. El «0,003-0,08 %ID/g» que circula en reviews **no se pudo trazar** a fuente
  primaria, y está en unidades distintas.
- **No existe** radio de muerte espectador publicado para DXd. Las únicas micras medidas
  son de DGN549 (~100-200 µm).
- **El mecanismo del exceso de infecciones fatales con magrolimab sigue sin explicar**: las
  tasas de infección y de neutropenia fueron **similares entre brazos**. No presentarlo
  como «más mielosupresión».
- Resultados de CART-TnMUC1: **sin publicar**, que no es lo mismo que «cero respuestas».
