# Bitácora de investigación

Registro cronológico, entrada a entrada, de todo lo que se ha probado en el simulador
desde el 19 de octubre de 2026: hallazgos, **errores propios corregidos**, expectativas
que fallaron y decisiones de diseño. Se publica tal cual se escribió, con los tropiezos,
porque en investigación el camino importa tanto como el resultado.

> Notas de lectura
> - **Quién escribe:** las entradas son de **Claude** (Anthropic), la IA que desarrolla el
>   proyecto. Cuando dicen «mi error» o «mi expectativa», hablan de Claude. «El autor» es
>   **Ismael Ben**, que dirige el proyecto y propone las ideas.
> - Muchas entradas se escribieron sin tildes (era un cuaderno de trabajo).
> - «Peligrosas» = células de carcinoma in situ + invasivo + células madre tumorales.
> - Algunas entradas citan scripts exploratorios de un solo uso que no forman parte del
>   repositorio; los resultados que importan están en los scripts principales y en
>   `validacion.py`.
> - Las entradas más recientes están al final. El resumen de cada avance está en
>   [`docs/avances/`](avances/).

### 2026-10-19 · inicio
Proyecto creado el 2026-10-19. Conducto mamario 2D: genes driver reales (PIK3CA, ERBB2, MYC, TP53, PTEN, CDH1, BRCA1/2, HLA, CD274), O2 con necrosis, inmunes T/NK, terapias "IA" (cazador, car_her2, reparador, hibrido, centinela, quimio, mi_celula). main.py = visor con Play; experimento.py = comparativa sin ventana.

### 2026-10-19
 añadida jerarquía reina/obreras (células madre tumorales + plasticidad); la estrategia del autor "arreglar la reina y extirpar obreras" (terapia reina) es la mejor en experimento.py; el autor rechaza la quimio por dañina.

### 2026-10-19 (2)
 guardianas (reina curada que fagocita obreras) con puerta Y = 0 sanas comidas; modos de vasos incl. "reconducir" O2 hacia el tumor (el autor NO quiere cortar oxígeno, quiere reconducirlo).

### 2026-10-19 (3)
 los restos tras la guardiana son lesiones NUEVAS; guardiana+memoria (15 vigilantes residentes) = 0 tumorales a 500 y 700 días. El autor no quiere cirugía ni daño al sujeto.

### 2026-10-19 (4)
 el autor quiere recorrer el camino hasta la realidad; hecho A1-A3 (datos_reales.py con TCGA vía cBioPortal, calibración --neoantigenos 0.1); HOJA_DE_RUTA.md con fases A4-A7 pendientes (PhysiCell 3D, gemelo digital por paciente). Se le explicó que laboratorio y ensayos no se reducen a meses.

### 2026-10-19 (5)
 A7 sensibilidad: la estrategia aguanta subtipo, BRCA1 y tratamiento tardío, pero se hunde con plasticidad x10 (312 peligrosas vs 2,7). La pregunta decisiva es biológica: facilidad de desdiferenciación in vivo.

### 2026-10-19 (6)
 A8 reproducción: añadido TERT/telomerasa y contador de divisiones (Hayflick). TERT emerge al 100%. "Castración" (antitelomerasa) sola casi inutil (8251->7199, como los ensayos reales fallidos) pero como refuerzo de guardiana+memoria mejora 6x (264->42). Propagador oncolitico amplifica pero es a cara o cruz.

### 2026-10-19 (7)
 A5 iniciado sin PhysiCell: tejido3d.py (cilindro, 26 vecinos, vías vena/intratumoral/aguja, --matriz, --tiodg 6-tio-dG). El 3D destapa el problema de ACCESO (interior de la esfera inalcanzable). identificadores.py -> INFORME_IDENTIFICADORES.md: ningún marcador tiene cobertura+exclusividad; salida = puertas AND (telomerasa Y aneuploidía, o fuera del nicho). Siguiente natural: meter esa puerta AND en el simulador y atacar el acceso en 3D.


### 2026-10-20
 A5b. Puertas AND implementadas (nicho, tert, tert+nicho) + puertas.py (sensib./especif. por puerta): la aneuploidia llega TARDE (0% a dia 120, 79% a dia 300), asi que no sirve en fase precoz. Sin barrera, guardiana_nicho+memoria fue la mejor (487 peligrosas vs 3104 de guardiana).

### 2026-10-20 (2)
 anadida PRESION INTERSTICIAL (campo sim.P, solo frena a quien ENTRA; baja con degradar_matriz y normalizar vasos) + barrera.py. Anadidos MACROFAGOS, CD47 y calreticulina ("comeme"): la fagocitosis es una puerta Y natural. bloquear_cd47(dirigido=True) = biespecifico (evita anemia).

### 2026-10-20 (3)
 HALLAZGO. La secuencia del autor (1 anti-CD47, 2 castrar, 3 que se lo coman los macrofagos) falla SOLO por el acceso: con presion 3 fagocitosis, sin presion 920 (x300), con 620 celulas ya exponiendo "comeme". La puerta Y protege de la anemia (27 sanas sin freno -> 0 comidas). PERO comerse 920 celulas no curo (957 -> 932 peligrosas): limpiar no es cerrar la fuente. El autor NO quiere fase precoz: quiere curar en cualquier fase.


### 2026-10-20 (4)
 cadena.py = secuencia completa del autor (abrir matriz+vasos / anti-CD47 dirigido / 6-tio-dG / macrofagos / cerrar fuente con guardiana_nicho+memoria). Resultados 3 semillas, 500 dias: nada 7286, solo abrir 4302, solo guardiana+memoria 7230 (!), cadena SIN abrir 9258, cadena SIN cerrar 10828 (el peor), CADENA COMPLETA 5658 pero semillas [2, 10922, 6050]. CONCLUSIONES: (a) la presion intersticial TUMBA la estrategia de la guardiana que funcionaba antes (487 -> 7230): todo el trabajo previo asumia un tumor sin barrera; (b) los dos controles validan que hacen falta abrir Y cerrar la fuente; (c) la varianza entre semillas domina al tratamiento: usar medianas y 5-7 semillas, nunca medias de 3.


### 2026-10-20 (5)
 TECHO de la cadena (7 semillas, correlacion carga/desenlace): por encima de ~4000 celulas el dia del tratamiento NO funciona nunca (3/3 fracaso); por debajo de ~1000 funciona 3 de 4. Lo que decide es el tamano al empezar, no el tratamiento. Es decir: es una herramienta de fase precoz, justo lo que el autor NO queria. Todo escrito en INFORME_BARRERA.md.


### 2026-10-20 (6)
 LA MALLA. Implementados CAFs (myCAF que se activan por TGF-b del tumor y por rigidez = bucle YAP/TAZ), campo de matriz sim.M/Ms, LOX, y tratar_matriz() con 5 modos separados (hialuronidasa / anti_lox / losartan / anti_tgfb / matar_caf). CLAVE: la malla tambien CONTIENE (contencion_matriz) -> matar_caf sale PEOR que no hacer nada, reproduciendo Ozdemir 2014 y Rhim 2014. Ninguna intervencion sobre la malla abre el acceso (0 fagocitosis en las tres probadas). En 3D, degradar matriz 50% no movio nada: 26% accesible y 1% de nucleo alcanzado antes y despues (ese 26% es GEOMETRIA pura, no presion).

### 2026-10-20 (7)
 CABALLO DE TROYA (idea del autor). Agente con troyano=0.9 (la presion no le frena: entra invitado como los monocitos reclutados por CCL2) + entrega= (suelta la terapia real al verse rodeado) + oxigeno= (bombona) + castra= (esteriliza en vez de matar). Resultados 3 semillas: el transporte FUNCIONA (32-473 entregas vs 3 fagocitosis nadando). La dosis importa y satura: 40 -> [940], 120 -> [9], 300 -> [5]. La bombona de O2 NO aporta nada (7058 sin O2 vs 7053 con O2). CASTRAR ENRIQUECE EN REINAS (1896 -> 2301): al esterilizar a las que se dividen seleccionas a las quiescentes = reinventamos el enriquecimiento en CSC de la quimio. Anadir el componente de memoria lo frena (2282 -> 1979).

### 2026-10-20 (8)
 INVESTIGACION de fracasos reales. Magrolimab (anti-CD47, Gilead, 4.900M$): ENHANCE-3 parado por AUMENTO de mortalidad (infecciones, fallo respiratorio), MDS y AML por futilidad. PEGPH20 fase III fallo y empeoro supervivencia. Nanoparticulas: 0,7% de la dosis llega (Wilhelm 2016, sin mejorar en 10 anos). 82% de los farmacos oncologicos que fracasan en fase tardia es por falta de EFICACIA; 62% de los fase III no alcanzan su objetivo. CONCLUSION DE LA SESION: cada fracaso resolvio UNA barrera de varias; el cancer es un proceso evolutivo (toda presion unica selecciona resistentes); y las terapias se validan en el regimen donde funcionan -- el salto raton->humano es nuestro salto semilla 0 -> semilla 1.


### 2026-10-20 (9)
 CORRECCION IMPORTANTE DEL RELATO DEL PROYECTO. La "guardiana" (reina curada y reprogramada) NO hace el trabajo: en 450 dias se crean 23, queda 1 viva y fagocita 5 celulas de 10.029 peligrosas. Lo que funciona son 82.102 muertes por ataque directo + 3.520 reparaciones genicas. Por eso el barrido de reversion salio plano (6651/6601/6652/6598/6635 para estado permanente -> vida media 5 dias): el estado reprogramado es IRRELEVANTE, no es que la herencia lo proteja. Consecuencia BUENA: la estrategia no depende de la pieza inventada, depende de dos cosas que si existen -- celulas terapeuticas persistentes que matan (CAR-T) y edicion genica in vivo (CRISPR). README.md y HOJA_DE_RUTA.md describen un proyecto que ya no es el que hay: PENDIENTE reescribirlos. Tambien pendiente: barrido de reversion sin herencia (5 filas) -- lo mato el gestor por falta de memoria del sistema, no relanzar sin permiso.

### 2026-10-20 (10)
 mejor configuracion actual = troyano+guardiana_nicho+memoria a dosis 300 (6651 peligrosas, 826 reinas, semillas [3, 9922, 10029]). El vehiculo troyano SI aporta: reduce reinas de 1540 a 826 y es lo unico que ha movido las semillas grandes en todo el dia (~4%, posible ruido). El castrador hay que QUITARLO: enriquece en reinas (826 -> ~2300).


### 2026-10-21
 E0 APLICADO (instrumentacion). tejido.py: 5 contadores nuevos; las tres ramas de _reparar se cuentan por separado (reparaciones_que_curan / _p53_letales / _esteriles) y el farmaco p53 (tejido.py:665) y el cebo (tejido.py:714) salen del contador compartido a reparaciones_farmaco/_cebo; visor.py muestra el desglose. HALLAZGO: 91,2% de las reparaciones son ESTERILES (gen corregido, celula sigue tumoral), 8,3% curan, 0,5% matan via p53, mas 1738 cargas gastadas en genes ya sanos (8x las reparaciones efectivas).

### 2026-10-21 (2)
 BUG DE DISENO encontrado por el agente investigador. TERT cuenta como driver (n_drivers, tejido.py:85) pero NO esta en TODOS_LOS_GENES (politicas.py:154). Lo tiene el 91,5% de las reinas -> el 91,5% es INCURABLE POR CONSTRUCCION, nunca llega a n_drivers==0 ni puede volverse guardiana. Las 23 guardianas de 450 dias eran exactamente el ~9% de reinas sin TERT. Fuga silenciosa adicional en politicas.py:107-109: si quedan solo genes no reparables, la reina se marca {} = "sana, no tocar" y ni se repara ni se ataca.

### 2026-10-21 (3)
 EXP-2 (comer vs castrar) con bajas EMPAREJADAS EXACTAS (tope de bajas, 5 semillas, 401 vs 400 bajas): matar 5485 peligrosas / 408 reinas; castrar 6625 / 838. La HIPOTESIS DEL HUECO QUEDA FALSADA: matar no sale peor por dejar la celda libre, el resultado en peligrosas es ruido (matar gana 2, castrar gana 2, empate 1). La interpretacion previa del proyecto SE MANTIENE: el dano de limpiar es no cerrar la fuente, no liberar espacio. Robusto si: castrar duplica las reinas (enriquecimiento en CSC, 3a vez que aparece). METODO: no se pueden emparejar bajas bajando la dosis (castrar es un sumidero por tejido.py:1100 'if c.senescente: return'); hay que emparejar por tope de bajas.

### 2026-10-21 (4)
 INVESTIGACION real del brazo fagocitosis. Magrolimab muerto en las 3 indicaciones (ENHANCE-3: 19,0% eventos fatales vs 11,4%, por INFECCIONES no por anemia); el fase 2 en triple negativo (ELEVATE-TNBC) terminado sin publicar; 0 ensayos activos 2026. En MAMA no existe ni una evidencia aleatorizada del bloqueo CD47. Evorpacept (Fc inactivo) es lo mas vivo pero fallo ASPEN-03/04 y renuncio a fase 3 en EEUU. CAR-M CT-0508 (Nature Medicine 2025, 14 pacientes, 8 de mama): CERO respuestas objetivas, Carisma en liquidacion. Agonistas STING: fracaso en bloque, y STING MATA las propias celulas T (Larkin 2017). "Limpiar sin cerrar la fuente" tiene 3 ensayos aleatorizados en contra: LION (NEJM 2019, quitar ganglios con metastasis ocultas en 56% no cambio nada), CARMENA (NEJM 2018, cirugia salio PEOR: 13,9 vs 18,4 meses), GOG-0213. La excepcion DESKTOP III solo fue positiva donde se lograba reseccion COMPLETA = cerrar la fuente.

### 2026-10-21 (5)
 el puente innata->adaptativa NO EXISTE en el modelo (fagocitar solo borra la celda y suma contador, tejido.py:1109-1118: sin presentacion de antigeno, sin dendriticas, sin reclutamiento secundario). En la realidad esta demostrado en raton, disputado (macrofago vs cDC1) y NUNCA traducido a clinica. Lo unico que mantiene un cancer oculto indefinidamente (equilibrio, Koebel 2007) es inmunidad ADAPTATIVA.

### 2026-10-21 (6)
 calibradores humanos para no venderse humo. Rezatapopt (reactivador p53 Y220C, PYNNACLE fase 2, corte ago-2025): ORR 33%, duracion mediana 6,2 meses, MAMA 18%. Restaurar un driver da respuestas y meses, no curaciones. No hay NINGUN ensayo de edicion in vivo en tumores solidos (lo clinico es hepatico); lo mejor en solido es preclinico (LNP intratumoral anti-SOX2, ~17% edicion -> 50% regresiones, Masarwy 2025) y entra POR AGUJA, no por vena. La correccion genica NO tiene efecto espectador (el "10% basta" de Freeman 1993 es de MATAR con HSV-TK via uniones comunicantes), asi que hay que editarlas casi todas salvo que se encienda propagacion_correccion (biologia.py:136-145, _corregir_vecinas) -- es el unico mecanismo del repertorio cuyo coste no crece con el numero de celulas.

### 2026-10-21 (7)
 tambien es dudosa la rama optimista n_drivers==0 -> celula sana (tejido.py:1144): la celula "curada" conserva la aneuploidia (tejido.py:398), presente en ~95% de tumores. La rama BIEN avalada es la de p53 restaurado -> 50% apoptosis (Ventura y Xue, Nature 2007).


### 2026-10-21 (8)
 E1 (TERT en los vectores, parcheado en memoria, NO consolidado en el default). Confirma el bug: guardianas y reparaciones_que_curan x2,65 (20 -> 53), esteriles bajan (184 -> 138), reparaciones totales planas (212 -> 196, el cuello no eran las reparaciones). PERO el desenlace NO mejora: peligrosas por semilla [326,1,145,374,1111] sin TERT vs [273,1,261,451,2100] con TERT = empeora en 3 de 5 y casi DUPLICA en la semilla grande. Decision: no consolidar (ademas en mama TERT se activa por amplificacion, no por la mutacion de promotor editable).

### 2026-10-21 (9)
 BARRIDO DEL EFECTO ESPECTADOR (propagacion_correccion x alcance_correccion x error_edicion, 5 brazos x 5 semillas). El mecanismo propaga de verdad (266 correcciones, 125 celulas corregidas por completo) pero (a) NO es exponencial: x10 en probabilidad y x2 en alcance solo lleva de 110 a 266 propagadas = SUBLINEAL, la cadena se ahoga porque corrige un alelo por intento y cada corregida queda rodeada de tumorales que se dividen mas rapido; (b) en el tumor grande EMPEORA SIEMPRE: semilla 5 control 1111 vs [2839, 2471, 3070, 3057] = 4 de 4 brazos, entre 2,2x y 2,8x. El error_edicion al 5% no fue lo que lo estropeo (solo 15 errores).

### 2026-10-21 (10)
 CONCLUSION DEL DIA, y es un cierre. TRES experimentos independientes apuntan al mismo mecanismo: (1) barrido de reversion plano (ayer), (2) E1 TERT, (3) efecto espectador. Todos comparten la causa, que esta en el codigo: convertir una tumoral en "curada" la SACA DE LA LISTA DE DIANAS pero conserva la aneuploidia (tejido.py:398) y puede re-mutar sin proteccion -> se fabrican celulas que nadie vigila y vuelven a romperse, y cuanto mas grande el tumor mas caro sale. EL BRAZO DE "CURARLO REPROGRAMANDO" ESTA MUERTO EN EL SIMULADOR, y no por falta de potencia: se le dio toda, incluido un mecanismo que en la realidad NO EXISTE (la correccion genica no tiene efecto espectador), y sale peor. Lo que funciona es MATAR. Encaja con rezatapopt: restaurar p53 da 33% de respuestas y 6,2 meses, y la resistencia vuelve por drivers alternativos generados por inestabilidad cromosomica = la misma aneuploidia que nuestras celulas "curadas" conservan.



### 2026-10-21 (11)
 MATAR Y CERRAR. La fuente de reinas tiene TRES grifos (tejido.py:262-278): division asimetrica (se cierra matando), plasticidad x3 CON HIPOXIA (se cierra oxigenando = "reconducir", justo lo que el autor siempre quiso), y SASP de las senescentes (farmaco_anti_sasp, JAK/ruxolitinib). A 250 dias: reconducir baja desdiferenciaciones 141->83 y reinas 90->38 (4 de 5 semillas). "Solo cerrar sin matar" NO controla (8121->5673 pero por semilla es ruido): hacen falta las dos cosas.

### 2026-10-21 (12)
 EL SASP ES CONDICIONAL, y asi se explica por fin el enriquecimiento en reinas. Con el brazo MATADOR (18 senescentes) el anti-SASP no hace NADA: resultado identico digito a digito, el termino no se activa ni una vez. Con el brazo CASTRADOR (132 senescentes) SI: desdiferenciaciones 900->431 y reinas 875->407, bajando en 5 DE 5 SEMILLAS. O sea: las celulas que jubilas GRITAN y el grito fabrica reinas nuevas (senescencia inducida por tratamiento + SASP promoviendo stemness, literatura establecida, antidoto real = inhibidores JAK). PERO peligrosas no se mueve (6393->6370): entiendes el mecanismo, lo neutralizas, y no sirve solo. Matar sin castrar sigue siendo 6x mejor (981).

### 2026-10-21 (13)
 BIFURCACION a 500 dias, 7 semillas (el barrido mas caro del dia, ~22 min). solo matar pel=[9913,1,8692,5608,9906,9718,9996] vs matar+reconducir pel=[12,1,8268,10392,10428,10270,9]. NO es ruido y NO es victoria: 2 curaciones casi totales (9913->12, 9996->9), 2 empates, 3 EMPEORAMIENTOS. Mediana enganosa (9718->8268). Mecanismo propuesto: reconducir O2 hace dos cosas OPUESTAS -- quita la hipoxia que triplica la fabricacion de reinas (bueno) pero riega mejor el tumor y lo deja crecer (malo); segun la semilla gana una u otra, y cuando gana la buena gana por goleada. PENDIENTE: comprobar si las semillas que curan son las que llegaban pequenas al dia 120 (encaja con el techo de ~4000 celulas de 2026-10-20 (5)).



### 2026-10-21 (14)
 la bifurcacion NO la explica la carga al dia 120 (resultado negativo, pero matiza el techo de ayer). Cruce de las 7 semillas (carga/reinas/hipoxicas al dia 120 vs desenlace a 500): no es monotonico en NINGUNA variable -- s1 CURA con 665 celulas y s3 NO cura con 395; cura con 33 reinas y no cura con 20. Todas las semillas llegan pequenas al dia 120 (<1700) y aun asi se parten en dos. Conclusion: el tamano inicial explica el FRACASO de los tumores grandes (techo ~4000) pero NO explica el EXITO entre los pequenos; lo que decide la bifurcacion ocurre DESPUES de intervenir, no antes.



### 2026-10-21 (15)
 el SASP NO se abre solo con el tiempo. A 500 dias el brazo matador da resultado IDENTICO digito a digito con y sin anti-SASP (5 semillas), y las senescentes DISMINUYEN con el tiempo (22 d120 -> 20 -> 13 -> 10 -> 6 final): el tumor no acumula senescencia, la diluye creciendo. CERRADO: el SASP es exclusivamente un efecto secundario de CASTRAR. Leccion practica = "no castres", no "castra y anade ruxolitinib".

### 2026-10-21 (16)
 HALLAZGO INESPERADO de la serie temporal (lo mejor del dia). Hay RESPUESTA Y RECAIDA, emergente, sin programarlo: peligrosas 673 (d120, tratamiento) -> 298 (d200, MINIMO) -> 335 (d300) -> 1290 (d400) -> 8692 (d500). Es la curva clinica clasica respuesta/meseta/progresion. La causa esta en la columna de al lado: las desdiferenciaciones NO PARAN NUNCA (28 -> 109 -> 185 -> 367 -> 837). REENCUADRE: el problema NO es matar poco (de hecho el tumor baja a menos de la mitad), es DURABILIDAD -- el dia 200 vamos ganando y el dia 500 hemos perdido, porque el grifo de la fuente sigue abierto mientras matamos. Coincide con el problema real del campo (recurrencias en mama repartidas entre los anos 5 y 20, EBCTCG 2017). CONSECUENCIA PARA LA BIFURCACION: mirar si en las semillas que colapsan las REINAS llegan a CERO en ese minimo del dia 200 = punto de no retorno.

### 2026-10-21 (17)
 METODO CONTRA EL ECO. Al implementar lecciones de fracasos ajenos NO meterlas como parametros (eso produjo el falso "hallazgo" de contencion_matriz=0.8 que solo repetia a Ozdemir/Rhim). Los fracasos entran como BANCO DE PRUEBAS: el simulador debe reproducirlos SIN que se los contemos. Ya hay 2 aciertos no programados: "limpiar sin cerrar la fuente no sirve" (coincide con LION/CARMENA/GOG-0213) y la curva respuesta-recaida. PENDIENTE: el simulador apila intervenciones GRATIS y eso es una mentira -- las toxicidades se suman en la realidad (magrolimab solo era tolerable; lo que mato fue combinarlo con azacitidina+venetoclax). Hace falta un coste por intervencion o el modelo dira siempre "dalo todo a la vez".



### 2026-10-21 (18)
 LA BIFURCACION ES CAOS DETERMINISTA. El 9913->12 NO es un efecto terapeutico. Prueba: misma semilla, mismo dia, cambiando la dosis UNA SOLA CELULA -> s1: 299=10102, 300=12, 301=9725; s7: 299=7, 300=9, 301=659. Y variando solo el dia de intervencion el estado del tumor empeora MONOTONAMENTE (mas celulas, mas reinas, menos O2) mientras el desenlace salta sin orden: d118=9811, d120=12, d122=8122 (en d118 el tumor estaba MEJOR y acabo mil veces peor). Eso demuestra que NINGUNA funcion monotona de esas variables puede predecir el desenlace: no es que falte el biomarcador, es que no puede existir. Validacion fuera de muestra (semillas 8-14): CERO rescates; en 6 de 7 el tumor nunca fue hipoxico y reconducir no hizo literalmente nada (trayectorias identicas). REGLA: reportar siempre distribuciones sobre muchas semillas y NUNCA presentar un 9913->12 como efecto terapeutico.

### 2026-10-21 (19)
 que es "reconducir" en realidad, cuantificado: x6,7 velocidad de division (tejido.py:321, pdiv *= 0.15 si hipoxica) a cambio de /3 fabricacion de reinas (tejido.py:273) y de quitarle al tumor el motivo para angiogenesis (tejido.py:314). Dos fuerzas grandes que se deciden cerca de la frontera, y cerca de la frontera el sistema es caotico.

### 2026-10-21 (20)
 DOS PIEZAS SOLIDAS QUE INVIERTEN LA ESTRATEGIA. (a) La enfermedad residual controlada es COMPARTIMENTO MADRE PURO: en los 9 finales controlados peligrosas == reinas exactamente (12=12, 9=9, 19=19, 4=4...). (b) Llegar a CERO REINAS NO ES ABSORBENTE: una corrida estuvo 90 DIAS SEGUIDOS con 0 reinas y 0 peligrosas, con ~90 obreras vivas, y se rehizo (las desdiferenciaciones siguieron subiendo durante esa ventana: cada reina nueva nacia de una obrera y la mataban, en equilibrio, hasta que una escapo). CONSECUENCIA: durante tres dias la idea fue "mata la reina, que es la fuente"; pero mientras quede UNA OBRERA la plasticidad rehace la reina. El punto de no retorno, si existe, es VACIAR EL DEPOSITO DE OBRERAS, no erradicar reinas.

### 2026-10-21 (21)
 por que hay respuesta y recaida, mecanismo verificado: la oleada terapeutica es TRANSITORIA (troyano vida 400, guardiana_nicho vida 500; a partir del dia ~340 solo quedan los 300 agentes permanentes de memoria, que tiene expansion=0) y justo entonces el tumor SE ARREGLA SU PROPIA HIPOXIA construyendo vasos: O2 medio 0,236 -> 0,275 -> 0,430 -> 0,49 y las peligrosas detras 324 -> 335 -> 672 -> 1606. El O2 sube ANTES del despegue. La tasa de muerte se queda clavada en 20-45 muertes/dia aunque el tumor pase de 427 a 2761 celulas. PENDIENTE sin explorar: las reinas tienen neoantigeno efectivo 0 SIEMPRE (_observar hace neo //= 3 para madres), asi que el TCR nunca las ve -- las mata la puerta de nicho, no el antigeno.



### 2026-10-21 (22)
 MAPA DEL CEMENTERIO -> INFORME_CEMENTERIO.md (nuevo). De ONCE barreras en 3 decadas, solo DOS resueltas (B4 heterogeneidad antigenica y B5 indice terapeutico), las dos por la MISMA tecnologia (ADC con payload permeable) y el MISMO mecanismo (derramar el citotoxico ~100 um desde donde se unio): brentuximab 9/14 respuestas con diana escasa/ausente; T-DM1 en EMILIA MENOS EA grado >=3 que el comparador (41% vs 57%) CON mejor eficacia. MATIZ CRITICO: los ADC NO resolvieron "llegar" -- entregan 1-2% de dosis inyectada, el mismo orden que el 0,7% de las nanoparticulas; subieron la potencia del payload 3 ordenes de magnitud. Es solucion de POTENCIA, no de TRANSPORTE (distincion clave para un simulador que modela transporte).

### 2026-10-21 (23)
 LA BARRERA QUE NO HA RESUELTO NADIE = la CONJUNCION de tres propiedades en la misma celula: LEJOS del vaso + QUIESCENTE + en OTRO ORGANO. Covarian, y por eso es irreducible (doxorrubicina cae a la mitad a 40-50 um del vaso, la hipoxia empieza a 90-140 um, y esa misma celula es la que no se divide). >80% de las pacientes de mama albergan SOLO DTCs Ki67-negativas; el log-kill de Skipper se derivo de una leucemia murina exponencial. Los profarmacos de hipoxia (el intento mas directo) murieron en 3 fase III por "pacientes no seleccionados por hipoxia" = nadie sabe MEDIR la barrera. Excepcion que lo confirma: el tamoxifeno es lo unico con efecto a 15 anos y NO mata la celula dormida, suprime la senal que la despierta (y solo en ER+). Lo mas nuevo: CLEVER (Nature Medicine 2025, n=51 SIN control) redujo DTCs 78-87% -- primera demostracion en humanos, no barrera resuelta.

### 2026-10-21 (24)
 ATACAR "LA REINA" YA SE INTENTO Y ESTA MUERTO. Notch/Wnt/Hedgehog murieron por B1: la diana es COMPARTIDA CON LA CELULA MADRE NORMAL (vantictumab e ipafricept dieron FRACTURAS OSEAS). Es el mismo muro que encontro identificadores.py por su cuenta. CONVERGENCIA: nosotros medimos que matar reinas no cierra nada mientras queden obreras y que la residual es compartimento madre puro; el mapa dice que la barrera abierta es la celula quiescente lejos del vaso y en otro organo. Son la misma cosa vista desde dos lados.

### 2026-10-21 (25)
 PRECIO DE APILAR (el simulador hoy lo regala = mentira). Serie Kurzrock, % de ensayos en que TODOS los farmacos pudieron darse a dosis plena: 2 dirigidos 51% (n=8.568), dirigido+citotoxico 38% (n=24.326), 3 farmacos 28% (n=37.763). Heuristica: 2 dirigidos de clases NO solapantes ~70% de dosis cada uno; si SOLAPAN ~30% cada uno. Y la toxicidad de combinacion es MAS RAPIDA: mediana hasta muerte toxica 40 dias (monoterapia) -> 14,5 dias (combinacion), con mortalidad 0,36% -> 1,23% (JAMA Oncol 2018, 112 ensayos, 19.217 pacientes). OJO: la asimetria "eficacia independiente / toxicidad aditiva" NO la afirma ningun trabajo unico -- hay que ensamblarla con Palmer&Sorger 2017 (eficacia; dice literalmente que la toxicidad no la estudia) + Kucukosmanoglu 2024 (toxicidad, Bliss como nulo). NO existe formula publicada de "cuantos agentes antes de que la toxicidad limite".

### 2026-10-21 (26)
 PRIORIDAD SIGUIENTE = A15, la BATERIA DE VALIDACION (12 tests falsables en INFORME_CEMENTERIO.md §5). Razon: mientras el simulador no reproduzca fracasos conocidos que NO se le han programado, sus predicciones no valen fuera de la carpeta. Ya acierta TRES sin ayuda: (1) "limpiar sin cerrar la fuente no sirve" = LION/CARMENA/GOG-0213; (2) la curva respuesta-recaida; (3) "matar la reina no termina la partida" (su T7), medido aqui ANTES de leerlo en la literatura. El test que mas probablemente falle hoy es T2 (gradiente de penetracion del farmaco); el que corrige la mentira actual es T9 (apilar debe costar).



### 2026-10-21 (27)
 A17 RESPONDIDO: el limite NO es VER, y tampoco es MATAR. Parcheando tejido.py:972 (neo //= 3 para madres) a divisor 1 = reina tan visible como obrera: peligrosas [9913,1,8692,5608,9906] -> [9918,1,8492,10136,10079] = CERO mejoras en 5 semillas (4 planas, 1 empeora mucho: s4 5608->10136 y sus reinas 424->1433; la mediana la domina esa sola semilla, asi que NO se puede afirmar "empeora"). Lo robusto: al hacerlas visibles se MATA MAS (inmunes 1676->2895, terapia 9473->11879, ~25% mas) y hay igual o mas tumor. DATO ADICIONAL: las reinas NO son uniformemente invisibles -- visibles/vivas por semilla era (3,1138),(0,1),(416,611),(0,424),(67,1717): en una semilla 416 de 611 ya eran visibles y el tumor gano igual. EXPLICACION que encaja con todo lo del dia: es problema de ATENCION, no de deteccion -- al ver las reinas los agentes gastan sus pasos persiguiendolas (estan en el nicho, son caras de matar) en vez de atacar OBRERAS, y el deposito de obreras es lo que regenera reinas. Dar vision al sistema inmune lo concentra en la diana equivocada. Coincide con KEYNOTE-119 (ORR 9,6% vs 10,6% del control: el checkpoint en mama metastasica rindio PEOR que la quimio). CONCLUSION: el limite esta en A QUE SE DEDICA lo que mata, y en que la oleada se agota antes que el deposito.



### 2026-10-21 (28)
 A15 EMPEZADO -> validacion.py (nuevo). Bateria de 5 tests implementados (T1, T2, T7, T10, T12) + 7 listados como pendientes con lo que necesita cada uno. Cada test lleva su fuente al lado y hay TRES veredictos (PASA / FALLA / no concluyente). Regla escrita en la cabecera: un fracaso clinico NO se programa como parametro, se usa como TEST; y NO se toca el modelo para que un test pase. Escala: 1 celda = 15 um (estaba en un comentario de biologia.py:59, no como constante). Nota tecnica: numpy devuelve np.bool_ y rompe los `is True`; y hay que forzar sys.stdout a UTF-8 porque Windows usa cp1252 y revienta al redirigir.

### 2026-10-21 (29)
 RESULTADO DE LA BATERIA: 4 PASAN, 1 no concluyente. T1 (geometria del O2) PASA y bien: hipoxia a 165 um frente a los 150-180 um de Thomlinson & Gray 1955, necrosis a 240 um -- la escala espacial del simulador esta calibrada sin que nadie la calibrara, y era el test base (si fallaba, todas las conclusiones sobre acceso y barrera fisica serian humo). T7 (matar la reina no termina la partida) PASA: reinas minimo 31 -> final 1125. T10 (citorreduccion incompleta no ayuda) PASA: 0/3 semillas mejoran. T12 (encogimiento != supervivencia) PASA: minimo 112 -> final 6957, recaida >10x en 2/3.

### 2026-10-21 (30)
 LECCION DE METODO con T2. Primero salio FALLA con letalidad 12% cerca del vaso y 45% lejos (al reves de lo esperado). Pero el test estaba MAL: contaba como "muerta por farmaco" cualquier celula no viva al final, y a 180-210 um esta el nucleo necrotico (lo confirma el propio T1, necrosis a 240 um). Anadido brazo CONTROL sin farmaco y restada la mortalidad basal -> la letalidad ATRIBUIBLE es ~0% a TODAS las distancias (0/9/0/-0/1%). El hallazgo real no es que el farmaco llegue uniformemente: es que el agente "quimio" es CASI INERTE. No se puede medir el gradiente de algo que no mata, asi que T2 queda en NO CONCLUYENTE hasta darle efecto a quimio. Moraleja: el primer instinto fue "el modelo falla"; el fallo estaba en el test.



### 2026-10-21 (31)
 T13 ANADIDO a validacion.py = el REQUISITO DEL AUTOR como test no negociable ("eliminar lo tumoral sin danar al sujeto"). PASA limpio: 0 dano a sanas sobre 14.906 muertes totales, 3/3 semillas (colateral, sanas comidas, madres sanas tocadas, castradas sanas, extirpadas sanas = todo 0). T2 reformulado para medir con NUESTRA terapia (guardiana_nicho) en vez de quimio -- el autor rechaza la quimio y ademas una celula terapeutica penetra PEOR que una molecula pequena, asi que el test queda mas duro. Bateria actual: 6 tests, 5 pasan, T2 pasa con reserva (el perfil no es un decaimiento monotono: el maximo esta a 90-120 um, no pegado al vaso).

### 2026-10-21 (32)
 HALLAZGO QUE CIERRA EL DIA -- LA TERAPIA SELECCIONA EL REFUGIO HIPOXICO. Mi hipotesis era que matar la periferia oxigena el nucleo (por un -21% raro en T2). REFUTADA: los absolutos enganan (necrosis 660->59, hipoxicas 865->212) porque el tumor es 20x mas pequeno (8208->424 tumorales). NORMALIZADO es al reves: fraccion hipoxica 0,094 -> 0,486 (x5), O2 medio del tumor 0,402 -> 0,241 (-40%), necrosis por celula tumoral 0,078 -> 0,139. La terapia mata donde LLEGA (periferia oxigenada) y lo que sobrevive es el nucleo hipoxico: no rescata al tumor, lo DESTILA.

### 2026-10-21 (33)
 EL BUCLE QUE EXPLICA LA RECAIDA, y conecta todo lo del dia. (1) la terapia mata la periferia oxigenada -> (2) sobrevive poblacion hipoxica, del 9% al 49% -> (3) las hipoxicas tienen PLASTICIDAD x3 (tejido.py:273) -> (4) fabrican reinas al triple de velocidad (desdiferenciaciones 28->109->185->367->837) -> (5) recaida. EL TRATAMIENTO SE FABRICA SU PROPIA RECAIDA, y no por fallar sino por FUNCIONAR: al eliminar lo accesible enriquece el tumor en la poblacion que mas reinas produce. Presion selectiva clasica, emergente, nadie la programo. Es literalmente la barrera que el cementerio dice que no ha resuelto nadie (celula lejos del vaso + quiescente), vista aparecer sola y con el mecanismo de enriquecimiento incluido. REINTERPRETA "reconducir": quitar la hipoxia rompe el paso 3 del bucle, de ahi los colapsos totales en algunas semillas (aunque cerca de la frontera decide el azar).



### 2026-10-22
 EL ORDEN IMPORTA PARA EL GRIFO PERO NO PARA EL DESENLACE. Cerrar la hipoxia ANTES (d90), a la vez (d120) o DESPUES (d150) de matar, 7 semillas, 400 dias: reconducir elimina la hipoxia por completo (fraccion 0,063 -> 0,000) en los tres. Victorias en DESDIFERENCIACIONES: antes 5/7 > a la vez 4/7 > despues 2/7 = dosis-respuesta TEMPORAL monotona en la direccion predicha (firma de causalidad real). Desdif medianas: A matar solo 455, B a la vez 122, C antes 182, D despues 345. PERO en peligrosas C gana solo 3/7 y pierde 1/7 = por debajo del umbral prefijado (5/7). TERCERA VEZ EN EL DIA que se mueve el mecanismo y NO se mueve el desenlace (anti-SASP, visibilidad, grifo).

### 2026-10-22 (2)
 **EL HALLAZGO MAS GRANDE DEL PROYECTO: NO HAY UMBRAL DE DEPOSITO Y EL TUMOR SE REGENERA DE NOVO.** Sonda de medida (borrar instantaneamente una fraccion del tumor SIN tocar nada sano; no es terapia, es sonda): dejar 100% -> 10.411 peligrosas, 50% -> 10.388, 20% -> 10.411, 10% -> 10.422, 5% -> 10.380, 1% (8 celulas) -> 10.413. BORRAR EL 80% O EL 99% DEL TUMOR NO CAMBIA NADA. Y con el 0% (erradicacion TOTAL verificada: 0 tumorales inmediatamente despues) 4 de 5 semillas lo REGENERAN. Parte (b): borrando TODAS las obreras y dejando las reinas, vuelve completo (10.483); borrando el 100% de obreras, igual.

### 2026-10-22 (3)
 VERIFICACION del hallazgo anterior (3 comprobaciones, las 3 pasan). (1) la sonda borra de verdad: 0 tumorales inmediatamente despues en 5/5. (2) La reaparicion es LENTA, no expansion de residuo: d150=4-6 celulas, d300=25-1365; el tumor original tardo 120 dias en llegar a 740 y el regenerado ~180 en llegar a ~1300 = velocidad comparable -> carcinogenesis DE NOVO. (3) Hay rastro del campo: las mutaciones acumuladas en celulas SANAS suben de 0,0 a 0,4-0,7. COHERENCIA INTERNA que nadie programo: la semilla 5 es la UNICA con 0,0 mutaciones en sanas y la UNICA que no reincide. (Bug propio detectado y anotado: la columna d400 del script era artefacto, el bucle llega a 1599 pasos y el dia 400 exacto nunca ocurre; datos validos hasta d300.)

### 2026-10-22 (4)
 CONSECUENCIA. EL PROBLEMA NO ES EL TUMOR, ES EL TEJIDO PREMUTADO QUE LO PRODUCE. Vaciarlo al 100% no cura: borras el producto y dejas la fabrica. Es CANCERIZACION DE CAMPO (Slaughter 1953). Explica de golpe: por que siempre hay respuesta y recaida; por que nada funciona por encima de ~4000 celulas (mas tiempo de campo activo, no mas volumen); por que las recurrencias en mama se reparten entre los anos 5 y 20; y por que el TAMOXIFENO es lo unico del cementerio con efecto a 15 anos -- es el unico que NO MATA CELULAS, apaga la senal del TEJIDO (EBCTCG 2011, ~21.000 mujeres, RR 0,53/0,68/0,97 por quinquenios).

### 2026-10-22 (5)
 ESTRATEGIA DE TRES CAPAS en marcha (idea del autor: "si no arreglamos lo que lo causa volvera a aparecer"). capa 1 vaciar (troyano+guardiana_nicho+memoria, PUNTUAL) + capa 2 cerrar grifo (reconducir, PUNTUAL) + capa 3 APAGAR EL CAMPO (p.estrogeno 1.0 -> 0.1, MANTENIDO). La palanca del campo YA EXISTIA en el modelo (experimento.py --hormonal) y llevabamos tres dias sin usarla. Diferencia clave: las capas 1 y 2 se agotan, la 3 no -- y el problema de todo el dia ha sido que la oleada se agota. Se mide tambien la RETIRADA del campo al dia 350: si al retirarlo vuelve el riesgo, el modelo reproduciria el RR 0,97 de los anos 10-14 sin que nadie se lo haya ensenado (seria la 5a coincidencia no programada).



### 2026-10-22 (6)
 TRES CAPAS, veredicto final con 20 SEMILLAS (13 nunca usadas para decidir nada). 2 capas (vaciar+grifo): 12/20 colapsos (<100), 3 intermedios, 5 fracasos (>8000), desv 4242. 3 capas (+campo): 12/20 colapsos, 6 intermedios, 2 FRACASOS, desv 3408. Metrica corregida: gana 5 / pierde 3 / EMPATA 12. LECTURA: anadir el campo NO mejora los casos buenos (12 vs 12 colapsos identicos) pero RECORTA LA COLA MALA (fracasos 5 -> 2). No sacrifica nada. Honestidad: 5 vs 3 con n=20 NO es prueba, es tendencia consistente con un mecanismo que entendemos; harian falta 50-100 semillas. EL NUMERO DEL PROYECTO: vaciar+grifo controla el tumor en 12 DE 20 SEMILLAS (60%), con validacion fuera de muestra.

### 2026-10-22 (7)
 BUG DE METRICA CORREGIDO (importante, afecta a todo lo medido antes). Un umbral puramente RELATIVO miente en los extremos: pasar de 1 a 8 peligrosas es +700% y no significa nada (ambas son control total). Con la metrica vieja el veredicto era 9 gana / 8 pierde = ruido indistinguible; con suelo absoluto de 100 aparece la senal real (5/3/12 empates). Anadidas mejor_que() y peor_que() en validacion.py con el suelo, documentadas con el caso que las motivo, y aplicadas tambien a T10. REGLA PARA EL FUTURO: nunca comparar con umbral relativo solo cuando los valores pueden acercarse a cero.



### 2026-10-22 (8)
 PREVENCION -- EL MODELO FALLA, y es un limite importante que ya esta documentado. Primer intento mal planteado (dos errores mios): (a) el escenario luminal SIEMBRA una celula alterada en tejido.py:189, asi que no media prevencion sino retraso de una lesion ya iniciada; (b) use estrogeno 1.0 como control cuando el valor natural del escenario es 0.6 (tejido.py:116 lo sobreescribe), lo que exageraba el efecto. Rehecho con TEJIDO VIRGEN (semilla={}) y control 0.6, 8 semillas, 700 dias: INCIDENCIA 8/8 en AMBOS brazos = 0% de reduccion; carga final 10.474 vs 10.478; el dia de aparicion pasa de 107 a 151 pero 2 de 8 semillas ADELANTAN.

### 2026-10-22 (9)
 DOS FALLOS DEL MODELO frente a la realidad, misma causa. (1) el tejido virgen desarrolla cancer en el 100% de los casos (riesgo real de por vida ~12%); (2) apagar el campo no previene (tamoxifeno en prevencion primaria reduce incidencia ~50%). CAUSA: las tasas de mutacion estan ACELERADAS para que pase algo en una ventana corta -- razonable para tratamiento, INVALIDA el modelo para prevencion, porque si todo muta rapido cualquier tejido acaba produciendo tumor. No se arregla tocando el estrogeno: habria que recalibrar la tasa de mutacion contra la incidencia real. Anadido como T14 en validacion.py, FALLANDO A LA VISTA a proposito. MATIZ PARA LA CAPA 3: de las tres propiedades reales de la hormonoterapia el modelo reproduce DOS (inutil en solitario, el riesgo vuelve al retirarla) y FALLA la tercera (prevencion). La estrategia de tres capas sigue valiendo porque se evaluo en TRATAMIENTO, que es el regimen donde el modelo funciona.

### 2026-10-22 (10)
 CONFIRMADO que la carcinogenesis de novo es real y no era artefacto: el tejido VIRGEN (sin ninguna celula sembrada) desarrolla carcinoma solo, en ~100 dias, 8/8 semillas. La celula sembrada del escenario luminal solo ACELERA lo que iba a pasar igual.

### 2026-10-22 (11)
 CUATRO ERRORES DE DISENO EXPERIMENTAL MIOS EN UN DIA, todos del mismo tipo (medir una cosa creyendo medir otra): (1) T2 atribuia al farmaco las muertes por necrosis del nucleo -- se arreglo con brazo control; (2) la metrica de victorias usaba umbral solo relativo y contaba 1->8 como derrota; (3) el experimento de prevencion partia de una celula ya sembrada; (4) ese mismo experimento usaba estrogeno 1.0 como control cuando el natural es 0.6. LECCION ADICIONAL DE INFRAESTRUCTURA: cuando un experimento puede interrumpirse (el gestor mato dos barridos por memoria), hay que (a) imprimir CADA semilla al terminarla y (b) ALTERNAR los brazos, para que cualquier prefijo del resultado siga siendo un experimento valido.



### 2026-10-22 (12)
 **EL MEJOR RESULTADO DEL PROYECTO: EXTIRPAR PRIMERO LO CAMBIA TODO.** El autor levanto su veto a la cirugia ("solo inyectar o extirpar, sin quimio, sin debilitar al sujeto") y eso abrio resecar(), que nunca habiamos usado en serio. Protocolo: resecar(cobertura 0.95, margen 2) el dia 120 -> troyano+guardiana_nicho+memoria dosis 300 -> reconducir -> campo apagado. 8 semillas, 500 dias. A. inyectar SIN extirpar: 0/8 erradica, 1/8 controlado, pico 4698, 0/380 dias controlado, finales [3276,117,2751,8849,10062,6121,6354,86]. B. EXTIRPAR+inyectar: 4/8 ERRADICA, 7/8 controlado, pico 28, 380/380 dias. C. protocolo completo: 4/8 erradica, 8/8 controlado, finales [16,0,42,28,0,0,22,0]. PRECIO: 210 celulas sanas extirpadas (~1,7% de un conducto de ~12.000) y CERO colateral de las terapias (T13 sigue pasando).

### 2026-10-22 (13)
 POR QUE funciona extirpar, y conecta con todo lo medido antes: el techo del proyecto era >4000 celulas no funciona NADA / <1000 funciona casi siempre. Las terapias celulares no fallaban por malas sino por no alcanzar el nucleo de una masa grande. Resecar deja un residuo de 0-15 celulas y ahi la MISMA terapia que fracasaba funciona. Convierte un problema de VOLUMEN en uno de SUPERFICIE -- literalmente lo que dice el docstring de resecar(), escrito antes de que supieramos que ese era el muro. SEPTIMA COINCIDENCIA con la realidad y la mas importante: el simulador redescubrio solo el esquema clinico estandar de los ultimos 100 anos (cirugia primero, sistemico sobre el residuo despues; nunca sistemico solo contra una masa grande).

### 2026-10-22 (14)
 HONESTIDAD sobre ese resultado: (a) la VIGILANCIA NO APORTO NADA -- retratamientos=0 en las 8 semillas, nunca se disparo el umbral (>=5 peligrosas en biopsia mensual), asi que B y C son el mismo brazo y la unica diferencia (110 vs 42 en una semilla) es ruido; hay que ajustar umbral o frecuencia; (b) son solo 8 semillas, repetir con mas; (c) el protocolo NO resuelve el campo: 4 de 8 no erradican y el tejido premutado sigue ahi.

### 2026-10-22 (15)
 EL JARDINERO (idea del autor: "curamos cancer y luego reparamos tejido"). Montado como MEDIDA DE REQUISITOS MINIMOS, no como terapia: barrido cobertura x eficacia x frecuencia, incluido un brazo deliberadamente imposible (100% de las celulas, eficacia 1.0, cada 10 dias) porque si ni eso evita el rebrote la via esta muerta y no hay farmaco que buscar. Diana = celula luminal con n_drivers==0 pero con golpes PARCIALES (un supresor necesita 2 golpes: con 1 parece sana y esta a un paso). CLAVE: el jardinero NO PUEDE SELECCIONAR porque una premutada no se distingue por fuera (el golpe esta en el ADN), asi que trata a ciegas a todo el que alcanza. Y el tejido sigue mutando mientras lo reparas: tasa_mutacion = 6e-3 por gen y division, ACELERADA ~6000x sobre la real (~1e-6), confesado en el comentario de biologia.py:61.

### 2026-10-22 (16)
 AVISO SOBRE "COSAS NATURALES" en quimioprevencion, por si vuelve a surgir: la via se probo en serio y salio MAL con dano medible. CARET (beta-caroteno + vitamina A en fumadores) se paro porque el cancer de pulmon AUMENTO 28%. SELECT (vitamina E, >35.000 hombres) aumento el cancer de prostata 17%. Los dos eran grandes, aleatorizados y bienintencionados, con la logica razonable de "antioxidante = menos dano al ADN = menos cancer". Lo que SI tiene evidencia en quimioprevencion de mama son los SERM (tamoxifeno, raloxifeno) e inhibidores de aromatasa, y todos actuan reduciendo DIVISIONES, no reparando ADN.



### 2026-10-22 (17)
 EL JARDINERO NO TIENE A QUIEN REPARAR, y obliga a CORREGIR el hallazgo del campo. Barrido cobertura x eficacia x frecuencia: peligrosas ~10.400 en TODOS los brazos incluido el imposible (100% de celulas, eficacia 1.0, cada 10 dias) -- y en ese brazo solo hubo 13 EDICIONES en 500 dias, con premutadas=0 en todos. CAUSA: con la tasa acelerada, la celula que recibe el primer golpe en un supresor recibe el segundo casi inmediatamente; el estado intermedio (parece sana, esta a un paso) dura un suspiro. CORRECCION IMPORTANTE: lo que medimos el 2026-10-22 (2-4) como "tejido premutado" (mutaciones en sanas 0,0 -> 0,4-0,7) eran mutaciones PASAJERAS (c.mutaciones, las que alimentan neoantigenos), NO golpes en genes driver (c.genoma). El reservorio de drivers a medias esta VACIO. Asi que el tumor que reaparecia tras borrar el 100% no venia de un tejido cargado sino de MUTACION DE NOVO a velocidad desbocada: no es cancerizacion de campo, es artefacto de calibracion.

### 2026-10-23
 CALIBRACION HECHA Y OPERATIVA. Barrido sobre tejido virgen (6 semillas, 700 dias, umbral 100 celulas): 6e-3 (6000x real) -> 6/6 incidencia, dia 110, carga 10.483; 2e-3 -> 4/6, dia 370, carga 584; 6e-4 -> 3/6, dia 470; 2e-4 -> 1/6 (~17%, cerca del 12% real), dia 500, carga 2; 6e-5 -> 0/6. NO HAY UNA TASA CORRECTA UNICA: depende de la pregunta. Anadido biologia.CALIBRACIONES {tratamiento 6e-3, intermedia 2e-3, prevencion 6e-4, realista 2e-4} con la tabla completa y el razonamiento en comentario, y expuesto como `experimento.py --calibracion`. T14 actualizado: ya no dice "falla y punto" sino que apunta a la salida.

### 2026-10-23 (2)
 OBSERVACION VALIOSA que emergio de la calibracion: cuanto mas realista haces la incidencia, MAS CARO es medir si algo la previene. A 2e-4 la incidencia (17%) se acerca a la real pero con 1 evento de 6 no se puede detectar una reduccion del 50%. Es exactamente por lo que los ensayos de prevencion reales necesitan decenas de miles de participantes y decadas. El simulador reprodujo sin querer la dificultad METODOLOGICA del campo, no solo la biologica.

### 2026-10-23 (3)
 Y el reservorio PREMUTADO nunca pasa de 14 celulas a NINGUNA tasa del barrido. El estado intermedio no se sostiene como poblacion, asi que "reparar el tejido" (el jardinero) no tiene a quien reparar y eso NO se arregla calibrando: haria falta cambiar la relacion entre tasa de mutacion y ritmo de division para que el primer y el segundo golpe se separen en el tiempo.



### 2026-10-23 (4)
 EPIDEMIOLOGIA (investigada a fondo). Lo que tienen en comun TODAS las pacientes: ser mujer + haber envejecido + un epitelio con decadas de ciclos hormonales. Riesgo x8 de los 30 a los 70 (ningun factor de estilo de vida se acerca). A IGUAL EDAD, premenopausicas RR 1,43 vs postmenopausicas = no es la edad, es el numero de ciclos. ~56% de los casos NO se explica por ningun factor clasico. Solo el 23% es prevenible (CRUK/Brown 2018: alcohol 8%, obesidad 8%, no lactancia 5%, THS 2%, radiacion 1%). Migrantes japonesas a EEUU: alcanzan la incidencia americana en 1-2 generaciones y tras ajustar por factores medidos el RR es 1,11 = el gradiente es AMBIENTAL y el "ambiente" es el patron reproductivo/metabolico.

### 2026-10-23 (5)
 CORRECCION CLAVE A LA IDEA DE "QUITAR ESTROGENOS". En el WHI el ESTROGENO SOLO **REDUJO** el cancer de mama (HR 0,78, 20 anos de seguimiento). El villano es la COMBINACION con progestageno: E+P diario RR 2,30, E+P intermitente 1,93, E solo 1,33. La progesterona actua via RANKL sobre las CELULAS MADRE mamarias -- justo el compartimento que el simulador ya tiene. PREVENCION QUE FUNCIONA: inhibidores de aromatasa -65% (MAP.3) con efecto que PERSISTE tras dejar el farmaco (IBIS-II, HR 0,51 a 131 meses); SERM -38% (solo RE+); mastectomia reductora en BRCA -90% o mas. LO QUE NO FUNCIONO en mama: dieta baja en grasa (WHI, n=48.835, HR 0,91 NS), vitamina D y omega-3 (VITAL), fenretinida. NADA basado en dieta o suplementos ha superado un ensayo aleatorizado.

### 2026-10-23 (6)
 LO QUE LE FALTA AL SIMULADOR segun la literatura (11 puntos; los 4 mejores). (1) el estrogeno tambien es GENOTOXICO DIRECTO (quinonas catecol-estrogeno, aductos depurinantes), no solo mitogeno: deberia multiplicar tambien la TASA DE MUTACION, no solo las divisiones -- y eso explicaria por que los IA (bloquean sintesis) superan a los SERM (bloquean receptor). (2) falta el PROGESTAGENO, dominante en THS, via RANKL sobre celulas madre. (3) la PARADOJA DEL EMBARAZO: HR 1,80 a los 5 anos post-parto y 0,77 a los 34 (Nichols, Ann Intern Med 2019) -- no se puede producir con un modelo que solo cuenta ciclos; requiere diferenciacion terminal + involucion post-parto (y el modelo YA tiene matriz, CAFs y macrofagos para implementarlo). (4) LA INVERSION DEL IMC (protector premenopausia RR 0,92, danino postmenopausia RR 1,12) es EL MEJOR TEST DE VALIDACION disponible: si al subir la grasa el simulador no cambia de signo al cruzar la menopausia, el acoplamiento hormonal esta mal. Ademas: el modelo de Pike 1983 ("edad del tejido mamario") es basicamente lo que el simulador reimplementa, y explica la inflexion de la curva de incidencia a los 50.

### 2026-10-23 (7)
 TNBC -- EL PANORAMA CAMBIO EN 18 MESES. TRES aprobaciones nuevas en 1a linea metastasica, todas ADC anti-Trop-2: TROPION-Breast02 (PFS 10,8 vs 5,6; OS 23,7 vs 18,7; ORR 63% vs 29%), ASCENT-03/04, PANKU-Breast02 (PFS 8,5 vs 3,1 HR 0,29). El TNBC ya NO esta desamparado. Y la noticia de 2026 son los ADC BIESPECIFICOS: iza-bren (EGFRxHER3) con fase 3 doblemente positiva (ORR 51,7%) y JSKN016 (TROP2xHER3) fase I con ORR 64,5%, mPFS 8,5 m, TRAE >=G3 solo 24,6% y NINGUNA EPI. Dos antigenos malos por separado que combinados amplian cobertura y reparten avidez (lo que ademas deberia mejorar la penetracion).

### 2026-10-23 (8)
 NO EXISTE DIANA EXCLUSIVA, y hay tres pruebas duras. (a) El MISMO scFv de sacituzumab (farmaco aprobado y seguro) montado en un CAR-T fue LETAL en ratones humanizados para TROP2, con dano pulmonar severo (Sun, JITC 2025) -- EL INDICE TERAPEUTICO VIENE DEL FORMATO, NO DE LA DIANA (el ADC no se amplifica, se dosifica intermitente, y necesita una concentracion intracelular que el tejido sano con expresion baja rara vez alcanza). (b) UPTIDER (autopsia rapida, 909 muestras, 30 pacientes, Nat Commun 2025): comparo 72 dianas de ADC -> TROP2 es LA MEJOR de las 72 por diferencial tumor/normal, B7-H4 la segunda, MUC1 clasificada de ALTO RIESGO. No hay exclusividad: hay un RANKING de diferenciales. (c) las mutaciones de perdida de funcion de TACSTD2 causan distrofia corneal = Trop-2 tiene funcion en epitelio normal. CONSECUENCIA: la industria y la FDA han ABANDONADO el modelo "diana = llave" -- sacituzumab y Dato-DXd se aprobaron en TNBC SIN exigir expresion de Trop-2, y T-DXd se aprobo en HER2-ultralow. Modelo nuevo: "la diana solo enriquece la entrega, el payload hace el resto en 50-200 um".

### 2026-10-23 (9)
 LOS DOS PATRONES QUE EXPLICAN CASI TODOS LOS FRACASOS DE ADC. (1) ES EL PAYLOAD, NO LA DIANA: los ADC con MMAE/auristatina clasica han caido uno tras otro (LIV-1, GPNMB, enfortumab en mama, B7-H4 felmetatug, NaPi2b, ROR1); los de inhibidor de TOPOISOMERASA I funcionan (Dato-DXd, sac-TMT, AZD8205, iza-bren, SYS6043). Seis fracasos y cinco exitos explicados por la quimica del payload. (2) "AUSENCIA EN TEJIDO NORMAL" ERA MENTIRA EN CASI TODOS: ROR1 (el dogma era FALSO -- los anticuerpos antiguos no eran sensibles; con el clon 6D4 aparece en PARATIROIDES, islotes pancreaticos, esofago, estomago, duodeno; siete farmacos muertos), GPNMB, MUC1, LIV-1, CEACAM5. El campo se construyo sobre mapas de tejido normal incompletos.

### 2026-10-23 (10)
 LA LECCION DE FRalfa -- EL RIESGO MAS CARO Y MENOS MODELADO. La literatura repetia "~80% de los TNBC sobreexpresan FRalfa". Ensayo real (Yam, Invest New Drugs 2021): 96 pacientes cribadas, IHC en 80 -> TASA REAL 10,0% (8/80), solo 2 recibieron el farmaco, ORR 0%, terminado precozmente. Se repite: CEACAM5 24% real, TROP2 56% con los cortes dicotomicos de los ensayos frente al 82-88% publicado, NaPi2b con ARNm alto pero PROTEINA infraexpresada vs normal. REGLA: medir la prevalencia con EL ENSAYO Y EL CORTE QUE SE VAN A USAR, no con los de la literatura. Debe entrar en el simulador como parametro propio: prevalencia efectiva en cribado prospectivo != prevalencia publicada.



### 2026-10-23 (11)
 ADC IMPLEMENTADO en el modelo -> tejido.liberar_adc(diana, afinidad, prevalencia, dosis, radio_espectador, depende_ciclo) + biologia.PREVALENCIA_DIANA. NO es un agente que camina: es una MOLECULA que difunde desde el vaso, con tres propiedades medidas en humanos: (1) penetracion limitada -- alta afinidad ~20 um del vaso, baja ~85 um ("barrera del sitio de union": agarrarse fuerte a la primera celula impide llegar al fondo); (2) el payload se DERRAMA 50-200 um y mata vecinas SIN la diana (es lo que resolvio la heterogeneidad y por lo que la FDA aprueba sin exigir expresion); (3) depende del CICLO CELULAR. PREVALENCIA_DIANA usa las cifras MEDIDAS en cribado prospectivo, no las publicadas (trop2 0.56 no 0.82-0.88; folato 0.10 no 0.80).

### 2026-10-23 (12)
 DOS HALLAZGOS GEOMETRICOS DEL ADC, no buscados. (a) EL CARCINOMA IN SITU ES INALCANZABLE POR GEOMETRIA: la celula tumoral mas cercana a un vaso esta a 60 um, la mediana a 210 um, la maxima a 405; con alta afinidad (20 um) es alcanzable el 0% y con baja (90 um) el 8,2%. El tumor dentro del conducto esta separado del estroma vascularizado y se nutre por difusion. ESO EXPLICA POR QUE LOS ADC SE USAN EN ENFERMEDAD AVANZADA Y NO EN IN SITU, y el modelo lo produjo solo. (b) La ANGIOGENESIS tardo 162 dias en aparecer y en cuanto aparecio el tumor paso de 1.179 a 10.255 celulas: el acceso vascular es lo que lo desata. IRONIA: el ADC solo puede llegar cuando el tumor ya esta en el estado (grande y mayoritariamente quiescente) en que el ADC funciona peor.

### 2026-10-23 (13)
 LA SUPOSICION DEL CICLO CELULAR NO DETERMINA EL DESENLACE (control metodologico importante, hecho A PROPOSITO antes de concluir nada, para no repetir el error de contencion_matriz=0.8). 4 ciclos de ADC cada 21 dias desde el dia 250, 5 semillas, 500 dias: sin ADC 10.420 peligrosas; ADC que NO mata quiescentes 10.419 (mato 978); ADC que SI mata quiescentes 10.472 (mato 13.173); alta afinidad 10.447. EL BRAZO QUE MATA 13.173 CELULAS DEJA EL MISMO TUMOR QUE EL QUE MATA 978. Es decir: el simulador NO devuelve mi suposicion, y ademas es la SEXTA vez que se mueve el mecanismo sin mover el desenlace. OCTAVA COINCIDENCIA CON LA REALIDAD: los ADC tampoco curan en clinica -- las tres aprobaciones nuevas en TNBC metastasico dan PFS 10,8 meses y OS 23,7 meses. Compran tiempo, no curan, y el modelo lo reprodujo sin que se lo ensenaramos.


### 2026-10-28
 FASE 1 DISEMINACION (diseminacion.py --d4, 32 semillas): D4 FALLA -> fase 1 NO LISTA. D2 (dormidas) es el supuesto devuelto (prob_entra_dormida 0,70 -> D2 58 %); D3 ~ 1-prob_sobrevive_transito; D1 base 38 % (en diana) pero recorre 25-50 % con prob_sobrevive_transito; competencia_temprana casi no pesa. Siguiente: que la dormancia salga de un mecanismo (nicho Ghajar: microvaso estable/brote) y no de un parametro.

### 2026-10-28 (2)
 T15 en validacion.py (anti-PD-1, inspirado en la aprobacion FDA del dostarlimab en recto dMMR, AZUR-1) FALLA: 0/8 mejoras y 0 RC en dMMR (prob_neoantigeno 1.0) y proficiente (0.1); el farmaco quita 35.000 frenos pero no hay expansion clonal T (tope 40 T+NK). Anadido tejido.bloquear_pd1() sin cambiar defaults. Falta: expansion clonal, PD-L1 adaptativo por IFN-g, reservorio T exhausto.

### 2026-10-28 (3)
 INMUNIDAD ADAPTATIVA (expansion clonal Kaech/Ahmed, PD-L1 adaptativo Taube/Ribas, Tex progenitores revigorizables Siddiqui/Miller) tras flag Parametros(inmunidad_adaptativa=True), APAGADO por defecto (9 SUPUESTOS). T15 sigue FALLA 0/8 (el farmaco frena durante el bloqueo pero el dia 500 esta saturado; RC imposible porque las reinas presentan neo//3 -> invisibles). Bateria sigue pasando, T12 a una semilla de fallar. Posible artefacto: los T se quedan pegados a dianas PD-L1+ (~200.000 frenos).

### 2026-10-28 (4)
 DORMANCIA POR NICHO VASCULAR (v2 diseminacion.py; fraccion vasos tipo H 0,16 medida, Wang 2017; estado sigue al vaso, reversible). D2 agregado 75 % (FALLA, dentro del ruido), D4 NO CONCLUYENTE: solo ~24 DTC -> cada una pesa 4 puntos. Despertar irreversible -> D2 37 %. Fase 1 NO LISTA; falta resolucion (128-256 semillas o medir todas las DTC que llegan). Copia v1 en un script exploratorio, no incluido.

### 2026-10-28 (5)
 ARTEFACTO "T PEGADOS" CONFIRMADO Y CORREGIDO: cazar() sin memoria volvia siempre a la misma diana PD-L1+. Sin flag adaptativo era menor (5 frenos en 30 h, veredictos T1-T15 iguales tras corregir). Con flag los T quedaban anclados hasta 41 dias (94 % frenos repetidos) -> todo lo medido con inmunidad_adaptativa=True antes de la correccion queda en duda. Fix: biologia.dias_desenganche_pd1=2.0 (SUPUESTO), el T no reintenta la celula que le freno. Aun con fix, con flag los T matan casi nada: lo limita el halo de PD-L1 inducido por IFN-g. Con flag on T13 FALLA (1 sana comida por guardiana puerta nicho, no por T): falso positivo posicional pendiente. tejido3d.py tiene el mismo artefacto sin corregir.

### 2026-10-28 (6)
 T13 ERA FRAGIL TAMBIEN SIN FLAG: con 24 semillas fallaba 3/24 (5 sanas comidas por la guardiana puerta nicho); pasaba solo porque las semillas 1-3 no caian en el caso. Causa: madres SANAS vivian fuera de la capa basal con marcador de madre. Fix en tejido.py _madre_sana_sin_nicho (identidad de madre sana depende del nicho basal: division orientada, diferenciacion al perder nicho). T13 0/24 con y sin flag; eficacia sin cambio (12/16 -> 11/16 controlados). OJO: T1 paso de 240 a 330 um (limite 360) sin explicar. Limite: en mama real hay progenitores luminales ALDH1+ suprabasales que una puerta real tocaria. LECCION: 3 semillas no bastan para el test de seguridad.

### 2026-10-28 (7)
 T1 SIGUE VALIDADO: 240->330 um era ruido de 3 semillas (solo 1 con necrosis) + metrica dependiente del tamano. Con 48 semillas: necrosis 255 vs 270 um, hipoxia 180 vs 180; a igual tamano identico; consumo O2 no cambia. Propuesta metrica robusta por anillos (frente hipoxia 150 um; necrosis umbral 25 %) en (script exploratorio /t1/CRITERIO_T1_ROBUSTO.txt, no incluido), NO implementada. T13 ahora con minimo 24 semillas (T13_SEMILLAS_MIN), PASA 0/24. EFECTO LATERAL SIN INVESTIGAR: tras _madre_sana_sin_nicho los tumores a dia 250 son la MITAD (mediana 2327 -> 1172 celdas).

### 2026-10-28 (8)
 T1 CAMBIADO A METRICA POR ANILLOS (validacion.py, minimo 12 semillas, paralelo, ~20 s). PASA: frente de hipoxia 11 celdas = 165 um (diana 150-180, clavado); inicio de necrosis (>=50 % del anillo) 23 celdas = 345 um, a 15 um del limite 360. Perfil limpio y monotono: necrosis aparece a 14 celdas (210 um), 25 % a ~19 (285), 50 % a 23 en anillos ya escasos (n=55). La transicion necrotica del modelo es gradual (~135 um) y mas lejana que la realidad; el umbral 50 % cae en la cola con poca muestra -> fragil. No se cambio el criterio tras ver el dato.

### 2026-10-28 (9)
 **DEFECTO: LAS HIPOXICAS SON CASI INMORTALES.** Solo mueren con O2<0,04 (p=0,08/paso, vida 3,1 d); entre 0,04 y 0,15 riesgo ~0 (vida ~1500-2500 d, menos que las oxigenadas). Real: hipoxicas humanas viven 4-10 d en xenoinjerto (Durand & Sham IJROBP 1998); a 0,5 % O2 no hay efecto (Papandreou 2005). Sondas sin cambiar el modelo (script exploratorio /o2, no incluido): A vida 7 d bajo 0,15; B riesgo lineal 0->0,08 entre 0,15 y 0,04. Con A o B la fraccion hipoxica cae a 0-0,13. EN DUDA: "la terapia selecciona el refugio hipoxico (9->49 %)" (a igual tamano un tumor sin tratar ya esta al 54-58 %: efecto pre-angiogenico, no seleccion) y "el tratamiento se fabrica su propia recaida" (solo 2/8 semillas fallan y saltan con cada variante). T1 pasa por Manhattan + retirada de necroticas (euclidea sin retirada: necrosis 50 % a 240 um). Pendiente decision del autor: A o B, metrica T1 euclidea, re-medir refugio a igual tamano con 24 semillas.

### 2026-10-28 (10)
 VARIANTE B APLICADA (eleccion del autor): necrosis por paso lineal 0 en O2=0,15 -> 0,08 en 0,04 (Durand & Sham 1998, Papandreou 2005), tejido.py ~269-290. T1 ahora con distancia EUCLIDEA (validacion._distancia_euclidea_a_vaso): hipoxia 150 um, necrosis 210 um, PASA. Bateria igual (T2,T7,T10,T12,T13 0/24 PASAN; T15 FALLA). Refugio con 24 semillas: modelo viejo cura 16/24, B 17/24; recaidas B 9-53 peligrosas (antes 296-8189). VEREDICTO: "la terapia selecciona el refugio" -> exceso real +0,16-0,40 sobre igual tamano solo en el modelo viejo, desaparece con B; "el tratamiento se fabrica su recaida" REFUTADO: las recaidas las predice al 100 % un tumor grande e hipoxico YA en el dia 120 (antes de tratar). Pendiente: tejido3d.py:146 sigue con la regla vieja; DISEMINACION.md hay que re-medirla (tumores algo menores). Salidas en (script exploratorio /o2/., no incluido)

### 2026-10-28 (11)
 CONTRASTES DE OPTIMISMO (T16-T18 nuevos en validacion.py; criterios en (script exploratorio /optimismo/CRITERIOS.md, no incluido). NO hay optimismo general tras la variante B. T16 duplicacion sin tratar: 21 d vs 100-260 d reales -> FALLA, tiempo comprimido ~10x (pesimista). T17 quimio neoadyuvante luminal: pCR 0/16, respuesta 0/16 aunque mata ~2x el tumor basal (repuebla a la vez) -> FALLA pesimista. T18 cirugia sola (calibracion 2e-3, no se ejecuta por defecto): recaida atribuible +8 % vs 35 % EBCTCG -> FALLA OPTIMISTA: resecar() corta un disco sobre foco unico compacto = R0 perfecta (residuo mediana 2 celulas); sin multifocalidad, extension intraductal, margenes positivos. Residuo->recaida: 2->17 %, 32->50 %, 80->88 %, 161->100 %. CONSECUENCIA: "extirpar primero lo cambia todo" (2026-10-22 (12)) se apoyaba en esa R0 irreal -> en duda. Opciones planteadas al autor: declarar limite / anadir residuo realista (focos satelite) / NO calibrar cobertura (seria programar el resultado). El modelo no guarda linaje madre-hija.


### 2026-10-29
 **LA CAUSA DEL FALLO DE T18 ERA FALSA; resecar() NUNCA FUE UNA TUMORECTOMIA.** Tres sondas que no miran el desenlace (script exploratorio cirugia/, no incluido): (a) geometria.py -- el tumor ocupa el 90 % del anillo ductal en 5 componentes conexas, o sea el modelo SI tiene multifocalidad y no habia que anadirla (lo escrito el 28 era falso); (b) que_corta.py -- resecar(0.95,2) extirpa el 93 % del CONDUCTO y el 49 % del campo, porque el tumor rodea la luz, el centroide cae en el centro del conducto y el disco (radio ~28) engloba el conducto entero (radio 22); (c) sonda3 -- lo mismo el dia 120 con tasa por defecto: 96 %. Censo del campo: 1.251 celulas = 188 luminales + 224 mioepiteliales + 839 fibroblastos; el conducto son 412, no las "~12.000" que decia la nota del 22-10 (comparaba celulas del modelo con un conducto humano real). LIMITE DE ESCALA: campo 1,65 mm, conducto 330 um de radio, margen quirurgico clinico 1-2 mm = 67-133 celdas > semilado del campo (55): ninguna cirugia con margen realista cabe, la TASA absoluta de recaida no es contrastable.


### 2026-10-29 (2)
 T18 REESCRITO COMO RAZON, NO TASA (eleccion del autor entre 3 opciones). Criterio cerrado antes (script exploratorio cirugia/CRITERIOS_T18_v2.md, no incluido) y medido con semillas 101-148 NUNCA usadas, para que fuera ciego: brazo resecar(0.95,2), margen positivo (residuo>0) n=19 recae 32 %, margen NEGATIVO (residuo=0) n=23 recae 0 % -> RR = infinito (15,6 con correccion de continuidad) frente al OR 2,44 de Houssami 2010 (Eur J Cancer 46:3219, 28.162 pacientes). Control de novo (vaciado total) 7 %: la senal no esta tapada. VEREDICTO: FALLA "el margen lo decide todo", con causa ESTRUCTURAL nombrable: en el modelo margen negativo cura siempre, en la clinica una conservadora sin RT con margen negativo recae 20-25 % porque queda enfermedad en el resto de la mama. FUERA DEL CONDUCTO NO HAY MAMA: el campo tiene UN sistema ductal, la mama real 15-20. El modelo no es optimista ni pesimista: LE FALTA EL ORGANO. Lo que SI reprodujo sin ensenarselo: el signo y la dosis-respuesta (residuo 0 -> 0 %, 31 -> 57 %) y 45 % de margenes positivos vs 20-25 % clinico (mismo orden). NO se corrige calibrando cobertura (seria programar el resultado) ni sembrando focos fuera del conducto (exige campo mayor = cambio de modelo).


### 2026-10-29 (3)
 "EXTIRPAR PRIMERO LO CAMBIA TODO" (22-10) REINTERPRETADO, NO RETIRADO (eleccion del autor). Los numeros siguen (4/8 erradica con extirpar+inyectar vs 0/8 sin extirpar; volumen -> superficie). Lo que cambia es QUE operacion es: no "una cirugia conservadora basta" sino "quitar el organo entero del campo basta". El precio real no es 1,7 % de nada: son ~210 sanas sobre un epitelio luminal de 188 celulas. Escrito en README.md (seccion nueva "Extirpar primero rompe el techo -- y que es en realidad esa cirugia"), en el docstring de tejido.resecar() (aviso de interpretacion) y en validacion.T18. PENDIENTE: los scripts de las sondas eran exploratorios y no se incluyen.


### 2026-10-29 (4)
 DEUDA DE tejido3d.py PAGADA (eleccion del autor). (a) NECROSIS VARIANTE B llevada al 3D (antes solo moria la celula con O2 < o2_necrosis a p=0,08; ahora rampa lineal igual que tejido.py). Comprobado: 120 d, semilla 1 -> 8.326 necrosis y 2.291 necroticas en la BANDA hipoxica, que con la regla vieja era territorio de celulas casi inmortales. Las dos mallas ya envejecen igual. (b) DESENGANCHE PD-1 escrito tambien en tejido3d._atacar (agentes.politica_linfocito_t YA lo consultaba pero el 3D no lo escribia nunca). OJO, HONESTIDAD: en 3D NO HAY RECLUTAMIENTO INMUNE -- los agentes solo entran por liberar_terapia y son tipo "IA", asi que no existe ningun agente tipo "T" y el artefacto NUNCA pudo manifestarse en 3D. La correccion es correcta y necesaria para cuando el 3D recluta inmunes, pero NO cambia nada de lo medido en 3D hasta hoy. No vender el fix como si arreglara una medida.


### 2026-10-29 (5)
 **AGUJERO EN EL FIX DEL 28: EL CAR-HER2 SIGUE PEGADO.** El desenganche del 28 solo se escribe si ag.tipo == "T", y los agentes de terapia son tipo "IA". De los que pueden ser frenados por PD-L1, car_her2 es el UNICO sin resistente_pdl1 (cazador, hibrido y reina si lo llevan) y ademas politicas.car_her2() no consulta ag.memoria["frenado_por"]. MEDIDO (2D, luminal, car_her2 dosis 300 el dia 120, 220 d, semillas 1-3): semilla 2 -> 95 frenos PD-L1, el 80 % en pasos CONSECUTIVOS sobre el MISMO par (agente, celula), racha maxima 25 pasos = 6,25 dias; semillas 1 y 3 sin frenos. Magnitud pequena frente a los ~114.000 frenos del caso con flag adaptativo, pero es el MISMO defecto. PENDIENTE DE DECISION DEL AUTOR: arreglarlo (que el desenganche valga para cualquier agente frenado y que car_her2() filtre por frenado_por) cambiaria resultados ya publicados con car_her2.


### 2026-10-29 (6)
 AGUJERO DEL CAR-HER2 CORREGIDO (el autor dijo que si). El desenganche de tejido._atacar ya no depende de ag.tipo == "T": vale para CUALQUIER agente frenado por PD-L1 (en la practica solo car_her2, el unico agente de terapia con receptor y sin resistente_pdl1; confirmado con censo de politicas.TERAPIAS). politicas.car_her2() filtra ahora por frenado_por, y tejido3d.py lleva la misma generalizacion. ALCANCE REAL: car_her2 NO se usa en ningun test de validacion.py ni en ningun informe .md -- solo existe en politicas.py -- asi que NO toca ningun resultado publicado (mi aviso previo de que si lo tocaba era exagerado y se corrigio). ANTES/DESPUES, 12 semillas, dosis 300 el dia 120, 220 d: consecutivos sobre el mismo par 703 (45 %) -> 0; racha maxima 25 pasos -> 0; muertes terapia 467 -> 450; peligrosas mediana 613 -> 620.


### 2026-10-29 (7)
 **MI EXPECTATIVA FALLO Y SE DEJA ESCRITA SIN REESCRIBIRLA** (EXPECTATIVA_CAR.md del (script exploratorio cirugia/, no incluido). Predije MENOS frenos totales tras el fix; hubo MAS: 1.556 -> 1.606 (+3 %). Causa: el CAR-HER2 desenganchado SE MUEVE y encuentra OTRAS celulas PD-L1+ que lo frenan a su vez -- menos frenos repetidos pero mas dianas distintas. Efecto sobre el tumor NULO y sin direccion: la semilla 4 empeora (667 -> 872) y la 6 mejora (6.096 -> 5.616). SEPTIMA VEZ en el proyecto que se corrige un MECANISMO y el DESENLACE no se mueve (el mismo patron que el ADC y el ciclo celular). La leccion se repite: en este modelo el desenlace lo fija el tamano/estado del tumor, no el detalle del mecanismo.



### 2026-11-05
 **SIN BISTURI (avances 2025-2026: Kuerer JAMA Oncol 2025 0/31 recaidas sin operar tras pCR por biopsia; FIRST crioablacion 89,6 % completa, CDIS en borde 12,5 %; RT ablativa+HT ESTRO 2026 53 % pCR; RESPONDER FNR biopsia 17,8 %).** Anadido tejido.ablacionar(cobertura, margen, borde, supervivencia_borde) (misma geometria que resecar, deja necrosis en vez de retirar) y tejido.biopsia(cilindros=6, celulas_por_cilindro=12) (FNR emergente, no programado). Banco sin_cirugia.py (luminal, dia 120, arsenal troyano+guardiana_nicho+memoria 300, 500 d, 24 semillas, 6 procesos, 43 min). Control <100 / erradicacion =0: nada 0/24; solo arsenal 20/24 control, 0/24 erradica; cirugia+arsenal 24/24, 13 a 0; ABLACION borde letal = CIRUGIA (24/24, 13 a 0, mismas 324 sanas destruidas); ablacion borde 30 % vivo 24/24 pero solo 7 a 0; via Kuerer (arsenal -> biopsia d180 -> ablacion solo si positiva) 23/24 control, 10/24 SIN tocar la mama, sanas mediana 247. PERO la biopsia fallo en las 10 negativas (todas tenian tumor vivo: FNR 10/24 = 42 % vs 17,8 % real); 9 de esas 10 se controlaron igual (la memoria residente limpia el resto), 1 recayo (1497). Conclusion: en el modelo el bisturi no es lo que cura, lo es destruir el volumen; la puerta para no tocar nada es una prueba de residuo mejor que la biopsia.


### 2026-11-05 (2)
 PRUEBA DE RESIDUO (prueba_residuo.py + tejido.ctdna). CORREGIDO un fallo mio: biopsia() contaba como positivo/tumor cualquier celula con driver (iniciadas invisibles al microscopio); ahora solo carcinoma (IN_SITU/INVASIVA/REINA). Foto dia 180 (60 d tras arsenal), 24 semillas: 18/24 con carcinoma residual, 11 de ellas con 1-4 celulas. FNR biopsia por vacio: 6 cil 74 %, 12 66 %, 24 52 %, 48 37 %, 96 28 %. ctDNA una extraccion (calibrado a sensibilidad pre-tratamiento 50/75/90 % -> 2,48e-2 / 9,91e-2 / 2,28e-1 fragmentos/celula): FNR 69/51/37 %, 0 falsos positivos. Protocolos 500 d (sin_cirugia.py, brazos nuevos): Kuerer biopsia corregida 20/24 sin tocar, 23/24 control pero residuos persistentes 1-39 y 1 recaida (1497), FNR 12/16; Kuerer+ctDNA mensual 24/24 control, 14 a 0, pero 0/24 sin tocar (todas acaban ablacionadas, rescates del dia 180 al 480); solo ctDNA mensual (sin pinchar) 24/24 control, 15 a 0, 2/24 sin tocar, sanas 292. LECTURA: el ctDNA seriado acaba encontrando algo en todas porque el campo genera lesiones de novo (tasa de mutacion acelerada 6000x): los positivos tardios son probablemente lesiones NUEVAS, artefacto a medir con la tasa calibrada. PENDIENTE: ablacionar() y resecar() apuntan al centroide de TODAS las celulas con driver (iniciadas incluidas) -> el disco engloba el conducto entero; una ablacion guiada por imagen deberia apuntar solo al carcinoma visible y destruiria mucho menos sano.


### 2026-11-05 (3)
 ABLACION GUIADA (ablacionar(guiada=True): solo carcinoma visible, union de discos de margen 2, borde 1; SUPUESTO OPTIMISTA: la imagen ve focos de pocas celulas) + TASA REALISTA tras tratar (6e-3 hasta el dia 120, CALIBRACIONES["realista"] 2e-4 despues). 24 semillas, 500 d. Control/erradica/sanas destruidas/sin tocar: ctDNA+disco 24/15/292/2; ctDNA+GUIADA 24/3/14/2 (10-20x menos dano pero el campo invisible sobrevive y deja restos de 1-23); REAL solo arsenal 23/4/0/24 (1 fallo 354); REAL ctDNA+disco 24/17/285/6; **REAL ctDNA+GUIADA 24/24 control, 12/24 a 0, 8 sanas destruidas (vs 324 de la cirugia), 6/24 sin tocar nunca** = mejor equilibrio del proyecto. Con tasa realista los sin tocar suben de 2 a 6: parte de los positivos tardios eran lesiones de novo (confirmado). Mismas semillas = mismos dias de rescate en disco y guiada (diseno emparejado limpio). Disyuntiva que queda: disco erradica mas (mata el campo invisible) vs guiada casi no dana.

### 2026-11-05 (7)
UMBRAL DE IMAGEN (issue #1). La imagen solo ve focos de >=N celulas (tejido.focos_carcinoma, componentes conexas); si el ctDNA da positivo y no se ve nada, solo se repite el arsenal. Tasa realista tras tratar, 24 semillas. Control/erradica/sanas/sin tocar: >=5 24/5/0/18; >=25 24/4/0/20; >=100 24/4/0/22; ciega 23/4/0/24. El control NO depende de ver lo invisible. La ciega no mejora al arsenal solo (23/24 en los dos): repetir el arsenal no aporta; la paciente que falta (semilla 6, tumor de ~900) se rescata solo si se ve y se ablaciona su foco grande. Con imagen perfecta se erradica a 12; con realista a 4-5, el resto conviven con restos de 2-11 celulas contenidos por la memoria inmune (control, no curacion; 380 dias de seguimiento no bastan para saber si se sostiene).

### 2026-11-06
SIN BISTURI EN 3D (sin_cirugia_3d.py, criterio escrito antes de ejecutar; luminal, tratamiento dia 60, 240 d, 8 semillas, tasa realista desde el 60, ctDNA calibrado en 3D 2,25e-3 a sensibilidad 75 %). Control = <5 % de la mediana sin tratar (4.855). Control/a 0/sanas mediana: nada 0/0/0; solo arsenal 0/0/0 (96.756 vs 97.102: los farmacos solos son inutiles en 3D, muro de superficie medido sobre el desenlace); disco+arsenal dia 60 4/0/8.763; **ctDNA+GUIADA 8/8, 2 a 0, 3.800 sanas**; ctDNA+disco 1/0/5.980. Todas las pacientes dan ctDNA positivo el dia 90 (0/8 sin tocar). VEREDICTO: AGUANTA. La ventaja no es de calendario: el disco el mismo dia 90 solo controla 1/8; lo que cuenta es apuntar a cada foco (el disco al centroide deja focos fuera y en las semillas 4 y 8 el tumor rellena el conducto por encima del brazo sin tratar). En 3D la guiada dana mas que en 2D (3.800 vs 8; incluye mioepiteliales y fibroblastos). Contratiempo: el banco se corto 3 veces por falta de memoria del ordenador; anadida cache en disco (resultados/cache_3d/, fuera del repo) para retomarlo. Limitaciones: 8 semillas (<12), sin reclutamiento inmune en 3D, imagen que ve focos de 1 celula. Captura 07 (paciente 4 el dia 90: disco 6.045 sanas y 365 peligrosas fuera; guiada 4.199 y 0).

### 2026-11-06 (2)
UMBRAL DE IMAGEN EN 3D (sin_cirugia_3d.py --umbrales, criterio escrito antes; mismas 8 semillas; ctDNA+ sin nada visible -> solo se repite el arsenal). Control/a 0/sanas mediana: >=1 8/2/3.800; >=5 8/3/3.778; >=25 8/3/3.774; >=100 8/3/3.774 (identico a >=25 celula a celula: el dia 90 los focos 3D son grandes); CIEGA 0/0/0 (96.630 vs 97.102 sin tratar). VEREDICTO: el control NO depende de ver focos diminutos (>=100 da 8/8, minimo 7/8). Contrario a 2D (ciega 23/24): en 3D la ablacion es imprescindible porque los farmacos no llegan al nucleo. Semilla 5 3.705 con >=5 y 320 con >=25: azar, no efecto. El banco se corto otra vez por memoria (otros proyectos abiertos); terminado con 1 proceso (12.923 s).
