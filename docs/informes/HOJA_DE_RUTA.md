# Hoja de ruta: de la hipótesis a la paciente

> **Actualizado el 2026-10-21. La estrategia ha cambiado.**
>
> La estrategia original era: *curar la célula madre tumoral ("reina") y reprogramarla
> como guardiana que elimina a las obreras con una puerta lógica Y*. Se ha medido a fondo
> y **ese brazo no funciona**: la guardiana fagocita 5 células de 10.029; arreglar el bug
> que hacía incurable al 91,5 % de las reinas multiplicó las guardianas por 2,65 y dejó el
> tumor igual o peor; y el efecto espectador de la corrección, a toda potencia, triplicó
> el tumor grande. Detalle en el README, sección *Curar reprogramando: acta de defunción*.
>
> **ACTUALIZACIÓN 2026-10-22 — la estrategia cambia otra vez, y es el cambio grande.**
>
> Se midió lo que pasa al borrar el **100 %** del tumor sin tocar una célula sana: **4 de 5
> semillas lo regeneran de novo**, desde un tejido premutado (mutaciones en células sanas
> de 0,0 a 0,4-0,7). No hay umbral de depósito: dejar 8 células da el mismo resultado que
> dejar 740. **El problema no es el tumor, es el tejido que lo produce** — cancerización de
> campo (Slaughter 1953).
>
> Eso añade una **tercera capa** a la estrategia, y es la única que no se agota:
> **apagar el campo** (`p.estrogeno`, es decir tamoxifeno o inhibidor de aromatasa),
> **mantenido en el tiempo**. El modelo reprodujo solo dos propiedades reales de esa capa:
> sola no sirve (es adyuvante) y retirarla devuelve el riesgo.
>
> Número del proyecto: **vaciar + cerrar el grifo controla el tumor en 12 de 20 semillas
> (60 %)**, con 13 semillas fuera de muestra. Añadir el campo recorta la cola mala
> (fracasos de 5 a 2) sin tocar los casos buenos.
>
> **Estrategia actual**: (1) **matar** con células terapéuticas persistentes que entran
> con vehículo (troyano) y usan una puerta Y que no toca tejido sano; (2) **cerrar la
> fuente** de células madre tumorales, que tiene tres grifos: división de las reinas,
> plasticidad por hipoxia, y SASP de las senescentes; (3) dejar memoria residente.
> Sin cirugía ni quimioterapia.
>
> Lo bueno del cambio: la estrategia ya **no depende de ninguna pieza inventada**. CAR-T
> persistentes, normalización vascular e inhibidores JAK existen los tres.

## Fase A: computacional (lo que se puede hacer desde un ordenador)

| Paso | Estado | Resultado |
|---|---|---|
| A1. Modelo por agentes 2D con genes reales | ✅ hecho | `main.py`, `experimento.py` |
| A2. Contraste con ~1.000 pacientes reales (TCGA) | ✅ hecho | `datos_reales.py` → `INFORME_DATOS_REALES.md` |
| A3. Calibrar neoantígenos con datos reales | ✅ hecho | con `--neoantigenos 0.1` la estrategia sigue controlando lo peligroso (749 → 2,7) |
| A4. Calibrar el resto de parámetros con literatura (ciclo celular, plasticidad, fracción de células madre CD44⁺/CD24⁻) | pendiente | semanas |
| A5. Pasar a 3D con PhysiCell + PhysiBoSS (redes de señalización reales por célula) | pendiente | 1-3 meses |
| A6. Simular cada paciente con SU genoma (gemelo digital) y su propio juego de vectores | pendiente | meses |
| A7. Análisis de sensibilidad: ¿qué parámetro, si estuviéramos equivocados, rompe la estrategia? | ✅ hecho | ver abajo |
| A8. Reproducción: telomerasa y contador de divisiones | ✅ hecho | ver abajo |
| A9. Barrera física: presión intersticial, malla, CAFs | ✅ hecho | `INFORME_BARRERA.md`; tumbó los resultados previos |
| A10. Instrumentar en qué acaba cada reparación | ✅ hecho (2026-10-21) | 91,2 % estériles, 8,3 % curan |
| A11. ¿Sirve curar reprogramando? | ✅ **respondido: no** | tres vías independientes, todas iguales o peores |
| A12. Cerrar la fuente: los tres grifos de reinas | ✅ hecho (2026-10-21) | el SASP fabrica reinas; antídoto real (JAK) |
| A13. ¿Por qué `reconducir` cura a unos y empeora a otros? | ✅ **respondido: caos** | ±1 célula de dosis cambia 12 en 10.102; no puede existir biomarcador |
| A14. Mapa del cementerio: qué se intentó y contra qué murió | ✅ hecho (2026-10-21) | `INFORME_CEMENTERIO.md`; 11 barreras, 2 resueltas |
| A15. Batería de validación: reproducir 12 fracasos conocidos | **siguiente** | es la prioridad: convierte el simulador en algo falsable |
| A16. Coste de apilar intervenciones (hoy es gratis = mentira) | pendiente | serie de Kurzrock: 51 % / 38 % / 28 % a dosis plena |
| A17. ¿El límite es matar o VER? Las reinas tienen neoantígeno efectivo 0 | en curso | `tejido.py:972` hace `neo //= 3` para madres |

**La prioridad ahora es A15**, y no por capricho: mientras el simulador no reproduzca
fracasos conocidos que no se le han programado, sus predicciones no valen para nada fuera
de esta carpeta. Ya acierta tres sin ayuda (ver `INFORME_CEMENTERIO.md` §5).

## Análisis de sensibilidad (A7): qué rompe la estrategia

`guardiana+memoria`, 500 días, neoantígenos calibrados (0,1), 3 semillas.
*Peligrosas* = in situ + invasivas + reinas.

| Suposición que cambiamos | Peligrosas | ¿Aguanta? |
|---|---|---|
| caso base | 2,7 | ✅ |
| escenario BRCA1 (mutación heredada) | 4,3 | ✅ |
| escenario HER2 | 6,3 | ✅ |
| tratar muy tarde (día 400 en vez de 250) | 5,0 | ✅ |
| **plasticidad ×10** (las obreras se vuelven reina 10 veces más fácil) | **312,7** | ❌ **se hunde** |

**Conclusión**: la estrategia es robusta al subtipo de tumor, a la mutación heredada y
al momento del tratamiento. Solo la rompe una cosa: **la plasticidad**. Si cualquier
célula tumoral puede convertirse en célula madre con facilidad, la jerarquía "reina /
obreras" deja de existir, no hay una fuente que reparar y toda la idea se cae.

**Por tanto, la pregunta que decide si esta idea vale algo no es de ingeniería, sino de
biología básica**: ¿con qué facilidad una célula tumoral de mama se desdiferencia en
célula madre *in vivo*? Es un debate abierto en la literatura (Gupta et al., Cell 2011;
Chaffer et al., PNAS 2011; Batlle & Clevers, Nat Med 2017). Ese sería el primer
experimento que pedir a un laboratorio, y es barato comparado con todo lo demás.

## Lo que dicen los datos reales (A2)

1. **Señal exclusiva del tumor** (cuello de botella nº 1): TMB mediana 1,3 mut/Mb y solo un
   2,7 % de tumores con muchas mutaciones. La puerta Y necesita otra señal: marcadores de
   célula madre + secuenciación (como hace la guardiana), HER2 en el 73 % de los HER2+, o
   péptidos de mutaciones recurrentes (PIK3CA H1047R).
2. **No hay un único bug**: solo el 22 % de las pacientes tiene todos sus drivers en nuestros
   8 genes; el 43 % tiene 3 o más. Hace falta **reparación personalizada por paciente**.
3. **Subtipos distintos, bugs distintos**: luminal A → PIK3CA (48 %); basal → TP53 (90 %);
   HER2 → ERBB2 (73 %). La estrategia debe adaptarse a cada subtipo.

## Cómo se reproduce (A8): el contador de divisiones

Una célula tumoral **no tiene sexo**: se divide en dos por mitosis, clones idénticos.
No hay "machos y hembras" que esterilizar. Pero sí hay un cuello de botella universal:

- Toda célula lleva un **contador de divisiones** (los telómeros se acortan en cada copia).
  Al llegar a cero, la célula se jubila: **senescencia**, el límite de Hayflick.
- Para ser cáncer hay que **rebobinar ese contador**, y eso lo hace la **telomerasa (TERT)**.
  Activa en el 73-90 % de los cánceres; en mama la vía alternativa (ALT) es solo ~3 %,
  así que casi todos dependen de ella.

Al añadir TERT al simulador, **el 100 % del tumor acabó teniéndolo sin programarlo**
(emerge por selección, como en los datos reales), y el tumor pasó de ~680 a ~8.300 células:
la inmortalidad es lo que lo vuelve peligroso.

### ¿Funciona la "castración"? (día 250 → 500, neoantígenos 0,1, 3 semillas)

| Estrategia | Peligrosas | Senescentes | Madres sanas afectadas |
|---|---|---|---|
| ninguna | 8.251 | 1 | 0 |
| **solo bloquear la telomerasa** | **7.199** | 99 | 2,7 |
| propagador (virus oncolítico) solo | 8.747 | 1 | 0 |
| guardiana + memoria | 264 | 5 | 0 |
| **guardiana + memoria + bloquear telomerasa** | **42,7** | 8 | 3,3 |

**Bloquear la telomerasa sola casi no sirve** (8.251 → 7.199): las células tienen reserva
de divisiones suficiente y una sola dosis no alcanza a las que nacen después. Esto coincide
con la realidad: los inhibidores de telomerasa (imetelstat y otros) **fracasaron en ensayos
como tratamiento único**, por este mismo motivo.

**Pero como refuerzo es excelente**: combinado con guardiana+memoria mejora 6 veces
(264 → 42,7) tocando solo ~3 células madre sanas.

### Propagación

El propagador se replica de verdad: de 4 unidades iniciales llegó a un pico de 4.539
saltando de célula en célula. Pero el desenlace es **a cara o cruz** (una semilla acabó
en 952 peligrosas y otra en 10.421): es una carrera entre la propagación y el crecimiento
del tumor. Coincide con la variabilidad real de la terapia con virus oncolíticos.

## Fase B: laboratorio (años; necesita un grupo de investigación)

| Pieza | ¿Existe hoy? | Qué falta |
|---|---|---|
| Células que reconocen dos señales a la vez (puerta Y) | Sí: synNotch CAR-T (Roybal 2016), en ensayos | adaptarlas a células madre de mama |
| Encontrar células madre tumorales | Sí: CD44⁺/CD24⁻, ALDH1 (Al-Hajj 2003) | marcador que NO lleven las madres sanas |
| Vehículo que entra invitado (troyano) | Sí: CAR-M (CT-0508), células mesenquimales | CT-0508 dio **cero respuestas** en 14 pacientes, 8 de mama |
| Normalizar / reconducir la vasculatura | Sí, en clínica | saber **a quién** ayuda: en el modelo cura a unos y empeora a otros |
| Frenar el SASP de las senescentes | Sí: inhibidores JAK (ruxolitinib), aprobados | probar que frena la stemness inducida por el tratamiento |
| Memoria de por vida | Sí: linfocitos T residentes, vacunas ARNm personalizadas (en ensayo) | dirigirlas contra lesiones nuevas |
| ~~Reprogramar una reina como guardiana~~ | No | **ya no hace falta**: el brazo se midió y no funciona |
| ~~Editar genes dentro del cuerpo~~ | Parcial: CRISPR in vivo en hígado/sangre | **despriorizado**: ver el acta de defunción en el README |

Nótese lo que ha cambiado: la lista ya no contiene **ninguna pieza que no exista**. Antes
la idea dependía de reprogramar una célula madre tumoral, que no se ha hecho nunca. Ahora
depende de cinco cosas que están en clínica o aprobadas — lo cual no la hace buena, pero
sí la hace *discutible con alguien que trabaje en esto*.

Orden lógico: organoides de mama de pacientes (mini-tumores en placa, semanas) →
ratones con tumores humanos (meses) → seguridad en animales grandes.

## Fase C: ensayos clínicos (años; obligatorio por seguridad)

Fase I (¿es seguro?), fase II (¿funciona?), fase III (¿mejor que lo actual?). No se puede
saltar: es lo que impide que un tratamiento prometedor dañe a las personas.

## Quién podría hacer la fase B en España

CNIO (Madrid), VHIO (Barcelona), IRB Barcelona, CRG, Instituto de Investigación Biomédica
de Salamanca / Centro de Investigación del Cáncer (CIC). El material de este proyecto
(modelo + informe con datos reales) es lo que se lleva a una primera conversación.
