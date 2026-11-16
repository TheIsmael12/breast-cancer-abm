# Cómo contribuir

Gracias por querer ayudar. Este es un proyecto de **investigación abierta**: lo que más
falta no es código, sino **criterio biológico**. Si trabajas en oncología, biología
celular, inmunología o modelado, tus correcciones valen oro.

## Cómo se trabaja aquí

El proyecto lo desarrolla una IA (Claude, de Anthropic) dirigida por una persona que no
es del sector sanitario. Eso significa dos cosas: el trabajo es rápido y está muy
documentado, pero **le falta la mirada de quien conoce la biología por dentro**. Las
revisiones de los PR también se hacen con ayuda de Claude, y cualquier corrección con
fuente se toma en serio aunque contradiga resultados ya publicados.

## Qué se busca

| Tipo | Ejemplo | Cómo empezar |
|---|---|---|
| **Corregir un supuesto** | «Las células hipóxicas viven más de lo que dice el modelo» | Abre un *issue* «Supuesto biológico» con la fuente |
| **Proponer un experimento** | «¿Y si el arsenal se da antes de la ablación?» | *Issue* «Nuevo experimento», con el criterio de éxito |
| **Añadir un test de validación** | Un dato clínico que el modelo deba reproducir | PR con un `Tn` nuevo en `studies/validation.py` |
| **Mejorar el código** | Rendimiento, 3D, visor | PR directo |
| **Datos** | Series de pacientes, parámetros medidos | *Issue* «Dato o fuente» |

## Las normas de la casa

Son las que han sostenido el proyecto desde el principio. Los PR que no las cumplan se
devolverán con cariño.

1. **El criterio de éxito se escribe ANTES de ejecutar.** Si cambias el criterio después
   de ver el resultado, no es un experimento, es una justificación.
2. **Varias semillas, siempre.** Como mínimo 12 para un resultado y 24 para una
   conclusión. Una semilla no demuestra nada.
3. **Un resultado que contradice la expectativa se publica igual.** La
   [bitácora](docs/LOGBOOK.md) está llena de ellos y es lo más valioso del proyecto.
4. **Cada supuesto inventado se marca como `SUPUESTO`** en el código, con por qué se
   eligió ese valor. Cada valor medido lleva su fuente.
5. **Los agentes solo perciben la superficie celular.** Para ver el genoma tienen que
   *Analizar*. No se le da a una terapia información que no podría tener.
6. **No se cambia un valor por defecto sin medir antes y después** con las mismas
   semillas. Si cambia un resultado ya publicado, se dice en el PR.

## Flujo de trabajo

```bash
git clone https://github.com/TheIsmael12/breast-cancer-abm
cd breast-cancer-abm
python -m venv .venv && .venv/Scripts/activate      # en Linux/macOS: source .venv/bin/activate
pip install -e ".[app,dev]"
git checkout -b mi-cambio

# antes de abrir el PR, lo mismo que hace la CI:
ruff check src tests
ruff format --check src tests
pytest
```

`tests/test_fingerprint.py` comprueba que el modelo da **exactamente** los mismos números
que antes. Si tu cambio la rompe a propósito (un cambio de modelo, no de estilo), regenera
la huella con `python tests/fingerprint.py` y explícalo en el PR y en la bitácora.

### Dónde va cada cosa

| Carpeta | Qué contiene |
|---|---|
| `src/breast_cancer_abm/model/` | el motor: biología, agentes, terapias, tejido 2D y 3D |
| `src/breast_cancer_abm/studies/` | un banco de pruebas por pregunta (`bcabm study --list`) |
| `src/breast_cancer_abm/viz/` | visor de matplotlib y figuras de `docs/figures/` |
| `src/breast_cancer_abm/app/` | la aplicación (`bcabm app`) |
| `tests/` | prueba de humo y huella del modelo |
| `data/`, `results/`, `docs/` | datos de entrada, salidas brutas y documentación |

El código (nombres de módulos, funciones, variables, opciones) va en **inglés**; los
comentarios, la documentación y los textos que ve el usuario, en **español**.

Si tu cambio es un **avance** (un resultado nuevo), añade también:

- `docs/progress/AAAA-MM-DD-titulo.md`: pregunta, fuentes, qué cambiaste, resultados,
  errores encontrados y limitaciones (usa [el de 2026-11-05](docs/progress/2026-11-05-treatment-without-surgery.md)
  como plantilla), y una fila en [`docs/progress/README.md`](docs/progress/README.md).
- `results/AAAA-MM-DD/`: la salida bruta de cada ejecución.
- Capturas en `docs/figures/` si ayudan (`bcabm figures`).

## Cuidado con el ordenador

Los bancos de pruebas usan `multiprocessing`. Pasa siempre `--processes` con menos núcleos
de los que tienes; con todos, el equipo se vuelve inusable durante la ejecución.

## Lo que este proyecto NO es

No es una herramienta clínica. Ningún resultado de aquí justifica una decisión médica,
y ningún PR debe presentarlo así.
