# Cómo contribuir

Gracias por querer ayudar. Este es un proyecto de **investigación abierta**: lo que más
falta no es código, sino **criterio biológico**. Si trabajas en oncología, biología
celular, inmunología o modelado, tus correcciones valen oro.

## Qué se busca

| Tipo | Ejemplo | Cómo empezar |
|---|---|---|
| **Corregir un supuesto** | «Las células hipóxicas viven más de lo que dice el modelo» | Abre un *issue* «Supuesto biológico» con la fuente |
| **Proponer un experimento** | «¿Y si el arsenal se da antes de la ablación?» | *Issue* «Nuevo experimento», con el criterio de éxito |
| **Añadir un test de validación** | Un dato clínico que el modelo deba reproducir | PR con un `Tn` nuevo en `validacion.py` |
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
   [bitácora](docs/BITACORA.md) está llena de ellos y es lo más valioso del proyecto.
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
pip install -r requirements.txt pyflakes pytest
git checkout -b mi-cambio

# antes de abrir el PR, lo mismo que hace la CI:
python -m pyflakes *.py tests/*.py
python -m pytest -q tests
```

Si tu cambio es un **avance** (un resultado nuevo), añade también:

- `docs/avances/AAAA-MM-DD-titulo.md`: pregunta, fuentes, qué cambiaste, resultados,
  errores encontrados y limitaciones (usa [el de 2026-11-05](docs/avances/2026-11-05-tratar-sin-operar.md)
  como plantilla), y una fila en [`docs/avances/README.md`](docs/avances/README.md).
- `resultados/AAAA-MM-DD/`: la salida bruta de cada ejecución.
- Capturas en `docs/capturas/` si ayudan (`python capturas.py`).

## Cuidado con el ordenador

Los bancos de pruebas usan `multiprocessing`. Pasa siempre `--procesos` con menos núcleos
de los que tienes; con todos, el equipo se vuelve inusable durante la ejecución.

## Lo que este proyecto NO es

No es una herramienta clínica. Ningún resultado de aquí justifica una decisión médica,
y ningún PR debe presentarlo así.
