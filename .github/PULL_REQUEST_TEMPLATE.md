## Qué cambia

<!-- Una o dos frases. -->

## Tipo

- [ ] Corrección de un supuesto biológico
- [ ] Experimento / avance nuevo
- [ ] Test de validación
- [ ] Código (rendimiento, visor, 3D…)
- [ ] Documentación

## Si es un experimento o cambia la biología

- **Criterio de éxito, escrito antes de ejecutar:**
- **Semillas usadas:** <!-- mínimo 12 -->
- **Antes / después con las mismas semillas:**
- **Fuentes** (artículo, DOI o enlace):
- [ ] He marcado como `SUPUESTO` cada valor que no sale de una fuente
- [ ] Si cambia un resultado ya publicado, lo digo aquí y en `docs/`

## Comprobaciones

- [ ] `python -m pyflakes *.py tests/*.py` sin avisos
- [ ] `python -m pytest -q tests` pasa
- [ ] Si es un avance: documento en `docs/avances/`, salidas en `resultados/`
