"""Estudios del proyecto: cada uno es un banco de simulaciones con su criterio escrito antes.

``STUDIES`` es el registro que usan el comando ``bcabm study`` y la aplicación.
"""

from __future__ import annotations

import importlib
import inspect
from dataclasses import dataclass
from functools import cached_property


@dataclass(frozen=True)
class Study:
    name: str  # nombre para la línea de comandos
    module: str  # módulo con main(argv)
    topic: str  # agrupación en la aplicación

    @cached_property
    def doc(self) -> str:
        """Docstring del módulo (importarlo es barato: solo trae el modelo, no matplotlib)."""
        return inspect.cleandoc(importlib.import_module(self.module).__doc__ or "")

    @property
    def title(self) -> str:
        return self.doc.splitlines()[0] if self.doc else self.name

    def run(self, argv: list[str] | None = None):
        return importlib.import_module(self.module).main(argv)


_S = "breast_cancer_abm.studies"
STUDIES: dict[str, Study] = {
    s.name: s
    for s in (
        Study("no-surgery", f"{_S}.no_surgery", "Tratar sin operar"),
        Study("no-surgery-3d", f"{_S}.no_surgery_3d", "Tratar sin operar"),
        Study("residual-test", f"{_S}.residual_test", "Tratar sin operar"),
        Study("experiment", f"{_S}.experiment", "Terapias"),
        Study("chain", f"{_S}.chain", "Terapias"),
        Study("propagation", f"{_S}.propagation", "Terapias"),
        Study("reversion", f"{_S}.reversion", "Terapias"),
        Study("residual", f"{_S}.residual", "Terapias"),
        Study("gates", f"{_S}.gates", "Identificar el tumor"),
        Study("identifiers", f"{_S}.identifiers", "Identificar el tumor"),
        Study("real-data", f"{_S}.real_data", "Identificar el tumor"),
        Study("barrier", f"{_S}.barrier", "Acceso al tumor"),
        Study("mesh", f"{_S}.mesh", "Acceso al tumor"),
        Study("tissue3d", "breast_cancer_abm.model.tissue3d", "Acceso al tumor"),
        Study("dissemination", "breast_cancer_abm.model.dissemination", "Diseminación"),
        Study("validation", f"{_S}.validation", "Validación"),
    )
}
