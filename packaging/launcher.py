"""Punto de entrada del ejecutable: abre la aplicación en su propia ventana."""

import multiprocessing
import sys

if __name__ == "__main__":
    multiprocessing.freeze_support()  # pywebview y los estudios usan procesos hijos
    from breast_cancer_abm.app.main import main

    sys.exit(main(sys.argv[1:]))
