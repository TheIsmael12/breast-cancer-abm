"""Hace importable tests/fingerprint.py desde las pruebas."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
