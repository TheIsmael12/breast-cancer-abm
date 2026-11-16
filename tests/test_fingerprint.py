"""El modelo sigue dando exactamente los mismos números (ver fingerprint.py)."""

from __future__ import annotations

import json

from fingerprint import FILE, compute


def test_model_fingerprint_unchanged():
    expected = json.loads(FILE.read_text(encoding="utf-8"))
    got = compute()
    for key in expected:
        assert got[key] == expected[key], f"la huella '{key}' ha cambiado"
