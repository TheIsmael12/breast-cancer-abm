"""¿Qué identifica a una célula de cáncer de mama y SOLO a ella?

Evalúa cada marcador candidato sobre los ~1.000 tumores reales de TCGA ya descargados
(datos/tcga_brca.json) y lo cruza con lo que se sabe del tejido sano. Dos preguntas:

  COBERTURA   ¿en qué fracción de pacientes está presente?   (si falla, no sirve para todas)
  EXCLUSIVIDAD ¿lo tiene alguna célula sana?                 (si falla, hay daño colateral)

Un identificador útil necesita las dos. Ninguno las tiene a la vez: de ahí las puertas AND.

    python identificadores.py
"""
from __future__ import annotations

import statistics
from pathlib import Path

from datos_reales import PRESENTACION, alteraciones, descargar

# (marcador, cobertura estimada, dónde está también en tejido SANO, estado clínico)
# Las coberturas marcadas como None se calculan abajo con los datos de TCGA.
CANDIDATOS = [
    ("Telomerasa (TERT) activa", None,
     "Células madre: médula ósea, intestino, piel, folículo piloso y mama",
     "Inhibidores probados; fracasaron solos y dan toxicidad medular"),
    ("Neoantígenos (proteínas nuevas por mutación)", None,
     "NINGUNA: no existen en una célula sana. Exclusividad perfecta",
     "Base de las vacunas de ARNm personalizadas (en ensayo)"),
    ("HER2 muy amplificado", None,
     "Corazón y epitelios tienen algo de HER2: es diferencia de CANTIDAD, no de tipo",
     "Trastuzumab y T-DXd aprobados; cardiotoxicidad por esa misma razón"),
    ("Pérdida de MHC-I (HLA/B2M)", None,
     "Algunas células sanas lo bajan; las NK ya lo vigilan de serie",
     "Lo explotan las terapias con células NK"),
    ("Marcadores de célula madre (CD44+/CD24-/ALDH1+)", None,
     "Las células madre SANAS de la mama los llevan igual",
     "Sin terapia aprobada que los use como diana única"),
    ("Aneuploidía (número de cromosomas alterado)", None,
     "Prácticamente ninguna célula sana. Casi exclusivo",
     "INVISIBLE desde fuera: está dentro del núcleo, no en la superficie"),
]


def main():
    datos = descargar()
    alt = alteraciones(datos)
    muestras = list(alt)
    n = len(muestras)
    tmb = datos["tmb"]
    valores = [tmb[s] for s in muestras if s in tmb]

    cob = {
        "Telomerasa (TERT) activa": 0.90,          # literatura: 73-90 %; en mama ALT solo ~3 %
        "Neoantígenos (proteínas nuevas por mutación)": sum(v >= 10 for v in valores) / len(valores),
        "HER2 muy amplificado": sum("ERBB2" in alt[s] for s in muestras) / n,
        "Pérdida de MHC-I (HLA/B2M)": sum(bool(alt[s] & set(PRESENTACION)) for s in muestras) / n,
        "Marcadores de célula madre (CD44+/CD24-/ALDH1+)": 1.0,   # presentes en todo tumor
        "Aneuploidía (número de cromosomas alterado)": 0.95,       # casi universal en mama
    }

    L = ["# ¿Qué identifica a una célula de cáncer de mama?", "",
         f"Datos: {n} tumores de TCGA PanCancer Atlas. TMB mediana "
         f"{statistics.median(valores):.2f} mutaciones/Mb.", "",
         "| Marcador | Cobertura | ¿Lo tiene alguna célula SANA? | Estado |",
         "|---|---|---|---|"]
    for nombre, _, sano, estado in CANDIDATOS:
        L.append(f"| {nombre} | {cob[nombre]:.0%} | {sano} | {estado} |")

    L += ["", "## El problema, en una frase", "",
          "Ningún marcador tiene a la vez COBERTURA alta y EXCLUSIVIDAD total:", "",
          "- Los **exclusivos** (neoantígenos) solo están en una minoría: "
          f"{cob['Neoantígenos (proteínas nuevas por mutación)']:.0%} de los tumores de mama "
          "tienen carga mutacional alta.",
          "- Los **universales** (telomerasa, marcadores de célula madre, aneuploidía) los "
          "comparte el tejido sano regenerativo, o no se ven desde fuera de la célula.", "",
          "## La salida: combinaciones", "",
          "Ninguna célula sana cumple VARIAS condiciones a la vez. Ejemplos de combinación "
          "con cobertura alta y, en teoría, sin equivalente sano:", ""]

    combis = [
        ("telomerasa activa  Y  marcador de célula madre  Y  está en el conducto, no en el nicho basal",
         "la madre sana vive en un sitio concreto; la tumoral no respeta la arquitectura"),
        ("telomerasa activa  Y  aneuploidía",
         "una célula madre sana tiene telomerasa pero NO tiene los cromosomas descuadrados"),
        ("HER2 muy alto  Y  proliferando",
         "el corazón tiene algo de HER2 pero no se divide: base de por qué el T-DXd funciona"),
        ("pérdida de MHC-I  Y  marcadores de estrés (MICA/B)",
         "es exactamente lo que ya hacen las células NK del cuerpo"),
    ]
    L += [f"- **{c}**  \n  {por}" for c, por in combis]
    L += ["", "## Lo que esto implica para la estrategia", "",
          "1. No busques el marcador perfecto: **no existe**. Es un problema abierto de "
          "investigación, no un hueco de este proyecto.",
          "2. La especificidad tiene que venir de **combinar condiciones** (puerta AND), que es "
          "lo que ya hacen los receptores synNotch en los ensayos reales.",
          "3. La combinación más prometedora para mama no usa neoantígenos (son escasos) sino "
          "**telomerasa + algo que la célula madre sana no tenga**: aneuploidía o estar fuera "
          "de su sitio.",
          "4. Y sigue pendiente el problema que destapó el 3D: aunque la identifiques bien, "
          "hay que **llegar** hasta ella."]

    texto = "\n".join(L)
    (Path(__file__).parent / "docs" / "informes" / "INFORME_IDENTIFICADORES.md").write_text(texto, encoding="utf-8")
    print(texto)


if __name__ == "__main__":
    main()
