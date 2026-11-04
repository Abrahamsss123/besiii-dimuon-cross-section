# Cómo se extrajeron los 133 puntos (paso a paso)

**Idea central.** La Fig. 3 del artículo de BESIII (PRD 102, 112009, página 8 del PDF publicado) está incrustada como **gráfico
vectorial**, no como imagen. Eso significa que el PDF guarda, con precisión de décimas de milésima de punto, la posición exacta de cada
marcador y de cada barra de error. No hay que "adivinar" píxeles ni hacer clic sobre la figura: se leen esos números y se convierten
a unidades físicas. El script es `src/extract_fig3.py`; aquí se explica qué hace en cada etapa.

## 1. Calibración de los ejes
En el PDF las coordenadas están en puntos tipográficos (1 pt = 1/72 in). Se leyeron del PDF las posiciones de las marcas mayores de los
ejes (9 en x: 3.8, 3.9 … 4.6 GeV; 11 en y: 6.2, 6.0 … 4.2 nb) y se ajustó una recta posición → valor para cada eje
(`src/figcal.py`). Resultado: 174.767 pt/GeV y 46.722 pt/nb, con residuos máximos de 2.2×10⁻⁴ GeV y 8.6×10⁻⁴ nb.
Esos residuos fijan la incertidumbre aleatoria de la extracción (0.001 nb por punto).

## 2. Los marcadores (la posición de cada punto)
Con PyMuPDF (`page.get_drawings()`) se listan todos los trazos del panel izquierdo. Cada marcador es un dibujo de 4 segmentos con menos de
2 pt de lado. El color distingue la luminosidad: **relleno rojo = alta** luminosidad (29 puntos), **contorno negro hueco = baja** (104).
El centro del marcador da (E_cm, σ) tras aplicar la calibración. El script comprueba que salen exactamente 133 marcadores (29 + 104),
el mismo número de energías que reporta BESIII; si no, se detiene.

## 3. Las barras de error
PyMuPDF no entrega las barras como trazos separados. Por eso la página se convierte a SVG con `pdftocairo -svg` y se leen los segmentos
verticales de grosor 0.484 o 0.726 pt, aplicando la matriz de transformación (`matrix(a,b,c,d,e,f)`) de cada trazo. Cada marcador se
asocia con los segmentos a menos de 0.06 pt de su abscisa: uno que termina en el centro por arriba y otro que empieza en el centro por
abajo. La barra es el promedio de las dos semi-longitudes, convertido a nb con la escala vertical. Se recuperaron 239 segmentos y la barra
de 128 de los 133 puntos.

**Barras ocultas.** En 5 puntos de alta luminosidad la barra es más pequeña que el marcador y no se ve. El CSV deja el campo vacío y el
análisis usa como cota superior el radio del marcador (0.645 pt = 0.0138 nb, `common.py`).

## 4. Qué significa "σ_BCDR"
Los puntos de esta figura no son la sección eficaz observada, sino la **BCDR** (Born-continuum-dressed-resonance): la observada
dividida por el factor de corrección de radiación del estado inicial (Fig. 2 = observada; Fig. 3 = BCDR). Por eso la ISR ya está
removida en estos datos. El manuscrito explica la identificación (Sec. 3.2).

## 5. Validación (src/validate_extraction.py)
La misma figura contiene las curvas que BESIII ajustó (azul: ajuste completo; cian: continuo f_cnt/s). Se leen sus vértices, se calibran
igual y se comparan con `src/besiii_model.py`:

| Comparación | rms | desfase medio |
|---|---|---|
| curva cian vs. f_cnt/s | 0.012 % | −0.14 % |
| curva azul vs. modelo completo (9 resonancias + S) | 0.020 % | −0.13 % |
| curva azul vs. modelo sin las 6 resonancias bajas | 0.137 % (máx. 0.96 %) | +0.50 % |

El acuerdo valida a la vez la calibración y la implementación del modelo. El desfase constante de 0.13 % se adopta como incertidumbre
sistemática de la extracción (`SYS_EXTR`). `docs/extraccion_overlay.png` superpone los puntos extraídos sobre la figura original para
inspección visual.

## 6. Limitaciones del método
* Solo recupera lo que está dibujado: el error estadístico de la barra, no la matriz de covarianza ni la sistemática de 2.91 %
  (que se trata como normalización común).
* Depende de la versión del PDF (página y grosores de línea); con otra versión el script falla de forma explícita.
* No sustituye a una tabla numérica oficial de BESIII; si la colaboración la publica, debería usarse en su lugar.

## Reproducir el paso
```bash
python3 src/extract_fig3.py --pdf PhysRevD_102_112009.pdf --out data/besiii_bcdr_133.csv
python3 src/validate_extraction.py --pdf PhysRevD_102_112009.pdf
```
Requiere `pymupdf`, `numpy`, `matplotlib` y `pdftocairo` (poppler-utils). El CSV generado es idéntico byte a byte al incluido.
