#!/usr/bin/env bash
# Reproduce todo el análisis. Uso:
#   ./run_all.sh [ruta/a/PhysRevD_102_112009.pdf]
# Sin argumento se omiten los pasos 1-2 (extracción y validación) y se usa data/besiii_bcdr_133.csv ya incluido.
set -euo pipefail
cd "$(dirname "$0")"
PDF="${1:-}"
if [ -n "$PDF" ]; then
  echo "== Paso 1: extracción de los 133 puntos desde el PDF vectorial =="
  python3 src/extract_fig3.py --pdf "$PDF" --out data/besiii_bcdr_133.csv
  echo "== Paso 2: validación contra las curvas de BESIII =="
  python3 src/validate_extraction.py --pdf "$PDF"
else
  echo "(sin PDF: se usa el CSV incluido; pasos 1-2 omitidos)"
fi
echo "== Paso 3: análisis =="        ; python3 src/analysis.py
echo "== Paso 4: Monte Carlo =="     ; python3 src/montecarlo.py
echo "== Paso 5: figuras =="         ; python3 src/figures.py
echo "== Paso 6: macros y tablas =="  ; python3 src/make_tex.py
echo "== Comprobación de resultados =="; python3 tests/check_results.py
cp figs/*.pdf paper/
if command -v pdflatex >/dev/null; then
  echo "== Paso 7: compilar el manuscrito =="
  (cd paper && for i in 1 2 3; do pdflatex -interaction=nonstopmode Manuscrito_RIF_AbrahamG_corregido.tex >/dev/null; done)
  echo "PDF: paper/Manuscrito_RIF_AbrahamG_corregido.pdf"
fi
