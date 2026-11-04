# Dependencia energética de σ(e⁺e⁻ → μ⁺μ⁻) con datos de BESIII (3.8–4.6 GeV)

Código y datos del manuscrito *"Dependencia energética de la sección eficaz dimuónica con datos de BESIII en 3.8–4.6 GeV"*
(Abraham H. Gutierrez). Todo número, tabla y figura del manuscrito se genera con los scripts de este repositorio.

*English summary: code and data to reproduce the power-law (σ ∝ E⁻ⁿ) analysis of the BESIII e⁺e⁻→μ⁺μ⁻ cross section between 3.8 and
4.6 GeV: vector-graphics extraction of 133 data points from the published figure, log–log fit, BESIII resonance model, Monte Carlo.*

## Qué hay aquí

| Archivo | Qué hace |
|---|---|
| `data/besiii_bcdr_133.csv` | Los **133 puntos** de la Fig. 3 de BESIII (E_cm, σ_BCDR, error estadístico, tipo de luminosidad). Es el resultado del paso 1. |
| `src/extract_fig3.py` | **Paso 1 – Extracción de los datos** desde el PDF vectorial (ver `docs/EXTRACCION.md`). |
| `src/figcal.py` | Calibración de los ejes de la figura (posiciones de las marcas → GeV y nb). Lo usan los pasos 1 y 2. |
| `src/validate_extraction.py` | **Paso 2 – Validación**: compara las curvas que BESIII dibuja en la figura con nuestro modelo (acuerdo ≈ 0.02 %); genera `docs/extraccion_overlay.png`. |
| `src/besiii_model.py` | Modelo de BESIII (Ec. 3 de su artículo): continuo + 9 resonancias + S(4220); QED a orden árbol; factor de masa del muón; polarización del vacío leptónica; FSR. |
| `src/fits.py` | Ajuste log–log ponderado con incertidumbres **no circulares**, energía pivote, χ², normalización libre. |
| `src/common.py` | Rutas, constantes del análisis (errores de extracción, sistemática de 2.91 %, semilla) y lectura del CSV. |
| `src/analysis.py` | **Paso 3 – Análisis**: ajuste de n, comparación con el modelo, escalera de modelos, contribución de cada resonancia, ISR/FSR/masa/polarización del vacío. Escribe `results/results.json`. |
| `src/montecarlo.py` | **Paso 4 – Monte Carlo** (hipótesis nula n = 2 exacto): compara estimar α con el intercepto vs. con el pivote. Escribe `results/mc_toys.npz`. |
| `src/figures.py` | **Paso 5** – Figuras 2–5 del manuscrito (`figs/`). |
| `src/make_tex.py` | **Paso 6** – Genera `paper/numeros.tex` (todas las cifras del texto) y las tablas LaTeX. |
| `tests/check_results.py` | Comprueba que los resultados clave coinciden con los del manuscrito. |
| `paper/` | Fuente LaTeX del manuscrito (`Manuscrito_RIF_AbrahamG_corregido.tex`) y los archivos que genera `make_tex.py`. |
| `run_all.sh` | Ejecuta todo en orden. |

## Cómo reproducir

```bash
pip install -r requirements.txt          # numpy, scipy, matplotlib, pymupdf
# además: poppler-utils (pdftocairo) y, para compilar el manuscrito, una instalación de LaTeX (pdflatex)

# (a) Solo el análisis, con el CSV incluido (no necesita el PDF de BESIII):
./run_all.sh

# (b) Desde cero, incluyendo la extracción de los datos (necesita el PDF publicado de BESIII):
./run_all.sh /ruta/a/PhysRevD_102_112009.pdf
```

El PDF de BESIII (M. Ablikim et al., *Phys. Rev. D* **102**, 112009 (2020), arXiv:2007.12872) **no se redistribuye** aquí:
descárguelo de la revista. Se usó la versión publicada, cuyo SHA-256 es
`2cf8be0918ba694f6ee428aaf9c3d917f7886c54a4341f48253cb090c1415ec4` (el script lo imprime). Con otra versión del PDF el script se
detiene si no encuentra exactamente 133 marcadores.

Resultado esperado, con n como exponente efectivo (se verifica con `tests/check_results.py`): n = 2.036 ± 0.010, χ²/ndf = 222.5/131, E₀ = 4.206 GeV;
el modelo completo de BESIII da n = 2.026 y un χ² solo descriptivo de 119.6 (sus parámetros los fijó BESIII con estos mismos datos); el Monte Carlo muestra que el intercepto amplifica ≈ 31× el error de α.

## El Monte Carlo en una frase

Se simulan 40 000 conjuntos de 133 puntos que siguen **exactamente** σ ∝ E⁻² (con las energías y errores reales) y se estima α de dos
maneras; el intercepto fluctúa 0.73 % y el pivote 0.024 %, aunque no haya ningún efecto físico: la causa es la correlación
ρ(n, log k) ≈ −0.9995 entre pendiente e intercepto. Semilla fija en `src/common.py` (`SEED = 20261103`).

## Limitaciones (también en el manuscrito)

* Los puntos BCDR se obtienen con el propio modelo de BESIII: el acuerdo es una prueba de consistencia, no una validación independiente.
* En 5 puntos de alta luminosidad la barra de error queda tras el marcador: se usa su radio como cota superior.
* No se calcula la polarización del vacío hadrónica ni la FSR a órdenes superiores.

## Licencia y créditos

Código: MIT (`LICENSE`). Los datos extraídos (`data/`) provienen de la figura publicada por la colaboración BESIII; cítese
Phys. Rev. D 102, 112009 (2020) al usarlos. Cite este trabajo con `CITATION.cff` (agregue el DOI de Zenodo al publicar).
