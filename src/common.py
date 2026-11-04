"""Rutas, constantes del análisis y lectura de los datos. Lo importan analysis.py, montecarlo.py y figures.py."""
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA_CSV = ROOT / "data" / "besiii_bcdr_133.csv"
RESULTS = ROOT / "results"
FIGS = ROOT / "figs"
LATEX = RESULTS / "latex"          # macros y tablas LaTeX generadas (src/make_tex.py)

# --- incertidumbres de los datos --------------------------------------------------------------------------
PT_PER_NB = 46.722                       # escala vertical de la figura de BESIII (pt por nb), de figcal.py
MARKER_R_PT = 0.645                      # radio del marcador en pt
MARKER_R_NB = MARKER_R_PT / PT_PER_NB    # cota superior de la barra cuando queda oculta tras el marcador (nb)
EXTR_NB = 0.0010                         # incertidumbre aleatoria de la extracción (nb): residuo de calibración de ejes
SYS = 0.0291                             # sistemática total de BESIII (2.91 %), totalmente correlacionada
SYS_EXTR = 0.0013                        # desfase constante (0.13 %) hallado en validate_extraction.py
SEED = 20261103                          # semilla del Monte Carlo


def load_data(csv=DATA_CSV):
    """Devuelve dict con E, S (nb), err_stat, dS (error total por punto), cls ('alta'/'baja') y hidden (barra oculta)."""
    d = np.genfromtxt(csv, delimiter=",", names=True, dtype=None, encoding="utf8")
    E = d["Ecm_GeV"].astype(float)
    S = d["sigma_BCDR_nb"].astype(float)
    err = d["err_stat_nb"].astype(float)          # campo vacío -> nan
    cls = d["lumi_class"].astype(str)
    hidden = np.isnan(err)
    err_stat = np.where(hidden, MARKER_R_NB, err)
    dS = np.sqrt(err_stat**2 + EXTR_NB**2)
    return dict(E=E, S=S, err_stat=err_stat, dS=dS, cls=cls, hidden=hidden)
