"""
Modelo de seccion eficaz 'dressed' / BCDR de BESIII, Ec. (3) de PRD 102, 112009 (2020), con x = 0:

    sigma_D(s) = | A_cnt + sum_{k=1..9} e^{i phi_k} A_Rk + e^{i phi_S} A_S |^2

    A_cnt = sqrt(f_cnt / s)
    A_R   = sqrt(12 pi Gee Gmm) / [ (s - M^2) + i Gtot M ]      (unidades naturales; x C -> nb)

Unidades: GeV, nb.   C = (hbar c)^2 = 389379 nb GeV^2.
Parametros de las 6 resonancias bajas fijos a PDG con fase 0 (como en BESIII); lepton universality Gee = Gmm.
"""
import numpy as np

C = 389379.0            # nb GeV^2
ALPHA0 = 1 / 137.035999
F_CNT = 88.51           # nb GeV^2 (Solucion I, BESIII Tabla/texto)

# (M, Gtot, Gee=Gmm) en GeV.  PDG (listados 2017-2020; BESIII usa PDG 2018).
LOW = {
    "rho(770)":  (0.77526, 0.1491,    7.04e-6),
    "omega(782)": (0.78265, 0.00849,  0.60e-6),
    "phi(1020)": (1.019461, 0.004249, 2.973e-4 * 4.249e-3),   # Gee = B(ee) * Gtot
    "J/psi":     (3.096900, 92.9e-6,  5.53e-6),
    "psi(2S)":   (3.686097, 296e-6,   2.34e-6),
    "psi(3770)": (3.77313, 0.0272,    0.262e-6),
}
# psi(4040), psi(4160), psi(4415): M y Gtot fijos a PDG (BESIII); Gmm y fase libres (Tabla I)
HEAVY = {"psi(4040)": (4.039, 0.080), "psi(4160)": (4.191, 0.070), "psi(4415)": (4.421, 0.062)}

# Tabla I de BESIII: Gmm [keV], fase [grados]; S(4220): M, Gtot [MeV], Gmm [keV], fase
_T = {
 "I":    dict(g=(0.73, 2.45, 1.25), p=(283, 153, 334), S=(4216.7, 47.2, 1.53, 20)),
 "II":   dict(g=(0.62, 2.36, 1.26), p=(285, 138, 332), S=(4213.6, 39.9, 1.28, 0)),
 "III":  dict(g=(0.58, 2.28, 1.27), p=(286, 136, 332), S=(4213.7, 38.5, 1.20, 359)),
 "IV":   dict(g=(0.71, 2.41, 1.24), p=(283, 150, 333), S=(4216.2, 45.5, 1.46, 17)),
 "V":    dict(g=(0.74, 2.28, 1.24), p=(282, 157, 335), S=(4219.4, 49.6, 1.50, 31)),
 "VI":   dict(g=(0.58, 2.22, 1.27), p=(287, 132, 331), S=(4212.8, 36.4, 1.12, 352)),
 "VII":  dict(g=(0.66, 2.08, 1.24), p=(291, 143, 332), S=(4216.1, 37.8, 1.09, 10)),
 "VIII": dict(g=(0.80, 2.31, 1.25), p=(284, 154, 332), S=(4217.3, 45.5, 1.40, 22)),
}
SOLUTIONS = _T


def _A(s, M, Gt, Gee, Gmm):
    return np.sqrt(C * 12 * np.pi * Gee * Gmm) / ((s - M**2) + 1j * M * Gt)


def amplitudes(E, sol="I", f_cnt=F_CNT):
    """Devuelve dict nombre -> amplitud compleja (con fase incluida), en sqrt(nb)."""
    E = np.asarray(E, float)
    s = E**2
    T = _T[sol]
    out = {"continuo": np.sqrt(f_cnt / s) + 0j}
    for k, (M, Gt, Ge) in LOW.items():
        out[k] = _A(s, M, Gt, Ge, Ge)
    for i, (k, (M, Gt)) in enumerate(HEAVY.items()):
        g = T["g"][i] * 1e-6
        out[k] = np.exp(1j * np.radians(T["p"][i])) * _A(s, M, Gt, g, g)
    M4, G4, g4, p4 = T["S"]
    g4 *= 1e-6
    out["S(4220)"] = np.exp(1j * np.radians(p4)) * _A(s, M4 * 1e-3, G4 * 1e-3, g4, g4)
    return out


def sigma_D(E, sol="I", f_cnt=F_CNT, drop=()):
    amp = amplitudes(E, sol, f_cnt)
    tot = sum(v for k, v in amp.items() if k not in drop)
    return np.abs(tot) ** 2


def interference_percent(E, sol="I", f_cnt=F_CNT):
    """Contribucion de cada resonancia al sigma total (interferencia lineal con el continuo, en % de f_cnt/s)."""
    amp = amplitudes(E, sol, f_cnt)
    cnt = amp["continuo"]
    return {k: 200 * np.real(np.conj(cnt) * v) / np.abs(cnt) ** 2 for k, v in amp.items() if k != "continuo"}


def sigma_tree(E, alpha=ALPHA0):
    """QED a orden arbol, sin masa:  sigma = 4 pi alpha^2 / (3 s)   [nb]."""
    return 4 * np.pi * alpha**2 / 3 * C / np.asarray(E, float) ** 2


def beta_mu_factor(E, m_mu=0.1056583755):
    """sigma / sigma_(m=0) = beta (3 - beta^2)/2  (seccion eficaz total, mu con masa)."""
    b = np.sqrt(1 - 4 * m_mu**2 / np.asarray(E, float) ** 2)
    return b * (3 - b**2) / 2


# ----- polarizacion del vacio leptonica (1 lazo, exacta) -----
_MASS = {"e": 0.51099895e-3, "mu": 0.1056583755, "tau": 1.77686}


def delta_alpha_lep(E):
    s = np.asarray(E, float) ** 2
    tot = 0.0
    for m in _MASS.values():
        b = np.sqrt(1 - 4 * m**2 / s)
        tot = tot + ALPHA0 / (3 * np.pi) * (-5 / 3 - 4 * m**2 / s + (1 + 2 * m**2 / s) * b * np.log((1 + b) / (1 - b)))
    return tot


def alpha_running(E):
    return ALPHA0 / (1 - delta_alpha_lep(E))
