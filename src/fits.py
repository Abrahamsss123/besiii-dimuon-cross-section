"""Ajustes log-log ponderados (varianza efectiva autoconsistente) y utilidades estadisticas."""
import numpy as np
from scipy import stats

LN10 = np.log(10.0)


def wls_loglog(E, S, dS, dE=0.0008, tol=1e-9, itmax=200):
    """
    Ajuste  log10(sigma) = a*log10(E) + b   (n = -a)  por minimos cuadrados ponderados.
    Incertidumbre de la ordenada Y:  dY^2 = (dS/(S ln10))^2 + (a dE /(E ln10))^2
    - el segundo termino usa la pendiente AJUSTADA (autoconsistente); la semilla es a = 0 (sin teoria).
    """
    E, S, dS = map(np.asarray, (E, S, dS))
    X, Y = np.log10(E), np.log10(S)
    dY0 = dS / (S * LN10)
    a = 0.0
    for it in range(itmax):
        dY2 = dY0**2 + (a * dE / (E * LN10)) ** 2
        w = 1 / dY2
        Sw, Sx, Sy = w.sum(), (w * X).sum(), (w * Y).sum()
        Sxx, Sxy = (w * X * X).sum(), (w * X * Y).sum()
        D = Sw * Sxx - Sx**2
        a_new = (Sw * Sxy - Sx * Sy) / D
        b = (Sxx * Sy - Sx * Sxy) / D
        if abs(a_new - a) < tol:
            a = a_new
            break
        a = a_new
    dY2 = dY0**2 + (a * dE / (E * LN10)) ** 2
    w = 1 / dY2
    Sw, Sx, Sy = w.sum(), (w * X).sum(), (w * Y).sum()
    Sxx, Sxy = (w * X * X).sum(), (w * X * Y).sum()
    D = Sw * Sxx - Sx**2
    a = (Sw * Sxy - Sx * Sy) / D
    b = (Sxx * Sy - Sx * Sxy) / D
    var_a, var_b, cov_ab = Sw / D, Sxx / D, -Sx / D
    chi2 = float((w * (Y - a * X - b) ** 2).sum())
    ndf = len(E) - 2
    X0 = Sx / Sw                                  # pivote (media ponderada de log10 E)
    Y0 = a * X0 + b
    return dict(
        n=-a, sig_n=np.sqrt(var_a), b=b, sig_b=np.sqrt(var_b), cov_ab=cov_ab,
        rho=cov_ab / np.sqrt(var_a * var_b), chi2=chi2, ndf=ndf, p=float(stats.chi2.sf(chi2, ndf)),
        iters=it + 1, X0=X0, E0=10**X0, sigma0=10**Y0, rel_sigma0=LN10 / np.sqrt(Sw),
        Sw=Sw, Sx=Sx, D=D, a=a,
    )


def band_log10(fit, Xgrid):
    """Desviacion estandar (1 sigma) de la recta ajustada en log10(sigma), con la covarianza completa."""
    v = fit["sig_b"] ** 2 + Xgrid**2 * fit["sig_n"] ** 2 + 2 * Xgrid * fit["cov_ab"]
    return np.sqrt(v)


def chi2_norm_free(S, dS, model):
    """chi2 de  S = N*model  con N libre (solucion analitica)."""
    w = 1 / dS**2
    N = (w * S * model).sum() / (w * model**2).sum()
    sigN = 1 / np.sqrt((w * model**2).sum())
    chi2 = float((w * (S - N * model) ** 2).sum())
    ndf = len(S) - 1
    return dict(N=N, sig_N=sigN, chi2=chi2, ndf=ndf, p=float(stats.chi2.sf(chi2, ndf)))


def chi2_fixed(S, dS, model):
    chi2 = float((((S - model) / dS) ** 2).sum())
    ndf = len(S)
    return dict(chi2=chi2, ndf=ndf, p=float(stats.chi2.sf(chi2, ndf)))
