"""
PASO 4. Monte Carlo bajo la hipótesis nula: ¿por qué la normalización se mide en el pivote y no en el intercepto?

Se generan N_toy conjuntos de 133 puntos con la ley EXACTA sigma ~ E^-2 (sin ningún efecto físico adicional) usando las
energías reales y los errores relativos reales de los datos, con ruido gaussiano. En cada conjunto se ajusta la recta
log-log y se estima alpha de dos maneras:
   (a) con el intercepto k (extrapolación a E = 1 GeV):      alpha_k  = sqrt(3 k /(4 pi C))      [solo válido si n = 2]
   (b) con la normalización en el pivote E0 (zona de datos):  alpha_p  = sqrt(3 sigma0 E0^2/(4 pi C))
Se compara cuánto fluctúa cada estimador respecto del valor verdadero. Cualquier fluctuación es puro ruido estadístico.

Uso: python3 src/montecarlo.py [--ntoy 40000]      (usa la semilla common.SEED; requiere haber corrido analysis.py)
Escribe: results/mc_toys.npz y añade la clave "mc" a results/results.json
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import besiii_model as bm, common

ap = argparse.ArgumentParser(); ap.add_argument("--ntoy", type=int, default=40000)
Ntoy = ap.parse_args().ntoy

rng = np.random.default_rng(common.SEED)
D = common.load_data()
E, S, dS = D["E"], D["S"], D["dS"]
res = json.load(open(common.RESULTS / "results.json"))
f = res["fits"]["all133"]

X = np.log10(E)
sigY = dS / (S * np.log(10))                    # errores relativos reales
w = 1 / sigY**2
Sw, Sx, Sxx = w.sum(), (w * X).sum(), (w * X * X).sum()
Dn = Sw * Sxx - Sx**2
X0 = Sx / Sw                                    # pivote
ytrue = np.log10(f["sigma0"]) - 2 * (X - X0)    # verdad: n = 2 exacto, con la normalización medida
Y = ytrue[None, :] + rng.normal(0, 1, (Ntoy, len(E))) * sigY[None, :]
Sy, Sxy = (Y * w).sum(1), (Y * w * X).sum(1)
a_hat = (Sw * Sxy - Sx * Sy) / Dn               # pendiente ajustada en cada toy
b_hat = (Sxx * Sy - Sx * Sxy) / Dn              # intercepto en cada toy
k_true = 10 ** (np.log10(f["sigma0"]) + 2 * X0)
alpha_true = np.sqrt(3 * k_true / (4 * np.pi * bm.C))
alpha_k = np.sqrt(3 * 10 ** b_hat / (4 * np.pi * bm.C))
alpha_p = np.sqrt(3 * 10 ** (a_hat * X0 + b_hat) * (10 ** X0) ** 2 / (4 * np.pi * bm.C))
dk, dp = alpha_k / alpha_true - 1, alpha_p / alpha_true - 1

res["mc"] = dict(Ntoy=Ntoy, std_k=float(dk.std() * 100), std_pivot=float(dp.std() * 100), amplif=float(dk.std() / dp.std()),
                 frac_k_gt_1pct=float((np.abs(dk) > 0.01).mean() * 100), frac_p_gt_1pct=float((np.abs(dp) > 0.01).mean() * 100),
                 sig_n_toy=float(a_hat.std()), rho_toy=float(np.corrcoef(-a_hat, b_hat)[0, 1]))
np.savez(common.RESULTS / "mc_toys.npz", dk=dk[:20000], dp=dp[:20000], n=-a_hat[:20000])
json.dump(res, open(common.RESULTS / "results.json", "w"), indent=1, default=float)
m = res["mc"]
print(f"N_toy={Ntoy}: desv. est. alpha_k = {m['std_k']:.3f} %, alpha_pivote = {m['std_pivot']:.3f} %  (amplificación {m['amplif']:.1f}x)")
print(f"fracción con |error| > 1 %: intercepto {m['frac_k_gt_1pct']:.1f} %, pivote {m['frac_p_gt_1pct']:.1f} %;  rho(n, log k) = {m['rho_toy']:.4f}")
