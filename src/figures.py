"""PASO 5. Figuras del paper (Okabe-Ito, apta para daltonismo). Uso: python3 src/figures.py   -> figs/*.pdf, *.png"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import besiii_model as bm, fits, common

OUT = str(common.FIGS); common.FIGS.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 9.5, "axes.labelsize": 10, "legend.fontsize": 8.3, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
    "figure.dpi": 150, "savefig.dpi": 300, "pdf.fonttype": 42, "font.family": "DejaVu Sans",
})
# paleta Okabe-Ito (apta para daltonismo)
C_LOW, C_HIGH = "#4d4d4d", "#E69F00"
C_BES, C_TREE, C_RUN, C_FIT = "#0072B2", "#000000", "#CC79A7", "#009E73"

res = json.load(open(common.RESULTS / "results.json"))
D = common.load_data()
E, S, dS, cls = D["E"], D["S"], D["dS"], D["cls"]
lo, hi = cls == "baja", cls == "alta"
Eg = np.linspace(3.80, 4.60, 800)


def pts(ax, y, dy, scale=1.0):
    ax.errorbar(E[lo], y[lo] * scale, yerr=dy[lo] * scale, fmt="o", ms=3.4, mfc="white", mec=C_LOW, ecolor=C_LOW, elinewidth=0.9,
                capsize=0, label="BESIII, baja luminosidad (104)", zorder=3)
    ax.errorbar(E[hi], y[hi] * scale, yerr=dy[hi] * scale, fmt="o", ms=3.4, mfc=C_HIGH, mec="k", mew=0.4, ecolor=C_HIGH, elinewidth=0.9,
                capsize=0, label="BESIII, alta luminosidad (29)", zorder=4)


# ============================================================ Fig. 2: lineal + razon a QED arbol
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.6, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[2.3, 1.2], hspace=0.07))
pts(a1, S, dS)
a1.plot(Eg, bm.sigma_tree(Eg), "--", color=C_TREE, lw=1.4, label=r"QED árbol, $\alpha(0)$: $4\pi\alpha^2/3s$")
a1.plot(Eg, bm.sigma_tree(Eg, 1) * bm.alpha_running(Eg) ** 2, "-.", color=C_RUN, lw=1.4, label=r"QED árbol con $\alpha(s)$ (polarización del vacío)")
a1.plot(Eg, bm.sigma_D(Eg, "I"), "-", color=C_BES, lw=1.5, label="Modelo completo de BESIII (Sol. I)")
a1.set_ylabel(r"$\sigma$ [nb]"); a1.legend(frameon=False, loc="upper right")
a1.set_ylim(4.0, 6.5)
ref = bm.sigma_tree(E)
pts(a2, S / ref - 1, dS / ref, 100)
a2.plot(Eg, 100 * (bm.sigma_D(Eg, "I") / bm.sigma_tree(Eg) - 1), "-", color=C_BES, lw=1.5)
a2.plot(Eg, 100 * ((bm.sigma_tree(Eg, 1) * bm.alpha_running(Eg) ** 2) / bm.sigma_tree(Eg) - 1), "-.", color=C_RUN, lw=1.4)
a2.axhline(0, color=C_TREE, ls="--", lw=1.2)
a2.set_ylabel(r"$\sigma/\sigma_{\rm árbol}-1$ [%]"); a2.set_xlabel(r"$E_{\rm cm}$ [GeV]"); a2.set_xlim(3.79, 4.61)
fig.savefig(f"{OUT}/fig_lineal.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/fig_lineal.png", bbox_inches="tight"); plt.close(fig)

# ============================================================ Fig. 3: log-log + residuos con banda de covarianza completa
F = fits.wls_loglog(E, S, dS)
X = np.log10(E); Xg = np.log10(Eg)
fit_g = 10 ** (F["a"] * Xg + F["b"])
band = np.log(10) * fits.band_log10(F, Xg) * 100             # % (1 sigma)
k0 = F["sigma0"] * F["E0"] ** 2
fit_pts = 10 ** (F["a"] * X + F["b"])

# --- version limpia: ambos paneles comparten el eje x en GeV con escala logaritmica (log-log real)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.6, 6.4), sharex=True, gridspec_kw=dict(height_ratios=[2.0, 1.5], hspace=0.07))
a1.set_xscale("log"); a1.set_yscale("log")
a1.errorbar(E[lo], S[lo], yerr=dS[lo], fmt="o", ms=3.0, mfc="white", mec=C_LOW, ecolor=C_LOW, elinewidth=0.8, label="BESIII, baja lum. (104)")
a1.errorbar(E[hi], S[hi], yerr=dS[hi], fmt="o", ms=3.0, mfc=C_HIGH, mec="k", mew=0.4, ecolor=C_HIGH, elinewidth=0.8, label="BESIII, alta lum. (29)")
a1.plot(Eg, fit_g, "-", color=C_FIT, lw=1.6, label=rf"ajuste libre: $n={F['n']:.3f}\pm{F['sig_n']:.3f}$")
a1.plot(Eg, k0 / Eg**2, "--", color=C_TREE, lw=1.2, label=r"pendiente $-2$ (misma normalización en $E_0$)")
a1.plot(F["E0"], F["sigma0"], marker="D", ms=6, mfc="none", mec=C_FIT, mew=1.4, ls="none", label=rf"pivote $E_0={F['E0']:.3f}$ GeV")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
h, l = a1.get_legend_handles_labels()
h += [Patch(facecolor=C_FIT, alpha=0.22, label=r"banda $\pm1\sigma$ del ajuste (covarianza completa $n$–$\log k$)"),
      Line2D([], [], color=C_BES, lw=1.5, label="modelo completo de BESIII (Sol. I), ver panel inferior")]
l += [x.get_label() for x in h[-2:]]
a1.set_ylabel(r"$\sigma$ [nb]"); a1.legend(h, l, frameon=False, loc="lower left", fontsize=7.8)
a1.set_yticks([4.5, 5, 5.5, 6]); a1.set_yticklabels(["4.5", "5.0", "5.5", "6.0"]); a1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
a2.fill_between(Eg, -band, band, color=C_FIT, alpha=0.22, lw=0, label=r"banda $\pm1\sigma$ del ajuste (covarianza completa $n$–$\log k$)")
a2.errorbar(E[lo], 100 * (S[lo] / fit_pts[lo] - 1), yerr=100 * dS[lo] / fit_pts[lo], fmt="o", ms=3.0, mfc="white", mec=C_LOW, ecolor=C_LOW, elinewidth=0.8)
a2.errorbar(E[hi], 100 * (S[hi] / fit_pts[hi] - 1), yerr=100 * dS[hi] / fit_pts[hi], fmt="o", ms=3.0, mfc=C_HIGH, mec="k", mew=0.4, ecolor=C_HIGH, elinewidth=0.8)
a2.plot(Eg, 100 * (bm.sigma_D(Eg, "I") / fit_g - 1), "-", color=C_BES, lw=1.5, label="modelo completo de BESIII (Sol. I)")
a2.plot(Eg, 100 * ((k0 / Eg**2) / fit_g - 1), "--", color=C_TREE, lw=1.2, label=r"pendiente $-2$")
a2.axhline(0, color="gray", lw=0.6)
a2.set_ylabel(r"$\sigma/\sigma_{\rm ajuste}-1$ [%]"); a2.set_ylim(-3.2, 3.2)
a2.set_xlabel(r"$E_{\rm cm}$ [GeV] (escala logarítmica)")
a2.set_xlim(3.79, 4.62); a2.set_xticks([3.8, 4.0, 4.2, 4.4, 4.6]); a2.set_xticklabels(["3.8", "4.0", "4.2", "4.4", "4.6"]); a2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
fig.savefig(f"{OUT}/fig_loglog.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/fig_loglog.png", bbox_inches="tight"); plt.close(fig)

# ============================================================ Fig. 4: anatomia de las desviaciones (interferencias)
ipg = bm.interference_percent(Eg, "I")
low_sum = sum(ipg[k] for k in bm.LOW)
tot = 100 * (bm.sigma_D(Eg, "I") / (bm.F_CNT / Eg**2) - 1)
fig, ax = plt.subplots(figsize=(6.6, 4.4))
cnt = bm.F_CNT / E**2
ax.errorbar(E[lo], 100 * (S[lo] / cnt[lo] - 1), yerr=100 * dS[lo] / cnt[lo], fmt="o", ms=3.0, mfc="white", mec=C_LOW, ecolor="#9a9a9a", elinewidth=0.8, label="datos BESIII", zorder=2)
ax.errorbar(E[hi], 100 * (S[hi] / cnt[hi] - 1), yerr=100 * dS[hi] / cnt[hi], fmt="o", ms=3.0, mfc=C_HIGH, mec="k", mew=0.4, ecolor="#9a9a9a", elinewidth=0.8, zorder=3)
ax.plot(Eg, tot, "-", color=C_BES, lw=2.0, label="suma (modelo completo)", zorder=5)
for name, col, ls in (("psi(4040)", "#D55E00", "-"), ("psi(4160)", "#009E73", "-"), ("psi(4415)", "#CC79A7", "-"), ("S(4220)", "#E69F00", "--")):
    ax.plot(Eg, ipg[name], ls, color=col, lw=1.2, label=name.replace("psi", r"$\psi$") if name != "S(4220)" else r"$S(4220)$", zorder=4)
ax.plot(Eg, low_sum, "-", color="#56B4E9", lw=1.4, label=r"6 estados ligeros y $J/\psi,\psi(2S),\psi(3770)$", zorder=4)
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel(r"$E_{\rm cm}$ [GeV]"); ax.set_ylabel(r"$\sigma/(f_{\rm cnt}/s)-1$ [%]")
ax.set_xlim(3.79, 4.61); ax.set_ylim(-2.7, 3.0)
ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=7.6, columnspacing=1.2, handlelength=1.6)
fig.savefig(f"{OUT}/fig_resonancias.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/fig_resonancias.png", bbox_inches="tight"); plt.close(fig)

# ============================================================ Fig. 5: Monte Carlo (hipótesis nula n = 2 exacto)
m = np.load(common.RESULTS / "mc_toys.npz")
sk, sp, nn = m["dk"] * 100, m["dp"] * 100, m["n"]
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4))
ax = axs[0]
lim = max(np.abs(sk).max(), 1.0) * 1.02
ax.hist(sk, bins=np.linspace(-lim, lim, 90), color="#D55E00", alpha=0.75, density=True, label=rf"intercepto $k$ (desv. est. = {sk.std():.2f} %)")
ax.hist(sp, bins=np.linspace(-lim, lim, 600), color=C_BES, alpha=0.9, density=True, label=rf"pivote $E_0$ (desv. est. = {sp.std():.2f} %)")
ax.set_yscale("log"); ax.set_ylim(1e-3, 3e3); ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel(r"$\alpha_{\rm estimado}/\alpha_{\rm verdadero}-1$ [%]"); ax.set_ylabel("densidad (escala log)")
ax.set_title("(a) distribución del error", fontsize=9.5); ax.legend(frameon=False, fontsize=7.6, loc="upper right")
ax = axs[1]
sel = slice(0, 2500)
ax.plot(nn[sel] - 2, sk[sel], ".", ms=2.2, color="#D55E00", alpha=0.45, label=r"desde el intercepto $k$")
ax.plot(nn[sel] - 2, sp[sel], ".", ms=2.2, color=C_BES, alpha=0.6, label=r"desde el pivote $E_0$")
ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel(r"pendiente ajustada $n-2$ en cada simulación"); ax.set_ylabel(r"$\alpha_{\rm estimado}/\alpha_{\rm verdadero}-1$ [%]")
ax.set_title("(b) el error del intercepto sigue a la pendiente", fontsize=9.5); ax.legend(frameon=False, fontsize=7.6, loc="upper left")
fig.tight_layout()
fig.savefig(f"{OUT}/fig_mc.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/fig_mc.png", bbox_inches="tight"); plt.close(fig)
print("figuras ok")
