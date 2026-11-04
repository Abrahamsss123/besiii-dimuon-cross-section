"""
PASO 3. Análisis: ajuste log-log, modelo completo de BESIII, escalera de modelos, resonancias, polarización del vacío.
(El Monte Carlo está en montecarlo.py.)   Uso: python3 src/analysis.py      Escribe results/results.json
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0, str(Path(__file__).resolve().parent))
import besiii_model as bm, fits, common
from common import SYS, SYS_EXTR

res = {}
D = common.load_data()
E, S, dS, cls, hidden = D["E"], D["S"], D["dS"], D["cls"], D["hidden"]
rel = dS / S
res["data"] = dict(N=len(E), n_alta=int((cls == "alta").sum()), n_baja=int((cls == "baja").sum()),
                   n_hidden=int(hidden.sum()), marker_bound_nb=common.MARKER_R_NB,
                   rel_med_baja=float(np.median(rel[(cls == "baja")]) * 100),
                   rel_med_alta=float(np.median(rel[(cls == "alta") & ~hidden]) * 100),
                   rel_min=float(rel.min() * 100), rel_max=float(rel.max() * 100),
                   Emin=float(E.min()), Emax=float(E.max()))


# ---- 1) ley de potencia  sigma = k E^-n  (incertidumbre no circular: la pendiente de la propagación es la AJUSTADA)
def fitres(mask, label):
    f = fits.wls_loglog(E[mask], S[mask], dS[mask])
    f["label"] = label
    f["N"] = int(mask.sum())
    f["birge"] = float(np.sqrt(f["chi2"] / f["ndf"]))
    f["sig_n_scaled"] = f["sig_n"] * max(1.0, f["birge"])
    f["k"] = 10 ** f["b"]
    f["sig_logk"] = f["sig_b"]
    return f


F = {"all133": fitres(np.ones(len(E), bool), "133 puntos"),
     "baja104": fitres(cls == "baja", "104 baja lum."),
     "alta29": fitres(cls == "alta", "29 alta lum.")}
res["fits"] = {k: {kk: (float(vv) if isinstance(vv, (np.floating, float, int, np.integer)) else vv) for kk, vv in v.items()} for k, v in F.items()}

# ---- 1b) efecto del tratamiento de incertidumbres sobre n (barra del gráfico, + extracción, + energía)
Ef = np.log(10)
variants = {"barra": fits.wls_loglog(E, S, D["err_stat"], dE=0.0),
            "barra_extr": fits.wls_loglog(E, S, dS, dE=0.0),
            "completo": F["all133"]}
dE_term = (F["all133"]["n"] * 0.0008 / (E * Ef)) ** 2
dY_tot = (dS / (S * Ef)) ** 2 + dE_term
res["unc"] = {k: dict(n=float(v["n"]), sig_n=float(v["sig_n"]), chi2=float(v["chi2"]), ndf=int(v["ndf"])) for k, v in variants.items()}
res["unc"]["energy_var_pct"] = float(np.median(dE_term / dY_tot) * 100)

# ---- 2) normalización alpha_eff en el pivote; comparación con alpha(0) y con el running leptónico
f = F["all133"]
alpha_eff = lambda sig0, E0: np.sqrt(3 * sig0 * E0**2 / (4 * np.pi * bm.C))
E0 = f["E0"]
a_eff = alpha_eff(f["sigma0"], E0)
rel_s0 = np.sqrt(f["rel_sigma0"] ** 2 + SYS**2 + SYS_EXTR**2)
res["alpha"] = dict(
    E0=float(E0), sigma0=float(f["sigma0"]), rel_sigma0_fit=float(f["rel_sigma0"]), rel_sigma0_tot=float(rel_s0),
    alpha_eff=float(a_eff), alpha_eff_sig=float(a_eff * rel_s0 / 2),
    alpha_PDG=bm.ALPHA0, dev_PDG=float(a_eff / bm.ALPHA0 - 1), dev_PDG_sigma=float((a_eff / bm.ALPHA0 - 1) / (rel_s0 / 2)),
    dal_E0=float(bm.delta_alpha_lep(E0)), alpha_lep_E0=float(bm.alpha_running(E0)),
    dev_lep=float(bm.alpha_running(E0) / bm.ALPHA0 - 1),
    sigma_ratio_lep=float((bm.alpha_running(E0) / bm.ALPHA0) ** 2 - 1),
    f_qed=float(4 * np.pi * bm.ALPHA0**2 / 3 * bm.C),
)
a = res["alpha"]
a["tension_lep_sigma"] = float((a["alpha_eff"] / bm.ALPHA0 - a["dev_lep"] - 1) / (rel_s0 / 2))
f_ratio = bm.F_CNT / a["f_qed"]
a["fcnt"] = bm.F_CNT
a["fcnt_ratio"] = float(f_ratio)
a["alpha_cnt_dev"] = float(np.sqrt(f_ratio) - 1)
a["alpha_cnt_stat"] = float(0.5 * 0.11 / bm.F_CNT)
a["alpha_cnt_syst"] = float(0.5 * SYS)

# ---- 3) escalera de modelos (normalización libre y fija)
models = {
    "arbol": bm.sigma_tree(E),
    "arbol_mass": bm.sigma_tree(E) * bm.beta_mu_factor(E),
    "vp_lep": bm.sigma_tree(E, 1.0) * bm.alpha_running(E) ** 2 * bm.beta_mu_factor(E),
    "vp_lep_fsr": bm.sigma_tree(E, 1.0) * bm.alpha_running(E) ** 2 * bm.beta_mu_factor(E) * (1 + 3 * bm.ALPHA0 / (4 * np.pi)),
    "continuo_BESIII": bm.F_CNT / E**2,
    "BESIII_completo_I": bm.sigma_D(E, "I"),
    "BESIII_sin_bajas": bm.sigma_D(E, "I", drop=tuple(bm.LOW.keys())),
}
lad = {}
for k, m in models.items():
    cf, cx = fits.chi2_norm_free(S, dS, m), fits.chi2_fixed(S, dS, m)
    r = S / m - 1
    lad[k] = dict(N=float(cf["N"]), sig_N=float(cf["sig_N"]), chi2_free=cf["chi2"], ndf_free=cf["ndf"], p_free=cf["p"],
                  chi2_fixed=cx["chi2"], ndf_fixed=cx["ndf"], p_fixed=cx["p"],
                  mean_resid=float(r.mean() * 100), rms_resid=float(r.std() * 100))
res["ladder"] = lad

# ---- 4) las 8 soluciones (casi degeneradas) de BESIII
sols, allm_ = {}, []
for s_ in bm.SOLUTIONS:
    m = bm.sigma_D(E, s_)
    allm_.append(m)
    cf, cx = fits.chi2_norm_free(S, dS, m), fits.chi2_fixed(S, dS, m)
    sols[s_] = dict(chi2_free=cf["chi2"], N=float(cf["N"]), chi2_fixed=cx["chi2"])
allm_ = np.array(allm_)
spread = (allm_.max(0) - allm_.min(0)) / allm_.mean(0) * 100
res["solutions"] = dict(per=sols, spread_max=float(spread.max()), spread_med=float(np.median(spread)), E_at_max=float(E[np.argmax(spread)]))

# ---- 5) contribución de cada resonancia (interferencia, % de sigma)
ip = bm.interference_percent(E, "I")
res["interference"] = {k: dict(at_3p85=float(bm.interference_percent(3.85, "I")[k]), at_4p59=float(bm.interference_percent(4.59, "I")[k]),
                               max_abs=float(np.max(np.abs(ip[k])))) for k in ip}
tot_int = bm.sigma_D(E, "I") / (bm.F_CNT / E**2) - 1
res["interference_total"] = dict(min=float(tot_int.min() * 100), max=float(tot_int.max() * 100), rms=float(tot_int.std() * 100),
                                 low_only_at_3p85=float(sum(bm.interference_percent(3.85, "I")[k] for k in bm.LOW)),
                                 low_only_at_4p59=float(sum(bm.interference_percent(4.59, "I")[k] for k in bm.LOW)))

# ---- 6) ISR / FSR / masa / polarización del vacío
sigD442 = float(bm.sigma_D(4.42, "I"))
res["isr"] = dict(sigma_obs_4420=5.818, sigma_BCDR_model_4420=sigD442, f_ISR=5.818 / sigD442)
res["fsr"] = dict(schwinger=float(3 * bm.ALPHA0 / (4 * np.pi)))
res["mass"] = dict(at_3p85=float(bm.beta_mu_factor(3.85) - 1), at_4p6=float(bm.beta_mu_factor(4.6) - 1),
                   naive_wrong=float(4 * 0.1056583755**2 / 4.0**2))
res["vp"] = dict(MZ=float(bm.delta_alpha_lep(91.1876)), at_3p85=float(bm.delta_alpha_lep(3.85)), at_4p6=float(bm.delta_alpha_lep(4.6)))
m_run = bm.sigma_tree(E, 1.0) * bm.alpha_running(E) ** 2
res["vp"]["n_lep"] = float(fits.wls_loglog(E, m_run, 0.007 * m_run)["n"])

# ---- 6b) pendiente que un ajuste de ley de potencia VERÍA si la física fuera el modelo completo de BESIII (mismos puntos y pesos)
nm = {s_: float(fits.wls_loglog(E, bm.sigma_D(E, s_), dS)["n"]) for s_ in bm.SOLUTIONS}
fm_low = fits.wls_loglog(E, bm.sigma_D(E, "I", drop=tuple(bm.LOW.keys())), dS)
fm_cnt = fits.wls_loglog(E, bm.F_CNT / E**2, dS)
Fa = F["all133"]
res["n_model"] = dict(per_solution=nm, I=nm["I"], min=min(nm.values()), max=max(nm.values()),
                      sin_bajas=float(fm_low["n"]), continuo=float(fm_cnt["n"]),
                      n_data=float(Fa["n"]), sig_n=float(Fa["sig_n"]),
                      pull_I=float((Fa["n"] - nm["I"]) / Fa["sig_n"]), pull_I_scaled=float((Fa["n"] - nm["I"]) / Fa["sig_n_scaled"]),
                      pull_n2=float((Fa["n"] - 2) / Fa["sig_n"]), pull_n2_scaled=float((Fa["n"] - 2) / Fa["sig_n_scaled"]),
                      pull_nlep=float((Fa["n"] - res["vp"]["n_lep"]) / Fa["sig_n"]))

# ---- 6c) prueba del running de alpha DENTRO del modelo completo
amp = bm.amplitudes(E, "I")
res_amp = sum(v for k, v in amp.items() if k != "continuo")
r_run = bm.alpha_running(E) / bm.alpha_running(f["E0"])


def _chi2(c, run):
    cont = amp["continuo"] * c * (r_run if run else 1.0)
    return float((((S - np.abs(cont + res_amp) ** 2) / dS) ** 2).sum())


o0 = minimize_scalar(lambda c: _chi2(c, False), bounds=(0.9, 1.1), method="bounded", options=dict(xatol=1e-10))
o1 = minimize_scalar(lambda c: _chi2(c, True), bounds=(0.9, 1.1), method="bounded", options=dict(xatol=1e-10))
res["running_test"] = dict(chi2_const=o0.fun, c_const=float(o0.x), fcnt_const=float(bm.F_CNT * o0.x**2),
                           chi2_run=o1.fun, c_run=float(o1.x), fcnt_run_at_E0=float(bm.F_CNT * o1.x**2),
                           delta_chi2=float(o0.fun - o1.fun), ndf=len(E) - 1)

common.RESULTS.mkdir(exist_ok=True)
json.dump(res, open(common.RESULTS / "results.json", "w"), indent=1, default=float)
print(f"n = {Fa['n']:.3f} ± {Fa['sig_n']:.3f}   chi2/ndf = {Fa['chi2']:.1f}/{int(Fa['ndf'])}   E0 = {Fa['E0']:.3f} GeV")
print("results/results.json escrito")
