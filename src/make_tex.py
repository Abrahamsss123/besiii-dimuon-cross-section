"""PASO 6. Genera results/latex/numeros.tex (macros con las cifras del análisis) y las tablas .tex a partir de results.json y del CSV.
Uso: python3 src/make_tex.py"""
import sys, json
import numpy as np

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import besiii_model as bm, common

OUT = str(common.LATEX); common.LATEX.mkdir(parents=True, exist_ok=True)
r = json.load(open(common.RESULTS / "results.json"))
f = r["fits"]["all133"]


def fx(x, nd):
    return f"{x:.{nd}f}"


def pfmt(p):
    if p < 1e-10:
        return r"\ensuremath{<10^{-10}}"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return rf"\ensuremath{{{float(m):.1f}\times10^{{{int(e)}}}}}"
    if p < 0.01:
        return r"\ensuremath{<0.01}"
    return f"{p:.2f}"


M = {}  # macros
d = r["data"]
M.update(NDat=d["N"], NAlta=d["n_alta"], NBaja=d["n_baja"], NHid=d["n_hidden"], NAltaBar=d["n_alta"] - d["n_hidden"],
         RelBaja=fx(d["rel_med_baja"], 2), RelAlta=fx(d["rel_med_alta"], 2), RelMin=fx(d["rel_min"], 2), RelMax=fx(d["rel_max"], 2),
         MarkerBound=fx(d["marker_bound_nb"], 4), MarkerBoundPct=fx(100 * d["marker_bound_nb"] / 5.0, 2),
         Emin=fx(d["Emin"], 2), Emax=fx(d["Emax"], 2))
M.update(nAll=fx(f["n"], 3), snAll=fx(f["sig_n"], 3), snAllSc=fx(f["sig_n_scaled"], 3), chiAll=fx(f["chi2"], 1), ndfAll=int(f["ndf"]),
         chindfAll=fx(f["chi2"] / f["ndf"], 2), pAll=pfmt(f["p"]), Ezero=fx(f["E0"], 3), sigmaZero=fx(f["sigma0"], 3),
         sigmaZeroRel=fx(100 * f["rel_sigma0"], 3), rhoAll=fx(f["rho"], 4), kAll=fx(f["k"], 1), logkAll=fx(f["b"], 3), slogkAll=fx(f["sig_b"], 3),
         birgeAll=fx(f["birge"], 2), pullTwo=fx(r["n_model"]["pull_n2"], 1), pullTwoSc=fx(r["n_model"]["pull_n2_scaled"], 1),
         devNTwo=fx(100 * (f["n"] / 2 - 1), 1))
nm = r["n_model"]
M.update(nModel=fx(nm["I"], 3), nModelMin=fx(nm["min"], 3), nModelMax=fx(nm["max"], 3), nModelNoLow=fx(nm["sin_bajas"], 3),
         pullModel=fx(nm["pull_I"], 1), pullModelSc=fx(nm["pull_I_scaled"], 1), nLep=fx(r["vp"]["n_lep"], 3), pullLep=fx(nm["pull_nlep"], 1),
         dnLow=fx(nm["I"] - nm["sin_bajas"], 3))
a = r["alpha"]
M.update(alphaEff=fx(a["alpha_eff"] * 1e3, 3), alphaEffSig=fx(a["alpha_eff_sig"] * 1e3, 3), devPDG=fx(100 * a["dev_PDG"], 2),
         devPDGsig=fx(a["dev_PDG_sigma"], 1), alphaLepEz=fx(a["alpha_lep_E0"] * 1e3, 3), devLep=fx(100 * a["dev_lep"], 2),
         tensLep=fx(abs(a["tension_lep_sigma"]), 1), fcnt=fx(a["fcnt"], 2), fqed=fx(a["f_qed"], 2), fcntRatio=fx(a["fcnt_ratio"], 4),
         fcntExcess=fx(100 * (a["fcnt_ratio"] - 1), 1), alphaCntDev=fx(100 * a["alpha_cnt_dev"], 2), alphaCntStat=fx(100 * a["alpha_cnt_stat"], 2),
         alphaCntSyst=fx(100 * a["alpha_cnt_syst"], 2), sigRatioLep=fx(100 * a["sigma_ratio_lep"], 2), dalEz=fx(100 * a["dal_E0"], 2),
         relSigmaTot=fx(100 * a["rel_sigma0_tot"], 2))
M["alphaPDG"] = fx(bm.ALPHA0 * 1e3, 4)
rt = r["running_test"]
M.update(chiConst=fx(rt["chi2_const"], 1), chiRun=fx(rt["chi2_run"], 1), dChiRun=fx(-rt["delta_chi2"], 1), fcntFit=fx(rt["fcnt_const"], 2),
         ndfRun=rt["ndf"])
mc = r["mc"]
M.update(mcStdK=fx(mc["std_k"], 2), mcStdP=fx(mc["std_pivot"], 3), mcAmp=fx(mc["amplif"], 0), mcFracK=fx(mc["frac_k_gt_1pct"], 0),
         mcNtoy=f"{mc['Ntoy']:,}".replace(",", r"\,"), mcRho=fx(mc["rho_toy"], 4), mcSigN=fx(mc["sig_n_toy"], 3))
u = r["unc"]
M.update(uncEnergyPct=fx(u["energy_var_pct"], 1), vpChi=fx(r["ladder"]["vp_lep"]["chi2_free"], 1))
it = r["interference_total"]
M.update(intMin=fx(it["min"], 2), intMax=fx(it["max"], 2), intRms=fx(it["rms"], 2), lowAtLo=fx(it["low_only_at_3p85"], 2),
         lowAtHi=fx(it["low_only_at_4p59"], 2))
so = r["solutions"]
chis = [v["chi2_free"] for v in so["per"].values()]
M.update(solSpreadMed=fx(so["spread_med"], 2), solSpreadMax=fx(so["spread_max"], 2), solChiMin=fx(min(chis), 1), solChiMax=fx(max(chis), 1),
         solEmax=fx(so["E_at_max"], 2))
M.update(fISR=fx(r["isr"]["f_ISR"], 3), fISRpct=fx(100 * (r["isr"]["f_ISR"] - 1), 0), sigBCDRmodel=fx(r["isr"]["sigma_BCDR_model_4420"], 3),
         fsr=fx(100 * r["fsr"]["schwinger"], 2), massLo=fx(abs(r["mass"]["at_3p85"]) * 1e6, 1), massHi=fx(abs(r["mass"]["at_4p6"]) * 1e6, 1),
         dalMZ=fx(r["vp"]["MZ"], 4), dalLo=fx(100 * r["vp"]["at_3p85"], 2), dalHi=fx(100 * r["vp"]["at_4p6"], 2))
lad = r["ladder"]
full = lad["BESIII_completo_I"]
M.update(fullChi=fx(full["chi2_free"], 1), fullNdf=full["ndf_free"], fullChindf=fx(full["chi2_free"] / full["ndf_free"], 2), fullP=fx(full["p_free"], 2),
         fullChiFix=fx(full["chi2_fixed"], 1), fullNdfFix=full["ndf_fixed"], fullPFix=fx(full["p_fixed"], 2),
         fullN=fx(full["N"], 4), fullNsig=fx(full["sig_N"], 4), fullRms=fx(full["rms_resid"], 2),
         noLowChi=fx(lad["BESIII_sin_bajas"]["chi2_free"], 1), noLowChiFix=fx(lad["BESIII_sin_bajas"]["chi2_fixed"], 1),
         treeChi=fx(lad["arbol"]["chi2_free"], 1), treeN=fx(lad["arbol"]["N"], 4), treeMean=fx(lad["arbol"]["mean_resid"], 2))

with open(f"{OUT}/numeros.tex", "w") as fh:
    fh.write("% Generado automaticamente por make_tex.py a partir de results.json -- NO editar a mano\n")
    for k, v in M.items():
        fh.write(rf"\newcommand{{\{k}}}{{{v}}}" + "\n")

# ------------------------------------------------------------ tabla de ajustes (5 columnas: cabe en el ancho de página)
rows = [("all133", "Todos los puntos"), ("baja104", "Solo baja luminosidad"), ("alta29", "Solo alta luminosidad")]
with open(f"{OUT}/tab_ajustes.tex", "w") as fh:
    fh.write("\\begin{tabular}{lcccc}\n\\toprule\n")
    fh.write(r"Muestra & Puntos & Exponente efectivo $n$ & $\chi^2/\mathrm{ndf}$ & Valor $p$ \\" + "\n\\midrule\n")
    for k, lab in rows:
        x = r["fits"][k]
        fh.write(f"{lab} & {int(x['N'])} & ${x['n']:.3f}\\pm{x['sig_n']:.3f}$ & {x['chi2']:.1f}/{int(x['ndf'])} = {x['chi2']/x['ndf']:.2f} & {pfmt(x['p'])} \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# ------------------------------------------------------------ tratamiento de incertidumbres
with open(f"{OUT}/tab_incert.tex", "w") as fh:
    fh.write("\\begin{tabular}{lccc}\n\\toprule\n")
    fh.write(r"Incertidumbres incluidas en cada punto & Exponente efectivo $n$ & $\sigma_n$ & $\chi^2/\mathrm{ndf}$ \\" + "\n\\midrule\n")
    for k, lab in (("barra", "Solo la barra del gráfico"), ("barra_extr", "Barra del gráfico y error de extracción"),
                   ("completo", "Barra, error de extracción y error de energía")):
        x = r["unc"][k]
        fh.write(f"{lab} & {x['n']:.4f} & {x['sig_n']:.4f} & {x['chi2']/x['ndf']:.3f} \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# ------------------------------------------------------------ escalera de modelos
lab = [("arbol", r"QED a orden árbol, con $\alpha(0)$"),
       ("vp_lep", r"QED árbol con polarización del vacío leptónica"),
       ("vp_lep_fsr", r"Lo anterior más radiación del estado final"),
       ("continuo_BESIII", r"Solo el continuo de BESIII$^\dagger$"),
       ("BESIII_sin_bajas", r"Modelo de BESIII sin las seis resonancias ligeras$^\dagger$"),
       ("BESIII_completo_I", r"Modelo completo de BESIII (Solución I)$^\dagger$")]
with open(f"{OUT}/tab_escalera.tex", "w") as fh:
    fh.write("\\begin{tabular}{>{\\raggedright\\arraybackslash}p{6.4cm}ccc}\n\\toprule\n")
    fh.write(r"Modelo & Normalización $N$ & $\chi^2$ ($N$ libre) & $\chi^2$ ($N=1$) \\" + "\n\\midrule\n")
    for k, t_ in lab:
        x = lad[k]
        fh.write(f"{t_} & ${x['N']:.4f}\\pm{x['sig_N']:.4f}$ & {x['chi2_free']:.1f} & {x['chi2_fixed']:.1f} \\\\\n")
        if k == "vp_lep_fsr":
            fh.write("\\midrule\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# ------------------------------------------------------------ interferencias
names = [("rho(770)", r"$\rho(770)$"), ("omega(782)", r"$\omega(782)$"), ("phi(1020)", r"$\phi(1020)$"), ("J/psi", r"$J/\psi$"),
         ("psi(2S)", r"$\psi(2S)$"), ("psi(3770)", r"$\psi(3770)$"), ("psi(4040)", r"$\psi(4040)$"), ("psi(4160)", r"$\psi(4160)$"),
         ("psi(4415)", r"$\psi(4415)$"), ("S(4220)", r"$S(4220)$")]
with open(f"{OUT}/tab_interf.tex", "w") as fh:
    fh.write("\\begin{tabular}{lccc}\n\\toprule\n")
    fh.write(r"Estado & en $3.85$ GeV [\%] & en $4.59$ GeV [\%] & máx.\ $|\cdot|$ en el rango [\%] \\" + "\n\\midrule\n")
    for k, t in names:
        x = r["interference"][k]
        fh.write(f"{t} & ${x['at_3p85']:+.2f}$ & ${x['at_4p59']:+.2f}$ & {x['max_abs']:.2f} \\\\\n")
        if k == "psi(3770)":
            fh.write("\\midrule\n")
    fh.write("\\midrule\n")
    fh.write(f"Suma de los 6 estados ligeros & ${it['low_only_at_3p85']:+.2f}$ & ${it['low_only_at_4p59']:+.2f}$ & --- \\\\\n")
    fh.write(f"Suma total (modelo completo) & --- & --- & de ${it['min']:+.2f}$ a ${it['max']:+.2f}$ \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")

# ------------------------------------------------------------ tabla de los 133 puntos
dd = np.genfromtxt(common.DATA_CSV, delimiter=",", names=True, dtype=None, encoding="utf8")
E = dd["Ecm_GeV"].astype(float); S = dd["sigma_BCDR_nb"].astype(float); cl = dd["lumi_class"].astype(str)
er = np.array([np.nan if (x is None or np.isnan(x)) else x for x in dd["err_stat_nb"]], float)
cells = []
for i in range(len(E)):
    e_txt = r"$<$0.014" if np.isnan(er[i]) else f"{er[i]:.3f}"
    cells.append(f"{E[i]:.4f} & {S[i]:.3f} & {e_txt} & {'A' if cl[i]=='alta' else 'B'}")
nrow = int(np.ceil(len(cells) / 3))
with open(f"{OUT}/tab_datos.tex", "w") as fh:
    fh.write("\\begin{tabular}{rrrc|rrrc|rrrc}\n\\toprule\n")
    hdr = r"$E_{\rm cm}$ & $\sigma$ & $\delta\sigma$ & L"
    fh.write(" & ".join([hdr] * 3) + " \\\\\n\\midrule\n")
    for j in range(nrow):
        parts = []
        for c in range(3):
            idx = j + c * nrow
            parts.append(cells[idx] if idx < len(cells) else " & & & ")
        fh.write(" & ".join(parts) + " \\\\\n")
    fh.write("\\bottomrule\n\\end{tabular}\n")
print("tex ok:", len(M), "macros")
