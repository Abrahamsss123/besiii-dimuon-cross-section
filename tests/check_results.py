"""Comprobación rápida: los resultados clave coinciden con los valores publicados en el manuscrito (tolerancias estrechas)."""
import json, sys
from pathlib import Path
R = json.load(open(Path(__file__).resolve().parents[1] / "results" / "results.json"))
f = R["fits"]["all133"]
checks = [
    ("N de puntos", R["data"]["N"], 133, 0),
    ("n (133 puntos)", f["n"], 2.036, 0.0006),
    ("sigma_n", f["sig_n"], 0.010, 0.0006),
    ("chi2/ndf", f["chi2"] / f["ndf"], 1.70, 0.006),
    ("E0 [GeV]", f["E0"], 4.206, 0.0006),
    ("n modelo BESIII (Sol. I)", R["n_model"]["I"], 2.026, 0.0006),
    ("chi2 modelo completo (N libre)", R["ladder"]["BESIII_completo_I"]["chi2_free"], 119.6, 0.06),
    ("amplificación MC", R["mc"]["amplif"], 31.0, 0.6),
]
bad = 0
for name, got, exp, tol in checks:
    ok = abs(got - exp) <= tol
    bad += not ok
    print(f"{'OK ' if ok else 'FALLA'}  {name:34s} obtenido {got:10.4f}   esperado {exp:10.4f} ± {tol}")
sys.exit(1 if bad else 0)
