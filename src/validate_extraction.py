"""
PASO 2. Validación de la extracción y del modelo (contra las curvas vectoriales que BESIII dibuja en la misma figura).

La Fig. 3 de BESIII contiene, además de los puntos, dos curvas: la azul (ajuste completo, Ec. 3 del artículo) y la cian
(continuo f_cnt/s). Se leen sus vértices del PDF, se calibran igual que los puntos y se comparan con nuestro modelo
(besiii_model.py). Un acuerdo a ~0.02 % valida a la vez: (i) la calibración de ejes y (ii) la implementación del modelo.
Además escribe docs/extraccion_overlay.png: los puntos extraídos (cruces) sobre la figura original, para inspección visual.

Uso:  python3 src/validate_extraction.py --pdf PhysRevD_102_112009.pdf
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pymupdf as fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import besiii_model as bm, figcal, common

ap = argparse.ArgumentParser()
ap.add_argument("--pdf", required=True)
a = ap.parse_args()

c = figcal.calibrate()
page = fitz.open(a.pdf)[figcal.PAGE_INDEX]
zone = figcal.panel_zone(page)
blue, cyan = [], []
for d in page.get_drawings():
    if not zone.intersects(d["rect"]) or d["type"] != "s":
        continue
    col = d.get("color")
    pts = []
    for it in d["items"]:
        if it[0] == "l":
            pts += [it[1], it[2]]
        elif it[0] == "c":
            pts += [it[1], it[4]]
    if col == (0.0, 0.0, 1.0):
        blue = pts
    elif col == (0.0, 1.0, 1.0):
        cyan = pts
Eb = np.array([c["ax"] * p.x + c["bx"] for p in blue]); Sb = np.array([c["ay"] * p.y + c["by"] for p in blue])
Ec = np.array([c["ax"] * p.x + c["bx"] for p in cyan]); Sc = np.array([c["ay"] * p.y + c["by"] for p in cyan])
print(f"curva azul: {len(Eb)} vértices;  curva cian: {len(Ec)} vértices")
print(f"calibración: residuo máx x = {c['res_x']:.1e} GeV, y = {c['res_y']:.1e} nb")

mc = (Ec > 3.805) & (Ec < 4.598)
rc = Sc[mc] / (bm.F_CNT / Ec[mc] ** 2) - 1
print(f"[continuo] cian vs f_cnt/s : media {100*rc.mean():+.3f}%  rms {100*rc.std():.3f}%  max|.| {100*np.abs(rc).max():.3f}%")
mb = (Eb > 3.805) & (Eb < 4.598)
for label, drop in (("completo (9 res + S)", ()), ("sin las 6 resonancias bajas", tuple(bm.LOW.keys()))):
    rb = Sb[mb] / bm.sigma_D(Eb[mb], "I", drop=drop) - 1
    print(f"[azul] vs modelo {label:30s}: media {100*rb.mean():+.3f}%  rms {100*rb.std():.3f}%  max|.| {100*np.abs(rb).max():.3f}%")
print(f"\nDelta_alpha_lep(M_Z=91.1876) = {float(bm.delta_alpha_lep(91.1876)):.5f}   (valor conocido ~0.0314)")
print(f"beta_mu factor a 4 GeV      = {float(bm.beta_mu_factor(4.0))-1:.2e}")

# ---- imagen de control: puntos extraídos sobre la figura original ----
D = common.load_data()
clip = fitz.Rect(140, 40, 305, 160)
zoom = 6
pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), clip=clip)
img = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)
x_pt = (D["E"] - c["bx"]) / c["ax"]; y_pt = (D["S"] - c["by"]) / c["ay"]
fig, ax = plt.subplots(figsize=(9, 6.5))
ax.imshow(img, extent=(clip.x0, clip.x1, clip.y1, clip.y0))
lo = D["cls"] == "baja"
ax.plot(x_pt[lo], y_pt[lo], "+", color="#D55E00", ms=7, mew=0.9, label="extraído (baja luminosidad)")
ax.plot(x_pt[~lo], y_pt[~lo], "x", color="#0072B2", ms=5, mew=0.9, label="extraído (alta luminosidad)")
ax.legend(loc="upper right", fontsize=8); ax.set_axis_off()
out = common.ROOT / "docs" / "extraccion_overlay.png"; out.parent.mkdir(exist_ok=True)
fig.savefig(out, dpi=170, bbox_inches="tight"); print("imagen de control ->", out)
