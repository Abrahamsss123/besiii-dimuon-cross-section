"""
PASO 1. Extracción de los 133 puntos BCDR de la Fig. 3 (panel izquierdo) de BESIII, PRD 102, 112009 (2020).

No es digitalización manual: la figura del PDF es gráfico VECTORIAL, así que las coordenadas de cada marcador y de
cada barra de error están guardadas exactamente. El script las lee y las convierte a (E_cm [GeV], sigma [nb]).

  1. Calibración: ajuste lineal posición(pt) -> valor con las marcas mayores de los ejes (figcal.py).
  2. Marcadores: dibujos de 4 trazos y < 2 pt de lado en la zona del panel.
       - relleno rojo (1,0,0)   -> muestras de ALTA luminosidad
       - contorno negro, hueco  -> muestras de BAJA luminosidad
     Centro del marcador -> (E_cm, sigma).
  3. Barras de error: PyMuPDF no entrega las barras como trazos sueltos, así que la página se convierte a SVG con
     `pdftocairo -svg` y se leen los segmentos verticales (grosor 0.484 o 0.726 pt), aplicando la matriz de
     transformación de cada <path>. Cada marcador se asocia con los segmentos situados a < 0.06 pt de su abscisa;
     la barra es la media de la semi-longitud superior e inferior.
  4. Si la barra es menor que el marcador queda oculta (5 puntos de alta luminosidad): el CSV deja el campo vacío
     y el análisis usa como cota superior el radio del marcador (common.MARKER_R_NB).

Uso:   python3 src/extract_fig3.py --pdf PhysRevD_102_112009.pdf --out data/besiii_bcdr_133.csv
Requiere: PyMuPDF, NumPy y la utilidad `pdftocairo` (paquete poppler-utils).
Resultado esperado: 133 marcadores (104 baja, 29 alta). Si el PDF es otra versión (p. ej. arXiv) el número de página
o los grosores de línea pueden cambiar y el script lo detecta y se detiene.
"""
import argparse, hashlib, re, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
import pymupdf as fitz
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figcal

ap = argparse.ArgumentParser()
ap.add_argument("--pdf", required=True, help="PhysRevD_102_112009.pdf (versión publicada)")
ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "data" / "besiii_bcdr_133.csv"))
a = ap.parse_args()

print("[pdf] sha256 =", hashlib.sha256(open(a.pdf, "rb").read()).hexdigest())

# ---------- 1. calibración ----------
c = figcal.calibrate()
print(f"[calibración] x: {1/c['ax']:.3f} pt/GeV (residuo máx {c['res_x']:.1e} GeV) | "
      f"y: {c['pt_per_nb']:.3f} pt/nb (residuo máx {c['res_y']:.1e} nb)")

# ---------- 2. marcadores ----------
doc = fitz.open(a.pdf)
page = doc[figcal.PAGE_INDEX]
zone = figcal.panel_zone(page)
marks = []
for d in page.get_drawings():
    if not zone.intersects(d["rect"]):
        continue
    r = d["rect"]
    if r.width < 2 and r.height < 2 and len(d["items"]) == 4:
        if d["type"] == "f" and d.get("fill") == (1.0, 0.0, 0.0):
            marks.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2, "alta"))
        elif d["type"] == "s" and d.get("color") == (0.0, 0.0, 0.0):
            marks.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2, "baja"))
n_alta = sum(m[2] == "alta" for m in marks); n_baja = sum(m[2] == "baja" for m in marks)
print(f"[marcadores] total={len(marks)}  alta={n_alta}  baja={n_baja}")
if (len(marks), n_alta, n_baja) != (133, 29, 104):
    sys.exit("ERROR: se esperaban 133 marcadores (29 alta, 104 baja); ¿es el PDF publicado de PRD 102, 112009?")

# ---------- 3. barras de error (SVG) ----------
with tempfile.TemporaryDirectory() as tmp:
    svg = Path(tmp) / "p.svg"
    pg = str(figcal.PAGE_INDEX + 1)
    subprocess.run(["pdftocairo", "-svg", "-f", pg, "-l", pg, a.pdf, str(svg)], check=True)
    s = svg.read_text(encoding="utf8")
segs = []
for p in re.findall(r"<path[^>]*>", s):
    dm = re.search(r' d="([^"]+)"', p)
    if not dm:
        continue
    sw = re.search(r'stroke-width="([\d.]+)"', p)
    if not sw or sw.group(1) not in ("0.484", "0.726"):
        continue
    tr = re.search(r"matrix\(([^)]*)\)", p)
    A, B, C_, D, E_, F_ = map(float, tr.group(1).split(",")) if tr else (1, 0, 0, 1, 0, 0)
    for x1, y1, x2, y2 in re.findall(r"M\s*(-?[\d.e-]+)\s+(-?[\d.e-]+)\s+L\s*(-?[\d.e-]+)\s+(-?[\d.e-]+)", dm.group(1)):
        x1, y1, x2, y2 = map(float, (x1, y1, x2, y2))
        X1, Y1 = A * x1 + C_ * y1 + E_, B * x1 + D * y1 + F_
        X2, Y2 = A * x2 + C_ * y2 + E_, B * x2 + D * y2 + F_
        if abs(X1 - X2) < 1e-3 and 154 < X1 < 296 and 48 < min(Y1, Y2) and max(Y1, Y2) < 152:
            segs.append((X1, min(Y1, Y2), max(Y1, Y2)))
segs = np.array(segs)
print(f"[barras] segmentos verticales de barras de error: {len(segs)}")

rows = []
for cx, cy, kind in marks:
    near = segs[np.abs(segs[:, 0] - cx) < 0.06]
    up = near[(near[:, 2] <= cy + 0.15) & (near[:, 2] >= cy - 1.2)]    # segmento superior: termina cerca del centro
    dn = near[(near[:, 1] >= cy - 0.15) & (near[:, 1] <= cy + 1.2)]    # segmento inferior: empieza cerca del centro
    h_up = (cy - up[:, 1].min()) if len(up) else np.nan
    h_dn = (dn[:, 2].max() - cy) if len(dn) else np.nan
    h = np.nanmean([h_up, h_dn]) if not (np.isnan(h_up) and np.isnan(h_dn)) else np.nan
    rows.append((c["ax"] * cx + c["bx"], c["ay"] * cy + c["by"], h / c["pt_per_nb"] if not np.isnan(h) else np.nan, kind))

rows.sort()
out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
with open(out, "w") as fh:
    fh.write("Ecm_GeV,sigma_BCDR_nb,err_stat_nb,lumi_class\n")
    for E, sig, err, kind in rows:
        fh.write(f"{E:.5f},{sig:.5f},{'' if np.isnan(err) else f'{err:.5f}'},{kind}\n")
n_err = sum(not np.isnan(r[2]) for r in rows)
print(f"[salida] {len(rows)} puntos -> {out}  (con barra recuperada: {n_err}; barra oculta: {len(rows)-n_err})")
