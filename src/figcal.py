"""Calibración de ejes de la Fig. 3 (panel izquierdo) de BESIII, PRD 102, 112009, página 8 del PDF publicado.

Las posiciones (en puntos PDF) de las marcas mayores de los ejes se leyeron del propio PDF; los valores son los
rótulos impresos en la figura. Un ajuste lineal posición -> valor da la conversión de cualquier coordenada del
dibujo a GeV y nb.
"""
import numpy as np

PAGE_INDEX = 7           # página 8 (índice 7)
X_TICKS_PT = np.array([154.98, 172.48, 189.98, 207.477, 224.898, 242.395, 259.895, 277.316, 294.816])
X_TICKS_VAL = np.array([3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6])                       # GeV
Y_TICKS_PT = np.array([55.34, 64.613, 73.969, 83.32, 92.676, 102.031, 111.387, 120.742, 130.098, 139.371, 148.723])
Y_TICKS_VAL = np.array([6.2, 6.0, 5.8, 5.6, 5.4, 5.2, 5.0, 4.8, 4.6, 4.4, 4.2])             # nb


def calibrate():
    ax, bx = np.polyfit(X_TICKS_PT, X_TICKS_VAL, 1)
    ay, by = np.polyfit(Y_TICKS_PT, Y_TICKS_VAL, 1)
    res_x = float(np.max(np.abs(ax * X_TICKS_PT + bx - X_TICKS_VAL)))
    res_y = float(np.max(np.abs(ay * Y_TICKS_PT + by - Y_TICKS_VAL)))
    return dict(ax=ax, bx=bx, ay=ay, by=by, res_x=res_x, res_y=res_y, pt_per_nb=-1.0 / ay)


def panel_zone(page):
    """Zona de la página que contiene el panel izquierdo de la figura (mitad izquierda, 20 % superior)."""
    import pymupdf as fitz
    R = page.rect
    return fitz.Rect(0, 0, R.width * 0.5, R.height * 0.2)
