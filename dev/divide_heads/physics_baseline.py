"""Constructive counterpart to the (negative) divide-feature search: what DOES
locate channel heads, and does it survive coarsening? Tests the Montgomery-
Dietrich physics -- drainage area and local slope -- at Clubb's heads.

For each resolution: AUC of log10(area) for heads vs ANY-area channel cells
(does area separate heads? trivially yes), and AUC of local slope for heads vs
SAME-area channel cells (does slope add, at fixed area?). If slope discriminates
at 1 m but fades by 12 m (a derivative, like curvature), that is the argument for
chi-z (an integral) -- divides for robust structure, chi-z for robust location.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 physics_baseline.py
"""
import numpy as np
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
import divides_lib as L

ALG = '/tmp/dreich_algorithm'; NR, NC = 7107, 6266; ND = -9999.0
WEST, NORTH = L.WEST, L.NORTH
z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
z1 = np.where(z1 == ND, np.nan, z1)
clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
P = 260; R0, R1, C0, C1 = cr.min() - P, cr.max() + P, cc.min() - P, cc.max() + P
z1 = z1[R0:R1, C0:C1]
clubb_xy = np.column_stack([clubb[:, 0] - (WEST + C0), (NORTH - R0) - clubb[:, 1]])


def block_mean(z, n):
    nr2, nc2 = z.shape[0] // n, z.shape[1] // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def auc(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    return np.mean([np.mean(b < x) + 0.5 * np.mean(b == x) for x in a])


def run(res):
    agg = block_mean(z1, res) if res > 1 else z1
    dem = np.where(np.isfinite(agg), agg, ND).astype(np.float32)
    filled = D.fill(dem, ND, 0.0001, float(res))
    fi = D.build_flowinfo(filled, ND, float(res)); D.contributing_area(fi)
    row, col, ncon = fi['row_of'], fi['col_of'], fi['ncontrib']
    f = np.where(filled == np.float32(ND), np.nan, filled).astype(float)
    gy, gx = np.gradient(f, res)
    slope = np.hypot(gx, gy)                    # local gradient magnitude
    amin, amax = max(2, int(50 / res / res)), int(8000 / res / res)
    chan = np.where((ncon >= amin) & (ncon <= amax))[0]
    tree = cKDTree(np.column_stack([(col[chan] + 0.5) * res, (row[chan] + 0.5) * res]))
    _, idx = tree.query(clubb_xy); head = chan[idx]
    rng = np.random.default_rng(0)
    lo, hi = np.percentile(ncon[head], [10, 90])
    sameA = chan[(ncon[chan] >= lo) & (ncon[chan] <= hi)]
    sameA = sameA[rng.integers(0, len(sameA), size=400)]
    anyA = chan[rng.integers(0, len(chan), size=400)]
    sl_h = slope[row[head], col[head]]; sl_c = slope[row[sameA], col[sameA]]
    a_h = np.log10(ncon[head].astype(float)); a_any = np.log10(ncon[anyA].astype(float))
    return auc(a_h, a_any), auc(sl_h, sl_c), np.median(sl_h), np.median(sl_c)


print('res   area_AUC(vs any)   slope_AUC(vs same-area)   slope_head  slope_ctrl')
for res in (1, 2, 3, 5, 8, 12):
    aa, sa, sh, sc = run(res)
    print('%3d        %.3f                %.3f              %.3f      %.3f'
          % (res, aa, sa, sh, sc), flush=True)
