"""The decisive test: do the divide head-signals survive coarsening? For each
resolution we aggregate the raw 1 m MBR DEM, route it, and recompute the AUC
(head-apex vs same-area control-apex) of each divide feature. A signal that holds
its AUC from 1 m to 12 m -- where cross-valley curvature collapses to ~0 -- is the
resolution-robust head detector we are after.

Features: order_contrast (|Δ Strahler| across the divide), area_contrast
(|Δlog10 area|), multiplicity (# basins within R, expect inverse), relief
(apex - head elevation).

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 resolution_screen.py
"""
import numpy as np
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx.dreich import _find_farthest_upslope
import divides_lib as L

ALG = '/tmp/dreich_algorithm'; NR, NC = 7107, 6266; ND = -9999.0
WEST, NORTH = L.WEST, L.NORTH; T_M2 = 200.0
RES = [1, 2, 3, 5, 8, 12]

z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
z1 = np.where(z1 == ND, np.nan, z1)
clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
P = 260
R0, R1, C0, C1 = cr.min() - P, cr.max() + P, cc.min() - P, cc.max() + P
z1 = z1[R0:R1, C0:C1]
cx0, cy0 = WEST + C0, NORTH - R0           # clip origin (m)
clubb_xy = np.column_stack([clubb[:, 0] - cx0, cy0 - clubb[:, 1]])   # metres from clip origin


def block_mean(z, n):
    nr2, nc2 = z.shape[0] // n, z.shape[1] // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def auc(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    return np.mean([np.mean(b < x) + 0.5 * np.mean(b == x) for x in a]) if len(a) else np.nan


def run(res):
    agg = block_mean(z1, res) if res > 1 else z1
    dem = np.where(np.isfinite(agg), agg, ND).astype(np.float32)
    filled = D.fill(dem, ND, 0.0001, float(res))
    fi = D.build_flowinfo(filled, ND, float(res)); D.contributing_area(fi)
    D.build_svector(fi); D.distance_from_outlet(fi)
    row, col, ncon, NI = fi['row_of'], fi['col_of'], fi['ncontrib'], fi['NodeIndex']
    nr, nc = dem.shape
    T = max(2, int(round(T_M2 / (res * res))))
    sources = D.get_sources(fi, T); D.junction_network(fi, sources); SO = fi['SO']
    _, lab = D.drainage_divides(fi, threshold=T)
    segs = D.channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
    area, order = {}, {}
    for s in segs:
        dr, dj = s['cells'][-1]; nd = NI[dr, dj]
        area[s['cat']] = float(ncon[nd]); order[s['cat']] = int(SO[nd]) if SO[nd] > 0 else 1
    ok = lab >= 0
    sigg = np.full((nr, nc), np.nan); org = np.full((nr, nc), np.nan)
    sigg[ok] = np.log10([area.get(int(c), np.nan) for c in lab[ok]])
    org[ok] = [order.get(int(c), 1) for c in lab[ok]]

    def cross(field):
        out = np.zeros((nr, nc))
        for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
            sh = np.roll(np.roll(field, di, 0), dj, 1); shl = np.roll(np.roll(lab, di, 0), dj, 1)
            dd = np.abs(field - sh); v = (lab >= 0) & (shl >= 0) & (lab != shl) & np.isfinite(dd)
            out = np.where(v & (dd > out), dd, out)
        return out
    ac, oc = cross(sigg), cross(org)
    elev = np.where(dem == np.float32(ND), np.nan, dem)

    # snap heads (metres -> cells at this res) to head-scale channel cells
    amin, amax = max(2, int(50 / res / res)), int(8000 / res / res)
    chan = np.where((ncon >= amin) & (ncon <= amax))[0]
    if len(chan) < 50:
        return None
    tree = cKDTree(np.column_stack([(col[chan] + 0.5) * res, (row[chan] + 0.5) * res]))
    _, idx = tree.query(clubb_xy); head_nodes = chan[idx]
    rng = np.random.default_rng(0)
    lo, hi = np.percentile(ncon[head_nodes], [10, 90])
    pool = chan[(ncon[chan] >= lo) & (ncon[chan] <= hi)]
    ctrl = pool[rng.integers(0, len(pool), size=400)]

    Rm = max(1, int(15 / res))

    def feats(nodes):
        ap = np.array([_find_farthest_upslope(fi, int(n)) for n in nodes])
        ai, aj = row[ap], col[ap]
        m = np.array([len(np.unique((lambda w: w[w >= 0])(
            lab[max(0, i - Rm):i + Rm + 1, max(0, j - Rm):j + Rm + 1]))) for i, j in zip(ai, aj)])
        return dict(order_contrast=oc[ai, aj], area_contrast=ac[ai, aj],
                    mult=m.astype(float), relief=elev[ai, aj] - elev[row[nodes], col[nodes]])
    H, C = feats(head_nodes), feats(ctrl)
    return {k: auc(H[k], C[k]) for k in H}


print('res   order_contrast  area_contrast  multiplicity  relief')
for res in RES:
    r = run(res)
    if r is None:
        print('%3d   (too few channel cells)' % res); continue
    print('%3d        %.3f          %.3f         %.3f      %.3f'
          % (res, r['order_contrast'], r['area_contrast'], r['mult'], r['relief']), flush=True)
