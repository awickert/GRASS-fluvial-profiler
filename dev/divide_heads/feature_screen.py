"""Feature screen (Mid Bailey Run, 1 m): which divide-network feature, read at the
apex above a channel head, actually distinguishes heads from same-area random
points? Tests Andy's "across from a significantly different / touching catchment"
family of signals against the control that killed the curvature.

For each feature we report AUC = P(feature at a head-apex > feature at a control
apex): 0.5 = no signal, ->1 heads higher, ->0 heads lower. Same-area control
throughout (the apex of a random head-scale channel cell).

Features at the apex (catchment top above each head):
  area_contrast   |Δlog10 area|  across the divide (discharge/scale difference)
  order_contrast  |Δ Strahler|   across the divide
  multiplicity_R  # distinct basins within R cells (how many catchments meet)
  relief          apex elev - head elev   (hillslope relief above the head)
  wrap            divide turning angle    (the curvature signal, for reference)

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 feature_screen.py
"""
import pickle
import numpy as np
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx import divides as DV
from rivernetworkx.dreich import _find_farthest_upslope
import divides_lib as L

WEST, NORTH = L.WEST, L.NORTH; ND = -9999.0; THRESH = 200
print('loading cache...', flush=True)
d = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))
full = d['filled']; clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
M = 180; r0, r1 = cr.min() - M, cr.max() + M; c0, c1 = cc.min() - M, cc.max() + M
dem = full[r0:r1, c0:c1].astype(np.float32); nr, nc = dem.shape
fi = D.build_flowinfo(dem, ND, 1.0)
D.contributing_area(fi); D.build_svector(fi); D.distance_from_outlet(fi)
row, col, ncon, NI = fi['row_of'], fi['col_of'], fi['ncontrib'], fi['NodeIndex']
sources = D.get_sources(fi, THRESH)
D.junction_network(fi, sources); SO = fi['SO']
_, lab = D.drainage_divides(fi, threshold=THRESH)
segs = D.channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
area, order = {}, {}
for s in segs:
    dr, dcj = s['cells'][-1]; nd = NI[dr, dcj]
    area[s['cat']] = float(ncon[nd]); order[s['cat']] = int(SO[nd]) if SO[nd] > 0 else 1

ok = lab >= 0
sigg = np.full((nr, nc), np.nan); org = np.full((nr, nc), np.nan)
sigg[ok] = np.log10([area.get(int(c), np.nan) for c in lab[ok]])
org[ok] = [order.get(int(c), 1) for c in lab[ok]]


def cross(field):
    out = np.zeros((nr, nc))
    for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        sh = np.roll(np.roll(field, di, 0), dj, 1); shl = np.roll(np.roll(lab, di, 0), dj, 1)
        diff = np.abs(field - sh)
        v = (lab >= 0) & (shl >= 0) & (lab != shl) & np.isfinite(diff)
        out = np.where(v & (diff > out), diff, out)
    return out


area_contrast = cross(sigg); order_contrast = cross(org)
elev = np.where(dem == np.float32(ND), np.nan, dem)


def multiplicity(ai, aj, R):
    out = np.empty(len(ai))
    for k in range(len(ai)):
        i, j = ai[k], aj[k]
        w = lab[max(0, i - R):i + R + 1, max(0, j - R):j + R + 1]
        out[k] = len(np.unique(w[w >= 0]))
    return out


# apexes for Clubb heads and same-area control
chan = np.where((ncon >= 50) & (ncon <= 5000))[0]
# clip-local (row,col) -> absolute coords: add the clip offset (c0, r0)
htree = cKDTree(np.column_stack([WEST + c0 + col[chan] + 0.5, NORTH - r0 - row[chan] - 0.5]))
sn, idx = htree.query(clubb); head_nodes = chan[idx]
print('head snap distance: med=%.1f max=%.1f m' % (np.median(sn), sn.max()))
rng = np.random.default_rng(0)
lo, hi = np.percentile(ncon[head_nodes], [10, 90])
pool = chan[(ncon[chan] >= lo) & (ncon[chan] <= hi)]
ctrl_nodes = pool[rng.integers(0, len(pool), size=400)]


def apex_feats(nodes):
    ap = np.array([_find_farthest_upslope(fi, int(n)) for n in nodes])
    ai, aj = row[ap], col[ap]
    return dict(
        area_contrast=area_contrast[ai, aj], order_contrast=order_contrast[ai, aj],
        mult15=multiplicity(ai, aj, 15), mult30=multiplicity(ai, aj, 30),
        relief=elev[ai, aj] - elev[row[nodes], col[nodes]],
    )


H = apex_feats(head_nodes); C = apex_feats(ctrl_nodes)


def auc(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    n = 0; s = 0.0
    for x in a:
        s += np.mean(b < x) + 0.5 * np.mean(b == x); n += 1
    return s / n


print('\nfeature           head_med  ctrl_med   AUC(head>ctrl)')
for f in ('area_contrast', 'order_contrast', 'mult15', 'mult30', 'relief'):
    print('%-16s  %7.2f  %7.2f    %.3f' % (f, np.median(H[f]), np.median(C[f]), auc(H[f], C[f])))
print('\n(AUC 0.5 = no signal; >0.6 worth developing; <0.4 heads lower)')
