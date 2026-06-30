"""Resolve the discrepancy: feature_screen (cache-filled clip) gave order_contrast
AUC 0.778; resolution_screen (raw-refilled clip) gave 0.444 at the SAME 1 m. Which
is real? Compute the feature both ways on the SAME clip extent + SAME heads, and
look at whether the apexes/features actually differ.
"""
import pickle
import numpy as np
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx.dreich import _find_farthest_upslope
import divides_lib as L

ALG = '/tmp/dreich_algorithm'; NR, NC = 7107, 6266; ND = -9999.0
WEST, NORTH = L.WEST, L.NORTH; THRESH = 200
clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
P = 260; R0, R1, C0, C1 = cr.min() - P, cr.max() + P, cc.min() - P, cc.max() + P

raw = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
raw = np.where(raw == ND, np.nan, raw)[R0:R1, C0:C1]
cache = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))['filled'][R0:R1, C0:C1].astype(np.float32)


def auc(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    return np.mean([np.mean(b < x) + 0.5 * np.mean(b == x) for x in a])


def setup(dem):
    fi = D.build_flowinfo(dem.astype(np.float32), ND, 1.0)
    D.contributing_area(fi); D.build_svector(fi); D.distance_from_outlet(fi)
    row, col, ncon, NI = fi['row_of'], fi['col_of'], fi['ncontrib'], fi['NodeIndex']
    nr, nc = dem.shape
    sources = D.get_sources(fi, THRESH); D.junction_network(fi, sources); SO = fi['SO']
    _, lab = D.drainage_divides(fi, threshold=THRESH)
    segs = D.channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
    order = {}
    for s in segs:
        dr, dj = s['cells'][-1]; order[s['cat']] = int(SO[NI[dr, dj]]) if SO[NI[dr, dj]] > 0 else 1
    org = np.full((nr, nc), np.nan); ok = lab >= 0
    org[ok] = [order.get(int(c), 1) for c in lab[ok]]
    oc = np.zeros((nr, nc))
    for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        sh = np.roll(np.roll(org, di, 0), dj, 1); shl = np.roll(np.roll(lab, di, 0), dj, 1)
        dd = np.abs(org - sh); v = (lab >= 0) & (shl >= 0) & (lab != shl) & np.isfinite(dd)
        oc = np.where(v & (dd > oc), dd, oc)
    chan = np.where((ncon >= 50) & (ncon <= 5000))[0]
    tree = cKDTree(np.column_stack([WEST + C0 + col[chan] + 0.5, NORTH - R0 - row[chan] - 0.5]))
    _, idx = tree.query(np.column_stack([clubb[:, 0], clubb[:, 1]]))
    head_nodes = chan[idx]
    rng = np.random.default_rng(0)
    lo, hi = np.percentile(ncon[head_nodes], [10, 90])
    pool = chan[(ncon[chan] >= lo) & (ncon[chan] <= hi)]
    ctrl = pool[rng.integers(0, len(pool), size=400)]
    aph = np.array([_find_farthest_upslope(fi, int(n)) for n in head_nodes])
    apc = np.array([_find_farthest_upslope(fi, int(n)) for n in ctrl])
    return oc[row[aph], col[aph]], oc[row[apc], col[apc]], fi


for name, dem in (('cache-filled', cache), ('raw-refilled', D.fill(np.where(np.isfinite(raw), raw, ND).astype(np.float32), ND, 0.0001, 1.0))):
    h, c, fi = setup(dem)
    print('%-14s order_contrast AUC=%.3f  head_med=%.2f ctrl_med=%.2f  head>0:%.2f ctrl>0:%.2f'
          % (name, auc(h, c), np.median(h), np.median(c), np.mean(h > 0), np.mean(c > 0)))
