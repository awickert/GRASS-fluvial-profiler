"""Apex detector test (Mid Bailey Run, 1 m): does the divide directly above a
channel head carry the geometric signature of an APEX -- a tight ridge wrap where
the valley's bounding divides converge -- detected from divide geometry ALONE?

Non-circular: the wrap signal is a property of the divide line (independent of the
heads). For each Clubb field head we take the above-head divide point (the
farthest-upslope point of its catchment -- on the ridge) and read off the divide
wrap there, comparing it to the wrap over the whole divide network and a random-
divide baseline. If above-head points sit in the high-wrap tail, the ridge-
convergence apex is a real, detectable head signature.

wrap(k) = 1 - chord(k-W, k+W) / (2W) along each divide polyline: 0 = straight
ridge, -> 1 at a U-turn (the valley-head pinch).

Outputs figures/apex_detect_*.png + stats. Clips around Clubb to iterate fast.
Run:  PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 apex_detect.py
"""
import os
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx import divides as DV
from rivernetworkx.dreich import _find_farthest_upslope
import divides_lib as L

HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, 'figures')
os.makedirs(FIG, exist_ok=True)
WEST, NORTH = L.WEST, L.NORTH
ND = -9999.0; THRESH = 200; W = 12        # wrap window (corners)

print('loading cache...', flush=True)
d = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))
full = d['filled']
clubb = L.load_clubb()
cr = ((NORTH - clubb[:, 1])).astype(int); cc = ((clubb[:, 0] - WEST)).astype(int)
M = 180
r0, r1 = cr.min() - M, cr.max() + M; c0, c1 = cc.min() - M, cc.max() + M
dem = full[r0:r1, c0:c1].astype(np.float32)
nr, nc = dem.shape; cw = nc + 1
print('clip %dx%d around %d Clubb heads' % (nr, nc, len(clubb)))

fi = D.build_flowinfo(dem, ND, 1.0)
D.contributing_area(fi); D.build_svector(fi); D.distance_from_outlet(fi)
row, col = fi['row_of'], fi['col_of']
_, lab = D.drainage_divides(fi, threshold=THRESH)
sources = D.get_sources(fi, THRESH)
segs = D.channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
outlet = np.zeros((nr, nc), bool)
for s in segs:
    rr, ccc = s['cells'][-1]; outlet[rr, ccc] = True
edges, meta = DV.extract_divide_edges(lab, stream=outlet)
fx = DV.diagonal_flow_corners(fi['recv'], row, col, nr, nc)
edges = DV.split_diagonal_crossings(edges, fx, cw)
g = DV.build_divide_graph(edges)
print('divide network: %d segments' % len(g['segments']))

# ---- wrap signal along each divide polyline: turning angle on a SMOOTHED line
# (smoothing kills the unit-step staircase that otherwise dominates) ----
def smooth(a, win):
    if len(a) < win:
        return a
    k = np.ones(win) / win
    return np.convolve(a, k, mode='same')

dc_i, dc_j, dc_w = [], [], []      # divide corner (i,j) and its wrap (|turn angle|, rad)
for p in g['segments']:
    ii, jj = DV.decode_corner(np.asarray(p) // 2, cw)
    n = len(ii)
    if n < 2 * W + 1:
        continue
    si = smooth(ii.astype(float), W); sj = smooth(jj.astype(float), W)
    for k in range(W, n - W):
        vin = np.array([si[k] - si[k - W], sj[k] - sj[k - W]])
        vout = np.array([si[k + W] - si[k], sj[k + W] - sj[k]])
        nin = np.hypot(*vin); nout = np.hypot(*vout)
        if nin < 1e-6 or nout < 1e-6:
            continue
        cosang = np.clip(np.dot(vin, vout) / (nin * nout), -1, 1)
        dc_i.append(ii[k]); dc_j.append(jj[k]); dc_w.append(np.arccos(cosang))
dc_i = np.array(dc_i); dc_j = np.array(dc_j); dc_w = np.array(dc_w)
dtree = cKDTree(np.column_stack([dc_i, dc_j]).astype(float))
print('divide-with-wrap points: %d (wrap med=%.2f p90=%.2f)'
      % (len(dc_w), np.median(dc_w), np.percentile(dc_w, 90)))

# ---- above-head divide points: catchment top of each Clubb head ----
# NOTE (corrected): the head snap must add the clip offset (c0, r0); the original
# version omitted it, snapping Clubb heads to the wrong cells and producing a
# SPURIOUS wrap signal. The corrected, controlled screen (feature_screen.py,
# resolution_screen.py) shows no divide feature distinguishes heads -- see README.
chan = np.where(fi['ncontrib'] >= 100)[0]
hx = WEST + c0 + (col[chan] + 0.5); hy = NORTH - r0 - (row[chan] + 0.5)
htree = cKDTree(np.column_stack([hx, hy]))
_, idx = htree.query(np.column_stack([clubb[:, 0], clubb[:, 1]]))
head_nodes = chan[idx]
ah_i = np.array([row[_find_farthest_upslope(fi, int(n))] for n in head_nodes])
ah_j = np.array([col[_find_farthest_upslope(fi, int(n))] for n in head_nodes])
dist, di = dtree.query(np.column_stack([ah_i, ah_j]).astype(float))
ah_wrap = dc_w[di]                  # wrap at the divide point nearest each above-head point
print('above-head point -> nearest divide-wrap point: dist med=%.1f p90=%.1f cells'
      % (np.median(dist), np.percentile(dist, 90)))

# random divide baseline
rng = np.random.default_rng(0)
base = dc_w[rng.integers(0, len(dc_w), size=2000)]
# DECISIVE CONTROL: catchment tops of random head-scale channel cells (same area)
lo, hi = np.percentile(fi['ncontrib'][head_nodes], [10, 90])
pool = chan[(fi['ncontrib'][chan] >= lo) & (fi['ncontrib'][chan] <= hi)]
nn = pool[rng.integers(0, len(pool), size=300)]
ct_i = np.array([row[_find_farthest_upslope(fi, int(n))] for n in nn])
ct_j = np.array([col[_find_farthest_upslope(fi, int(n))] for n in nn])
_, cdi = dtree.query(np.column_stack([ct_i, ct_j]).astype(float))
ctrl_wrap = dc_w[cdi]
print('\nwrap at above-head divide points:      med=%.2f mean=%.2f' % (np.median(ah_wrap), ah_wrap.mean()))
print('wrap at same-area control apexes:      med=%.2f mean=%.2f' % (np.median(ctrl_wrap), ctrl_wrap.mean()))
print('wrap over divide network (all points): med=%.2f mean=%.2f' % (np.median(base), base.mean()))
pct = np.array([(dc_w < w).mean() for w in ah_wrap])
cpct = np.array([(dc_w < w).mean() for w in ctrl_wrap])
print('percentile vs network: heads med=%.0f%%  control med=%.0f%%' % (100 * np.median(pct), 100 * np.median(cpct)))
# precision context: how common are high-wrap "apex" points?
for thr in (1.0, 1.5):
    print('  divide points with wrap>%.1f rad: %.1f%% of network (=%d points)'
          % (thr, 100 * np.mean(dc_w > thr), int(np.sum(dc_w > thr))))

# ---- figures ----
fig, ax = plt.subplots(figsize=(7, 4.5))
b = np.linspace(0, max(ah_wrap.max(), np.percentile(base, 99)), 25)
ax.hist(base, bins=b, density=True, alpha=0.45, color='0.6', label='divide network')
ax.hist(ah_wrap, bins=b, density=True, alpha=0.6, color='C3', label='above Clubb heads')
ax.set_xlabel('divide wrap (0 straight -> 1 U-turn)'); ax.set_ylabel('density')
ax.set_title('Divide wrap: above heads vs network (above-head median pct=%.0f%%)' % (100 * np.median(pct)))
ax.legend(); fig.savefig(os.path.join(FIG, 'apex_detect_hist.png'), dpi=130, bbox_inches='tight')

sub = np.where(dem == np.float32(ND), np.nan, dem)
fig, ax = plt.subplots(figsize=(11, 11))
ax.imshow(sub, cmap='gray', extent=[WEST + c0, WEST + c0 + nc, NORTH - r0 - nr, NORTH - r0], origin='upper')
sc = ax.scatter(WEST + c0 + dc_j, NORTH - r0 - dc_i, c=dc_w, s=2, cmap='viridis', vmin=0, vmax=np.percentile(dc_w, 98))
ax.plot(clubb[:, 0], clubb[:, 1], 'o', mfc='none', mec='r', ms=7, mew=1.4, label='Clubb head')
ax.plot(WEST + c0 + ah_j, NORTH - r0 - ah_i, '^', color='C1', ms=5, label='above-head divide pt')
plt.colorbar(sc, ax=ax, label='divide wrap', shrink=0.6); ax.legend()
ax.set_title('Divide network coloured by wrap; Clubb heads + their above-head points')
fig.savefig(os.path.join(FIG, 'apex_detect_map.png'), dpi=130, bbox_inches='tight')
print('saved figures/apex_detect_hist.png, apex_detect_map.png')
