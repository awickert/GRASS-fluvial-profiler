"""Cross-divide catchment-contrast test (Mid Bailey Run, 1 m). Andy's idea: a
channel head is where you sit across the divide from a *significantly different*
catchment -- not the lateral side-by-side contact of two first-order valleys, but
where their UPSTREAM tips touch. So we read the contrast at the APEX (the upstream
point above each head), where the wrapping divide separates the small first-order
valley from whatever is behind it.

contrast(cell) = max over differently-labelled neighbours of |log10 A_here -
log10 A_there|, where A_x is the drainage area of the channel link that side
drains to (a discharge/scale proxy). Tested at above-head apexes vs same-area
random apexes (the same-area control that killed the curvature signal).

Outputs figures/cross_divide_*.png + stats. Run on a clip around Clubb's heads.
PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 cross_divide.py
"""
import os
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx.dreich import _find_farthest_upslope
import divides_lib as L

HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, 'figures')
os.makedirs(FIG, exist_ok=True)
WEST, NORTH = L.WEST, L.NORTH; ND = -9999.0; THRESH = 200

print('loading cache...', flush=True)
d = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))
full = d['filled']; clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
M = 180
r0, r1 = cr.min() - M, cr.max() + M; c0, c1 = cc.min() - M, cc.max() + M
dem = full[r0:r1, c0:c1].astype(np.float32); nr, nc = dem.shape
fi = D.build_flowinfo(dem, ND, 1.0)
D.contributing_area(fi); D.build_svector(fi); D.distance_from_outlet(fi)
row, col, ncon, NI = fi['row_of'], fi['col_of'], fi['ncontrib'], fi['NodeIndex']
_, lab = D.drainage_divides(fi, threshold=THRESH)

# per-link (basin) significance = drainage area at the link's downstream end
sources = D.get_sources(fi, THRESH)
segs = D.channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
sig = {}
for s in segs:
    dr, dcc = s['cells'][-1]; sig[s['cat']] = float(ncon[NI[dr, dcc]])
sigg = np.full((nr, nc), np.nan)             # log10 area of the link each cell drains to
ok = lab >= 0
cats = lab[ok]
sigg[ok] = np.log10([sig.get(int(ca), np.nan) for ca in cats])

# cross-divide contrast per cell = max |Δ log area| with a differently-labelled
# neighbour (np.roll wraps the clip edges; the interior test points are unaffected)
contrast = np.zeros((nr, nc))
for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
    sh = np.roll(np.roll(sigg, di, 0), dj, 1)
    shl = np.roll(np.roll(lab, di, 0), dj, 1)
    diff = np.abs(sigg - sh)
    valid = (lab >= 0) & (shl >= 0) & (lab != shl) & np.isfinite(diff)
    contrast = np.where(valid & (diff > contrast), diff, contrast)
print('contrast field: nonzero cells=%d, med(nonzero)=%.2f p90=%.2f'
      % (np.sum(contrast > 0), np.median(contrast[contrast > 0]), np.percentile(contrast[contrast > 0], 90)))

# snap Clubb heads to head-scale channel cells; above-head apex = catchment top.
# NOTE (corrected): must add the clip offset (c0, r0); the original omitted it and
# the head result was a snapping artifact -- see README and feature_screen.py.
chan = np.where((ncon >= 50) & (ncon <= 5000))[0]
htree = cKDTree(np.column_stack([WEST + c0 + col[chan] + 0.5, NORTH - r0 - row[chan] - 0.5]))
_, idx = htree.query(clubb); head_nodes = chan[idx]
ap = np.array([_find_farthest_upslope(fi, int(n)) for n in head_nodes])
ah_contrast = contrast[row[ap], col[ap]]

rng = np.random.default_rng(0)
lo, hi = np.percentile(ncon[head_nodes], [10, 90])
pool = chan[(ncon[chan] >= lo) & (ncon[chan] <= hi)]
nn = pool[rng.integers(0, len(pool), size=300)]
apn = np.array([_find_farthest_upslope(fi, int(n)) for n in nn])
ctrl_contrast = contrast[row[apn], col[apn]]

print('\ncross-divide contrast (|Δlog10 area|) at apex:')
print('  above Clubb heads:        med=%.2f mean=%.2f frac>1=%.2f' % (np.median(ah_contrast), ah_contrast.mean(), np.mean(ah_contrast > 1)))
print('  same-area control apexes: med=%.2f mean=%.2f frac>1=%.2f' % (np.median(ctrl_contrast), ctrl_contrast.mean(), np.mean(ctrl_contrast > 1)))
allc = contrast[contrast > 0]
print('  divide network at large:  med=%.2f mean=%.2f' % (np.median(allc), allc.mean()))

# ---- figures ----
fig, ax = plt.subplots(figsize=(7, 4.5))
b = np.linspace(0, max(np.percentile(ah_contrast, 98), np.percentile(ctrl_contrast, 98), 1), 25)
ax.hist(ctrl_contrast, bins=b, density=True, alpha=0.5, color='0.6', label='same-area control apexes')
ax.hist(ah_contrast, bins=b, density=True, alpha=0.6, color='C3', label='above Clubb heads')
ax.set_xlabel('cross-divide contrast  |Δ log10 catchment area|'); ax.set_ylabel('density')
ax.set_title('Cross-divide contrast at apexes: heads (med=%.2f) vs control (med=%.2f)'
             % (np.median(ah_contrast), np.median(ctrl_contrast)))
ax.legend(); fig.savefig(os.path.join(FIG, 'cross_divide_hist.png'), dpi=130, bbox_inches='tight')

sub = np.where(dem == np.float32(ND), np.nan, dem)
fig, ax = plt.subplots(figsize=(11, 11))
ax.imshow(sub, cmap='gray', extent=[WEST + c0, WEST + c0 + nc, NORTH - r0 - nr, NORTH - r0], origin='upper')
yy, xx = np.where(contrast > 0)
sc = ax.scatter(WEST + c0 + xx, NORTH - r0 - yy, c=contrast[yy, xx], s=2, cmap='magma',
                vmin=0, vmax=np.percentile(allc, 97))
ax.plot(clubb[:, 0], clubb[:, 1], 'o', mfc='none', mec='c', ms=7, mew=1.4, label='Clubb head')
ax.plot(WEST + c0 + col[ap], NORTH - r0 - row[ap], '^', color='lime', ms=4, label='above-head apex')
plt.colorbar(sc, ax=ax, label='cross-divide contrast', shrink=0.6); ax.legend()
ax.set_title('Divides coloured by cross-divide catchment contrast; Clubb heads + apexes')
fig.savefig(os.path.join(FIG, 'cross_divide_map.png'), dpi=130, bbox_inches='tight')
print('saved figures/cross_divide_hist.png, cross_divide_map.png')
