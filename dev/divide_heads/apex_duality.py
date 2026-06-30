"""Duality test (Mid Bailey Run, 1 m): is a channel head a characteristic
hillslope length L_H below its divide apex?

For each Clubb (2014) field head we take the apex above it to be the
farthest-upslope point of its own catchment (the back of the hollow), and measure
the head->apex flow distance L_H, the relief, and the straight-line distance. If
L_H is a TIGHT distribution (low coefficient of variation) there is a real
characteristic length and the head = apex - L_H duality holds; a scattered L_H
would weaken it. Compared against a null of random channel cells at the same
drainage-area scale.

Outputs: figures/apex_duality_hist.png, figures/apex_duality_map.png, and printed
statistics. Uses the cached 1 m MBR FlowInfo (/tmp/divides_proto_fi.pkl, rebuilt
by ../divides_setup.py) and Clubb's heads (/tmp/clubb_channel_heads.xlsx).

Run:  PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 apex_duality.py
"""
import os
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from rivernetworkx.dreich import _find_farthest_upslope
import divides_lib as L

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figures')
os.makedirs(FIG, exist_ok=True)
WEST, NORTH = L.WEST, L.NORTH

print('loading cache...', flush=True)
d = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))
fi = d['fi']; filled = d['filled']
row, col, fd, ncon = fi['row_of'], fi['col_of'], fi['fd'], fi['ncontrib']

# snap each Clubb head to the nearest HEAD-SCALE channel cell (avoid deep trunks)
AMIN, AMAX = 50, 5000                       # head-scale drainage area (cells); Clubb med ~769
hs = np.where((ncon >= AMIN) & (ncon <= AMAX))[0]
hx = WEST + (col[hs] + 0.5); hy = NORTH - (row[hs] + 0.5)
tree = cKDTree(np.column_stack([hx, hy]))
clubb = L.load_clubb()
snap, idx = tree.query(clubb)
head_nodes = hs[idx]
print('Clubb heads: %d; snap dist median=%.1f m (max %.1f m); snapped area med=%d cells'
      % (len(clubb), np.median(snap), snap.max(), np.median(ncon[head_nodes])))


def apex_stats(nodes):
    LH = np.empty(len(nodes)); relief = np.empty(len(nodes)); straight = np.empty(len(nodes))
    apex = np.empty(len(nodes), dtype=np.int64)
    for i, nd in enumerate(nodes):
        ap = _find_farthest_upslope(fi, int(nd)); apex[i] = ap
        LH[i] = fd[ap] - fd[nd]
        relief[i] = float(filled[row[ap], col[ap]] - filled[row[nd], col[nd]])
        straight[i] = np.hypot(row[ap] - row[nd], col[ap] - col[nd])
    return LH, relief, straight, apex


LH, relief, straight, apex = apex_stats(head_nodes)


def desc(x):
    return 'med=%.0f mean=%.0f std=%.0f CV=%.2f p10=%.0f p90=%.0f' % (
        np.median(x), x.mean(), x.std(), x.std() / x.mean(),
        np.percentile(x, 10), np.percentile(x, 90))


print('\nClubb heads -> apex:')
print('  L_H flow dist (m): %s' % desc(LH))
print('  relief (m):        %s' % desc(relief))
print('  straight (m):      %s' % desc(straight))

# null: random channel cells at the same area scale
lo, hi = np.percentile(ncon[head_nodes], [10, 90])
pool = hs[(ncon[hs] >= lo) & (ncon[hs] <= hi)]
rng = np.random.default_rng(0)
nn = pool[rng.integers(0, len(pool), size=min(500, len(pool)))]
nLH, _, _, _ = apex_stats(nn)
print('\nNULL (random channel cells, same area scale): L_H %s' % desc(nLH))
print('CV ratio (null / Clubb) = %.2f  (>1 means heads are tighter -> duality)'
      % ((nLH.std() / nLH.mean()) / (LH.std() / LH.mean())))

# ---- figure 1: L_H histograms, heads vs null ----
fig, ax = plt.subplots(figsize=(7, 4.5))
bins = np.linspace(0, np.percentile(np.r_[LH, nLH], 98), 30)
ax.hist(nLH, bins=bins, density=True, alpha=0.45, color='0.6', label='null channel cells')
ax.hist(LH, bins=bins, density=True, alpha=0.6, color='C3', label='Clubb field heads')
ax.axvline(np.median(LH), color='C3', ls='--', lw=1)
ax.set_xlabel('head -> apex flow distance  L_H  (m)'); ax.set_ylabel('density')
ax.set_title('MBR: hillslope length above channel heads (CV=%.2f) vs null (CV=%.2f)'
             % (LH.std() / LH.mean(), nLH.std() / nLH.mean()))
ax.legend()
fig.savefig(os.path.join(FIG, 'apex_duality_hist.png'), dpi=130, bbox_inches='tight')

# ---- figure 2: heads, apexes, and the connecting hillslope, over the DEM ----
rr = np.r_[row[head_nodes], row[apex]]; cc = np.r_[col[head_nodes], col[apex]]
r0, r1 = rr.min() - 60, rr.max() + 60; c0, c1 = cc.min() - 60, cc.max() + 60
sub = np.where(filled[r0:r1, c0:c1] == np.float32(-9999), np.nan, filled[r0:r1, c0:c1])
fig, ax = plt.subplots(figsize=(11, 11))
ax.imshow(sub, cmap='gray', extent=[WEST + c0, WEST + c1, NORTH - r1, NORTH - r0], origin='upper')
hxh = WEST + col[head_nodes]; hyh = NORTH - row[head_nodes]
axx = WEST + col[apex]; ayy = NORTH - row[apex]
for k in range(len(head_nodes)):
    ax.plot([hxh[k], axx[k]], [hyh[k], ayy[k]], '-', color='C0', lw=0.8, alpha=0.8)
ax.plot(hxh, hyh, 'o', color='C3', ms=4, label='Clubb head')
ax.plot(axx, ayy, '^', color='C0', ms=4, label='divide apex (farthest upslope)')
ax.legend(); ax.set_title('MBR: each Clubb head and its divide apex (hillslope above the head)')
fig.savefig(os.path.join(FIG, 'apex_duality_map.png'), dpi=130, bbox_inches='tight')
print('\nsaved figures/apex_duality_hist.png and figures/apex_duality_map.png')
