"""Map the first-order valley catchments, bounded downstream by the first tributary
AREA-STEP (the universal signal). For each T-source, mouth = trunk node just above
its first significant tributary junction; valley catchment = upstream subtree(mouth)
(contiguous in the Braun-Willett SVector). Colour each valley over a hillshade;
green = field heads, red = first tributary junctions.
Run: PYTHONPATH=<repo>:<repo>/dev/divide_heads /usr/bin/python3 plot_catchments.py <mbr|feather_s> [T]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, ListedColormap
from rivernetworkx import dreich as D
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
T = float(sys.argv[2]) if len(sys.argv) > 2 else 5000.0
TRIB_FRAC = 0.20
cfg = SITES[key]
z, west, north, res, nd = read_flt(cfg['flt']); z = np.where(z == nd, np.nan, z)
field = load_field_heads(cfg['key'], cfg.get('nfilter'))
nr, nc = z.shape
cr = ((north - field[:, 1]) / res).astype(int); cc = ((field[:, 0] - west) / res).astype(int)
ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < nc)
R0, R1 = max(0, cr[ins].min() - 300), min(nr, cr[ins].max() + 300)
C0, C1 = max(0, cc[ins].min() - 300), min(nc, cc[ins].max() + 300)
zc = z[R0:R1, C0:C1]; X0, Y0 = west + C0 * res, north - R0 * res
H, W = zc.shape
dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32); field = field[ins]
print('routing %s (T=%.0f)...' % (cfg['name'], T), flush=True)
filled = D.fill(dem, nd, cellsize=res); fi = D.build_flowinfo(filled, nd, res)
D.contributing_area(fi); D.build_svector(fi)
ncon, recv, ro, co = fi['ncontrib'], fi['recv'], fi['row_of'], fi['col_of']
SV, SVI = fi['SVector'], fi['SVectorIndex']; pxa = res * res
Tc = T / pxa


def first_junction_mouth(src):
    """Return (mouth, junction) for the first tributary >= TRIB_FRAC of trunk below
    src, or (last-trunk-node, None) if none (valley runs to the domain edge)."""
    cur = src
    while True:
        nx = int(recv[cur])
        if nx == cur:
            return cur, None
        if ncon[nx] - ncon[cur] - 1 >= TRIB_FRAC * ncon[cur]:
            return cur, nx
        cur = nx


label = np.full(fi['N'], -1, dtype=np.int64)
junc_nodes = []
for lid, s in enumerate(sorted(D.get_sources(fi, Tc).tolist(), key=lambda s: -ncon[s])):
    mouth, jn = first_junction_mouth(int(s))
    block = SV[SVI[mouth]:SVI[mouth] + ncon[mouth]]
    label[block[label[block] == -1]] = lid                # first (largest) valley claims a cell
    if jn is not None:
        junc_nodes.append(jn)
lab = label.reshape(H, W)

ls = LightSource(azdeg=315, altdeg=45)
hs = ls.hillshade(np.where(np.isfinite(zc), zc, np.nanmin(zc)), vert_exag=2, dx=res, dy=res)
nlab = int(lab.max()) + 1
rng = np.arange(nlab); cols = plt.cm.tab20(rng % 20 / 20.0)
cmap = ListedColormap(cols)
masked = np.ma.masked_where(lab < 0, lab % 20)

fig, ax = plt.subplots(figsize=(12, 11))
ax.imshow(hs, cmap='gray', extent=[X0, X0 + W * res, Y0 - H * res, Y0], origin='upper')
ax.imshow(np.ma.masked_where(lab < 0, lab % 20), cmap='tab20', alpha=0.45, vmin=0, vmax=19,
          extent=[X0, X0 + W * res, Y0 - H * res, Y0], origin='upper')
ax.plot(field[:, 0], field[:, 1], 'o', mfc='none', mec='lime', mew=1.6, ms=9, label='field heads')
jr = np.array([ro[n] for n in junc_nodes]); jc = np.array([co[n] for n in junc_nodes])
ax.plot(X0 + (jc + .5) * res, Y0 - (jr + .5) * res, 'r.', ms=6, label='first tributary junctions')
ax.set_title('%s: first-order valley catchments (bounded by first tributary area-step, T=%.0f m$^2$), n=%d'
             % (cfg['name'], T, nlab))
ax.legend(loc='upper right', fontsize=9); ax.set_xlabel('E (m)'); ax.set_ylabel('N (m)')
plt.tight_layout(); plt.savefig('dev/divide_heads/figures/catchments_%s.png' % key, dpi=130)
print('wrote dev/divide_heads/figures/catchments_%s.png  (%d valleys)' % (key, nlab))
