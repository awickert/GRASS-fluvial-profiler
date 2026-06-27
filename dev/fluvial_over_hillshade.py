import numpy as np
import networkx as nx
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from rivernetworkx.core import build_graph
from rivernetworkx.grass_io import _read_raster, sample_raster
from grass.pygrass.vector import VectorTopo
from grass.script import run_command, region as _region


def fit(la, ls, nbins=18, minper=8):
    if len(la) < 2 * minper:
        return None
    edges = np.linspace(la.min(), la.max(), nbins + 1); idx = np.digitize(la, edges); bc, bm = [], []
    for i in range(1, nbins + 1):
        sel = idx == i
        if sel.sum() >= minper:
            bc.append(np.median(la[sel])); bm.append(np.median(ls[sel]))
    bc, bm = np.array(bc), np.array(bm)
    if len(bc) < 6:
        return None
    cl = np.polyfit(bc, bm, 1); rl = float(np.sum((bm - np.polyval(cl, bc)) ** 2)); best = None
    for k in np.linspace(bc[1], bc[-2], 50):
        nL = int((bc <= k).sum())
        if nL < 2 or len(bc) - nL < 2:
            continue
        X = np.column_stack([np.ones_like(bc), np.maximum(0., bc - k)]); co, *_ = np.linalg.lstsq(X, bm, rcond=None)
        rss = float(np.sum((bm - X @ co) ** 2))
        if best is None or rss < best[0]:
            best = (rss, k, co[1])
    rss, k, sl = best; rng = bc.max() - bc.min()
    R = 1 - rss / rl if rl > 0 else 0.0; edge = min(k - bc.min(), bc.max() - k) / rng
    return dict(A_star=10 ** k, theta=-sl, R=R, edge=edge, ok=(R >= 0.3 and edge >= 0.1 and 0.2 <= -sl <= 1.5))


run_command('g.region', raster='dem', quiet=True)
r = _region(); cn = (float(r['n']) + float(r['s'])) / 2; ce = (float(r['e']) + float(r['w'])) / 2
ns = float(r['nsres']); ew = float(r['ewres'])
run_command('g.region', n=cn + 600 * ns, s=cn - 600 * ns, e=ce + 600 * ew, w=ce - 600 * ew, align='dem', quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='acc_h', drainage='ddir_h', overwrite=True, quiet=True)
run_command('r.slope.aspect', elevation='dem', slope='slp_h', format='percent', overwrite=True, quiet=True)
run_command('r.relief', input='dem', output='hs_h', overwrite=True, quiet=True)
acc, b = _read_raster('acc_h'); slp, _ = _read_raster('slp_h'); dem, _ = _read_raster('dem'); hs, _ = _read_raster('hs_h')
nrows, ncols = acc.shape
west, north, nsr, ewr = b['west'], b['north'], b['nsres'], b['ewres']
east = west + ncols * ewr; south = north - nrows * nsr
Aabs = np.abs(acc); Sg = slp / 100.0
ev = dem.ravel(); sv = Sg.ravel(); fnm = np.isfinite(ev) & np.isfinite(sv)
qs = np.percentile(ev[fnm], np.arange(2, 60, 2)); cut = qs[0]
for q in qs:
    band = fnm & (ev <= q) & (ev > q - (qs[1] - qs[0]))
    if band.sum() and (sv[band] < 1e-3).mean() < 0.4:
        cut = q; break
bluff = dem > cut
valid = np.isfinite(Aabs) & np.isfinite(Sg) & (Aabs > 0) & (Sg > 1e-3) & bluff
LA = np.where(valid, np.log10(np.where(Aabs > 0, Aabs, 1)), np.nan).ravel()
LS = np.where(valid, np.log10(np.where(Sg > 0, Sg, 1)), np.nan).ravel()
gA = fit(LA[np.isfinite(LA)], LS[np.isfinite(LS)], 30, 50); T = max(5, int(round(gA['A_star'] / 3.0)))
run_command('r.stream.extract', elevation='dem', accumulation='acc_h', threshold=T,
            stream_vector='str_h', stream_raster='srast_h', direction='sdir_h', d8cut=0, overwrite=True, quiet=True)
run_command('r.stream.basins', direction='sdir_h', stream_rast='srast_h', basins='bas_h', overwrite=True, quiet=True)
basin = _read_raster('bas_h')[0]; basin_r = basin.ravel()
srast = _read_raster('srast_h')[0]; srast_r = srast.ravel()
vt = VectorTopo('str_h'); vt.open('r'); recs = []
for ln in vt.viter('lines'):
    if ln.cat is None:
        continue
    en = ln.to_array(); x, y = en[:, 0], en[:, 1]
    keep = np.concatenate(([True], (np.diff(x) != 0) | (np.diff(y) != 0))); x, y = x[keep], y[keep]
    if len(x) < 2:
        continue
    Av = np.abs(sample_raster(acc, x, y, **b))
    if Av[0] > Av[-1]:
        x, y, Av = x[::-1], y[::-1], Av[::-1]
    recs.append({'cat': int(ln.cat), 'x': x, 'y': y, 'A': Av})
vt.close()
def key(x, y):
    return (round(float(x), 2), round(float(y), 2))
upi = {key(r['x'][0], r['y'][0]): r['cat'] for r in recs}
for r in recs:
    t = upi.get(key(r['x'][-1], r['y'][-1]), 0); r['tostream'] = 0 if t == r['cat'] else t
rbc = {r['cat']: r for r in recs}
G = build_graph(recs)
o = np.argsort(basin_r, kind='stable'); bs = basin_r[o]; uq, st = np.unique(bs, return_index=True); cells_of = {}
for i, u in enumerate(uq):
    if np.isfinite(u):
        e = st[i + 1] if i + 1 < len(st) else len(o); cells_of[int(u)] = o[st[i]:e]

seg2astar = {}; hx = []; hy = []; resolved = set()
for c in sorted(rbc, key=lambda c: rbc[c]['A'].max()):
    if c in resolved:
        continue
    anc = nx.ancestors(G, c) | {c}; idxs = [cells_of[a] for a in anc if a in cells_of]
    if not idxs:
        continue
    cell = np.concatenate(idxs); la = LA[cell]; ls = LS[cell]; ok = np.isfinite(la) & np.isfinite(ls)
    if ok.sum() < 60:
        continue
    f = fit(la[ok], ls[ok])
    if f and f['ok']:
        As = f['A_star']
        for s in anc:
            if s not in seg2astar:                       # finest (most-upstream) hollow wins
                seg2astar[s] = As
        rr = rbc[c]; Av = rr['A']
        if Av[0] <= As <= Av[-1]:
            j = int(np.argmax(Av >= As)); hx.append(rr['x'][j]); hy.append(rr['y'][j])
        else:
            hx.append(rr['x'][0] if As < Av[0] else rr['x'][-1]); hy.append(rr['y'][0] if As < Av[0] else rr['y'][-1])
        resolved |= anc
print('hollows=%d ; segments with a hollow A*: %d' % (len(hx), len(seg2astar)))

# fluvial network = stream cells downstream of their hollow (A >= local A*)
maxcat = int(np.nanmax(basin_r)) + 1
lut = np.full(maxcat + 1, np.nan)
for s, As in seg2astar.items():
    if 0 <= s <= maxcat:
        lut[s] = As
bint = np.where(np.isfinite(basin_r), basin_r, 0).astype(int)
astar_cell = lut[np.clip(bint, 0, maxcat)]
stream_cell = np.isfinite(srast_r) & (srast_r > 0)        # actual extracted channels
fluvial_cell = stream_cell & np.isfinite(astar_cell) & (Aabs.ravel() >= astar_cell)
print('channel cells=%d ; fluvial (downstream of hollow)=%d (%.0f%% of channel)'
      % (stream_cell.sum(), fluvial_cell.sum(), 100.0 * fluvial_cell.sum() / stream_cell.sum()))

stream2d = stream_cell.reshape(nrows, ncols); fluv2d = fluvial_cell.reshape(nrows, ncols)
ext = [west, east, south, north]
fig, ax = plt.subplots(figsize=(11, 9))
ax.imshow(hs, cmap='gray', extent=ext, origin='upper', interpolation='none')
ax.imshow(np.where(stream2d, 1.0, np.nan), cmap=matplotlib.colors.ListedColormap(['#cc9966']),
          extent=ext, origin='upper', alpha=0.45, interpolation='none')      # full extracted network (faint)
ax.imshow(np.where(fluv2d, 1.0, np.nan), cmap=matplotlib.colors.ListedColormap(['#1f77ff']),
          extent=ext, origin='upper', alpha=0.9, interpolation='none')        # fluvial (kept)
ax.scatter(hx, hy, s=10, c='red', edgecolors='k', linewidths=0.2, zorder=5, label='hollows (%d)' % len(hx))
ax.set_title('Trempealeau: fluvial network pruned at hollows (blue) over hillshade\n'
             'tan = full over-extracted network; red = hollow / colluvial-fluvial transition')
ax.legend(loc='upper right', fontsize=8); ax.set_xlabel('Easting'); ax.set_ylabel('Northing')
plt.tight_layout(); plt.savefig('/tmp/tremp_fluvial_hillshade.png', dpi=130)
print('saved /tmp/tremp_fluvial_hillshade.png')
run_command('g.remove', type='raster', name='acc_h,ddir_h,slp_h,srast_h,sdir_h,bas_h,hs_h', flags='f', quiet=True)
run_command('g.remove', type='vector', name='str_h', flags='f', quiet=True)
