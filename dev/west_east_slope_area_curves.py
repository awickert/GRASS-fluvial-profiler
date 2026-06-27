import numpy as np
import networkx as nx
from rivernetworkx.core import build_graph
from rivernetworkx.grass_io import _read_raster, sample_raster
from grass.pygrass.vector import VectorTopo
from grass.script import run_command, region as _region

run_command('g.region', raster='dem', quiet=True)
r = _region(); cn = (float(r['n']) + float(r['s'])) / 2; ce = (float(r['e']) + float(r['w'])) / 2
ns = float(r['nsres']); ew = float(r['ewres'])
run_command('g.region', n=cn + 600 * ns, s=cn - 600 * ns, e=ce + 600 * ew, w=ce - 600 * ew, align='dem', quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='acc_h', drainage='ddir_h', overwrite=True, quiet=True)
run_command('r.slope.aspect', elevation='dem', slope='slp_h', format='percent', overwrite=True, quiet=True)
acc, b = _read_raster('acc_h'); slp, _ = _read_raster('slp_h'); dem, _ = _read_raster('dem')
Aabs = np.abs(acc); Sg = slp / 100.0
ev = dem.ravel(); sv = Sg.ravel(); fnm = np.isfinite(ev) & np.isfinite(sv)
qs = np.percentile(ev[fnm], np.arange(2, 60, 2)); cut = qs[0]
for q in qs:
    band = fnm & (ev <= q) & (ev > q - (qs[1] - qs[0]))
    if band.sum() and (sv[band] < 1e-3).mean() < 0.4:
        cut = q; break
bluff = (dem > cut).ravel()
valid = np.isfinite(Aabs).ravel() & np.isfinite(Sg).ravel() & (Aabs.ravel() > 0) & (Sg.ravel() > 1e-4) & bluff
LA = np.where(valid, np.log10(np.where(Aabs.ravel() > 0, Aabs.ravel(), 1)), np.nan)
LS = np.where(valid, np.log10(np.where(Sg.ravel() > 0, Sg.ravel(), 1)), np.nan)
run_command('r.stream.extract', elevation='dem', accumulation='acc_h', threshold=133,
            stream_vector='str_h', stream_raster='srast_h', direction='sdir_h', d8cut=0, overwrite=True, quiet=True)
run_command('r.stream.basins', direction='sdir_h', stream_rast='srast_h', basins='bas_h', overwrite=True, quiet=True)
basin_r = _read_raster('bas_h')[0].ravel()
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

def curve(c, label):
    anc = nx.ancestors(G, c) | {c}; idxs = [cells_of[a] for a in anc if a in cells_of]
    cell = np.concatenate(idxs); la = LA[cell]; ls = LS[cell]; ok = np.isfinite(la) & np.isfinite(ls)
    la, ls = la[ok], ls[ok]
    edges = np.linspace(la.min(), la.max(), 16); idx = np.digitize(la, edges)
    print('\n%s (cat %d, x=%.0f, n=%d): binned logA -> median logS' % (label, c, float(rbc[c]['x'][0]), len(la)))
    s = ''
    for i in range(1, 16):
        sel = idx == i
        if sel.sum() > 15:
            s += ' (%.2f,%+.2f)' % (np.median(la[sel]), np.median(ls[sel]))
    print('  ' + s)

# a large WEST subwatershed (failed) and a large EAST one
west = sorted([c for c in rbc if float(rbc[c]['x'][0]) < 402400 and rbc[c]['A'].max() < 50000],
              key=lambda c: rbc[c]['A'].max(), reverse=True)
east = sorted([c for c in rbc if float(rbc[c]['x'][0]) >= 402400 and rbc[c]['A'].max() < 50000],
              key=lambda c: rbc[c]['A'].max(), reverse=True)
curve(west[0], 'WEST bowl'); curve(west[1], 'WEST bowl 2')
curve(east[0], 'EAST bluff'); curve(east[1], 'EAST bluff 2')
run_command('g.remove', type='raster', name='acc_h,ddir_h,slp_h,srast_h,sdir_h,bas_h', flags='f', quiet=True)
run_command('g.remove', type='vector', name='str_h', flags='f', quiet=True)
