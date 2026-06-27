import numpy as np
import networkx as nx
from collections import Counter
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
acc, b = _read_raster('acc_h'); slp, _ = _read_raster('slp_h'); dem, _ = _read_raster('dem')
Aabs = np.abs(acc); Sg = slp / 100.0
ev = dem.ravel(); sv = Sg.ravel(); fnm = np.isfinite(ev) & np.isfinite(sv)
qs = np.percentile(ev[fnm], np.arange(2, 60, 2)); cut = qs[0]
for q in qs:
    band = fnm & (ev <= q) & (ev > q - (qs[1] - qs[0]))
    if band.sum() and (sv[band] < 1e-3).mean() < 0.4:
        cut = q; break
bluff = (dem > cut).ravel()
Ar = Aabs.ravel(); Sr = Sg.ravel(); logA_all = np.log10(np.where(Ar > 0, Ar, np.nan))
T = max(5, int(round(fit(logA_all[np.isfinite(logA_all) & bluff & (Sr > 1e-3)],
                         np.log10(np.where((Sr > 1e-3) & bluff & np.isfinite(logA_all), Sr, np.nan))[
                             np.isfinite(logA_all) & bluff & (Sr > 1e-3)], 30, 50)['A_star'] / 3.0)))
run_command('r.stream.extract', elevation='dem', accumulation='acc_h', threshold=T,
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
cats_byA = sorted(rbc, key=lambda c: rbc[c]['A'].max())
XMID = 402400.0
print('T=%d' % T)

def run_pass(floor):
    valid = np.isfinite(Aabs).ravel() & np.isfinite(Sg).ravel() & (Ar > 0) & (Sr > floor) & bluff
    LA = np.where(valid, np.log10(np.where(Ar > 0, Ar, 1)), np.nan)
    LS = np.where(valid, np.log10(np.where(Sr > 0, Sr, 1)), np.nan)
    resolved = set(); wp = ep = 0
    for c in cats_byA:
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
            resolved |= anc
            if float(rbc[c]['x'][0]) < XMID:
                wp += 1
            else:
                ep += 1
    return wp, ep

print('floor      WESTpass  EASTpass')
for fl in (1e-5, 1e-4, 3e-4, 1e-3):
    wp, ep = run_pass(fl)
    print('  %.0e      %5d     %5d' % (fl, wp, ep))
