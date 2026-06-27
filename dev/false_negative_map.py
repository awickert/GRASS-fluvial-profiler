import numpy as np
import networkx as nx
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
            best = (rss, k, co[0], co[1])
    rss, k, h, sl = best; rng = bc.max() - bc.min()
    R = 1 - rss / rl if rl > 0 else 0.0; edge = min(k - bc.min(), bc.max() - k) / rng
    return dict(A_star=10 ** k, logA_star=k, logS_star=h, theta=-sl, R=R, edge=edge,
                ok=(R >= 0.3 and edge >= 0.1 and 0.2 <= -sl <= 1.5))


run_command('g.region', raster='dem', quiet=True)
r = _region(); cn = (float(r['n']) + float(r['s'])) / 2; ce = (float(r['e']) + float(r['w'])) / 2
ns = float(r['nsres']); ew = float(r['ewres'])
run_command('g.region', n=cn + 600 * ns, s=cn - 600 * ns, e=ce + 600 * ew, w=ce - 600 * ew, align='dem', quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='acc_h', drainage='ddir_h', overwrite=True, quiet=True)
run_command('r.slope.aspect', elevation='dem', slope='slp_h', format='percent', overwrite=True, quiet=True)
acc, b = _read_raster('acc_h'); slp, _ = _read_raster('slp_h'); dem, _ = _read_raster('dem')
nrows, ncols = acc.shape
Aabs = np.abs(acc); Sg = slp / 100.0
ev = dem.ravel(); sv = Sg.ravel(); fn = np.isfinite(ev) & np.isfinite(sv)
qs = np.percentile(ev[fn], np.arange(2, 60, 2)); cut = qs[0]
for q in qs:
    band = fn & (ev <= q) & (ev > q - (qs[1] - qs[0]))
    if band.sum() and (sv[band] < 1e-3).mean() < 0.4:
        cut = q; break
bluff = dem > cut
valid = np.isfinite(Aabs) & np.isfinite(Sg) & (Aabs > 0) & (Sg > 1e-3) & bluff
LA = np.where(valid, np.log10(np.where(Aabs > 0, Aabs, 1)), np.nan).ravel()
LS = np.where(valid, np.log10(np.where(Sg > 0, Sg, 1)), np.nan).ravel()
validr = valid.ravel()
gA = fit(LA[np.isfinite(LA)], LS[np.isfinite(LS)], 30, 50); A_star = gA['A_star']; T = max(5, int(round(A_star / 3.0)))

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
o = np.argsort(basin_r, kind='stable'); bs = basin_r[o]
uniq, st = np.unique(bs, return_index=True); cells_of = {}
for i, u in enumerate(uniq):
    if np.isfinite(u):
        e = st[i + 1] if i + 1 < len(st) else len(o)
        cells_of[int(u)] = o[st[i]:e]

cats_byA = sorted(rbc, key=lambda c: rbc[c]['A'].max())
resolved = set(); hsig = []
for c in cats_byA:
    if c in resolved:
        continue
    anc = nx.ancestors(G, c) | {c}
    idxs = [cells_of[a] for a in anc if a in cells_of]
    if not idxs:
        continue
    cell = np.concatenate(idxs); la = LA[cell]; ls = LS[cell]; ok = np.isfinite(la) & np.isfinite(ls)
    la, ls = la[ok], ls[ok]
    if len(la) < 60:
        continue
    f = fit(la, ls)
    if f and f['ok']:
        rr = rbc[c]; Av = rr['A']; As = f['A_star']
        if Av[0] <= As <= Av[-1]:
            j = int(np.argmax(Av >= As)); xa, ya = rr['x'][j], rr['y'][j]
        else:
            xa, ya = (rr['x'][0], rr['y'][0]) if As < Av[0] else (rr['x'][-1], rr['y'][-1])
        hsig.append((float(xa), float(ya), f['logA_star'], f['logS_star']))
        resolved |= anc
print('hollows=%d  segments pruned=%d/%d' % (len(hsig), len(resolved), len(rbc)))

# ---- false-negative map ----
res_cells = np.zeros(LA.shape, bool)
for c in resolved:
    if c in cells_of:
        res_cells[cells_of[c]] = True
unres = validr & ~res_cells & np.isfinite(LA) & np.isfinite(LS)
ui = np.where(unres)[0]
rngs = np.random.default_rng(0)
samp = ui if len(ui) <= 30000 else rngs.choice(ui, 30000, replace=False)
ucol = (samp % ncols).astype(float); urow = (samp // ncols).astype(float)
ux = b['west'] + (ucol + 0.5) * b['ewres']; uy = b['north'] - (urow + 0.5) * b['nsres']
ula = LA[samp]; uls = LS[samp]
H = np.array(hsig); HX, HY, HlA, HlS = H[:, 0], H[:, 1], H[:, 2], H[:, 3]
tolA, tolS, Rsp = 0.30, 0.30, 300.0
best = np.full(len(samp), np.inf)
for j in range(len(H)):
    d2 = ((ula - HlA[j]) / tolA) ** 2 + ((uls - HlS[j]) / tolS) ** 2 \
        + (((ux - HX[j]) ** 2 + (uy - HY[j]) ** 2) / Rsp ** 2)
    best = np.minimum(best, d2)
fnmask = best < 1.0
nfn = int(fnmask.sum())
print('unresolved sampled cells: %d ; FN candidates: %d (%.1f%%)' % (len(samp), nfn, 100.0 * nfn / len(samp)))

# concentration: which unresolved branches carry the FN cells?
from collections import defaultdict
seg_of = basin_r[samp]
tot = defaultdict(int); fnc = defaultdict(int)
for s, ff in zip(seg_of, fnmask):
    if np.isfinite(s) and int(s) in rbc and int(s) not in resolved:
        tot[int(s)] += 1
        if ff:
            fnc[int(s)] += 1
fracs = sorted(((fnc[s] / n, s, n, fnc[s]) for s, n in tot.items() if n >= 8), reverse=True)
hi = [f for f in fracs if f[0] > 0.5]
totfn = sum(f[3] for f in fracs)
print('\nunresolved branches with >=8 sampled cells: %d' % len(fracs))
print('  median FN-fraction per branch: %.2f' % (np.median([f[0] for f in fracs]) if fracs else 0))
print('  branches with FN-fraction >0.5 (candidate missed hollows): %d' % len(hi))
print('  share of all FN cells sitting on those branches: %.0f%%' % (100.0 * sum(f[3] for f in hi) / max(1, totfn)))
print('\n  top FN-heavy branches, re-fit their subwatershed (R: ~0.25-0.30 = NARROW miss=real; <<0.3 = no structure=noise):')
for frac, s, n, nf in fracs[:15]:
    anc = nx.ancestors(G, s) | {s}; idxs = [cells_of[a] for a in anc if a in cells_of]
    cell = np.concatenate(idxs); la = LA[cell]; ls = LS[cell]; ok = np.isfinite(la) & np.isfinite(ls)
    f = fit(la[ok], ls[ok])
    print('     seg %5d  FNfrac=%.2f  nbasin=%6d  R=%s edge=%s theta=%s'
          % (s, frac, ok.sum(),
             ('%.2f' % f['R']) if f else 'none', ('%.2f' % f['edge']) if f else '-',
             ('%.2f' % f['theta']) if f else '-'))
run_command('g.remove', type='raster', name='acc_h,ddir_h,slp_h,srast_h,sdir_h,bas_h', flags='f', quiet=True)
run_command('g.remove', type='vector', name='str_h', flags='f', quiet=True)
