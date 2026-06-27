import numpy as np
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command, region as _region


def fit_binned(la, ls, nbins, minper):
    edges = np.linspace(la.min(), la.max(), nbins + 1); idx = np.digitize(la, edges); bc, bm = [], []
    for i in range(1, nbins + 1):
        sel = idx == i
        if sel.sum() >= minper:
            bc.append(np.median(la[sel])); bm.append(np.median(ls[sel]))
    bc, bm = np.array(bc), np.array(bm)
    if len(bc) < 6:
        return None
    cl = np.polyfit(bc, bm, 1); rl = float(np.sum((bm - np.polyval(cl, bc)) ** 2)); best = None
    for k in np.linspace(bc[1], bc[-2], 60):
        nL = int((bc <= k).sum())
        if nL < 2 or len(bc) - nL < 2:
            continue
        X = np.column_stack([np.ones_like(bc), np.maximum(0., bc - k)]); co, *_ = np.linalg.lstsq(X, bm, rcond=None)
        rss = float(np.sum((bm - X @ co) ** 2))
        if best is None or rss < best[0]:
            best = (rss, k, co[1])
    rss, k, sl = best; rng = bc.max() - bc.min()
    R = 1 - rss / rl if rl > 0 else 0.0; edge = min(k - bc.min(), bc.max() - k) / rng
    ok = (R >= 0.3 and edge >= 0.1 and 0.2 <= -sl <= 1.5)
    return dict(A_star=10 ** k, theta=-sl, R=R, edge=edge, ok=ok)


run_command('g.region', raster='dem', quiet=True)
r = _region(); cn = (float(r['n']) + float(r['s'])) / 2; ce = (float(r['e']) + float(r['w'])) / 2
ns = float(r['nsres']); ew = float(r['ewres'])
run_command('g.region', n=cn + 600 * ns, s=cn - 600 * ns, e=ce + 600 * ew, w=ce - 600 * ew, align='dem', quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='acc_w', drainage='ddir_w', overwrite=True, quiet=True)
run_command('r.slope.aspect', elevation='dem', slope='slp_w', format='percent', overwrite=True, quiet=True)
acc, b = _read_raster('acc_w'); slp, _ = _read_raster('slp_w')
Aabs = np.abs(acc); Sg = slp / 100.0
valid = np.isfinite(Aabs) & np.isfinite(Sg) & (Aabs > 0) & (Sg > 1e-4)
LA = np.where(valid, np.log10(np.where(Aabs > 0, Aabs, 1)), np.nan)
LS = np.where(valid, np.log10(np.where(Sg > 0, Sg, 1)), np.nan)
g = fit_binned(LA[valid], LS[valid], 30, 50)
A_star = g['A_star']; T = max(5, int(round(A_star / 3.0)))
print('global all-cell A*=%.0f R=%.2f -> T=%d' % (A_star, g['R'], T))

rows, cols = Aabs.shape; rng = np.random.default_rng(0)
def rc_to_en(rr, cc):
    return b['west'] + (cc + 0.5) * b['ewres'], b['north'] - (rr + 0.5) * b['nsres']

bands = [(T, 200), (200, 800), (800, 3000), (3000, 12000), (12000, 10 ** 12)]
print('\nall-cell subwatershed detection (r.water.outlet basins), by outlet drainage area:')
print('  band(cells)        nTested  passRate  medN     med A*')
for lo, hi in bands:
    rr, cc = np.where(valid & (Aabs >= lo) & (Aabs < hi))
    if len(rr) == 0:
        continue
    pick = rng.choice(len(rr), min(8, len(rr)), replace=False)
    npass = 0; ns_pts = []; astars = []
    for p in pick:
        E, N = rc_to_en(rr[p], cc[p])
        run_command('r.water.outlet', input='ddir_w', output='bas_w',
                    coordinates='%f,%f' % (E, N), overwrite=True, quiet=True)
        bas, _ = _read_raster('bas_w'); mask = np.isfinite(bas) & (bas > 0) & valid
        la = LA[mask]; ls = LS[mask]
        if len(la) < 60:
            continue
        ns_pts.append(len(la))
        f = fit_binned(la, ls, 20, 10)
        if f is not None:
            astars.append(f['A_star'])
            if f['ok']:
                npass += 1
    print('  %6d-%-9d %6d   %5.0f%%   %6d   %6.0f'
          % (lo, hi, len(pick), 100.0 * npass / len(pick),
             int(np.median(ns_pts)) if ns_pts else 0,
             np.median(astars) if astars else -1))
run_command('g.remove', type='raster', name='acc_w,ddir_w,slp_w,bas_w', flags='f', quiet=True)
