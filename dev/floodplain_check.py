import numpy as np
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command


def fit_binned(la, ls, nbins=30, minper=50):
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
    return dict(A_star=10 ** k, theta=-sl, R=1 - rss / rl if rl > 0 else 0.0,
                edge=min(k - bc.min(), bc.max() - k) / rng)


run_command('g.region', raster='dem', quiet=True)
run_command('r.slope.aspect', elevation='dem', slope='slp_tmp', format='percent', overwrite=True, quiet=True)
acc, b = _read_raster('accum'); slp, _ = _read_raster('slp_tmp')
A = np.abs(acc).ravel(); S = (slp / 100.0).ravel()
fin = np.isfinite(A) & np.isfinite(S) & (A > 0)
A, S = A[fin], S[fin]
print('cells: %d' % len(A))
for thr in (1e-4, 1e-3, 1e-2):
    flat = S < thr
    if flat.sum():
        aq = np.percentile(np.log10(A[flat]), [50, 90])
        print('  S < %.0e : %6.2f%% of cells; their logA median/p90 = %.2f / %.2f'
              % (thr, 100.0 * flat.mean(), aq[0], aq[1]))

# binned median logS across logA, to see the high-A (floodplain) behavior
la = np.log10(A); ls = np.log10(np.where(S > 0, S, np.nan))
edges = np.linspace(la.min(), la.max(), 24); idx = np.digitize(la, edges)
print('\nbinned median logS by logA (watch high-logA for floodplain crash):')
for i in range(1, 24):
    sel = (idx == i) & np.isfinite(ls)
    if sel.sum() > 100:
        print('   logA~%.2f  medlogS=%+.2f  fracS<1e-3=%4.1f%%'
              % (np.median(la[sel]), np.median(ls[sel]), 100.0 * (S[(idx == i)] < 1e-3).mean()))

# fit under three cleaning regimes
def fit_regime(mask, label):
    a, s = A[mask], S[mask]
    ok = s > 0
    f = fit_binned(np.log10(a[ok]), np.log10(s[ok]))
    print('  %-28s A*=%8.0f theta=%.2f R=%.2f edge=%.2f  (n=%d)'
          % (label, f['A_star'], f['theta'], f['R'], f['edge'], ok.sum()))

print('\nfit sensitivity to floodplain cleaning:')
fit_regime(np.ones(len(A), bool), 'no filter')
fit_regime(S > 1e-3, 'S > 1e-3 (drop flats)')
fit_regime((S > 1e-3) & (A < np.percentile(A, 99)), 'S>1e-3 & A<p99 (drop trunk)')
run_command('g.remove', type='raster', name='slp_tmp', flags='f', quiet=True)
