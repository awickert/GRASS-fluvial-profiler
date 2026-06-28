"""Calibrate a resolution-adaptive tangential-curvature threshold for DrEICH on
the Ohio test case (Mid Bailey Run 1 m LiDAR).

Problem: the standard tan_curv_threshold = 0.1 is 1 m-specific; averaged coarse
DEMs have far smaller curvature magnitudes, so 0.1 flags ~nothing beyond 1 m.

Scheme: hold the FRACTION of cells flagged "valley-like" constant. Anchor at the
fraction f that curvature > 0.1 flags at 1 m, then for each resolution set the
threshold to the (1 - f) quantile of THAT resolution's own curvature, so the same
fraction f of cells clears it. Curvature window held at 7 cells (window=7*res),
matching the sweep. DEM coarsening = block mean (matches r.resamp.stats average).

Writes "<res> <threshold>" lines to /tmp/dreich_curv_pctl.txt for
dev/dreich_resolution_sweep.sh adapt.
"""
import numpy as np
from rivernetworkx import dreich as D

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
RES = [1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 30]
OUT = '/tmp/dreich_curv_pctl.txt'


def block_mean(z, n):
    """Mean-aggregate a (NaN-aware) array into n x n blocks."""
    nr, nc = z.shape
    nr2, nc2 = nr // n, nc // n
    zt = z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n)
    return np.nanmean(zt, axis=(1, 3))


def curvature(zc, res):
    zf = np.where(np.isfinite(zc), zc, ND).astype(np.float32)
    filled = D.fill(zf, ND, 0.0001, float(res))
    cu = D.tangential_curvature(filled, ND, float(res), window_radius=7 * res)
    ok = (filled != np.float32(ND)) & np.isfinite(cu)
    return cu[ok].astype(np.float64)


def main():
    # 1 m anchor: use the validated faithful curvature field (no refit/fill).
    c1 = np.fromfile('%s/bailey_run_dem_tan_curv.flt' % ALG, '<f4').reshape(NR, NC)
    c1 = c1[c1 != np.float32(ND)].astype(np.float64)
    f = float(np.mean(c1 > 0.1))                      # fraction valley-like at 1 m
    print('1 m anchor: frac(curv>0.1) = %.5f  (threshold = %.1f pctile)'
          % (f, 100 * (1 - f)))

    z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
    z1 = np.where(z1 == ND, np.nan, z1)

    rows = []
    for n in RES:
        if n == 1:
            thr = 0.1
        else:
            cu = curvature(block_mean(z1, n), n)
            thr = float(np.quantile(cu, 1 - f))
        rows.append((n, thr))
        print('res=%2d m   tan_curv_threshold = %.6g' % (n, thr))

    with open(OUT, 'w') as fh:
        for n, thr in rows:
            fh.write('%d %.6g\n' % (n, thr))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
