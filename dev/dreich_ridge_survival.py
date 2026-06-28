"""Ridge-survival test (Ohio / Mid Bailey Run): do drainage divides (ridgelines)
hold their position under coarsening, where DrEICH channel heads collapse?

For each resolution we mean-aggregate the 1 m DEM, route it with build_flowinfo,
and trace ridges with drainage_divides at a FIXED PHYSICAL basin area (so the
ridge scale is constant). We compare each resolution's ridge cells to the 1 m
ridge on a common 50 m reference grid (recall/precision of ridge-bearing blocks),
exactly the metric we used for heads -- where head recall fell 0.49 (2 m) -> ~0.

If ridges survive where heads do not, the ridge-bounded-valley strategy is viable.
"""
import numpy as np
from rivernetworkx import dreich as D

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
WEST, NORTH = 398767.32685636, 4369656.1523925
RESLIST = [1, 2, 3, 5, 8, 10, 15, 20, 30]
THRESH_M2 = 50000.0      # first-order-valley basin scale (~ DrEICH head scale)
REF = 50.0               # reference block size (m) for the position comparison


def block_mean(z, n):
    nr, nc = z.shape
    nr2, nc2 = nr // n, nc // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def ref_blocks(divide, n):
    """Set of 50 m reference-grid blocks that contain a ridge cell."""
    rr, cc = np.where(divide)
    x = WEST + (cc + 0.5) * n
    y = NORTH - (rr + 0.5) * n
    bx = ((x - WEST) / REF).astype(np.int64)
    by = ((NORTH - y) / REF).astype(np.int64)
    return set(zip(bx.tolist(), by.tolist()))


def main():
    z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
    z1 = np.where(z1 == ND, np.nan, z1)
    prefill = np.fromfile('%s/bailey_run_dem_fill.flt' % ALG, '<f4').reshape(NR, NC)

    print('  THRESH=%g m2 basins; %g m reference grid' % (THRESH_M2, REF))
    print('res  thr_cells  ridge_%%valid  ridge_blocks  recall  precision')
    base_blocks = None
    for n in RESLIST:
        filled = prefill if n == 1 else D.fill(
            np.where(np.isfinite(block_mean(z1, n)), block_mean(z1, n), ND).astype(np.float32),
            ND, 0.0001, float(n))
        fi = D.build_flowinfo(filled, ND, float(n))
        D.contributing_area(fi)
        thr = max(2, int(round(THRESH_M2 / (n * n))))
        divide, _ = D.drainage_divides(fi, threshold=thr)
        valid = np.isfinite(filled) & (filled != np.float32(ND))
        dens = 100.0 * divide.sum() / valid.sum()
        blk = ref_blocks(divide, n)
        if n == 1:
            base_blocks = blk
            rec = prec = 1.0
        else:
            inter = len(blk & base_blocks)
            rec = inter / len(base_blocks)
            prec = inter / max(1, len(blk))
        print('%3d  %8d  %10.1f  %12d  %.2f    %.2f'
              % (n, thr, dens, len(blk), rec, prec))


if __name__ == '__main__':
    main()
