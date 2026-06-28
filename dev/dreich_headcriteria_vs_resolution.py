"""Supervised test (Ohio / Mid Bailey Run): given the KNOWN 1 m channel heads,
which coarse-DEM quantity is resolution-invariant at the heads -- i.e. which
thresholded criterion would reproduce the head positions on coarse data?

For each resolution we mean-aggregate the 1 m DEM, route it with the SAME
build_flowinfo routing that the DrEICH pipeline uses (internal consistency), and
read drainage area A, local slope S, and the slope-area product A*S^theta at each
known head (snapped to the local channel = max-accumulation cell in a 3x3 window).
Resolution-invariance of the median value at the heads = the criterion that
transfers; its DRIFT with resolution is the test ("the degradation is the test").
"""
import numpy as np
from rivernetworkx import dreich as D

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
WEST, NORTH, RES1 = 398767.32685636, 4369656.1523925, 1.0
RESLIST = [1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 30]
THETA = 0.5


def block_mean(z, n):
    nr, nc = z.shape
    nr2, nc2 = nr // n, nc // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def fields(filled, n):
    """Per-cell drainage area (m^2) and steepest-descent slope on a filled DEM."""
    fi = D.build_flowinfo(filled, ND, float(n))
    D.contributing_area(fi)
    nr, nc = filled.shape
    row, col, recv, flc = fi['row_of'], fi['col_of'], fi['recv'], fi['flc']
    z = filled.astype(np.float64)
    steplen = np.where(flc == 1, n, np.where(flc == 2, n * np.sqrt(2.0), np.nan))
    slope = (z[row, col] - z[fi['row_of'][recv], fi['col_of'][recv]]) / steplen
    A = np.full((nr, nc), np.nan); S = np.full((nr, nc), np.nan)
    A[row, col] = fi['ncontrib'].astype(float) * (n * n)
    S[row, col] = slope
    return A, S


def main():
    z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
    z1 = np.where(z1 == ND, np.nan, z1)
    ref = np.load('/tmp/dreich_ref_heads.npy')                 # known heads, E/N
    hcol = ((ref[:, 0] - WEST) / RES1).astype(int)
    hrow = ((NORTH - ref[:, 1]) / RES1).astype(int)
    prefill = np.fromfile('%s/bailey_run_dem_fill.flt' % ALG, '<f4').reshape(NR, NC)

    print('res  nH   medA(m2)  driftA   medS    driftS   med(A*S^.5)  drift   cvA  cvS  cvAS')
    base = {}
    for n in RESLIST:
        filled = prefill if n == 1 else D.fill(
            np.where(np.isfinite(block_mean(z1, n)), block_mean(z1, n), ND).astype(np.float32),
            ND, 0.0001, float(n))
        A, S = fields(filled, n)
        nr, nc = filled.shape
        hr, hc = hrow // n, hcol // n
        vals = []
        for r, c in zip(hr, hc):                               # snap to local channel
            if not (0 <= r < nr and 0 <= c < nc):
                continue
            r0, r1, c0, c1 = max(0, r - 1), min(nr, r + 2), max(0, c - 1), min(nc, c + 2)
            sub = A[r0:r1, c0:c1]
            if not np.isfinite(sub).any():
                continue
            dr, dc = np.unravel_index(np.nanargmax(sub), sub.shape)
            rr, cc = r0 + dr, c0 + dc
            a, s = A[rr, cc], S[rr, cc]
            if np.isfinite(a) and np.isfinite(s) and s > 0:
                vals.append((a, s, a * s ** THETA))
        v = np.array(vals)
        mA, mS, mAS = np.median(v, axis=0)
        cv = v.std(0) / v.mean(0)
        if n == 1:
            base = dict(A=mA, S=mS, AS=mAS)
        print('%3d %4d  %8.0f  %+5.0f%%  %.4f  %+5.0f%%  %10.2f  %+5.0f%%   %.2f %.2f %.2f'
              % (n, len(v), mA, 100 * (mA / base['A'] - 1), mS, 100 * (mS / base['S'] - 1),
                 mAS, 100 * (mAS / base['AS'] - 1), cv[0], cv[1], cv[2]))


if __name__ == '__main__':
    main()
