"""Supervised resolution test, v2 (Ohio / Mid Bailey Run): at the KNOWN 1 m heads,
how do candidate criteria DRIFT with resolution? Fixes the v1 snapping artifact
(fixed PHYSICAL snap radius, not a fixed cell count) and adds the focus variable:

  chi(head) = integral_outlet^head (A_0 / A)^theta dx   along the flow path.

chi is NONLOCAL (an upstream integral of area, not a point value), so the
hypothesis is that it transfers across resolution better than point slope (which
we saw degrade ~39% by 30 m). Routing is the consistent build_flowinfo D8.
"""
import numpy as np
from rivernetworkx import dreich as D

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
WEST, NORTH = 398767.32685636, 4369656.1523925
RESLIST = [1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 30]
THETA, A_0 = 0.525, 1000.0       # DrEICH chi parameters
SNAP_M = 15.0                    # fixed physical snap radius (metres)


def block_mean(z, n):
    nr, nc = z.shape
    nr2, nc2 = nr // n, nc // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def main():
    z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
    z1 = np.where(z1 == ND, np.nan, z1)
    ref = np.load('/tmp/dreich_ref_heads.npy')
    hcol = ((ref[:, 0] - WEST)).astype(int)
    hrow = ((NORTH - ref[:, 1])).astype(int)
    prefill = np.fromfile('%s/bailey_run_dem_fill.flt' % ALG, '<f4').reshape(NR, NC)

    print('  R=%g m fixed snap; theta=%.3f A_0=%g' % (SNAP_M, THETA, A_0))
    print('res  nH   medA(m2) dA    medS    dS     med_chi   dchi    cvA  cvS  cvChi')
    base = {}
    for n in RESLIST:
        filled = prefill if n == 1 else D.fill(
            np.where(np.isfinite(block_mean(z1, n)), block_mean(z1, n), ND).astype(np.float32),
            ND, 0.0001, float(n))
        fi = D.build_flowinfo(filled, ND, float(n))
        D.contributing_area(fi)
        nr, nc = filled.shape
        row, col, recv, flc, NodeIndex = (fi['row_of'], fi['col_of'], fi['recv'],
                                          fi['flc'], fi['NodeIndex'])
        A_m2 = fi['ncontrib'].astype(float) * (n * n)
        z = filled.astype(np.float64)
        steplen = np.where(flc == 1, n, np.where(flc == 2, n * np.sqrt(2.0), 0.0))
        slope_nd = (z[row, col] - z[fi['row_of'][recv], fi['col_of'][recv]]) / np.where(steplen > 0, steplen, np.nan)
        incr = np.where(A_m2 > 0, (A_0 / A_m2) ** THETA * steplen, 0.0)   # chi integrand*dx
        Agrid = np.full((nr, nc), np.nan); Agrid[row, col] = A_m2

        w = int(round(SNAP_M / n))
        Avals, Svals, Cvals = [], [], []
        for r0, c0 in zip(hrow // n, hcol // n):
            if not (0 <= r0 < nr and 0 <= c0 < nc):
                continue
            ra, rb = max(0, r0 - w), min(nr, r0 + w + 1)
            ca, cb = max(0, c0 - w), min(nc, c0 + w + 1)
            sub = Agrid[ra:rb, ca:cb]
            if not np.isfinite(sub).any():
                continue
            dr, dc = np.unravel_index(np.nanargmax(sub), sub.shape)
            r, c = ra + dr, ca + dc
            nd = int(NodeIndex[r, c])
            a = A_m2[nd]; s = slope_nd[nd]
            if not (np.isfinite(a) and np.isfinite(s) and s > 0):
                continue
            # chi: walk downstream summing the integrand
            chi = 0.0; node = nd
            while recv[node] != node:
                chi += incr[node]; node = recv[node]
            Avals.append(a); Svals.append(s); Cvals.append(chi)
        A, S, C = np.array(Avals), np.array(Svals), np.array(Cvals)
        mA, mS, mC = np.median(A), np.median(S), np.median(C)
        if n == 1:
            base = dict(A=mA, S=mS, C=mC)
        print('%3d %4d  %8.0f %+4.0f%% %.4f %+5.0f%% %9.1f %+5.0f%%  %.2f %.2f %.2f'
              % (n, len(A), mA, 100 * (mA / base['A'] - 1), mS, 100 * (mS / base['S'] - 1),
                 mC, 100 * (mC / base['C'] - 1),
                 A.std() / A.mean(), S.std() / S.mean(), C.std() / C.mean()))


if __name__ == '__main__':
    main()
