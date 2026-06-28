"""Estimate the hillslope length scale of the Mid Bailey Run 1 m DEM, to set a
geomorphically-motivated constant curvature window (Grieve et al., 2016; Roering
et al., 2010; Hurst et al., 2012): compute the spread (std and IQR) of tangential
curvature over a range of window radii and find the scaling break. Below the
break curvature is noise/microtopography-dominated (steep decline with window);
above it, hillslope-form-dominated (gentler). The break radius ~ hillslope scale.
"""
import numpy as np
from rivernetworkx import dreich as D

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
WINDOWS = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30]   # metres (cellsize = 1 m)


def main():
    z = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
    z = np.where(z == ND, np.nan, z)
    # a NaN-free block (curvature spread is distributional; a block suffices)
    S = 2048
    blk = None
    for r0 in (1500, 2500, 3500):
        for c0 in (3200, 4000, 2400):
            b = z[r0:r0 + S, c0:c0 + S]
            if b.shape == (S, S) and not np.isnan(b).any():
                blk = b.astype(np.float32)
                break
        if blk is not None:
            break
    filled = D.fill(blk, ND, 0.0001, 1.0)

    print(' window(m)   std(C_tan)     IQR(C_tan)   dlog(std)/dlog(w)')
    stds, iqrs = [], []
    for w in WINDOWS:
        cu = D.tangential_curvature(filled, ND, 1.0, window_radius=w)
        v = cu[cu != np.float32(ND)].astype(np.float64)
        stds.append(v.std())
        iqrs.append(np.percentile(v, 75) - np.percentile(v, 25))
    stds, iqrs = np.array(stds), np.array(iqrs)
    w = np.array(WINDOWS, float)
    slope = np.gradient(np.log(stds), np.log(w))           # local log-log slope
    for i, ww in enumerate(WINDOWS):
        print('%7d   %11.5g   %11.5g   %12.3f' % (ww, stds[i], iqrs[i], slope[i]))
    # the break: where the (negative) slope is steepest -> flattens just after
    brk = WINDOWS[int(np.argmin(slope))]
    print('\nsteepest-decline window (slope min) ~ %d m; '
          'hillslope scale is at/just above the break.' % brk)


if __name__ == '__main__':
    main()
