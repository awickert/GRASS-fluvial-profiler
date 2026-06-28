"""Faithful Python port of LSDRaster::calculate_polyfit_surface_metrics tangential
curvature (raster_selection[6]), window_radius=7, on LSDTT's filled DEM.

Transcribed from /tmp/LSDTopoTools2/src/LSDRaster.cpp:
  - calculate_polyfit_surface_metrics @3930 (mask radial_dist<=window_radius, no floor)
  - 2nd-order surface z = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
  - x = (krow-kr)*res (along ROWS), y = (kcol-kr)*res (along COLS)
  - A (6x6) constant Gram matrix; bb summed over the CIRCULAR mask; solve A*coeffs=bb
  - fx=d fy=e fxx=2a fyy=2b fxy=c ; p=fx^2+fy^2 ; q=p+1
  - tan_curv = (fxx*fy^2 - 2*fxy*fx*fy + fyy*fx^2)/(p*sqrt(q))   when q>0 and p*sqrt(q)!=0
A is constant, so each bb field is a correlation of the DEM with a fixed 15x15 kernel.
Validate cell-for-cell vs the instrumented C++ raster bailey_run_dem_tan_curv.flt.
"""
import numpy as np
from scipy import ndimage

NR, NC = 7107, 6266
ND = -9999.0
RES = 1.0
WINDOW_RADIUS = 7.0

z = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_fill.flt', dtype='<f4').reshape(NR, NC).astype(np.float64)
valid_data = z != ND

kr = int(np.ceil(WINDOW_RADIUS / RES))            # 7
kw = 2 * kr + 1                                   # 15
ii = np.arange(kw)
xk = (ii[:, None] - kr) * RES * np.ones((1, kw))  # x varies along kernel rows
yk = np.ones((kw, 1)) * (ii[None, :] - kr) * RES  # y varies along kernel cols
# dreich_algorithm/LSDRaster.cpp:1070 uses floor(radial_dist) <= window_radius
# (the standalone repo the driver #includes; NOT LSDTopoTools2's radial_dist<=wr).
mask = (np.floor(np.sqrt(xk * xk + yk * yk)) <= WINDOW_RADIUS).astype(np.float64)

# Constant A (6x6), accumulated over masked positions exactly as the C++ does.
x = xk[mask == 1]; y = yk[mask == 1]
terms = [x**2, y**2, x * y, x, y, np.ones_like(x)]     # basis: x^2 y^2 xy x y 1
A = np.array([[np.sum(terms[r] * terms[cc]) for cc in range(6)] for r in range(6)])
Ainv = np.linalg.inv(A)

# bb kernels (centered): bb_k[i,j] = sum_window K_k * z  -> correlation with K_k.
# K_k = basis_k * mask  (zero outside the circle).
K = [(x**2), (y**2), (x * y), x, y, np.ones_like(x)]   # values at masked cells
kernels = []
for kvals in K:
    ker = np.zeros((kw, kw))
    ker[mask == 1] = kvals
    kernels.append(ker)

# Only coeffs 0..4 (a b c d e) are needed for tangential curvature.
bb = [ndimage.correlate(z, ker, mode='constant', cval=0.0) for ker in kernels]
bb = np.stack(bb, axis=0)                              # (6, NR, NC)
coef = np.tensordot(Ainv[:5], bb, axes=([1], [0]))     # (5, NR, NC) = a,b,c,d,e

a, b, c, d, e = coef
fx, fy, fxx, fyy, fxy = d, e, 2 * a, 2 * b, c
p = fx * fx + fy * fy
q = p + 1.0
denom = p * np.sqrt(q)
with np.errstate(invalid='ignore', divide='ignore'):
    tanc = (fxx * fy * fy - 2 * fxy * fx * fy + fyy * fx * fx) / denom

# Validity: not edge, center is data, and NO NoData anywhere in the kw x kw square.
interior = np.zeros((NR, NC), bool)
interior[kr:NR - kr, kr:NC - kr] = True
nodata_in_square = ndimage.maximum_filter((~valid_data).astype(np.uint8), size=kw) > 0
good = interior & valid_data & (~nodata_in_square)
# C++ sets NoData where q<=0 or p*sqrt(q)==0 as well
good &= (q > 0) & (denom != 0)

tan_out = np.full((NR, NC), ND, dtype=np.float32)
tan_out[good] = tanc[good].astype(np.float32)

# ---- validate vs C++ reference ----
ref = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(NR, NC)
ref_valid = ref != np.float32(ND)
mine_valid = tan_out != np.float32(ND)

print('valid cells  mine=%d  ref=%d  (set diff: mine-only=%d ref-only=%d)' %
      (mine_valid.sum(), ref_valid.sum(),
       int((mine_valid & ~ref_valid).sum()), int((~mine_valid & ref_valid).sum())))

both = mine_valid & ref_valid
diff = tan_out[both].astype(np.float64) - ref[both].astype(np.float64)
ad = np.abs(diff)
mag = np.abs(ref[both].astype(np.float64))
print('on %d common cells: max|diff|=%.3e  median|diff|=%.3e  99.99pct=%.3e' %
      (both.sum(), ad.max(), np.median(ad), np.percentile(ad, 99.99)))
print('relative (|diff|/(|ref|+1e-6)): median=%.2e  max=%.2e' %
      (np.median(ad / (mag + 1e-6)), (ad / (mag + 1e-6)).max()))

# Downstream is a threshold at 0.1 (find_valleys). Count cells where mine vs ref
# land on opposite sides of 0.1 -- those are the only ones that can change valleys.
for thr in (0.1,):
    flip = both & ((tan_out > thr) != (ref > thr))
    near = both & (np.abs(ref.astype(np.float64) - thr) < 1e-3)
    print('threshold %.2f: cells flipping side = %d ; ref within 1e-3 of thr = %d'
          % (thr, int(flip.sum()), int(near.sum())))

np.save('/tmp/dreich_tancurv_mine.npy', tan_out)
print('saved /tmp/dreich_tancurv_mine.npy')
