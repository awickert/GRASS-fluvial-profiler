import time
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import openpyxl
import rivernetworkx as rnx
from rivernetworkx.grass_io import _read_raster, sample_raster, assemble_records
from grass.pygrass.vector import VectorTopo

# 1. Clubb Mid Bailey Run ground-truth heads (UTM 17N E/N)
wb = openpyxl.load_workbook('/tmp/clubb_channel_heads.xlsx', data_only=True, read_only=True)
ws = wb['Sheet1']
hE, hN = [], []
for row in ws.iter_rows(min_row=3, values_only=True):
    if row[0] and 'Bailey' in str(row[0]):
        hN.append(float(row[2])); hE.append(float(row[3]))
hE, hN = np.array(hE), np.array(hN)
print('MBR ground-truth heads: %d' % len(hE))

# 2. rasters over the current (MBR window) region
acc, bounds = _read_raster('acc')
dem, _ = _read_raster('dem')
print('region array %s' % (acc.shape,))

# 3. drainage area at each mapped head. The GPS points are a few m off the 1 m
# flowlines, so raw sampling lands on hillslope cells (A~1). Snap to the channel
# the head sits on: take the max |accumulation| in a small neighborhood.
def snap_max_accum(arr, E, N, bnds, radius_cells):
    nr, nc = arr.shape
    out = []
    for e, n in zip(E, N):
        col = int(np.floor((e - bnds['west']) / bnds['ewres']))
        row = int(np.floor((bnds['north'] - n) / bnds['nsres']))
        r0, r1 = max(0, row - radius_cells), min(nr, row + radius_cells + 1)
        c0, c1 = max(0, col - radius_cells), min(nc, col + radius_cells + 1)
        sub = np.abs(arr[r0:r1, c0:c1]); sub = sub[np.isfinite(sub)]
        out.append(np.nan if sub.size == 0 else sub.max())
    return np.array(out)

rawA = np.abs(sample_raster(acc, hE, hN, **bounds))
for label, hA in (('raw (at point)', rawA),
                  ('snap r=5m', snap_max_accum(acc, hE, hN, bounds, 5)),
                  ('snap r=10m', snap_max_accum(acc, hE, hN, bounds, 10))):
    hA = hA[np.isfinite(hA)]
    ps = [np.percentile(hA, p) for p in (25, 50, 75)]
    print('head A %-16s n=%d  p25/50/75 = %.0f / %.0f / %.0f m2  (log10 med=%.2f)'
          % (label, len(hA), ps[0], ps[1], ps[2], np.log10(ps[1])))
headA = snap_max_accum(acc, hE, hN, bounds, 5)
headA = headA[np.isfinite(headA)]

# 4. dense-network geometry (one pass; drop coincident vertices)
t = time.time()
vt = VectorTopo('streams_dense'); vt.open('r')
geom = {}
for ln in vt.viter('lines'):
    if ln.cat is None:
        continue
    en = ln.to_array()
    x, y = en[:, 0], en[:, 1]
    keep = np.concatenate(([True], (np.diff(x) != 0) | (np.diff(y) != 0)))
    if keep.sum() >= 2:
        geom[ln.cat] = (x[keep], y[keep])
vt.close()
print('read %d segment geometries in %.1fs' % (len(geom), time.time() - t))

# 5. sample z, A along vertices; slope-area (topology not needed -> dummy tostream)
cats = list(geom)
zsmp = {c: sample_raster(dem, gx, gy, **bounds) for c, (gx, gy) in geom.items()}
Asmp = {c: np.abs(sample_raster(acc, gx, gy, **bounds)) for c, (gx, gy) in geom.items()}
recs = assemble_records(cats, {c: 0 for c in cats}, geom, z=zsmp, A=Asmp)
t = time.time()
logA, logS = rnx.slope_area(recs, window=30.0, log=True)
m = np.isfinite(logA) & np.isfinite(logS)
logA, logS = logA[m], logS[m]
print('slope_area: %d finite pts in %.1fs; logA %.2f..%.2f' %
      (len(logA), time.time() - t, logA.min(), logA.max()))

# 6. plot cloud + binned median, with the Clubb head-area band overlaid
plt.figure(figsize=(7.5, 6.5))
plt.scatter(logA, logS, s=2, alpha=0.05, color='0.6')
bins = np.linspace(logA.min(), logA.max(), 30)
ix = np.digitize(logA, bins); mA, mS = [], []
for i in range(1, len(bins)):
    sel = ix == i
    if sel.sum() > 40:
        mA.append(np.median(logA[sel])); mS.append(np.median(logS[sel]))
plt.plot(mA, mS, 'r.-', lw=2, label='binned median')
hlo, hmed, hhi = (np.log10(np.percentile(headA, 25)),
                  np.log10(np.median(headA)),
                  np.log10(np.percentile(headA, 75)))
plt.axvspan(hlo, hhi, color='C0', alpha=0.15, label='Clubb head A (IQR)')
plt.axvline(hmed, color='C0', ls='--', lw=1.6, label='Clubb head A (median)')
plt.xlabel('log10 drainage area [cells; 1 cell = 1 m2]')
plt.ylabel('log10 channel slope (-dz/ds, flowline-smoothed 30 m)')
plt.title('Mid Bailey Run, OH (extract thr=100): slope-area vs Clubb heads')
plt.legend(loc='lower left'); plt.tight_layout()
plt.savefig('/tmp/midbaileyrun_slope_area.png', dpi=115)
print('saved /tmp/midbaileyrun_slope_area.png')

# --- broken-stick fit: flat hillslope + power-law fluvial, continuous knot ---
# Drop S~0 flat-valley artifacts (logS < -6) and the high-A trunk noise.
fm = (logA >= 2.1) & (logA <= 5.0) & (logS > -6)
la, ls = logA[fm], logS[fm]
best = None
for k in np.linspace(2.4, 4.2, 181):
    X = np.column_stack([np.ones_like(la), np.maximum(0.0, la - k)])
    coef, *_ = np.linalg.lstsq(X, ls, rcond=None)
    sse = float(np.sum((ls - X @ coef) ** 2))
    if best is None or sse < best[0]:
        best = (sse, k, coef)
sse, k, (h, sl) = best
print('BROKEN-STICK pooled: head A* = 10^%.2f = %.0f m2 ; theta=%.2f ; hillslope logS=%.2f'
      % (k, 10 ** k, -sl, h))
print('Clubb head A median (snap r=5m) = 10^%.2f = %.0f m2'
      % (np.log10(np.median(headA)), np.median(headA)))

plt.figure(figsize=(7.5, 6.5))
plt.scatter(logA, logS, s=2, alpha=0.05, color='0.6')
plt.plot(mA, mS, 'k.-', lw=1, alpha=0.5, label='binned median')
xx = np.linspace(la.min(), la.max(), 200)
yy = h + np.minimum(0.0, sl * (xx - k)) * 0 + sl * np.maximum(0.0, xx - k)
plt.plot(xx, yy, 'r-', lw=2.5, label='broken-stick fit')
plt.axvline(k, color='r', ls=':', lw=1.5, label='predicted head A*=%.0f m2' % 10 ** k)
plt.axvspan(hlo, hhi, color='C0', alpha=0.15, label='Clubb head A (IQR)')
plt.axvline(hmed, color='C0', ls='--', lw=1.6, label='Clubb head A (median)')
plt.xlabel('log10 drainage area [cells; 1 cell = 1 m2]')
plt.ylabel('log10 channel slope (-dz/ds, flowline-smoothed 30 m)')
plt.title('Mid Bailey Run: broken-stick S-A fit vs Clubb heads')
plt.legend(loc='lower left', fontsize=8); plt.tight_layout()
plt.savefig('/tmp/midbaileyrun_break.png', dpi=115)
print('saved /tmp/midbaileyrun_break.png')
