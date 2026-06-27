import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import openpyxl
import rivernetworkx as rnx
from rivernetworkx.grass_io import _read_raster

# Clubb heads
wb = openpyxl.load_workbook('/tmp/clubb_channel_heads.xlsx', data_only=True, read_only=True)
ws = wb['Sheet1']
hE, hN = [], []
for row in ws.iter_rows(min_row=3, values_only=True):
    if row[0] and 'Bailey' in str(row[0]):
        hN.append(float(row[2])); hE.append(float(row[3]))
hE, hN = np.array(hE), np.array(hN)

acc, b = _read_raster('acc')
# A* from the module's own pipeline (consistency)
recs = rnx.read_stream_segments('streams_dense', elevation='dem',
                                accumulation='acc', assume_complete=True)
logA, logS = rnx.slope_area(recs, window=30.0, log=True)
keep = (logS >= np.log10(1e-4)) & (logA <= np.log10(1e5))
fit = rnx.fit_sa_break(logA[keep], logS[keep])
A_star = fit['A_star']

# Clubb head area: snap to nearest extracted stream cell
strm, _ = _read_raster('streams_rast')
sr, sc = np.where(np.isfinite(strm) & (strm != 0))
sx = b['west'] + (sc + 0.5) * b['ewres']
sy = b['north'] - (sr + 0.5) * b['nsres']
sA = np.abs(acc[sr, sc])
A_clubb = np.empty(len(hE))
for i, (e, n) in enumerate(zip(hE, hN)):
    j = int(np.argmin(np.hypot(sx - e, sy - n)))
    A_clubb[i] = sA[j]

ratio = A_clubb / A_star
print('A* (fluvial initiation) = %.0f m2' % A_star)
print('Clubb head area / A*  -- is the relationship consistent?')
for p in (10, 25, 50, 75, 90):
    print('  p%-2d ratio = %.2f' % (p, np.percentile(ratio, p)))
print('  log10(ratio): mean=%.2f  std=%.2f  (std is the spread = (in)consistency)'
      % (np.mean(np.log10(ratio)), np.std(np.log10(ratio))))
print('  fraction of Clubb heads upstream of fluvial onset (A_clubb < A*): %.0f%%'
      % (100.0 * np.mean(A_clubb < A_star)))

plt.figure(figsize=(7, 5))
plt.hist(np.log10(ratio), bins=16, color='0.6', edgecolor='k')
plt.axvline(0, color='r', lw=2, label='A_clubb = A* (fluvial onset)')
plt.axvline(np.median(np.log10(ratio)), color='b', ls='--', lw=2,
            label='median ratio = %.2f' % np.median(ratio))
plt.xlabel('log10( Clubb head area / A* )'); plt.ylabel('count')
plt.title('Mid Bailey Run: Clubb heads relative to fluvial initiation A*')
plt.legend(); plt.tight_layout()
plt.savefig('/tmp/mbr_relationship.png', dpi=120)
print('saved /tmp/mbr_relationship.png')
