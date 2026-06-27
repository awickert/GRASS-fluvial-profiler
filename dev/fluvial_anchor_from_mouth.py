import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from rivernetworkx.core import moving_average
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command, region as _region

run_command('g.region', raster='dem', quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='acc_f', drainage='dir_f', overwrite=True, quiet=True)
run_command('r.slope.aspect', elevation='dem', slope='slp_f', format='percent', overwrite=True, quiet=True)
run_command('r.stream.extract', elevation='dem', accumulation='acc_f', threshold=50,
            stream_raster='srast_f', direction='sdir_f', d8cut=0, overwrite=True, quiet=True)
acc, b = _read_raster('acc_f'); dem, _ = _read_raster('dem'); slp, _ = _read_raster('slp_f'); sr, _ = _read_raster('srast_f')
nr, nc = dem.shape; res = b['ewres']
A = np.abs(acc); S = slp / 100.0
# floodplain elevation cut (for marking bluff vs floodplain; NOT for trimming)
ev = dem.ravel(); sv = S.ravel(); fn = np.isfinite(ev) & np.isfinite(sv)
qs = np.percentile(ev[fn], np.arange(2, 60, 2)); cut = qs[0]
for q in qs:
    band = fn & (ev <= q) & (ev > q - (qs[1] - qs[0]))
    if band.sum() and (sv[band] < 1e-3).mean() < 0.4:
        cut = q; break
print('floodplain cut = %.1f m' % cut)

move = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}
def trace_down(r0, c0):
    path = [(r0, c0)]; r, c = r0, c0
    for _ in range(200000):
        d = sr_dir[r, c]
        if not np.isfinite(d) or int(d) <= 0 or int(d) not in move:
            break
        dr, dc = move[int(d)]; r2, c2 = r + dr, c + dc
        if r2 < 0 or r2 >= nr or c2 < 0 or c2 >= nc:
            break
        r, c = r2, c2; path.append((r, c))
    return path
sr_dir, _ = _read_raster('dir_f')

# pick the longest flowline among the highest-elevation channel heads
chan = np.isfinite(sr) & (sr > 0)
ce = np.where(chan, dem, -np.inf)
heads = np.dstack(np.unravel_index(np.argsort(ce.ravel())[::-1][:40], dem.shape))[0]
best = []
for r0, c0 in heads:
    p = trace_down(int(r0), int(c0))
    if len(p) > len(best):
        best = p
path = best
rr = np.array([p[0] for p in path]); cc = np.array([p[1] for p in path])  # head -> mouth
print('flowline length: %d cells (%.0f m)' % (len(path), len(path) * res))

# along-flowline arrays (head -> mouth)
zf = dem[rr, cc]; Af = A[rr, cc]; elevf = dem[rr, cc]
step = np.hypot(np.diff(rr.astype(float)), np.diff(cc.astype(float))) * res
dist = np.concatenate([[0], np.cumsum(step)])
zs = moving_average(dist, zf, 30.0)
Sf = -np.gradient(zs, dist); Sf[~np.isfinite(Sf)] = np.nan
ok = np.isfinite(Af) & np.isfinite(Sf) & (Af > 0) & (Sf > 0)
la = np.log10(Af); ls = np.log10(np.where(Sf > 0, Sf, np.nan))
bluff = elevf > cut

# fit fluvial S~A^-theta on the bluff fluvial reach (upper half by A), walk up to departure
m = ok & bluff
order = np.argsort(la[m])[::-1]                 # high A -> low A (mouth -> head, bluff)
idxs = np.where(m)[0][order]
laB = la[idxs]; lsB = ls[idxs]
nfit = max(8, int(0.45 * len(laB)))             # developed fluvial reach = highest-A ~45%
th, intc = np.polyfit(laB[:nfit], lsB[:nfit], 1)
resid = lsB - (intc + th * laB)
tol = 2.5 * np.median(np.abs(resid[:nfit] - np.median(resid[:nfit])))
trans = None
for i in range(nfit, len(laB)):
    if abs(resid[i]) > tol and (i + 1 >= len(laB) or abs(resid[i + 1]) > tol):
        trans = i; break
print('fluvial fit theta=%.2f ; departure (transition) at logA=%.2f (A=%.0f cells)'
      % (-th, laB[trans] if trans else np.nan, 10 ** laB[trans] if trans else -1))

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].scatter(la[ok & ~bluff], ls[ok & ~bluff], s=8, c='0.7', label='floodplain (kept, not fit)')
ax[0].scatter(la[ok & bluff], ls[ok & bluff], s=10, c='#1f77ff', label='bluff channel')
xx = np.linspace(np.nanmin(laB), np.nanmax(laB), 50)
ax[0].plot(xx, intc + th * xx, 'k--', lw=2, label='fluvial fit S~A^-%.2f' % (-th))
if trans is not None:
    ax[0].scatter([laB[trans]], [lsB[trans]], s=120, c='red', marker='*', zorder=6, label='transition (departure)')
ax[0].set_xlabel('log10 A (cells)'); ax[0].set_ylabel('log10 S'); ax[0].legend(fontsize=8)
ax[0].set_title('walk up from mouth: fit fluvial, find upslope departure')
ax[1].plot(dist, zf, '0.6', lw=1); ax[1].plot(dist, zs, 'b', lw=1.5, label='smoothed z')
if trans is not None:
    dtr = dist[idxs[trans]]
    ax[1].axvline(dtr, color='red', lw=2, label='transition')
ax[1].set_xlabel('distance from head (m)'); ax[1].set_ylabel('elevation (m)'); ax[1].legend(fontsize=8)
ax[1].set_title('long profile (head -> mouth)')
plt.tight_layout(); plt.savefig('/tmp/tremp_fluvanchor.png', dpi=120)
print('saved /tmp/tremp_fluvanchor.png')
run_command('g.remove', type='raster', name='acc_f,dir_f,slp_f,srast_f,sdir_f', flags='f', quiet=True)
