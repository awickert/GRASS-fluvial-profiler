import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from rivernetworkx.core import moving_average
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command, region as _region

run_command('g.region', raster='dem', quiet=True)
reg = _region(); north = float(reg['n']); west = float(reg['w']); south = float(reg['s']); east = float(reg['e'])
nsres = float(reg['nsres']); ewres = float(reg['ewres'])
run_command('r.watershed', flags='s', elevation='dem', accumulation='acc_p', drainage='dir_p', overwrite=True, quiet=True)
run_command('r.stream.extract', elevation='dem', accumulation='acc_p', threshold=30,
            stream_raster='srast_p', direction='sdir_p', d8cut=0, overwrite=True, quiet=True)
run_command('r.relief', input='dem', output='relief_p', overwrite=True, quiet=True)
acc, b = _read_raster('acc_p'); dem, _ = _read_raster('dem'); sr, _ = _read_raster('srast_p')
sr_dir, _ = _read_raster('dir_p'); relief, _ = _read_raster('relief_p')
nr, nc = dem.shape; res = ewres; A = np.abs(acc)

move = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}
chan = np.isfinite(sr) & (sr > 0)
rs, cs = np.where(chan); dvals = np.where(np.isfinite(sr_dir[rs, cs]), sr_dir[rs, cs], 0).astype(int)
inflow = np.zeros(dem.shape, int)
for d, (dr, dc) in move.items():
    sel = dvals == d; tr = rs[sel] + dr; tc = cs[sel] + dc
    inb = (tr >= 0) & (tr < nr) & (tc >= 0) & (tc < nc)
    np.add.at(inflow, (tr[inb], tc[inb]), 1)
sources = np.argwhere(chan & (inflow == 0))
print('channel heads (sources): %d' % len(sources))

def trace_down(r0, c0):
    path = [(r0, c0)]; r, c = r0, c0
    for _ in range(50000):
        d = sr_dir[r, c]
        if not np.isfinite(d) or int(d) <= 0 or int(d) not in move:
            break
        dr, dc = move[int(d)]; r2, c2 = r + dr, c + dc
        if r2 < 0 or r2 >= nr or c2 < 0 or c2 >= nc:
            break
        r, c = r2, c2; path.append((r, c))
    return path

W = max(3, int(round(40.0 / res)))    # confirm-decline window (~40 m)
def peak_transition(rr, cc):
    if len(rr) < W + 3:
        return None
    z = dem[rr, cc]; Av = A[rr, cc]
    step = np.hypot(np.diff(rr.astype(float)), np.diff(cc.astype(float))) * res
    dist = np.concatenate([[0], np.cumsum(step)])
    zs = moving_average(dist, z, 40.0)
    S = -np.gradient(zs, dist)
    Ss = moving_average(dist, S, 40.0)
    runmax = -np.inf; argmax = None; conf = 0
    for i in range(len(Ss)):
        if not (np.isfinite(Ss[i]) and np.isfinite(Av[i]) and Av[i] > 0 and Ss[i] > 0):
            continue
        if Ss[i] > runmax:
            runmax = Ss[i]; argmax = i; conf = 0
        elif Ss[i] < 0.8 * runmax:
            conf += 1
            if conf >= W and argmax is not None:
                return argmax, Av[argmax], Ss[argmax]
        else:
            conf = 0
    return None

txy = []  # (x, y, A, S)
for r0, c0 in sources:
    p = trace_down(int(r0), int(c0))
    rr = np.array([q[0] for q in p]); cc = np.array([q[1] for q in p])
    pk = peak_transition(rr, cc)
    if pk is None:
        continue
    i, Apk, Spk = pk
    x = west + (cc[i] + 0.5) * ewres; y = north - (rr[i] + 0.5) * nsres
    txy.append((x, y, float(Apk), float(Spk)))
txy = np.array(txy)
print('transitions found: %d' % len(txy))
areas = txy[:, 2] * ewres * nsres
print('transition area (m^2): p10=%.0f  median=%.0f  p90=%.0f' %
      (np.percentile(areas, 10), np.median(areas), np.percentile(areas, 90)))
print('transition A (cells): p10=%.0f  median=%.0f  p90=%.0f' %
      (np.percentile(txy[:, 2], 10), np.median(txy[:, 2]), np.percentile(txy[:, 2], 90)))
xmid = (west + east) / 2
print('WEST (x<%.0f): %d   EAST: %d' % (xmid, (txy[:, 0] < xmid).sum(), (txy[:, 0] >= xmid).sum()))

fig, ax = plt.subplots(figsize=(11, 10))
ax.imshow(relief, cmap='gray', extent=[west, east, south, north], origin='upper')
ax.scatter(txy[:, 0], txy[:, 1], s=14, c='red', edgecolors='k', linewidths=0.3, label='peak transition')
ax.axvline(xmid, color='cyan', lw=0.7, ls=':')
ax.set_title('S-A peak transitions (down from each head) — %d heads, %d transitions' % (len(sources), len(txy)))
ax.legend(loc='upper right'); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); plt.savefig('/tmp/tremp_peak.png', dpi=120)
print('saved /tmp/tremp_peak.png')
run_command('g.remove', type='raster', name='acc_p,dir_p,srast_p,sdir_p,relief_p', flags='f', quiet=True)
