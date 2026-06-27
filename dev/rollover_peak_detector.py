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

def trace_down(r0, c0, maxlen=500, maxA=6000):
    path = [(r0, c0)]; r, c = r0, c0
    for _ in range(maxlen):
        d = sr_dir[r, c]
        if not np.isfinite(d) or int(d) <= 0 or int(d) not in move:
            break
        dr, dc = move[int(d)]; r2, c2 = r + dr, c + dc
        if r2 < 0 or r2 >= nr or c2 < 0 or c2 >= nc:
            break
        r, c = r2, c2; path.append((r, c))
        if A[r, c] > maxA:
            break
    return path

def lin(x, y):
    n = len(x); sx = x.sum(); sy = y.sum(); sxx = (x * x).sum(); sxy = (x * y).sum()
    den = n * sxx - sx * sx
    if den == 0:
        return 0.0, y.mean()
    m = (n * sxy - sx * sy) / den; bq = (sy - m * sx) / n
    return m, bq

def rollover(la, ls):
    """Look up AND down: tent fit with a rising (colluvial) then falling (fluvial)
    limb. Returns (peak_logA, R, theta_down) or None."""
    o = np.argsort(la); la = la[o]; ls = ls[o]
    nb = 12; edges = np.linspace(la.min(), la.max(), nb + 1)
    idx = np.clip(np.digitize(la, edges) - 1, 0, nb - 1)
    bx, by = [], []
    for i in range(nb):
        sel = idx == i
        if sel.sum() >= 2:
            bx.append(np.median(la[sel])); by.append(np.median(ls[sel]))
    bx, by = np.array(bx), np.array(by)
    if len(bx) < 7:
        return None
    ml, bl = lin(bx, by); rss_line = float(np.sum((by - (ml * bx + bl)) ** 2))
    best = None
    for p in range(2, len(bx) - 2):
        au, bu = lin(bx[:p + 1], by[:p + 1]); ad, bd = lin(bx[p:], by[p:])
        if au > 0.05 and ad < -0.10:
            pred = np.where(np.arange(len(bx)) <= p, au * bx + bu, ad * bx + bd)
            rss = float(np.sum((by - pred) ** 2))
            if best is None or rss < best[0]:
                best = (rss, bx[p], -ad)
    if best is None:
        return None
    rss, peak_la, theta = best
    R = 1 - rss / rss_line if rss_line > 0 else 0.0
    if R < 0.5:
        return None
    return peak_la, R, theta

txy = []
for r0, c0 in sources:
    p = trace_down(int(r0), int(c0))
    if len(p) < 12:
        continue
    rr = np.array([q[0] for q in p]); cc = np.array([q[1] for q in p])
    z = dem[rr, cc]; Av = A[rr, cc]
    step = np.hypot(np.diff(rr.astype(float)), np.diff(cc.astype(float))) * res
    dist = np.concatenate([[0], np.cumsum(step)])
    zs = moving_average(dist, z, 40.0); S = -np.gradient(zs, dist); Ss = moving_average(dist, S, 40.0)
    ok = np.isfinite(Av) & np.isfinite(Ss) & (Av > 0) & (Ss > 0)
    if ok.sum() < 12:
        continue
    la = np.log10(Av[ok]); ls = np.log10(Ss[ok])
    ro = rollover(la, ls)
    if ro is None:
        continue
    peak_la, R, theta = ro
    j = np.argmin(np.abs(np.log10(np.where(Av > 0, Av, np.nan)) - peak_la))
    x = west + (cc[j] + 0.5) * ewres; y = north - (rr[j] + 0.5) * nsres
    txy.append((x, y, float(Av[j]), R, theta))
txy = np.array(txy)
print('heads with a genuine rollover: %d' % len(txy))

# dedup: snap to ~12 m grid, keep one (highest R) per cell
snap = 12.0
keys = (np.round(txy[:, 0] / snap).astype(int), np.round(txy[:, 1] / snap).astype(int))
seen = {}
for i in range(len(txy)):
    k = (keys[0][i], keys[1][i])
    if k not in seen or txy[i, 3] > txy[seen[k], 3]:
        seen[k] = i
uq = txy[sorted(seen.values())]
print('unique transitions (deduped @ %.0f m): %d' % (snap, len(uq)))
areas = uq[:, 2] * ewres * nsres
print('area (m^2): p10=%.0f median=%.0f p90=%.0f' %
      (np.percentile(areas, 10), np.median(areas), np.percentile(areas, 90)))
print('theta: median=%.2f   R: median=%.2f' % (np.median(uq[:, 4]), np.median(uq[:, 3])))
xmid = (west + east) / 2
print('WEST: %d   EAST: %d' % ((uq[:, 0] < xmid).sum(), (uq[:, 0] >= xmid).sum()))
np.save('/tmp/tremp_transitions.npy', uq)            # for instant replot

aspect = (east - west) / (north - south)
H = 9.0; fig, ax = plt.subplots(figsize=(H * aspect + 1.5, H))
ax.imshow(relief, cmap='gray', extent=[west, east, south, north], origin='upper')
ax.set_aspect('equal')
ax.scatter(uq[:, 0], uq[:, 1], s=10, c='red', edgecolors='none', alpha=0.85)
ax.set_title('Trempealeau Hills — colluvial-fluvial transitions (robust S-A rollover)\n'
             '%d transitions   median area %.0f m$^2$   $\\theta$=%.2f' %
             (len(uq), np.median(areas), np.median(uq[:, 4])), fontsize=11)
ax.set_xlabel('Easting (m)'); ax.set_ylabel('Northing (m)')
plt.tight_layout(); plt.savefig('/tmp/tremp_peak2.png', dpi=150)
print('saved /tmp/tremp_peak2.png')
run_command('g.remove', type='raster', name='acc_p,dir_p,srast_p,sdir_p,relief_p', flags='f', quiet=True)
