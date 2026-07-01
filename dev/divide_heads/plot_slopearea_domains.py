"""Process-domain diagnostic in SLOPE-AREA space (the area-based lens, not chi).

Read the full hilltop->downstream profile for field-head valleys and look for the
THREE transitions Andy named, each with its own operator:
  (1) hillslope->colluvial (valley head)  = slope-area PEAK/rollover (S rises with A
      on the hillslope limb, then turns over; Montgomery & Foufoula-Georgiou 1993);
  (2) colluvial->fluvial (channel head)    = CONCAVITY (theta) change in the declining
      limb -- colluvial/debris-flow low theta, fluvial theta~0.4-0.6 (Stock &
      Dietrich 2003);
  (3) lateral fluvial input                 = discrete STEP in area (tributary junction).

Per valley: log-log slope vs area (top), and area vs downstream distance (bottom,
to see tributary steps). Green line = field head area. Slope smoothed by regression
of elev on distance over a moving window. Just VISUALIZE -- no detector yet.
Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 plot_slopearea_domains.py <mbr|feather_s>
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx.dreich import _find_farthest_upslope
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
DOWN_TO = 3e5                                             # follow flow down to this area (m2)
WIN = 15.0                                                # slope regression window (m along flow)
cfg = SITES[key]
z, west, north, res, nd = read_flt(cfg['flt']); z = np.where(z == nd, np.nan, z)
field = load_field_heads(cfg['key'], cfg.get('nfilter'))
nr, nc = z.shape
cr = ((north - field[:, 1]) / res).astype(int); cc = ((field[:, 0] - west) / res).astype(int)
ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < nc)
R0, R1 = max(0, cr[ins].min() - 400), min(nr, cr[ins].max() + 400)
C0, C1 = max(0, cc[ins].min() - 400), min(nc, cc[ins].max() + 400)
zc = z[R0:R1, C0:C1]; X0, Y0 = west + C0 * res, north - R0 * res
dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32); field = field[ins]
print('routing %s...' % cfg['name'], flush=True)
filled = D.fill(dem, nd, cellsize=res); fi = D.build_flowinfo(filled, nd, res)
D.contributing_area(fi); D.build_svector(fi); D.distance_from_outlet(fi)
ncon, recv, ro, co, fd = fi['ncontrib'], fi['recv'], fi['row_of'], fi['col_of'], fi['fd']
pxa = res * res
zf = filled.ravel()
chan = ncon * pxa >= 2000; cn = np.where(chan)[0]
cx = X0 + (co[chan] + .5) * res; cy = Y0 - (ro[chan] + .5) * res
dsnap, cidx = cKDTree(np.column_stack([cx, cy])).query(field)
good = [k for k in range(len(field)) if dsnap[k] < 20][:6]


def profile(fhc):
    """hilltop -> downstream node sequence through fhc; return (area, dist, elev)."""
    hilltop = _find_farthest_upslope(fi, fhc)
    seq = [hilltop]
    cur = hilltop
    while ncon[cur] * pxa < DOWN_TO:
        nx = int(recv[cur])
        if nx == cur:
            break
        seq.append(nx); cur = nx
    seq = np.array(seq)
    return ncon[seq] * pxa, fd[seq], zf[seq], seq


def smooth_slope(dist, elev):
    """Local channel slope by regression of elev on distance over +/-WIN/2 window."""
    s = np.full(len(dist), np.nan)
    for i in range(len(dist)):
        m = np.abs(dist - dist[i]) <= WIN / 2
        if m.sum() >= 3:
            A = np.vstack([dist[m], np.ones(m.sum())]).T
            s[i] = abs(np.linalg.lstsq(A, elev[m], rcond=None)[0][0])
    return s


fig, axes = plt.subplots(2, 3, figsize=(17, 8))
for ax, w in zip(axes.T, good):
    a1, a2 = ax
    fhc = int(cn[cidx[w]]); fa = ncon[fhc] * pxa
    area, dist, elev, seq = profile(fhc)
    slope = smooth_slope(dist, elev)
    ddown = dist.max() - dist                             # 0 at hilltop, grows downstream
    ok = np.isfinite(slope) & (slope > 0) & (area > 0)
    a1.loglog(area[ok], slope[ok], 'k.', ms=3)
    a1.axvline(fa, color='g', lw=2, label='field head %.0f m$^2$' % fa)
    a1.set_xlabel('drainage area (m$^2$)'); a1.set_ylabel('slope'); a1.legend(fontsize=7)
    a1.set_title(cfg['name'], fontsize=8)
    a2.plot(ddown, area, 'b-'); a2.axhline(fa, color='g', lw=1.5)
    a2.set_xlabel('distance downstream (m)'); a2.set_ylabel('area (m$^2$)'); a2.set_yscale('log')
plt.suptitle('%s: slope-area (top) & area-distance (bottom) -- hillslope|colluvial|fluvial + tributary steps' % cfg['name'])
plt.tight_layout(); plt.savefig('dev/divide_heads/figures/slopearea_%s.png' % key, dpi=125)
print('wrote dev/divide_heads/figures/slopearea_%s.png' % key)
