"""One MBR profile in detail: chi-elevation (top) and its slope d(elev)/d(chi)
= channel steepness proxy (bottom), to answer (1) why a downstream channel->channel
break exists, and (2) what the biggest rollover ABOVE the channel head is.
Marks: green=field head, blue=source-end pick, red=junction-end pick, purple=chi-z
test stat, black=max chi-curvature (biggest rollover)."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx.dreich import _find_farthest_upslope, _build_channel_chi, _calculate_channel_head, _reg32
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
which = int(sys.argv[2]) if len(sys.argv) > 2 else 0
T = 5000; MS = 10
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
filled = D.fill(dem, nd, cellsize=res); fi = D.build_flowinfo(filled, nd, res)
D.contributing_area(fi); D.junction_network(fi, D.get_sources(fi, T)); D.build_svector(fi); D.distance_from_outlet(fi)
ncon, recv, ro, co, JIdx = fi['ncontrib'], fi['recv'], fi['row_of'], fi['col_of'], fi['JIdx']; pxa = res * res
chan = ncon * pxa >= 2000; cn = np.where(chan)[0]
cx = X0 + (co[chan] + .5) * res; cy = Y0 - (ro[chan] + .5) * res
dsnap, cidx = cKDTree(np.column_stack([cx, cy])).query(field)
good = [k for k in range(len(field)) if dsnap[k] < 20]
i = good[which]
fhc = int(cn[cidx[i]]); fa = ncon[fhc] * pxa
s = fhc
while ncon[s] * pxa < T:
    nx = int(recv[s])
    if nx == s: break
    s = nx
hilltop = _find_farthest_upslope(fi, s); jn = s
while True:
    nx = int(recv[jn])
    if nx == jn: break
    jn = nx
    if JIdx[jn] != -9999: break
nsq, chi, elev = _build_channel_chi(fi, filled, hilltop, jn, 1000., .525)
chi = np.array(chi); elev = np.array(elev); areas = np.array([ncon[x] * pxa for x in nsq])
n = len(chi); tst = np.full(n, np.nan)
for h in range(MS, n - MS + 1):
    r2, _ = _reg32(chi[h:], elev[h:]); _, dw = _reg32(chi[:h], elev[:h]); tst[h] = r2 - (dw - 2) / 2
ns_s, ch_s, el_s = _build_channel_chi(fi, filled, hilltop, s, 1000., .525)
hs = _calculate_channel_head(ns_s, ch_s, el_s, MS); hj = _calculate_channel_head(nsq, chi, elev, MS)
slope = np.gradient(elev, chi)                       # d elev / d chi = channel steepness
curv = np.gradient(slope, chi)                       # chi-space curvature (rollover)
chi_of = {int(nn): float(c) for nn, c in zip(nsq, chi)}
roll = int(np.nanargmax(np.abs(curv[MS:-MS])) + MS)  # biggest rollover
mk = lambda nnode: chi_of.get(int(nnode), np.nan)

fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax.plot(chi, elev, 'k-'); ax.set_ylabel('elevation (m)')
for nnode, col, lab in [(fhc, 'g', 'field head %.0f m2' % fa), (hs, 'b', 'source-end %.0f' % (ncon[hs] * pxa)),
                        (hj, 'r', 'junction-end %.0f' % (ncon[hj] * pxa))]:
    ax.axvline(mk(nnode), color=col, lw=2, label=lab)
ax.axvline(chi[roll], color='k', ls='--', lw=1.5, label='max rollover %.0f m2' % areas[roll])
ax.axvline(mk(s), color='orange', ls=':', label='source %.0f' % (ncon[s] * pxa)); ax.legend(fontsize=8)
ax.set_title('%s profile #%d (chi 0 = junction, right = hilltop)' % (cfg['name'], which))
a3 = ax.twinx(); a3.plot(chi, tst, color='purple', alpha=0.4); a3.set_ylabel('chi-z test stat', color='purple')
ax2.plot(chi, slope, 'm-'); ax2.set_ylabel('slope d(elev)/d(chi) = ksn'); ax2.set_xlabel('chi (upstream ->)')
for nnode, col in [(fhc, 'g'), (hs, 'b'), (hj, 'r')]:
    ax2.axvline(mk(nnode), color=col, lw=2)
ax2.axvline(chi[roll], color='k', ls='--', lw=1.5)
plt.tight_layout(); plt.savefig('dev/divide_heads/figures/one_profile_%s.png' % key, dpi=130)
print('field %.0f | source %.0f | src-end %.0f | jct-end %.0f | max-rollover %.0f m2'
      % (fa, ncon[s] * pxa, ncon[hs] * pxa, ncon[hj] * pxa, areas[roll]))
print('wrote dev/divide_heads/figures/one_profile_%s.png' % key)
