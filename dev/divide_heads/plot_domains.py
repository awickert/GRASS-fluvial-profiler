"""Process-domain diagnostic: read the FULL hilltop->downstream profile (no T
truncation) for field-head valleys and look for the TWO transitions --
hillslope->colluvial (valley head, chi rollover) and colluvial->fluvial (channel
head, ksn step / fluvial onset). Does MBR show a thin colluvial package and
Feather a thick one, and does the LOWER break sit on the field head at both?
Panels per valley: chi-elevation (top) and slope=d(elev)/d(chi) i.e. ksn (bottom).
Green = field head. Run for <mbr|feather_s>."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
from rivernetworkx.dreich import _find_farthest_upslope, _build_channel_chi
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
DOWN_TO = 3e5                                              # follow flow down to this area (m2)
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
ncon, recv, ro, co = fi['ncontrib'], fi['recv'], fi['row_of'], fi['col_of']; pxa = res * res
chan = ncon * pxa >= 2000; cn = np.where(chan)[0]
cx = X0 + (co[chan] + .5) * res; cy = Y0 - (ro[chan] + .5) * res
dsnap, cidx = cKDTree(np.column_stack([cx, cy])).query(field)
good = [k for k in range(len(field)) if dsnap[k] < 20][:6]

fig, axes = plt.subplots(2, 3, figsize=(17, 8))
for ax, w in zip(axes.T, good):                           # each column = one valley (2 stacked panels)
    a1, a2 = ax
    fhc = int(cn[cidx[w]]); fa = ncon[fhc] * pxa
    hilltop = _find_farthest_upslope(fi, fhc)
    end = fhc                                             # follow the full flow path far downstream
    while ncon[end] * pxa < DOWN_TO:
        nx = int(recv[end])
        if nx == end:
            break
        end = nx
    nsq, chi, elev = _build_channel_chi(fi, filled, hilltop, end, 1000.0, 0.525)
    chi = np.array(chi); elev = np.array(elev); areas = np.array([ncon[x] * pxa for x in nsq])
    slope = np.gradient(elev, chi)
    fh = int(np.argmin(np.abs(areas - fa)))
    a1.plot(chi, elev, 'k-'); a1.axvline(chi[fh], color='g', lw=2, label='field head %.0f m$^2$' % fa)
    a1.set_ylabel('elev (m)'); a1.legend(fontsize=7); a1.set_title('%s (chi 0=down, right=hilltop)' % cfg['name'], fontsize=8)
    a2.semilogy(chi, np.clip(slope, 1e-4, None), 'm-'); a2.axvline(chi[fh], color='g', lw=2)
    a2.set_ylabel('slope d(elev)/d(chi)=ksn'); a2.set_xlabel('$\\chi$')
plt.suptitle('%s: full hilltop->downstream profiles -- looking for hillslope|colluvial|fluvial (two breaks)' % cfg['name'])
plt.tight_layout(); plt.savefig('dev/divide_heads/figures/domains_%s.png' % key, dpi=125)
print('wrote dev/divide_heads/figures/domains_%s.png' % key)
