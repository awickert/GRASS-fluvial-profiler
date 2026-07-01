"""Plot the METHOD'S OWN chi-elevation profiles (no preconditioning on the field
head). For each field head we find its valley the way the method does -- walk up
the main stem to the divide source (area>=T), take that source's hilltop->tributary-
junction profile -- then OVERLAY where the field head falls and where the detector
picks under source-end vs junction-end. This shows, honestly, whether the head
sits at an upper hillslope->channel break and whether a downstream break steals
the global-max pick.
Run: PYTHONPATH=<repo>:<repo>/dev:<repo>/dev/divide_heads /usr/bin/python3 plot_profiles.py <mbr|feather_s>
"""
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
dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32)
field = field[ins]
print('routing %s...' % cfg['name'], flush=True)
filled = D.fill(dem, nd, cellsize=res); fi = D.build_flowinfo(filled, nd, res)
D.contributing_area(fi)
D.junction_network(fi, D.get_sources(fi, T)); D.build_svector(fi); D.distance_from_outlet(fi)
ncon, recv, ro, co, NI, JIdx = fi['ncontrib'], fi['recv'], fi['row_of'], fi['col_of'], fi['NodeIndex'], fi['JIdx']
delta, ds = fi['delta'], fi['donorstack']; pxa = res * res


def source_of(fhc):                              # walk DOWN to first cell with area>=T (the divide source)
    cur = fhc
    while ncon[cur] * pxa < T:
        nx = int(recv[cur])
        if nx == cur:
            break
        cur = nx
    return cur


def profile(a, b):                               # nodeseq a(top)->b(bottom) + chi + elev
    return _build_channel_chi(fi, filled, a, b, 1000.0, 0.525)


# snap field heads to channel cells, keep well-snapped ones
chan = ncon * pxa >= 2000; cn = np.where(chan)[0]
cx = X0 + (co[chan] + .5) * res; cy = Y0 - (ro[chan] + .5) * res
dsnap, cidx = cKDTree(np.column_stack([cx, cy])).query(field)
sel = [i for i in range(len(field)) if dsnap[i] < 20][:6]

fig, axes = plt.subplots(2, 3, figsize=(16, 8))
for ax, i in zip(axes.ravel(), sel):
    fhc = int(cn[cidx[i]]); fa = ncon[fhc] * pxa
    s = source_of(fhc)                           # the method's source for this valley
    hilltop = _find_farthest_upslope(fi, s)
    jn = s                                        # walk source -> tributary junction
    while True:
        nx = int(recv[jn])
        if nx == jn:
            break
        jn = nx
        if JIdx[jn] != -9999:
            break
    nsq, chi, elev = profile(hilltop, jn)         # full method profile
    chi_of = {int(n): float(c) for n, c in zip(nsq, chi)}
    n = len(chi); tst = np.full(n, np.nan)
    for h in range(MS, n - MS + 1):
        r2, _ = _reg32(chi[h:], elev[h:]); _, dw = _reg32(chi[:h], elev[:h]); tst[h] = r2 - (dw - 2) / 2
    # detector picks: source-end (hilltop->source) vs junction-end (hilltop->junction)
    ns_s, ch_s, el_s = profile(hilltop, s)
    h_src = _calculate_channel_head(ns_s, ch_s, el_s, MS)
    h_jun = _calculate_channel_head(nsq, chi, elev, MS)
    ax.plot(chi, elev, '-', color='0.4', lw=1.5)
    if s in chi_of: ax.axvline(chi_of[s], color='orange', ls=':', lw=1.5, label='source (T=%d)' % T)
    ax.axvline(chi_of.get(fhc, np.nan), color='g', lw=2, label='field head %.0f m$^2$' % fa)
    if h_src in chi_of: ax.plot(chi_of[h_src], elev[nsq.index(h_src)] if h_src in nsq else np.nan, 'b^', ms=11, label='source-end pick', zorder=5)
    if h_jun in chi_of: ax.plot(chi_of[h_jun], elev[nsq.index(h_jun)], 'r*', ms=15, label='junction-end pick', zorder=6)
    a2 = ax.twinx(); a2.plot(chi, tst, color='purple', alpha=0.4, lw=1); a2.set_ylabel('test stat', color='purple', fontsize=8)
    ax.set_xlabel('$\\chi$ (upstream $\\rightarrow$)'); ax.set_ylabel('elev (m)'); ax.legend(fontsize=6.5, loc='upper left')
plt.suptitle('%s: method-own chi-elevation profiles (hilltop->junction). green=field head, orange=source, blue=source-end pick, red=junction-end pick' % cfg['name'], fontsize=10)
plt.tight_layout(); plt.savefig('dev/divide_heads/figures/profiles_%s.png' % key, dpi=130)
print('wrote dev/divide_heads/figures/profiles_%s.png' % key)
