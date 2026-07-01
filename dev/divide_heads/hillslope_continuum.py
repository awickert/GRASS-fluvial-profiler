"""Keystone test of the threshold-hillslope continuum (Roering et al. 1999, 2007).

Prediction: the hillslope-> colluvial transition is a CONTINUUM set by nondimensional
erosion rate E* = E L_H/(K S_c).
  - low E* (gentle, soil-mantled MBR): slope rises ~linearly with distance from the
    crest to a PEAK at the valley head (parabolic hilltop), S << S_c.
  - high E* (steep Feather): slope SATURATES fast to a near-constant threshold S_c;
    no rising limb / peak, a plateau instead.

Measures, per site:
  * slope S vs distance-from-crest  (crest = flow-tree leaf, ncon==1; distance via
    farthest-upslope hilltop and flow distance fd)  -- binned median +/- IQR;
  * S vs drainage area over the hillslope+colluvial range;
  * hilltop curvature C_HT = Laplacian(z) at ridge cells (convex => negative;
    ~ -E/K, Roering 2007);
  * S_c estimate = high percentile of hillslope gradient; hillslope length L_H.
DEM lightly Gaussian-smoothed (sigma px) before slope/curvature to tame 1 m / TIN noise.
Run: PYTHONPATH=<repo>:<repo>/dev/divide_heads /usr/bin/python3 hillslope_continuum.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from rivernetworkx import dreich as D
from rivernetworkx.dreich import _find_farthest_upslope
from site_experiment import read_flt, load_field_heads, SITES

SIGMA = 2.0                                               # smoothing (px) for slope/curvature
A_MAX = 5000.0                                            # sample cells up to this area (hillslope+colluvial)
NSAMP = 25000
DMAX = 160.0                                              # distance-from-crest axis limit (m)


def analyse(key):
    cfg = SITES[key]
    z, west, north, res, nd = read_flt(cfg['flt']); z = np.where(z == nd, np.nan, z)
    field = load_field_heads(cfg['key'], cfg.get('nfilter'))
    nr, nc = z.shape
    cr = ((north - field[:, 1]) / res).astype(int); cc = ((field[:, 0] - west) / res).astype(int)
    ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < nc)
    R0, R1 = max(0, cr[ins].min() - 300), min(nr, cr[ins].max() + 300)
    C0, C1 = max(0, cc[ins].min() - 300), min(nc, cc[ins].max() + 300)
    zc = z[R0:R1, C0:C1]; H, W = zc.shape
    dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32)
    print('routing %s...' % cfg['name'], flush=True)
    filled = D.fill(dem, nd, cellsize=res); fi = D.build_flowinfo(filled, nd, res)
    D.contributing_area(fi); D.build_svector(fi); D.distance_from_outlet(fi)
    ncon, fd = fi['ncontrib'], fi['fd']; pxa = res * res
    # slope + curvature on lightly smoothed RAW topography (faithful hillslope form)
    zs = zc.copy(); m = ~np.isfinite(zs); zs[m] = np.nanmean(zc)
    zs = gaussian_filter(zs, SIGMA)
    gy, gx = np.gradient(zs, res); slope = np.sqrt(gx ** 2 + gy ** 2)      # tan(theta)
    lap = (np.gradient(np.gradient(zs, res, axis=0), res, axis=0)
           + np.gradient(np.gradient(zs, res, axis=1), res, axis=1))       # Laplacian curvature
    sflat = slope.ravel(); lapflat = lap.ravel()
    # ridge cells = flow-tree leaves (no donors) => ncon==1, and finite topo
    ridge = (ncon == 1) & np.isfinite(zc.ravel())
    C_HT = np.nanmedian(lapflat[ridge])
    # sample hillslope+colluvial cells; distance-from-crest = fd[hilltop]-fd[cell]
    cand = np.where((ncon * pxa <= A_MAX) & (ncon >= 1) & np.isfinite(zc.ravel()))[0]
    rs = np.random.default_rng(0).choice(cand, min(NSAMP, cand.size), replace=False)
    dist = np.empty(rs.size); Sv = sflat[rs]; Av = ncon[rs] * pxa
    for i, cnode in enumerate(rs):
        ht = _find_farthest_upslope(fi, int(cnode))
        dist[i] = fd[ht] - fd[cnode]
    ok = np.isfinite(dist) & np.isfinite(Sv) & (dist >= 0)
    dist, Sv, Av = dist[ok], Sv[ok], Av[ok]
    # Physical S_c = the value slope SATURATES to on the mid/lower hillslope (robust
    # median away from the crest), NOT a high percentile (that is cliff/TIN noise).
    Sc = np.nanmedian(Sv[(dist > 25) & (Av < 800)])
    # hillslope length: distance from crest at which area first ~ valley scale (per crest sample)
    LH = np.nanmedian(dist[(Av > 300) & (Av < 800)]) if np.any((Av > 300) & (Av < 800)) else np.nan
    print('  %-22s C_HT(med)=%.4f 1/m  Sc(plateau)=%.2f (%.0f deg)  L_H~%.0f m  n=%d'
          % (cfg['name'], C_HT, Sc, np.degrees(np.arctan(Sc)), LH, dist.size))
    return cfg['name'], dist, Sv, Av, C_HT, Sc, lapflat[ridge]


def binmed(x, y, edges):
    idx = np.digitize(x, edges); lo, md, hi, ctr = [], [], [], []
    for b in range(1, len(edges)):
        yy = y[idx == b]
        if yy.size >= 20:
            lo.append(np.percentile(yy, 25)); md.append(np.median(yy)); hi.append(np.percentile(yy, 75))
            ctr.append(0.5 * (edges[b - 1] + edges[b]))
    return np.array(ctr), np.array(lo), np.array(md), np.array(hi)


res = [analyse(k) for k in ('mbr', 'feather_s')]
edges = np.linspace(0, DMAX, 33)
fig, ax = plt.subplots(2, 2, figsize=(15, 11))
cols = {'Mid Bailey Run, OH': 'tab:blue', 'Feather River South, CA': 'tab:red'}
for name, dist, Sv, Av, C_HT, Sc, ridgelap in res:
    c = cols[name]
    ctr, lo, md, hi = binmed(dist, Sv, edges)
    ax[0, 0].plot(ctr, md, '-', color=c, lw=2, label=name)
    ax[0, 0].fill_between(ctr, lo, hi, color=c, alpha=0.15)
    ax[0, 0].axhline(Sc, color=c, ls=':', lw=1)
    aedges = np.logspace(0, np.log10(A_MAX), 25)
    ac, alo, amd, ahi = binmed(Av, Sv, aedges)
    ax[0, 1].plot(ac, amd, '-', color=c, lw=2, label=name)
    ax[0, 1].fill_between(ac, alo, ahi, color=c, alpha=0.15)
    ax[1, 1].hist(ridgelap, bins=np.linspace(-0.4, 0.4, 60), color=c, alpha=0.5,
                  density=True, label='%s  C_HT med=%.3f' % (name.split(',')[0], C_HT))
    ax[1, 1].axvline(C_HT, color=c, ls='--')
ax[0, 0].set_xlabel('distance from crest (m)'); ax[0, 0].set_ylabel('slope tan($\\theta$)')
ax[0, 0].set_ylim(0, 1.5); ax[0, 0].set_title('slope vs distance from crest (median, IQR; dotted = $S_c$ plateau)')
ax[0, 0].legend(fontsize=9)
ax[0, 1].set_xscale('log'); ax[0, 1].set_xlabel('drainage area (m$^2$)'); ax[0, 1].set_ylabel('slope')
ax[0, 1].set_ylim(0, 1.5); ax[0, 1].set_title('slope-area, hillslope+colluvial range'); ax[0, 1].legend(fontsize=9)
ax[1, 1].set_xlabel('ridge-cell curvature $\\nabla^2 z$ (1/m); convex = negative')
ax[1, 1].set_title('hilltop curvature $C_{HT}\\approx -E/K$'); ax[1, 1].legend(fontsize=8)
ax[1, 0].axis('off')
ax[1, 0].text(0.02, 0.5,
    'Prediction: MBR (low E*) slope RISES ~linearly to a peak, S<<S_c;\n'
    'Feather (high E*) slope SATURATES fast to near-constant S_c (plateau),\n'
    'no rising limb. => hillslope->colluvial break is a CONTINUUM;\n'
    'redefine it as ONSET OF NEGATIVE S-A SCALING (peak OR plateau-end).',
    fontsize=11, va='center')
plt.suptitle('Threshold-hillslope continuum test: MBR (gentle) vs Feather S (steep)')
plt.tight_layout(); plt.savefig('dev/divide_heads/figures/hillslope_continuum.png', dpi=125)
print('wrote dev/divide_heads/figures/hillslope_continuum.png')
