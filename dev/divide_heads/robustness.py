"""Resolution robustness of the channel-head methods against the FIELD heads.

Headline question (Andy): as the DEM coarsens, does method=divides keep recovering
Clubb (2014)'s field-mapped channel heads -- where curvature-DrEICH collapses?
Primary metric: field-head recall/precision within Clubb's survey hull
(density-controlled), vs resolution, for divides vs curvature. Self-consistency
(coarse heads vs 1 m heads) is a secondary diagnostic only.

Validity guards: the valley scale is held PHYSICALLY constant
(threshold_cells = T_m2/res^2); min_segment_length is held at a fixed physical
length (so a coarse profile is not starved into a false collapse); and the
head->absolute-coordinate mapping carries the clip offset + res explicitly.

Uses the productionised method (extract_channel_heads(valleys=...)).
Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 robustness.py [quick]
"""
import sys
import numpy as np
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from rivernetworkx import dreich as D
import divides_lib as L

ALG = '/tmp/dreich_algorithm'; NR, NC = 7107, 6266; ND = -9999.0
WEST, NORTH = L.WEST, L.NORTH
T_M2 = 10000.0           # divides: valley scale (matched head count/area at 1 m)
T_SRC_M2 = 100.0         # curvature: fine source threshold (DrEICH default 100 @1 m)
SEG_M = 10.0             # chi-z min segment length, physical (10 nodes @1 m)
TOLS = (30.0, 50.0)      # field-head match tolerances (m)

z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
z1 = np.where(z1 == ND, np.nan, z1)
clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
P = 400
R0, R1, C0, C1 = cr.min() - P, cr.max() + P, cc.min() - P, cc.max() + P
z1 = z1[R0:R1, C0:C1]
X0, Y0 = WEST + C0, NORTH - R0                    # clip origin in absolute coords

# Clubb survey hull (+50 m buffer) for the density-controlled, in-hull comparison
hull = ConvexHull(clubb); cen = clubb.mean(0); hp = clubb[hull.vertices]
hb = cen + (hp - cen) * (1.0 + 50.0 / np.linalg.norm(hp - cen, axis=1, keepdims=True))
tri = Delaunay(hb)
in_hull = lambda xy: tri.find_simplex(xy) >= 0
clubb_in = clubb[in_hull(clubb)]


def block_mean(z, n):
    nr2, nc2 = z.shape[0] // n, z.shape[1] // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def heads_xy(heads, res):
    return np.array([[X0 + (c + 0.5) * res, Y0 - (r + 0.5) * res] for (r, c) in heads])


def run_method(dem, res, valleys):
    T = max(2, int(round((T_M2 if valleys == 'divides' else T_SRC_M2) / (res * res))))
    seg = max(3, int(round(SEG_M / res)))
    win = max(7.0, 1.6 * res)
    out = D.extract_channel_heads(dem, nodata=ND, cellsize=float(res),
                                  valleys=valleys, threshold=T,
                                  min_segment_length=seg, window_radius=win,
                                  return_flowinfo=True)
    heads, fi = out
    nsrc = len(D.get_sources(fi, T)) if valleys == 'divides' else -1   # # first-order valleys
    return heads_xy(heads, res), T, nsrc


def metrics(mxy, ref):
    """In-hull recall/precision of mxy against reference points ref, at TOLS."""
    mh = mxy[in_hull(mxy)] if len(mxy) else mxy
    out = {}
    if len(mh) == 0 or len(ref) == 0:
        return {t: (0.0, 0.0) for t in TOLS}, len(mh)
    d_ref, _ = cKDTree(mh).query(ref)
    d_h, _ = cKDTree(ref).query(mh)
    for t in TOLS:
        out[t] = (float(np.mean(d_ref <= t)), float(np.mean(d_h <= t)))
    return out, len(mh)


RES = ([1, 30] if len(sys.argv) > 1 and sys.argv[1] == 'quick'
       else [1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
# first-order valley spacing on MBR ~235 m; print the resolution/valley-scale ratios
print('clip %dx%d; %d Clubb heads (%d in hull); valley scale=%g m2 (~%.0f m wide)' % (
    z1.shape[0], z1.shape[1], len(clubb), len(clubb_in), T_M2, np.sqrt(T_M2)))
print('\nres Tcell nVal nHd | DIV recall@30/50 prec@50 | CRV @50 | selfcon | cells/valley')
heads1 = None
rows = []
for res in RES:
    agg = block_mean(z1, res) if res > 1 else z1
    dem = np.where(np.isfinite(agg), agg, ND).astype(np.float32)
    dxy, T, nsrc = run_method(dem, res, 'divides')
    cxy, _, _ = run_method(dem, res, 'curvature')
    dm, dn = metrics(dxy, clubb_in); cm, cn = metrics(cxy, clubb_in)
    if res == 1:
        heads1 = dxy[in_hull(dxy)]
    sc = metrics(dxy, heads1)[0][50.0][0] if heads1 is not None and len(heads1) else float('nan')
    rows.append((res, dm, cm, nsrc, len(dxy)))
    # cells across a first-order valley (~sqrt(valley area) / res), a Nyquist gauge
    cells_across = np.sqrt(T_M2) / res
    print('%3d %5d %4d %3d | %.2f/%.2f %.2f | %.2f | %.2f | %.1f'
          % (res, T, nsrc, len(dxy), dm[30.0][0], dm[50.0][0], dm[50.0][1],
             cm[50.0][0], sc, cells_across), flush=True)

# two panels: field-head recall (top) and the failure mechanism (valley & head
# counts, bottom) vs resolution
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os
rr = [x[0] for x in rows]
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.5, 8), sharex=True)
for t, ls in ((30.0, '-'), (50.0, '--')):
    ax.plot(rr, [x[1][t][0] for x in rows], 'o' + ls, color='C0', label='divides @%dm' % t)
    ax.plot(rr, [x[2][t][0] for x in rows], 's' + ls, color='C3', label='curvature @%dm' % t)
ax.set_ylabel('field-head recall (in hull)'); ax.set_ylim(0, 1); ax.legend(fontsize=8)
ax.set_title('MBR: field-head recovery vs resolution')
ax2.plot(rr, [x[3] for x in rows], 'o-', color='C2', label='first-order valleys found')
ax2.plot(rr, [x[4] for x in rows], 's-', color='C4', label='divide heads (total)')
ax2.set_xlabel('DEM resolution (m)'); ax2.set_ylabel('count'); ax2.legend(fontsize=8)
ax2.set_title('Failure mechanism: valleys resolved vs heads located')
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'robustness.png')
fig.savefig(out, dpi=130, bbox_inches='tight'); print('saved', out)
