"""Coarse-resolution sweep to pin VALLEY-STRUCTURE loss directly (Andy: test the
45-90 m prediction head-on, not via the 2-cell gauge).

Two losses, separated, vs DEM resolution:
  * VALLEY structure  -- positional self-consistency of the first-order valleys
    (T-sources) against the 1 m valleys, in Clubb's hull:
       recall    = fraction of 1 m valleys with a coarse valley within VTOL
                   (drops when real valleys disappear)
       precision = fraction of coarse valleys near a 1 m valley
                   (craters when the threshold floors and the area floods with
                    2-cell 'sources' -- the artifact we must not mistake for signal)
       n_valleys = in-hull count (explodes under flooding)
  * CHANNEL-HEAD location -- field-head recall@50 in hull (the chi-z deliverable)

The valley scale is held physically constant (T_cells = T_M2/res^2, floored at 2);
the floor (where 10000 m^2 < 2 cells) is flagged. min_segment_length is held at a
fixed physical length, so the production downstream-extension fix is active.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 valley_structure_sweep.py
"""
import numpy as np
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from rivernetworkx import dreich as D
import divides_lib as L

ALG = '/tmp/dreich_algorithm'; NR, NC = 7107, 6266; ND = -9999.0
WEST, NORTH = L.WEST, L.NORTH
T_M2 = 10000.0            # valley scale
SEG_M = 10.0             # chi-z min segment length, physical
VTOL = 100.0            # valley-match tolerance (~ valley radius; spacing ~235 m)
HTOL = 50.0            # field-head match tolerance

z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
z1 = np.where(z1 == ND, np.nan, z1)
clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
P = 400
R0, R1, C0, C1 = cr.min() - P, cr.max() + P, cc.min() - P, cc.max() + P
z1 = z1[R0:R1, C0:C1]
X0, Y0 = WEST + C0, NORTH - R0

hull = ConvexHull(clubb); cen = clubb.mean(0); hp = clubb[hull.vertices]
hb = cen + (hp - cen) * (1.0 + 50.0 / np.linalg.norm(hp - cen, axis=1, keepdims=True))
tri = Delaunay(hb)
in_hull = lambda xy: tri.find_simplex(xy) >= 0
clubb_in = clubb[in_hull(clubb)]


def block_mean(z, n):
    nr2, nc2 = z.shape[0] // n, z.shape[1] // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def run(dem, res):
    T = max(2, int(round(T_M2 / (res * res))))
    floored = T_M2 / (res * res) < 2
    seg = max(3, int(round(SEG_M / res)))
    heads, fi = D.extract_channel_heads(dem, nodata=ND, cellsize=float(res),
                                        valleys='divides', threshold=T,
                                        min_segment_length=seg, return_flowinfo=True)
    src = D.get_sources(fi, T)
    vxy = np.array([[X0 + (fi['col_of'][s] + 0.5) * res, Y0 - (fi['row_of'][s] + 0.5) * res]
                    for s in src]) if len(src) else np.empty((0, 2))
    hxy = np.array([[X0 + (c + 0.5) * res, Y0 - (r + 0.5) * res] for (r, c) in heads]) \
        if len(heads) else np.empty((0, 2))
    return vxy, hxy, T, floored


RES = [1, 20, 30, 40, 45, 50, 60, 70, 80, 90]
print('clip %dx%d; %d Clubb heads (%d in hull); valley scale %g m2 (~%.0f m wide), spacing ~235 m'
      % (z1.shape[0], z1.shape[1], len(clubb), len(clubb_in), T_M2, np.sqrt(T_M2)))
print('VTOL %g m (valley match), HTOL %g m (head match)\n' % (VTOL, HTOL))
print('res Tcl flr | nVal_inhull | VALLEY rec/prec | HEAD recall@50 | cells/valley')
v1 = None; rows = []
for res in RES:
    agg = block_mean(z1, res) if res > 1 else z1
    dem = np.where(np.isfinite(agg), agg, ND).astype(np.float32)
    vxy, hxy, T, floored = run(dem, res)
    v_in = vxy[in_hull(vxy)] if len(vxy) else vxy
    h_in = hxy[in_hull(hxy)] if len(hxy) else hxy
    if res == 1:
        v1 = v_in.copy()
    # valley structure self-consistency vs 1 m valleys, at a tight + loose tol
    def rp(tol):
        if not (len(v1) and len(v_in)):
            return 0.0, 0.0
        return (float(np.mean(cKDTree(v_in).query(v1)[0] <= tol)),
                float(np.mean(cKDTree(v1).query(v_in)[0] <= tol)))
    v_rec, v_prec = rp(VTOL); v_rec50, v_prec50 = rp(50.0)
    # channel-head recall vs Clubb field heads
    h_rec = float(np.mean(cKDTree(h_in).query(clubb_in)[0] <= HTOL)) if len(h_in) else 0.0
    cells = np.sqrt(T_M2) / res
    rows.append((res, len(v_in), v_rec, v_prec, h_rec, cells, v_rec50))
    print('%3d %3d %3s | %5d | t100 %.2f/%.2f  t50 %.2f/%.2f | %.2f | %.1f'
          % (res, T, 'YES' if floored else '.', len(v_in), v_rec, v_prec,
             v_rec50, v_prec50, h_rec, cells), flush=True)

import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os
rr = [x[0] for x in rows]
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.5, 8), sharex=True)
ax.fill_between(rr, [x[6] for x in rows], [x[2] for x in rows], color='C0', alpha=0.18)
ax.plot(rr, [x[2] for x in rows], 'o-', color='C0', label='valley recall (tol 100 m)')
ax.plot(rr, [x[6] for x in rows], 'o:', color='C0', alpha=0.7, label='valley recall (tol 50 m)')
ax.plot(rr, [x[4] for x in rows], 's-', color='C3', label='channel-head recall@50 (field)')
ax.axhline(0.5, color='grey', lw=0.6, ls=':')
ax.set_ylabel('recall / precision'); ax.set_ylim(0, 1.05); ax.legend(fontsize=8)
ax.set_title('MBR: valley structure vs channel-head location, by resolution')
ax2.plot(rr, [x[1] for x in rows], 'o-', color='C2', label='valleys in hull (count)')
ax2.set_xlabel('DEM resolution (m)'); ax2.set_ylabel('count'); ax2.legend(fontsize=8)
ax2.set_title('valley count: explosion = threshold-floor flooding, not real valleys')
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'valley_structure_sweep.png')
fig.savefig(out, dpi=130, bbox_inches='tight'); print('\nsaved', out)
