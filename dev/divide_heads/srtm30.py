"""Is 30 m (SRTM) a hard wall for method=divides, or just the wrong valley scale?
The 30 m collapse is chi-z profile starvation (valleys still resolved, profiles
too short). A COARSER valley scale -> bigger valleys -> longer profiles should
recover heads -- but bigger valleys give coarser heads that may miss the fine
field heads. Sweep the valley scale at 30 m (and 20 m for reference) and watch
both: heads-recovered (nHd vs nVal) and field-head recall.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 srtm30.py
"""
import numpy as np
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from rivernetworkx import dreich as D
import divides_lib as L

ALG = '/tmp/dreich_algorithm'; NR, NC = 7107, 6266; ND = -9999.0
WEST, NORTH = L.WEST, L.NORTH
z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
z1 = np.where(z1 == ND, np.nan, z1)
clubb = L.load_clubb()
cr = (NORTH - clubb[:, 1]).astype(int); cc = (clubb[:, 0] - WEST).astype(int)
P = 400; R0, C0 = cr.min() - P, cc.min() - P
z1 = z1[R0:cr.max() + P, C0:cc.max() + P]
X0, Y0 = WEST + C0, NORTH - R0
hull = ConvexHull(clubb); cen = clubb.mean(0); hp = clubb[hull.vertices]
hb = cen + (hp - cen) * (1.0 + 50.0 / np.linalg.norm(hp - cen, axis=1, keepdims=True))
tri = Delaunay(hb); in_hull = lambda xy: tri.find_simplex(xy) >= 0
clubb_in = clubb[in_hull(clubb)]


def block_mean(z, n):
    a, b = z.shape[0] // n, z.shape[1] // n
    return np.nanmean(z[:a * n, :b * n].reshape(a, n, b, n), axis=(1, 3))


print('res  T_m2   T_cell cells/val nVal nHd | recall@50 prec@50  (in hull, vs field)')
for res in (20, 30):
    dem = np.where(np.isfinite(block_mean(z1, res)), block_mean(z1, res), ND).astype(np.float32)
    for T_m2 in (10000, 20000, 40000, 80000, 160000):
        T = max(2, int(round(T_m2 / (res * res))))
        seg = max(3, int(round(10.0 / res)))
        heads, fi = D.extract_channel_heads(dem, nodata=ND, cellsize=float(res),
                                            valleys='divides', threshold=T,
                                            min_segment_length=seg, return_flowinfo=True)
        nval = len(D.get_sources(fi, T))
        xy = np.array([[X0 + (c + 0.5) * res, Y0 - (r + 0.5) * res] for (r, c) in heads])
        mh = xy[in_hull(xy)] if len(xy) else xy
        if len(mh):
            dref, _ = cKDTree(mh).query(clubb_in); dh, _ = cKDTree(clubb_in).query(mh)
            rec, prec = np.mean(dref <= 50), np.mean(dh <= 50)
        else:
            rec = prec = 0.0
        print('%3d %6d %5d   %5.1f   %4d %3d | %.2f      %.2f'
              % (res, T_m2, T, np.sqrt(T_m2) / res, nval, len(mh), rec, prec), flush=True)
    print()
