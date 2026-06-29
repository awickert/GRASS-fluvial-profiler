"""Set the valley scale T against Clubb HONESTLY: restrict to the area Clubb
surveyed (convex hull of the 53 field heads, lightly buffered) and report
precision as well as recall inside it. Basin-wide Clubb recall is confounded by
head density (finer T trivially wins on a sparse target); within the surveyed
area, a too-fine T shows up as low precision (many heads, no field head nearby).

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_clubb_extent.py
"""
import time
import numpy as np
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from rivernetworkx import dreich as D
import divides_lib as L

T0 = time.time()
def log(m): print('[%6.1fs] %s' % (time.time() - T0, m), flush=True)

C = L.load_cache()
fi, filled = C['fi'], C['filled']
clubb = L.load_clubb()

# Clubb survey extent: convex hull of the field heads, buffered outward by BUF m
BUF = 50.0
hull = ConvexHull(clubb)
cen = clubb.mean(axis=0)
hull_pts = clubb[hull.vertices]
hull_buf = cen + (hull_pts - cen) * (1.0 + BUF / np.linalg.norm(hull_pts - cen, axis=1, keepdims=True))
tri = Delaunay(hull_buf)
def inside(xy):
    return tri.find_simplex(xy) >= 0
log('Clubb hull area ~ %.2f km2; %d field heads' % (hull.volume / 1e6, len(clubb)))


def heads_for_T(T):
    heads = []
    for s in D.get_sources(fi, T):
        h, _s, _n = L.head_and_score(fi, filled, int(s), min_segment_length=10)
        heads.append(h)
    return [h for h in L.dedup_furthest_upstream(fi, filled, heads) if h != 0]


# curvature reference too, for context
methods = [('curvature', None)] + [('T=%d' % T, T) for T in (3000, 5000, 8000, 10000)]
print('\nWITHIN Clubb survey hull (buffer %gm):' % BUF)
print('method      heads_in   recall  precision   med_dist   (tol=30m)')
ckc = cKDTree(clubb)
for name, T in methods:
    nodes = C['ref_heads'] if T is None else heads_for_T(T)
    xy = L.nodes_xy(fi, nodes)
    xin = xy[inside(xy)]
    if len(xin) == 0:
        print('%-10s  %7d   ----     ----' % (name, 0)); continue
    mt = cKDTree(xin)
    d_field, _ = mt.query(clubb)           # nearest in-hull head to each field head
    d_head, _ = ckc.query(xin)             # nearest field head to each in-hull head
    tol = 30.0
    recall = np.mean(d_field <= tol)       # frac of 53 field heads matched
    prec = np.mean(d_head <= tol)          # frac of in-hull heads that are real
    print('%-10s  %7d   %.2f      %.2f       %4.0fm'
          % (name, len(xin), recall, prec, np.median(d_field)), flush=True)
log('done')
