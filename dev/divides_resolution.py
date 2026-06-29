"""The test that actually matters: does divides->chi-z heads survive coarsening
where curvature-DrEICH collapsed to 0 heads (recall 0.49 @2m -> 0 @>=10m)?

Mean-aggregate the 1 m MBR DEM to 2/3/5/8/10/12 m, run the full divides->heads
pipeline at each (fill, route, contributing area, T-scale valley sources, chi-z
head per valley, dedup) at a FIXED PHYSICAL valley scale T_m2. Compare each
resolution's heads to the 1 m divides heads (self-consistency) on the same
tol = max(50 m, 2*res) grid the curvature resolution study used.

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_resolution.py
"""
import time
import numpy as np
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
import divides_lib as L

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
WEST, NORTH = 398767.32685636, 4369656.1523925
RESLIST = [1, 2, 3, 5, 8, 10, 12]
T_M2 = 10000.0            # physical valley scale (matched head-area at 1 m)

T0 = time.time()
def log(m): print('[%6.1fs] %s' % (time.time() - T0, m), flush=True)


def block_mean(z, n):
    nr, nc = z.shape
    nr2, nc2 = nr // n, nc // n
    return np.nanmean(z[:nr2 * n, :nc2 * n].reshape(nr2, n, nc2, n), axis=(1, 3))


def heads_xy(fi, nodes, res):
    r = fi['row_of'][np.asarray(nodes, dtype=np.int64)]
    c = fi['col_of'][np.asarray(nodes, dtype=np.int64)]
    return np.column_stack([WEST + (c + 0.5) * res, NORTH - (r + 0.5) * res])


def divide_heads(filled, res):
    """Full divides->chi-z head pipeline on a filled DEM at given resolution."""
    fi = D.build_flowinfo(filled, ND, float(res))
    D.contributing_area(fi)
    D.build_svector(fi)
    D.distance_from_outlet(fi)
    T = max(2, int(round(T_M2 / (res * res))))
    tsrc = D.get_sources(fi, T)
    heads = []
    for s in tsrc:
        h, _s, _n = L.head_and_score(fi, filled, int(s), min_segment_length=10)
        heads.append(h)
    final = L.dedup_furthest_upstream(fi, filled, heads)
    final = [h for h in final if h != 0]
    return fi, final


# 1 m: reuse the cached filled DEM + reference curvature heads
C = L.load_cache()
fi1, filled1 = C['fi'], C['filled']
ref_xy = heads_xy(fi1, C['ref_heads'], 1)        # 634 curvature heads (1 m)
log('1 m divides heads (T=%g m2)' % T_M2)
_, heads1 = divide_heads(filled1, 1)
xy1 = heads_xy(fi1, heads1, 1)
log('1 m: %d divide heads' % len(heads1))

z1 = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
z1 = np.where(z1 == ND, np.nan, z1)

print('\nres  heads  selfR  selfP   curvR  curvP    (tol=max(50,2res))')
base = cKDTree(xy1)
refck = cKDTree(ref_xy)
for res in RESLIST:
    if res == 1:
        xy = xy1
    else:
        agg = block_mean(z1, res)
        filled = D.fill(np.where(np.isfinite(agg), agg, ND).astype(np.float32),
                        ND, 0.0001, float(res))
        fi, heads = divide_heads(filled, res)
        xy = heads_xy(fi, heads, res)
    tol = max(50.0, 2.0 * res)
    if len(xy) == 0:
        print('%3d  %5d   ----   ----    ----  ----' % (res, 0)); continue
    ck = cKDTree(xy)
    dself, _ = ck.query(xy1)                 # nearest res-head to each 1m-head
    dself_p, _ = base.query(xy)              # nearest 1m-head to each res-head
    dcurv, _ = ck.query(ref_xy)              # nearest res-head to each curv head
    dcurv_p, _ = refck.query(xy)
    print('%3d  %5d   %.2f   %.2f    %.2f  %.2f'
          % (res, len(xy), np.mean(dself <= tol), np.mean(dself_p <= tol),
             np.mean(dcurv <= tol), np.mean(dcurv_p <= tol)), flush=True)
    log('res=%d done' % res)
log('done')
