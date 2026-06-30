"""Head-to-head: divides+chi-z vs curvature-DrEICH (Clubb's method, our faithful
reimplementation) against the field heads, SAME harness, density-controlled.

Divides numbers come from the site_experiment npz (T sweep). Curvature is swept
over its source threshold here. Both reported at their field-density-matched
operating point, in-hull, recall@30/@50 + precision + median.

NOTE on "Clubb's recall": the curvature path is our DrEICH port, validated to
reproduce Clubb's C++ head LOCATIONS (dreich-port-status), so it is a fair proxy
for Clubb's METHOD. It is NOT Clubb et al. (2014)'s published field-head numbers
-- that literature comparison needs their paper.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 compare_methods.py <mbr|indian>
"""
import sys
import numpy as np
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from rivernetworkx import dreich as D
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
cfg = SITES[key]
TSRC_LIST = (50, 100, 200, 400, 800)       # curvature source threshold (cells @1 m)
buf, hullbuf = 400, 50.0

z, west, north, res, nd = read_flt(cfg['flt'])
z = np.where(z == nd, np.nan, z)
field = load_field_heads(cfg['key'])
nr, ncc = z.shape
cr = ((north - field[:, 1]) / res).astype(int); cc = ((field[:, 0] - west) / res).astype(int)
ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < ncc); field = field[ins]
R0, R1 = max(0, cr[ins].min() - buf), min(nr, cr[ins].max() + buf)
C0, C1 = max(0, cc[ins].min() - buf), min(ncc, cc[ins].max() + buf)
zc = z[R0:R1, C0:C1]; X0, Y0 = west + C0 * res, north - R0 * res
dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32)

hull = ConvexHull(field); cen = field.mean(0); hp = field[hull.vertices]
hb = cen + (hp - cen) * (1.0 + hullbuf / np.linalg.norm(hp - cen, axis=1, keepdims=True))
tri = Delaunay(hb); in_hull = lambda xy: tri.find_simplex(xy) >= 0
hull_km2 = ConvexHull(hb).volume / 1e6
field_in = field[in_hull(field)]; nf = len(field_in); field_dens = nf / hull_km2
print('%s: %d field heads in hull, %.1f heads/km2' % (cfg['name'], nf, field_dens))


def metrics(heads):
    hxy = np.array([[X0 + (c + 0.5) * res, Y0 - (r + 0.5) * res] for (r, c) in heads]) \
        if len(heads) else np.empty((0, 2))
    hin = hxy[in_hull(hxy)] if len(hxy) else hxy
    if not (len(hin) and nf):
        return len(hin), len(hin) / hull_km2, 0., 0., float('nan'), float('nan')
    d_ref = cKDTree(hin).query(field_in)[0]; d_h = cKDTree(field_in).query(hin)[0]
    return (len(hin), len(hin) / hull_km2, float(np.mean(d_ref <= 30)),
            float(np.mean(d_ref <= 50)), float(np.mean(d_h <= 50)), float(np.median(d_ref)))


print('\nCURVATURE-DrEICH (Clubb method, our port) -- source-threshold sweep:')
print(' Tsrc  nHd  dens | recall@30 recall@50 prec@50 median_m')
crows = []
for ts in TSRC_LIST:
    heads = D.extract_channel_heads(dem, nodata=nd, cellsize=float(res),
                                    valleys='curvature', threshold=ts)
    n, dn, r30, r50, pr, med = metrics(heads)
    crows.append((dn, r30, r50, pr, med))
    print('  %4d %4d %5.1f | %.2f      %.2f      %.2f    %.1f' % (ts, n, dn, r30, r50, pr, med), flush=True)

# interpolate curvature metrics to field density (if reachable); else report its best
A = np.array(crows); o = np.argsort(A[:, 0])
def interp(col):
    return float(np.interp(field_dens, A[o, 0], A[o, col]))
cv_reach = A[:, 0].min() <= field_dens <= A[:, 0].max()

d = np.load('/tmp/site_%s.npz' % key, allow_pickle=True); m = d['matched']  # divides (r30,r50,prc,med)

print('\n==== density-matched comparison @ %.1f heads/km2 ====' % field_dens)
print('             recall@30  recall@50  prec@50  median_m')
print('  DIVIDES      %.2f       %.2f       %.2f     %.1f' % (m[0], m[1], m[2], m[3]))
if cv_reach:
    print('  CURVATURE    %.2f       %.2f       %.2f     %.1f'
          % (interp(1), interp(2), interp(3), interp(4)))
else:
    best = max(crows, key=lambda r: r[2])   # curvature can't reach field density; report its best recall@50
    print('  CURVATURE    %.2f       %.2f       %.2f     %.1f   (max density %.1f < field; best-recall row)'
          % (best[1], best[2], best[3], best[4], A[:, 0].max()))
