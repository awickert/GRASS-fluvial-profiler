"""Arm B -- resolution robustness of method=divides, per site.

Coarsens the DEM (block-mean) over a fixed grid, holding the valley scale T
PHYSICALLY constant, and measures vs resolution:
  * field-head recall@30/@50 + median (in-hull, density-controlled by fixed T)
  * valley-structure self-consistency (first-order valleys vs the 1 m valleys)
  * curvature-DrEICH field recall@50 baseline (expected to collapse early)

Three axes for the breakdown (the novel part):
  (1) absolute resolution (m);
  (2) dimensionless res/(2R) using the site median valley radius R -- so sites
      with different valley sizes are comparable and Andy's "full breakdown when
      one pixel = the whole valley-head width = 2R" is a universal res/(2R)=1 law;
  (3) PER-VALLEY normalized: each field head carries its own local 2R_i (from the
      nearest 1 m divide head's measured radius). If the aggregate breakdown is
      gradual only because valleys vary in size, plotting detection vs res/(2R_i)
      SHARPENS it to a cliff at ~1; if it stays gradual, the METHOD is graceful.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 arm_b.py <mbr|indian> [quick]
"""
import sys
import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from rivernetworkx import dreich as D
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
quick = len(sys.argv) > 2 and sys.argv[2] == 'quick'
cfg = SITES[key]
T_M2 = 5000.0                      # density-appropriate valley scale (Arm A)
SEG_M = 10.0
HTOL, VTOL = 50.0, 100.0
buf, hullbuf = 500, 50.0
RES_LIST = ([1, 30, 90, 180] if quick else
            [1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90,
             100, 120, 140, 160, 180, 200])
STR8 = np.ones((3, 3), bool)

z, west, north, res0, nd = read_flt(cfg['flt'])
z = np.where(z == nd, np.nan, z)
field = load_field_heads(cfg['key'])
nr, ncc = z.shape
cr = ((north - field[:, 1]) / res0).astype(int); cc = ((field[:, 0] - west) / res0).astype(int)
ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < ncc); field = field[ins]
R0, R1 = max(0, cr[ins].min() - buf), min(nr, cr[ins].max() + buf)
C0, C1 = max(0, cc[ins].min() - buf), min(ncc, cc[ins].max() + buf)
zc = z[R0:R1, C0:C1]; X0, Y0 = west + C0 * res0, north - R0 * res0

hull = ConvexHull(field); cen = field.mean(0); hp = field[hull.vertices]
hb = cen + (hp - cen) * (1.0 + hullbuf / np.linalg.norm(hp - cen, axis=1, keepdims=True))
tri = Delaunay(hb); in_hull = lambda xy: tri.find_simplex(xy) >= 0
field_in = field[in_hull(field)]; nf = len(field_in)
print('%s: clip %s, %d field heads in hull' % (cfg['name'], zc.shape, nf), flush=True)


def block_mean(a, n):
    if n == 1:
        return a
    h, w = a.shape[0] // n, a.shape[1] // n
    return np.nanmean(a[:h * n, :w * n].reshape(h, n, w, n), axis=(1, 3))


def heads_xy(heads, res):
    return (np.array([[X0 + (c + 0.5) * res, Y0 - (r + 0.5) * res] for (r, c) in heads])
            if len(heads) else np.empty((0, 2)))


def head_radius(fi, h):
    s = int(fi['SVectorIndex'][h]); n = int(fi['ncontrib'][h])
    U = fi['SVector'][s:s + n]; rows = fi['row_of'][U]; cols = fi['col_of'][U]
    r0, c0 = rows.min(), cols.min()
    H = rows.max() - r0 + 3; W = cols.max() - c0 + 3
    mask = np.zeros((H, W), bool); mask[rows - r0 + 1, cols - c0 + 1] = True
    locn = np.full((H, W), -1, np.int64); locn[rows - r0 + 1, cols - c0 + 1] = U
    rim = mask & ~ndimage.binary_erosion(mask, STR8, border_value=0)
    rn = locn[rim]
    d = np.hypot(fi['row_of'][rn] - fi['row_of'][h], fi['col_of'][rn] - fi['col_of'][h])
    return float(np.median(d)) * res0

# ---- 1 m reference: divide heads, valleys, and per-head radius -----------------
dem1 = np.where(np.isfinite(zc), zc, nd).astype(np.float32)
T1 = max(2, int(round(T_M2 / res0 ** 2)))
heads1, fi1 = D.extract_channel_heads(dem1, nodata=nd, cellsize=float(res0),
                                      valleys='divides', threshold=T1, return_flowinfo=True)
h1xy = heads_xy(heads1, res0)
R1arr = np.array([head_radius(fi1, fi1['NodeIndex'][r, c]) for (r, c) in heads1])
in1 = in_hull(h1xy)
v1xy = heads_xy([(int(fi1['row_of'][s]), int(fi1['col_of'][s])) for s in D.get_sources(fi1, T1)], res0)
v1_in = v1xy[in_hull(v1xy)] if len(v1xy) else v1xy
R_site = float(np.median(R1arr[in1])) if in1.any() else float('nan')
# each field head's local valley radius = nearest 1 m divide head's radius
fld_R = R1arr[cKDTree(h1xy).query(field_in)[1]] if len(h1xy) else np.full(nf, np.nan)
print('  1 m: %d divide heads (%d in hull), median valley radius R=%.0f m -> 2R=%.0f m'
      % (len(heads1), in1.sum(), R_site, 2 * R_site), flush=True)

# ---- resolution sweep ---------------------------------------------------------
rows = []; det = []          # det[res] = bool per field head (detected within HTOL)
print('  res Tcell | divR@30 divR@50 med | curvR@50 | valley rec/prec', flush=True)
for res in RES_LIST:
    agg = block_mean(zc, res)
    dem = np.where(np.isfinite(agg), agg, nd).astype(np.float32)
    T = max(2, int(round(T_M2 / (res * res))))
    seg = max(3, int(round(SEG_M / res)))
    heads, fi = D.extract_channel_heads(dem, nodata=nd, cellsize=float(res), valleys='divides',
                                        threshold=T, min_segment_length=seg, return_flowinfo=True)
    hin = (lambda h: h[in_hull(h)] if len(h) else h)(heads_xy(heads, res))
    if len(hin):
        dref = cKDTree(hin).query(field_in)[0]
        r30, r50, med = float(np.mean(dref <= 30)), float(np.mean(dref <= 50)), float(np.median(dref))
        detected = dref <= HTOL
    else:
        r30 = r50 = 0.0; med = float('nan'); detected = np.zeros(nf, bool)
    # valley structure vs 1 m
    vxy = heads_xy([(int(fi['row_of'][s]), int(fi['col_of'][s])) for s in D.get_sources(fi, T)], res)
    v_in = vxy[in_hull(vxy)] if len(vxy) else vxy
    if len(v1_in) and len(v_in):
        vrec = float(np.mean(cKDTree(v_in).query(v1_in)[0] <= VTOL))
        vprec = float(np.mean(cKDTree(v1_in).query(v_in)[0] <= VTOL))
    else:
        vrec = vprec = 0.0
    # curvature baseline (field recall@50)
    ch = D.extract_channel_heads(dem, nodata=nd, cellsize=float(res), valleys='curvature',
                                 threshold=max(2, int(round(100 / (res * res)))), min_segment_length=seg,
                                 window_radius=max(7.0, 1.6 * res))
    chin = (lambda h: h[in_hull(h)] if len(h) else h)(heads_xy(ch, res))
    cr50 = float(np.mean(cKDTree(chin).query(field_in)[0] <= 50)) if len(chin) else 0.0
    rows.append((res, T, r30, r50, med, cr50, vrec, vprec, len(hin)))
    det.append(detected)
    print('  %3d %5d | %.2f    %.2f    %4.0f | %.2f     | %.2f / %.2f'
          % (res, T, r30, r50, med if np.isfinite(med) else -1, cr50, vrec, vprec), flush=True)

rows = np.array(rows); det = np.array(det)            # det: (n_res, nf)
np.savez('/tmp/arm_b_%s.npz' % key, rows=rows, det=det, fld_R=fld_R, R_site=R_site,
         res_list=np.array(RES_LIST), site=cfg['name'], nf=nf)
# usable resolution = where divides recall@50 crosses 0.5
r50 = rows[:, 3]
passing = np.where(r50 >= 0.5)[0]
ur = RES_LIST[passing[-1]] if len(passing) else 0
vrec = rows[:, 6]
vpass = np.where(vrec >= 0.5)[0]
vur = RES_LIST[vpass[-1]] if len(vpass) else 0
print('  usable resolution: channel-head recall@50>=0.5 to ~%g m; valley structure to ~%g m  (2R=%.0f m)'
      % (ur, vur, 2 * R_site))
print('saved /tmp/arm_b_%s.npz' % key)
