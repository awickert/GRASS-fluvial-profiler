"""Standard channel-head validation experiment, identical across field sites.

For a site (DEM .flt + .hdr, and its Clubb et al. field heads), at native
resolution: sweep the valley scale T (m^2), run method=divides, and measure the
density-controlled field-head recall/precision INSIDE the survey convex hull
(+50 m), exactly as for Mid Bailey Run. The density-matched row (detected in-hull
count ~ field in-hull count) is the apples-to-apples cross-site number.

Site-agnostic: georef + dims come from the .hdr; field heads from the 3-site
spreadsheet by site name. Anchor-checks: (1) the new multi-site head loader must
reproduce the old MBR loader; (2) MBR through this harness must reproduce our
established density-controlled numbers.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 site_experiment.py <mbr|indian|all>
"""
import sys
import numpy as np
from scipy.spatial import cKDTree, ConvexHull, Delaunay
from rivernetworkx import dreich as D

XLSX = '/home/awickert/Downloads/Channel_head_coords.xlsx'
SITES = {
    'mbr':    dict(name='Mid Bailey Run, OH', flt='/tmp/dreich_algorithm/bailey_run_dem.flt', key='Bailey'),
    'indian': dict(name='Indian Creek, OH',   flt='/tmp/dreich_algorithm/indian_creek_dem.flt', key='Indian'),
    # Feather River, CA: 15 heads in two clusters split at N~4,394,000 (10 S, 5 N)
    'feather_s': dict(name='Feather River South, CA', flt='/tmp/feather_south/output.tin.tif',
                      key='Feather', nfilter=('<', 4394000.0)),
    'feather_n': dict(name='Feather River North, CA', flt='/tmp/feather_north/output.tin.tif',
                      key='Feather', nfilter=('>', 4394000.0)),
}
T_LIST = (2000, 3000, 4000, 5000, 7000, 10000, 15000, 20000, 40000)
TOLS = (30.0, 50.0)


def read_flt(path):
    """Read a DEM as (z, west, north_top, res, nodata). Accepts an ESRI .flt (+.hdr)
    or a GeoTIFF (.tif/.tiff, via rasterio) so OpenTopography DEMs ingest directly."""
    if path.lower().endswith(('.tif', '.tiff')):
        import rasterio
        with rasterio.open(path) as src:
            z = src.read(1).astype(np.float64)
            b = src.bounds
            nd = -9999.0 if src.nodata is None else float(src.nodata)
            return z, b.left, b.top, float(src.res[0]), nd
    h = {}
    for line in open(path[:-4] + '.hdr'):
        p = line.split()
        if len(p) == 2:
            h[p[0].lower()] = p[1]
    nr, nc, res = int(h['nrows']), int(h['ncols']), float(h['cellsize'])
    nd = float(h['nodata_value'])
    z = np.fromfile(path, '<f4').reshape(nr, nc).astype(np.float64)
    return z, float(h['xllcorner']), float(h['yllcorner']) + nr * res, res, nd   # west, north(top)


def load_field_heads(site_key, nfilter=None):
    """Field heads (E, N) for a site; ``nfilter=('<'|'>', N)`` splits a site by
    northing (used to separate the two Feather River clusters)."""
    import openpyxl
    wb = openpyxl.load_workbook(XLSX, data_only=True, read_only=True)
    ws = wb['Sheet1']
    EN = [(float(r[3]), float(r[2])) for r in ws.iter_rows(min_row=3, values_only=True)
          if r[0] and site_key in str(r[0])]
    arr = np.array(EN)
    if nfilter is not None and len(arr):
        op, val = nfilter
        arr = arr[arr[:, 1] < val] if op == '<' else arr[arr[:, 1] > val]
    return arr


def anchor_check_loader():
    """New multi-site MBR heads must match the old single-site loader."""
    try:
        import divides_lib as L
        old = L.load_clubb()
        new = load_field_heads('Bailey')
        d = cKDTree(new).query(old)[0]
        ok = len(old) == len(new) and d.max() < 1e-3
        print('anchor (head loader): old %d vs new %d MBR heads, max coord diff %.4g m -> %s'
              % (len(old), len(new), d.max(), 'OK' if ok else 'MISMATCH'))
    except Exception as e:
        print('anchor (head loader): skipped (%s)' % e)


def run_site(key, buf=400, hullbuf=50.0):
    cfg = SITES[key]
    z, west, north, res, nd = read_flt(cfg['flt'])
    z = np.where(z == nd, np.nan, z)
    field = load_field_heads(cfg['key'], cfg.get('nfilter'))
    nr, nc = z.shape
    cr = ((north - field[:, 1]) / res).astype(int)
    cc = ((field[:, 0] - west) / res).astype(int)
    inside = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < nc)
    print('%s: DEM %dx%d @ %g m; %d field heads (%d inside DEM extent)'
          % (cfg['name'], nr, nc, res, len(field), inside.sum()))
    if inside.sum() < len(field):
        print('  WARNING: %d field heads fall OUTSIDE the DEM extent (dropped)' % (len(field) - inside.sum()))
    field = field[inside]
    R0, R1 = max(0, cr[inside].min() - buf), min(nr, cr[inside].max() + buf)
    C0, C1 = max(0, cc[inside].min() - buf), min(nc, cc[inside].max() + buf)
    zc = z[R0:R1, C0:C1]
    X0, Y0 = west + C0 * res, north - R0 * res
    dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32)

    hull = ConvexHull(field); cen = field.mean(0); hp = field[hull.vertices]
    hb = cen + (hp - cen) * (1.0 + hullbuf / np.linalg.norm(hp - cen, axis=1, keepdims=True))
    tri = Delaunay(hb); in_hull = lambda xy: tri.find_simplex(xy) >= 0
    hull_km2 = ConvexHull(hb).volume / 1e6                 # buffered-hull area (km^2)
    field_in = field[in_hull(field)]
    nf = len(field_in); field_dens = nf / hull_km2
    print('  clip %dx%d; %d field heads in hull (%.3f km2 -> %.1f heads/km2)'
          % (zc.shape[0], zc.shape[1], nf, hull_km2, field_dens))
    print('  T_m2  Tcell  nHd  dens | recall@30 recall@50 prec@50  median_m')
    rows = []
    for T_M2 in T_LIST:
        T = max(2, int(round(T_M2 / (res * res))))
        heads = D.extract_channel_heads(dem, nodata=nd, cellsize=float(res),
                                        valleys='divides', threshold=T)
        hxy = np.array([[X0 + (c + 0.5) * res, Y0 - (r + 0.5) * res] for (r, c) in heads]) \
            if len(heads) else np.empty((0, 2))
        hin = hxy[in_hull(hxy)] if len(hxy) else hxy
        if len(hin) and len(field_in):
            d_ref = cKDTree(hin).query(field_in)[0]
            d_h = cKDTree(field_in).query(hin)[0]
            rec = {t: float(np.mean(d_ref <= t)) for t in TOLS}
            prec = float(np.mean(d_h <= 50.0)); med = float(np.median(d_ref))
        else:
            rec = {t: 0.0 for t in TOLS}; prec = med = float('nan')
        rows.append((T_M2, T, len(hin), len(hin) / hull_km2, rec[30.0], rec[50.0], prec, med))
        print('  %5d %5d %4d %5.1f | %.2f      %.2f      %.2f    %.1f'
              % (T_M2, T, len(hin), len(hin) / hull_km2, rec[30.0], rec[50.0], prec, med), flush=True)
    # density-matched: interpolate each metric to the EXACT field density
    # (defined, repeatable tie-break -- not snap-to-nearest-swept-T)
    A = np.array(rows)
    dens = A[:, 3]; o = np.argsort(dens)                    # ascending density
    interp = lambda col: float(np.interp(field_dens, dens[o], A[o, col]))
    r30, r50, prc, med = interp(4), interp(5), interp(6), interp(7)
    print('  density-matched (interp @ %.1f heads/km2): recall@30 %.2f recall@50 %.2f prec@50 %.2f median %.1f m'
          % (field_dens, r30, r50, prc, med))
    return dict(site=cfg['name'], n_field_in=nf, hull_km2=hull_km2, field_dens=field_dens,
                rows=A, matched=(r30, r50, prc, med))


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    anchor_check_loader()
    print()
    keys = ['mbr', 'indian'] if which == 'all' else [which]
    out = {}
    for k in keys:
        res = run_site(k)
        out[k] = res
        np.savez('/tmp/site_%s.npz' % k, rows=res['rows'], n_field_in=res['n_field_in'],
                 hull_km2=res['hull_km2'], field_dens=res['field_dens'],
                 matched=np.array(res['matched']), site=res['site'])
        print()
    print('saved per-site npz to /tmp/site_<key>.npz')
