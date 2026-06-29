"""Localize the divide residual (routing is already proven identical). For the
edges TopoToolbox has but I miss (recall gap) and the edges I have but it lacks
(precision gap), report distance-to-stream and Topo order -- telling apart a
too-aggressive stream trim (near-stream gaps) from a basin-partition difference
(ridge gaps).

Run:  PYTHONPATH=/path/to/r.fluvial /usr/bin/python3 residual_analysis.py
"""
import numpy as np
import rasterio
from scipy import ndimage
from rivernetworkx import dreich as D

TIF = '/tmp/topotoolbox/DEMdata/srtm_bigtujunga30m_utm11.tif'
ND = -9999.0; THRESH = 200
R0, R1, C0, C1 = 150, 450, 300, 700

with rasterio.open(TIF) as src:
    dem = src.read(1).astype(np.float64)[R0:R1, C0:C1]
    cs = float(src.transform.a)
    if src.nodata is not None:
        dem[dem == src.nodata] = np.nan
nr, nc = dem.shape
filled = D.fill(dem.astype(np.float32), ND, 0.0001, cs)
fi = D.build_flowinfo(filled, ND, cs); D.contributing_area(fi)
stream = np.zeros((nr, nc), bool)
stream[fi['row_of'], fi['col_of']] = fi['ncontrib'] >= THRESH
# distance (in cells) from each corner to the nearest channel cell
cdist = ndimage.distance_transform_edt(~stream)
cw = nc + 1
cornerdist = np.full((nr + 1, nc + 1), np.inf)
for di in (0, 1):
    for dj in (0, 1):
        sub = cornerdist[di:di + nr, dj:dj + nc]
        np.minimum(sub, cdist, out=sub)


def edge_dist(key):
    return min(cornerdist[i, j] for (i, j) in key)


# rebuild both edge sets (same as compare.py)
ndr = nr + 1
IX = np.genfromtxt('/tmp/ttb_out_IX.txt'); od = np.genfromtxt('/tmp/ttb_out_order.txt')
ttb = {}
for i in range(len(IX) - 1):
    a, b = IX[i], IX[i + 1]
    if np.isnan(a) or np.isnan(b):
        continue
    ka = (int((a - 1) % ndr), int((a - 1) // ndr)); kb = (int((b - 1) % ndr), int((b - 1) // ndr))
    key = frozenset([ka, kb]); o = np.nanmin([od[i], od[i + 1]])
    if key not in ttb or o < ttb[key]:
        ttb[key] = o
mine = np.load('/tmp/ttb_mine_edges.npy')
my = {}
for i1, j1, i2, j2, o in mine:
    key = frozenset([(int(i1), int(j1)), (int(i2), int(j2))])
    if key not in my or o < my[key]:
        my[key] = o

T, Mn = set(ttb), set(my)
miss = T - Mn        # TTB has, I lack (recall gap)
extra = Mn - T       # I have, TTB lacks (precision gap)


def summary(name, edges, omap):
    ds = np.array([edge_dist(k) for k in edges])
    os = np.array([omap[k] for k in edges])
    print('%s: %d edges' % (name, len(edges)))
    print('   dist-to-stream: <=1: %.0f%%  <=2: %.0f%%  >=4: %.0f%%  median=%.1f'
          % (100*np.mean(ds <= 1), 100*np.mean(ds <= 2), 100*np.mean(ds >= 4), np.median(ds)))
    print('   order: ==1: %.0f%%  median=%.0f  max=%.0f'
          % (100*np.mean(os == 1), np.median(os), np.nanmax(os)))


summary('MISS  (TTB has, I lack)', miss, ttb)
summary('EXTRA (I have, TTB lacks)', extra, my)
shared = T & Mn
ds = np.array([edge_dist(k) for k in shared])
print('SHARED: %d edges; dist-to-stream <=1: %.0f%% median=%.1f'
      % (len(shared), 100*np.mean(ds <= 1), np.median(ds)))
