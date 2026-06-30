"""Is the order-tail divergence metric fragility, or arbitrary network difference?
Compare Topo/Strahler/Shreve order between my outlet-trim network and TopoToolbox
on the SHARED (geometrically agreed) edges. If the robust schemes (Strahler,
Shreve) agree well where Topo does not, the divergence is the metric's
sensitivity, not non-robust networks.

Needs /tmp/ttb_o_*.txt from ttb_orders.m. Run with PYTHONPATH=/path/to/r.fluvial.
"""
import numpy as np
import rasterio
from rivernetworkx import dreich as D
from rivernetworkx import divides as DV

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
row, col = fi['row_of'], fi['col_of']
_, lab = D.drainage_divides(fi, threshold=THRESH)
sources = D.get_sources(fi, THRESH)
segs = D.channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
outlet = np.zeros((nr, nc), bool)
for s in segs:
    r, c = s['cells'][-1]; outlet[r, c] = True

edges, meta = DV.extract_divide_edges(lab, stream=outlet)
fx = DV.diagonal_flow_corners(fi['recv'], row, col, nr, nc)
edges = DV.split_diagonal_crossings(edges, fx, meta['cw'])
g = DV.build_divide_graph(edges)
cw = meta['cw']


def my_edgeorder(scheme):
    order, _ = DV.order_divides(g, scheme=scheme)
    E = {}
    for k, p in enumerate(g['segments']):
        ii, jj = DV.decode_corner(np.asarray(p) // 2, cw)
        for a in range(len(p) - 1):
            key = frozenset([(int(ii[a]), int(jj[a])), (int(ii[a + 1]), int(jj[a + 1]))])
            E[key] = min(E.get(key, np.inf), order[k])
    return E


ndr = nr + 1
IX = np.genfromtxt('/tmp/ttb_o_IX.txt')


def ttb_edgeorder(scheme):
    od = np.genfromtxt('/tmp/ttb_o_%s.txt' % scheme)
    E = {}
    for i in range(len(IX) - 1):
        a, b = IX[i], IX[i + 1]
        if np.isnan(a) or np.isnan(b):
            continue
        ka = (int((a-1) % ndr), int((a-1)//ndr)); kb = (int((b-1) % ndr), int((b-1)//ndr))
        key = frozenset([ka, kb]); o = np.nanmin([od[i], od[i+1]])
        E[key] = min(E.get(key, np.inf), o)
    return E


print('scheme     my_max ttb_max  shared  Pearson_r  exact%%  within1%%  median|diff|')
for sch in ('topo', 'strahler', 'shreve'):
    me = my_edgeorder(sch); tt = ttb_edgeorder(sch)
    shared = set(me) & set(tt)
    mo = np.array([me[k] for k in shared]); to = np.array([tt[k] for k in shared])
    ok = np.isfinite(mo) & np.isfinite(to)
    mo, to = mo[ok], to[ok]
    r = np.corrcoef(mo, to)[0, 1]
    print('%-9s  %6g %7g  %6d   %.3f     %4.0f    %5.0f      %.0f'
          % (sch, mo.max(), to.max(), len(shared), r,
             100*np.mean(mo == to), 100*np.mean(np.abs(mo-to) <= 1), np.median(np.abs(mo-to))))
