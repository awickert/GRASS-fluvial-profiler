"""Faithful getdivide: open each per-link basin loop at its OUTLET only (the
link's downstream-most cell), keeping the rest of the near-stream perimeter --
vs the earlier all-channel-adjacent trim. Compare both to the TopoToolbox
reference (/tmp/ttb_out_IX.txt from run_divideobj.m) on the same window/routing.

Run:  PYTHONPATH=/path/to/r.fluvial /usr/bin/python3 outlet_trim_experiment.py
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

channel = np.zeros((nr, nc), bool)
channel[row, col] = fi['ncontrib'] >= THRESH

# link outlets = downstream-most cell of each channel link (confluences + outlets)
sources = D.get_sources(fi, THRESH)
segs = D.channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
outlet = np.zeros((nr, nc), bool)
for s in segs:
    r, c = s['cells'][-1]
    outlet[r, c] = True
print('links=%d, outlet cells=%d, channel cells=%d' % (len(segs), outlet.sum(), channel.sum()))

# TopoToolbox reference edges
ndr = nr + 1
IX = np.genfromtxt('/tmp/ttb_out_IX.txt')
ttb = set()
for i in range(len(IX) - 1):
    a, b = IX[i], IX[i + 1]
    if np.isnan(a) or np.isnan(b):
        continue
    ka = (int((a - 1) % ndr), int((a - 1) // ndr)); kb = (int((b - 1) % ndr), int((b - 1) // ndr))
    ttb.add(frozenset([ka, kb]))


def my_edges(stream_mask):
    edges, meta = DV.extract_divide_edges(lab, stream=stream_mask)
    fx = DV.diagonal_flow_corners(fi['recv'], row, col, nr, nc)
    edges = DV.split_diagonal_crossings(edges, fx, meta['cw'])
    g = DV.build_divide_graph(edges)
    order, dd = DV.order_divides(g, scheme='topo')
    E = set()
    for p in g['segments']:
        ii, jj = DV.decode_corner(np.asarray(p) // 2, meta['cw'])
        for a in range(len(p) - 1):
            E.add(frozenset([(int(ii[a]), int(jj[a])), (int(ii[a + 1]), int(jj[a + 1]))]))
    return E, len(g['segments']), int((order == 0).sum()), int(order.max())


M = 12
def interior(k):
    return all(M <= i <= nr - M and M <= j <= nc - M for (i, j) in k)
ttb_i = {k for k in ttb if interior(k)}

for name, mask in [('all-channel trim', channel), ('OUTLET-only trim', outlet)]:
    E, nseg, nun, omax = my_edges(mask)
    inter = E & ttb
    print('\n%-18s edges=%d segs=%d unordered=%d maxorder=%d' % (name, len(E), nseg, nun, omax))
    print('   ALL      recall=%.3f  precision=%.3f  Jaccard=%.3f'
          % (len(inter) / len(ttb), len(inter) / len(E), len(inter) / len(E | ttb)))
    Ei = {k for k in E if interior(k)}; ii = Ei & ttb_i
    print('   INTERIOR recall=%.3f  precision=%.3f  Jaccard=%.3f'
          % (len(ii) / len(ttb_i), len(ii) / len(Ei), len(ii) / len(Ei | ttb_i)))
