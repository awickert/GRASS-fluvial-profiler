"""Where are the unordered (order-0) divide segments? Concentrated in the -1
drains-nowhere (endorheic) regions, or spread through a genuine cyclic core?
This decides the faithful fix (clean the -1 inflation vs. change basin construction).
Uses the cached lab; renders order-0 (red) vs ordered (blue) over the DEM.
"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from rivernetworkx import dreich as D
from rivernetworkx import divides as DV
import divides_lib as L

STREAM_T = 10000
lab = np.load('/tmp/divides_lab_T%d.npy' % STREAM_T)
C = L.load_cache(); fi = C['fi']; filled = C['filled']
edges, meta = DV.extract_divide_edges(lab)
g = DV.build_divide_graph(edges)
nr, nc, cw = meta['nr'], meta['nc'], meta['cw']

chan = np.zeros((nr, nc), bool); chan[fi['row_of'], fi['col_of']] = fi['ncontrib'] >= STREAM_T
cs = np.zeros((nr + 1, nc + 1), bool)
for di in (0, 1):
    for dj in (0, 1):
        cs[di:di + nr, dj:dj + nc] |= chan
sf = cs.ravel()
leaf = np.array([bool(sf[np.asarray(p)].any()) for p in g['segments']])
order, dd = DV.order_divides(g, leaf_segments=leaf, scheme='topo')

# nodata-adjacency per segment: does it border a -1 cell?
nd = np.zeros((nr + 1, nc + 1), bool)
ndcell = (lab == -1)
for di in (0, 1):
    for dj in (0, 1):
        nd[di:di + nr, dj:dj + nc] |= ndcell
ndflag = nd.ravel()
seg_touches_nd = np.array([bool(ndflag[np.asarray(p)].any()) for p in g['segments']])
seglen = np.array([len(p) - 1 for p in g['segments']])

un = order == 0
print('segments: %d  unordered(order0): %d (%.0f%%)' % (len(order), un.sum(), 100 * un.mean()))
print('unordered touching a -1 region: %.2f   ordered touching -1: %.2f'
      % (seg_touches_nd[un].mean(), seg_touches_nd[~un].mean()))
print('unordered seg length: med=%.0f p90=%.0f   ordered: med=%.0f'
      % (np.median(seglen[un]), np.percentile(seglen[un], 90), np.median(seglen[~un])))

fig, ax = plt.subplots(figsize=(13, 15))
ax.imshow(np.where(filled == np.float32(-9999), np.nan, filled), cmap='gray',
          extent=[L.WEST, L.WEST + nc, L.NORTH - nr, L.NORTH], origin='upper')
ax.imshow(np.where(ndcell, 1.0, np.nan), cmap='autumn', alpha=0.25,
          extent=[L.WEST, L.WEST + nc, L.NORTH - nr, L.NORTH], origin='upper')
lines_u, lines_o = [], []
for i, p in enumerate(g['segments']):
    ii, jj = DV.decode_corner(np.array(p), cw)
    xy = np.column_stack([L.WEST + jj, L.NORTH - ii])
    (lines_u if order[i] == 0 else lines_o).append(xy)
ax.add_collection(LineCollection(lines_o, colors='blue', linewidths=0.4))
ax.add_collection(LineCollection(lines_u, colors='red', linewidths=0.4))
ax.set_title('MBR divides: ordered (blue) vs unordered/order-0 (red); -1 regions shaded')
out = '/tmp/mbr_divides_unordered.png'
fig.savefig(out, dpi=130, bbox_inches='tight')
print('saved', out)
