"""Test the getdivide effect: trim divide edges adjacent to channels so divides
terminate at streams (do not run along them). Compare the ordered network to the
untrimmed one -- does it open the loops (fewer order-0) and move high order off
the valleys onto the interfluves? Uses the cached lab.

Run:  PYTHONPATH=.:dev /usr/bin/python3 dev/divides_getdivide_test.py
"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from rivernetworkx import divides as DV
import divides_lib as L

STREAM_T = 10000
lab = np.load('/tmp/divides_lab_T%d.npy' % STREAM_T)
C = L.load_cache(); fi = C['fi']; filled = C['filled']
nr, nc = lab.shape; cw = nc + 1
stream = np.zeros((nr, nc), bool)
stream[fi['row_of'], fi['col_of']] = fi['ncontrib'] >= STREAM_T


def run(trim):
    edges, meta = DV.extract_divide_edges(lab, stream=stream if trim else None)
    fx = DV.diagonal_flow_corners(fi['recv'], fi['row_of'], fi['col_of'], nr, nc)
    edges = DV.split_diagonal_crossings(edges, fx, cw)
    g = DV.build_divide_graph(edges)
    order, dd = DV.order_divides(g, scheme='topo')   # default: degree-1 endpoints
    return g, order, dd


for trim in (False, True):
    g, order, dd = run(trim)
    seglen = np.array([len(p) - 1 for p in g['segments']])
    un = (order == 0)
    print('%-10s segments=%5d  endpoints=%5d  junctions=%5d  unordered=%4d (%.0f%%)  maxorder=%d'
          % ('TRIM' if trim else 'no-trim', len(g['segments']), len(g['endpoints']),
             len(g['junctions']), un.sum(), 100 * un.mean(), order.max()))

# render the trimmed, ordered network
g, order, dd = run(True)
fig, ax = plt.subplots(figsize=(13, 15))
ax.imshow(np.where(filled == np.float32(-9999), np.nan, filled), cmap='gray',
          extent=[L.WEST, L.WEST + nc, L.NORTH - nr, L.NORTH], origin='upper')
omax = max(1, int(order.max()))
cmap = plt.get_cmap('viridis')
lines, colors, lws = [], [], []
for i, p in enumerate(g['segments']):
    if order[i] == 0:
        continue
    ii, jj = DV.decode_corner(np.array(p) // 2, cw)
    lines.append(np.column_stack([L.WEST + jj, L.NORTH - ii]))
    colors.append(cmap(np.log1p(order[i]) / np.log1p(omax)))   # log colour
    lws.append(0.25 + 0.45 * np.log1p(order[i]))               # log width, modest
ax.add_collection(LineCollection(lines, colors=colors, linewidths=lws))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(1, omax))
fig.colorbar(sm, ax=ax, label='Topo divide order', shrink=0.6)
ax.set_title('MBR divides, getdivide-trimmed + crossing-resolved, Topo order')
fig.savefig('/tmp/mbr_divides_trimmed.png', dpi=130, bbox_inches='tight')
print('saved /tmp/mbr_divides_trimmed.png')
