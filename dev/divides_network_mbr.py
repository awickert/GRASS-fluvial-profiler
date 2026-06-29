"""Run the Scherler & Schwanghart divide network on 1 m Mid Bailey Run: extract
the edge divides, build the graph, order it (Topo), and render the ordered network
+ report the order / length distributions. This is the artifact to look at
together: the deterministic divide hierarchy, where short low-order divides are
candidate noise (coarsening / messiness) and long high-order divides are signal.

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_network_mbr.py [stream_threshold]
"""
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from rivernetworkx import dreich as D
from rivernetworkx import divides as DV
import divides_lib as L

STREAM_T = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
WEST, NORTH = L.WEST, L.NORTH
RES = 1.0

T0 = time.time()
def log(m): print('[%6.1fs] %s' % (time.time() - T0, m), flush=True)

C = L.load_cache()
fi, filled = C['fi'], C['filled']
log('drainage_divides(threshold=%d) -> per-cell basin labels' % STREAM_T)
_, lab = D.drainage_divides(fi, threshold=STREAM_T)

log('extract_divide_edges')
edges, meta = DV.extract_divide_edges(lab)
log('  %d divide edges' % len(edges))

# Stream-coincident corners = divide endpoints (S&S): a corner touching a channel
# cell (contributing area >= STREAM_T). Divides terminate where they reach a river.
log('stream-coincident terminal corners')
chan = np.zeros((meta['nr'], meta['nc']), dtype=bool)
chan[fi['row_of'], fi['col_of']] = fi['ncontrib'] >= STREAM_T
cw = meta['cw']
corner_stream = np.zeros((meta['nr'] + 1, meta['nc'] + 1), dtype=bool)
for di in (0, 1):
    for dj in (0, 1):                          # corner (i,j) touches cells (i-di, j-dj)
        corner_stream[di:di + meta['nr'], dj:dj + meta['nc']] |= chan
streamflag = corner_stream.ravel()            # indexed by REAL corner id (i*cw + j)
log('  %d stream-coincident corners' % int(streamflag.sum()))

# resolve D8 diagonal-flow crossings (S&S FX): split such corners into through
# nodes so the divide passes through rather than forming a spurious junction
log('diagonal_flow_corners + split_diagonal_crossings')
fx = DV.diagonal_flow_corners(fi['recv'], fi['row_of'], fi['col_of'],
                              meta['nr'], meta['nc'])
log('  %d diagonal-crossing corners flagged' % int((fx > 0).sum()))
edges = DV.split_diagonal_crossings(edges, fx, cw)     # now in doubled id space

log('build_divide_graph')
g = DV.build_divide_graph(edges)
log('  %d endpoints, %d junctions, %d segments'
    % (len(g['endpoints']), len(g['junctions']), len(g['segments'])))

# leaf = a divide segment that reaches a river (any corner stream-coincident).
# nodes are in doubled space now, so the real corner id is node // 2.
leaf = np.array([bool(streamflag[np.asarray(p) // 2].any()) for p in g['segments']])
log('order_divides (Topo); %d/%d segments are stream-terminating leaves'
    % (leaf.sum(), len(leaf)))
order, dd = DV.order_divides(g, leaf_segments=leaf, scheme='topo')
seglen = np.array([len(p) - 1 for p in g['segments']])
log('  order: max=%d  dist over orders: %s'
    % (order.max(), np.bincount(order).tolist()))
print('segment length (corner steps): med=%.0f p90=%.0f max=%.0f'
      % (np.median(seglen), np.percentile(seglen, 90), seglen.max()))
print('order-1 segments: %d (%.0f%%); their median length=%.0f vs order>=2 median=%.0f'
      % ((order == 1).sum(), 100.0 * (order == 1).mean(),
         np.median(seglen[order == 1]),
         np.median(seglen[order >= 2]) if (order >= 2).any() else -1))

# --- figure: ordered divide network over the DEM ---
log('rendering figure')
cw = meta['cw']
fig, ax = plt.subplots(figsize=(14, 16))
ax.imshow(np.where(filled == np.float32(-9999), np.nan, filled),
          cmap='gray', extent=[WEST, WEST + meta['nc'] * RES,
                               NORTH - meta['nr'] * RES, NORTH], origin='upper')
omax = int(order.max())
lines, colors, lws = [], [], []
cmap = plt.get_cmap('viridis')
for i, p in enumerate(g['segments']):
    if order[i] == 0:
        continue
    ii, jj = DV.decode_corner(np.array(p) // 2, cw)
    xy = np.column_stack([WEST + jj * RES, NORTH - ii * RES])
    lines.append(xy)
    colors.append(cmap((order[i] - 1) / max(1, omax - 1)))
    lws.append(0.3 + 0.5 * order[i])
ax.add_collection(LineCollection(lines, colors=colors, linewidths=lws))
ax.set_title('MBR divide network, Topo order (stream T=%d) -- %d segments'
             % (STREAM_T, len(g['segments'])))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(1, omax))
fig.colorbar(sm, ax=ax, label='Topo divide order', shrink=0.6)
out = '/tmp/mbr_divide_network_T%d.png' % STREAM_T
fig.savefig(out, dpi=130, bbox_inches='tight')
log('saved %s' % out)
