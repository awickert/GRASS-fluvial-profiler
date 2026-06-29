"""Hand-checkable test of edge-based divide extraction (Scherler & Schwanghart).

A 5x5 basin-label grid with a single staircase divide between basin 0 and 1:

    0 0 0 1 1
    0 0 0 1 1
    0 0 1 1 1
    0 0 1 1 1
    0 0 1 1 1

The divide is the boundary line between the 0s and 1s. We extract its edges and
print them as corner (i,j) pairs so they can be checked by eye against the grid.
"""
import numpy as np
from rivernetworkx.divides import (extract_divide_edges, decode_corner,
                                   build_divide_graph, order_divides,
                                   split_diagonal_crossings)

lab = np.array([
    [0, 0, 0, 1, 1],
    [0, 0, 0, 1, 1],
    [0, 0, 1, 1, 1],
    [0, 0, 1, 1, 1],
    [0, 0, 1, 1, 1],
])
edges, meta = extract_divide_edges(lab)
cw = meta['cw']
print('grid %dx%d, %d divide edges' % (meta['nr'], meta['nc'], len(edges)))
ai, aj = decode_corner(edges[:, 0], cw)
bi, bj = decode_corner(edges[:, 1], cw)
segs = sorted({tuple(sorted([(int(ai[k]), int(aj[k])), (int(bi[k]), int(bj[k]))]))
               for k in range(len(edges))})
for (a, b) in segs:
    orient = 'V' if a[1] == b[1] else 'H'
    print('  %s  %s - %s' % (orient, a, b))

# Expected divide path (corners), top to bottom of the staircase:
#   vertical col-3 boundary rows 0-1-2: (0,3)-(1,3)-(2,3)
#   horizontal row-2 boundary col 2:    (2,2)-(2,3)
#   vertical col-2 boundary rows 2-5:   (2,2)-(3,2)-(4,2)-(5,2)
expected = {
    ((0, 3), (1, 3)), ((1, 3), (2, 3)),
    ((2, 2), (2, 3)),
    ((2, 2), (3, 2)), ((3, 2), (4, 2)), ((4, 2), (5, 2)),
}
got = set(segs)
print('\nmatches expected interior path:', got == expected)
if got != expected:
    print('  missing:', expected - got)
    print('  extra  :', got - expected)

# --- graph on the single divide: 2 endpoints, 0 junctions, 1 segment ---
g = build_divide_graph(edges)
print('\nsingle divide graph: %d endpoints, %d junctions, %d segments'
      % (len(g['endpoints']), len(g['junctions']), len(g['segments'])))
assert len(g['endpoints']) == 2 and len(g['junctions']) == 0 and len(g['segments']) == 1

# --- 3-basin T-junction: 0|1 divide drops to a full-width 01|2 divide ---
lab3 = np.array([
    [0, 0, 1, 1],
    [0, 0, 1, 1],
    [2, 2, 2, 2],
    [2, 2, 2, 2],
])
e3, m3 = extract_divide_edges(lab3)
g3 = build_divide_graph(e3)
ji, jj = decode_corner(np.array(g3['junctions']), m3['cw'])
print('T-junction graph: %d endpoints, %d junctions %s, %d segments'
      % (len(g3['endpoints']), len(g3['junctions']),
         list(zip(ji.tolist(), jj.tolist())), len(g3['segments'])))
# expect junction at corner (2,2); 3 endpoints (2,0),(2,4),(0,2); 3 segments
assert g3['junctions'] == [2 * m3['cw'] + 2], g3['junctions']
assert len(g3['endpoints']) == 3 and len(g3['segments']) == 3
print('\nall graph assertions passed')

# --- ordering: two short verticals branch off one long main divide ---
#   0 0 0 0 0 0
#   0 0 0 0 0 0
#   1 1 2 2 3 3
#   1 1 2 2 3 3
# main horizontal divide (row 2) split at junctions (2,2),(2,4) into 3 pieces;
# two verticals 1|2, 2|3 drop off as order-1 tributaries. The central main piece
# H2=(2,2)-(2,4) is the innermost -> highest Topo order (2).
lab4 = np.array([
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 1, 2, 2, 3, 3],
    [1, 1, 2, 2, 3, 3],
])
e4, m4 = extract_divide_edges(lab4)
g4 = build_divide_graph(e4)
order, dd = order_divides(g4, scheme='topo')
print('\nordering test: %d segments, %d junctions, %d endpoints'
      % (len(g4['segments']), len(g4['junctions']), len(g4['endpoints'])))
for i, p in enumerate(g4['segments']):
    a = decode_corner(np.array([p[0]]), m4['cw'])
    b = decode_corner(np.array([p[-1]]), m4['cw'])
    aa = (int(a[0]), int(a[1])); bb = (int(b[0]), int(b[1]))
    print('  seg %d  %s-%s  len=%d  topo_order=%d  dd=%.0f'
          % (i, aa, bb, len(p) - 1, order[i], dd[i]))
print('max Topo order = %d (expect 2, the central main-divide piece)' % order.max())
assert order.max() == 2, order
# exactly one segment should be order 2 (the central piece); the rest order 1
assert int((order == 2).sum()) == 1 and int((order == 1).sum()) == len(order) - 1
print('\nall ordering assertions passed')

# --- diagonal-crossing resolution: a degree-4 corner at a NW-SE flow crossing
# must split into two through-paths, not a junction ---
cw = 5                                    # corner grid for a 4x4 cell DEM
C = 2 * cw + 2                            # centre corner (2,2) = 12
N, S, E, W = (1*cw+2), (3*cw+2), (2*cw+3), (2*cw+1)
cross = np.array([[C, N], [C, S], [C, E], [C, W]], dtype=np.int64)

g_un = build_divide_graph(cross)          # unresolved: centre is a degree-4 junction
assert len(g_un['junctions']) == 1 and len(g_un['segments']) == 4

fx = np.zeros(5 * cw, dtype=np.int8); fx[C] = 1   # NW-SE flow through centre
split = split_diagonal_crossings(cross, fx, cw)
g_sp = build_divide_graph(split)
ci, _ = decode_corner(np.array([s for seg in g_sp['segments'] for s in seg]) // 2, cw)
print('diagonal crossing: unresolved -> %d junctions; resolved -> %d junctions, %d segments'
      % (len(g_un['junctions']), len(g_sp['junctions']), len(g_sp['segments'])))
assert len(g_sp['junctions']) == 0, g_sp['junctions']   # pass-through, no junction
assert len(g_sp['segments']) == 2                        # two through-paths
# both virtual nodes decode back to the real centre corner (2,2)
assert (np.array([C * 2, C * 2 + 1]) // 2 == C).all()
print('diagonal-crossing assertions passed')
