"""Isolate where Design A goes wrong: feed find_valleys_gated the CURVATURE signal
and check it reproduces the 634 reference heads. If yes -> walk+compare are fine,
the divide signal is the problem. Then inspect divide-method head locations."""
import numpy as np
from rivernetworkx import dreich as D
import divides_lib as L

C = L.load_cache()
fi, filled, tcurv, ref_heads = C['fi'], C['filled'], C['tcurv'], C['ref_heads']
sources = C['sources']
ref_rc = L.heads_rc(fi, ref_heads)

# sanity: compare reference to itself
print('self-compare ref vs ref:', L.compare_heads(ref_rc, ref_rc))

# 1) curvature signal through find_valleys_gated
gate_c = (tcurv > 0.1)
gnode = L.gate_from_grid(fi, gate_c)
val_c = L.find_valleys_gated(fi, gnode, sources, n_connecting_nodes=10)
nv = len(np.unique(val_c[val_c != L.ND]))
_, heads_c = D.channel_heads_from_valleys(fi, filled, val_c, 10, 1000.0, 0.525)
print('curv-gate via find_valleys_gated: %d valleys, %d heads' % (nv, len(heads_c)))
print('  vs ref:', L.compare_heads(L.heads_rc(fi, heads_c), ref_rc))

# how many heads are node 0 (the no-positive-test degenerate)?
heads_c = np.asarray(heads_c)
print('  heads at node 0:', int((heads_c == 0).sum()))
print('  sample head (r,c):', list(zip(fi['row_of'][heads_c[:5]].tolist(),
                                       fi['col_of'][heads_c[:5]].tolist())))

# 2) one divide-gate run, inspect head locations
divide, _ = D.drainage_divides(fi, threshold=50000)
grid = L.valley_flank_grid(divide, 120)
val_d = L.find_valleys_gated(fi, L.gate_from_grid(fi, grid), sources, 10)
_, heads_d = D.channel_heads_from_valleys(fi, filled, val_d, 10, 1000.0, 0.525)
heads_d = np.asarray(heads_d)
print('\ndivide-gate(50000,W120): %d heads, %d at node0'
      % (len(heads_d), int((heads_d == 0).sum())))
nz = heads_d[heads_d != 0]
print('  head areas: med=%d p10=%d p90=%d'
      % (np.median(fi['ncontrib'][nz]), np.percentile(fi['ncontrib'][nz], 10),
         np.percentile(fi['ncontrib'][nz], 90)))
print('  ref  head areas: med=%d' % np.median(fi['ncontrib'][ref_heads]))
print('  vs ref:', L.compare_heads(L.heads_rc(fi, heads_d), ref_rc))
