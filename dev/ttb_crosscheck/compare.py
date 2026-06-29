"""Cross-check, comparison side. TopoToolbox DIVIDEobj (run in Octave on my exact
routing) vs my divide pipeline on the same window. Compares the divide EDGE SET
(undirected corner-pair) and the Topo order on matched edges. Both use the corner
grid (nr+1)x(nc+1); corner (i,j) = NW corner of cell (i,j), 0-based.

Run:  /usr/bin/python3 dev/ttb_compare.py
"""
import numpy as np

nr, nc, cs = np.loadtxt('/tmp/ttb_out_meta.txt', delimiter=',')
nr = int(nr); nc = int(nc)
ndr = nr + 1                                   # divide-grid rows (MATLAB col-major)
IX = np.genfromtxt('/tmp/ttb_out_IX.txt')      # nan-separated 1-based linear idx
od = np.genfromtxt('/tmp/ttb_out_order.txt')   # per-node order


def lin_to_ij(L):                              # 1-based col-major -> 0-based (i,j)
    k = L.astype(np.int64) - 1
    return k % ndr, k // ndr


# TopoToolbox edges: consecutive non-nan nodes within a segment
ttb_edges = {}
i = 0
while i < len(IX) - 1:
    a, b = IX[i], IX[i + 1]
    if not (np.isnan(a) or np.isnan(b)):
        ai, aj = lin_to_ij(np.array([a])); bi, bj = lin_to_ij(np.array([b]))
        key = frozenset([(int(ai[0]), int(aj[0])), (int(bi[0]), int(bj[0]))])
        o = np.nanmin([od[i], od[i + 1]])
        if key not in ttb_edges or o < ttb_edges[key]:
            ttb_edges[key] = o
    i += 1

# my edges
mine = np.load('/tmp/ttb_mine_edges.npy')
my_edges = {}
for i1, j1, i2, j2, o in mine:
    key = frozenset([(int(i1), int(j1)), (int(i2), int(j2))])
    if key not in my_edges or o < my_edges[key]:
        my_edges[key] = o

T = set(ttb_edges); Mn = set(my_edges)
inter = T & Mn
print('TopoToolbox edges: %d   mine: %d   shared: %d' % (len(T), len(Mn), len(inter)))
print('recall (TTB edges I reproduce):   %.3f' % (len(inter) / len(T)))
print('precision (my edges TTB also has): %.3f' % (len(inter) / len(Mn)))
print('Jaccard: %.3f' % (len(inter) / len(T | Mn)))

# order agreement on shared edges
if inter:
    to = np.array([ttb_edges[k] for k in inter])
    mo = np.array([my_edges[k] for k in inter])
    good = np.isfinite(to) & np.isfinite(mo)
    r = np.corrcoef(to[good], mo[good])[0, 1]
    print('\norder on shared edges: TTB max=%g, mine max=%g, Pearson r=%.3f'
          % (np.nanmax(to), mo.max(), r))
    print('order-1 edges: TTB %d, mine %d, shared-as-both-order1 %d'
          % ((to == 1).sum(), (mo == 1).sum(), ((to == 1) & (mo == 1)).sum()))

# interior-only (exclude a boundary margin where the Octave shims + outlets=false
# and window-clipping act) -- the fair test of the core construction
M = 12
def interior(key):
    return all(M <= i <= nr - M and M <= j <= nc - M for (i, j) in key)
Ti = {k for k in T if interior(k)}; Mi = {k for k in Mn if interior(k)}
ii = Ti & Mi
print('\nINTERIOR ONLY (>=%d from edge): TTB %d, mine %d, shared %d' % (M, len(Ti), len(Mi), len(ii)))
print('  recall %.3f   precision %.3f   Jaccard %.3f'
      % (len(ii) / len(Ti), len(ii) / len(Mi), len(ii) / len(Ti | Mi)))
