"""Junction-network foundation: build the segment graph from r.stream.extract's
stream raster + drainage direction (reliable raster topology), run
rivernetworkx.strahler_order on the FULL Bailey network. Verifies the
cell->segment->junction plumbing and benchmarks strahler_order vs the
r.stream.order 3h non-finish."""
import time
import numpy as np
import networkx as nx
import rivernetworkx as rnx
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command

run_command('g.region', n=4369656.1523925, s=4362549.1523925,
            w=398767.32685636, e=405033.32685636, nsres=1, ewres=1, quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='accd',
            drainage='dird', overwrite=True, quiet=True)
run_command('r.stream.extract', elevation='dem', accumulation='accd', threshold=100,
            stream_raster='srast', direction='sdird', d8cut=0, overwrite=True, quiet=True)
srast, b = _read_raster('srast'); dird, _ = _read_raster('dird'); accr, _ = _read_raster('accd')
nr, nc = srast.shape
cat = np.where(np.isfinite(srast) & (srast > 0), srast, 0).astype(np.int64)
A = np.abs(accr)
move = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}

t0 = time.time()
# segment topology from the raster: a stream cell is its segment's downstream end
# when its flow-receiver carries a different cat; tostream = that cat (0 = off-map).
rr, cc = np.where(cat > 0)
dv = np.where(np.isfinite(dird[rr, cc]), dird[rr, cc], 0).astype(int)
best = {}                                     # cat -> (accum, tostream) at its true outlet
for d, (dr, dc) in move.items():
    sel = dv == d
    sr, sc = rr[sel], cc[sel]
    tr, tc = sr + dr, sc + dc
    inb = (tr >= 0) & (tr < nr) & (tc >= 0) & (tc < nc)
    own = cat[sr, sc]
    tgt = np.zeros(len(sr), np.int64)
    tgt[inb] = cat[tr[inb], tc[inb]]
    av = A[sr, sc]
    change = own != tgt                       # a boundary cell of segment `own`
    for k, ts, a in zip(own[change], tgt[change], av[change]):
        k = int(k)
        if k not in best or a > best[k][0]:   # true outlet = highest-accumulation boundary cell
            best[k] = (a, int(ts))
edges = {k: v[1] for k, v in best.items()}
t_topo = time.time() - t0
G = nx.DiGraph()
for k, ts in edges.items():
    G.add_edge(k, ts)
print('segments: %d ; graph nodes: %d ; edges: %d ; topology built in %.1fs'
      % (len(edges), G.number_of_nodes(), G.number_of_edges(), t_topo))

isdag = nx.is_directed_acyclic_graph(G)
print('is DAG: %s ; self-loops: %d' % (isdag, nx.number_of_selfloops(G)))
if not isdag:
    ncyc = 0; exlen = None
    for cyc in nx.simple_cycles(G):
        ncyc += 1
        if exlen is None:
            exlen = len(cyc)
        if ncyc >= 50:
            break
    print('cycles found (capped at 50): %d ; example cycle length: %s' % (ncyc, exlen))

t1 = time.time()
so = rnx.strahler_order(G)
t_so = time.time() - t1
vals = np.array([v for v in so.values()])
print('strahler_order on full network: %d segments ordered in %.2fs ; max order %d'
      % (len(so), t_so, int(vals.max())))
import collections
print('order histogram:', dict(sorted(collections.Counter(vals.tolist()).items())))
run_command('g.remove', type='raster', name='accd,dird,srast,sdird', flags='f', quiet=True)
