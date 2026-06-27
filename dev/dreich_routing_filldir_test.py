"""Compatibility diagnostic: run the faithful DrEICH pipeline on r.fill.dir
STEEPEST-DESCENT D8 routing (vs r.watershed least-cost) to see how much of the
~40 m head displacement vs the C++ reference is routing-algorithm choice.
Accumulation + flow-distance derived in Python from the r.fill.dir direction."""
import numpy as np
from collections import deque
import rivernetworkx as rnx
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command

run_command('g.region', n=4369656.1523925, s=4362549.1523925,
            w=398767.32685636, e=405033.32685636, nsres=1, ewres=1, quiet=True)
run_command('r.fill.dir', input='dem', output='filled', direction='fdir',
            overwrite=True, quiet=True)
dem, b = _read_raster('filled'); deg, _ = _read_raster('fdir')
nr, nc = dem.shape; res = b['ewres']; north = b['north']; west = b['west']; nsr = b['nsres']
A0 = 1000.0; MN = 0.525; MINSEG = 10; CURV_THRESH = 0.1; NCONNECT = 10

tcurv = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(7107, 6266)
tcurv = np.where(tcurv == -9999, np.nan, tcurv)

# direction degrees -> 1..8 (45->1 NE ... 360->8 E), same convention as r.watershed
code = np.where(np.isfinite(deg), np.round(deg / 45.0).astype(int), 0)
code[code == 0] = 8 if False else code[code == 0]    # keep 0 as nodata sentinel
move = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}
valid = np.isfinite(dem)
N = nr * nc
recv = np.full(N, -1, dtype=np.int64)
steplen = np.zeros(N, dtype=np.float32)
for d, (dr, dc) in move.items():
    m = (code == d) & valid
    rs, cs = np.where(m)
    tr, tc = rs + dr, cs + dc
    inb = (tr >= 0) & (tr < nr) & (tc >= 0) & (tc < nc)
    src_idx = (rs * nc + cs)
    tgt_idx = np.where(inb, (tr * nc + tc), -1)
    recv[src_idx] = tgt_idx
    steplen[src_idx] = (res if (dr == 0 or dc == 0) else res * 1.41421356)
# Kahn topological accumulation + flow distance to outlet
indeg = np.zeros(N, dtype=np.int32)
has = recv >= 0
np.add.at(indeg, recv[has], 1)
acc = np.where(valid.ravel(), 1.0, 0.0).astype(np.float64)
q = deque(np.where((indeg == 0) & valid.ravel())[0].tolist())
order = []
while q:
    u = q.popleft(); order.append(u); v = recv[u]
    if v >= 0:
        acc[v] += acc[u]; indeg[v] -= 1
        if indeg[v] == 0:
            q.append(v)
fd = np.zeros(N, dtype=np.float64)
for u in reversed(order):
    v = recv[u]
    if v >= 0:
        fd[u] = fd[v] + steplen[u]
A = acc.reshape(nr, nc)
FD = fd.reshape(nr, nc)
print('routing built (steepest-descent D8 from r.fill.dir)')

def dn(r, c):
    v = recv[r * nc + c]
    return (v // nc, v % nc) if v >= 0 else None
def up_neighbors(r, c):
    out = []
    for d, (dr, dc) in move.items():
        u, w = r - dr, c - dc
        if 0 <= u < nr and 0 <= w < nc and code[u, w] == d and valid[u, w]:
            out.append((u, w))
    return out

# coarse sources: accumulation >= 100 with no upstream cell >= 100
chan = (A >= 100) & valid
rr, cc = np.where(chan)
src = [(int(r), int(c)) for r, c in zip(rr, cc) if not any(chan[u] for u in up_neighbors(r, c))]
print('coarse sources: %d' % len(src))

visited = np.zeros(dem.shape, bool); valley_nodes = []
for (r0, c0) in src:
    r, c = r0, c0; consec = 0
    while True:
        if not np.isfinite(tcurv[r, c]):
            break
        visited[r, c] = True
        consec = consec + 1 if tcurv[r, c] > CURV_THRESH else 0
        if consec > NCONNECT:
            valley_nodes.append((r, c)); break
        nb = dn(r, c)
        if nb is None or visited[nb]:
            break
        r, c = nb
print('valley nodes: %d' % len(valley_nodes))

def farthest_upslope(r, c, maxlen=6000):
    path = [(r, c)]
    for _ in range(maxlen):
        ups = up_neighbors(r, c)
        if not ups:
            break
        nb = max(ups, key=lambda p: FD[p])
        if FD[nb] <= FD[(r, c)]:
            break
        r, c = nb; path.append((r, c))
    return path[::-1]
def downstream_reach(r, c, node_A, maxlen=800):
    path = []
    for _ in range(maxlen):
        nb = dn(r, c)
        if nb is None:
            break
        r, c = nb; path.append((r, c))
        if A[r, c] > 8 * node_A:
            break
    return path

head_cells = []
for (r, c) in valley_nodes:
    pr = farthest_upslope(r, c) + downstream_reach(r, c, A[r, c])
    if len(pr) < 2 * MINSEG + 2:
        continue
    rs = np.array([p[0] for p in pr]); cs = np.array([p[1] for p in pr])
    z = dem[rs, cs]; av = A[rs, cs]
    step = np.hypot(np.diff(rs.astype(float)), np.diff(cs.astype(float))) * res
    dist = np.concatenate([[0], np.cumsum(step)])
    ok = np.isfinite(z) & np.isfinite(av) & (av > 0)
    if ok.sum() < 2 * MINSEG + 2:
        continue
    rs, cs, z, av, dist = rs[ok], cs[ok], z[ok], av[ok], dist[ok]
    ch = rnx.chi(av, dist, theta=MN, ref_area=A0)
    if ch[0] < ch[-1]:
        ch = ch[::-1]; z = z[::-1]; rs = rs[::-1]; cs = cs[::-1]
    idx = rnx.channel_head_chi_split(ch, z, min_segment_length=MINSEG)
    if idx >= 0:
        head_cells.append((int(rs[idx]), int(cs[idx])))
head_set = set(head_cells); remove = set()
for (r, c) in head_cells:
    rr2, cc2 = r, c
    for _ in range(2000):
        nb = dn(rr2, cc2)
        if nb is None:
            break
        rr2, cc2 = nb
        if (rr2, cc2) in head_set and (rr2, cc2) != (r, c):
            remove.add((rr2, cc2))
heads = [h for h in dict.fromkeys(head_cells) if h not in remove]
E = np.array([west + (c + 0.5) * res for (r, c) in heads])
Nn = np.array([north - (r + 0.5) * nsr for (r, c) in heads])
print('fill.dir port heads: %d  (C++ reference 634)' % len(heads))
ref = np.load('/tmp/dreich_ref_heads.npy')
d = np.array([np.min(np.hypot(ref[:, 0] - e, ref[:, 1] - n)) for e, n in zip(E, Nn)])
print('fill.dir-vs-reference (Euclidean): median %.1f m | <=5m %.0f%% | <=10m %.0f%% | <=20m %.0f%%'
      % (np.median(d), 100 * np.mean(d <= 5), 100 * np.mean(d <= 10), 100 * np.mean(d <= 20)))
run_command('g.remove', type='raster', name='filled,fdir', flags='f', quiet=True)
