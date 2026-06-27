"""Definitive reproduction test: route on LSDTT's OWN filled DEM
(bailey_run_dem_fill.flt) with steepest-descent D8 computed in Python, then run
the faithful DrEICH pipeline. Removes both the fill difference and the routing-
algorithm difference, so a close match to the 634 reference would confirm the
methodology port is faithful and isolate the ~40 m to fill+routing."""
import numpy as np
from collections import deque
import rivernetworkx as rnx
from grass.script import run_command

run_command('g.region', n=4369656.1523925, s=4362549.1523925,
            w=398767.32685636, e=405033.32685636, nsres=1, ewres=1, quiet=True)
nr, nc = 7107, 6266
res = 1.0; north = 4369656.1523925; west = 398767.32685636; nsr = 1.0
A0 = 1000.0; MN = 0.525; MINSEG = 10; CURV_THRESH = 0.1; NCONNECT = 10

dem = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_fill.flt', dtype='<f4').reshape(nr, nc)
dem = np.where(dem == -9999, np.nan, dem).astype(np.float64)
tcurv = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(nr, nc)
tcurv = np.where(tcurv == -9999, np.nan, tcurv)
valid = np.isfinite(dem)
N = nr * nc
move = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}

# steepest-descent D8 receiver on the LSDTT-filled surface
recv = np.full(N, -1, dtype=np.int64)
steplen = np.zeros(N, dtype=np.float64)
best = np.zeros((nr, nc), dtype=np.float64)        # require slope > 0 (downhill)
for d, (dr, dc) in move.items():
    dist = res if (dr == 0 or dc == 0) else res * 1.41421356
    r0, r1 = max(0, -dr), nr - max(0, dr)
    c0, c1 = max(0, -dc), nc - max(0, dc)
    src = dem[r0:r1, c0:c1]; nbr = dem[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
    sl = (src - nbr) / dist
    better = np.isfinite(sl) & (sl > best[r0:r1, c0:c1])
    best[r0:r1, c0:c1][better] = sl[better]
    rr, cc = np.where(better)
    sidx = (rr + r0) * nc + (cc + c0); tidx = (rr + r0 + dr) * nc + (cc + c0 + dc)
    recv[sidx] = tidx; steplen[sidx] = dist
print('steepest-descent D8 built on LSDTT filled DEM')

indeg = np.zeros(N, dtype=np.int32); has = recv >= 0
np.add.at(indeg, recv[has], 1)
acc = np.where(valid.ravel(), 1.0, 0.0)
q = deque(np.where((indeg == 0) & valid.ravel())[0].tolist()); order = []
while q:
    u = q.popleft(); order.append(u); v = recv[u]
    if v >= 0:
        acc[v] += acc[u]; indeg[v] -= 1
        if indeg[v] == 0:
            q.append(v)
fd = np.zeros(N)
for u in reversed(order):
    v = recv[u]
    if v >= 0:
        fd[u] = fd[v] + steplen[u]
A = acc.reshape(nr, nc); FD = fd.reshape(nr, nc)
print('accumulation + flow distance built')

def dn(r, c):
    v = recv[r * nc + c]
    return (v // nc, v % nc) if v >= 0 else None
def up_neighbors(r, c):
    idx = r * nc + c; out = []
    for d, (dr, dc) in move.items():
        u, w = r - dr, c - dc
        if 0 <= u < nr and 0 <= w < nc and recv[u * nc + w] == idx:
            out.append((u, w))
    return out

chan = (A >= 100) & valid
rr, cc = np.where(chan)
src = [(int(r), int(c)) for r, c in zip(rr, cc) if not any(chan[u] for u in up_neighbors(r, c))]
print('coarse sources: %d' % len(src))
visited = np.zeros((nr, nc), bool); valley_nodes = []
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
print('LSDTT-fill port heads: %d  (C++ reference 634)' % len(heads))
ref = np.load('/tmp/dreich_ref_heads.npy')
d = np.array([np.min(np.hypot(ref[:, 0] - e, ref[:, 1] - n)) for e, n in zip(E, Nn)])
print('LSDTT-fill-vs-reference (Euclidean): median %.1f m | <=5m %.0f%% | <=10m %.0f%% | <=20m %.0f%%'
      % (np.median(d), 100 * np.mean(d <= 5), 100 * np.mean(d <= 10), 100 * np.mean(d <= 20)))
np.save('/tmp/dreich_lsdttfill_heads.npy', np.column_stack([E, Nn]))
