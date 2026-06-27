"""Faithful DrEICH reproduction (LSDTT methodology) on the full Bailey window:
find_valleys (curvature) -> per valley: max-flow-distance hilltop + downstream
channel reach -> chi -> chi-z split -> furthest-upstream dedup.
Routing from GRASS; tangential curvature from the C++ reference raster.
Validate against the 634 C++ reference heads."""
import numpy as np
import rivernetworkx as rnx
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command

run_command('g.region', n=4369656.1523925, s=4362549.1523925,
            w=398767.32685636, e=405033.32685636, nsres=1, ewres=1, quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='accd',
            drainage='dird', overwrite=True, quiet=True)
run_command('r.stream.extract', elevation='dem', accumulation='accd', threshold=100,
            stream_raster='strd', direction='sdird', d8cut=0, overwrite=True, quiet=True)
run_command('r.stream.distance', flags='o', stream_rast='strd', direction='dird',
            elevation='dem', method='downstream', distance='fdout', overwrite=True, quiet=True)
acc, b = _read_raster('accd'); dem, _ = _read_raster('dem'); dird, _ = _read_raster('dird')
strd, _ = _read_raster('strd'); fd, _ = _read_raster('fdout')
nr, nc = dem.shape; res = b['ewres']; north = b['north']; west = b['west']; nsr = b['nsres']
A = np.abs(acc)
A0 = 1000.0; MN = 0.525; MINSEG = 10; CURV_THRESH = 0.1; NCONNECT = 10

tcurv = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(7107, 6266)
assert tcurv.shape == dem.shape
tcurv = np.where(tcurv == -9999, np.nan, tcurv)
fd = np.where(np.isfinite(fd), fd, -1.0)

move = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}
def dn(r, c):
    d = dird[r, c]
    if not np.isfinite(d) or int(d) <= 0 or int(d) not in move:
        return None
    dr, dc = move[int(d)]; r2, c2 = r + dr, c + dc
    return (r2, c2) if (0 <= r2 < nr and 0 <= c2 < nc) else None
def up_neighbors(r, c):
    out = []
    for d, (dr, dc) in move.items():
        u, v = r - dr, c - dc
        if 0 <= u < nr and 0 <= v < nc and np.isfinite(dird[u, v]) and int(dird[u, v]) == d:
            out.append((u, v))
    return out

chan = np.isfinite(strd) & (strd > 0)
rr, cc = np.where(chan)
src = [(int(r), int(c)) for r, c in zip(rr, cc) if not any(chan[u] for u in up_neighbors(r, c))]
print('coarse sources: %d' % len(src))

# find_valleys
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
    """find_farthest_upslope_node: climb by max flow-distance-to-outlet."""
    path = [(r, c)]
    for _ in range(maxlen):
        ups = up_neighbors(r, c)
        if not ups:
            break
        nb = max(ups, key=lambda p: fd[p])
        if fd[nb] <= fd[(r, c)]:
            break
        r, c = nb; path.append((r, c))
    return path[::-1]                                   # hilltop ... valley_node
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
    up = farthest_upslope(r, c); dnp = downstream_reach(r, c, A[r, c])
    pr = up + dnp
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
    if idx < 0:
        continue
    head_cells.append((int(rs[idx]), int(cs[idx])))

# furthest-upstream dedup: drop a head if another head drains through it
head_set = set(head_cells); remove = set()
for (r, c) in head_cells:
    rr2, cc2 = r, c
    for _ in range(2000):
        nb = dn(rr2, cc2)
        if nb is None:
            break
        rr2, cc2 = nb
        if (rr2, cc2) in head_set and (rr2, cc2) != (r, c):
            remove.add((rr2, cc2))                       # downstream head has one above it
heads = [h for h in dict.fromkeys(head_cells) if h not in remove]
E = np.array([west + (c + 0.5) * res for (r, c) in heads])
N = np.array([north - (r + 0.5) * nsr for (r, c) in heads])
print('faithful port heads: %d   (C++ reference: 634)' % len(heads))

ref = np.load('/tmp/dreich_ref_heads.npy')
d = np.array([np.min(np.hypot(ref[:, 0] - e, ref[:, 1] - n)) for e, n in zip(E, N)])
print('port-vs-reference (Euclidean): median %.1f m | <=5m %.0f%% | <=10m %.0f%% | <=20m %.0f%%'
      % (np.median(d), 100 * np.mean(d <= 5), 100 * np.mean(d <= 10), 100 * np.mean(d <= 20)))
dr = np.array([np.min(np.hypot(E - e, N - n)) for e, n in zip(ref[:, 0], ref[:, 1])])
print('reference-recall (each ref head -> nearest port): median %.1f m | <=5m %.0f%% | <=10m %.0f%%'
      % (np.median(dr), 100 * np.mean(dr <= 5), 100 * np.mean(dr <= 10)))
np.save('/tmp/dreich_faithful_heads.npy', np.column_stack([E, N]))
run_command('g.remove', type='raster', name='accd,dird,strd,sdird,fdout', flags='f', quiet=True)
