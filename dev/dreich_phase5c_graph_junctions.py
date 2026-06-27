"""Phase 5c: junction-network-faithful DrEICH using the rivernetworkx segment
graph (Strahler order + 2nd-order-junction anchoring via graph hops). Validate
vs the 634 LSDTT-on-our-DEM reference."""
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
run_command('r.stream.distance', flags='o', stream_rast='srast', direction='dird',
            elevation='dem', method='downstream', distance='fdout', overwrite=True, quiet=True)
srast, b = _read_raster('srast'); dird, _ = _read_raster('dird'); accr, _ = _read_raster('accd')
dem, _ = _read_raster('dem'); fd, _ = _read_raster('fdout')
nr, nc = srast.shape; res = b['ewres']; north = b['north']; west = b['west']; nsr = b['nsres']
cat = np.where(np.isfinite(srast) & (srast > 0), srast, 0).astype(np.int64)
A = np.abs(accr); fd = np.where(np.isfinite(fd), fd, -1.0)
A0 = 1000.0; MN = 0.525; MINSEG = 10; CURV_THRESH = 0.1; NCONNECT = 10
tcurv = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(7107, 6266)
tcurv = np.where(tcurv == -9999, np.nan, tcurv)
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

# segment topology (true outlet = highest-accumulation boundary cell) + graph + Strahler
rr, cc = np.where(cat > 0)
dv = np.where(np.isfinite(dird[rr, cc]), dird[rr, cc], 0).astype(int)
best = {}
for d, (dr, dc) in move.items():
    sel = dv == d; sr, sc = rr[sel], cc[sel]
    tr, tc = sr + dr, sc + dc; inb = (tr >= 0) & (tr < nr) & (tc >= 0) & (tc < nc)
    own = cat[sr, sc]; tgt = np.zeros(len(sr), np.int64); tgt[inb] = cat[tr[inb], tc[inb]]
    av = A[sr, sc]; ch = own != tgt
    for k, ts, a, oc, orow in zip(own[ch], tgt[ch], av[ch], sc[ch], sr[ch]):
        k = int(k)
        if k not in best or a > best[k][0]:
            best[k] = (a, int(ts), (int(orow), int(oc)))
tostream = {k: v[1] for k, v in best.items()}
outlet_cell = {k: v[2] for k, v in best.items()}
G = nx.DiGraph()
for k, ts in tostream.items():
    G.add_edge(k, ts)
SOseg = rnx.strahler_order(G)
SO = np.zeros(nr * nc, np.int16)
flat = cat.ravel()
order_arr = np.array([SOseg.get(int(k), 0) for k in range(flat.max() + 1)], dtype=np.int16)
SO = order_arr[flat].reshape(nr, nc)
print('segments %d ; max Strahler %d' % (len(tostream), int(max(SOseg.values()))))

chan = cat > 0
rrc, ccc = np.where(chan)
src = [(int(r), int(c)) for r, c in zip(rrc, ccc) if not any(chan[u] for u in up_neighbors(r, c))]
print('coarse sources: %d' % len(src))

def second_order_anchor(kv):                  # 2 segment-hops downstream -> that outlet cell
    k1 = tostream.get(kv, 0)
    if k1 == 0:
        return outlet_cell.get(kv)
    k2 = tostream.get(k1, 0)
    return outlet_cell.get(k2, outlet_cell.get(k1, outlet_cell.get(kv)))

visited = np.zeros(dem.shape, bool); vjuncs = {}
for (r0, c0) in src:
    r, c = r0, c0; consec = 0
    while True:
        if not np.isfinite(tcurv[r, c]):
            break
        visited[r, c] = True
        consec = consec + 1 if tcurv[r, c] > CURV_THRESH else 0
        if consec > NCONNECT:
            # the C++ valley junction is the UPSTREAM junction of the valley's
            # link = the valley's own source, kept distinct per valley (not the
            # shared downstream confluence).
            vjuncs[(r0, c0)] = int(cat[r0, c0]); break
        nb = dn(r, c)
        if nb is None or visited[nb]:
            break
        r, c = nb
print('valley junctions: %d  (C++ reference ~768)' % len(vjuncs))

def farthest_upslope(r, c, maxlen=6000):
    path = [(r, c)]
    for _ in range(maxlen):
        ups = up_neighbors(r, c)
        if not ups:
            break
        nb = max(ups, key=lambda p: fd[p])
        if fd[nb] <= fd[(r, c)]:
            break
        r, c = nb; path.append((r, c))
    return path[::-1]
def trace_to(r, c, target, maxlen=3000):
    path = []
    for _ in range(maxlen):
        nb = dn(r, c)
        if nb is None:
            break
        path.append(nb); r, c = nb
        if (r, c) == target:
            break
    return path

head_cells = []
for (vj, kv) in vjuncs.items():
    anchor = second_order_anchor(kv)
    if anchor is None:
        continue
    pr = farthest_upslope(*vj) + trace_to(vj[0], vj[1], anchor)
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
    for _ in range(3000):
        nb = dn(rr2, cc2)
        if nb is None:
            break
        rr2, cc2 = nb
        if (rr2, cc2) in head_set and (rr2, cc2) != (r, c):
            remove.add((rr2, cc2))
heads = [h for h in dict.fromkeys(head_cells) if h not in remove]
E = np.array([west + (c + 0.5) * res for (r, c) in heads])
Nn = np.array([north - (r + 0.5) * nsr for (r, c) in heads])
print('Phase5c heads: %d  (C++ reference 634)  heads/junction %.2f'
      % (len(heads), len(heads) / max(1, len(vjuncs))))
ref = np.load('/tmp/dreich_ref_heads.npy')
d = np.array([np.min(np.hypot(ref[:, 0] - e, ref[:, 1] - n)) for e, n in zip(E, Nn)])
print('Phase5c-vs-reference: median %.1f m | <=5m %.0f%% | <=10m %.0f%% | <=20m %.0f%%'
      % (np.median(d), 100 * np.mean(d <= 5), 100 * np.mean(d <= 10), 100 * np.mean(d <= 20)))
np.save('/tmp/dreich_phase5c_heads.npy', np.column_stack([E, Nn]))
run_command('g.remove', type='raster', name='accd,dird,srast,sdird,fdout', flags='f', quiet=True)
