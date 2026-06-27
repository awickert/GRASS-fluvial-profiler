"""DrEICH port prototype v1: chi-z head-finder on real Bailey channels.
Validate against the C++ reference (634 heads) and Clubb's 53 field heads."""
import numpy as np
import openpyxl
import rivernetworkx as rnx
from rivernetworkx.grass_io import _read_raster
from grass.script import run_command, region as _region

# Bailey window matching the reference export
run_command('g.region', n=4369656.1523925, s=4362549.1523925,
            w=398767.32685636, e=405033.32685636, nsres=1, ewres=1, quiet=True)
run_command('r.watershed', flags='s', elevation='dem', accumulation='accd',
            drainage='dird', overwrite=True, quiet=True)
run_command('r.stream.extract', elevation='dem', accumulation='accd', threshold=100,
            stream_raster='strd', direction='sdird', d8cut=0, overwrite=True, quiet=True)
acc, b = _read_raster('accd'); dem, _ = _read_raster('dem')
dird, _ = _read_raster('dird'); strd, _ = _read_raster('strd')
nr, nc = dem.shape; res = b['ewres']; north = b['north']; west = b['west']
A = np.abs(acc)
A0 = 1000.0; MN = 0.525; MINSEG = 10

move = {1: (-1, 1), 2: (-1, 0), 3: (-1, -1), 4: (0, -1), 5: (1, -1), 6: (1, 0), 7: (1, 1), 8: (0, 1)}

def down_neighbor(r, c):
    d = dird[r, c]
    if not np.isfinite(d) or int(d) <= 0 or int(d) not in move:
        return None
    dr, dc = move[int(d)]; r2, c2 = r + dr, c + dc
    if 0 <= r2 < nr and 0 <= c2 < nc:
        return r2, c2
    return None

def up_neighbors(r, c):
    """neighbors that drain INTO (r,c): their down-neighbor is (r,c)."""
    out = []
    for d, (dr, dc) in move.items():
        nr_, nc_ = r - dr, c - dc          # neighbor at -offset flows in via dir d
        if 0 <= nr_ < nr and 0 <= nc_ < nc and np.isfinite(dird[nr_, nc_]) \
                and int(dird[nr_, nc_]) == d:
            out.append((nr_, nc_))
    return out

# coarse channel sources: stream cells with no upstream stream inflow
chan = np.isfinite(strd) & (strd > 0)
src = []
rr, cc = np.where(chan)
for r, c in zip(rr, cc):
    if not any(chan[u] for u in up_neighbors(r, c)):
        src.append((int(r), int(c)))
print('coarse sources (threshold=100): %d' % len(src))

def trace_to_hilltop(r, c, maxlen=4000):
    """Follow max-accumulation upslope contributor to the ridge."""
    path = [(r, c)]
    for _ in range(maxlen):
        ups = up_neighbors(r, c)
        if not ups:
            break
        r, c = max(ups, key=lambda p: A[p])
        path.append((r, c))
    return path[::-1]                      # ridge -> source order? path built source->ridge; reverse

def trace_down(r, c, src_A, maxlen=600):
    path = []
    for _ in range(maxlen):
        nb = down_neighbor(r, c)
        if nb is None:
            break
        r, c = nb; path.append((r, c))
        if A[r, c] > 25 * src_A:           # well past head into an established channel
            break
    return path

heads = []
for (r, c) in src:
    src_A = A[r, c]
    up = trace_to_hilltop(r, c)            # hilltop ... source
    dn = trace_down(r, c, src_A)           # below source
    pr = up + dn                           # hilltop -> downstream
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
    ch = rnx.chi(av, dist, theta=MN, ref_area=A0)          # oriented by area
    # order hilltop(0)->downstream: ensure index 0 is the largest-chi (upslope) end
    if ch[0] < ch[-1]:
        ch = ch[::-1]; z = z[::-1]; rs = rs[::-1]; cs = cs[::-1]
    idx = rnx.channel_head_chi_split(ch, z, min_segment_length=MINSEG)
    if idx < 0:
        continue
    heads.append((west + (cs[idx] + 0.5) * res, north - (rs[idx] + 0.5) * b['nsres']))
heads = np.array(heads)
# spatial dedup at 5 m
snap = 5.0; seen = set(); uq = []
for h in heads:
    k = (round(h[0] / snap), round(h[1] / snap))
    if k not in seen:
        seen.add(k); uq.append(h)
uq = np.array(uq)
print('port v1 heads (deduped): %d' % len(uq))

ref = np.load('/tmp/dreich_ref_heads.npy')
wb = openpyxl.load_workbook('/tmp/clubb_channel_heads.xlsx', data_only=True, read_only=True)
ws = wb['Sheet1']; hE, hN = [], []
for row in ws.iter_rows(min_row=3, values_only=True):
    if row[0] and 'Bailey' in str(row[0]):
        hN.append(float(row[2])); hE.append(float(row[3]))
hE, hN = np.array(hE), np.array(hN)
dref = np.array([np.min(np.hypot(ref[:, 0] - e, ref[:, 1] - n)) for e, n in zip(uq[:, 0], uq[:, 1])])
print('port vs C++ reference: median nearest %.1f m | <=10m %.0f%% | <=20m %.0f%%'
      % (np.median(dref), 100 * np.mean(dref <= 10), 100 * np.mean(dref <= 20)))
dcl = np.array([np.min(np.hypot(uq[:, 0] - e, uq[:, 1] - n)) for e, n in zip(hE, hN)])
print('Clubb53 -> nearest port head: median %.1f m | <=20m %.0f%% | <=50m %.0f%%'
      % (np.median(dcl), 100 * np.mean(dcl <= 20), 100 * np.mean(dcl <= 50)))
np.save('/tmp/dreich_port_v1_heads.npy', uq)
run_command('g.remove', type='raster', name='accd,dird,strd,sdird', flags='f', quiet=True)
