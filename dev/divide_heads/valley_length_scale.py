"""Valley head length scale (hillslope length, divide -> channel head) on 1 m
Mid Bailey Run.

Andy's definition: for each divide-defined channel head, the half-width of the
valley head is the distance from the *head-bounding (upstream) divide arc* to the
head. We realise the "head-bounding arc" as the boundary of the head's own
contributing area: every cell that drains through the head is upstream of it by
definition, so that catchment's perimeter ridge IS the upstream arc -- there is
no downstream portion to exclude. The rim-to-head distance is then the hillslope
length L_h (the amphitheater radius), whose companion anchor is drainage density
Dd ~ 1/(2 L_h).

Two distances are stored per rim cell (Andy: more diagnostic info is better):
  * linear   -- straight-line head->rim (sets the grid resolution needed to
                resolve the valley head: ~half-width)
  * flowpath -- along-flow head->rim = fd[rim] - fd[head] (the true hillslope
                flow length; flowpath >> linear flags an elongated/sinuous bowl)

Per head we store the full distribution AND characteristic (median/mean/...) of
both, plus the head's contributing area and an area-equivalent radius sqrt(A/pi)
as a per-head self-consistency anchor.

ANCHOR-CHECK FIRST (per the harness-bug lesson): reproduce the validated 943
heads at T=10000 before trusting any distance number; if the count is off, the
harness is wrong and the lengths are meaningless.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 valley_length_scale.py
"""
import pickle
import time
import numpy as np
from scipy import ndimage
from rivernetworkx import dreich as D
import divides_lib as L

T_VALLEY = 10000
RES = 1.0
WEST, NORTH = L.WEST, L.NORTH

t0 = time.time()
print('loading 1 m cache...', flush=True)
filled = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))['filled']

print('extracting divide-based heads (T=%d) + flowinfo...' % T_VALLEY, flush=True)
heads, fi = D.extract_channel_heads(filled, nodata=-9999.0, cellsize=RES,
                                    fill_dem=False, filled=filled,
                                    threshold=T_VALLEY, valleys='divides',
                                    return_flowinfo=True)
print('  %d heads in %.0fs' % (len(heads), time.time() - t0), flush=True)

# ---- ANCHOR-CHECK: validated pipeline gives 943 heads at this scale ----------
assert 900 <= len(heads) <= 980, \
    'head count %d off the validated 943 -- harness suspect, aborting' % len(heads)
print('  anchor OK (count consistent with validated 943)\n', flush=True)

NodeIndex = fi['NodeIndex']; SVector = fi['SVector']; SVI = fi['SVectorIndex']
ncontrib = fi['ncontrib']; row_of = fi['row_of']; col_of = fi['col_of']; fd = fi['fd']
head_nodes = np.array([NodeIndex[r, c] for (r, c) in heads], dtype=np.int64)

STR8 = np.ones((3, 3), bool)


def rim_distances(h):
    """Return (linear_m, flowpath_m, rim_nodes, area_cells) for head node h.
    Rim = interior boundary cells of h's contributing area (the upstream arc)."""
    s = int(SVI[h]); n = int(ncontrib[h])
    U = SVector[s:s + n]
    rows = row_of[U]; cols = col_of[U]
    r0, c0 = rows.min(), cols.min()
    H = rows.max() - r0 + 3; Wd = cols.max() - c0 + 3
    lr = rows - r0 + 1; lc = cols - c0 + 1
    mask = np.zeros((H, Wd), bool); mask[lr, lc] = True
    locnode = np.full((H, Wd), -1, dtype=np.int64); locnode[lr, lc] = U
    er = ndimage.binary_erosion(mask, structure=STR8, border_value=0)
    rim = mask & ~er
    rim_nodes = locnode[rim]
    rh, ch = row_of[h], col_of[h]
    dr = (row_of[rim_nodes] - rh).astype(np.float64)
    dc = (col_of[rim_nodes] - ch).astype(np.float64)
    linear = RES * np.sqrt(dr * dr + dc * dc)
    flow = (fd[rim_nodes] - fd[h]).astype(np.float64)
    return linear, flow, rim_nodes, n


def stats(a):
    return dict(median=float(np.median(a)), mean=float(np.mean(a)),
                p25=float(np.percentile(a, 25)), p75=float(np.percentile(a, 75)),
                min=float(a.min()), max=float(a.max()), n=int(a.size))


print('computing rim distances per head...', flush=True)
recs = []
lin_all, flow_all = [], []           # per-head object arrays (full distributions)
med_lin, med_flow, areas = [], [], []
tA = time.time()
for h in head_nodes:
    lin, flow, rim_nodes, n = rim_distances(int(h))
    flow = np.clip(flow, 0, None)    # fd float32 round-off can dip <0 by ~um
    recs.append((int(row_of[h]), int(col_of[h]),
                 WEST + (col_of[h] + 0.5) * RES, NORTH - (row_of[h] + 0.5) * RES,
                 n, n * RES * RES, np.sqrt(n * RES * RES / np.pi),
                 stats(lin)['median'], stats(lin)['mean'], lin.min(), lin.max(),
                 stats(flow)['median'], stats(flow)['mean'], flow.max(), rim_nodes.size))
    lin_all.append(lin.astype(np.float32)); flow_all.append(flow.astype(np.float32))
    med_lin.append(np.median(lin)); med_flow.append(np.median(flow)); areas.append(n * RES * RES)
print('  done in %.0fs\n' % (time.time() - tA), flush=True)

med_lin = np.array(med_lin); med_flow = np.array(med_flow); areas = np.array(areas)
r_area = np.sqrt(areas / np.pi)          # area-equivalent radius (full disk)

# ---- independent anchor: drainage density Dd -> 1/(2 Dd) ----------------------
segs = D.channel_network_segments(fi, heads)
chan_cells = sum(len(s['cells']) for s in segs)
chan_len_m = chan_cells * RES            # approx (cardinal-equiv) channel length
valid_area_m2 = fi['N'] * RES * RES
Dd = chan_len_m / valid_area_m2
Lh_from_Dd = 1.0 / (2.0 * Dd)

print('=== valley head length scale (1 m MBR, %d heads) ===' % len(heads))
print('LINEAR  half-width (per-head median, m):  median %.1f  mean %.1f  p25 %.1f  p75 %.1f'
      % (np.median(med_lin), np.mean(med_lin), np.percentile(med_lin, 25), np.percentile(med_lin, 75)))
print('FLOWPATH length    (per-head median, m):  median %.1f  mean %.1f  p25 %.1f  p75 %.1f'
      % (np.median(med_flow), np.mean(med_flow), np.percentile(med_flow, 25), np.percentile(med_flow, 75)))
print('flowpath/linear ratio (per head):         median %.2f' % np.median(med_flow / np.maximum(med_lin, 1e-9)))
print('--- anchors ---')
print('area-equiv radius sqrt(A/pi) (m):         median %.1f   (self-consistency vs linear half-width)'
      % np.median(r_area))
print('drainage density Dd:                      %.4f m/m^2  -> 1/(2Dd) = %.1f m  (independent, basin-avg)'
      % (Dd, Lh_from_Dd))

# ---- persist per-head data ---------------------------------------------------
cols = ['row', 'col', 'x', 'y', 'area_cells', 'area_m2', 'r_area_m',
        'lin_median', 'lin_mean', 'lin_min', 'lin_max',
        'flow_median', 'flow_mean', 'flow_max', 'n_rim']
arr = np.array(recs, dtype=np.float64)
np.savez('/tmp/valley_length_scale_1m.npz',
         columns=np.array(cols), table=arr,
         lin_dist=np.array(lin_all, dtype=object),
         flow_dist=np.array(flow_all, dtype=object),
         Dd=Dd, Lh_from_Dd=Lh_from_Dd)
print('\nper-head table + full distributions -> /tmp/valley_length_scale_1m.npz '
      '(%d heads, columns=%s)' % (len(recs), cols))
print('total %.0fs' % (time.time() - t0))
