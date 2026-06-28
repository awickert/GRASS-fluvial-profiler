"""Faithful Python port of LSDFlowInfo::create() receivers + NContributingNodes +
get_sources_index_threshold, on LSDTT's exact filled DEM. EXACT validation vs the
instrumented C++ rasters _CP (contributing pixels) and _SRC (sources)."""
import numpy as np
from collections import deque

NR, NC = 7107, 6266
ND = -9999
ROOT2INV = np.float32(0.707106781)
THRESH = 100

z = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_fill.flt', dtype='<f4').reshape(NR, NC)
valid = z != np.float32(ND)

# 8 neighbours in LSDTT order: N NE E SE S SW W NW ; even=cardinal, odd=diagonal
OFF = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]
NEG = np.float32(-1e30)
slopes = np.full((8, NR, NC), NEG, dtype=np.float32)
for d, (dr, dc) in enumerate(OFF):
    fac = np.float32(1.0) if d % 2 == 0 else ROOT2INV
    # neighbour elevation, shifted; out-of-bounds stays NoData-equivalent
    nb = np.full((NR, NC), np.float32(ND), dtype=np.float32)
    r0, r1 = max(0, -dr), NR - max(0, dr)
    c0, c1 = max(0, -dc), NC - max(0, dc)
    nb[r0:r1, c0:c1] = z[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
    drop = (z - nb).astype(np.float32) * fac
    ok = valid & (nb != np.float32(ND))
    slopes[d] = np.where(ok, drop, NEG)

best_dir = np.argmax(slopes, axis=0)              # first (lowest-index) max
best_slope = np.max(slopes, axis=0)
# receiver cell: neighbour at best_dir if best_slope > 0 (strict), else self
recv_r = np.arange(NR)[:, None] * np.ones((1, NC), int)
recv_c = np.ones((NR, 1), int) * np.arange(NC)[None, :]
has = valid & (best_slope > np.float32(0.0))
drv = np.array([o[0] for o in OFF]); dcv = np.array([o[1] for o in OFF])
recv_r = np.where(has, recv_r + drv[best_dir], recv_r)
recv_c = np.where(has, recv_c + dcv[best_dir], recv_c)
print('receivers computed')

# contributing nodes (drainage area, pixels) via topological (Kahn) accumulation
rflat = (recv_r * NC + recv_c).ravel()
vflat = valid.ravel()
idx = np.arange(NR * NC)
is_self = (rflat == idx) | (~vflat)
indeg = np.zeros(NR * NC, np.int32)
don = idx[vflat & ~is_self]
np.add.at(indeg, rflat[don], 1)
area = np.where(vflat, 1, 0).astype(np.int64)
q = deque(idx[vflat & (indeg == 0)].tolist())
while q:
    u = q.popleft()
    if is_self[u]:
        continue
    v = rflat[u]
    area[v] += area[u]
    indeg[v] -= 1
    if indeg[v] == 0:
        q.append(v)
CP = area.reshape(NR, NC)

# validate vs C++ _CP
cp_ref = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_CP.flt', dtype='<f4').reshape(NR, NC)
cmp = valid & (cp_ref != ND)
mism = np.sum(CP[cmp] != cp_ref[cmp].astype(np.int64))
print('_CP exact match: %d / %d cells differ (%.4f%%)' % (mism, cmp.sum(), 100 * mism / cmp.sum()))
if mism:
    dr_, dc_ = np.where(cmp & (CP != cp_ref.astype(np.int64)))
    print('  example diffs (mine vs ref):', [(int(CP[r, c]), int(cp_ref[r, c])) for r, c in zip(dr_[:5], dc_[:5])])

# sources: area>=thresh AND all donors < thresh (donor = cell whose receiver is this cell)
chan = valid & (CP >= THRESH)
# a channel cell is a source unless one of its donors is also a channel
donor_is_chan = np.zeros((NR, NC), bool)
cf = chan.ravel()
src_cells = idx[cf]
np.add.at(donor_is_chan.reshape(-1), rflat[idx[vflat & ~is_self & cf]], True)  # mark receivers that have a channel donor
sources = chan & (~donor_is_chan)
nsrc = int(sources.sum())
print('my sources: %d  (C++ N_SOURCES 112757)' % nsrc)

src_ref = np.fromfile('/tmp/dreich_algorithm/bailey_run_dem_SRC.flt', dtype='<f4').reshape(NR, NC)
src_ref_mask = (src_ref != ND)
print('C++ _SRC cells: %d ; exact-cell overlap with mine: %d' %
      (src_ref_mask.sum(), int((sources & src_ref_mask).sum())))

# --- characterize the _CP mismatches ---
dr_, dc_ = np.where(cmp & (CP != cp_ref.astype(np.int64)))
diff = CP[dr_, dc_] - cp_ref[dr_, dc_].astype(np.int64)
print('MISMATCH chars: n=%d ; row range %d-%d ; col range %d-%d' % (len(dr_), dr_.min(), dr_.max(), dc_.min(), dc_.max()))
print('  diff values (mine-ref): min %d max %d ; |diff| median %d' % (diff.min(), diff.max(), int(np.median(np.abs(diff)))))
print('  min CP value among mismatches: %d (high=trunk) ; any mismatch with CP<100 (headwater)? %d' % (CP[dr_,dc_].min(), int(np.sum(CP[dr_,dc_]<100))))
print('  mismatches adjacent to NoData edge: %d of %d' % (int(np.sum(~valid[np.clip(dr_-1,0,NR-1),dc_] | ~valid[np.clip(dr_+1,0,NR-1),dc_] | ~valid[dr_,np.clip(dc_-1,0,NC-1)] | ~valid[dr_,np.clip(dc_+1,0,NC-1)])), len(dr_)))
