"""Faithful port of LSDRaster::fill(MinSlope) -- priority-flood depression filling
with an imposed minimum gradient (Wang & Liu / Barnes, +epsilon).
Transcribed from dreich_algorithm/LSDRaster.cpp:3241.

Validate vs the reference filled DEM: fill(raw bailey_run_dem) == bailey_run_dem_fill.
Tie-breaking in flats follows insertion order (the C++ std::priority_queue tie order
is unspecified), so flat interiors may differ by the epsilon pattern -- these do not
affect headwater channel heads (same class as the 1196 trunk _CP float ties).
"""
import heapq
import numpy as np

NR, NC = 7107, 6266
ND = np.float32(-9999.0)
ALG = '/tmp/dreich_algorithm'
ONE_OVER_ROOT2 = np.float32(0.707106781)

# neighbour order N NE E SE S SW W NW (even cardinal, odd diagonal) -- as fill()'s kernal
OFF = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]


def fill(z, min_slope=0.0001, res=1.0):
    """Return the depression-filled DEM (float32). z: float32 array, NoData == ND."""
    z = z.astype(np.float32)
    nr, nc = z.shape
    filled = z.copy()
    data = z != ND
    # FillIndex: 0 = data unprocessed, 1 = in queue, 2 = processed, -1 = nodata
    state = np.where(data, 0, -1).astype(np.int8)
    inc_card = np.float32(min_slope * res)
    inc_diag = np.float32(min_slope * res * float(ONE_OVER_ROOT2))

    # seed the queue with boundary / NoData-adjacent data cells
    pad = np.ones((nr, nc), bool)
    pad[1:-1, 1:-1] = False                              # array edge
    ndmask = ~data
    nbr_nd = np.zeros((nr, nc), bool)
    for dr, dc in OFF:
        sl = np.zeros((nr, nc), bool)
        r0, r1 = max(0, -dr), nr - max(0, dr)
        c0, c1 = max(0, -dc), nc - max(0, dc)
        sl[r0:r1, c0:c1] = ndmask[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        nbr_nd |= sl
    seed = data & (pad | nbr_nd)

    heap = []
    push = heapq.heappush
    cnt = 0
    sr, sc = np.where(seed)
    for r, c in zip(sr.tolist(), sc.tolist()):
        push(heap, (float(filled[r, c]), cnt, r, c)); cnt += 1
        state[r, c] = 1

    while heap:
        zc, _, row, col = heapq.heappop(heap)
        state[row, col] = 2
        zc = np.float32(zc)
        for k, (dr, dc) in enumerate(OFF):
            nr_, nc_ = row + dr, col + dc
            if nr_ < 0 or nr_ >= nr or nc_ < 0 or nc_ >= nc:
                continue
            if state[nr_, nc_] != 0:                     # nodata / in-queue / processed
                continue
            if filled[nr_, nc_] <= zc:
                if min_slope > 0:
                    filled[nr_, nc_] = np.float32(zc + (inc_card if k % 2 == 0 else inc_diag))
                else:
                    filled[nr_, nc_] = zc
            push(heap, (float(filled[nr_, nc_]), cnt, nr_, nc_)); cnt += 1
            state[nr_, nc_] = 1
    return filled


if __name__ == '__main__':
    import time
    raw = np.fromfile(f'{ALG}/bailey_run_dem.flt', dtype='<f4').reshape(NR, NC)
    ref = np.fromfile(f'{ALG}/bailey_run_dem_fill.flt', dtype='<f4').reshape(NR, NC)
    t = time.time()
    mine = fill(raw, min_slope=0.0001, res=1.0)
    print('fill %.1fs' % (time.time() - t), flush=True)
    data = (raw != ND) & (ref != ND)
    d = mine[data].astype(np.float64) - ref[data].astype(np.float64)
    ad = np.abs(d)
    print('vs reference fill: cells=%d  exact=%d (%.4f%%)  max|d|=%.4g median|d|=%.4g'
          % (data.sum(), int((ad == 0).sum()), 100 * (ad == 0).mean(), ad.max(), np.median(ad)))
    # where do differences sit? raised-above-raw cells (depressions) are where fill acts
    raised = data & (ref.astype(np.float64) > raw.astype(np.float64) + 1e-6)
    print('reference raised %d cells above raw; of those mine differs from ref at %d'
          % (int(raised.sum()), int((raised & (ad > 1e-4)).sum())))
    np.save('/tmp/dreich_fill_mine.npy', mine)
