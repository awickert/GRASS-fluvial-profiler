"""
DrEICH channel-head extraction (Clubb et al., 2014), a faithful Python/numpy port
of the LSDTopoTools pipeline used by ``channel_heads_driver`` (the standalone
``dreich_algorithm`` repo).

Pure numpy, GRASS-free. The GRASS module ``r.fluvial.channelheads`` (method
``dreich``) is a thin wrapper that reads a DEM, calls :func:`extract_channel_heads`,
and writes the head points.

Pipeline (all on a depression-filled DEM):

    fill -> D8 FlowInfo (receivers, donor stack, ordered stack, contributing area,
    flow distance) -> sources (area threshold) -> junction network (stream order,
    junctions) -> tangential curvature (2nd-order polyfit) -> find_valleys
    (sustained high curvature) -> per-valley hilltop->junction chi profile ->
    chi-z split (R2 / Durbin-Watson) -> furthest-upstream dedup.

Validation (Bailey 7107x6266, vs an instrumented C++ run): every stage is
bit-exact -- routing/sources (112757), valley junctions (768), hilltops (768/768),
head count (634), and 630/634 head locations. The few residual heads are float32-
degenerate chi-z split optima (the test statistic is flat to <1e-6 on near-singular
2x2 regressions); this code runs the algorithm in float where LSDTT does, but is a
*faithful port*, not a bit-for-bit float emulator. See the project notes.

Faithfulness details preserved from the C++:
  * D8 receivers: neighbour order N NE E SE S SW W NW; slope = drop (cardinal) or
    0.707106781*drop (diagonal); max_slope init 0 strict, first-max tie-break;
    no-flux edges (a cell with no downhill neighbour is base level).
  * tangential curvature: 2nd-order polyfit, window radius 7,
    floor(radial_dist) <= window_radius mask, formula
    (fxx*fy^2 - 2*fxy*fx*fy + fyy*fx^2)/(p*sqrt(q)).
  * chi: dx*(A_0/(NContributingNodes*pixel_area))^(m/n), accumulated upstream.
  * chi-z split: float32 simple_linear_regression (sequential summation + JAMA-LU
    2x2 solve) + Durbin-Watson; test = R2_channel - (DW_hillslope - 2)/2.
"""
import heapq
import numpy as np

ROOT2INV = np.float32(0.707106781)        # LSDTT's literal 1/sqrt(2)
# neighbour order N NE E SE S SW W NW; even index cardinal, odd diagonal
_OFF = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]

# r.watershed drainage-direction encoding: |dir| in 1..8, CCW from NE
#   1 NE  2 N  3 NW  4 W  5 SW  6 S  7 SE  8 E   (negative = drains off-region)
# |dir| -> (drow, dcol), indexed by |dir| (slot 0 unused):
_RWS_DROW = np.array([0, -1, -1, -1, 0, 1, 1, 1, 0])
_RWS_DCOL = np.array([0, 1, 0, -1, -1, -1, 0, 1, 1])
_RWS_CARDINAL = frozenset({2, 4, 6, 8})
# inverse: (drow, dcol) -> r.watershed code (for encoding receivers to directions)
_RWS_CODE = {(-1, 1): 1, (-1, 0): 2, (-1, -1): 3, (0, -1): 4,
             (1, -1): 5, (1, 0): 6, (1, 1): 7, (0, 1): 8}


# ---------------------------------------------------------------- depression fill
def fill(dem, nodata=-9999.0, min_slope=0.0001, cellsize=1.0):
    """Priority-flood depression fill with an imposed minimum gradient
    (LSDRaster::fill; Wang & Liu / Barnes +epsilon). Returns a float32 array."""
    ndf = np.float32(nodata)
    z = dem.astype(np.float32)
    nr, nc = z.shape
    filled = z.copy()
    data = z != ndf
    state = np.where(data, 0, -1).astype(np.int8)         # 0 unproc, 1 queued, 2 done, -1 nodata
    inc_card = np.float32(min_slope * cellsize)
    inc_diag = np.float32(min_slope * cellsize * float(ROOT2INV))

    pad = np.ones((nr, nc), bool); pad[1:-1, 1:-1] = False
    ndmask = ~data
    nbr_nd = np.zeros((nr, nc), bool)
    for dr, dc in _OFF:
        sl = np.zeros((nr, nc), bool)
        r0, r1 = max(0, -dr), nr - max(0, dr)
        c0, c1 = max(0, -dc), nc - max(0, dc)
        sl[r0:r1, c0:c1] = ndmask[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        nbr_nd |= sl
    seed = data & (pad | nbr_nd)

    heap = []
    cnt = 0
    sr, sc = np.where(seed)
    for r, c in zip(sr.tolist(), sc.tolist()):
        heapq.heappush(heap, (float(filled[r, c]), cnt, r, c)); cnt += 1
        state[r, c] = 1
    while heap:
        zc, _, row, col = heapq.heappop(heap)
        state[row, col] = 2
        zc = np.float32(zc)
        for k, (dr, dc) in enumerate(_OFF):
            i, j = row + dr, col + dc
            if i < 0 or i >= nr or j < 0 or j >= nc or state[i, j] != 0:
                continue
            if filled[i, j] <= zc:
                filled[i, j] = np.float32(zc + (inc_card if k % 2 == 0 else inc_diag)) if min_slope > 0 else zc
            heapq.heappush(heap, (float(filled[i, j]), cnt, i, j)); cnt += 1
            state[i, j] = 1
    return filled


# ------------------------------------------------------- tangential curvature
def tangential_curvature(dem, nodata=-9999.0, cellsize=1.0, window_radius=7):
    """2nd-order polynomial-fit tangential curvature on a (filled) DEM
    (LSDRaster::calculate_polyfit_surface_metrics, raster_selection[6]).
    Vectorised: the constant normal-equation matrix makes each coefficient field a
    correlation of the DEM with a fixed kernel."""
    from scipy import ndimage
    ndf = np.float32(nodata)
    z = dem.astype(np.float64)
    nr, nc = z.shape
    valid = z != ndf
    kr = int(np.ceil(window_radius / cellsize))
    kw = 2 * kr + 1
    ii = np.arange(kw)
    xk = (ii[:, None] - kr) * cellsize * np.ones((1, kw))     # x along rows
    yk = np.ones((kw, 1)) * (ii[None, :] - kr) * cellsize     # y along cols
    mask = (np.floor(np.sqrt(xk * xk + yk * yk)) <= window_radius).astype(np.float64)

    x = xk[mask == 1]; y = yk[mask == 1]
    basis = [x ** 2, y ** 2, x * y, x, y, np.ones_like(x)]
    A = np.array([[np.sum(basis[r] * basis[c]) for c in range(6)] for r in range(6)])
    if int(mask.sum()) < 6 or np.linalg.matrix_rank(A) < 6:
        raise ValueError(
            "tangential_curvature: curvature window too small for the cell size "
            "(window_radius=%g, cellsize=%g -> only %d cells inside the window; "
            "the 2nd-order polyfit needs >= 6). The window must reach the "
            "neighbours: use window_radius >= ~1.5*cellsize."
            % (window_radius, cellsize, int(mask.sum())))
    Ainv = np.linalg.inv(A)
    kernels = []
    for kvals in basis:
        ker = np.zeros((kw, kw)); ker[mask == 1] = kvals; kernels.append(ker)
    bb = np.stack([ndimage.correlate(z, ker, mode='constant', cval=0.0) for ker in kernels])
    a, b, c, d, e = np.tensordot(Ainv[:5], bb, axes=([1], [0]))   # x^2 y^2 xy x y coeffs

    fx, fy, fxx, fyy, fxy = d, e, 2 * a, 2 * b, c
    p = fx * fx + fy * fy
    q = p + 1.0
    denom = p * np.sqrt(q)
    with np.errstate(invalid='ignore', divide='ignore'):
        tanc = (fxx * fy * fy - 2 * fxy * fx * fy + fyy * fx * fx) / denom

    interior = np.zeros((nr, nc), bool); interior[kr:nr - kr, kr:nc - kr] = True
    nodata_in_square = ndimage.maximum_filter((~valid).astype(np.uint8), size=kw) > 0
    good = interior & valid & ~nodata_in_square & (q > 0) & (denom != 0)
    out = np.full((nr, nc), ndf, dtype=np.float32)
    out[good] = tanc[good].astype(np.float32)
    return out


# ---------------------------------------------------------------- D8 FlowInfo
def build_flowinfo(dem, nodata=-9999.0, cellsize=1.0):
    """D8 LSDFlowInfo::create. Returns a dict of node-indexed structures."""
    ndf = np.float32(nodata)
    z = dem.astype(np.float32)
    nr, nc = z.shape
    valid = z != ndf
    NEG = np.float32(-1e30)
    slopes = np.full((8, nr, nc), NEG, dtype=np.float32)
    for dn, (dr, dc) in enumerate(_OFF):
        fac = np.float32(1.0) if dn % 2 == 0 else ROOT2INV
        nb = np.full((nr, nc), ndf, dtype=np.float32)
        r0, r1 = max(0, -dr), nr - max(0, dr)
        c0, c1 = max(0, -dc), nc - max(0, dc)
        nb[r0:r1, c0:c1] = z[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        drop = (z - nb).astype(np.float32) * fac
        slopes[dn] = np.where(valid & (nb != ndf), drop, NEG)
    best_dir = np.argmax(slopes, axis=0)
    best_slope = np.max(slopes, axis=0)
    has = valid & (best_slope > np.float32(0.0))

    NodeIndex = np.full((nr, nc), -1, dtype=np.int64)
    rr, cc = np.where(valid)
    N = rr.size
    NodeIndex[rr, cc] = np.arange(N)
    drv = np.array([o[0] for o in _OFF]); dcv = np.array([o[1] for o in _OFF])
    recv_r = np.where(has, np.arange(nr)[:, None] + drv[best_dir],
                      np.arange(nr)[:, None] * np.ones((1, nc), int))
    recv_c = np.where(has, np.arange(nc)[None, :] + dcv[best_dir],
                      np.ones((nr, 1), int) * np.arange(nc)[None, :])
    recv = NodeIndex[recv_r[rr, cc], recv_c[rr, cc]]
    even = (best_dir % 2 == 0)
    flc = np.where(has[rr, cc], np.where(even[rr, cc], 1, 2), 0).astype(np.int8)

    ndon = np.zeros(N, dtype=np.int64)
    np.add.at(ndon, recv, 1)
    delta = np.zeros(N + 1, dtype=np.int64); delta[1:] = np.cumsum(ndon)
    donorstack = np.argsort(recv, kind='stable')

    return dict(nr=nr, nc=nc, ndf=ndf, res=float(cellsize), N=N,
                NodeIndex=NodeIndex, row_of=rr.astype(np.int64), col_of=cc.astype(np.int64),
                recv=recv, flc=flc, ndon=ndon, delta=delta, donorstack=donorstack,
                baselevel=np.where(recv == np.arange(N))[0])


def build_flowinfo_from_directions(direction, dem, nodata=-9999.0, cellsize=1.0):
    """Build a FlowInfo from an external D8 drainage-direction raster instead of
    computing steepest descent internally. ``direction`` uses the r.watershed
    encoding (1..8, CCW from NE; negative = drains off the region edge; 0/NaN =
    no flow). Only ``recv`` / ``flc`` come from the directions; node indexing and
    the donor stack are built exactly as :func:`build_flowinfo`, so the result is
    a drop-in FlowInfo for the rest of the DrEICH pipeline.

    A node's receiver is the in-region, non-no-data neighbour its direction points
    to; cells whose direction is off-grid, off-region, no-data, or unset become
    base levels (self-receiver). ``dem`` supplies the grid shape and the valid /
    no-data mask (use the same surface the rest of the pipeline runs on).
    """
    z = dem.astype(np.float32)
    nr, nc = z.shape
    valid = z != np.float32(nodata)
    NodeIndex = np.full((nr, nc), -1, dtype=np.int64)
    rr, cc = np.where(valid)
    N = rr.size
    NodeIndex[rr, cc] = np.arange(N)

    d = np.asarray(direction)[rr, cc]
    d = np.where(np.isfinite(d), d, 0)                  # NaN/NULL -> no flow
    a = np.abs(np.rint(d)).astype(np.int64)
    a = np.where((a >= 1) & (a <= 8), a, 0)             # out-of-range -> no flow
    tr = rr + _RWS_DROW[a]
    tc = cc + _RWS_DCOL[a]
    ingrid = (a > 0) & (tr >= 0) & (tr < nr) & (tc >= 0) & (tc < nc)
    trc = np.clip(tr, 0, nr - 1)
    tcc = np.clip(tc, 0, nc - 1)
    routes = ingrid & valid[trc, tcc]                  # receiver is a real cell
    recv = np.where(routes, NodeIndex[trc, tcc], np.arange(N))   # else self (outlet)
    flc = np.where(routes, np.where(np.isin(a, list(_RWS_CARDINAL)), 1, 2),
                   0).astype(np.int8)

    ndon = np.zeros(N, dtype=np.int64)
    np.add.at(ndon, recv, 1)
    delta = np.zeros(N + 1, dtype=np.int64)
    delta[1:] = np.cumsum(ndon)
    donorstack = np.argsort(recv, kind='stable')
    return dict(nr=nr, nc=nc, ndf=np.float32(nodata), res=float(cellsize), N=N,
                NodeIndex=NodeIndex, row_of=rr.astype(np.int64),
                col_of=cc.astype(np.int64), recv=recv, flc=flc, ndon=ndon,
                delta=delta, donorstack=donorstack,
                baselevel=np.where(recv == np.arange(N))[0])


def directions_from_flowinfo(fi):
    """Encode FlowInfo receivers as an r.watershed drainage-direction raster
    (1..8 CCW from NE; 0 = baselevel/outlet). Inverse of
    :func:`build_flowinfo_from_directions` -- for emitting a routing raster
    (e.g. ``r.fluvial.fastscape direction=``) that ``r.fluvial.channelheads
    direction=`` or ``r.stream.distance`` can consume.

    Returns an int32 (nr, nc) array: 1..8 where flow leaves the cell, 0 at base
    levels, and -1 at non-node (no-data) cells (use -1 as the raster NULL value)."""
    nr, nc = fi['nr'], fi['nc']
    row, col, recv = fi['row_of'], fi['col_of'], fi['recv']
    dr = row[recv] - row
    dc = col[recv] - col
    codes = np.zeros(fi['N'], dtype=np.int32)        # 0 = baselevel (recv == self)
    for (ddr, ddc), cd in _RWS_CODE.items():
        codes[(dr == ddr) & (dc == ddc)] = cd
    out = np.full((nr, nc), -1, dtype=np.int32)      # -1 = no-data / NULL
    out[row, col] = codes
    return out


def contributing_area(fi):
    """NContributingNodes via vectorised reverse-topological waves."""
    N = fi['N']; recv = fi['recv']
    self_loop = recv == np.arange(N)
    indeg = np.zeros(N, dtype=np.int64)
    nonself = np.where(~self_loop)[0]
    np.add.at(indeg, recv[nonself], 1)
    area = np.ones(N, dtype=np.int64)
    active = np.where((indeg == 0) & ~self_loop)[0]
    while active.size:
        r = recv[active]
        np.add.at(area, r, area[active])
        np.add.at(indeg, r, -1)
        nxt = r[indeg[r] == 0]
        active = np.unique(nxt[recv[nxt] != nxt])
    fi['ncontrib'] = area
    return area


def build_svector(fi):
    """Braun & Willett ordered stack (DFS from base level, donors in DonorStack
    order) + SVectorIndex. Iterative; hot arrays as Python lists for speed."""
    N = fi['N']; delta = fi['delta']
    ds = fi['donorstack'].copy()
    for k in fi['baselevel']:
        b0 = delta[k]
        if ds[b0] != k:
            seg = ds[delta[k]:delta[k + 1]]
            pos = np.where(seg == k)[0]
            if pos.size:
                p = b0 + int(pos[0]); ds[p] = ds[b0]; ds[b0] = k
    ds = ds.tolist(); delta_l = delta.tolist()
    SVector = np.empty(N, dtype=np.int64)
    j = 0
    stack = [int(x) for x in fi['baselevel'][::-1]]
    while stack:
        node = stack.pop()
        SVector[j] = node; j += 1
        for m in range(delta_l[node + 1] - 1, delta_l[node] - 1, -1):
            d = ds[m]
            if d != node:
                stack.append(d)
    SVectorIndex = np.empty(N, dtype=np.int64); SVectorIndex[SVector] = np.arange(N)
    fi['SVector'] = SVector; fi['SVectorIndex'] = SVectorIndex
    return fi


def distance_from_outlet(fi):
    """Flow distance to outlet (float32, as LSDTT Array2D<float>), forward BFS waves."""
    N = fi['N']; recv = fi['recv']; flc = fi['flc']
    delta = fi['delta']; donorstack = fi['donorstack']
    diag = np.float32(np.sqrt(2.0) * fi['res'])
    fd = np.full(N, -1.0, dtype=np.float32)
    fd[fi['baselevel']] = np.float32(0.0)
    step = np.where(flc == 2, diag, np.float32(fi['res'])).astype(np.float32)
    active = fi['baselevel']
    while active.size:
        lens = delta[active + 1] - delta[active]
        if lens.sum() == 0:
            break
        idx = np.concatenate([np.arange(s, s + l) for s, l in zip(delta[active], lens)])
        donors = donorstack[idx]
        donors = donors[donors != recv[donors]]
        fd[donors] = fd[recv[donors]] + step[donors]
        active = donors
    fi['fd'] = fd
    return fd


def get_sources(fi, threshold):
    """Channel cells (ncontrib >= threshold) with no channel donor are sources."""
    N = fi['N']; recv = fi['recv']
    chan = fi['ncontrib'] >= threshold
    donor_is_chan = np.zeros(N, dtype=bool)
    cf = np.where(chan)[0]
    nonself = cf[recv[cf] != cf]
    np.add.at(donor_is_chan, recv[nonself], True)
    return np.where(chan & ~donor_is_chan)[0]


# ------------------------------------------------------------- junction network
def junction_network(fi, sources):
    """LSDJunctionNetwork::create -> StreamOrderArray, JunctionArray (confluence
    markers), JunctionIndexArray (node->junction), JunctionVector (junction->node),
    all node-indexed."""
    N = fi['N']; recv = fi['recv']; delta = fi['delta']; donorstack = fi['donorstack']
    ND = -9999
    SO = np.full(N, ND, dtype=np.int64)
    JArr = np.full(N, ND, dtype=np.int64)
    sources = list(sources)

    for s in sources:                                       # first loop: stream order + confluences
        baselevel_switch = 0; junction_switch = 0
        cur = s; receiver = recv[cur]; cso = 1
        if cur == receiver:
            baselevel_switch = 1
        while baselevel_switch < 2:
            if SO[cur] == ND:
                SO[cur] = cso
            else:
                if junction_switch == 0:
                    junction_switch = 1; JArr[cur] = 1
                    if SO[cur] == cso:
                        cso += 1; SO[cur] = cso
                    else:
                        baselevel_switch = 2
                else:
                    if JArr[cur] != 1:
                        if cso > SO[cur]:
                            SO[cur] = cso
                        else:
                            baselevel_switch = 2
                    else:
                        if SO[cur] > cso:
                            baselevel_switch = 2
                        elif SO[cur] == cso:
                            ndon_same = 0
                            for m in range(delta[cur], delta[cur + 1]):
                                dn = donorstack[m]
                                if dn != cur and SO[dn] == cso:
                                    ndon_same += 1
                            if ndon_same >= 2:
                                cso += 1; SO[cur] = cso
                            else:
                                baselevel_switch = 2
                        else:
                            SO[cur] = cso
            cur = recv[cur]; receiver = recv[cur]
            if cur == receiver:
                baselevel_switch += 1

    JIdx = np.full(N, ND, dtype=np.int64)
    JunctionVector = []
    this_junction = -1
    for s in sources:                                       # second loop: junction indices
        this_junction += 1
        baselevel_switch = 0; junction_switch = 0
        cur = s; receiver = recv[cur]
        JunctionVector.append(cur); JIdx[cur] = this_junction
        if receiver == cur:
            baselevel_switch = 1
        while baselevel_switch == 0 and junction_switch < 2:
            cur = receiver; receiver = recv[cur]
            if cur == receiver:
                if JIdx[cur] == ND:
                    this_junction += 1; JIdx[cur] = this_junction; JunctionVector.append(cur)
                junction_switch = 2; baselevel_switch = 1
            else:
                if JArr[cur] != ND:
                    junction_switch = JArr[cur]; JArr[cur] += 1
                    if JIdx[cur] == ND:
                        this_junction += 1; JIdx[cur] = this_junction; JunctionVector.append(cur)
    fi['SO'] = SO; fi['JArr'] = JArr; fi['JIdx'] = JIdx
    fi['JunctionVector'] = np.array(JunctionVector, dtype=np.int64)
    return fi


def _find_upstream_junction(fi, channel_node):
    """find_upstream_junction_from_channel_nodeindex: walk up same-order donors to
    the first junction cell."""
    recv = fi['recv']; delta = fi['delta']; donorstack = fi['donorstack']
    SO = fi['SO']; JIdx = fi['JIdx']; ND = -9999
    cur = channel_node; tco = SO[cur]
    if tco == ND:
        return ND
    if JIdx[cur] != ND:
        return JIdx[cur]
    while True:
        b, e = delta[cur], delta[cur + 1]; n_don = e - b
        this_donor = 0; cr = cur
        while True:
            cr = donorstack[b + this_donor]
            this_donor += 1
            if not (SO[cr] != tco and this_donor < n_don):
                break
        if JIdx[cr] != ND:
            return JIdx[cr]
        cur = donorstack[b + this_donor - 1]


def find_valleys(fi, tan_curv, sources, n_connecting_nodes=10, tan_curv_threshold=0.1):
    """LSDJunctionNetwork::find_valleys. Node-indexed valley-junction array."""
    recv = fi['recv']; SO = fi['SO']; row_of = fi['row_of']; col_of = fi['col_of']
    N = fi['N']; ND = -9999; ndf = fi['ndf']
    visited = np.zeros(N, dtype=bool)
    valley = np.full(N, ND, dtype=np.int64)
    for s in sources:
        end = False; consec = 0; cur = s
        while not end:
            rcv = recv[cur]
            cv = tan_curv[row_of[cur], col_of[cur]]
            if cv != ndf:
                visited[cur] = True
                consec = consec + 1 if cv > tan_curv_threshold else 0
                if consec > n_connecting_nodes:
                    end = True
                    this_node = cur; seen = set(); outlet = False
                    while not outlet:
                        dn = recv[this_node]
                        if SO[dn] > SO[this_node]:
                            valley[this_node] = _find_upstream_junction(fi, this_node)
                            outlet = True
                        seen.add(this_node)
                        if dn in seen:
                            outlet = True
                        else:
                            this_node = dn
                if visited[rcv]:
                    end = True
                else:
                    cur = rcv
            else:
                end = True
    return valley


# ----------------------------------------------------------- chi-z head finder
def _find_farthest_upslope(fi, node):
    s = int(fi['SVectorIndex'][node]); n = int(fi['ncontrib'][node])
    upl = fi['SVector'][s:s + n]
    fdv = fi['fd'][upl]
    k = int(np.argmax(fdv))
    return int(upl[k]) if fdv[k] > 0.0 else node


def _build_channel_chi(fi, dem, hilltop, end_node, A_0, m_over_n):
    recv = fi['recv']; flc = fi['flc']; ncontrib = fi['ncontrib']
    row_of = fi['row_of']; col_of = fi['col_of']
    diag = np.float32(np.sqrt(2.0) * fi['res']); res = np.float32(fi['res'])
    pixel_area = fi['res'] * fi['res']
    nodeseq = [hilltop]; cur = hilltop
    while cur != end_node:
        r = recv[cur]; nodeseq.append(int(r))
        if r == cur:
            break
        cur = int(r)
    n = len(nodeseq)
    elev = np.empty(n, dtype=np.float32); chi = np.zeros(n, dtype=np.float32)
    for i, nd in enumerate(nodeseq):
        elev[i] = dem[row_of[nd], col_of[nd]]
    A_0 = np.float32(A_0)
    for i in range(n - 2, -1, -1):
        nd = nodeseq[i]
        dx = diag if flc[nd] == 2 else res
        area = np.float32(np.float32(ncontrib[nd]) * np.float32(pixel_area))
        ratio = np.float32(A_0 / area)
        chi[i] = np.float32(float(dx) * (float(ratio) ** m_over_n) + float(chi[i + 1]))
    return nodeseq, chi, elev


def _lu2_f32(a00, a01, a10, a11, b0, b1):
    """TNT JAMA LU solve of a 2x2 system in float32 (partial pivoting on column 0)."""
    f32 = np.float32
    if abs(a10) > abs(a00):
        a00, a10 = a10, a00; a01, a11 = a11, a01; b0, b1 = b1, b0
    mult = f32(a10 / a00)
    u11 = f32(a11 - f32(mult * a01))
    y1 = f32(b1 - f32(b0 * mult))
    x1 = f32(y1 / u11)
    x0 = f32(f32(b0 - f32(x1 * a01)) / a00)
    return x0, x1


def _ssum(arr):
    """Sequential left-to-right float32 sum (== C++ TNT matmult / stat-loop order)."""
    if arr.size == 0:
        return np.float32(0.0)
    return np.cumsum(arr, dtype=np.float32)[-1]


def _reg32(x, y):
    """float32 simple_linear_regression -> (R2, Durbin-Watson)."""
    f32 = np.float32
    x = x.astype(np.float32, copy=False); y = y.astype(np.float32, copy=False)
    n = f32(x.size)
    Sxx = _ssum(x * x); Sx = _ssum(x); Sxy = _ssum(x * y); Sy = _ssum(y)
    m, b = _lu2_f32(Sxx, Sx, Sx, n, Sxy, Sy)
    mean = f32(Sy / n)
    SST = _ssum(((y - mean) ** 2).astype(np.float32))
    pred = (m * x + b).astype(np.float32)
    resid = (pred - y).astype(np.float32)
    resid[np.abs(resid) < f32(1e-12)] = f32(0.0)
    SSerr = _ssum((resid * resid).astype(np.float32))
    r2 = f32(f32(1.0) - SSerr / SST)
    dresid = np.diff(resid).astype(np.float32)
    top = _ssum((dresid * dresid).astype(np.float32))
    bottom = _ssum((resid * resid).astype(np.float32))
    if bottom == f32(0.0):
        bottom = f32(1e-10)
    return r2, f32(top / bottom)


def _calculate_channel_head(nodeseq, chi, elev, min_segment_length):
    """LSDChannel::calculate_channel_heads chi-z split -> head node."""
    f32 = np.float32
    end = len(chi); max_test = f32(0.0); elev_inter = f32(0.0)
    for hill in range(min_segment_length, end - min_segment_length + 1):
        r2_chan, _ = _reg32(chi[hill:], elev[hill:])
        _, dw_hill = _reg32(chi[:hill], elev[:hill])
        test = f32(r2_chan - f32((dw_hill - f32(2.0)) / f32(2.0)))
        if test > max_test:
            max_test = test; elev_inter = elev[hill]
    node_index = 0
    for i in range(len(elev)):
        if elev[i] == elev_inter:
            node_index = nodeseq[i]
    return node_index


def channel_heads_from_valleys(fi, dem, valley, min_segment_length=10,
                               A_0=1000.0, m_over_n=0.525):
    """GetChannelHeadsChiMethodFromValleys: per valley junction -> head, then
    furthest-upstream dedup. Returns (heads_per_junction, final_heads) node lists."""
    row_of = fi['row_of']; col_of = fi['col_of']; recv = fi['recv']
    JunctionVector = fi['JunctionVector']; ND = -9999
    nz = np.where(valley != ND)[0]
    order = np.lexsort((col_of[nz], row_of[nz]))            # row-major junction_list
    junction_list = valley[nz[order]]
    heads_temp = []
    for jn in junction_list:
        dn = int(JunctionVector[jn])
        hilltop = _find_farthest_upslope(fi, dn)
        nodeseq, chi, elev = _build_channel_chi(fi, dem, hilltop, dn, A_0, m_over_n)
        heads_temp.append(_calculate_channel_head(nodeseq, chi, elev, min_segment_length))
    heads_temp = np.array(heads_temp, dtype=np.int64)

    elevs = dem[row_of[heads_temp], col_of[heads_temp]].astype(np.float64)
    sorted_heads = heads_temp[np.argsort(-elevs, kind='stable')]   # highest elevation first
    binary = {int(h): 1 for h in heads_temp}
    visited = set()
    for h in sorted_heads:
        cur = int(h); end = False
        while not end:
            visited.add(cur)
            r = int(recv[cur])
            if binary.get(r, 0) == 1:
                binary[r] = 0
            cur = r
            if r in visited:
                end = True
    final = [h for h in binary if binary[h] == 1]
    return heads_temp, final


def drainage_divides(fi, threshold=100):
    """Trace drainage divides (ridgelines) as basin boundaries from a FlowInfo.

    Channels are the cells whose contributing area >= ``threshold``; each cell is
    labelled by the channel LINK it ultimately drains into (the network split at
    confluences, as in :func:`channel_network_segments`); a divide cell is one
    that is 8-adjacent to a cell draining to a *different* link. Divides are the
    topological dual of the channel network -- the ridges that wrap each
    first-order valley -- and (being basin boundaries) are far more robust to
    coarsening than the convergent valley-head feature itself.

    Returns ``(divide, basin_label)``, both (nr, nc): ``divide`` is boolean
    (ridge cells), ``basin_label`` is the per-cell draining-link id (-1 where the
    cell drains nowhere on the network)."""
    if 'ncontrib' not in fi:
        contributing_area(fi)
    nr, nc, N = fi['nr'], fi['nc'], fi['N']
    recv, row, col, NodeIndex = fi['recv'], fi['row_of'], fi['col_of'], fi['NodeIndex']

    sources = get_sources(fi, threshold)
    segs = channel_network_segments(fi, [(int(row[s]), int(col[s])) for s in sources])
    link = np.full(N, -1, dtype=np.int64)
    for seg in segs:                                  # label channel cells by link
        cells = np.asarray(seg['cells'])
        link[NodeIndex[cells[:, 0], cells[:, 1]]] = seg['cat']

    # propagate each hillslope cell's draining-link upstream (downstream-first wave)
    label = link.copy()
    known = label >= 0
    idx = np.arange(N)
    while True:
        ready = (~known) & (recv != idx) & known[recv]
        if not ready.any():
            break
        label[ready] = label[recv[ready]]
        known[ready] = True

    lab = np.full((nr, nc), -1, dtype=np.int64)
    lab[row, col] = label
    divide = np.zeros((nr, nc), dtype=bool)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            r0, r1 = max(0, -dr), nr - max(0, dr)
            c0, c1 = max(0, -dc), nc - max(0, dc)
            a = lab[r0:r1, c0:c1]
            b = lab[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
            divide[r0:r1, c0:c1] |= (a != b) & (a >= 0) & (b >= 0)
    return divide, lab


# ---------------------------------------------------------------- top-level API
def extract_channel_heads(dem, nodata=-9999.0, cellsize=1.0, *, fill_dem=True,
                          threshold=100, min_slope=0.0001, A_0=1000.0, m_over_n=0.525,
                          n_connecting_nodes=10, min_segment_length=10,
                          window_radius=7, tan_curv_threshold=0.1, return_filled=False,
                          return_flowinfo=False, filled=None, curvature=None,
                          direction=None):
    """Extract DrEICH channel heads from a DEM array.

    Parameters
    ----------
    dem : 2-D ndarray
        Elevation, row 0 = north. ``nodata`` marks no-data cells.
    nodata, cellsize : float
        No-data sentinel and (square) cell size in map units.
    fill_dem : bool
        If True, depression-fill the DEM first (LSDTT priority-flood). Set False if
        ``dem`` is already filled.
    threshold : int
        Source drainage-area threshold (cells).
    min_slope, A_0, m_over_n, n_connecting_nodes, min_segment_length,
    window_radius, tan_curv_threshold :
        DrEICH parameters (LSDTT defaults: 0.0001, 1000, 0.525, 10, 10, 7, 0.1).
    filled, curvature, direction : optional ndarrays
        Inject a pre-filled DEM, a curvature field, and/or an external D8
        drainage-direction raster (r.watershed encoding) to replace the
        corresponding internal step. ``direction`` lets the routing come from
        another module (e.g. ``r.watershed -s`` or ``r.fluvial.fastscape``).

    Returns
    -------
    heads : list of (row, col)
        Channel-head cells. If ``return_filled``, the filled DEM is appended; if
        ``return_flowinfo``, the FlowInfo dict (D8 routing, for tracing the
        channel network downstream of the heads) is appended. The return is a
        tuple in the order ``(heads[, filled][, flowinfo])`` when either is set.
    """
    # ``filled``, ``curvature`` and ``direction`` may be injected to swap in
    # alternative components for comparison or to consume routing from another
    # module: ``direction`` is an external D8 drainage-direction raster
    # (r.watershed encoding) used instead of internal steepest descent. Curvature
    # is always faithful (still computed on the filled surface).
    if filled is None:
        filled = fill(dem, nodata, min_slope, cellsize) if fill_dem else dem.astype(np.float32)
    else:
        filled = filled.astype(np.float32)
    tcurv = curvature if curvature is not None else \
        tangential_curvature(filled, nodata, cellsize, window_radius)
    fi = (build_flowinfo_from_directions(direction, filled, nodata, cellsize)
          if direction is not None else build_flowinfo(filled, nodata, cellsize))
    contributing_area(fi)
    sources = get_sources(fi, threshold)
    junction_network(fi, sources)
    valley = find_valleys(fi, tcurv, sources, n_connecting_nodes, tan_curv_threshold)
    build_svector(fi)
    distance_from_outlet(fi)
    _, final = channel_heads_from_valleys(fi, filled, valley, min_segment_length, A_0, m_over_n)
    heads = [(int(fi['row_of'][n]), int(fi['col_of'][n])) for n in final]
    if not (return_filled or return_flowinfo):
        return heads
    out = [heads]
    if return_filled:
        out.append(filled)
    if return_flowinfo:
        out.append(fi)
    return tuple(out)


def channel_network_segments(fi, heads):
    """Trace the D8 channel network downstream of the channel ``heads`` and split
    it into stream links at confluences -- the v.stream.network directed-graph
    representation.

    The network is every cell on a flow path from a channel head down to its
    outlet (region edge or internal pit). It is cut into links at junctions:
    a link runs from a source or a confluence down to the next confluence or the
    outlet, with the confluence cell shared between the incoming links and the
    single outgoing link. Each link points to exactly one downstream link
    (``tostream``), so the network is a converging tree (``OFFMAP = 0`` at the
    outlet), matching what r.stream.extract + v.stream.network produce.

    Parameters
    ----------
    fi : dict
        FlowInfo from :func:`build_flowinfo` (uses ``recv``, ``NodeIndex``,
        ``row_of``, ``col_of``, ``N``, ``nr``, ``nc``).
    heads : list of (row, col)
        Channel-head cells (e.g. from :func:`extract_channel_heads`).

    Returns
    -------
    segments : list of dict
        One per link, ordered and ``cat``-numbered (1-based) by the (row, col)
        of its upstream end. Each dict has ``cat`` (int), ``cells`` (list of
        (row, col) from upstream junction to downstream junction, inclusive),
        and ``tostream`` (cat of the downstream link, or 0 at the outlet).
    """
    recv = fi['recv']; NodeIndex = fi['NodeIndex']
    row_of = fi['row_of']; col_of = fi['col_of']; N = fi['N']

    # 1. mark every node on a flow path downstream of a head. Stop as soon as a
    #    path merges into an already-marked one (the rest is already on).
    on = np.zeros(N, dtype=bool)
    for (r, c) in heads:
        node = int(NodeIndex[r, c])
        if node < 0:
            continue
        while not on[node]:
            on[node] = True
            nxt = int(recv[node])
            if nxt == node:                     # baselevel / internal pit
                break
            node = nxt
    chan = np.where(on)[0]

    # 2. channel in-degree: how many channel nodes drain straight into each node.
    rv = recv[chan]
    feeds = on[rv] & (rv != chan)               # exclude self-loops and off-network
    indeg = np.zeros(N, dtype=np.int64)
    np.add.at(indeg, rv[feeds], 1)

    def is_outlet(nd):
        nxt = int(recv[nd])
        return nxt == nd or not on[nxt]

    # 3. links begin at sources (indeg 0) and confluences (indeg >= 2).
    starts = [int(nd) for nd in chan if indeg[nd] == 0 or indeg[nd] >= 2]

    # 4. walk each start downstream to the next confluence or the outlet.
    raw = []
    for s in starts:
        cells = [s]; nd = s
        while True:
            if is_outlet(nd):
                end = nd; break
            nxt = int(recv[nd])
            cells.append(nxt); nd = nxt
            if indeg[nd] >= 2:                   # confluence: link ends (shared cell)
                end = nd; break
        if len(cells) >= 2:                      # drop a degenerate 1-cell outlet head
            raw.append((s, end, cells))

    # 5. cat-number by upstream (row, col), then resolve tostream from the end
    #    node: the link that *starts* at this link's end confluence, else 0.
    raw.sort(key=lambda t: (int(row_of[t[0]]), int(col_of[t[0]])))
    seg_of_start = {s: i + 1 for i, (s, _e, _c) in enumerate(raw)}
    segments = []
    for i, (_s, end, cells) in enumerate(raw):
        segments.append(dict(
            cat=i + 1,
            cells=[(int(row_of[n]), int(col_of[n])) for n in cells],
            tostream=seg_of_start.get(end, 0)))
    return segments
