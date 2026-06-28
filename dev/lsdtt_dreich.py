"""Faithful direct port of the LSDTT/DrEICH channel-head pipeline (Clubb et al.
2014), operating on LSDTT's filled DEM, validated stage-by-stage against the
instrumented C++ rasters in /tmp/dreich_algorithm.

Ports (from the standalone dreich_algorithm repo the driver #includes):
  LSDFlowInfo::create        -> node index, receivers, FlowLengthCode, donor stack
                                (NDonors/Delta/DonorStack), ordered stack SVector,
                                NContributingNodes
  LSDFlowInfo helpers        -> distance_from_outlet, find_farthest_upslope_node,
                                is_node_upstream (via SVector contiguity)
  LSDJunctionNetwork::create -> StreamOrderArray, JunctionArray, JunctionIndexArray,
                                JunctionVector
  find_upstream_junction_from_channel_nodeindex
  find_valleys               -> valley junctions          (validate vs _VJ = 768)
  LSDChannel + calculate_channel_heads (chi-z split)
  GetChannelHeadsChiMethodFromValleys -> channel heads     (validate vs _CH = 634)

Receiver/slope computation matches the already-EXACT lsdtt_replica_flowinfo.py
(float32, neighbour order N NE E SE S SW W NW, max_slope init 0 strict,
one_ov_root2 diagonal). Boundary conditions: No / no flux / no flux / No flux
(no base-level boundary), so a node is base level iff it has no downhill neighbour
(FlowLengthCode == 0).
"""
import sys
import numpy as np

from rivernetworkx.core import _linfit_r2_dw   # validated chi-z regression (R2, DW)

NR, NC = 7107, 6266
ND = -9999
NDF = np.float32(ND)
RES = 1.0
ROOT2INV = np.float32(0.707106781)
DIAG = np.float32(np.sqrt(2.0) * RES)
ALG = '/tmp/dreich_algorithm'

# neighbour order N NE E SE S SW W NW; even index = cardinal, odd = diagonal
OFF = [(-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)]


def build_flowinfo(z):
    """LSDFlowInfo::create on a float32 filled DEM. Returns a dict of node-indexed
    structures (everything the junction network + head finder need)."""
    valid = z != NDF
    NEG = np.float32(-1e30)
    slopes = np.full((8, NR, NC), NEG, dtype=np.float32)
    for d, (dr, dc) in enumerate(OFF):
        fac = np.float32(1.0) if d % 2 == 0 else ROOT2INV
        nb = np.full((NR, NC), NDF, dtype=np.float32)
        r0, r1 = max(0, -dr), NR - max(0, dr)
        c0, c1 = max(0, -dc), NC - max(0, dc)
        nb[r0:r1, c0:c1] = z[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        drop = (z - nb).astype(np.float32) * fac
        ok = valid & (nb != NDF)
        slopes[d] = np.where(ok, drop, NEG)
    best_dir = np.argmax(slopes, axis=0)              # first (lowest-index) max
    best_slope = np.max(slopes, axis=0)
    has = valid & (best_slope > np.float32(0.0))      # strict > 0 -> has a receiver

    # node indexing: row-major over valid cells (== LSDFlowInfo NodeIndex)
    NodeIndex = np.full((NR, NC), -1, dtype=np.int64)
    rr, cc = np.where(valid)
    N = rr.size
    NodeIndex[rr, cc] = np.arange(N)
    row_of = rr.astype(np.int64)
    col_of = cc.astype(np.int64)

    # receiver grid -> receiver node
    drv = np.array([o[0] for o in OFF]); dcv = np.array([o[1] for o in OFF])
    recv_r = np.where(has, np.arange(NR)[:, None] + drv[best_dir], np.arange(NR)[:, None] * np.ones((1, NC), int))
    recv_c = np.where(has, np.arange(NC)[None, :] + dcv[best_dir], np.ones((NR, 1), int) * np.arange(NC)[None, :])
    recv = NodeIndex[recv_r[rr, cc], recv_c[rr, cc]]    # node index of receiver
    # FlowLengthCode: 0 base level, 1 cardinal, 2 diagonal
    flc_grid_even = (best_dir % 2 == 0)
    flc = np.where(has[rr, cc], np.where(flc_grid_even[rr, cc], 1, 2), 0).astype(np.int8)

    # donor stack (Braun & Willett): NDonors, Delta, DonorStack
    ndon = np.zeros(N, dtype=np.int64)
    np.add.at(ndon, recv, 1)                            # every node donates once (self for base level)
    delta = np.zeros(N + 1, dtype=np.int64)
    delta[1:] = np.cumsum(ndon)
    donorstack = np.argsort(recv, kind='stable')        # nodes bucketed by receiver, ascending node within bucket

    return dict(valid=valid, N=N, NodeIndex=NodeIndex, row_of=row_of, col_of=col_of,
                recv=recv, flc=flc, ndon=ndon, delta=delta, donorstack=donorstack,
                baselevel=np.where(recv == np.arange(N))[0])


def contributing_area(fi):
    """NContributingNodes (subtree size in the receiver forest) via vectorized
    reverse-topological waves -- O(N) total, no per-node Python loop. Equivalent
    to the SVector accumulation; used where the exact SVector order is not needed
    (e.g. selecting sources)."""
    N = fi['N']; recv = fi['recv']
    self_loop = recv == np.arange(N)
    indeg = np.zeros(N, dtype=np.int64)            # number of non-self donors
    nonself = np.where(~self_loop)[0]
    np.add.at(indeg, recv[nonself], 1)
    area = np.ones(N, dtype=np.int64)
    active = np.where((indeg == 0) & ~self_loop)[0]   # leaves (exclude base-level self-loops handled below)
    # base-level leaves (self-loop with no donors) just keep area 1; they never push
    while active.size:
        r = recv[active]
        np.add.at(area, r, area[active])
        np.add.at(indeg, r, -1)
        nxt = r[indeg[r] == 0]
        nxt = nxt[recv[nxt] != nxt]                   # don't propagate through base level
        active = np.unique(nxt)
    fi['ncontrib'] = area
    return area


def build_svector(fi):
    """Ordered stack SVector (DFS from base level, donors in DonorStack order),
    SVectorIndex. Iterative explicit-stack DFS; hot arrays as Python lists for
    speed over ~35M nodes."""
    N = fi['N']; delta = fi['delta']
    # base-level-first swap inside each base level node's donor bucket (create() l.633)
    ds = fi['donorstack'].copy()
    for k in fi['baselevel']:
        b0 = delta[k]
        if ds[b0] != k:
            seg = ds[delta[k]:delta[k + 1]]
            pos = np.where(seg == k)[0]
            if pos.size:
                p = b0 + int(pos[0])
                ds[p] = ds[b0]; ds[b0] = k
    ds = ds.tolist()
    delta_l = delta.tolist()
    SVector = np.empty(N, dtype=np.int64)
    sv = SVector
    j = 0
    stack = [int(x) for x in fi['baselevel'][::-1]]
    push = stack.append; pop = stack.pop
    while stack:
        node = pop()
        sv[j] = node; j += 1
        b = delta_l[node]; e = delta_l[node + 1]
        for m in range(e - 1, b - 1, -1):           # reverse so donors pop in stack order
            d = ds[m]
            if d != node:                            # skip base-level self-donor
                push(d)
    assert j == N, (j, N)
    SVectorIndex = np.empty(N, dtype=np.int64)
    SVectorIndex[SVector] = np.arange(N)
    fi['SVector'] = SVector; fi['SVectorIndex'] = SVectorIndex
    return fi


def distance_from_outlet(fi):
    """Flow distance to outlet, forward BFS waves from base level over the donor
    structure (parent computed before child). O(N), no per-node Python loop."""
    # C++ flow_distance is Array2D<float> (float32); mirror that exactly so
    # find_farthest_upslope_node's argmax breaks near-ties the same way.
    N = fi['N']; recv = fi['recv']; flc = fi['flc']
    delta = fi['delta']; donorstack = fi['donorstack']
    fd = np.full(N, -1.0, dtype=np.float32)
    fd[fi['baselevel']] = np.float32(0.0)
    step = np.where(flc == 2, DIAG, np.float32(RES)).astype(np.float32)
    active = fi['baselevel']
    while active.size:
        # donors of all active nodes (exclude self-donors)
        lens = delta[active + 1] - delta[active]
        if lens.sum() == 0:
            break
        starts = delta[active]
        idx = np.concatenate([np.arange(s, s + l) for s, l in zip(starts, lens)]) if active.size else np.array([], int)
        donors = donorstack[idx]
        donors = donors[donors != recv[donors]]      # drop base-level self-donors
        # set fd for donors (each donor has exactly one receiver, already fd-set)
        fd[donors] = fd[recv[donors]] + step[donors]
        active = donors
    fi['fd'] = fd
    return fd


def get_sources_index_threshold(fi, threshold):
    """Channel cells = ncontrib >= threshold; a source is a channel cell with no
    channel donor (matches lsdtt_replica_flowinfo, validated EXACT = 112757)."""
    N = fi['N']; recv = fi['recv']; ncontrib = fi['ncontrib']
    chan = ncontrib >= threshold
    donor_is_chan = np.zeros(N, dtype=bool)
    cf = np.where(chan)[0]
    # mark receivers that have a channel donor (exclude base-level self)
    nonself = cf[recv[cf] != cf]
    np.add.at(donor_is_chan, recv[nonself], True)
    src = np.where(chan & ~donor_is_chan)[0]
    return src


def junction_network(fi, sources):
    """LSDJunctionNetwork::create -> StreamOrderArray (node), JunctionArray (node),
    JunctionIndexArray (node-> junction or ND), JunctionVector (junction-> node)."""
    N = fi['N']; recv = fi['recv']; delta = fi['delta']; donorstack = fi['donorstack']
    SO = np.full(N, ND, dtype=np.int64)        # StreamOrderArray (node-indexed)
    JArr = np.full(N, ND, dtype=np.int64)      # JunctionArray (confluence markers)
    sources = list(sources)

    # ---- first loop: stream orders + junction (confluence) markers ----
    for s in sources:
        baselevel_switch = 0
        junction_switch = 0
        cur = s
        receiver = recv[cur]
        cso = 1
        if cur == receiver:
            baselevel_switch = 1
        while baselevel_switch < 2:
            if SO[cur] == ND:
                SO[cur] = cso
            else:
                if junction_switch == 0:
                    junction_switch = 1
                    JArr[cur] = 1
                    if SO[cur] == cso:
                        cso += 1
                        SO[cur] = cso
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
                                cso += 1
                                SO[cur] = cso
                            else:
                                baselevel_switch = 2
                        else:   # SO[cur] < cso
                            SO[cur] = cso
            cur = recv[cur]
            receiver = recv[cur]
            if cur == receiver:
                baselevel_switch += 1

    # ---- second loop: junction indices + JunctionVector ----
    JIdx = np.full(N, ND, dtype=np.int64)      # JunctionIndexArray
    JunctionVector = []
    this_junction = -1
    for s in sources:
        this_junction += 1
        baselevel_switch = 0
        junction_switch = 0
        cur = s
        receiver = recv[cur]
        JunctionVector.append(cur)             # each source is a junction
        JIdx[cur] = this_junction
        if receiver == cur:
            baselevel_switch = 1
        while baselevel_switch == 0 and junction_switch < 2:
            cur = receiver
            receiver = recv[cur]
            if cur == receiver:                # base level
                if JIdx[cur] == ND:
                    this_junction += 1
                    JIdx[cur] = this_junction
                    JunctionVector.append(cur)
                junction_switch = 2
                baselevel_switch = 1
            else:
                if JArr[cur] != ND:            # this is a confluence
                    junction_switch = JArr[cur]
                    JArr[cur] += 1
                    if JIdx[cur] != ND:
                        pass                   # already-visited junction: stop
                    else:
                        this_junction += 1
                        JIdx[cur] = this_junction
                        JunctionVector.append(cur)
    fi['SO'] = SO; fi['JArr'] = JArr; fi['JIdx'] = JIdx
    fi['JunctionVector'] = np.array(JunctionVector, dtype=np.int64)
    return fi


def find_upstream_junction(fi, channel_node):
    """find_upstream_junction_from_channel_nodeindex: walk up same-order donors to
    the first junction cell."""
    recv = fi['recv']; delta = fi['delta']; donorstack = fi['donorstack']
    SO = fi['SO']; JIdx = fi['JIdx']
    cur = channel_node
    tco = SO[cur]
    if tco == ND:
        return ND
    if JIdx[cur] != ND:
        return JIdx[cur]
    while True:
        b, e = delta[cur], delta[cur + 1]
        n_don = e - b
        this_donor = 0
        dco = ND
        cr = cur
        while True:
            cr = donorstack[b + this_donor]
            dco = SO[cr]
            this_donor += 1
            if not (dco != tco and this_donor < n_don):
                break
        if JIdx[cr] != ND:
            return JIdx[cr]
        cur = donorstack[b + this_donor - 1]


def find_valleys(fi, tan_curv_grid, sources, no_connecting_nodes=10, tan_curv_threshold=0.1):
    """LSDJunctionNetwork::find_valleys. Returns a node-indexed array: valley
    junction index at the valley-outlet cell, ND elsewhere."""
    N = fi['N']; recv = fi['recv']; SO = fi['SO']
    row_of = fi['row_of']; col_of = fi['col_of']
    tcurv = tan_curv_grid
    visited = np.zeros(N, dtype=bool)
    valley_junc = np.full(N, ND, dtype=np.int64)

    for s in sources:
        EndofReach = False
        consec = 0
        cur = s
        while not EndofReach:
            rcv = recv[cur]
            cv = tcurv[row_of[cur], col_of[cur]]
            if cv != NDF:
                visited[cur] = True
                if cv > tan_curv_threshold:
                    consec += 1
                else:
                    consec = 0
                if consec > no_connecting_nodes:
                    EndofReach = True
                    this_node = cur
                    visited_tmp = set()
                    reached_outlet = False
                    while not reached_outlet:
                        dn = recv[this_node]
                        cur_SO = SO[this_node]
                        dn_SO = SO[dn]
                        visited_tmp.add(this_node)
                        if dn_SO > cur_SO:
                            vj = find_upstream_junction(fi, this_node)
                            valley_junc[this_node] = vj
                            reached_outlet = True
                        if dn in visited_tmp:
                            reached_outlet = True
                        else:
                            this_node = dn
                # receiver-visited test
                if visited[rcv]:
                    EndofReach = True
                else:
                    cur = rcv
            else:
                EndofReach = True
    return valley_junc


def find_farthest_upslope_node(fi, node):
    """find_farthest_upslope_node: the upslope node with greatest flow distance,
    first in SVector order on ties (C++ inits farthest=0.0, strict >)."""
    s = int(fi['SVectorIndex'][node]); n = int(fi['ncontrib'][node])
    upl = fi['SVector'][s:s + n]
    fdv = fi['fd'][upl]
    k = int(np.argmax(fdv))                          # first max in SVector order
    if fdv[k] > 0.0:
        return int(upl[k])
    return node


def build_channel_chi(fi, z, hilltop, end_node, A_0, m_over_n):
    """LSDChannel(hilltop, end_node): node path hilltop->end via receivers, with chi
    (downslope_chi=0) computed inline as in LSDChannel::create_LSDC."""
    recv = fi['recv']; flc = fi['flc']; ncontrib = fi['ncontrib']
    row_of = fi['row_of']; col_of = fi['col_of']
    pixel_area = RES * RES
    nodeseq = [hilltop]
    cur = hilltop
    while cur != end_node:
        r = recv[cur]
        nodeseq.append(int(r))
        if r == cur:
            break
        cur = int(r)
    n = len(nodeseq)
    # elevation is read from a float32 raster; chi is stored float32 like C++.
    elev = np.empty(n, dtype=np.float32); chi = np.zeros(n, dtype=np.float32)
    for i, nd in enumerate(nodeseq):
        elev[i] = z[row_of[nd], col_of[nd]]
    A_0 = np.float32(A_0)
    # chi from downstream (last) upward, mirroring C++ float storage with double pow
    for i in range(n - 2, -1, -1):
        nd = nodeseq[i]
        dx = DIAG if flc[nd] == 2 else np.float32(RES)
        area = np.float32(np.float32(ncontrib[nd]) * np.float32(pixel_area))
        ratio = np.float32(A_0 / area)                 # float32 division
        val = float(dx) * (float(ratio) ** m_over_n) + float(chi[i + 1])
        chi[i] = np.float32(val)
    return nodeseq, chi, elev


def _lu2_f32(a00, a01, a10, a11, b0, b1):
    """TNT JAMA LU solve of a 2x2 system in float32 (partial pivoting on column 0),
    exactly as simple_linear_regression's LU<float>.solve. Returns (x0, x1)."""
    f32 = np.float32
    if abs(a10) > abs(a00):                            # pivot: swap rows
        a00, a10 = a10, a00
        a01, a11 = a11, a01
        b0, b1 = b1, b0
    mult = f32(a10 / a00)                              # L[1][0]
    u11 = f32(a11 - f32(mult * a01))                   # U[1][1]
    y1 = f32(b1 - f32(b0 * mult))                      # forward solve L y = Pb
    x1 = f32(y1 / u11)                                 # back solve U x = y
    x0 = f32(f32(b0 - f32(x1 * a01)) / a00)
    return x0, x1


def _ssum(arr):
    """Sequential left-to-right float32 sum (== C++ TNT matmult / stat-loop order).
    np.cumsum accumulates sequentially, unlike np.sum's pairwise summation."""
    if arr.size == 0:
        return np.float32(0.0)
    return np.cumsum(arr, dtype=np.float32)[-1]


def _reg32(x, y):
    """simple_linear_regression in float32, faithful to LSDStatsTools: sequential
    float32 summation (TNT matmult order) + JAMA LU 2x2 solve + R2 + Durbin-Watson."""
    f32 = np.float32
    x = x.astype(np.float32, copy=False); y = y.astype(np.float32, copy=False)
    n = f32(x.size)
    Sxx = _ssum(x * x)
    Sx = _ssum(x)
    Sxy = _ssum(x * y)
    Sy = _ssum(y)
    m, b = _lu2_f32(Sxx, Sx, Sx, n, Sxy, Sy)           # LHS = A^T A, RHS = A^T b
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
    dw = f32(top / bottom)
    return r2, dw


def calculate_channel_heads(nodeseq, chi, elev, min_seg):
    """LSDChannel::calculate_channel_heads chi-z split. Returns the head node.
    test = R2_channel - (DW_hillslope - 2)/2 over all splits; max (init 0, strict >);
    head node = NodeSequence at the (last) cell whose elevation == channel-front elev."""
    end = len(chi)
    f32 = np.float32
    max_test = f32(0.0)
    elev_inter = f32(0.0)
    for hill in range(min_seg, end - min_seg + 1):
        r2_chan, _ = _reg32(chi[hill:], elev[hill:])
        _, dw_hill = _reg32(chi[:hill], elev[:hill])
        test = f32(r2_chan - f32((dw_hill - f32(2.0)) / f32(2.0)))
        if test > max_test:                            # strict > -> first (lowest hill) wins on ties
            max_test = test
            elev_inter = elev[hill]                   # channel_elev.front()
    node_index = 0                                    # C++ init 0
    for i in range(len(elev)):
        if elev[i] == elev_inter:
            node_index = nodeseq[i]                   # last match wins (overwrite)
    return node_index


def get_channel_heads_from_valleys(fi, z, valley_junc, min_seg, A_0, m_over_n):
    """GetChannelHeadsChiMethodFromValleys: per valley junction -> hilltop ->
    chi-z head; then furthest-upstream dedup via the binary-map walk-down."""
    row_of = fi['row_of']; col_of = fi['col_of']; recv = fi['recv']
    JunctionVector = fi['JunctionVector']
    # junction_list in row-major order of the valley-outlet cells (C++ scans the grid)
    nz = np.where(valley_junc != ND)[0]
    order = np.lexsort((col_of[nz], row_of[nz]))      # sort by row, then col
    junction_list = valley_junc[nz[order]]

    heads_temp = []
    for jn in junction_list:
        downstream_node = int(JunctionVector[jn])
        hilltop = find_farthest_upslope_node(fi, downstream_node)
        nodeseq, chi, elev = build_channel_chi(fi, z, hilltop, downstream_node, A_0, m_over_n)
        head = calculate_channel_heads(nodeseq, chi, elev, min_seg)
        heads_temp.append(int(head))

    # ---- furthest-upstream dedup (binary_map walk-down, highest elevation first) ----
    heads_temp = np.array(heads_temp, dtype=np.int64)
    elevs = z[row_of[heads_temp], col_of[heads_temp]].astype(np.float64)
    # sort descending by elevation (stable -> ties keep input order, as matlab sort)
    order = np.argsort(-elevs, kind='stable')
    heads_sorted = heads_temp[order]
    binary = {int(h): 1 for h in heads_temp}          # binary_map at head cells
    visited = set()
    for h in heads_sorted:
        cur = int(h)
        EndOfReach = False
        while not EndOfReach:
            visited.add(cur)
            r = int(recv[cur])
            if binary.get(r, 0) == 1:
                binary[r] = 0
            cur = r
            if r in visited:
                EndOfReach = True
    final = [h for h in binary if binary[h] == 1]
    return heads_temp, final


# =========================================================================
if __name__ == '__main__':
    stage = sys.argv[1] if len(sys.argv) > 1 else 'vj'
    print('loading filled DEM + curvature ...', flush=True)
    z = np.fromfile(f'{ALG}/bailey_run_dem_fill.flt', dtype='<f4').reshape(NR, NC)
    tcurv = np.fromfile(f'{ALG}/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(NR, NC)

    fi = build_flowinfo(z)
    print('N data nodes: %d ; base level nodes: %d' % (fi['N'], fi['baselevel'].size), flush=True)
    contributing_area(fi)
    print('contributing area computed; max = %d' % fi['ncontrib'].max(), flush=True)

    src = get_sources_index_threshold(fi, 100)
    print('sources: %d  (C++ 112757)' % src.size, flush=True)

    print('building junction network ...', flush=True)
    fi = junction_network(fi, src)
    print('junctions: %d ; max stream order: %d' % (fi['JunctionVector'].size, fi['SO'].max()), flush=True)

    print('find_valleys ...', flush=True)
    vj = find_valleys(fi, tcurv, src)
    nvj = int(np.sum(vj != ND))
    print('valley junctions found: %d  (C++ _VJ = 768)' % nvj, flush=True)

    # validate vs _VJ raster (junction index at valley-outlet cells)
    vj_ref = np.fromfile(f'{ALG}/bailey_run_dem_VJ.flt', dtype='<f4').reshape(NR, NC)
    ref_mask = vj_ref != NDF
    mine_grid = np.full((NR, NC), ND, dtype=np.int64)
    nz = vj != ND
    mine_grid[fi['row_of'][nz], fi['col_of'][nz]] = vj[nz]
    mine_mask = mine_grid != ND
    print('_VJ cells: ref=%d mine=%d  overlap=%d  mine-only=%d ref-only=%d'
          % (ref_mask.sum(), mine_mask.sum(), int((ref_mask & mine_mask).sum()),
             int((mine_mask & ~ref_mask).sum()), int((ref_mask & ~mine_mask).sum())), flush=True)
    np.save('/tmp/dreich_vj_mine.npy', mine_grid)

    if stage == 'heads':
        import time
        t = time.time(); build_svector(fi); print('SVector %.1fs' % (time.time() - t), flush=True)
        t = time.time(); distance_from_outlet(fi); print('flow distance %.1fs' % (time.time() - t), flush=True)
        t = time.time()
        heads_temp, final = get_channel_heads_from_valleys(fi, z, vj, 10, 1000.0, 0.525)
        print('head finder %.1fs ; per-junction heads=%d ; after dedup=%d  (C++ _CH = 634)'
              % (time.time() - t, len(heads_temp), len(final)), flush=True)
        # validate vs _CH raster (channel-head node cells)
        ch_ref = np.fromfile(f'{ALG}/bailey_run_dem_CH.flt', dtype='<f4').reshape(NR, NC)
        ref_ch = ch_ref != NDF
        mine_ch = np.zeros((NR, NC), bool)
        fa = np.array(final, dtype=np.int64)
        mine_ch[fi['row_of'][fa], fi['col_of'][fa]] = True
        print('_CH cells: ref=%d mine=%d  overlap=%d  mine-only=%d ref-only=%d'
              % (ref_ch.sum(), mine_ch.sum(), int((ref_ch & mine_ch).sum()),
                 int((mine_ch & ~ref_ch).sum()), int((ref_ch & ~mine_ch).sum())), flush=True)
        np.save('/tmp/dreich_ch_mine.npy', mine_ch)

        # separate dedup vs per-junction-extraction differences
        ro_r, ro_c = np.where(ref_ch & ~mine_ch)
        ro_nodes = fi['NodeIndex'][ro_r, ro_c]
        temp_set = set(int(x) for x in heads_temp)
        in_temp = [int(n) in temp_set for n in ro_nodes]
        print('ref-only heads that ARE in my pre-dedup heads_temp (=> dedup diff): %d / %d'
              % (sum(in_temp), len(in_temp)), flush=True)
        np.save('/tmp/dreich_heads_temp.npy', heads_temp)
        np.save('/tmp/dreich_final.npy', np.array(final, dtype=np.int64))
