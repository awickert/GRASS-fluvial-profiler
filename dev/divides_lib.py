"""Shared helpers for the divides->channel-heads prototype (Mid Bailey Run, 1 m).

Loads the cached FlowInfo from divides_setup.py and provides the pieces both
Design A (divide-shape valley gate) and Design B (chi-z self-gating) build on:

  * load_cache()              -- unpickle fi, curvature, filled, reference heads
  * valley_flank_grid()       -- bool grid: cells flanked by divides on both
                                 sides (the divide analog of convergent curvature)
  * find_valleys_gated()      -- find_valleys walk param'd by a per-node bool gate
  * head_and_score()          -- chi-z head finder that ALSO returns the fit score
  * dedup_furthest_upstream() -- the DrEICH furthest-upstream head dedup
  * compare_heads()           -- recall/precision of test heads vs a reference

Nothing here imports GRASS; it works straight off the pickle.
"""
import pickle
import numpy as np
from scipy import ndimage
from rivernetworkx import dreich as D
from rivernetworkx.dreich import (_find_upstream_junction, _find_farthest_upslope,
                                   _build_channel_chi, _reg32)

CACHE = '/tmp/divides_proto_fi.pkl'
ND = -9999
# Mid Bailey Run 1 m geo-reference (raw DEM corner; res = 1 m)
WEST, NORTH = 398767.32685636, 4369656.1523925
CLUBB_XLSX = '/tmp/clubb_channel_heads.xlsx'


def load_cache():
    with open(CACHE, 'rb') as f:
        return pickle.load(f)


# --------------------------------------------------------- divide-flank signal
def _dir_within(divide, W, axis, forward):
    """Bool grid: a divide cell lies within W steps along (axis, forward).

    Uses a running 'index of last divide' via np.maximum.accumulate, so it is
    O(N) and fully vectorised."""
    nr, nc = divide.shape
    n = divide.shape[axis]
    idx = np.arange(n)
    shape = (1, n) if axis == 1 else (n, 1)
    coord = idx.reshape(shape)
    pos = np.where(divide, np.broadcast_to(coord, divide.shape), -10 * n)
    if forward:                       # divide at HIGHER index, within +W
        rev = np.flip(pos, axis=axis)
        nearest = np.flip(np.maximum.accumulate(rev, axis=axis), axis=axis)
        return (nearest - coord <= W) & (nearest >= 0)
    else:                             # divide at LOWER index, within -W
        nearest = np.maximum.accumulate(pos, axis=axis)
        return (coord - nearest <= W) & (nearest >= 0)


def valley_flank_grid(divide, W, diagonals=False):
    """Cells flanked by divide cells on BOTH sides within W, along H or V (and
    optionally the two diagonals). The divide-topology analog of a convergent
    (valley-bottom) cross-section: a valley bottom has ridges on both flanks; a
    planar hillslope or spur does not."""
    left = _dir_within(divide, W, axis=1, forward=False)
    right = _dir_within(divide, W, axis=1, forward=True)
    up = _dir_within(divide, W, axis=0, forward=False)
    down = _dir_within(divide, W, axis=0, forward=True)
    flanked = (left & right) | (up & down)
    if diagonals:
        # rotate 45deg via shear is messy; approximate diagonals by checking the
        # max-filtered divide in diagonal-offset windows using ndimage on a
        # rotated copy is overkill -- H/V suffices for a first cut (see notes).
        pass
    # a divide cell itself is a ridge, not a valley bottom
    return flanked & ~divide


# ---------------------------------------------- find_valleys with a generic gate
def find_valleys_gated(fi, gate, sources, n_connecting_nodes=10):
    """LSDJunctionNetwork::find_valleys with the per-cell curvature test replaced
    by a per-NODE boolean ``gate`` (True = valley-like cell). Returns the same
    node-indexed valley array find_valleys does, so channel_heads_from_valleys
    consumes it unchanged."""
    recv = fi['recv']; SO = fi['SO']; N = fi['N']
    visited = np.zeros(N, dtype=bool)
    valley = np.full(N, ND, dtype=np.int64)
    for s in sources:
        end = False; consec = 0; cur = int(s)
        while not end:
            rcv = recv[cur]
            visited[cur] = True
            consec = consec + 1 if gate[cur] else 0
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
                cur = int(rcv)
    return valley


def gate_from_grid(fi, grid):
    """Per-node bool gate from a (nr,nc) bool grid."""
    return grid[fi['row_of'], fi['col_of']]


# --------------------------------------------- chi-z head finder WITH the score
def head_and_score(fi, dem, dn, A_0=1000.0, m_over_n=0.525, min_segment_length=10):
    """Build the hilltop->dn chi-z profile and return (head_node, best_test,
    n_profile). best_test is the LSDChannel chi-z split score (R2_chan minus the
    hillslope autocorrelation term) -- the self-gating signal for Design B."""
    hilltop = _find_farthest_upslope(fi, int(dn))
    nodeseq, chi, elev = _build_channel_chi(fi, dem, hilltop, int(dn), A_0, m_over_n)
    n = len(chi)
    if n < 2 * min_segment_length:
        return int(dn), -np.inf, n
    f32 = np.float32
    best = f32(0.0); elev_inter = None
    for hill in range(min_segment_length, n - min_segment_length + 1):
        r2_chan, _ = _reg32(chi[hill:], elev[hill:])
        _, dw_hill = _reg32(chi[:hill], elev[:hill])
        test = f32(r2_chan - f32((dw_hill - f32(2.0)) / f32(2.0)))
        if test > best:
            best = test; elev_inter = elev[hill]
    if elev_inter is None:
        return int(nodeseq[0]), float(best), n
    head = int(nodeseq[0])
    for i in range(n):
        if elev[i] == elev_inter:
            head = int(nodeseq[i])
    return head, float(best), n


# --------------------------------------------------- furthest-upstream head dedup
def dedup_furthest_upstream(fi, dem, head_nodes):
    """DrEICH GetChannelHeadsChiMethodFromValleys dedup: of heads sharing a flow
    path, keep only the furthest upstream (highest-elevation) one."""
    recv = fi['recv']; row_of = fi['row_of']; col_of = fi['col_of']
    heads = np.asarray(list(head_nodes), dtype=np.int64)
    if heads.size == 0:
        return []
    elevs = dem[row_of[heads], col_of[heads]].astype(np.float64)
    sorted_heads = heads[np.argsort(-elevs, kind='stable')]
    binary = {int(h): 1 for h in heads}
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
    return [h for h in binary if binary[h] == 1]


# --------------------------------------------------------------- comparison metric
def load_clubb():
    """Clubb et al. (2014) Mid Bailey Run field heads as (E, N) projected coords."""
    import openpyxl
    wb = openpyxl.load_workbook(CLUBB_XLSX, data_only=True, read_only=True)
    ws = wb['Sheet1']
    E, N = [], []
    for row in ws.iter_rows(min_row=3, values_only=True):
        if row[0] and 'Bailey' in str(row[0]):
            N.append(float(row[2])); E.append(float(row[3]))
    return np.column_stack([E, N])


def nodes_xy(fi, nodes, res=1.0):
    """Map node ids -> (x, y) projected coords (cell centres)."""
    nodes = np.asarray(nodes, dtype=np.int64)
    return np.column_stack([WEST + (fi['col_of'][nodes] + 0.5) * res,
                            NORTH - (fi['row_of'][nodes] + 0.5) * res])


def recall_vs_points(test_xy, ref_xy, tols=(10, 30, 50, 100)):
    """Recall (frac of ref points with a test point within tol) + median nearest
    distance. For sparse ground-truth (Clubb) precision is not meaningful."""
    from scipy.spatial import cKDTree
    if len(test_xy) == 0:
        return {t: 0.0 for t in tols}, np.inf
    d, _ = cKDTree(test_xy).query(ref_xy)
    return {t: float(np.mean(d <= t)) for t in tols}, float(np.median(d))


def heads_rc(fi, nodes):
    return np.column_stack([fi['row_of'][np.asarray(nodes, dtype=np.int64)],
                            fi['col_of'][np.asarray(nodes, dtype=np.int64)]]).astype(float)


def compare_heads(test_rc, ref_rc, res=1.0, tols=(5.0, 10.0, 30.0)):
    """Recall (frac of ref heads with a test head within tol) and precision
    (frac of test heads within tol of a ref head). Distances in map units."""
    from scipy.spatial import cKDTree
    out = {}
    if len(test_rc) == 0 or len(ref_rc) == 0:
        for t in tols:
            out[t] = (0.0, 0.0)
        return out
    tt = cKDTree(test_rc * res); rt = cKDTree(ref_rc * res)
    d_ref, _ = tt.query(ref_rc * res)      # nearest test head to each ref head
    d_test, _ = rt.query(test_rc * res)    # nearest ref head to each test head
    for t in tols:
        recall = float(np.mean(d_ref <= t))
        prec = float(np.mean(d_test <= t))
        out[t] = (recall, prec)
    return out
