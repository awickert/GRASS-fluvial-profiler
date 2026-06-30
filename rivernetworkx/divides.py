"""Drainage-divide networks (Scherler & Schwanghart, 2020).

  Scherler, D. & Schwanghart, W. (2020). Drainage divide networks - Part 1:
  Identification and ordering in digital elevation models. Earth Surf. Dynam. 8,
  245-259. https://doi.org/10.5194/esurf-8-245-2020

A faithful, vectorised reimplementation: identify divides as the basin boundaries
(on pixel EDGES, infinitesimal width) of a D8 stream network, build the divide
network as a graph of segments between endpoints (at streams) and junctions,
resolve D8 diagonal-flow crossings (their ``divnet`` FX rule), open each basin
loop at its outlet (their ``getdivide``; pass per-basin outlet cells to
``extract_divide_edges(stream=)``), and order it by peeling segments from the
endpoints inward (Strahler / Shreve / Topo) -- identical to their ``sort.m``.

Built on the DrEICH FlowInfo (``rivernetworkx.dreich``): D8 receivers and
per-link drainage basins.

Status (June 2026): VERIFIED against the actual TopoToolbox ``DIVIDEobj``, run
under Octave on identical routing -- reproduces its divide network to ~98 % of
interior edges with Strahler order robustly matched (cross-check + reasoning in
``dev/ttb_crosscheck/``). Not yet wired into the package exports. Still to do:
per-divide relief metrics (HR / FD / DAI); then use the divides for channel
heads. The Topo/Shreve order tails are intentionally not matched (segmentation-
sensitive metrics; Strahler is the robust one).

Conventions
-----------
Cells are indexed ``(r, c)`` with ``r`` = 0 at the north edge. Divides live on the
``(nr+1) x (nc+1)`` grid of pixel CORNERS; a corner ``(i, j)`` is the NW corner of
cell ``(i, j)``. A divide EDGE is a unit segment between two adjacent corners,
lying on the shared boundary of two cells with different drainage-basin labels.
"""
import numpy as np


def extract_divide_edges(lab, nodata=-1, stream=None):
    """Divide edges = pixel boundaries between cells of different basin label.

    Parameters
    ----------
    lab : (nr, nc) int ndarray
        Per-cell drainage-basin label (e.g. the draining-link id from
        ``dreich.drainage_divides``); ``nodata`` marks cells off the network.
    nodata : int
        Label value to treat as "no basin" (still divides against real basins,
        as a domain/off-map boundary).
    stream : (nr, nc) bool ndarray, optional
        Cells at which to open the divide loops -- boundaries adjacent to a True
        cell are dropped, breaking each basin's boundary loop so the network is a
        tree that orders cleanly. For the **faithful** Scherler & Schwanghart
        ``getdivide`` (open each basin at its single outlet, keeping the rest of
        the perimeter), pass the per-basin **outlet** cells (each channel link's
        downstream-most cell). Cross-checked against TopoToolbox: outlet cells
        reproduce its divide network to ~98% (interior); a full *channel* mask
        also opens the loops but over-removes the near-stream perimeter that
        ``getdivide`` keeps (~92% interior). Divides then terminate at streams and
        do not run along them ("divide segments do not cross any rivers but their
        nodes may coincide with stream edges").

    Returns
    -------
    edges : (m, 2) int ndarray
        Each row is a divide edge as two corner ids ``(a, b)``, with a corner
        ``(i, j)`` encoded as ``i * (nc + 1) + j``. Vertical and horizontal edges
        are both included; the array is unordered.
    meta : dict
        ``nr, nc`` and the corner-grid width ``cw = nc + 1`` for decoding ids.
    """
    lab = np.asarray(lab)
    nr, nc = lab.shape
    cw = nc + 1

    def cid(i, j):
        return i * cw + j

    edges = []
    # vertical divide edges: boundary between cell (r,c) and (r,c+1), on the grid
    # line at column index c+1, from corner (r, c+1) to (r+1, c+1).
    vkeep = lab[:, :-1] != lab[:, 1:]
    if stream is not None:
        vkeep &= ~(stream[:, :-1] | stream[:, 1:])
    vr, vc = np.where(vkeep)                 # cell rows/cols of the LEFT cell
    a = cid(vr, vc + 1); b = cid(vr + 1, vc + 1)
    edges.append(np.column_stack([a, b]))

    # horizontal divide edges: boundary between cell (r,c) and (r+1,c), on the
    # grid line at row index r+1, from corner (r+1, c) to (r+1, c+1).
    hkeep = lab[:-1, :] != lab[1:, :]
    if stream is not None:
        hkeep &= ~(stream[:-1, :] | stream[1:, :])
    hr, hc = np.where(hkeep)                 # cell rows/cols of the UPPER cell
    a = cid(hr + 1, hc); b = cid(hr + 1, hc + 1)
    edges.append(np.column_stack([a, b]))

    edges = np.vstack(edges) if edges else np.empty((0, 2), dtype=np.int64)
    return edges.astype(np.int64), dict(nr=nr, nc=nc, cw=cw)


def decode_corner(cid, cw):
    """Corner id -> (i, j). With diagonal-crossing splitting, pass ``cid // 2``."""
    return np.asarray(cid) // cw, np.asarray(cid) % cw


def diagonal_flow_corners(recv, row_of, col_of, nr, nc):
    """Mark divide corners that sit at a D8 diagonal-flow crossing (Scherler &
    Schwanghart ``divnet`` ``FX`` grid). A cell flowing diagonally to its receiver
    passes through the shared corner; that corner is flagged 1 (NW-SE flow) or
    2 (NE-SW flow). Returns an int8 array indexed by corner id ``i*(nc+1)+j``."""
    r = np.asarray(row_of); c = np.asarray(col_of); recv = np.asarray(recv)
    dr = row_of[recv] - r
    dc = col_of[recv] - c
    diag = (np.abs(dr) == 1) & (np.abs(dc) == 1)
    ci = r + (dr > 0).astype(np.int64)            # shared corner of the diagonal step
    cj = c + (dc > 0).astype(np.int64)
    cw = nc + 1
    cid = ci * cw + cj
    fxval = np.where(dr * dc > 0, 1, 2).astype(np.int8)   # NW-SE vs NE-SW
    fx = np.zeros((nr + 1) * cw, dtype=np.int8)
    fx[cid[diag]] = fxval[diag]
    return fx


def split_diagonal_crossings(edges, fx, cw):
    """Resolve diagonal-flow crossings so the divide passes through rather than
    forming a spurious junction (S&S ``divnet``: a junction needs ``FX==0``).
    Returns edges in a DOUBLED corner-id space -- the real corner is ``new // 2``
    -- where a crossing corner becomes two independent through-nodes. The flow
    line separates the two through-paths: FX=1 (NW-SE flow) pairs the (N,E) edges
    and the (S,W) edges; FX=2 (NE-SW flow) pairs (N,W) and (E,S)."""
    pair0 = {1: ('N', 'E'), 2: ('N', 'W')}        # which directions take half 0

    def direction(frm, to):
        fi_, fj = divmod(frm, cw); ti_, tj = divmod(to, cw)
        if ti_ < fi_:
            return 'N'
        if ti_ > fi_:
            return 'S'
        return 'E' if tj > fj else 'W'

    def vid(corner, other):
        f = fx[corner]
        if f == 0:
            return corner * 2
        return corner * 2 + (0 if direction(corner, other) in pair0[f] else 1)
    new = np.empty_like(edges)
    for k in range(len(edges)):
        a = int(edges[k, 0]); b = int(edges[k, 1])
        new[k, 0] = vid(a, b); new[k, 1] = vid(b, a)
    return new


def build_divide_graph(edges):
    """Corner-node graph of the divide network: classify nodes by divide degree
    and trace maximal segments between breaks (endpoints/junctions).

    Endpoints (degree 1) are topological leaves (domain boundary). Junctions
    (degree >= 3) are where divide segments meet -- in the basin dual a stream
    confluence is naturally a degree-3 junction (two tributary basins + the trunk
    basin), so the network is NOT cut at streams here; stream-coincidence is used
    later as the *ordering seed* (see :func:`order_divides`), not as a cut.
    Degree-2 corners are interior to a segment. Pure cycles (an enclosed basin
    with no break) are traced as closed segments.

    Returns a dict with: ``adj`` (node -> neighbour list), ``degree`` (node ->
    int), ``endpoints``, ``junctions`` (sorted node lists), and ``segments`` (list
    of corner-id paths; break-to-break, or a cycle with the first node repeated at
    the end)."""
    adj = {}
    for a, b in edges:
        a = int(a); b = int(b)
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    degree = {n: len(ns) for n, ns in adj.items()}
    endpoints = sorted(n for n, d in degree.items() if d == 1)
    junctions = sorted(n for n, d in degree.items() if d >= 3)
    breaks = set(endpoints) | set(junctions)

    def canon(u, v):
        return (u, v) if u < v else (v, u)

    used = set()
    segments = []
    # 1. segments anchored at a break node, walking through degree-2 nodes
    for s in breaks:
        for nb in adj[s]:
            if canon(s, nb) in used:
                continue
            path = [s, nb]
            used.add(canon(s, nb))
            prev, cur = s, nb
            while cur not in breaks:
                nxts = [x for x in adj[cur] if x != prev]
                if not nxts:                       # dead end (shouldn't happen)
                    break
                nxt = nxts[0]
                if canon(cur, nxt) in used:
                    break
                used.add(canon(cur, nxt))
                path.append(nxt)
                prev, cur = cur, nxt
            segments.append(path)
    # 2. leftover edges form pure cycles (no break node)
    for a, b in edges:
        a = int(a); b = int(b)
        if canon(a, b) in used:
            continue
        path = [a, b]
        used.add(canon(a, b))
        prev, cur = a, b
        while cur != a:
            nxts = [x for x in adj[cur] if x != prev and canon(cur, x) not in used]
            if not nxts:
                break
            nxt = nxts[0]
            used.add(canon(cur, nxt))
            path.append(nxt)
            prev, cur = cur, nxt
        segments.append(path)
    return dict(adj=adj, degree=degree, endpoints=endpoints,
                junctions=junctions, segments=segments)


def _combine(orders, scheme):
    """Combine tributary orders at a junction (S&S 2020, Eq. set). An empty list
    (a leaf at an endpoint) seeds the base order 1."""
    if not orders:
        return 1
    m = max(orders)
    if scheme == 'topo':                       # +1 at every junction
        return m + 1
    if scheme == 'strahler':                   # +1 only if the max is shared
        return m + 1 if sum(o == m for o in orders) >= 2 else m
    if scheme == 'shreve':                      # additive magnitude
        return sum(orders)
    raise ValueError('scheme must be topo/strahler/shreve')


def order_divides(graph, leaf_segments=None, scheme='topo'):
    """Order the divide network by growing order from the stream-terminating
    segments inward (Scherler & Schwanghart, 2020). The **leaf segments** (those
    that reach a river -- where divides terminate) are order-1; order then grows
    away from the streams toward the main interfluves, combined at each junction
    via :func:`_combine`. ``leaf_segments`` is a per-segment boolean mask; it
    defaults to segments incident to a topological degree-1 endpoint (correct for
    small tree networks / tests). Also returns divide distance ``dd`` (max
    corner-steps from a leaf to the segment's inner end). True cycle segments (no
    break termini) stay order 0.

    Returns ``(order, dd)``: arrays indexed like ``graph['segments']``."""
    from collections import deque
    segs = graph['segments']
    nseg = len(segs)
    order = np.zeros(nseg, dtype=np.int64)
    dd = np.zeros(nseg, dtype=np.float64)
    breaks = set(graph['endpoints']) | set(graph['junctions'])
    incident = {}
    termini = []
    for i, p in enumerate(segs):
        a, b = p[0], p[-1]
        if a == b or a not in breaks or b not in breaks:
            termini.append(None)               # true cycle -> stays order 0
            continue
        termini.append((a, b))
        incident.setdefault(a, []).append(i)
        incident.setdefault(b, []).append(i)
    node_deg = {n: len(v) for n, v in incident.items()}
    sorted_seg = np.zeros(nseg, dtype=bool)
    incoming = {n: [] for n in incident}       # tributary orders arrived
    incoming_d = {n: [0.0] for n in incident}  # tributary distances arrived

    if leaf_segments is None:                   # default: touch a degree-1 endpoint
        ep = set(graph['endpoints'])
        leaf_segments = np.array([t is not None and (t[0] in ep or t[1] in ep)
                                  for t in termini])
    else:
        leaf_segments = np.asarray(leaf_segments, dtype=bool)

    q = deque()

    def push(o, d, n):
        node_deg[n] -= 1
        incoming[n].append(o)
        incoming_d[n].append(d)
        if node_deg[n] == 1:
            q.append(n)

    # seed: each leaf segment is order-1 and feeds BOTH its end junctions
    for i, t in enumerate(termini):
        if t is None or not leaf_segments[i]:
            continue
        order[i] = 1
        dd[i] = len(segs[i]) - 1
        sorted_seg[i] = True
        push(1, dd[i], t[0])
        push(1, dd[i], t[1])
    # peel inward: a junction with one unordered segment left releases it
    while q:
        n = q.popleft()
        if node_deg.get(n, 0) != 1:
            continue
        e = next((i for i in incident[n] if not sorted_seg[i]), None)
        if e is None:
            continue
        order[e] = _combine(incoming[n], scheme)
        dd[e] = max(incoming_d[n]) + (len(segs[e]) - 1)
        sorted_seg[e] = True
        a, b = termini[e]
        far = b if a == n else a
        node_deg[n] -= 1
        push(order[e], dd[e], far)
    return order, dd
