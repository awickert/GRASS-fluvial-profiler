############################################################################
#
# MODULE:       rivernetworkx.core
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Pure-NetworkX river-network construction, traversal, and I/O.
#               No GRASS dependencies: operates on plain edge records and
#               NetworkX graphs, so it is unit-testable anywhere and reusable
#               outside GRASS (e.g. the GRLP coupling).
#
# COPYRIGHT:    (c) 2025-2026 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v3). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################

import json
import math

import numpy as np
import networkx as nx
from networkx.readwrite import json_graph

# Node id of the off-map outlet, following GRASS GIS convention (tostream = 0).
OFFMAP = 0


#############
# GEOMETRY  #
#############

def segment_distances(x, y):
    """
    Along-segment distance from the upstream end (``s_downstream``) and from the
    downstream end (``s_upstream``), given vertex coordinates ordered
    upstream -> downstream.

    Returns (s_downstream, s_upstream), both 1-D arrays the length of x/y.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    along = np.hstack([0.0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    return along, along[-1] - along


#####################
# GRAPH CONSTRUCTION #
#####################

def build_graph(edge_records, outlet=OFFMAP):
    """
    Build a river-network DiGraph from edge records.

    Each record is a dict describing one stream segment:
        cat       (int)    segment id            -> source node
        tostream  (int)    downstream-neighbor   -> target node (``outlet`` if off-map)
        x, y      (arrays) vertex coords, ordered upstream -> downstream
        <other>   (arrays, optional) any further per-vertex quantities, e.g.
                  z (elevation), A (flow accumulation), slope -- carried through
                  generically.

    In the returned DiGraph, nodes are tributary junctions (plus the outlet),
    per-vertex arrays live on edges, values shared at junctions live on nodes,
    and ``node['s']`` is cumulative distance upstream of the outlet (computed by
    a branching-safe breadth-first sweep).
    """
    edge_records = list(edge_records)

    # Per-vertex attributes are any record keys beyond the ids and raw geometry.
    reserved = {'cat', 'tostream', 'x', 'y'}
    extra = []
    for rec in edge_records:
        for key, value in rec.items():
            if key not in reserved and value is not None and key not in extra:
                extra.append(key)
    attr_names = ['x', 'y', 's_upstream', 's_downstream'] + extra

    G = nx.DiGraph()
    for rec in edge_records:
        x = np.asarray(rec['x'], dtype=float)
        y = np.asarray(rec['y'], dtype=float)
        s_down, s_up = segment_distances(x, y)
        attrs = {'x': x, 'y': y, 's_upstream': s_up, 's_downstream': s_down}
        for key in extra:
            if rec.get(key) is not None:
                attrs[key] = np.asarray(rec[key], dtype=float)
        cat = int(rec['cat'])
        G.add_edge(cat, int(rec['tostream']), cat=cat, **attrs)

    _pull_first_to_parents(G, attr_names)
    _drop_downstream_endpoints(G, attr_names, outlet)
    _init_outlet(G, attr_names, outlet)
    _accumulate_distance(G, outlet)
    return G


def _pull_first_to_parents(G, attr_names):
    """
    Move the upstream-most entry of each edge array onto its upstream (parent)
    node, so that adjacent segments share values at tributary junctions.
    """
    for parent, child, data in G.edges(data=True):
        for attr in attr_names:
            arr = data.get(attr)
            if arr is None or len(arr) == 0:
                continue
            first = arr[0]
            node_attr = G.nodes[parent].get(attr, [])
            if not isinstance(node_attr, list):
                node_attr = [node_attr]
            node_attr.append(first)
            G.nodes[parent][attr] = node_attr
            if isinstance(arr, list):
                arr.pop(0)
            else:
                data[attr] = arr[1:]


def _drop_downstream_endpoints(G, attr_names, outlet):
    """
    Remove the downstream-most entry of each edge array, which now duplicates the
    value held on the downstream node. Edges that leave the map (child == outlet)
    keep their downstream point, since those elevations are real and may differ.
    """
    for parent, child, data in G.edges(data=True):
        if child == outlet:
            continue
        for attr in attr_names:
            arr = data.get(attr)
            if arr is None or len(arr) == 0:
                continue
            if isinstance(arr, list):
                arr.pop(-1)
            else:
                data[attr] = arr[:-1]


def _init_outlet(G, attr_names, outlet):
    """
    The outlet node is beyond the domain: give it defined placeholder values so
    the cumulative-distance sweep has a base case and NaNs do not propagate.
    """
    node = G.nodes[outlet]
    node['s_upstream'] = [0]
    node['s'] = [0]            # total distance upstream of the outlet starts at 0
    for attr in attr_names:
        if attr != 's_upstream':
            node.setdefault(attr, [np.nan])


def _accumulate_distance(G, outlet):
    """
    Cumulative distance upstream of the outlet, accumulated across the whole
    network by a breadth-first sweep from the outlet. Branching-safe.
    """
    for n in bfs_upward(G, outlet):
        for parent, child in G.in_edges(n):
            child_s = G.nodes[child]['s'][0]
            G.nodes[parent]['s'] = [child_s + G.nodes[parent]['s_upstream'][0]]
            G.edges[parent, child]['s'] = child_s + G.edges[parent, child]['s_upstream']


#############
# TRAVERSAL #
#############

def bfs_upward(G, start):
    """
    Yield nodes upstream of ``start`` in breadth-first order. Use this to sweep
    the network from an outlet, relying on downstream nodes being visited before
    the nodes above them.
    """
    R = G.reverse(copy=False)  # a view, no data copied
    for node in nx.bfs_tree(R, start):
        yield node


def upstream_subnetwork(G, node):
    """Return the sub-network of everything that drains to ``node`` (inclusive)."""
    return G.subgraph(nx.ancestors(G, node) | {node}).copy()


def downstream_path(G, node, outlet=OFFMAP):
    """Return the single downstream path of node ids from ``node`` to ``outlet``."""
    return nx.shortest_path(G, node, outlet)


###############################
# PROFILE ASSEMBLY / RESAMPLE #
###############################

def assemble_downstream_profile(records_by_cat, path, attrs=()):
    """
    Concatenate per-segment arrays into one continuous downstream profile.

    records_by_cat : {cat: record dict} with per-vertex ``x``/``y`` (and any
                     attrs), vertices ordered upstream -> downstream.
    path           : downstream-ordered cats [most upstream, ..., last on map];
                     e.g. ``downstream_path(G, start, outlet)`` with the outlet
                     dropped.
    attrs          : extra per-vertex attribute names to carry (e.g. ``'z'``,
                     ``'A'``, ``'slope'``).

    Returns a dict with ``s`` (cumulative distance upstream of the network
    mouth), ``cat`` (the source segment of each point), ``x``, ``y``, and each
    requested attribute, each a 1-D array running upstream -> downstream.
    Vertices shared at junctions are de-duplicated.
    """
    path = list(path)
    # Per-segment along-distance, and distance from the mouth at each segment's
    # downstream end (accumulated walking upstream from the mouth).
    s_down = {cat: segment_distances(records_by_cat[cat]['x'],
                                     records_by_cat[cat]['y'])[0]
              for cat in path}
    offset = {}
    running = 0.0
    for cat in reversed(path):                       # downstream -> upstream
        offset[cat] = running
        running += s_down[cat][-1]

    names = ['x', 'y'] + list(attrs)
    out = {name: [] for name in names}
    out['s'] = []
    out['cat'] = []                                  # source segment of each point
    prev = None
    for cat in path:                                 # upstream -> downstream
        rec = records_by_cat[cat]
        # distance from the mouth, decreasing downstream within the segment
        s_seg = offset[cat] + (s_down[cat][-1] - s_down[cat])
        start = 1 if prev is not None else 0         # drop shared junction point
        kept = s_seg[start:]
        out['s'].append(kept)
        out['cat'].append(np.full(len(kept), cat, dtype=int))
        for name in names:
            out[name].append(np.asarray(rec[name], dtype=float)[start:])
        prev = cat
    for name in list(out):
        out[name] = (np.concatenate(out[name]) if out[name]
                     else np.array([], dtype=float))
    return out


def densify(s, arrays, dx_target):
    """
    Resample 1-D arrays to ~uniform spacing ``dx_target`` along ``s``.

    s        : monotonic 1-D distance (either direction).
    arrays   : {name: 1-D array} sampled at ``s``.
    Returns (new_s ascending, {name: resampled array}); endpoints preserved,
    linear interpolation.
    """
    s = np.asarray(s, dtype=float)
    order = np.argsort(s)
    s_sorted = s[order]
    # Collapse coincident stations so np.interp has strictly increasing x.
    uniq = np.concatenate(([True], np.diff(s_sorted) > 0))
    s_u = s_sorted[uniq]
    if len(s_u) < 2:                       # zero-length / single-station segment
        return s_sorted, {name: np.asarray(arr, dtype=float)[order]
                          for name, arr in arrays.items()}
    s0, s1 = s_u[0], s_u[-1]
    n = max(int(np.ceil((s1 - s0) / float(dx_target))) + 1, 2)
    new_s = np.linspace(s0, s1, n)
    out = {name: np.interp(new_s, s_u,
                           np.asarray(arr, dtype=float)[order][uniq])
           for name, arr in arrays.items()}
    return new_s, out


def moving_average(s, y, window):
    """
    Centered moving average of ``y`` over a distance ``window`` in ``s`` units.

    ``s`` need not be evenly spaced; NaNs are ignored. Returns an array the same
    length as ``y``.
    """
    s = np.asarray(s, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n == 0:
        return np.empty(0, dtype=float)
    # Sort by s, then use prefix sums + searchsorted to get each point's window
    # sum and (non-NaN) count in O(1): O(n log n) overall, not the naive O(n^2)
    # per-point mask. NaNs are excluded from each window's mean.
    order = np.argsort(s, kind='stable')
    ss = s[order]
    yy = y[order]
    finite = ~np.isnan(yy)
    csum = np.concatenate(([0.0], np.cumsum(np.where(finite, yy, 0.0))))
    ccnt = np.concatenate(([0], np.cumsum(finite.astype(int))))
    half = window / 2.0
    lo = np.searchsorted(ss, ss - half, side='left')
    hi = np.searchsorted(ss, ss + half, side='right')
    cnt = ccnt[hi] - ccnt[lo]
    tot = csum[hi] - csum[lo]
    res = np.where(cnt > 0, tot / np.where(cnt > 0, cnt, 1), np.nan)
    out = np.empty(n, dtype=float)
    out[order] = res
    return out


def smooth_segment(record, attrs, window):
    """
    Per-vertex moving-average smoothing of each name in ``attrs`` within a
    single segment ``record``, over the segment's own along-distance.

    Segments are the reaches between tributary junctions, so smoothing here
    never crosses a junction: natural slope/area breaks at confluences are
    preserved, and a segment's smoothed values are independent of which
    flow path it belongs to. Returns ``{attr: smoothed array}``.
    """
    s_down, _ = segment_distances(record['x'], record['y'])
    return {a: moving_average(s_down, np.asarray(record[a], dtype=float), window)
            for a in attrs}


def channel_slope(s, z):
    """
    Downstream channel slope ``S = -dz/ds`` along ordered segment vertices.

    ``s`` is along-segment distance increasing downstream; ``z`` is elevation.
    S is positive where the channel descends downstream. Computed with a
    non-uniform-spacing gradient. NaN for segments shorter than two points.
    """
    s = np.asarray(s, dtype=float)
    z = np.asarray(z, dtype=float)
    if len(z) < 2:
        return np.full(len(z), np.nan)
    # Coincident vertices (zero spacing) make the gradient non-finite; report
    # those as NaN rather than +/-inf so they drop out of slope-area analysis.
    with np.errstate(divide='ignore', invalid='ignore'):
        S = -np.gradient(z, s)
    S[~np.isfinite(S)] = np.nan
    return S


def slope_area(records, window=None, log=False):
    """
    Paired channel slope and drainage area across the whole network, for
    slope--area analysis (e.g. locating the fluvial -> hillslope rollover).

    For each segment independently (junctions never crossed): optionally smooth
    elevation ``z`` and area ``A`` over ``window``, compute channel slope
    ``S = -dz/ds`` from the (smoothed) ``z``, and pair it with the (smoothed)
    ``A``. Each record must carry ``z`` and ``A``.

    Returns ``(A, S)`` 1-D arrays concatenated over all segments. With
    ``log=True`` returns ``(log10(A), log10(S))`` keeping only ``A > 0`` and
    ``S > 0`` (so the pairs are ready for a log--log fit / rollover search).
    """
    A_all, S_all = [], []
    for rec in records:
        s_down, _ = segment_distances(rec['x'], rec['y'])
        if window is not None:
            smoothed = smooth_segment(rec, ('z', 'A'), window)
            z, A = smoothed['z'], smoothed['A']
        else:
            z = np.asarray(rec['z'], dtype=float)
            A = np.asarray(rec['A'], dtype=float)
        A_all.append(A)
        S_all.append(channel_slope(s_down, z))
    A = np.concatenate(A_all) if A_all else np.array([], dtype=float)
    S = np.concatenate(S_all) if S_all else np.array([], dtype=float)
    if log:
        keep = (A > 0) & (S > 0)
        return np.log10(A[keep]), np.log10(S[keep])
    return A, S


def chi(area, flow_distance, theta=0.5, ref_area=1.0, base_chi=0.0):
    """
    Integral channel coordinate ``chi`` (Perron & Royden, 2013) along a single
    flow path:

        chi(x) = integral of (ref_area / A(x'))^theta  dx'

    integrated *upstream* from a downstream reference value ``base_chi``. ``chi``
    has units of length and linearises the steady-state detachment-limited long
    profile (elevation becomes ~linear in ``chi`` along a channel), so plotting
    elevation against ``chi`` separates channel reaches (linear) from hillslopes
    (non-linear) -- the basis of the DrEICH channel-head method (Clubb et al.,
    2014) -- and exposes profile disequilibrium.

    Parameters
    ----------
    area : (n,) array_like
        Drainage area at each node along one flow path (same units as
        ``ref_area``); must be strictly positive.
    flow_distance : (n,) array_like
        Cumulative flow distance along the same path (monotonic; only successive
        differences are used, so any consistent origin is fine).
    theta : float
        Concavity index ``m/n`` (Clubb et al. use ~0.45-0.5).
    ref_area : float
        Reference drainage area ``A_0``; sets the (arbitrary) scale of ``chi``.
    base_chi : float
        ``chi`` at the downstream end of the path.

    Returns
    -------
    chi : (n,) ndarray
        ``chi`` at each node, in the same order as the inputs.

    Notes
    -----
    Order-agnostic: the downstream end is taken to be the higher-drainage-area
    end, and ``chi`` increases upstream. Discretised as in LSDTopoTools
    (Clubb et al., 2014): each step adds ``(ref_area / A_upstream)^theta``
    times the step length.
    """
    area = np.asarray(area, dtype=float)
    flow_distance = np.asarray(flow_distance, dtype=float)
    n = len(area)
    if n == 0:
        return np.array([], dtype=float)
    if n == 1:
        return np.array([base_chi], dtype=float)
    # Orient upstream (index 0) -> downstream (index -1) by drainage area, so the
    # downstream (largest-area) node anchors at base_chi and chi grows upstream.
    flip = area[0] > area[-1]
    if flip:
        area = area[::-1]
        flow_distance = flow_distance[::-1]
    seg = np.abs(np.diff(flow_distance))            # length of each up->down step
    integrand = (ref_area / area) ** theta          # evaluated at the upstream node
    incr = integrand[:-1] * seg                     # chi added across each step
    # chi[i] = base_chi + sum(incr[i:]); downstream node (i = n-1) = base_chi
    chi_vals = np.concatenate([np.cumsum(incr[::-1])[::-1], [0.0]]) + base_chi
    return chi_vals[::-1] if flip else chi_vals


def _linfit_r2_dw(x, y):
    """Least-squares line ``y ~ x``; return ``(R2, durbin_watson)`` of the fit,
    matching LSDTopoTools ``simple_linear_regression``: ``R2 = 1 - SSE/SST`` and
    ``DW = sum((e_i - e_{i-1})^2) / sum(e_i^2)`` with residual = predicted -
    observed (DW and R2 are sign-invariant in the residual)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xm = x.mean()
    sxx = np.sum((x - xm) ** 2)
    if sxx == 0:
        return 0.0, 2.0
    m = np.sum((x - xm) * (y - y.mean())) / sxx
    resid = (m * x + (y.mean() - m * xm)) - y
    sst = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / sst if sst > 0 else 0.0
    bottom = np.sum(resid ** 2) or 1e-10
    dw = np.sum(np.diff(resid) ** 2) / bottom
    return r2, dw


def channel_head_chi_split(chi, elevation, min_segment_length=10):
    """
    DrEICH channel-head detector on a single hilltop -> junction profile
    (Clubb et al., 2014).

    Given ``chi`` and ``elevation`` ordered from the upslope hilltop (index 0,
    largest chi) down to a downstream junction (smallest chi), try every split
    into an upslope **hillslope** segment and a downslope **channel** segment,
    each at least ``min_segment_length`` nodes. For each split, regress
    elevation on ``chi`` within each segment and score

        test = R2_channel - (DurbinWatson_hillslope - 2) / 2

    which rewards a channel reach that is *linear in chi-z* (high R2) and a
    hillslope reach that is *not* (autocorrelated residuals -> DW far below 2).
    The channel head is the first node of the best-scoring channel segment.

    Returns the index (into the input arrays) of the channel head, or ``-1`` if
    the profile is shorter than ``2 * min_segment_length`` nodes.
    """
    chi = np.asarray(chi, dtype=float)
    elevation = np.asarray(elevation, dtype=float)
    n = len(chi)
    if n < 2 * min_segment_length:
        return -1
    best_test = -np.inf
    best_idx = -1
    for hill_len in range(min_segment_length, n - min_segment_length + 1):
        r2_chan, _ = _linfit_r2_dw(chi[hill_len:], elevation[hill_len:])
        _, dw_hill = _linfit_r2_dw(chi[:hill_len], elevation[:hill_len])
        test = r2_chan - (dw_hill - 2.0) / 2.0
        if test > best_test:
            best_test = test
            best_idx = hill_len
    return best_idx


def fit_sa_break(logA, logS, knots=None, min_side=10):
    """
    Locate the hillslope -> fluvial break in a slope--area cloud by a continuous,
    constrained broken-stick fit in log--log space:

        log S = h                              for log A <= log A*   (hillslope)
        log S = h - theta * (log A - log A*)   for log A >  log A*   (fluvial)

    The hillslope limb is flat (slope independent of area); the fluvial limb is a
    power law ``S ~ A^-theta``. The two limbs join at the knot ``log A*`` -- the
    colluvial-to-fluvial transition drainage area (where channels become
    fluvial). Fitting both limbs
    jointly -- scan the knot, least-squares ``(h, slope)`` at each, keep the knot
    with the smallest residual -- is the rigorous form of "fit the channel power
    law, fit the hillslope constant, take their intersection": the limbs meet at
    A* by construction.

    Operates on whatever ``(logA, logS)`` it is given -- the whole network
    (one basin-wide A*) or a single flowline (a local A*). Feed it the
    ``log=True`` output of :func:`slope_area`, ideally after dropping flat-valley
    artifacts (S ~ 0) and any unwanted high-area trunk reaches.

    Parameters
    ----------
    logA, logS : array-like
        Paired log10 drainage area and log10 channel slope.
    knots : array-like, optional
        Candidate ``log A*`` values to scan. Default: 200 points spanning the
        5th--95th percentile of ``logA`` (knots at the extremes leave too little
        data on one side to constrain a limb).
    min_side : int
        Require at least this many points on each side of the knot, so neither
        limb is fit from a handful of points.

    Returns
    -------
    dict or None
        ``{'logA_star', 'A_star', 'theta', 'hillslope_logS', 'rss', 'n'}``, or
        ``None`` if there are too few finite points to fit.
    """
    logA = np.asarray(logA, dtype=float)
    logS = np.asarray(logS, dtype=float)
    good = np.isfinite(logA) & np.isfinite(logS)
    logA, logS = logA[good], logS[good]
    if logA.size < 2 * min_side:
        return None
    if knots is None:
        lo, hi = np.percentile(logA, [5, 95])
        knots = np.linspace(lo, hi, 200)

    best = None
    for k in np.asarray(knots, dtype=float):
        n_left = int(np.count_nonzero(logA <= k))
        if n_left < min_side or (logA.size - n_left) < min_side:
            continue
        # Basis: a constant plus a hinge that is 0 on the hillslope side and
        # rises downstream of the knot, so the fitted left limb is exactly flat.
        X = np.column_stack([np.ones_like(logA), np.maximum(0.0, logA - k)])
        coef, *_ = np.linalg.lstsq(X, logS, rcond=None)
        rss = float(np.sum((logS - X @ coef) ** 2))
        if best is None or rss < best['rss']:
            h, slope = float(coef[0]), float(coef[1])
            best = {'logA_star': float(k), 'A_star': float(10.0 ** k),
                    'theta': -slope, 'hillslope_logS': h,
                    'rss': rss, 'n': int(logA.size)}
    return best


def colluvial_fluvial_transition(records, A_star):
    """
    Locate the colluvial-to-fluvial transition: the points where drainage area
    first reaches ``A_star`` going downstream -- i.e. where each channel becomes
    fluvial (the downslope limit of the colluvial hollow), with ``A_star`` the
    slope--area break from :func:`fit_sa_break`.

    Flow accumulation already integrates the upstream contributing area, so no
    network topology is needed: within a segment whose upstream end is below
    ``A_star`` and downstream end is at or above it, the crossing is exactly that
    channel's transition. The location is linearly interpolated in area between
    the two bracketing vertices.

    Each record needs per-vertex ``A`` (drainage area) and ``x``/``y``. Vertices
    may be ordered either way; the segment is oriented to ascending area first.

    Returns a list of ``(x, y, cat)`` -- one transition per channel that crosses
    into the fluvial regime within the network. Channels entirely below
    ``A_star`` (colluvial throughout the extracted network) or entirely above it
    (downstream of their transition) yield none.

    This is a PROCESS boundary, not a morphological field channel head -- the
    latter lies upslope, within the colluvial hollow (see r.fluvial.channelheads).
    """
    heads = []
    for rec in records:
        A = np.asarray(rec['A'], dtype=float)
        x = np.asarray(rec['x'], dtype=float)
        y = np.asarray(rec['y'], dtype=float)
        if A.size < 2:
            continue
        if A[0] > A[-1]:                       # orient upstream -> downstream
            A, x, y = A[::-1], x[::-1], y[::-1]
        if not (A[0] < A_star <= A[-1]):
            continue
        i = int(np.argmax(A >= A_star))        # first vertex at/above A_star
        if i == 0:
            xh, yh = x[0], y[0]
        else:
            f = (A_star - A[i - 1]) / (A[i] - A[i - 1])
            xh = x[i - 1] + f * (x[i] - x[i - 1])
            yh = y[i - 1] + f * (y[i] - y[i - 1])
        heads.append((float(xh), float(yh), int(rec['cat'])))
    return heads


#######
# I/O #
#######

def _to_jsonable(obj):
    """Recursively convert numpy arrays/scalars to plain Python for JSON.

    Non-finite floats (NaN/inf) become None, so the output is standards-compliant
    JSON (bare ``NaN`` is not valid JSON and strict parsers reject it);
    ``load_json`` restores them to NaN.
    """
    if isinstance(obj, np.ndarray):
        return _to_jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return _to_jsonable(obj.item())
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, (list, tuple)):
        return type(obj)(_to_jsonable(x) for x in obj)
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    return obj


def make_json_safe(G):
    """Return a deep copy of G with all numpy types converted to plain Python."""
    H = G.copy()
    for k in list(H.graph.keys()):
        H.graph[k] = _to_jsonable(H.graph[k])
    for _, data in H.nodes(data=True):
        for k in list(data.keys()):
            data[k] = _to_jsonable(data[k])
    for _, _, data in H.edges(data=True):
        for k in list(data.keys()):
            data[k] = _to_jsonable(data[k])
    return H


def to_json_dict(G):
    """Return the node-link dict for G (JSON-serializable)."""
    return json_graph.node_link_data(make_json_safe(G))


def export_json(G, path):
    """Write G to ``path`` as NetworkX node-link JSON."""
    with open(path, "w") as f:
        json.dump(to_json_dict(G), f, indent=2)


def _none_to_nan(obj):
    """Restore JSON null (written by _to_jsonable for non-finite floats) to NaN."""
    if isinstance(obj, list):
        return [_none_to_nan(x) for x in obj]
    if obj is None:
        return float('nan')
    return obj


def load_json(path):
    """Read a NetworkX node-link JSON file back into a DiGraph.

    JSON null inside attribute arrays (written for NaN by export_json) is
    restored to NaN.
    """
    with open(path) as f:
        G = json_graph.node_link_graph(json.load(f))
    for _, data in G.nodes(data=True):
        for k in list(data):
            data[k] = _none_to_nan(data[k])
    for _, _, data in G.edges(data=True):
        for k in list(data):
            data[k] = _none_to_nan(data[k])
    return G
