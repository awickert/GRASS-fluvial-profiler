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
            G.nodes[parent]['s'] = [G.nodes[child]['s'][0]
                                    + G.nodes[parent]['s_upstream'][0]]
            G.edges[parent, child]['s'] = (G.nodes[child]['s'][0]
                                           + G.edges[parent, child]['s_upstream'])


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
    mouth), ``x``, ``y``, and each requested attribute, each a 1-D array running
    upstream -> downstream. Vertices shared at junctions are de-duplicated.
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
    prev = None
    for cat in path:                                 # upstream -> downstream
        rec = records_by_cat[cat]
        # distance from the mouth, decreasing downstream within the segment
        s_seg = offset[cat] + (s_down[cat][-1] - s_down[cat])
        start = 1 if prev is not None else 0         # drop shared junction point
        out['s'].append(s_seg[start:])
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
    s0, s1 = s_sorted[0], s_sorted[-1]
    n = max(int(np.ceil((s1 - s0) / float(dx_target))) + 1, 2)
    new_s = np.linspace(s0, s1, n)
    out = {name: np.interp(new_s, s_sorted, np.asarray(arr, dtype=float)[order])
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
    half = window / 2.0
    out = np.empty(len(y), dtype=float)
    for i, si in enumerate(s):
        sel = (s >= si - half) & (s <= si + half)
        vals = y[sel]
        vals = vals[~np.isnan(vals)]
        out[i] = np.mean(vals) if len(vals) else np.nan
    return out


#######
# I/O #
#######

def _to_jsonable(obj):
    """Recursively convert numpy arrays/scalars to plain Python for JSON."""
    if isinstance(obj, np.ndarray):
        return _to_jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
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


def load_json(path):
    """Read a NetworkX node-link JSON file back into a DiGraph."""
    with open(path) as f:
        return json_graph.node_link_graph(json.load(f))
