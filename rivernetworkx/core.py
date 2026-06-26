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
        z, A      (arrays, optional) elevation / flow accumulation per vertex

    In the returned DiGraph, nodes are tributary junctions (plus the outlet),
    per-vertex arrays live on edges, values shared at junctions live on nodes,
    and ``node['s']`` is cumulative distance upstream of the outlet (computed by
    a branching-safe breadth-first sweep).
    """
    G = nx.DiGraph()
    attr_names = ['x', 'y', 's_upstream', 's_downstream']
    have_z = have_A = False

    for rec in edge_records:
        x = np.asarray(rec['x'], dtype=float)
        y = np.asarray(rec['y'], dtype=float)
        s_down, s_up = segment_distances(x, y)
        attrs = {'x': x, 'y': y, 's_upstream': s_up, 's_downstream': s_down}
        if rec.get('z') is not None:
            attrs['z'] = np.asarray(rec['z'], dtype=float)
            have_z = True
        if rec.get('A') is not None:
            attrs['A'] = np.asarray(rec['A'], dtype=float)
            have_A = True
        cat = int(rec['cat'])
        G.add_edge(cat, int(rec['tostream']), cat=cat, **attrs)

    if have_z:
        attr_names.append('z')
    if have_A:
        attr_names.append('A')

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
