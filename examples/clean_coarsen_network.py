#!/usr/bin/env python3
############################################################################
#
# SCRIPT:       clean_coarsen_network.py
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Post-process a river network exported by v.stream.network
#               (json= option): despike and smooth the elevations along each
#               segment, then coarsen (resample) the network, and re-export
#               it as JSON.
#
#               This is the "stage 2" step between v.stream.network (which
#               extracts and links the network from a DEM) and a downstream model such
#               as GRLP (which evolves its long profiles). Raw DEM-sampled
#               elevations are noisy and densely sampled; this cleans and
#               thins them while preserving the network topology.
#
# COPYRIGHT:    (c) 2025 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v3). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# Developed for the AGU 2025 presentation:
#   Wickert, A. D., and F. McNab (2025), Simulating Geomorphic Evolution
#   Through River Networks, EP23D-1702, AGU Fall Meeting, New Orleans, LA, USA.
#
# NOTE: bfs_upward(), _to_jsonable(), and make_json_safe_graph() are also
# provided by the importable `rivernetworkx` package (as bfs_upward,
# export_json/load_json, make_json_safe). They are kept vendored here so this
# example stays self-contained and runnable without installing rivernetworkx;
# import them from rivernetworkx instead if you prefer a single source.

import argparse
import copy
import json

import numpy as np
from scipy.interpolate import interp1d
import networkx as nx
from networkx.readwrite import json_graph


#############
# UTILITIES #
#############

def bfs_upward(G, start):
    """
    Breadth-first search going upwards.
    Use this to update overall distance upstream from all outlets
    in a single sweep through the network, relying on those network components
    closer to the outlet to be updated first
    """
    R = G.reverse(copy=False)  # just a view, no data duplication
    for node in nx.bfs_tree(R, start):
        yield node


def set_outlet_elevations(G, outlet_drop=4e-4):
    """
    v.stream.network exports outlet (off-map) nodes with z = nan, because
    they lie beyond the domain. If left as nan, the smoothing step propagates
    the nan upstream through the network. Replace each nan outlet elevation
    with the lowest channel elevation entering that node, dropped by
    `outlet_drop` so that a downstream gradient is preserved (e.g., a sand-bed
    slope times the cell size).
    """
    for node in G.nodes():
        z = G.nodes[node].get('z')
        if z is None:
            continue
        # Outlet nodes have no downstream segment and an undefined elevation
        if G.out_degree(node) == 0 and np.any(np.isnan(np.asarray(z, dtype=float))):
            in_edges = list(G.in_edges(node))
            if not in_edges:
                continue
            # Elevation at the outlet end (last array entry) of each inflow
            z_inflow = [G.edges[u, v]['z'][-1] for u, v in in_edges]
            G.nodes[node]['z'] = [min(z_inflow) - outlet_drop]


###############################################
# SMOOTH AND DESPIKE VALUES ON NETWORKX EDGES #
###############################################

def despike_1d(y, window=2, k=3.5):
    """
    Remove local outliers from a 1D array using a sliding median + MAD.

    y: 1D numpy array
    window: radius for local window (actual window size = 2*window + 1)
    k: threshold in units of local MAD
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy()

    x = y.copy()
    is_outlier = np.zeros(n, dtype=bool)

    for i in range(1, n - 1):  # typically keep endpoints as-is here
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        neighborhood = y[lo:hi]

        median = np.median(neighborhood)
        mad = np.median(np.abs(neighborhood - median))
        if mad == 0:
            continue

        z = 0.6745 * (y[i] - median) / mad
        if np.abs(z) > k:
            is_outlier[i] = True

    # Replace outliers with local median of non-outlier neighbors
    for i in np.where(is_outlier)[0]:
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        neighborhood = x[lo:hi]
        good = ~is_outlier[lo:hi]
        if np.any(good):
            x[i] = np.median(neighborhood[good])
        else:
            x[i] = np.median(neighborhood)  # fallback

    return x


def smooth_1d_fixed_ends(y, left_val, right_val, alpha=0.35, n_iters=15):
    """
    Diffusive smoothing of a 1D signal with fixed endpoints (Dirichlet BC).

    y: 1D numpy array (original data along the edge)
    left_val, right_val: fixed boundary values from the nodes
    alpha: smoothing strength per iteration (0 < alpha < 1)
    n_iters: number of iterations
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n == 0:
        return y.copy()
    if n == 1:
        # Degenerate edge; just pick something between the two nodes
        return np.array([(left_val + right_val) / 2.0])

    x = y.copy()
    x[0] = left_val
    x[-1] = right_val

    for _ in range(n_iters):
        x_new = x.copy()
        # update only interior points
        for i in range(1, n - 1):
            neighbor_mean = 0.5 * (x[i - 1] + x[i + 1])
            x_new[i] = (1 - alpha) * x[i] + alpha * neighbor_mean
        # re-enforce fixed endpoints each iteration
        x_new[0] = left_val
        x_new[-1] = right_val
        x = x_new

    return x


def clean_edge_values(
    G,
    node_attr="value",
    edge_attr="values",
    despike_window=2,
    despike_k=3.5,
    alpha=0.35,
    n_iters=15
):
    """
    For each edge in G:
      - take its values array
      - despike locally
      - smooth with fixed endpoints given by node_attr of endpoints
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        edge_vals = np.asarray(data[edge_attr], dtype=float)
        if edge_vals.size == 0:
            continue

        # node values (single-item lists)
        left_val = float(G.nodes[u][node_attr][0])
        right_val = float(G.nodes[v][node_attr][0])

        # 1) despike along the edge
        de_noised = despike_1d(edge_vals, window=despike_window, k=despike_k)

        # 2) smooth with fixed endpoint values from nodes
        smoothed = smooth_1d_fixed_ends(de_noised, left_val, right_val,
                                        alpha=alpha, n_iters=n_iters)

        # write back
        data[edge_attr] = smoothed.tolist()

    return G


###################
# COARSEN NETWORK #
###################

def downsample_edge_along_s_new(
    G,
    keep_fraction=0.1,
    attrs=("x", "y", "z", "A", "s"),
):
    """
    Return a NEW graph where each edge has had (x,y,z,A,s) resampled
    along s so that:

      - s is evenly spaced between min(s) and max(s)
      - endpoints are excluded (inset by one spacing each side)
      - ~ keep_fraction of original points remain
      - x, y, z, A are interpolated on s_new
      - Original graph G is left untouched
    """

    # Make a deep copy so we don't modify original
    H = copy.deepcopy(G)

    x_attr, y_attr, z_attr, A_attr, s_attr = attrs

    for u, v, data in H.edges(data=True):

        # Required geometry/elevation attrs; accumulation (A) is optional, since
        # v.stream.network json= may be run without an accumulation raster.
        if any(data.get(a) is None for a in (s_attr, x_attr, y_attr, z_attr)):
            continue
        s = np.asarray(data[s_attr])
        x = np.asarray(data[x_attr])
        y = np.asarray(data[y_attr])
        z = np.asarray(data[z_attr])
        has_A = data.get(A_attr) is not None
        A = np.asarray(data[A_attr]) if has_A else None

        # All present arrays must have the same length
        n = len(s)
        if not (len(x) == len(y) == len(z) == n and (not has_A or len(A) == n)):
            continue

        # Need >=3 points to create interior
        if n < 3:
            continue

        # Sort by s
        idx = np.argsort(s)
        s = s[idx]
        x = x[idx]
        y = y[idx]
        z = z[idx]
        if has_A:
            A = A[idx]

        # How many interior samples?
        n_new = int(round(keep_fraction * n))
        n_new = max(1, n_new)

        # Even spacing in s, excluding endpoints
        s_min, s_max = s[0], s[-1]
        ds = (s_max - s_min) / (n_new + 1)

        # Interior points only
        s_new = s_min + ds * np.arange(1, n_new + 1)
        s_new = s_new[::-1]  # Maintain original descending direction

        # Linear interpolators
        fx = interp1d(s, x, kind="linear")
        fy = interp1d(s, y, kind="linear")
        fz = interp1d(s, z, kind="linear")

        # Interpolate and store back on the NEW graph
        data[s_attr] = s_new.tolist()
        data[x_attr] = fx(s_new).tolist()
        data[y_attr] = fy(s_new).tolist()
        data[z_attr] = fz(s_new).tolist()
        if has_A:
            data[A_attr] = interp1d(s, A, kind="linear")(s_new).tolist()

        # Drop companion per-vertex arrays we did not resample (the export's
        # s_upstream/s_downstream, and any sampled slope): left untouched they
        # keep their original length and would disagree with x/y/z/A/s.
        for stale in ("s_upstream", "s_downstream", "slope"):
            data.pop(stale, None)

    return H


##########
# EXPORT #
##########

def _to_jsonable(obj):
    """
    Use this for export.
    Recursively convert objects into JSON-serializable forms:
    - numpy arrays -> lists
    - numpy scalars -> Python scalars
    - containers (dict/list/tuple) -> converted elementwise
    """

    # --- NumPy arrays ---
    if isinstance(obj, np.ndarray):
        # convert to list, then recurse (in case nested arrays, objects)
        return _to_jsonable(obj.tolist())

    # --- NumPy scalar types (float32, int64, etc.) ---
    if isinstance(obj, np.generic):
        return obj.item()  # returns a plain Python scalar

    # --- Containers: list/tuple ---
    if isinstance(obj, (list, tuple)):
        return type(obj)(_to_jsonable(x) for x in obj)

    # --- Dicts ---
    if isinstance(obj, dict):
        # Ensure values are converted; keys should already be JSON-OK
        return {k: _to_jsonable(v) for k, v in obj.items()}

    # Everything else is returned as-is; must already be JSON-serializable
    return obj


def make_json_safe_graph(G):
    """
    Use this for export
    Return a deep-copied version of G where:
    - all numpy arrays are converted to lists
    - all numpy scalar types are converted to Python scalars
    - node, edge, and graph attributes are processed
    """

    # deep copy so we don't mutate the original
    H = copy.deepcopy(G)

    # graph-level attributes
    for k in list(H.graph.keys()):
        H.graph[k] = _to_jsonable(H.graph[k])

    # node attributes
    for n, data in H.nodes(data=True):
        for k in list(data.keys()):
            data[k] = _to_jsonable(data[k])

    # edge attributes
    for u, v, data in H.edges(data=True):
        for k in list(data.keys()):
            data[k] = _to_jsonable(data[k])

    return H


def export_graph(G, outjson):
    """Write a NetworkX river-network graph to a JSON file."""
    H = make_json_safe_graph(G)
    data = json_graph.node_link_data(H)
    with open(outjson, "w") as f:
        json.dump(data, f, indent=2)


############
# PLOTTING #
############

def plot_long_profiles(G, outlet=0):
    """Plot every segment's long profile (elevation vs. distance upstream)."""
    from matplotlib import pyplot as plt
    plt.figure()
    for n in bfs_upward(G, outlet):
        for parent, child in G.in_edges(n):
            plt.plot(G.edges[parent, child]['s'], G.edges[parent, child]['z'],
                     'k-', linewidth=3, alpha=1)
            plt.plot(G.nodes[parent]['s'], G.nodes[parent]['z'], 'ko', alpha=1)
    plt.xlabel('Distance upstream of outlet')
    plt.ylabel('Elevation')
    plt.show()


########
# MAIN #
########

def main():
    parser = argparse.ArgumentParser(
        description="Despike, smooth, and coarsen a v.stream.network JSON "
                    "river network, then re-export it as JSON."
    )
    parser.add_argument("input", help="input JSON from v.stream.network json=")
    parser.add_argument("output", help="output (cleaned, coarsened) JSON")
    parser.add_argument("--outlet", type=int, default=0,
                        help="outlet node id (default: 0)")
    parser.add_argument("--outlet-drop", type=float, default=4e-4,
                        help="elevation drop applied at outlet nodes to keep a "
                             "downstream gradient (default: 4e-4)")
    parser.add_argument("--despike-window", type=int, default=100,
                        help="despiking half-window, in points (default: 100)")
    parser.add_argument("--despike-k", type=float, default=2.0,
                        help="despiking threshold in MADs (default: 2.0)")
    parser.add_argument("--alpha", type=float, default=0.85,
                        help="smoothing strength per iteration (default: 0.85)")
    parser.add_argument("--n-iters", type=int, default=2000,
                        help="number of smoothing iterations (default: 2000)")
    parser.add_argument("--keep-fraction", type=float, default=0.1,
                        help="fraction of points kept when coarsening "
                             "(default: 0.1)")
    parser.add_argument("--plot", action="store_true",
                        help="show long profiles of the cleaned network")
    args = parser.parse_args()

    # Load
    with open(args.input) as f:
        G = json_graph.node_link_graph(json.load(f))

    # Give outlet nodes real elevations so smoothing does not propagate nan
    set_outlet_elevations(G, outlet_drop=args.outlet_drop)

    # Despike + smooth elevations along each segment
    G = clean_edge_values(G, 'z', 'z',
                          despike_window=args.despike_window,
                          despike_k=args.despike_k,
                          alpha=args.alpha,
                          n_iters=args.n_iters)

    # Coarsen (resample) the network along s
    H = downsample_edge_along_s_new(G, keep_fraction=args.keep_fraction)

    if args.plot:
        plot_long_profiles(H, outlet=args.outlet)

    # Export
    export_graph(H, args.output)
    print("Wrote %s" % args.output)


if __name__ == "__main__":
    main()
