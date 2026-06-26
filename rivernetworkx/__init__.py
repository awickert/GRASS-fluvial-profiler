"""
rivernetworkx — river-network construction, traversal, and I/O on NetworkX.

The single in-memory representation shared across the GRASS fluvial toolkit:
``core`` is pure NetworkX (no GRASS, unit-testable, reusable by e.g. the GRLP
coupling); ``grass_io`` reads a GRASS stream-network vector into that
representation. GRASS modules are separate processes, so each rebuilds the graph
from the persisted vector at runtime.
"""

from .core import (
    OFFMAP,
    segment_distances,
    build_graph,
    bfs_upward,
    upstream_subnetwork,
    downstream_path,
    assemble_downstream_profile,
    densify,
    moving_average,
    smooth_segment,
    channel_slope,
    slope_area,
    fit_sa_break,
    channel_head_points,
    make_json_safe,
    to_json_dict,
    export_json,
    load_json,
)

# grass_io imports GRASS lazily (inside its reader functions), so this works
# without a GRASS session; build_network/read_stream_vector only need GRASS when
# actually called. sample_raster and assemble_records are pure.
from .grass_io import (  # noqa: F401
    sample_raster,
    offmap_inflow_cats,
    assemble_records,
    read_stream_vector,
    build_network,
)
