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
    strahler_order,
    assemble_downstream_profile,
    densify,
    moving_average,
    smooth_segment,
    channel_slope,
    slope_area,
    chi,
    channel_head_chi_split,
    fit_sa_break,
    colluvial_fluvial_transition,
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
    read_stream_segments,
    build_network,
)

# DrEICH channel-head extraction (Clubb et al., 2014): a faithful numpy port of the
# LSDTopoTools pipeline. Pure (no GRASS); scipy is imported lazily inside the
# curvature step. Used by r.fluvial.channelheads method=dreich.
from .dreich import (  # noqa: F401
    extract_channel_heads,
    channel_network_segments,
    build_flowinfo_from_directions,
    directions_from_flowinfo,
    fill as fill_dem,
    tangential_curvature,
)
