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
    make_json_safe,
    to_json_dict,
    export_json,
    load_json,
)

# grass_io requires a GRASS session; keep the package importable without one.
try:
    from .grass_io import read_stream_vector, build_network  # noqa: F401
except Exception:
    pass
