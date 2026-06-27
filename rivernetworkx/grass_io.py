############################################################################
#
# MODULE:       rivernetworkx.grass_io
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Read a GRASS stream-network vector (as linked by v.stream.network)
#               and the rasters sampled along it into rivernetworkx edge records,
#               then into a NetworkX graph.
#
# COPYRIGHT:    (c) 2025-2026 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v3). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# GRASS is imported lazily, inside the functions that need a live session, so
# this module (and its pure helpers below) import without GRASS -- which keeps
# the sampling and assembly logic unit-testable.

from contextlib import contextmanager

import numpy as np

from .core import OFFMAP, build_graph


########################
# PURE (no GRASS)      #
########################

def sample_raster(array, x, y, *, west, north, nsres, ewres):
    """
    Nearest-neighbor sample of a north-up raster ``array`` (row 0 = northmost)
    at map coordinates ``x`` (easting) and ``y`` (northing). Vectorized over
    arrays of points. Points outside the raster return NaN.

    west/north are the outer edges of the region; nsres/ewres the cell sizes.
    """
    array = np.asarray(array, dtype=float)
    nrows, ncols = array.shape
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    east = west + ncols * ewres
    south = north - nrows * nsres
    col = np.floor((x - west) / ewres).astype(int)
    row = np.floor((north - y) / nsres).astype(int)
    # Inside-ness is judged in map coordinates, inclusive of the outer edges; a
    # point exactly on the east/south edge belongs to the last cell, where
    # np.floor lands one index past the end, so clip rather than drop it.
    inside = (x >= west) & (x <= east) & (y >= south) & (y <= north)
    col = np.clip(col, 0, ncols - 1)
    row = np.clip(row, 0, nrows - 1)
    out = np.full(np.shape(x), np.nan, dtype=float)
    out[inside] = array[row[inside], col[inside]]
    return out


def offmap_inflow_cats(records):
    """
    Segment cats whose sampled accumulation is negative at the UPSTREAM-most
    vertex -- i.e. off-map contributing area reaches the channel head, which
    marks an incomplete catchment.

    GRASS ``r.watershed`` writes NEGATIVE flow accumulation for cells that may
    receive flow from outside the region, and the flag propagates downstream.
    In a COMPLETE catchment this is confined to the outlet's boundary cells (the
    downstream end), where the channel legitimately leaves the map -- that is
    normal and tolerated. An INCOMPLETE catchment instead has the negative reach
    up to a segment head, so we flag a segment only when its upstream-most
    accumulation is negative.

    Currently a hard error at build time; making off-map inflow queryable is
    enhancement issue #9. (Assumes off-map inflow enters at a segment head,
    which r.stream.extract segmentation gives; diffuse lateral mid-reach inflow,
    and a divide placed exactly on the region edge, are not handled -- issue #9.)
    """
    bad = []
    for rec in records:
        A = rec.get('A')
        if A is None:
            continue
        A = np.asarray(A, dtype=float)
        if A.size and A[0] < 0:          # negative at the channel head
            bad.append(rec['cat'])
    return bad


def assemble_records(cats, tostream, geometry, **vertex_attrs):
    """
    Build rivernetworkx edge records from already-read GRASS data.

    cats:       iterable of segment ids
    tostream:   {cat: downstream-neighbor cat}
    geometry:   {cat: (x_array, y_array)} vertices ordered upstream -> downstream
    vertex_attrs: name -> {cat: array}, e.g. z={...}, A={...}, slope={...};
                  a value of None is skipped.
    """
    records = []
    for cat in cats:
        x, y = geometry[cat]
        rec = {'cat': int(cat), 'tostream': int(tostream[cat]), 'x': x, 'y': y}
        for name, per_cat in vertex_attrs.items():
            if per_cat is not None:
                rec[name] = per_cat[cat]
        records.append(rec)
    return records


########################
# GRASS readers (lazy) #
########################

def _read_topology_table(streams):
    """Return (cats, tostream-dict) from the vector's attribute table."""
    import pandas as pd
    from grass.script import vector_db_select
    # !!!! GENERALIZE COLUMN NAME: 'tostream' is hardcoded here (and in
    # v.stream.network's tostream_cat_column option). Plumb the chosen column
    # name through instead of assuming the default. !!!!
    sel = vector_db_select(streams)
    df = pd.DataFrame(data=sel['values'].values(), columns=sel['columns'])
    if 'tostream' not in df.columns:
        from grass.script import fatal
        fatal("Column 'tostream' not found in <%s>. Run v.stream.network first."
              % streams)
    cats = [int(c) for c in df['cat']]
    tostream = {int(c): int(t) for c, t in zip(df['cat'], df['tostream'])}
    return cats, tostream


def _read_geometry(streams, cats):
    """Return {cat: (x_array, y_array)} of segment vertices.

    Consecutive coincident vertices are dropped: zero-length steps give a
    repeated along-distance, which makes np.interp (densify) ill-defined and
    np.gradient (channel_slope) divide by zero downstream.
    """
    from grass.pygrass.vector import VectorTopo
    vt = VectorTopo(streams)
    vt.open('r')
    geometry = {}
    for cat in cats:
        coords = vt.cat(cat_id=cat, vtype='lines')[0]
        en = coords.to_array()
        x, y = en[:, 0], en[:, 1]
        keep = np.concatenate(([True], (np.diff(x) != 0) | (np.diff(y) != 0)))
        geometry[cat] = (x[keep], y[keep])
    vt.close()
    return geometry


@contextmanager
def _region_over(geometry, align):
    """
    Temporarily limit the GRASS region to the network's extent, aligned to the
    ``align`` raster and grown one cell, so raster reads cover the whole network
    and no more (efficiency), and so a region that does not already span the
    network does not silently sample NaN. Restores the prior region on exit.

    Uses an explicit save/restore (not use_temp_region) so it nests safely
    inside a caller that already set its own temporary region (e.g. gunittest).
    """
    from grass.script import run_command, region as _region
    xs = np.concatenate([gx for gx, _ in geometry.values()])
    ys = np.concatenate([gy for _, gy in geometry.values()])
    saved = _region()
    pad_ns, pad_ew = float(saved['nsres']), float(saved['ewres'])
    try:
        # Pad the network bbox by a cell on each side: this both gives a
        # one-cell margin (so edge vertices sample cleanly) and guarantees a
        # non-degenerate extent for a perfectly straight, single-row or
        # single-column network -- a bare n==s / e==w would make g.region
        # reject the call. align= snaps the padded bounds back to the grid.
        run_command('g.region',
                    n=float(ys.max()) + pad_ns, s=float(ys.min()) - pad_ns,
                    e=float(xs.max()) + pad_ew, w=float(xs.min()) - pad_ew,
                    align=align, quiet=True)
        yield
    finally:
        run_command('g.region', n=saved['n'], s=saved['s'], e=saved['e'],
                    w=saved['w'], nsres=saved['nsres'], ewres=saved['ewres'],
                    quiet=True)


def _read_raster(rastname):
    """Return (north-up array, region-dict) for a raster over the current region.

    NULL cells are loaded as NaN (not garray's default fill, which would silently
    pass missing data through as 0); this keeps "no data" honest downstream.
    """
    from grass.script import array as garray
    from grass.pygrass.gis.region import Region
    arr = np.asarray(garray.array(mapname=rastname, null=np.nan))
    reg = Region()
    bounds = {'west': reg.west, 'north': reg.north,
              'nsres': reg.nsres, 'ewres': reg.ewres}
    return arr, bounds


def _read_geometry_all(streams):
    """Return {cat: (x_array, y_array)} for every line in the vector.

    Like _read_geometry but discovers the cats itself (one pass over the lines)
    instead of being handed a cat list, so callers that do not have / do not want
    the topology table (e.g. the channel-head detector reading a raw
    r.stream.extract network) need not run v.stream.network first. Consecutive
    coincident vertices are dropped (see _read_geometry).
    """
    from grass.pygrass.vector import VectorTopo
    vt = VectorTopo(streams)
    vt.open('r')
    geometry = {}
    for ln in vt.viter('lines'):
        if ln.cat is None:
            continue
        en = ln.to_array()
        x, y = en[:, 0], en[:, 1]
        keep = np.concatenate(([True], (np.diff(x) != 0) | (np.diff(y) != 0)))
        if keep.sum() >= 2:
            geometry[int(ln.cat)] = (x[keep], y[keep])
    vt.close()
    return geometry


def _read_records(geometry, tostream, elevation=None, accumulation=None,
                  slope=None, accum_mult=1.0, assume_complete=False):
    """
    Sample rasters along the given geometry and assemble validated edge records.

    Shared body of read_stream_vector (cats+tostream from the attribute table)
    and read_stream_segments (geometry only, no topology). ``geometry`` is
    {cat: (x, y)}; ``tostream`` is {cat: downstream cat} (all OFFMAP for the
    topology-free reader). Applies the NaN-coverage, off-map-inflow, and
    abs(accumulation) handling. Requires a live GRASS session.
    """
    cats = list(geometry)
    specs = (('z', elevation, 1.0),
             ('A', accumulation, accum_mult),
             ('slope', slope, 1.0))
    samples = {name: None for name, _, _ in specs}
    align = elevation or accumulation or slope
    if align:
        # Sample within a region clipped to the network (see _region_over).
        with _region_over(geometry, align):
            for name, rastname, mult in specs:
                if rastname:
                    arr, bounds = _read_raster(rastname)
                    samples[name] = {cat: mult * sample_raster(arr, gx, gy,
                                                               **bounds)
                                     for cat, (gx, gy) in geometry.items()}

    records = assemble_records(cats, tostream, geometry,
                               z=samples['z'], A=samples['A'],
                               slope=samples['slope'])

    from grass.script import fatal

    # Accumulation must be defined on every channel cell. A NaN means the
    # accumulation raster does not cover the network (NULL under the channel);
    # it would otherwise slip past the off-map head check below (NaN < 0 is
    # False) and be masked by abs(), so fail loudly.
    nan_cats = [rec['cat'] for rec in records
                if rec.get('A') is not None
                and np.isnan(np.asarray(rec['A'], dtype=float)).any()]
    if nan_cats:
        fatal("Accumulation raster has no data at channel cell(s) on "
              "segment(s) %s: it must cover the full stream network."
              % ', '.join(str(c) for c in nan_cats))

    # An incomplete catchment (off-map upstream contributing area) shows up as
    # negative flow accumulation reaching a channel head. Negative accumulation
    # only at the outlet boundary (the downstream end) is the normal map exit
    # and is fine. Off-map inflow is not yet supported: fail loudly rather than
    # build a misleading network -- unless the caller asserts a complete basin
    # (assume_complete), in which case a negative head is a boundary artifact
    # (a divide on the region edge). (issue #9)
    if not assume_complete:
        bad = offmap_inflow_cats(records)
        if bad:
            fatal("Negative flow accumulation reaches the head of segment(s) "
                  "%s: the catchment is incomplete (off-map contributing area "
                  "upstream; r.watershed marks such cells negative). Negative "
                  "accumulation only at the outlet boundary is fine. If the "
                  "region truly contains the full basin (a divide on the edge), "
                  "re-run asserting a complete catchment; otherwise use a "
                  "region/DEM that fully contains it, or omit accumulation. "
                  "(enhancement: issue #9)"
                  % ', '.join(str(c) for c in bad))

    # Elevation/slope with no data under a channel is less critical than
    # accumulation -- it stays NaN rather than masking a real value -- and slope
    # is legitimately NULL on region-boundary cells (r.slope.area needs
    # neighbours), so warn rather than error, but don't let the gap be silent.
    from grass.script import warning
    for name, label in (('z', 'elevation'), ('slope', 'slope')):
        gap = [rec['cat'] for rec in records
               if rec.get(name) is not None
               and np.isnan(np.asarray(rec[name], dtype=float)).any()]
        if gap:
            warning("%s raster has no data at channel cell(s) on segment(s) %s; "
                    "those points are NaN."
                    % (label, ', '.join(str(c) for c in gap)))

    # With off-map inflow ruled out, any remaining negative accumulation is
    # r.watershed's conservative boundary flag at the outlet, whose MAGNITUDE is
    # the true drainage area. Use the magnitude so the outlet (largest-area)
    # cells survive A>0 filtering in slope-area analysis, plots, and the export.
    for rec in records:
        if rec.get('A') is not None:
            rec['A'] = np.abs(np.asarray(rec['A'], dtype=float))
    return records


def read_stream_vector(streams, elevation=None, accumulation=None, slope=None,
                       accum_mult=1.0, assume_complete=False):
    """
    Read a v.stream.network-linked vector (and any rasters) into edge records.
    Requires a live GRASS session.

    ``assume_complete``: when True, the caller asserts the region contains the
    full basin, so negative accumulation at a channel head is treated as an
    r.watershed boundary artifact (a divide on the region edge) rather than
    off-map inflow, and the off-map error is skipped. The NaN-coverage check is
    always enforced.
    """
    cats, tostream = _read_topology_table(streams)
    geometry = _read_geometry(streams, cats)
    return _read_records(geometry, tostream, elevation=elevation,
                         accumulation=accumulation, slope=slope,
                         accum_mult=accum_mult, assume_complete=assume_complete)


def read_stream_segments(streams, elevation=None, accumulation=None, slope=None,
                         accum_mult=1.0, assume_complete=False):
    """
    Read a stream vector into edge records WITHOUT requiring topology.

    Like read_stream_vector but reads geometry straight from the lines and does
    not need a ``tostream`` column, so a raw r.stream.extract network can be used
    directly -- no v.stream.network step (and no O(N^2) linking) first. Each
    record carries ``tostream = OFFMAP`` (topology is absent); this is what the
    slope--area / fluvial-channel-head analysis needs, which is per-segment and
    uses flow accumulation rather than the network graph.
    """
    geometry = _read_geometry_all(streams)
    tostream = {cat: OFFMAP for cat in geometry}
    return _read_records(geometry, tostream, elevation=elevation,
                         accumulation=accumulation, slope=slope,
                         accum_mult=accum_mult, assume_complete=assume_complete)


def build_network(streams, elevation=None, accumulation=None, slope=None,
                  accum_mult=1.0, outlet=OFFMAP, assume_complete=False):
    """Read a GRASS stream-network vector and build the NetworkX DiGraph.

    ``assume_complete`` is forwarded to read_stream_vector (skip the off-map
    inflow error when the caller asserts the full basin is in the region).
    """
    records = read_stream_vector(streams, elevation=elevation,
                                 accumulation=accumulation, slope=slope,
                                 accum_mult=accum_mult,
                                 assume_complete=assume_complete)
    return build_graph(records, outlet=outlet)
