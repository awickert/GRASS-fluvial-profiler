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

from .core import build_graph


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
    try:
        run_command('g.region', n=float(ys.max()), s=float(ys.min()),
                    e=float(xs.max()), w=float(xs.min()), align=align,
                    quiet=True)
        run_command('g.region', grow=1, quiet=True)
        yield
    finally:
        run_command('g.region', n=saved['n'], s=saved['s'], e=saved['e'],
                    w=saved['w'], nsres=saved['nsres'], ewres=saved['ewres'],
                    quiet=True)


def _read_raster(rastname):
    """Return (north-up array, region-dict) for a raster over the current region."""
    from grass.script import array as garray
    from grass.pygrass.gis.region import Region
    arr = np.asarray(garray.array(mapname=rastname))
    reg = Region()
    bounds = {'west': reg.west, 'north': reg.north,
              'nsres': reg.nsres, 'ewres': reg.ewres}
    return arr, bounds


def read_stream_vector(streams, elevation=None, accumulation=None, slope=None,
                       accum_mult=1.0):
    """
    Read a v.stream.network-linked vector (and any rasters) into edge records.
    Requires a live GRASS session.
    """
    cats, tostream = _read_topology_table(streams)
    geometry = _read_geometry(streams, cats)

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

    # An incomplete catchment (off-map upstream contributing area) shows up as
    # negative flow accumulation reaching a channel head. Negative accumulation
    # only at the outlet boundary (the downstream end) is the normal map exit
    # and is fine. Off-map inflow is not yet supported: fail loudly rather than
    # build a misleading network. (issue #9)
    bad = offmap_inflow_cats(records)
    if bad:
        from grass.script import fatal
        fatal("Negative flow accumulation reaches the head of segment(s) %s: "
              "the catchment is incomplete (off-map contributing area upstream; "
              "r.watershed marks such cells negative). Negative accumulation "
              "only at the outlet boundary is fine, but off-map inflow is not "
              "yet supported (enhancement: issue #9). Use a region/DEM that "
              "fully contains the catchment, or omit accumulation."
              % ', '.join(str(c) for c in bad))
    return records


def build_network(streams, elevation=None, accumulation=None, slope=None,
                  accum_mult=1.0, outlet=0):
    """Read a GRASS stream-network vector and build the NetworkX DiGraph."""
    records = read_stream_vector(streams, elevation=elevation,
                                 accumulation=accumulation, slope=slope,
                                 accum_mult=accum_mult)
    return build_graph(records, outlet=outlet)
