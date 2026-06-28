#!/usr/bin/env python
############################################################################
#
# MODULE:       r.fluvial.channelheads
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Map the downslope limit of colluvial hollows -- the
#               colluvial-to-fluvial process transition -- from the slope-area
#               break: the drainage area at which channels become fluvial
#               (concave, S ~ A^-theta). Found from a deliberately
#               over-extracted stream network, placed where drainage area first
#               crosses that break.
#
# COPYRIGHT:    (c) 2026 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v3). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  a stream network from r.stream.extract, extracted with a SMALL
#         threshold so it over-shoots upslope past the colluvial hollows
#      -  flow accumulation (e.g. from r.watershed) covering that network
#      -  the rivernetworkx Python package
#
# This maps the COLLUVIAL-TO-FLUVIAL transition: the downslope limit of the
# colluvial hollow (the unchanneled, colluvium-mantled zero-order valley), where
# the landscape crosses from hillslope/colluvial to fluvial process dominance.
# It is a complete, reproducible PROCESS boundary -- NOT a field channel head.
# Morphological channel heads (e.g. Clubb et al. 2014's DrEICH) lie UPSLOPE of
# this transition, within the colluvial hollow.

#%module
#% description: Map channel heads / the colluvial-to-fluvial transition (DrEICH or slope-area break)
#% keyword: raster
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#% keyword: channel head
#% keyword: hollow
#%end

#%option
#%  key: method
#%  type: string
#%  label: Channel-head method
#%  description: dreich = DrEICH chi-z morphological channel heads (Clubb et al. 2014); slope_area = colluvial-to-fluvial (hollow) transition from the slope-area break
#%  options: slope_area,dreich
#%  answer: slope_area
#%  required: yes
#%end

#%option G_OPT_R_INPUT
#%  key: elevation
#%  label: Elevation (DEM)
#%  required: yes
#%end

#%option G_OPT_R_INPUT
#%  key: accumulation
#%  label: Flow-accumulation raster (e.g. from r.watershed), covering the network [method=slope_area]
#%  required: no
#%end

#%option G_OPT_V_INPUT
#%  key: streams
#%  label: Over-extracted stream network (r.stream.extract, small threshold) [method=slope_area]
#%  description: Extract with a small threshold so it over-shoots upslope past the hollows
#%  required: no
#%end

#%option G_OPT_R_INPUT
#%  key: direction
#%  label: External D8 drainage-direction raster for routing (r.watershed encoding) [method=dreich]
#%  description: Optional. If given, routing comes from this raster (e.g. r.watershed -s, or r.fluvial.fastscape); if omitted, the internal FastScape-style D8 routing is used (the canonical DrEICH pathway)
#%  required: no
#%end

#%option G_OPT_V_OUTPUT
#%  key: output
#%  label: Output channel-head points (dreich) or colluvial-to-fluvial transition points (slope_area)
#%  required: no
#%end

#%option G_OPT_V_OUTPUT
#%  key: network
#%  label: Output fluvial network (vector lines, channel heads and everything downstream) [method=dreich]
#%  description: Linked stream network in v.stream.network format (cat, x1, y1, x2, y2, tostream; tostream=0 = exits map), so it is a ready-to-use directed graph
#%  required: no
#%end

#%option G_OPT_R_OUTPUT
#%  key: network_raster
#%  label: Output fluvial network as a raster stream map (CELL; cell value = link cat, NULL off-network) [method=dreich]
#%  description: Raster form of the network for downstream GRASS modules; cell values match the network= vector cats
#%  required: no
#%end

#%option
#%  key: window
#%  type: double
#%  label: Flowline smoothing window (map units) for slope; 0 = none
#%  description: Per-segment moving average of elevation/area before computing slope; recommended (raw per-cell slope is noisy)
#%  answer: 0
#%  required: no
#%end

#%option
#%  key: accum_mult
#%  type: double
#%  label: Multiplier to convert flow accumulation to your chosen area unit
#%  answer: 1
#%  required: no
#%end

#%option
#%  key: min_slope
#%  type: double
#%  label: Minimum channel slope to include in the fit
#%  description: Drops near-flat-valley artifacts (slope ~ 0) that would bias the break
#%  answer: 0.0001
#%  required: no
#%end

#%option
#%  key: max_area
#%  type: double
#%  label: Maximum drainage area to include in the fit (optional)
#%  description: Excludes large downstream trunk reaches from the slope-area fit
#%  required: no
#%end

#%flag
#%  key: c
#%  description: Region contains the full basin: treat negative accumulation at hollow heads as a boundary artifact, not off-map inflow
#%end

#%option
#%  key: threshold
#%  type: integer
#%  label: Source drainage-area threshold, cells [method=dreich]
#%  answer: 100
#%  required: no
#%end

#%option
#%  key: a_0
#%  type: double
#%  label: Reference drainage area A_0 for chi [method=dreich]
#%  answer: 1000
#%  required: no
#%end

#%option
#%  key: m_over_n
#%  type: double
#%  label: Concavity m/n for chi [method=dreich]
#%  answer: 0.525
#%  required: no
#%end

#%option
#%  key: n_connecting_nodes
#%  type: integer
#%  label: Consecutive high-curvature nodes required to flag a valley [method=dreich]
#%  answer: 10
#%  required: no
#%end

#%option
#%  key: min_segment_length
#%  type: integer
#%  label: Minimum chi-z segment length for the channel-head split [method=dreich]
#%  answer: 10
#%  required: no
#%end

#%option
#%  key: window_radius
#%  type: double
#%  label: Polynomial-fit window radius (map units) for tangential curvature [method=dreich]
#%  answer: 7
#%  required: no
#%end

#%option
#%  key: tan_curv_threshold
#%  type: double
#%  label: Tangential-curvature threshold that marks valley cells [method=dreich]
#%  answer: 0.1
#%  required: no
#%end

import os
import tempfile

import numpy as np

from grass import script as gscript


def main():
    options, flags = gscript.parser()
    try:
        import rivernetworkx as rnx  # noqa: F401  (checked here for a clear message)
    except ImportError:
        gscript.fatal("r.fluvial.channelheads requires the 'rivernetworkx' package "
                      "(pip install -e . in your GRASS Python environment).")
    if options['method'] == 'dreich':
        _run_dreich(options, flags)
    else:
        _run_slope_area(options, flags)


def _run_dreich(options, flags):
    """DrEICH chi-z morphological channel heads (Clubb et al., 2014), via the
    faithful LSDTopoTools port in rivernetworkx.dreich. Needs only the DEM: fill,
    routing, tangential curvature, valleys and the chi-z split are all computed
    internally. Routing defaults to the internal FastScape-style D8 (the
    canonical DrEICH pathway) but can be supplied externally via direction=.
    Emits any of: channel-head points (output), the downstream fluvial network as
    vector lines (network) and/or as a raster stream map (network_raster)."""
    import numpy as np
    from rivernetworkx import dreich
    from rivernetworkx.grass_io import read_raster_gs, read_raster_int_gs

    elevation = options['elevation']
    output = options['output']
    network = options['network']
    network_raster = options['network_raster']
    if not (output or network or network_raster):
        gscript.fatal("method=dreich needs at least one output: 'output' "
                      "(channel-head points), 'network' (vector stream network) "
                      "and/or 'network_raster' (raster stream network).")
    want_fi = bool(network or network_raster)

    z, region = read_raster_gs(elevation)
    nsres, ewres = region['nsres'], region['ewres']
    if abs(nsres - ewres) > 1e-6 * max(nsres, ewres):
        gscript.warning("Non-square cells (nsres=%.6g, ewres=%.6g); DrEICH assumes "
                        "square cells, using nsres." % (nsres, ewres))
    cellsize = nsres
    nodata = -9999.0
    z = np.where(np.isfinite(z), z, nodata)

    direction = None
    if options['direction']:
        # CELL raster -> integer reader (r.out.bin float read would garble it)
        direction, dregion = read_raster_int_gs(options['direction'])
        if direction.shape != z.shape:
            gscript.fatal("direction raster <%s> (%dx%d) does not match the "
                          "elevation (%dx%d); align the region."
                          % (options['direction'], direction.shape[0],
                             direction.shape[1], z.shape[0], z.shape[1]))
        gscript.message("Routing from external direction raster <%s>."
                        % options['direction'])

    gscript.message("DrEICH: filling, routing, curvature, valleys, chi-z split.")
    result = dreich.extract_channel_heads(
        z, nodata=nodata, cellsize=cellsize,
        threshold=int(options['threshold']), min_slope=float(options['min_slope']),
        A_0=float(options['a_0']), m_over_n=float(options['m_over_n']),
        n_connecting_nodes=int(options['n_connecting_nodes']),
        min_segment_length=int(options['min_segment_length']),
        window_radius=float(options['window_radius']),
        tan_curv_threshold=float(options['tan_curv_threshold']),
        direction=direction, return_flowinfo=want_fi)
    heads, fi = (result if want_fi else (result, None))
    if not heads:
        gscript.fatal("No channel heads found; check the DEM, threshold, or "
                      "curvature parameters.")
    gscript.message("Found %d DrEICH channel heads." % len(heads))

    west, north = region['west'], region['north']

    if output:
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
        try:
            for row, col in heads:
                x = west + (col + 0.5) * ewres
                y = north - (row + 0.5) * nsres
                tmp.write("%.6f|%.6f\n" % (x, y))
            tmp.close()
            gscript.run_command(
                'v.in.ascii', input=tmp.name, output=output, format='point',
                separator='pipe', x=1, y=2, cat=0,
                columns="x double precision, y double precision",
                overwrite=gscript.overwrite(), quiet=True)
        finally:
            os.remove(tmp.name)
        gscript.message("Wrote %d channel-head points to vector <%s>."
                        % (len(heads), output))

    if want_fi:
        segments = dreich.channel_network_segments(fi, heads)
        if network:
            _write_network(segments, network, west, north, ewres, nsres)
            gscript.message("Wrote fluvial network (%d links) to vector <%s>."
                            % (len(segments), network))
        if network_raster:
            _write_network_raster(segments, network_raster, region)
            gscript.message("Wrote fluvial network (%d links) to raster <%s>."
                            % (len(segments), network_raster))


def _write_network(segments, network, west, north, ewres, nsres):
    """Write the DrEICH channel network as a vector line map in v.stream.network
    format: cat, upstream/downstream endpoints (x1,y1,x2,y2) and the cat of the
    single downstream link (tostream; 0 = exits the map). The result is a
    ready-to-use converging directed graph -- no v.stream.network run needed --
    and is also schema-compatible with it.

    Built without pygrass (headless-safe): GRASS ASCII 'standard' line primitives
    via v.in.ascii, then the attribute table via db.execute + v.db.connect."""
    def xy(row, col):
        return west + (col + 0.5) * ewres, north - (row + 0.5) * nsres

    ascii_tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    sql_tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.sql', delete=False)
    try:
        sql_tmp.write("CREATE TABLE %s (cat INTEGER, x1 DOUBLE PRECISION, "
                      "y1 DOUBLE PRECISION, x2 DOUBLE PRECISION, "
                      "y2 DOUBLE PRECISION, tostream INTEGER);\n" % network)
        for seg in segments:
            cells = seg['cells']
            pts = [xy(r, c) for r, c in cells]
            # GRASS ASCII standard line primitive: type, n_points, n_cats; the
            # vertices; then one 'layer cat' line.
            ascii_tmp.write("L %d 1\n" % len(pts))
            for x, y in pts:
                ascii_tmp.write(" %.6f %.6f\n" % (x, y))
            ascii_tmp.write(" 1 %d\n" % seg['cat'])
            x1, y1 = pts[0]
            x2, y2 = pts[-1]
            sql_tmp.write("INSERT INTO %s VALUES (%d, %.6f, %.6f, %.6f, %.6f, %d);\n"
                          % (network, seg['cat'], x1, y1, x2, y2, seg['tostream']))
        ascii_tmp.close()
        sql_tmp.close()
        gscript.run_command('v.in.ascii', flags='n', input=ascii_tmp.name,
                            output=network, format='standard',
                            overwrite=gscript.overwrite(), quiet=True)
        gscript.run_command('db.connect', flags='c', quiet=True)
        gscript.run_command('db.execute', input=sql_tmp.name)
        gscript.run_command('v.db.connect', map=network, table=network,
                            key='cat', layer=1, quiet=True)
    finally:
        os.remove(ascii_tmp.name)
        os.remove(sql_tmp.name)


def _write_network_raster(segments, network_raster, region):
    """Write the channel network as a CELL raster: each link's cells carry its
    cat (matching the network= vector), off-network cells become NULL. Headless
    via an int32 binary + r.in.bin (integer, bytes=4; the 0 background -> NULL)."""
    rows, cols = region['rows'], region['cols']
    grid = np.zeros((rows, cols), dtype='<i4')
    for seg in segments:
        cat = seg['cat']
        for (r, c) in seg['cells']:
            grid[r, c] = cat
    tmp = tempfile.NamedTemporaryFile(suffix='.bin', delete=False)
    tmp.close()
    try:
        grid.tofile(tmp.name)
        # integer data (no -f flag), bytes=4 -> CELL; the 0 background is NULL
        gscript.run_command('r.in.bin', input=tmp.name, output=network_raster,
                            bytes=4, anull=0,
                            north=region['north'], south=region['south'],
                            east=region['east'], west=region['west'],
                            rows=rows, cols=cols,
                            overwrite=gscript.overwrite(), quiet=True)
    finally:
        os.remove(tmp.name)


def _run_slope_area(options, flags):
    import rivernetworkx as rnx
    elevation = options['elevation']
    accumulation = options['accumulation']
    streams = options['streams']
    output = options['output']
    window = float(options['window']) or None
    accum_mult = float(options['accum_mult'])
    min_slope = float(options['min_slope'])
    max_area = float(options['max_area']) if options['max_area'] else None

    if not output:
        gscript.fatal("method=slope_area requires 'output' (transition points).")
    if not accumulation or not streams:
        gscript.fatal("method=slope_area requires both 'accumulation' and 'streams'.")
    if options['network'] or options['network_raster']:
        gscript.warning("network / network_raster outputs are only produced by "
                        "method=dreich; ignoring for method=slope_area.")

    # Read the over-extracted network with no topology requirement: the
    # slope-area break and the area-crossing placement are both per-segment and
    # use flow accumulation, so v.stream.network linking is not needed.
    gscript.message("Reading network and sampling elevation/accumulation.")
    records = rnx.read_stream_segments(streams, elevation=elevation,
                                       accumulation=accumulation,
                                       accum_mult=accum_mult,
                                       assume_complete=flags['c'])

    # Slope-area cloud (flowline-smoothed slope), then the constrained
    # broken-stick break A* = the colluvial-to-fluvial transition drainage area.
    logA, logS = rnx.slope_area(records, window=window, log=True)
    keep = logS >= np.log10(min_slope)
    if max_area:
        keep &= logA <= np.log10(max_area)
    logA, logS = logA[keep], logS[keep]
    if logA.size < 20:
        gscript.fatal("Too few channel points (%d) after filtering to fit a "
                      "slope-area break. Lower the extraction threshold, widen "
                      "max_area, or lower min_slope." % logA.size)

    fit = rnx.fit_sa_break(logA, logS)
    if fit is None:
        gscript.fatal("Could not fit a slope-area break (no clear "
                      "hillslope-to-fluvial transition in the data).")
    A_star = fit['A_star']
    gscript.message("Slope-area break: A* = %.1f (area units), theta = %.3f, "
                    "hillslope slope = %.4f"
                    % (A_star, fit['theta'], 10.0 ** fit['hillslope_logS']))

    # Place a transition point where each channel first reaches A* (the colluvial
    # hollow gives way to a fluvial channel).
    transitions = rnx.colluvial_fluvial_transition(records, A_star)
    if not transitions:
        gscript.fatal("No channel crossed A* = %.1f within the network; the "
                      "extraction may not reach upslope past the hollows."
                      % A_star)
    gscript.message("Found %d colluvial-to-fluvial transition points."
                    % len(transitions))

    # Write the points (with their source segment and the break area).
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    try:
        for x, y, cat in transitions:
            tmp.write("%.6f|%.6f|%d|%.6f\n" % (x, y, cat, A_star))
        tmp.close()
        gscript.run_command(
            'v.in.ascii', input=tmp.name, output=output, format='point',
            separator='pipe', x=1, y=2, cat=0,
            columns="x double precision, y double precision, "
                    "source_cat integer, a_star double precision",
            overwrite=gscript.overwrite(), quiet=True)
    finally:
        os.remove(tmp.name)

    gscript.message("Wrote %d transition points to vector <%s>."
                    % (len(transitions), output))


if __name__ == "__main__":
    main()
