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
#%  description: lsdtt = DrEICH chi-z morphological channel heads (Clubb et al. 2014); slope_area = colluvial-to-fluvial (hollow) transition from the slope-area break
#%  options: slope_area,lsdtt
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

#%option G_OPT_V_OUTPUT
#%  key: output
#%  label: Output channel-head points (lsdtt) or colluvial-to-fluvial transition points (slope_area)
#%  required: yes
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
#%  label: Source drainage-area threshold, cells [method=lsdtt]
#%  answer: 100
#%  required: no
#%end

#%option
#%  key: a_0
#%  type: double
#%  label: Reference drainage area A_0 for chi [method=lsdtt]
#%  answer: 1000
#%  required: no
#%end

#%option
#%  key: m_over_n
#%  type: double
#%  label: Concavity m/n for chi [method=lsdtt]
#%  answer: 0.525
#%  required: no
#%end

#%option
#%  key: n_connecting_nodes
#%  type: integer
#%  label: Consecutive high-curvature nodes required to flag a valley [method=lsdtt]
#%  answer: 10
#%  required: no
#%end

#%option
#%  key: min_segment_length
#%  type: integer
#%  label: Minimum chi-z segment length for the channel-head split [method=lsdtt]
#%  answer: 10
#%  required: no
#%end

#%option
#%  key: window_radius
#%  type: double
#%  label: Polynomial-fit window radius (map units) for tangential curvature [method=lsdtt]
#%  answer: 7
#%  required: no
#%end

#%option
#%  key: tan_curv_threshold
#%  type: double
#%  label: Tangential-curvature threshold that marks valley cells [method=lsdtt]
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
    if options['method'] == 'lsdtt':
        _run_lsdtt(options, flags)
    else:
        _run_slope_area(options, flags)


def _run_lsdtt(options, flags):
    """DrEICH chi-z morphological channel heads (Clubb et al., 2014), via the
    faithful LSDTopoTools port in rivernetworkx.dreich. Needs only the DEM; fill,
    flow routing, tangential curvature, valleys and the chi-z split are all
    computed internally."""
    import numpy as np
    from rivernetworkx import dreich
    from rivernetworkx.grass_io import read_raster_gs

    elevation = options['elevation']
    output = options['output']
    z, region = read_raster_gs(elevation)
    nsres, ewres = region['nsres'], region['ewres']
    if abs(nsres - ewres) > 1e-6 * max(nsres, ewres):
        gscript.warning("Non-square cells (nsres=%.6g, ewres=%.6g); DrEICH assumes "
                        "square cells, using nsres." % (nsres, ewres))
    cellsize = nsres
    nodata = -9999.0
    z = np.where(np.isfinite(z), z, nodata)

    gscript.message("DrEICH (lsdtt): filling, routing, curvature, valleys, chi-z split.")
    heads = dreich.extract_channel_heads(
        z, nodata=nodata, cellsize=cellsize,
        threshold=int(options['threshold']), min_slope=float(options['min_slope']),
        A_0=float(options['a_0']), m_over_n=float(options['m_over_n']),
        n_connecting_nodes=int(options['n_connecting_nodes']),
        min_segment_length=int(options['min_segment_length']),
        window_radius=float(options['window_radius']),
        tan_curv_threshold=float(options['tan_curv_threshold']))
    if not heads:
        gscript.fatal("No channel heads found; check the DEM, threshold, or "
                      "curvature parameters.")
    gscript.message("Found %d DrEICH channel heads." % len(heads))

    west, north = region['west'], region['north']
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
    gscript.message("Wrote %d channel-head points to vector <%s>." % (len(heads), output))


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

    if not accumulation or not streams:
        gscript.fatal("method=slope_area requires both 'accumulation' and 'streams'.")

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
