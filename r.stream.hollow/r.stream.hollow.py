#!/usr/bin/env python
############################################################################
#
# MODULE:       r.stream.hollow
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
#% description: Map the colluvial-to-fluvial (hollow) transition from the slope-area break
#% keyword: raster
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#% keyword: hollow
#%end

#%option G_OPT_R_INPUT
#%  key: elevation
#%  label: Elevation (DEM) sampled along the network
#%  required: yes
#%end

#%option G_OPT_R_INPUT
#%  key: accumulation
#%  label: Flow-accumulation raster (e.g. from r.watershed), covering the network
#%  required: yes
#%end

#%option G_OPT_V_INPUT
#%  key: streams
#%  label: Over-extracted stream network (r.stream.extract, small threshold)
#%  description: Extract with a small threshold so it over-shoots upslope past the hollows
#%  required: yes
#%end

#%option G_OPT_V_OUTPUT
#%  key: output
#%  label: Output colluvial-to-fluvial transition points (downslope limit of hollows)
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

import os
import tempfile

import numpy as np

from grass import script as gscript


def main():
    options, flags = gscript.parser()
    elevation = options['elevation']
    accumulation = options['accumulation']
    streams = options['streams']
    output = options['output']
    window = float(options['window']) or None
    accum_mult = float(options['accum_mult'])
    min_slope = float(options['min_slope'])
    max_area = float(options['max_area']) if options['max_area'] else None

    try:
        import rivernetworkx as rnx
    except ImportError:
        gscript.fatal("r.stream.hollow requires the 'rivernetworkx' package "
                      "(pip install -e . in your GRASS Python environment).")

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
    transitions = rnx.channel_head_points(records, A_star)
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
