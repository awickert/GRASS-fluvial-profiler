#!/usr/bin/env python
############################################################################
#
# MODULE:       r.stream.channelheads
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Locate fluvial channel heads as the slope-area break: the
#               drainage area at which channels become fluvial (concave,
#               S ~ A^-theta), found from a deliberately over-extracted stream
#               network and placed where drainage area first crosses that break.
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
#         threshold so it over-shoots upslope past the real channel heads
#      -  flow accumulation (e.g. from r.watershed) covering that network
#      -  the rivernetworkx Python package
#
# This finds FLUVIAL-INITIATION heads (where channels become fluvial), which is
# a complete, reproducible quantity -- it marks every fluvial channel head, so
# it legitimately returns more points than a sparse field-mapped set. It is NOT
# (yet) Clubb et al. (2014)'s DrEICH field-equivalent head, which sits somewhat
# further upslope; that is a planned second method.

#%module
#% description: Find fluvial channel heads from the slope-area break
#% keyword: raster
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#% keyword: channel head
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
#%  description: Extract with a small threshold so it over-shoots above the heads
#%  required: yes
#%end

#%option G_OPT_V_OUTPUT
#%  key: output
#%  label: Output channel-head points
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
#%  description: Region contains the full basin: treat negative accumulation at channel heads as a boundary artifact, not off-map inflow
#%end

import os
import sys
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
        gscript.fatal("r.stream.channelheads requires the 'rivernetworkx' "
                      "package (pip install -e . in your GRASS Python "
                      "environment).")

    # Read the over-extracted network with no topology requirement: the
    # slope-area break and the area-crossing head placement are both per-segment
    # and use flow accumulation, so v.stream.network linking is not needed.
    gscript.message("Reading network and sampling elevation/accumulation.")
    records = rnx.read_stream_segments(streams, elevation=elevation,
                                       accumulation=accumulation,
                                       accum_mult=accum_mult,
                                       assume_complete=flags['c'])

    # Slope-area cloud (flowline-smoothed slope), then the constrained
    # broken-stick break A* = the fluvial-initiation drainage area.
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

    # Place a head where each channel first reaches A* (becomes fluvial).
    heads = rnx.channel_head_points(records, A_star)
    if not heads:
        gscript.fatal("No channel crossed A* = %.1f within the network; the "
                      "extraction may not reach above the heads." % A_star)
    gscript.message("Found %d fluvial channel heads." % len(heads))

    # Write the points (with their source segment and the break area).
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    try:
        for x, y, cat in heads:
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

    gscript.message("Wrote %d channel heads to vector <%s>." % (len(heads), output))


if __name__ == "__main__":
    main()
