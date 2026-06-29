#!/usr/bin/env python
############################################################################
#
# MODULE:       r.fluvial.fastscape
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      D8 flow routing (depression fill, flow directions, accumulation)
#               on a DEM and, optionally (nsteps>0), landscape evolution under
#               rock uplift and detachment-limited stream-power incision
#               (dz/dt = U - K A^m S^n) with the implicit, unconditionally stable
#               FastScape algorithm of Braun & Willett (2013).
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
#      -  the rivernetworkx Python package
#
# The solver reuses the Braun & Willett D8 FlowInfo machinery (receivers, the
# ordered node stack, contributing area) shared with r.fluvial.channelheads.

#%module
#% description: D8 flow routing, with optional landscape evolution under stream-power incision (implicit FastScape; Braun & Willett 2013)
#% keyword: raster
#% keyword: hydrology
#% keyword: flow direction
#% keyword: flow accumulation
#% keyword: drainage
#% keyword: geomorphology
#% keyword: landscape evolution
#% keyword: fastscape
#%end

#%option G_OPT_R_INPUT
#%  key: input
#%  label: Initial elevation (DEM)
#%  required: yes
#%end

#%option G_OPT_R_OUTPUT
#%  key: evolved
#%  label: Output evolved-elevation DEM (landscape evolution; requires nsteps>0)
#%  required: no
#%  guisection: Landscape evolution
#%end

#%option G_OPT_R_OUTPUT
#%  key: direction
#%  label: Output D8 flow-direction raster (r.watershed encoding)
#%  description: The canonical FastScape/LSDFlowInfo routing, for r.fluvial.channelheads direction= or r.stream.distance (1-8 CCW from NE; 0 = pit/outlet)
#%  required: no
#%end

#%option G_OPT_R_OUTPUT
#%  key: accumulation
#%  label: Output flow-accumulation raster (number of upstream cells; always positive)
#%  required: no
#%end

#%option G_OPT_R_OUTPUT
#%  key: filled
#%  label: Output depression-filled DEM used for routing
#%  required: no
#%end

#%option
#%  key: nsteps
#%  type: integer
#%  label: Number of evolution timesteps (0 = route only, no landscape evolution)
#%  answer: 0
#%  required: no
#%end

#%option
#%  key: min_slope
#%  type: double
#%  label: Depression-fill minimum gradient used for flow routing
#%  answer: 1e-4
#%  required: no
#%end

#%option
#%  key: k
#%  type: double
#%  label: Stream-power erodibility coefficient K
#%  answer: 1e-5
#%  required: no
#%  guisection: Landscape evolution
#%end

#%option
#%  key: m
#%  type: double
#%  label: Drainage-area exponent m
#%  answer: 0.5
#%  required: no
#%  guisection: Landscape evolution
#%end

#%option
#%  key: n
#%  type: double
#%  label: Slope exponent n
#%  answer: 1.0
#%  required: no
#%  guisection: Landscape evolution
#%end

#%option
#%  key: uplift
#%  type: double
#%  label: Rock-uplift rate U (length / time), applied to the interior
#%  answer: 1e-3
#%  required: no
#%  guisection: Landscape evolution
#%end

#%option G_OPT_R_INPUT
#%  key: uplift_map
#%  label: Spatially-variable uplift raster (overrides 'uplift' if given)
#%  required: no
#%  guisection: Landscape evolution
#%end

#%option G_OPT_R_INPUT
#%  key: k_map
#%  label: Spatially-variable erodibility raster (overrides 'k' if given)
#%  required: no
#%  guisection: Landscape evolution
#%end

#%option
#%  key: dt
#%  type: double
#%  label: Timestep (consistent time units with K and uplift)
#%  answer: 1000.0
#%  required: no
#%  guisection: Landscape evolution
#%end

#%flag
#%  key: b
#%  description: Do NOT fix the domain edges as base level (default fixes the four edges)
#%  guisection: Landscape evolution
#%end

import numpy as np

from grass import script as gscript


def main():
    options, flags = gscript.parser()
    elevation = options['input']
    evolved = options['evolved']
    direction = options['direction']
    accumulation = options['accumulation']
    filled_out = options['filled']
    K = float(options['k'])
    m = float(options['m'])
    n = float(options['n'])
    uplift = float(options['uplift'])
    dt = float(options['dt'])
    nsteps = int(options['nsteps'])
    min_slope = float(options['min_slope'])

    if not (evolved or direction or accumulation or filled_out):
        gscript.fatal("No output requested: set at least one of evolved=, "
                      "direction=, accumulation=, filled=.")
    if evolved and nsteps == 0:
        gscript.warning("evolved= requested with nsteps=0: no evolution occurs, so "
                        "the output is just the input DEM. Set nsteps>0 to evolve.")

    try:
        from rivernetworkx import fastscape, dreich
        from rivernetworkx.grass_io import read_raster_gs, write_raster_gs
    except ImportError:
        gscript.fatal("r.fluvial.fastscape requires the 'rivernetworkx' package "
                      "(pip install -e . in your GRASS Python environment).")

    z, region = read_raster_gs(elevation)
    nsres, ewres = region['nsres'], region['ewres']
    if abs(nsres - ewres) > 1e-6 * max(nsres, ewres):
        gscript.warning("Non-square cells (nsres=%.6g, ewres=%.6g); FastScape "
                        "assumes square cells, using nsres." % (nsres, ewres))
    cellsize = nsres

    nodata = -9999.0
    z = np.where(np.isfinite(z), z, nodata).astype(np.float64)

    U = uplift
    if options['uplift_map']:
        um, _ = read_raster_gs(options['uplift_map'])
        U = np.where(np.isfinite(um), um, 0.0)
    Kp = K
    if options['k_map']:
        km, _ = read_raster_gs(options['k_map'])
        Kp = np.where(np.isfinite(km), km, K)

    fixed = None if flags['b'] else 'edges'

    if nsteps > 0:
        gscript.message("FastScape: %d steps, dt=%g, K=%g, m=%g, n=%g, U=%g, cell=%g"
                        % (nsteps, dt, K, m, n, uplift, cellsize))
    else:
        gscript.message("FastScape flow routing only (nsteps=0), cell=%g" % cellsize)
    # nsteps=0 returns the input unchanged: the routing outputs then describe the
    # input DEM (route-only); nsteps>0 routes the evolved landscape.
    zf = fastscape.evolve(z, nodata=nodata, cellsize=cellsize, K=Kp, m=m, n=n,
                          uplift=U, dt=dt, nsteps=nsteps, fixed_boundary=fixed,
                          min_slope=min_slope)

    if evolved:
        zf_out = np.where(zf == np.float32(nodata), np.nan, zf).astype(np.float64)
        write_raster_gs(zf_out, evolved, region, overwrite=gscript.overwrite())
        what = "evolved DEM" if nsteps > 0 else "input DEM unchanged (nsteps=0)"
        gscript.message("Wrote %s to <%s>." % (what, evolved))

    # Routing outputs from the final surface (== the input DEM when nsteps=0),
    # as r.watershed-encoded routing consumable by r.fluvial.channelheads direction=.
    if direction or accumulation or filled_out:
        filled = dreich.fill(zf, nodata, min_slope, cellsize)
        if filled_out:
            fout = np.where(filled == np.float32(nodata), np.nan, filled).astype(np.float64)
            write_raster_gs(fout, filled_out, region, overwrite=gscript.overwrite())
            gscript.message("Wrote filled DEM to <%s>." % filled_out)
        if direction or accumulation:                      # FlowInfo only needed for these
            fi = dreich.build_flowinfo(filled, nodata, cellsize)
            if direction:
                _write_direction(dreich.directions_from_flowinfo(fi), direction, region)
                gscript.message("Wrote D8 flow-direction raster to <%s>." % direction)
            if accumulation:
                dreich.contributing_area(fi)
                acc = np.full(z.shape, np.nan)
                acc[fi['row_of'], fi['col_of']] = fi['ncontrib'].astype(np.float64)
                write_raster_gs(acc, accumulation, region, overwrite=gscript.overwrite())
                gscript.message("Wrote flow-accumulation raster to <%s>." % accumulation)


def _write_direction(dirarr, name, region):
    """Write an int32 direction array as a CELL raster (-1 -> NULL), headless via
    r.in.bin (gscript only, no pygrass)."""
    import os
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix='.bin', delete=False)
    tmp.close()
    try:
        dirarr.astype('<i4').tofile(tmp.name)
        gscript.run_command('r.in.bin', input=tmp.name, output=name, bytes=4, anull=-1,
                            north=region['north'], south=region['south'],
                            east=region['east'], west=region['west'],
                            rows=region['rows'], cols=region['cols'],
                            overwrite=gscript.overwrite(), quiet=True)
    finally:
        os.remove(tmp.name)


if __name__ == "__main__":
    main()
