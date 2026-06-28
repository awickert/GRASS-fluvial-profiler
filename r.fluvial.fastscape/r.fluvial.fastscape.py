#!/usr/bin/env python
############################################################################
#
# MODULE:       r.fluvial.fastscape
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Evolve a landscape under rock uplift and detachment-limited
#               stream-power river incision (dz/dt = U - K A^m S^n) using the
#               implicit, unconditionally stable FastScape algorithm of
#               Braun & Willett (2013).
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
#% description: Evolve a landscape under uplift and stream-power incision (implicit FastScape)
#% keyword: raster
#% keyword: hydrology
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
#%  key: output
#%  label: Output evolved elevation (DEM)
#%  required: yes
#%end

#%option
#%  key: k
#%  type: double
#%  label: Stream-power erodibility coefficient K
#%  answer: 1e-5
#%  required: no
#%end

#%option
#%  key: m
#%  type: double
#%  label: Drainage-area exponent m
#%  answer: 0.5
#%  required: no
#%end

#%option
#%  key: n
#%  type: double
#%  label: Slope exponent n
#%  answer: 1.0
#%  required: no
#%end

#%option
#%  key: uplift
#%  type: double
#%  label: Rock-uplift rate U (length / time), applied to the interior
#%  answer: 1e-3
#%  required: no
#%end

#%option G_OPT_R_INPUT
#%  key: uplift_map
#%  label: Spatially-variable uplift raster (overrides 'uplift' if given)
#%  required: no
#%end

#%option G_OPT_R_INPUT
#%  key: k_map
#%  label: Spatially-variable erodibility raster (overrides 'k' if given)
#%  required: no
#%end

#%option
#%  key: dt
#%  type: double
#%  label: Timestep (consistent time units with K and uplift)
#%  answer: 1000.0
#%  required: no
#%end

#%option
#%  key: nsteps
#%  type: integer
#%  label: Number of timesteps
#%  answer: 100
#%  required: no
#%end

#%option
#%  key: min_slope
#%  type: double
#%  label: Depression-fill minimum gradient used for flow routing each step
#%  answer: 1e-4
#%  required: no
#%end

#%flag
#%  key: b
#%  description: Do NOT fix the domain edges as base level (default fixes the four edges)
#%end

import numpy as np

from grass import script as gscript


def main():
    options, flags = gscript.parser()
    elevation = options['input']
    output = options['output']
    K = float(options['k'])
    m = float(options['m'])
    n = float(options['n'])
    uplift = float(options['uplift'])
    dt = float(options['dt'])
    nsteps = int(options['nsteps'])
    min_slope = float(options['min_slope'])

    try:
        from rivernetworkx import fastscape
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

    gscript.message("FastScape: %d steps, dt=%g, K=%g, m=%g, n=%g, U=%g, cell=%g"
                    % (nsteps, dt, K, m, n, uplift, cellsize))
    zf = fastscape.evolve(z, nodata=nodata, cellsize=cellsize, K=Kp, m=m, n=n,
                          uplift=U, dt=dt, nsteps=nsteps, fixed_boundary=fixed,
                          min_slope=min_slope)

    # write the evolved DEM back, mapping nodata -> NULL
    zf_out = np.where(zf == np.float32(nodata), np.nan, zf).astype(np.float64)
    write_raster_gs(zf_out, output, region, overwrite=gscript.overwrite())
    gscript.message("Wrote evolved DEM to <%s> (relief %.2f)."
                    % (output, float(np.nanmax(zf_out))))


if __name__ == "__main__":
    main()
