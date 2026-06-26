#!/usr/bin/env python
############################################################################
#
# MODULE:       v.stream.profiler
#
# AUTHOR(S):    Andrew Wickert
#
# PURPOSE:      Build long profiles and slope--accumulation (e.g.,
#               slope--area) diagrams of a river network
#
# COPYRIGHT:    (c) 2016-2026 Andrew Wickert
#
#               This program is free software under the GNU General Public
#               License (>=v3). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################
#
# REQUIREMENTS:
#      -  uses inputs from r.stream.extract and v.stream.network
#      -  the rivernetworkx Python package (pip install -e . in the GRASS env)

# More information
# Started 14 October 2016; rebuilt on rivernetworkx 2026

#%module
#% description: Build river long profiles and slope-accumulation (e.g., slope-area) diagrams
#% keyword: vector
#% keyword: stream network
#% keyword: hydrology
#% keyword: geomorphology
#%end
#%option
#%  key: cat
#%  label: Starting line segment category
#%  required: yes
#%  guidependency: layer,column
#%end
#%option G_OPT_V_INPUT
#%  key: streams
#%  label: Vector stream network, linked by v.stream.network
#%  required: yes
#%end
#%option G_OPT_V_OUTPUT
#%  key: outstream
#%  label: Vector output of the extracted sub-network
#%  required: no
#%end
#%option
#%  key: direction
#%  type: string
#%  label: Which directon to march: up or down
#%  options: upstream,downstream
#%  answer: downstream
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: elevation
#%  label: Topography (DEM)
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: accumulation
#%  label: Flow accumulation raster
#%  required: no
#%end
#%option G_OPT_R_INPUT
#%  key: slope
#%  label: Map of slope created by r.slope.area
#%  required: no
#%end
#%option
#%  key: units
#%  type: string
#%  label: Flow accumulation units
#%  options: m2, km2, cumecs, cfs
#%  required: no
#%end
#%option
#%  key: accum_mult
#%  type: double
#%  label: Multiplier to convert flow accumulation to your chosen unit
#%  answer: 1
#%  required: no
#%end
#%option
#%  key: dx_target
#%  type: double
#%  label: Target distance between output stream points [map units]
#%  required: no
#%end
#%option
#%  key: window
#%  type: double
#%  label: Moving-average distance for smoothing [map units]
#%  required: no
#%end
#%option
#%  key: plots
#%  type: string
#%  label: Plots to generate
#%  options: LongProfile,SlopeAccum,SlopeDistance,AccumDistance
#%  required: no
#%  multiple: yes
#%end
#%option G_OPT_F_OUTPUT
#%  key: outfile
#%  label: Output text table of the long profile(s)
#%  required: no
#%end
#%option G_OPT_F_OUTPUT
#%  key: json
#%  label: Output NetworkX node-link JSON of the analyzed sub-network
#%  required: no
#%end

#%flag
#%  key: c
#%  description: Region contains the full basin: treat negative accumulation at channel heads as a boundary artifact, not off-map inflow
#%end

##################
# IMPORT MODULES #
##################
# PYTHON
import numpy as np
from matplotlib import pyplot as plt
# GRASS
from grass.pygrass.modules.shortcuts import vector as v
from grass import script as gscript
# rivernetworkx: the shared river-network library (graph, sampling, JSON I/O,
# and the profile-assembly / densify / smoothing helpers).
import rivernetworkx as rnx


###################
# UTILITY         #
###################

def accumulation_label(units):
    """Axis label for the flow-accumulation/discharge quantity."""
    return {'m2': 'Drainage area [m$^2$]',
            'km2': 'Drainage area [km$^2$]',
            'cumecs': 'Water discharge [m$^3$ s$^{-1}$]',
            'cfs': 'Water discharge [cfs]'}.get(units, 'Flow accumulation [$-$]')


def selected_paths(G, start, direction):
    """
    Cats making up the profile(s), as a list of downstream-ordered paths
    (each [most upstream, ..., most downstream], off-map outlet excluded).

    downstream: a single path from ``start`` to the map outlet.
    upstream:   one path per headwater draining to ``start``, ending at ``start``.
    """
    if direction == 'upstream':
        sub = rnx.upstream_subnetwork(G, start)
        headwaters = [n for n in sub.nodes if sub.in_degree(n) == 0]
        return [rnx.downstream_path(sub, hw, outlet=start) for hw in headwaters]
    # downstream
    return [rnx.downstream_path(G, start, outlet=rnx.OFFMAP)[:-1]]


def build_profile(records_by_cat, path, attrs, dx_target, window):
    """
    Assemble one continuous profile along ``path``. Densification and smoothing
    are done **per segment** (never across tributary junctions): each segment is
    resampled to ``dx_target`` and its quantities smoothed over ``window`` on
    its own along-distance, then the segments are concatenated. Returns a dict
    of 1-D arrays: ``s`` (distance from the path's downstream end), ``x``,
    ``y``, each attr, and ``<attr>_smoothed`` when ``window`` is set.

    In very dense networks some inter-junction reaches are too short to densify
    or smooth meaningfully; those reaches simply come through near-raw, which is
    accepted rather than smoothing across the junctions.
    """
    processed = {}
    for cat in path:
        rec = dict(records_by_cat[cat])
        if dx_target is not None:
            s_down, _ = rnx.segment_distances(rec['x'], rec['y'])
            dens_in = {name: rec[name] for name in (['x', 'y'] + list(attrs))}
            _, dens = rnx.densify(s_down, dens_in, dx_target)
            rec.update(dens)
        if window is not None:
            # !!!! TODO: switches to smooth z / slope / accumulation
            # individually (one window currently smooths them all) !!!!
            for name, arr in rnx.smooth_segment(rec, attrs, window).items():
                rec[name + '_smoothed'] = arr
        processed[cat] = rec
    out_attrs = list(attrs) + ([a + '_smoothed' for a in attrs] if window else [])
    return rnx.assemble_downstream_profile(processed, path, attrs=out_attrs)


###############
# MAIN MODULE #
###############

def main():
    options, flags = gscript.parser()

    start = int(options['cat'])
    streams = options['streams']
    direction = options['direction'] or 'downstream'
    elevation = options['elevation'] or None
    accumulation = options['accumulation'] or None
    slope = options['slope'] or None
    accum_mult = float(options['accum_mult'])
    outstream = options['outstream'] or None
    outfile = options['outfile'] or None
    outjson = options['json'] or None
    plots = [p for p in options['plots'].split(',') if p]
    accum_label = accumulation_label(options['units'])
    dx_target = float(options['dx_target']) if options['dx_target'] else None
    window = float(options['window']) if options['window'] else None

    # Read the linked network into edge records (with rasters sampled along each
    # segment) and build the graph. read_stream_vector errors helpfully if the
    # vector has not been linked by v.stream.network (no 'tostream' column).
    records = rnx.read_stream_vector(streams, elevation=elevation,
                                     accumulation=accumulation, slope=slope,
                                     accum_mult=accum_mult,
                                     assume_complete=flags['c'])
    records_by_cat = {rec['cat']: rec for rec in records}
    G = rnx.build_graph(records)
    if start not in G:
        gscript.fatal("Starting cat=%d is not a segment in <%s>."
                      % (start, streams))

    attrs = [name for name, on in (('z', elevation), ('A', accumulation),
                                   ('slope', slope)) if on]

    gscript.message("Extracting %s drainage pathway(s)..." % direction)
    paths = selected_paths(G, start, direction)
    profiles = [build_profile(records_by_cat, path, attrs, dx_target, window)
                for path in paths]
    selected_cats = sorted({cat for path in paths for cat in path})

    # Optional: extract the analyzed sub-network as its own vector
    if outstream:
        v.extract(input=streams, output=outstream,
                  cats=','.join(str(c) for c in selected_cats),
                  overwrite=gscript.overwrite())

    # Optional: export the analyzed sub-network as node-link JSON
    if outjson:
        if direction == 'upstream':
            H = rnx.upstream_subnetwork(G, start)
        else:
            H = G.subgraph(set(selected_cats) | {rnx.OFFMAP}).copy()
        rnx.export_json(H, outjson)
        gscript.message("Wrote JSON sub-network: %s" % outjson)

    # Warn (rather than silently skip) when a requested plot lacks its input.
    required = {'LongProfile': elevation, 'SlopeAccum': slope and accumulation,
                'SlopeDistance': slope, 'AccumDistance': accumulation}
    for p in plots:
        if not required.get(p, True):
            gscript.warning("Plot '%s' needs its input raster(s) "
                            "(elevation / accumulation / slope); skipping." % p)

    # Plots
    smooth = window is not None

    def pick(prof, name):
        """The smoothed series if a window was set, else the raw one."""
        return prof[name + '_smoothed'] if smooth else prof[name]

    if 'LongProfile' in plots and elevation:
        plt.figure()
        for prof in profiles:
            plt.plot(prof['s'] / 1000., pick(prof, 'z'), 'k-', linewidth=2)
        plt.xlabel('Distance from mouth [km]', fontsize=16)
        plt.ylabel('Elevation [m]', fontsize=16)
        plt.tight_layout()
    if 'SlopeAccum' in plots and slope and accumulation:
        plt.figure()
        for prof in profiles:
            S, A = pick(prof, 'slope'), pick(prof, 'A')
            keep = A > 0
            plt.loglog(A[keep], S[keep], 'k.', alpha=.5)
        plt.xlabel(accum_label, fontsize=16)
        plt.ylabel('Slope [$-$]', fontsize=16)
        plt.tight_layout()
    if 'SlopeDistance' in plots and slope:
        plt.figure()
        for prof in profiles:
            plt.plot(prof['s'] / 1000., pick(prof, 'slope'), 'k-', linewidth=2)
        plt.xlabel('Distance from mouth [km]', fontsize=16)
        plt.ylabel('Slope [$-$]', fontsize=16)
        plt.tight_layout()
    if 'AccumDistance' in plots and accumulation:
        plt.figure()
        for prof in profiles:
            A = pick(prof, 'A')
            keep = A > 0
            plt.plot(prof['s'][keep] / 1000., A[keep], 'k.', alpha=.5)
        plt.xlabel('Distance from mouth [km]', fontsize=16)
        plt.ylabel(accum_label, fontsize=16)
        plt.tight_layout()
    if plots:
        plt.show()

    # Text long-profile table: one row per point, labelled by its source
    # segment (cat). In upstream mode the headwater->start paths share their
    # trunk, so de-duplicate by cat to write each segment once (a segment's
    # distance from the start is path-independent in a converging network).
    if outfile:
        cols = ['cat', 's', 'x', 'y'] + attrs
        if smooth:
            cols += [a + '_smoothed' for a in attrs]
        seen = set()
        rows = []
        for prof in profiles:
            block = [prof['cat'], prof['s'], prof['x'], prof['y']]
            block += [prof[a] for a in attrs]
            if smooth:
                block += [prof[a + '_smoothed'] for a in attrs]
            block = np.column_stack(block)
            mask = np.array([int(c) not in seen for c in prof['cat']])
            seen.update(int(c) for c in prof['cat'])
            if mask.any():
                rows.append(block[mask])
        # cat is an integer segment id (the join key back to the vector); the
        # rest are floats. np.column_stack upcast cat to float, so format the
        # first column as %d to avoid writing it as e.g. "37.0".
        fmt = ['%d'] + ['%s'] * (len(cols) - 1)
        np.savetxt(outfile, np.vstack(rows), fmt=fmt,
                   header=' '.join(cols), comments='')
        gscript.message("Wrote long-profile table: %s" % outfile)


if __name__ == "__main__":
    main()
