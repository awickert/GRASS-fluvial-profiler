# r.fluvial: fluvial geomorphology tools for GRASS GIS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20848466.svg)](https://doi.org/10.5281/zenodo.20848466)

Fluvial geomorphology in [GRASS GIS](https://grass.osgeo.org/): river long
profiles, slope/accumulation (e.g., slope&ndash;area) relationships, river-network
graphs, channel steepness index (not yet), the hillslope-to-fluvial
(colluvial-to-fluvial) transition, and as much more as I add.

If you want a full-featured topographic-analysis suite, you should check out
[LSDTopoTools](https://lsdtopotools.github.io/) or
[TopoToolbox](https://topotoolbox.wordpress.com/).
If you want an easy way to perform simple analyses inside GRASS GIS, this is
your program.

## Modules

This repository is a small collection of [GRASS GIS addon
modules](https://grass.osgeo.org/grass-stable/manuals/addons/). They build on a
stream network extracted by
[`r.stream.extract`](https://grass.osgeo.org/grass-stable/manuals/r.stream.extract.html).

| Module | Purpose |
| --- | --- |
| **`v.stream.network`** | Adds topology columns to a stream-network vector: upstream/downstream node coordinates and `tostream`, the category of the next segment downstream (`0` if the stream leaves the map). Optionally (`json=`) exports the linked network as a [NetworkX](https://networkx.org/) node-link JSON graph, sampling elevation and flow accumulation along each segment and computing cumulative distance upstream of the outlet. |
| **`v.fluvial.profiler`** | Builds and plots river long profiles and slope&ndash;accumulation (e.g., slope&ndash;area) diagrams for a single downstream-directed channel. |
| **`r.fluvial.channelheads`** | Locates the upstream limit of the channel network, two ways. **`method=dreich`**: DrEICH morphological channel heads (Clubb et al. 2014) &mdash; a faithful port of the LSDTopoTools chi&ndash;elevation algorithm (fill &rarr; D8 routing &rarr; tangential curvature &rarr; valleys &rarr; chi&ndash;z split); needs only the DEM, optionally taking external routing (`direction=`, e.g. `r.watershed -s`) and emitting the downstream fluvial network as vector lines and/or a raster (v.stream.network directed-graph format). **`method=slope_area`** (default): the **colluvial-to-fluvial transition** (the downslope limit of colluvial hollows) from the slope&ndash;area break &mdash; a reproducible *process* boundary that lies downslope of the DrEICH head. |
| **`r.fluvial.fastscape`** | D8 flow routing (fill, direction, accumulation; route-only by default, `nsteps=0`) and, optionally, landscape evolution under rock uplift and detachment-limited stream-power incision (`dz/dt = U - K A^m S^n`) with the implicit, unconditionally stable FastScape algorithm (Braun &amp; Willett 2013). The routing feeds `r.fluvial.channelheads`. |

The graph construction, raster sampling, and JSON I/O shared by these modules
live in the [`rivernetworkx`](rivernetworkx/) Python package (pure NetworkX, no
GRASS), so the same network representation is reusable outside GRASS.

Each module has its own GRASS HTML manual page in its directory.

## Dependencies

- [GRASS GIS](https://grass.osgeo.org/) (a current release)
- Python 3 with `numpy`, `pandas`, `matplotlib`, `scipy`
- `networkx` (for the river-network graph and JSON export)
- The `rivernetworkx` package in this repository (`pip install -e .`), which the
  modules import for graph construction and JSON I/O

## Installation

Each module is a standard GRASS Python addon. From within a GRASS session you
can install one with `g.extension`, pointing its `url` at wherever the module
lives &mdash; a local path or this repository:

```
g.extension extension=v.stream.network url=<path-or-URL-to-the-module>
```

Alternatively, copy the module's `.py` file into your GRASS addons `scripts/`
directory and make it executable.

## Usage

A typical workflow extracts a stream network, links it, and then analyzes it:

```
g.region -p raster=DEM
r.mapcalc "cellArea_meters2 = nsres() * ewres()" --overwrite
r.watershed elevation=DEM flow=cellArea_meters2 accumulation=accumulation -s --overwrite
r.mapcalc "accumulation = if(isnull(DEM), null(), accumulation)" --overwrite
r.stream.extract elevation=DEM accumulation=accumulation \
    stream_raster=streams stream_vector=streams \
    threshold=30000000 direction=draindir d8cut=0 --overwrite

v.stream.network map=streams elevation=DEM accumulation=accumulation \
    json=network.json
```

See each module's manual page for full option lists and examples.

### Post-processing

`examples/clean_coarsen_network.py` is an optional "stage 2" for the JSON
exported by `v.stream.network json=`: it despikes and smooths the DEM-sampled
elevations along each segment and coarsens (resamples) the network, producing a
cleaner, thinner network suitable as input to a downstream long-profile model
such as [GRLP](https://github.com/awickert/GRLP). For example:

```
python examples/clean_coarsen_network.py network.json network_clean.json
```

Run with `--help` to see the despiking, smoothing, and coarsening options.

### Hillslope lengths to channels

Once you have a stream network and a flow-direction raster (both produced in the
workflow above),
[`r.stream.distance`](https://grass.osgeo.org/grass-stable/manuals/addons/r.stream.distance.html)
(from the `r.stream.*` addon suite) gives the **down-flow-path distance from each
cell to the channel it drains into** &mdash; the hillslope flow length &mdash;
and, with `elevation=`, the elevation drop along that path (Height Above Nearest
Drainage, HAND):

```
r.stream.distance stream_rast=streams direction=draindir elevation=DEM \
    method=downstream distance=hillslope_length difference=hand
```

`hillslope_length` is a useful input for **hydrological models** (overland-flow
travel time, runoff routing), **soil-erosion** estimates (the slope-length term
of the (R)USLE `LS` factor), and **hillslope geomorphology** (characteristic
hillslope length, drainage density). `hand` is a standard wetness/floodplain
predictor.

The channel network passed as `stream_rast` defines where hillslopes are taken to
end, so that choice propagates into every derived length. Thresholding
accumulation at `r.fluvial.channelheads`'s colluvial-to-fluvial transition area `A*`
(`accumulation >= A*`) gives a process-based channel definition rather than an
arbitrary extraction threshold.

## License

GNU General Public License, version 3 or later (GPL >= v3). See
[`LICENSE`](LICENSE).

## Author

Andrew D. Wickert
