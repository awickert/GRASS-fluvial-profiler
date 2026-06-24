# GRASS-fluvial-profiler

Fluvial geomorphology in [GRASS GIS](https://grass.osgeo.org/): river long
profiles, slope/accumulation (e.g., slope&ndash;area) relationships, river-network
graphs, channel steepness index (not yet), hillslope-to-fluvial transition
(not yet), and as much more as I add.

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
| **`v.stream.network`** | Adds topology columns to a stream-network vector: upstream/downstream node coordinates and `tostream`, the category of the next segment downstream (`0` if the stream leaves the map). |
| **`v.stream.networkx`** | Builds a [NetworkX](https://networkx.org/) directed graph of the network, samples elevation and flow accumulation along each segment, computes cumulative distance upstream of the outlet, and optionally exports the graph as JSON. |
| **`v.stream.profiler`** | Builds and plots river long profiles and slope&ndash;accumulation (e.g., slope&ndash;area) diagrams for a single downstream-directed channel. |

Each module has its own GRASS HTML manual page in its directory.

## Dependencies

- [GRASS GIS](https://grass.osgeo.org/) (a current release)
- Python 3 with `numpy`, `pandas`, `matplotlib`, `scipy`
- `networkx` (for `v.stream.networkx`)

## Installation

Each module is a standard GRASS Python addon. From within a GRASS session you
can install one directly from its source directory, for example:

```
g.extension extension=v.stream.networkx url=v.stream.networkx/
```

or point `g.extension` at this repository. Alternatively, copy the module's
`.py` file into your GRASS addons `scripts/` directory and make it executable.

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

v.stream.network map=streams
v.stream.networkx streams=streams elevation=DEM accumulation=accumulation \
    units=m2 outjson=network.json
```

See each module's manual page for full option lists and examples.

## License

GNU General Public License, version 2 or later (GPL >= v2). See
[`LICENSE`](LICENSE).

## Author

Andrew D. Wickert
