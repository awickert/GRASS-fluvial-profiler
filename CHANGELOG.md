# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **`rivernetworkx`**: a new importable Python package (pure NetworkX, no GRASS)
  that single-sources the river-network logic — graph construction
  (`build_network`/`build_graph`), raster sampling, breadth-first traversal,
  cumulative-distance accumulation, and node-link JSON I/O. The GRASS modules
  are now thin consumers of it, and the same representation is reusable outside
  GRASS (e.g. the GRLP coupling).

### Changed
- **`v.stream.network`** gains an optional `json=` export (with `elevation=`,
  `accumulation=`, `accum_mult=`): it builds the NetworkX river-network graph,
  samples rasters along each segment, computes cumulative distance upstream of
  the outlet, and writes node-link JSON. This is the capability formerly
  provided by the separate `v.stream.networkx` module. Diverging (distributary
  or braided) networks now stop with a clear error instead of a partial result.
- **`v.stream.profiler`** rebuilt on `rivernetworkx` (it was non-functional on
  current GRASS GIS and pandas). It now walks either downstream to the outlet or
  upstream through all tributaries, samples elevation/accumulation/slope along
  the channel, optionally densifies (`dx_target=`) and smooths (`window=`), and
  outputs a text long-profile table and/or a node-link JSON sub-network
  (`json=`). The previous branching bug is fixed by construction.

### Removed
- **`v.stream.networkx`**: retired. Its sole job (build graph &rarr; JSON) is now
  the `json=` option on `v.stream.network`, backed by `rivernetworkx`.
- **`v.stream.profiler/RiverNetwork.py`**: the bespoke network class is replaced
  by the shared `rivernetworkx` library.

## [0.2.0] - 2026-06-25

First release since 2017. It marks the modernization of the toolkit for current
GRASS GIS and the addition of NetworkX-based river-network analysis, developed
in collaboration with Fergus McNab (GFZ Potsdam) and presented at AGU 2025:

> Wickert, A. D., and F. McNab (2025), Simulating Geomorphic Evolution Through
> River Networks, EP23D-1702, *AGU Fall Meeting*, New Orleans, LA, USA.

The API is still pre-1.0 and may change.

### Added
- **`v.stream.networkx`**: a new module that builds a [NetworkX](https://networkx.org/)
  directed graph of a river network, samples elevation and flow accumulation
  along each segment, computes cumulative distance upstream of the outlet by a
  breadth-first sweep, and optionally exports the graph as JSON.
- `Makefile` and HTML manual page for `v.stream.networkx` so it builds and
  installs as a GRASS addon.
- `LICENSE` (GNU GPL v3 or later).
- `CITATION.cff` so the release is citable.
- `examples/clean_coarsen_network.py`: a post-processing ("stage 2") script
  that despikes, smooths, and coarsens a `v.stream.networkx` JSON export for
  use as input to a downstream long-profile model such as GRLP.
- Project `README` with module overview, dependencies, installation, and usage;
  this changelog.

### Changed
- **`v.stream.profiler`** substantially reworked and modernized for current
  GRASS GIS; adopts Pandas; produces river long profiles and
  slope&ndash;accumulation (e.g., slope&ndash;area) outputs.
- **`v.stream.network`** updated for current GRASS GIS; GRASS+NumPy raster
  reads modernized to a single step; `tostream` enforced as integer; the
  off-map outlet segment is coded as `0`.

### Removed
- Legacy scripts outside the fluvial-profiler scope (PRMS/MODFLOW/GSFLOW
  grid and parameter builders, an older standalone network extractor, and
  assorted one-off utilities). Standalone, multicore, and other Python-2
  development scripts were moved to `archive/`.

## [0.1.0] - 2017-09-16

Initial tagged release (code from December 2016): early GRASS GIS modules for
extracting and profiling river networks.

[0.2.0]: https://github.com/awickert/GRASS-fluvial-profiler/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/awickert/GRASS-fluvial-profiler/releases/tag/v0.1.0
