# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/).

## [0.2.0] - 2026-06-24

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
- `LICENSE` (GNU GPL v2 or later).
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
  assorted one-off utilities). One standalone profiler script was moved to
  `archive/`.

## [0.1.0] - 2017-09-16

Initial tagged release (code from December 2016): early GRASS GIS modules for
extracting and profiling river networks.

[0.2.0]: https://github.com/awickert/GRASS-fluvial-profiler/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/awickert/GRASS-fluvial-profiler/releases/tag/v0.1.0
