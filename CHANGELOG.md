# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **`r.fluvial.channelheads` `method=dreich` gains network outputs**: alongside
  the channel-head points (`points`, now optional), emit the **downstream fluvial
  network** (channel heads and everything below) as vector lines (`network=`) in
  **v.stream.network format** &mdash; each link carries `x1,y1,x2,y2` endpoints
  and a `tostream` cat (0 = exits the map), a ready-to-use converging directed
  graph with no separate `v.stream.network` run &mdash; and/or as a `CELL` raster
  stream map (`raster_network=`, cell value = link cat, NULL off-network) for
  downstream GRASS modules. At least one output is required. Backed by
  `rivernetworkx.channel_network_segments`, which traces the D8 network from the
  heads and splits it into links at confluences.
- **`r.fluvial.channelheads` `method=dreich` accepts external routing**: an
  optional `direction=` drainage-direction raster (r.watershed encoding) supplies
  the D8 routing from another module instead of the internal steepest descent.
  `r.watershed -s` and `r.fluvial.fastscape` are near-canonical sources
  (`r.watershed -s` reproduces the internal heads closely; the default
  multiple-flow-direction mode does not). Tangential curvature stays internal.
- **`r.fluvial.channelheads` gains `method=dreich`**: **DrEICH morphological channel
  heads** (Clubb et al., 2014), a faithful Python/numpy port of the LSDTopoTools
  chi&ndash;elevation pipeline &mdash; depression fill, D8 routing, second-order
  polynomial tangential curvature, valley detection (sustained high curvature),
  and the per-valley hilltop&rarr;junction chi&ndash;z regression split. Needs only
  the DEM. The existing slope&ndash;area method is now `method=slope_area` (the
  default). Validated stage-by-stage against an instrumented LSDTopoTools run on
  Mid Bailey Run: routing/sources, valley junctions (768), hilltops (768/768) and
  head count (634) are bit-exact; 630/634 head locations match (the few residual
  are float32-degenerate chi&ndash;z split optima).
- **`r.fluvial.fastscape`**: a new module for **D8 flow routing** and optional
  **landscape evolution**. By default (`nsteps=0`) it routes only &mdash;
  depression fill, `direction=` (D8, r.watershed encoding), and `accumulation=`
  (upstream cell count) &mdash; a fast O(n) drainage-routing solver (the
  Braun & Willett ordered-stack machinery shared with the DrEICH port) usable as
  a routing front-end to `r.fluvial.channelheads` alongside `r.watershed`,
  `r.fill.dir`, and `r.richdem`. With `nsteps>0` it evolves the landscape under
  uplift and detachment-limited stream-power incision (`dz/dt = U - K A^m S^n`)
  with the implicit, unconditionally stable FastScape algorithm
  (Braun & Willett, 2013), then routes the evolved surface. Outputs: `direction=`,
  `accumulation=`, `filled=` (the depression-filled DEM), and the evolved
  `output=` DEM; the routing is the canonical FastScape/LSDFlowInfo routing,
  internally consistent with `r.fluvial.channelheads direction=` and
  `r.stream.distance`.
- **`rivernetworkx`** gains the **χ–z channel-head primitives** in `core`:
  `chi` (the integral channel coordinate of Perron & Royden, 2013),
  `channel_head_chi_split` (the DrEICH χ–z regression head-finder &mdash; splits a
  hilltop&rarr;junction profile into an upslope hillslope and a downslope channel),
  and `strahler_order` (Strahler stream order on the segment graph).
- **`rivernetworkx`** gains `dreich` (`extract_channel_heads`,
  `channel_network_segments`, `drainage_divides`, `build_flowinfo_from_directions`,
  `directions_from_flowinfo`, `fill`,
  `tangential_curvature`) and `fastscape` (`evolve`), plus gscript-only raster I/O
  helpers `read_raster_gs` / `read_raster_int_gs` / `write_raster_gs` (via
  `r.out.bin` / `r.in.bin`, so modules run in headless `grass --exec` sessions
  where pygrass cannot load its shared libraries; `read_raster_int_gs` reads CELL
  maps such as a drainage-direction raster).
- **`r.fluvial.channelheads`**: a new module that maps the **colluvial-to-fluvial
  transition** (the downslope limit of colluvial hollows) from the
  slope&ndash;area break. It reads a deliberately over-extracted
  `r.stream.extract` network, fits a constrained broken-stick to the
  flowline-smoothed slope&ndash;area cloud (a flat hillslope/colluvial limb and a
  power-law fluvial limb), and places a transition point where each channel's
  drainage area first crosses the break A\*. This is a complete, reproducible
  *process* boundary &mdash; the appropriate upstream limit of the fluvial domain
  &mdash; not a field channel head (e.g. Clubb et al. 2014's DrEICH head, which
  lies upslope within the colluvial hollow). Characterized against Clubb's 53
  mapped heads at Mid Bailey Run, OH (field heads sit at ~0.5&times; A\*, mostly
  in the colluvial zone).
- **`rivernetworkx`** gains `fit_sa_break` (the constrained broken-stick
  slope&ndash;area break detector), `colluvial_fluvial_transition` (place points
  where drainage area first crosses A\*), and `read_stream_segments` (read a stream
  vector into edge records **without** a `tostream` column, so a raw
  `r.stream.extract` network is usable directly &mdash; no `v.stream.network`
  linking step required).

### Changed
- **`v.stream.profiler` renamed to `v.fluvial.profiler`.** The geomorphic-analysis
  modules form a `fluvial` family on the shared `rivernetworkx` core
  (`v.fluvial.profiler` plus the new `r.fluvial.channelheads`), kept distinct from
  **`v.stream.network`**, which retains its name as the network-topology builder
  (stream-network plumbing that composes with the `r.stream.*` ecosystem). The
  `r.`/`v.` prefixes follow data type. Breaking: the `v.stream.profiler` name no
  longer exists.
- **Repository renamed `GRASS-fluvial-profiler` &rarr; `r.fluvial`.** The project
  spans the whole `r.`/`v.fluvial` module family plus the `rivernetworkx` library,
  so the repo now carries the family name. GitHub redirects the old URLs; the
  `CITATION.cff`, `pyproject.toml`, and changelog links point at the new name.

### Fixed
- **`v.stream.network`** no longer crashes on a real network that has both
  tributary confluences and an off-map outlet. The downstream-cat values come
  back from `vector_db_select` as strings, which made the `to_cat` column
  string-typed; assigning the integer `0` off-map sentinel into it then raised a
  `TypeError` under pandas' string dtype. The linked cat is now cast to `int`.
  (Synthetic test fixtures have no confluence, so this only surfaced on real
  data — e.g. the Trempealeau network.)

## [0.3.0] - 2026-06-26

This release migrates the toolkit onto a single shared river-network library,
modernizes the modules, and makes several correctness fixes. The API is still
pre-1.0 and may change.

### Added
- **`rivernetworkx`**: a new importable Python package (pure NetworkX, no GRASS)
  that single-sources the river-network logic — graph construction
  (`build_network`/`build_graph`), raster sampling, breadth-first traversal,
  cumulative-distance accumulation, per-segment smoothing and densification,
  channel slope (&minus;d&#8202;z/d&#8202;s) and slope&ndash;area, and node-link
  JSON I/O. The GRASS modules are now thin consumers of it, and the same
  representation is reusable outside GRASS (e.g. the GRLP coupling).
- **`v.stream.network`** gains an optional `json=` export (with `elevation=`,
  `accumulation=`, `accum_mult=`): it builds the NetworkX river-network graph,
  samples rasters along each segment, computes cumulative distance upstream of
  the outlet, and writes node-link JSON. This is the capability formerly
  provided by the separate `v.stream.networkx` module.
- **`v.stream.network`** and **`v.stream.profiler`** gain a `-c` flag (the region
  contains the full basin): a negative-accumulation channel head is then treated
  as an `r.watershed` boundary artifact (a divide on the region edge) rather than
  off-map inflow, skipping the off-map error.

### Changed
- **`v.stream.profiler`** rebuilt on `rivernetworkx`. It walks either downstream
  to the outlet or upstream through all tributaries, samples
  elevation/accumulation/slope along the channel, optionally densifies
  (`dx_target=`) and smooths (`window=`) **per segment** (smoothing and
  densification never cross tributary junctions, preserving the natural breaks),
  and outputs a text long-profile table and/or a node-link JSON sub-network
  (`json=`).
- Drainage area sampled from `r.watershed` is returned as a positive magnitude
  (the outlet cells, which GRASS marks negative as a boundary flag, are
  recovered).
- `examples/clean_coarsen_network.py` now post-processes the
  `v.stream.network json=` export; its accumulation array is optional.

### Fixed
- **`v.stream.profiler`** was non-functional on current GRASS GIS and pandas (a
  stale `grass.raster` import and a `Series`&rarr;`int` coercion); it runs again,
  and its long-standing branching bug is fixed by construction.
- **`v.stream.network`** writes the off-map outlet `tostream` as `0` (an earlier
  partial migration wrote `-1`); diverging (distributary/braided) networks now
  stop with a clear error instead of writing a partial `tostream`.
- Raster sampling is limited to the network's extent (covers the whole network
  even when the current region does not, and avoids reading a full large
  region), includes points exactly on the east/south region edges, and tolerates
  coincident vertices.
- Incomplete catchments (off-map upstream contributing area, flagged by negative
  `r.watershed` accumulation reaching a channel head) are caught with a clear
  error rather than silently producing a misleading network; the `-c` flag
  overrides this when the region truly holds the full basin.
- Raster NULL cells are read as NaN rather than a silent `0`, so missing data
  stays honest; sampling accumulation with no data under a channel cell now
  errors clearly instead of slipping through as a fake `0`.
- `export_json` emits standards-compliant JSON (non-finite floats as `null`, not
  the bare `NaN` token that strict parsers reject); `load_json` restores them.

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

[0.3.0]: https://github.com/awickert/r.fluvial/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/awickert/r.fluvial/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/awickert/r.fluvial/releases/tag/v0.1.0
