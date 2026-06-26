# Design — channel network → hillslope geometry → hydrologic recession priors

*Thinking-piece, 2026-06-24. How GRASS-fluvial-profiler (GFP) can serve, from
one shared channel-network backbone, **both** geomorphic modeling and a
**hydrologic "recession-prior factory"** for lumped watershed models (MNiShed).
Aimed at framing a coherent **v1.0**. Not a spec.*

## The unifying idea

Channel-network **geometry sets the hydrologic response.** Two classic strands:

- **Surface/channel travel time** — the geomorphologic IUH (Rodríguez-Iturbe &
  Valdés 1979) and the **width function**: the distribution of flow distances is
  the travel-time distribution.
- **Subsurface/recession** — the hillslope-storage Boussinesq result: the
  **hillslope length `L`** (lateral distance to the channel) and slope `β`, with
  soil `K, f, H`, set the recession timescale **and** its nonlinear exponent.
  (See the MNiShed-side note `recession-from-hillslope-geometry`; advection
  `τ_adv ≈ L·f/(K·sinβ)`, diffusion `τ_diff ≈ L²·f/(K·H)`, ratio = hillslope
  Péclet; nonlinear Boussinesq → power-law recession `−dQ/dt = a·Q^b`, where `b`
  is MNiShed's `recession_exponent`: 1 = linear/sloping-kinematic, 3/2 = late
  Boussinesq, 3 = early.)

Both strands share one field: **distance to the nearest channel.** `L` is just
the inverse of **drainage density** (`d ≈ 1/(2·L̄)`). So a single off-channel
distance computation feeds the geomorphic and the hydrologic uses alike.

## What GFP already provides (the on-channel, geomorphic backbone)

- **`r.stream.extract`** (upstream of GFP) → the stream raster + vector +
  flow-direction (`draindir`).
- **`v.fluvial.network`** → topology (`tostream`, node coordinates).
- **`v.stream.networkx`** → a NetworkX directed graph with per-segment elevation,
  accumulation, and along-channel distances; cumulative distance upstream of the
  outlet; JSON export.
  *(v0.3.0: `v.stream.networkx` retired — this build-graph + JSON capability is now `v.fluvial.network json=`, backed by the `rivernetworkx` library.)*
- **`v.fluvial.profiler`** → long profiles and **slope–area** diagrams.

This is the geomorphic toolset (profiles, steepness/χ, knickpoints, slope–area)
**and** it is what *defines the channel network extent* — which is the geometry
both sides depend on.

## The gap for hydrologic priors (the off-channel, per-cell geometry)

The graph modules operate *along channels*; recession priors need *across
hillslopes*. Two pieces:

1. **Distance to nearest channel + hillslope slope — don't reinvent.** GRASS's
   **`r.stream.distance`** gives, per cell, the **flow distance to the stream
   `L`** and the **elevation drop**, from which `β = atan(drop/L)`. It uses the
   `draindir` already produced by `r.stream.extract`. This is mostly *wiring*,
   not new code.
2. **Where do channels begin? (the channel head / hillslope–fluvial
   transition).** `L` depends entirely on the channel-network *extent*, which
   should be a **process threshold**, not an arbitrary accumulation cutoff:
   a slope–area (`A·S^θ`) channel-initiation criterion (Montgomery & Dietrich
   1988/1992; Tarboton). GFP already computes slope–area in `v.fluvial.profiler`,
   so it is positioned to set this physically. **This is the README's
   "hillslope-to-fluvial transition (not yet)" item — and it is the shared
   linchpin: it defines the channel extent for the geomorphic analysis *and*
   sets `L` for the hydrologic priors.**

## The "recession-prior factory" workflow

```
DEM
 └─ r.watershed                         → accumulation, drainage direction
 └─ r.stream.extract                    → streams (raster+vector), draindir
       └─ [channel-head threshold]       ← slope–area (v.fluvial.profiler) — sets extent
 ├─ v.fluvial.network / v.stream.networkx → channel graph        (GEOMORPHIC track)
 └─ r.stream.distance                    → L (flow dist), drop → β   (per cell)
 + r.in.polaris                          → ksat, theta_s/theta_r, ... (per cell ± p5/p95)
       └─ (texture→Ksat pedotransfer where outside CONUS / SoilGrids)
 → r.mapcalc                             → per cell: τ_adv, τ_diff, Péclet, expected b
 → zonal stats over basin / sub-catchment polygons
                                          → per zone: τ, b PRIORS ± spread
 → export table (CSV / JSON)             → MNiShed config H0/bounds or run_and_score bounds
```

Two reasons this is more than a mean:

- **Heterogeneity is signal.** The *distribution* of `L` (and of `K`) across a
  zone predicts the recession **exponent `b`** — power-law recession can arise
  from heterogeneity of many linear stores (Harman, Sivapalan & Kumar 2009), not
  only from Boussinesq nonlinearity. So the **distance-to-channel histogram is a
  forward predictor of `b`**, computed for free from the same field.
- **Uncertainty propagates.** POLARIS p5/p50/p95 + the `L`/`β` distribution → a
  **range** on each prior (not a single number) — exactly what the calibration
  side wants for identifiability-aware priors.

## Architecture: a prior factory MNiShed consumes

`GFP + r.in.polaris + r.stream.distance` is the **producer**: DEM (+ zone
polygons) → per-zone recession priors `{τ, b} ± bounds`. **MNiShed consumes
scalars** — no GIS dependency inside the model. Same producer/consumer split
already chosen for the POLARIS importer and the open-water Penman package. The
eventual **handoff script** is then well-scoped:

> **`prior_factory`**: input a DEM (+ basin/sub-catchment polygons, optional
> substrate overrides); output a table of recession priors (`τ`, `b`, ± p5/p95)
> per zone, ready to drop into a MNiShed config or as `run_and_score` bounds.

## Serving both research lines toward v1.0

The same channel network and channel-head threshold serve two tracks:

| Track | Uses | Existing GFP pieces | New pieces |
|---|---|---|---|
| **Geomorphic** | profiles, k_sn / χ steepness, knickpoints, slope–area, drainage density, landscape-evolution setup | `v.fluvial.network(x)`, `v.fluvial.profiler` | (steepness/χ are on the README roadmap) |
| **Hydrologic** | hillslope `L`/`β` field → recession `τ, b` priors | slope–area (defines channels) | wire `r.stream.distance`; channel-head module; prior-export |

So bringing both "hydrologic + geomorphic research into the pipeline" for **v1.0**
is really: (1) finish the **channel-head / hillslope–fluvial-transition** module
(serves both); (2) wire **`r.stream.distance`** for the per-cell hillslope field;
(3) add a thin **prior-export** step (zonal stats → recession-prior table). The
geomorphic steepness/χ tools and the hydrologic prior factory then ride the same
backbone.

## Open questions

1. **A new module vs. a workflow script.** Channel-head determination wants to be
   a reusable module (e.g. `r.stream.channelheads` / `v.stream.channelheads`);
   the prior factory may be better as a documented workflow script that composes
   modules. Where to draw the line for v1.0?
2. **Per-cell-then-aggregate vs. aggregate-then-compute.** Computing `τ(L,β,K,H)`
   per cell and then taking the zonal distribution (mean + spread) is the honest
   path (captures heterogeneity → `b`); decide whether v1.0 does the full
   per-cell propagation or a first-cut mean ± simple spread.
3. **`H` (saturated thickness).** POLARIS/SoilGrids give the soil zone to ~2 m;
   the deeper regolith/aquifer thickness that matters for `τ_diff` needs geology
   or well data. v1.0 could expose `H` as a user input/override with the soil
   depth as default.
4. **Coupling to MNiShed.** Just emit a priors table (loose coupling), or provide
   a small reader on the MNiShed side? (Lean: emit a table; keep MNiShed
   GIS-free.)

## References

- Rodríguez-Iturbe & Valdés (1979) — geomorphologic IUH.
- Montgomery & Dietrich (1988, 1992) — channel initiation, slope–area thresholds.
- Brutsaert & Nieber (1977); Troch, van Loon & Hilberts (2003, HSB); Berne,
  Uijlenhoet & Troch (2005, hillslope Péclet); Rupp & Selker (2006, sloping
  aquifers); Harman, Sivapalan & Kumar (2009, power-law recession from
  heterogeneity). (See the MNiShed `recession-from-hillslope-geometry` note.)
- GRASS: `r.stream.extract`, `r.stream.distance`, `r.watershed`; the MNiMORPH
  `r.in.polaris` / `r.in.soilgrids` soil importers.
