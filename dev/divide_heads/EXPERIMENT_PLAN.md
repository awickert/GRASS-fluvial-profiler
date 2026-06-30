# Uniform channel-head validation experiment

One identical, scripted protocol applied to every field site, validating the
divides + chi-z channel-head method (`r.fluvial.channelheads method=divides`)
against Clubb et al. field-mapped channel heads.

**Sites** (Clubb et al. 3-site dataset, `Channel_head_coords.xlsx`):

| site | heads | UTM | DEM |
|------|------:|-----|-----|
| Mid Bailey Run, OH | 53 | 17N | `bailey_run_dem.flt` (1 m) — primary/anchor |
| Indian Creek, OH | 36 | 17N | `indian_creek_dem.flt` (1 m) |
| Feather River, CA | 15 | 10N | to source (OpenTopography / NCALM lidar) |

## Questions
- **A. Native-resolution accuracy** — does the method locate the field heads?
- **B. Resolution robustness** — how far can the DEM be coarsened before A fails?

## Fixed protocol (identical across sites)

**Inputs** (provenance recorded in the script header): DEM `.flt`+`.hdr` (georef,
res, UTM); field heads from the spreadsheet by site name.

**Evaluation domain:** convex hull of the field heads **+ 500 m buffer** (≈2× the
~235 m MBR valley spacing) so first-order catchments straddling the hull are not
truncated. All scoring is **in-hull only** — we measure where ground truth exists.

**Preprocessing:** clip → fill → D8 route, one code path, consistent nodata.

**Method:** `method=divides`, DrEICH defaults pinned (A₀=1000, m/n=0.525,
min_segment_length=10 physical, min_profile_nodes auto). One free parameter: the
valley scale **T (m²)**.

**Metrics** (all in-hull, at τ = 30 and 50 m): recall (field heads with a detected
head within τ), precision (detected heads within τ of a field head), F-score, and
the median + distribution of nearest distance.

## Decision 1 — valley scale T across sites  ✅ decided: density-matched
Drainage density differs between watersheds, so a single fixed T is unfair (it
conflates accuracy with how well one T fits each site). **Sweep T over a fixed
grid; report the density-matched operating point as the headline** — the metrics
**interpolated to the exact field-head density** (heads/km² in hull), a defined,
repeatable tie-break (not snap-to-nearest-swept-T). The full recall/precision-vs-
density curve is reported alongside. *Rationale, confirmed empirically:* fixing
T = 10000 m² understated MBR recall@50 by 0.34 (0.57 → 0.91) purely from an
off-density operating point; density control surfaces the true performance.
*(Implemented in `site_experiment.py`.)*

## Decision 2 — tolerance  ⏳ open
Absolute τ = 30/50 m is comparable across sites but the sites differ in scale.
Proposed addition: **also** report τ scaled to each site's valley radius (the
~80 m length scale from `valley_length_scale.py`, computed per site) — a
scale-free "head within X % of the local valley size." *Awaiting your call:
add it, or keep absolute-only.*

## Two arms

**Arm A — native-resolution accuracy** *(built, `site_experiment.py`)*: the
T-sweep + density-matched recall/precision/F/median, per site, with the
cross-site curve overlay (`plot_sites.py`).

**Arm B — resolution robustness** *(harness exists, `valley_structure_sweep.py`,
to generalize per site)*: repeat the density-matched validation over a coarsening
sweep (block-mean), T held **physically** constant, min_segment_length at fixed
physical length.

Grid (fine at the low end to catch the early curvature collapse, extended to
200 m to pass fully through breakdown):
```
1, 2, 3, 4, 5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90,
100, 120, 140, 160, 180, 200
```
Global-DEM landmarks flagged: 12 m (TanDEM-X), 30 m (SRTM/Copernicus), 90 m (SRTM3).

Per site report: (i) field-head recall@τ vs resolution; (ii) the
valley-structure-vs-channel-head separation (positional self-consistency of the
divides vs the 1 m valleys); (iii) a single **usable resolution** = where
recall@50 crosses 0.5. **Curvature-DrEICH baseline** at every resolution
(expected to collapse early — the contrast that motivates the method).

**Dimensionless axis (handles variable first-order valley size across landscapes):**
also plot every metric against `res / (2R)` — pixels per valley-head width, using
each site's own measured R — so the three sites are comparable and breakdown is
tested as a universal `res/(2R) = 1` law (prediction: full breakdown when one
pixel spans the whole valley head).

**Graceful-breakdown decomposition:** normalize *each field head* by its own
valley radius and re-plot detectability vs `res/(2R_i)`. If a gradual aggregate
tail **sharpens to a cliff** at `res/(2R_i) ≈ 1`, the gentleness was a
*distribution of valley sizes* (small valleys break first); if it **stays
gradual**, the *method* degrades gracefully. Separates the two causes explicitly.

## Controls and anchors (what makes it real)
- **Random-null control:** heads placed at random in-hull at matched density →
  the chance recall floor; report recall *above chance*. Fixed seed.
- **Curvature-DrEICH baseline** alongside divides everywhere.
- **Anchor-checks before trusting a new site:** reproduce the established MBR
  numbers through the harness (passing: loader 0 m diff; MBR T=10000 →
  recall@50 0.57, median 46 m), and verify a known field head maps to the
  correct cell (georef sanity).

## Pre-registration (predictions on the record, before running)
- **Indian Creek ≈ MBR** (same region, soil-mantled Ohio). → **confirmed**:
  density-matched recall@50 0.89–0.91, median ~21 m, both sites.
- **Feather River more different** (steep CA; expect a different density-matched
  T and a different usable resolution). → to test.

## Honest caveats (stated in outputs)
- Arm B coarsens by aggregating native 1 m DEMs — an *intrinsic-resolution* test,
  not real coarse sensors. Where a real coarse DEM exists (SRTM 30 m), add it as
  a gold check on ≥1 site.
- Negatives reported straight.

## Repeatability
One script + per-site config, pinned parameters, fixed seeds, versioned in
`dev/divide_heads/`. Outputs: per-site CSV/npz + figures + a combined cross-site
summary. DEM and field-head provenance in the script header.

## Status
- Arm A: built and run on MBR + Indian Creek (this session).
- Arm B: MBR done; generalize the harness to per-site.
- Feather River: awaiting DEM.
- Open: Decision 2 (tolerance); random-null control; the interpolation tie-break
  is implemented.
