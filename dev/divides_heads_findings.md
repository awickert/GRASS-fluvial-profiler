# Channel heads from drainage divides + chi-z (replacing cross-valley curvature)

Test site: 1 m USGS LiDAR, Mid Bailey Run, OH (the DrEICH validation DEM).
Reference: the faithful curvature-DrEICH port (634 heads at 1 m) and Clubb et al.
(2014)'s 53 field-mapped channel heads. Scripts: `divides_setup.py` (cache),
`divides_lib.py` (shared), `divides_valleys_v2.py` (1 m T-sweep),
`divides_resolution.py` (resolution sweep), `divides_vs_clubb.py` (field heads),
`dreich_ridge_survival.py` (premise).

## Idea

Cross-valley **tangential curvature** is the only resolution-fragile step in
DrEICH: it is the gate (`find_valleys`) that decides which junctions get a chi-z
head search. Curvature collapses under coarsening (window growth + DEM averaging
smooth away the few-metre convergent fingertip): head recall 0.49 @2 m -> ~0 @>=10 m
(`dreich_resolution_findings.md`). **Drainage divides** are basin boundaries -
topological, not differential - and survive coarsening. Replace the curvature
gate with the divide-defined valley set; keep the chi-z head finder unchanged.

## Premise: ridges survive where curvature heads do not (`dreich_ridge_survival.py`)

Ridge position recall vs the 1 m ridge, at a fixed 50000-cell basin scale:

| res (m) | 1 | 2 | 3 | 5 | 8 | 10 | 15 | 20 | 30 |
|--|--|--|--|--|--|--|--|--|--|
| ridge recall | 1.00 | 0.95 | 0.92 | 0.87 | 0.85 | 0.85 | 0.85 | 0.84 | 0.84 |

Ridge recall holds >=0.84 to 30 m; curvature head recall was ~0 by 10 m. The
premise holds. (Ridge *precision* falls 1.00->0.63 as divides fatten with
coarsening - density rises 1.8%->40% of cells - but position is preserved.)

## Realization (one false start, then the right one)

**Wrong:** a per-cell "flanked by divides on both sides" gate inside the
`find_valleys` walk. It fires high near the ridge crest and tags tiny upstream
junctions -> heads at median area 59 cells vs the reference's 6305, recall@30 0.01.
Abandoned.

**Right:** the divides define the first-order **valleys** directly. The channel
network at a valley-scale area threshold `T` partitions the terrain into
divide-bounded first-order basins; each basin's T-source is one valley, and one
chi-z head is found on its ridgetop->T-source profile. `T` is the valley-scale
knob (robust to coarsening); chi-z localizes the head within each valley. No
curvature anywhere.

## 1 m: valley scale T vs the curvature reference (634 heads, med area 6305)

| T (m^2) | heads | head med area | R@30 | R@50 | P@30 | P@50 |
|--:|--:|--:|--:|--:|--:|--:|
| 5000  | 1897 | 2990  | 0.42 | 0.64 | 0.14 | 0.22 |
| 10000 | 943  | 6206  | 0.33 | 0.51 | 0.22 | 0.31 |
| 20000 | 472  | 11501 | 0.21 | 0.31 | 0.27 | 0.36 |

Even at matched count + matched head area (T=10000), the divide method recovers
only ~half the curvature heads within 50 m. **It is a distinct head set, not a
curvature emulator.** Matching curvature at 1 m is the wrong success criterion -
see the field-head comparison.

## 1 m vs Clubb's 53 field heads (ground truth) - the decisive comparison

Recall = fraction of the 53 field heads with a method head within tol.

| method | heads | r@10 | r@30 | r@50 | r@100 | median dist |
|--|--:|--:|--:|--:|--:|--:|
| curvature-DrEICH | 634  | 0.06 | 0.28 | 0.45 | 0.79 | 54 m |
| divides T=10000  | 943  | 0.06 | 0.32 | 0.57 | 0.89 | 46 m |
| divides T=5000   | 1897 | 0.28 | 0.62 | 0.91 | 1.00 | 22 m |

The divide method matches the **field** heads at least as well as curvature-DrEICH
at comparable density (T=10000: 46 m vs 54 m), and better at finer T.
**CAVEAT (density confound):** T=5000's 22 m is partly inflated by head density -
1897 heads is 3x curvature's and 36x Clubb's, so more chances to land near a sparse
53-head set. The honest, density-controlled read is T=10000 ~ curvature, marginally
better. Neither method tightly matches field heads (median 22-54 m), consistent
with the known DrEICH-vs-Clubb offset (~150 m / ~0.5x area; see channelheads
DESIGN.md). Precision vs Clubb is not reported - Clubb mapped only part of the
basin, so most method heads lie in unmapped area.

### Density-controlled (within Clubb's survey hull) - the fair comparison

Basin-wide Clubb recall is confounded by head density. Restricting to the convex
hull of the 53 field heads (0.70 km^2, +50 m buffer) and reporting precision too
(`divides_clubb_extent.py`) removes the confound -- a too-fine T now shows up as
low precision, not free recall.

| method (in hull) | heads | recall | precision | med dist | (tol 30 m) |
|--|--:|--:|--:|--:|--|
| curvature-DrEICH | 30 | 0.28 | 0.43 | 57 m | |
| divides T=10000  | 25 | 0.32 | 0.52 | 46 m | |
| divides T=8000   | 33 | 0.45 | 0.52 | 34 m | |
| divides T=5000   | 45 | 0.62 | 0.51 | 22 m | |
| divides T=3000   | 68 | 0.79 | 0.50 | 18 m | |

**At matched head count (~30 vs 33), divides->chi-z beats curvature-DrEICH on every
metric:** recall 0.45 vs 0.28, precision 0.52 vs 0.43, median 34 vs 57 m. Divide
**precision is flat ~0.51 across all T** -- the finer-T recall gains are NOT a pure
density artifact; the method is ~half-accurate per head at any scale (vs curvature
0.43), and finer T simply finds more of the field heads at the same per-head
quality. The ~0.51 precision ceiling at 30 m (heads with no field head within 30 m)
reflects the known algorithm-vs-field offset (~0.5x area / ~150 m; DESIGN.md) plus
possibly real heads Clubb did not map; at 50 m tol it rises. This is the headline
1 m result: **on the field-head yardstick, density-controlled, the divide method is
better than curvature-DrEICH** -- and (below) far more resolution-robust.

## The other result that matters: resolution robustness (`divides_resolution.py`)

Fixed physical valley scale T = 10000 m^2; mean-aggregate 1 m -> coarser; full
divides->chi-z pipeline at each. Self-consistency = match to the 1 m divide heads
at tol = max(50 m, 2*res); curvR = match of the 634 1 m curvature heads to the
divide heads computed at that res.

| res (m) | heads | self recall | self prec | curvR | curvP |
|--:|--:|--:|--:|--:|--:|
| 1  | 943 | 1.00 | 1.00 | 0.51 | 0.31 |
| 2  | 936 | 0.85 | 0.86 | 0.51 | 0.32 |
| 3  | 903 | 0.78 | 0.81 | 0.49 | 0.31 |
| 5  | 891 | 0.70 | 0.75 | 0.37 | 0.25 |
| 8  | 875 | 0.61 | 0.66 | 0.52 | 0.35 |
| 10 | 854 | 0.67 | 0.73 | 0.56 | 0.39 |
| 12 | 831 | 0.66 | 0.74 | 0.56 | 0.39 |

**Divides->chi-z heads survive coarsening: self-recall 0.66 at 12 m vs curvature's
~0.01.** The method keeps producing field-plausible heads at 10-12 m, exactly
where curvature-DrEICH yields zero (std) or a handful of unreliable heads (adapt2).

**Two honest qualifications.**
1. Head *count* near-constancy (831-943 across 1-12 m) is partly **by
   construction** - T tracks a fixed physical area, so the number of ~10000 m^2
   valleys is ~fixed by the terrain. The meaningful robustness measure is the
   *self-consistency of position* (0.66 @12 m), not the count.
2. The 5 m row dips (0.70/0.37) below its neighbours - likely an aggregation
   grid-alignment artifact (extents snap to non-integer under `r.resamp` /
   block-mean; noted in `dreich_resolution_findings.md`). Not chased yet.

## Status / open items

- **Done:** premise (ridge survival); the valley-scale-T realization; 1 m vs
  curvature + vs Clubb; resolution robustness 1-12 m.
- **Design B (chi-z self-gating) -- BUILT, and it does not work** (`divides_designB.py`).
  Over-generate valleys at a fine T (4841 candidates @T=2000), keep those whose
  chi-z split score >= tau. The score is **saturated**: min 1.59, p25 1.96, median
  1.98, max 2.00 -- so every tau in [0, 0.99] keeps all 4841. The chi-z split
  ALWAYS finds a near-perfect channel-linear sub-segment (R2_chan~1) and an
  autocorrelated hillslope top (DW~0), for spurious candidates as readily as real
  ones, because chi linearises channels so well that any long lower segment is
  ~linear. **Verdict: chi-z fit quality cannot self-gate valleys.** Valley
  SELECTION must come from the divide scale T (Design A); chi-z only LOCALISES the
  head within a valley. Clean division of labour -- divides select, chi-z places.
- **Tuning T against Clubb is confounded by head density.** Clubb recall rises
  monotonically as T shrinks (T=5000: 1897 heads, r@50 0.91, 22 m; T=10000: 943,
  0.57, 46 m) purely because more heads cover a sparse 53-head target. Clubb
  mapped only PART of the basin, so a fair T needs Clubb's survey extent (compare
  counts/precision within it), not basin-wide recall. Owed: clip to Clubb's mapped
  sub-area before setting T.
- **Spatial adaptivity:** T is a single global scale. A spatially varying T (or
  the divide-shape selecting valleys locally) is the natural next step.
- **Second site + the 12 m TanDEM-X target:** MBR is one basin; the north-star is
  12 m TanDEM-X. Validate the robustness on real coarse data, not only aggregated
  1 m.
- **Performance:** the chi-z head loop is pure-Python O(n^2) per profile (~66 s for
  943 valleys at 1 m). Vectorize before productionizing.
- **Productionization:** none yet in `rivernetworkx`/the GRASS module - all dev/.
