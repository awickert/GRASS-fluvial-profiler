# DrEICH channel heads: resolution sensitivity (Ohio / Mid Bailey Run)

Test case: USGS 1 m LiDAR of Mid Bailey Run, OH (the DrEICH validation DEM),
mean-aggregated (`r.resamp.stats average`) to 1, 2, 3, 4, 5, 8, 10, 12, 15, 20,
30 m. DrEICH run with the internal FastScape-style D8 routing (the canonical
path). Curvature polyfit window held at 7-cell support (`window_radius = 7*res`).
**1 m result (637 heads) is taken as truth.** Recall = fraction of 1 m heads with
a detected head within `tol = max(50 m, 2*res)`; precision = fraction of detected
heads within `tol` of a 1 m head.

Scripts: `dreich_resolution_sweep.sh` (sweep), `dreich_curv_calibrate.py`
(adaptive threshold), `dreich_resolution_compare.py` (recall/precision). Outputs
saved in GRASS mapset `MidBaileyRun/dreich_restest` as
`heads_<scheme>_<N>m` / `net_<scheme>_<N>m` / `rnet_<scheme>_<N>m`.

## Three parameter schemes

- **std** — standard 1 m-tuned DrEICH: `tan_curv_threshold=0.1`,
  `n_connecting_nodes=10`, `threshold=100`, `min_segment_length=10`.
- **adapt** — std + a **resolution-adaptive curvature threshold**: hold the
  fraction of cells flagged valley-like constant (`0.1` flags the top 0.27 % of
  curvature at 1 m; each coarser res uses the matching percentile of its own
  curvature). Connectivity unchanged.
- **adapt2** — adapt threshold **+ relaxed connectivity**:
  `n_connecting_nodes=3`, `min_segment_length=5`, `threshold=50`.

## Note: curvature is in spatial units, which is why "adaptive" is a choice

The tangential curvature (`rivernetworkx.dreich.tangential_curvature`, faithful to
LSDTopoTools) is computed in **spatial units, 1/length (1/m)** — the polyfit
coordinates carry the cell size (`xk = (i-kr)*cellsize`), so the second
derivatives are 1/m. Therefore `tan_curv_threshold = 0.1` is a **fixed physical
threshold**, consistent across resolutions.

Two consequences for reading this study:
- The **std collapse to 0 heads beyond 1 m is a genuine physical signal loss**,
  not a unit artifact: an averaged/coarser surface really has less curvature at
  the 0.1 /m level (convergent fingertips smooth below the bar; the growing
  `7*res` m window smooths further).
- The **adapt / adapt2 threshold deliberately abandons that fixed physical bar**
  and instead flags the top-X% most-convergent cells at each resolution (a
  *relative* convergence criterion). So the coarse-resolution heads it recovers
  rest on genuinely weaker curvature than the 1 m ones — not the same physical
  feature class. The recovered counts (b, below) are bought by lowering the
  physical threshold.

So metric (b) "unavoidable loss" has two honest readings: under a fixed physical
threshold the loss is near-total past 1 m; under the relative criterion it is the
0.69 -> 0.06 curve below. (Units verified directly from the code; the "fixed-
threshold loss is near-total" reading is inferred from the std sweep, not a
separately controlled run.)

## Physical basis of the threshold scaling (and Grieve et al., 2016)

The calibrated thresholds follow a near-perfect power law
`tan_curv_threshold ≈ 0.118·res^(−1.12)` (R²_log = 0.993). For self-affine
topography, curvature measured at scale L scales as κ ~ L^(H−2) (H = Hurst
exponent); with the window = 7·res, that predicts threshold ~ res^(H−2), i.e.
H = β + 2 = 0.88. Measured **independently** from the 1 m DEM (2nd-order
structure function, isotropic), **H ≈ 0.93** — agreeing with 0.88 to ~0.05. So
the threshold reduction is governed by the terrain's self-affine scaling, not
fitting. (A power-spectrum estimate gave an unphysical H = 1.5 — Hanning/binning
bias; the structure function is the robust estimator and is trusted here.)
Code: `dev/dem_hurst.py`.

**Grieve, Mudd, Hurst, Milodowski & Furbish (2016), *Earth Surf. Dynam.* 4,
627–653** ("How does grid-resolution modulate the topographic expression of
geomorphic processes?") is essentially this experiment by the DrEICH group, and
corroborates it: identical tangential-curvature definition (their Eq. 4); the
curvature **window is set by the hillslope length scale, not grid resolution**,
and hillslope-scale curvature is **only reliable at resolutions ≤ ~10 m** (so our
12–30 m points are beyond the trustworthy range); their threshold is set from the
**curvature distribution** (deviation from normal on a Q–Q plot — same family as
our fraction-matching); and they find DrEICH heads degrade fast and "bear little
relation" to the 1 m heads at coarse resolution. Their **reliability** (Eq. 5)
and **sensitivity** (Eq. 6), with a 30 m coincidence radius, are exactly our
precision and recall.

## Head counts and correctness vs the 1 m truth

| res (m) | std | adapt | adapt2 | adapt2 recall | adapt2 precision | adapt2 med det→truth |
|--:|--:|--:|--:|--:|--:|--:|
| 1  | 637 | –   | 637 | 1.00 | 1.00 | 0 m  |
| 2  | 0   | 235 | 438 | 0.49 | 0.68 | 30 m |
| 3  | 0   | 107 | 262 | 0.28 | 0.61 | 39 m |
| 4  | 0   | 41  | 173 | 0.14 | 0.47 | 53 m |
| 5  | 0   | 26  | 119 | 0.12 | 0.53 | 46 m |
| 8  | 0   | 3   | 58  | 0.05 | 0.40 | 72 m |
| 10 | 0   | 0   | 38  | 0.03 | 0.34 | 67 m |
| 12 | 0   | 0   | 23  | 0.01 | 0.35 | 75 m |
| 15 | 0   | 0   | 10  | 0.01 | 0.50 | 48 m |
| 20 | 0   | 0   | 5   | 0.00 | 0.60 | 47 m |
| 30 | 0   | 0   | 4   | 0.00 | 0.00 | 98 m |

## Findings

1. **Standard 1 m parameters do not transfer at all** — heads only at 1 m, zero
   at every coarser step. Coarse-resolution DrEICH *requires* adapted parameters.
2. **Two levers, in order of importance.** (a) The **curvature threshold** must
   scale down ~20× from 1 m to 10 m — averaged DEMs have far smaller curvature
   magnitudes. (b) **Connectivity** (`n_connecting_nodes`): requiring 10 connected
   high-curvature cells is too strict on coarse grids; 3 recovers heads through
   30 m. (Window choice was secondary; `kr=1`, a 3×3 fit, is well-posed.)
3. **With both adapted, heads are recoverable at every resolution to 30 m** —
   relevant for global 12/30 m data.
4. **But recovery of the 1 m heads degrades steeply regardless of tuning**:
   recall 0.49 at 2 m, ~0.05 at 8 m, ~0 by 20–30 m. Coarsening *loses the fine
   headwater channels* — a physical resolution limit, separate from parameter
   tuning. The heads that survive are the larger, more-convergent channels, and
   they sit reasonably near real ones (precision ~0.34–0.68, median det→truth
   30–75 m) — **until 30 m, where precision collapses to 0** (the 4 heads do not
   correspond to any 1 m head).

## Caveats (read before trusting these numbers)

- **"1 m = truth" is an assumption.** Recall/precision measure recovery of the
  *1 m-detected* heads, not correctness against field-mapped channel heads. A
  coarse head is not "wrong" for failing to match a 1 m head — it detects what is
  resolvable at its scale. Whether the 1 m heads themselves are field-correct is a
  separate question (Clubb's mapped heads).
- **The adapt / adapt2 schemes are constructed here** (fraction-matched threshold,
  connectivity relaxation), not a published DrEICH-at-coarse-resolution standard.
  Defensible, but unproven against ground truth.
- Resolutions snap to non-integer extents under `g.region res=N -a`; edges differ
  by < one coarse cell between resolutions.
