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
`tan_curv_threshold ≈ 0.118·res^(−1.12)` (R²_log = 0.993). The terrain is
self-affine — Hurst exponent from the 1 m DEM (2nd-order structure function,
isotropic) is **H ≈ 0.93** (`dev/dem_hurst.py`; a power-spectrum estimate gave an
unphysical 1.5, a Hanning/binning artifact, so the structure function is trusted).

**CORRECTION (do not over-read the H-match).** I initially inferred H = β+2 = 0.88
from the threshold slope by assuming threshold ~ res^(H−2), i.e. that curvature
scales as the naive self-affine 2nd derivative κ ~ L^(H−2). That assumption is
**wrong here**: measured directly (`dev/dem_hillslope_scale.py`), the polyfit
curvature spread scales with *window* as only **W^(−0.6)**, not W^(H−2) = W^(−1.07).
The least-squares 2nd-order fit regularizes the curvature, so it is not a clean
point 2nd derivative. Consequently the threshold-vs-resolution slope (−1.12) is a
**combination** of window-growth (~−0.6, since window = 7·res) and DEM-averaging
(~−0.5), and its near-match to (H−2) = −1.07 is **partly coincidental**. What
stands: the terrain is self-affine (H≈0.93); the threshold follows a clean
empirical power law; curvature genuinely shrinks with coarsening (also Grieve et
al. Fig 4/5). What does NOT stand: "the threshold scaling directly measures the
Hurst exponent." The constant-window sweep disentangles the window vs averaging
contributions.

**Curvature measurement scale** (`dev/dem_hillslope_scale.py`): curvature spread
breaks marginally at ~4–6 m and flattens beyond (the break is weak because the
terrain is self-affine, so curvature-std vs window is near a featureless power law
in 2–30 m). This ~5–7 m is the curvature radius of ridge/valley-head features and
the appropriate curvature window — used as W_h = 7 m (the validated value).

**NOT the hillslope length.** The hillslope length L_H (divide to channel head)
for Bailey Run is ~100 m (from ~637 heads / ~36 km² → valleys ~235 m apart). The
DrEICH resolution limit is set by the few-metre **valley-head convergence
feature** (≈ the curvature window), NOT by L_H — you can resolve a 100 m hillslope
at 8 m, but the fine convergent fingertip that marks the channel head is smoothed
away. (Earlier text loosely called the ~5–7 m curvature scale the "hillslope
scale", following Grieve et al.'s loose usage; corrected here.)

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

## Constant hillslope-scale window + fixed threshold (disentangles the mechanism)

`cwfix`: window held at W_h = 7 m where valid (grows to ~1.5·res only past res≈5 m,
the validity limit), tan_curv = 0.1 fixed, standard other knobs. Heads:
1m=640, 2m=345, 3m=104, 4m=12, 5m=2, ≥8m=0. Compared to `std` (window 7·res,
same fixed 0.1) which gave 0 beyond 1 m, this isolates the cause:

- **Window-growth is the catastrophic lever.** At 2 m the ONLY difference between
  `std` (0 heads) and `cwfix` (345 heads) is the window (14 m vs 7 m). Oversmoothing
  by the growing window, not the averaging, is what zeroed the fixed-threshold runs.
- **DEM-averaging drives the count loss, even at a fixed window.** With window held
  at 7 m, heads still fall 640→345→104→12 over 1→4 m (count_kept 1.0→0.54→0.16→0.02).
  That steep loss is averaging removing the fine channels — unavoidable (metric b).
- **A fixed threshold is resolution-robust only while the window can stay constant**
  — here res ≤ ~4–5 m, set by the ~7 m curvature window (the valley-head feature
  scale, NOT the ~100 m hillslope length). Past that the window must grow (singular
  otherwise), and fixed 0.1 collapses to 0 by 8 m.

So: hold the window at the curvature/valley-head scale and a *fixed* threshold works
up to ~that resolution; beyond it, only an adaptive threshold recovers anything
(the `adapt2` row), and those heads are few and increasingly unreliable. Goodness
of fit of the `cwfix` detections (median distance to a 1 m head): 26 m @2m, 54 m
@3m; sensitivity/reliability @30 m: 0.31/0.57 @2m, 0.04/0.23 @3m.

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
