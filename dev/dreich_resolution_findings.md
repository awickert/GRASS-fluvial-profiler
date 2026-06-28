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
