# Cross-check: `rivernetworkx.divides` vs TopoToolbox `DIVIDEobj`

Development/testing harness that certifies the `rivernetworkx.divides`
drainage-divide network against the reference implementation, **TopoToolbox
`DIVIDEobj`**, by running both on *identical* flow routing and comparing the
divide edge set and Topo ordering.

> **Reference implementation — TopoToolbox** (GPL), by **Wolfgang Schwanghart**
> and **Dirk Scherler**: <https://github.com/wschwanghart/topotoolbox>
>
> - Schwanghart, W., Scherler, D. (2014). *TopoToolbox 2 — MATLAB-based software
>   for topographic analysis and modeling in Earth surface sciences.* Earth Surf.
>   Dynam. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014
> - Scherler, D., Schwanghart, W. (2020). *Drainage divide networks — Part 1:
>   Identification and ordering in digital elevation models.* Earth Surf. Dynam.
>   8, 245–259. https://doi.org/10.5194/esurf-8-245-2020
>
> This directory **does not vendor TopoToolbox**. You clone it yourself; the
> patches here only make a local clone run under Octave. Both TopoToolbox and
> r.fluvial are GPL.

## Why a controlled (shared-routing) comparison

A divide network is the dual of the drainage basins, which depend entirely on the
flow routing. To isolate the *divide construction* from any routing difference
(the lesson from the LSDTopoTools comparison — routing method dominates), we route
**once**, in the dreich pipeline, and feed those exact D8 receivers to TopoToolbox
as a sparse giver→receiver matrix `M` via `FLOWobj(M, …)`. Both tools then build
divides on the same basins.

## Result (Big Tujunga 30 m, 300×400 window, stream threshold 200 cells)

Opening each basin loop **at its outlet** (the faithful `getdivide`;
`outlet_trim_experiment.py`) reproduces TopoToolbox's divide network:

| construction | recall | precision | Jaccard | (interior, ≥12 from edge) |
|--------------|----:|----:|----:|---|
| all-channel trim (over-removes near-stream) | 0.886 | 0.960 | 0.854 | r 0.920 / p 0.974 |
| **outlet-only trim (= getdivide)** | **0.943** | 0.953 | 0.902 | **r 0.984 / p 0.967 / J 0.952** |

In the interior the outlet-trim matches TopoToolbox at **98.4 % recall / 96.7 %
precision** — a faithful reproduction. The remaining boundary gap is the Octave
shims + `outlets=false` + window clipping.

**Order schemes (`order_robustness.py`).** On the shared edges, every scheme has
median |Δorder| = 0 (the *majority* of edges match exactly), but the *tail*
divergence depends on the scheme's sensitivity to segmentation:

| scheme | my max | TTB max | exact | within ±1 |
|--------|---:|---:|---:|---:|
| Topo     | 25 | 38 | 78 % | 88 % |
| Strahler | 5  | 6  | 82 % | 95 % |
| Shreve   | 123| 266| 72 % | 83 % |

TopoToolbox segments the network more finely (its `divnet` dead-segment merging +
break-insertion). **Strahler is robust** to this (only increments when equal
orders meet → 5 vs 6); **Shreve is fragile** (additive → every extra segment
counts → 2×); Topo is intermediate. So the *geometry* is reproduced and the
*network* is robust; the deep-tail Topo/Shreve *values* are inherently
segmentation-sensitive descriptors. Matching them exactly needs bit-identical
segmentation; Strahler does not.

## Routing is identical (the residual is construction, not routing)

`routing_check.py` confirms TopoToolbox's `flowacc` on `FLOWobj(M)` equals
dreich's `ncontrib` for **100 % of cells** (max|diff| = 0). Since `flowacc` is a
strict function of the receiver tree, the routing is provably identical — the
`M` hand-off works exactly. So the divide residual is purely the *construction
algorithm*, not the flow routing.

## What the residual is (`residual_analysis.py`)

| set | edges | near a stream (≤1 cell) | order-1 |
|-----|------:|---:|---:|
| MISS (TopoToolbox has, I lack) | 1849 | **56 %** (median dist 0) | 66 % |
| EXTRA (I have, TopoToolbox lacks) | 604 | 19 % (median dist 3.6) | 96 % |
| SHARED (the ridge network) | 14342 | 9 % (median dist 7) | — |

The recall gap is dominated by **near-stream** edges: my `stream=` trim *removes
every channel-adjacent edge*, whereas `getdivide` *keeps* the basin perimeter and
opens each loop only at its single **outlet**. So the trim is a simplification of
`getdivide` that over-removes near-stream order-1 leaves; the per-link basins also
add a few deep spurs (the precision gap). The ridge network — what matters for
first-order valleys and channel heads — agrees.

## Reproduce

```bash
# 1. clone the reference (GPL; not vendored here) and make it Octave-runnable
git clone --depth 1 https://github.com/wschwanghart/topotoolbox.git /tmp/topotoolbox
TTB=/tmp/topotoolbox bash apply_octave_patches.sh

# 2. route + build M + compute my divides (writes /tmp/ttb_M.mat, /tmp/ttb_mine_edges.npy)
PYTHONPATH=/path/to/r.fluvial /usr/bin/python3 build_M.py

# 3. run TopoToolbox DIVIDEobj on that routing (writes /tmp/ttb_out_*.txt)
TTB=/tmp/topotoolbox SHIMS=$PWD/octave_shims octave --no-gui --quiet run_divideobj.m

# 4. compare
/usr/bin/python3 compare.py
```

Needs Octave with the `image` and `statistics` packages
(`apt install octave-image octave-statistics`).

## Octave compatibility shims/patches (all faithful or boundary-only)

TopoToolbox targets MATLAB; these make a clone run under Octave without changing
any algorithm:

- **`octave_shims/bwtraceboundary.m`** — Octave lacks `bwtraceboundary`; shimmed
  via `bwboundaries`. Returns the same boundary loop with a different start/
  direction, but the undirected boundary *edge set* (all that getdivide needs) is
  identical.
- **`octave_shims/verLessThan.m`** — returns `false` (assume modern MATLAB; pure
  version gate, no algorithm content).
- **`apply_octave_patches.sh`** edits the clone:
  - `polygon2GRIDobj` → force the non-`polyshape` `poly2mask` branch (Octave has
    no `polyshape`). Only used by `getdivide` for **edge-outlet** basins, so it
    affects window-boundary divides only.
  - `unique(…)` 3rd output (inverse map `ic`) computed via `ismember` — Octave
    does not implement it for `'stable'`/`'rows'`; `ismember` gives the exact
    MATLAB values.
  - `divorder.m`: rename variable/field `do` → `dvo` (`do` is an Octave keyword).
- **`run_divideobj.m`** uses `'outlets', false`: TopoToolbox's edge-outlet
  `getdivide` path indexes a 2nd neighbour that may not exist on a clipped window.

Because the `polyshape`/edge-outlet items touch the boundary, the comparison
reports the **interior** separately as the fair test of the core construction.

## Files

| file | role |
|------|------|
| `build_M.py`          | Python: DEM → dreich fill+route → sparse `M` (for FLOWobj) + my divides |
| `run_divideobj.m`     | Octave: `FLOWobj(M)` → `STREAMobj` → `DIVIDEobj` → `divorder`, export |
| `compare.py`          | Python: edge-set overlap + Topo-order agreement (overall + interior) |
| `routing_check.py`    | Python: prove TopoToolbox `flowacc` == dreich `ncontrib` (routing identity) |
| `residual_analysis.py`| Python: localize the residual by stream-distance + order |
| `outlet_trim_experiment.py` | Python: outlet-only trim (= getdivide) vs all-channel trim, vs the reference |
| `order_robustness.py` | Python: per-scheme order agreement -- Strahler robust, Shreve fragile |
| `apply_octave_patches.sh` | make a TopoToolbox clone Octave-runnable |
| `octave_shims/`       | `bwtraceboundary.m`, `verLessThan.m` |
