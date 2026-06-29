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

| metric | all edges | interior (≥12 from edge) |
|--------|----:|----:|
| precision (my edges TopoToolbox also has) | 0.960 | **0.974** |
| recall (TopoToolbox edges I reproduce)    | 0.886 | **0.920** |
| Jaccard                                   | 0.854 | 0.899 |
| first-order (order-1) divides matched     | **94 %** (6815/7247) | — |
| Topo order on shared edges                | Pearson r = 0.745 (max 38 vs 25) | — |

The edge-trim construction in `rivernetworkx.divides` reproduces TopoToolbox's
divide network at ~97 % precision / ~92 % interior recall, with first-order
divides agreeing 94 %. The Topo-order *tail* diverges (it is +1 per junction, so
small segmentation differences compound); low orders agree.

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
| `apply_octave_patches.sh` | make a TopoToolbox clone Octave-runnable |
| `octave_shims/`       | `bwtraceboundary.m`, `verLessThan.m` |
