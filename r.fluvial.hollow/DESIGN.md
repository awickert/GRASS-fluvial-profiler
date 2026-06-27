# r.fluvial.hollow — design and validation findings

`r.fluvial.hollow` maps the **colluvial-to-fluvial transition**: the downslope
limit of colluvial hollows, where channels become fluvial. This is a *process*
boundary (where fluvial incision begins to dominate hillslope/colluvial
transport), found from the **slope–area break** (the rollover from a flat
hillslope limb to a concave fluvial limb, `S ≈ k·A^−θ`).

It is **not** a morphological field channel head (e.g. Clubb et al. 2014's
DrEICH head) — see *Relation to field channel heads* below.

References: Montgomery & Foufoula-Georgiou (1993); Tarboton et al. (1991/92);
Stock & Dietrich (2003).

## Method

1. **Global slope–area, on *all cells*** (not an extracted network — see finding
   §1): per-cell slope (`r.slope.aspect` gradient) paired with flow accumulation
   over every cell, binned in log–log.
2. **Fit a constrained broken-stick to the binned medians** (`fit_sa_break`): a
   flat hillslope limb joined at a knot `log A*` to a power-law fluvial limb.
   `A*` is the transition drainage area.
3. **Validate** the rollover is real (§2) and **bracket** the extraction
   threshold below `A*` (§3); build a dense network once that over-shoots upslope
   past the hollows.
4. **Place transitions** where drainage area first crosses the (local) `A*`
   along the network (`colluvial_fluvial_transition`) — topology-free, since
   accumulation already integrates upstream area.

The analysis lives in `rivernetworkx` (`fit_sa_break`,
`colluvial_fluvial_transition`, `read_stream_segments`, plus the O(N)
endpoint-hash topology that replaces `v.stream.network`'s O(N²) linking — issue
#13); the module is a thin GRASS wrapper.

## Validation (Mid Bailey Run, OH — 1 m 3DEP; Clubb et al. 2014's 53 mapped heads)

- The slope–area rollover is clear and reproducible: binned-median broken-stick
  gives `R ≈ 0.95` (variance reduction vs a single line), `θ ≈ 0.42`.
- **Relation to field channel heads.** Clubb's 53 heads snap cleanly to the
  extracted channels (all within 10 m; median offset 1.6 m). Their drainage
  areas are **median ≈ 769 m²** (IQR 376–1754). Against the fluvial-transition
  `A*`, the ratio `A_clubb / A*` has **median ≈ 0.52, log-std ≈ 0.48 (×3 at 1σ),
  with 68 % of heads upslope of `A*`** (i.e. in the colluvial zone). So the field
  head and the fluvial transition are **related but distinct**: the field head
  sits roughly half an order of magnitude upslope, in the colluvial channel
  reach. The distance/area between them *is* the extent of that reach — a
  measurable quantity, not an error.

## Key findings

**§1 — The global slope–area must be all-cell (raster), not an extracted
network.** A thresholded `r.stream.extract` network is *all channel*, so its
slope–area data is a single declining fluvial band with **no flat hillslope
limb** — a broken-stick cannot beat one line, giving a *false* "no rollover."
Computing slope per cell over *every* cell restores the hillslope limb. On
Trempealeau this was decisive: extracted-network `R ≈ 0` (false negative) vs
all-cell `R ≈ 0.49` (real, if weak).

**§2 — Validate on binned medians, not raw points.** The rollover lives in the
central tendency; raw-point variance reduction is `R ≈ 0` *even for a genuine
rollover* (MBR global raw `R = 0.006`) because point scatter swamps it. Fit and
test the broken-stick on **binned medians**. Accept a rollover when: binned `R`
clears a **low bar (~0.3)**, the knot is **interior** (edge-margin ≳ 0.1, not
pinned at the data edge — this is what flags spurious edge-pinned fits), and `θ`
is plausible (~0.2–1). Reject only genuinely flat landscapes (`R ≈ 0`).

**§3 — Global `A*` is for bracketing only, not the output.** The `A*` *value* is
method-sensitive (MBR-extracted gave 1454 / 5079 / 803 from raw / binned /
windowed fits — the transition is genuinely gradual). That is fine: the global
step's only job is to **bracket the extraction threshold** — set the minimum
catchment area for `r.stream.extract` safely *below* `A*` so the built network
over-shoots upslope past the hollows. Precision is not needed; erring small is
safe. The actual transitions come from per-flowline detection on the built
network.

**§4 — Slope alone cannot localize individual heads.** At the heads, the
flowline-smoothed slope has log-std ≈ 3.4 — slope (a derivative of elevation) is
noise-limited there. The best slope–area threshold exponent is `p* = 0` (slope
adds only noise to area). So a single area threshold over-predicts head *count*
(20–40×) and mislocates (≥112 m), and per-single-flowline slope–area breaks are
unreliable. **χ–elevation** (an integral, hence smooth) is clean and is the path
to field heads — but Clubb's heads sit ~0.5× area / ~150 m upslope of the gross
χ–z break, so reproducing them needs DrEICH's statistical fluvial-linearity
selection. DrEICH is **deferred** as a complementary method (the planned route:
port the algorithm from LSDTopoTools, let GRASS supply the infrastructure).

## Caveats / edge cases

- **Trempealeau Hills** (the early prototype site) is a small, isolated-bluff
  clip — atypical terrain with a weak, noisy slope–area signal (all-cell
  `R ≈ 0.49`, `A* ≈ 48 cells`). It is a useful *"does the gate over-reject?"*
  sanity check but **not** a calibration target; tune to a proper catchment
  (MBR) and keep the accept bar low enough to admit weak-but-real rollovers.
- The transition is gradual, so head-area is not a sharp single value (>1 order
  of magnitude across a basin).

## Status and open items

- **Done:** building blocks in `rivernetworkx` (`fit_sa_break`,
  `colluvial_fluvial_transition`, `read_stream_segments`, O(N) topology); the v1
  module (over-extracted network → global break `A*` → transition points);
  validated against Clubb's MBR heads.
- **Next:** `slope_area_raster` + an all-cell *bracket-and-validate* helper
  (§1–§3), then **per-flowline local detection** of transitions on the built
  network (windowed/per-block `A*` was too noisy at uniform MBR to resolve real
  spatial variation; per-flowline is the route).
- **Future:** a DrEICH-equivalent for field channel heads (§4).
