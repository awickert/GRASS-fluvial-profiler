# Channel heads from drainage divides — exploratory analyses

Can the drainage-divide network *locate* channel heads (a resolution-robust
replacement for DrEICH's cross-valley curvature)? Reproducible R&D scripts +
figures; throwaway quality, not part of the addon.

Builds on the reproduced Scherler & Schwanghart divide network
(`rivernetworkx.divides`, see [`../ttb_crosscheck/`](../ttb_crosscheck/)) and the
earlier divides→chi-z head validation
([`../divides_heads_findings.md`](../divides_heads_findings.md)). Inputs: cached
1 m Mid Bailey Run FlowInfo `/tmp/divides_proto_fi.pkl` (rebuild with
`../divides_setup.py`), raw DEM `/tmp/dreich_algorithm/bailey_run_dem.flt`, and
Clubb (2014)'s 53 field heads `/tmp/clubb_channel_heads.xlsx`.

## The payoff: resolution robustness vs the FIELD heads (`robustness.py`)

The method `r.fluvial.channelheads method=divides` (divide-defined valleys + chi-z)
vs Clubb's field heads, in-hull (density-controlled), valley scale held physically
constant, `min_segment_length` scaled to physical length:

| res (m) | **divides** recall@50 | precision@50 | curvature recall@50 |
|--:|--:|--:|--:|
| 1  | 0.57 | 0.68 | 0.42 |
| 2  | 0.57 | 0.62 | 0.11 |
| 3  | 0.53 | 0.60 | **0.00** |
| 5  | 0.55 | 0.56 | 0.00 |
| 8  | 0.57 | 0.59 | 0.00 |
| 12 | **0.57** | 0.67 | 0.00 |

**The divide method recovers the field heads essentially as well at 12 m as at
1 m (recall@50 flat ~0.57) — resolution-INVARIANT — while curvature-DrEICH is dead
by 3 m, and divides beat curvature even at 1 m.** Both 1 m numbers reproduce the
earlier independent in-hull computations, so the harness is sound. Self-consistency
(coarse vs 1 m divide heads) holds at 0.92–0.96. (Caveat: this is *aggregated* 1 m
MBR, an intrinsic-resolution test; real 12 m TanDEM-X adds acquisition effects and
is the separate, north-star validation.) `figures/robustness.png`.

This is the constructive capstone to the negative exploration below: divide
*geometry* doesn't locate heads, but divide *structure* (the valleys) makes the
chi-z head finder resolution-robust against ground truth.

## Bottom line (2026-06-30)

**No divide-network feature we tried distinguishes channel heads from same-area
random points — at any resolution.** Curvature (ridge wrap), cross-divide area
contrast, cross-divide Strahler-order contrast, catchment multiplicity, and
hillslope relief all give AUC ≈ 0.5 against the decisive same-area control. The
geomorphic reason: **divide geometry is scale-free** — every nested hollow has a
curved, contrasting apex — while a channel head has an intrinsic *scale*, set by
discharge (∝ area) and slope (the Montgomery–Dietrich threshold). The divides
supply robust valley *structure*; they do **not** encode a head signature beyond
area. So the constructive path stays the validated **v2**: divides delineate the
valleys (robustly, to 12 m), area/discharge sets the scale, chi-z localizes the
transition.

> **A bug worth flagging.** First passes of the clipped scripts gave exciting
> AUCs (order-contrast 0.78, relief 0.88). They were a **head-snapping artifact**:
> the KD-tree used clip-local `col` without adding the clip offset, snapping Clubb
> heads to the wrong cells. Caught by (a) a second implementation
> (`resolution_screen`) disagreeing at 1 m and (b) the same-area control. Fixed;
> every signal then collapsed to null. The same-area control + cross-checking are
> what kept a false breakthrough out.

## Scripts and what each found

| script | question | result |
|--|--|--|
| `apex_duality.py` | is the head ~L_H below its catchment top? | **Sanity check only** (catchment-top-above-head is definitional). L_H ≈ 80 m, CV 0.42 — a useful scale, not a test. |
| `apex_detect.py` | does the divide *above* a head carry a tight-wrap apex signature? | **No** — head apexes are curved (83 %ile) but so are same-area control apexes, identically; "sharp bend" is 17 % of the network. Curvature is scale-free. |
| `cross_divide.py` | is the head across the divide from a *different* catchment (area)? | Weak/null once snapping fixed (AUC 0.43). |
| `feature_screen.py` | which of {area-, order-contrast, multiplicity, relief} discriminates heads? | **None** (AUC 0.43–0.58, corrected). The 0.78/0.88 first pass was the snapping bug. |
| `resolution_screen.py` | do any of these survive coarsening 1→12 m? | All null at all resolutions (the correctly-snapped reference). |
| `debug_disc.py` | isolate the cache-filled vs raw-refilled discrepancy | Showed both give AUC ≈ 0.45 → the earlier 0.78 was the bug, not the DEM. |

## What DOES locate heads (`physics_baseline.py`)

The constructive counterpart. At Clubb's heads, across 1→12 m:

| | 1 m | 12 m |
|--|--:|--:|
| drainage area AUC (heads vs random channel cells) | 0.74 | 0.75 — **resolution-stable** |
| local slope AUC (heads vs same-area cells) | 0.62 | 0.43 — **dies by 12 m** |

Area separates heads and survives coarsening; **slope adds a little at 1 m but is
gone by 12 m** — a derivative, fragile exactly like cross-valley curvature. So the
area–slope (Montgomery–Dietrich) threshold loses half its information at coarse
resolution. This is *why* v2 is built from **area (robust scale) + chi-z (an
integral transition, robust) + divide valley structure (robust)** — only
resolution-stable ingredients, no fragile derivatives. (area control is restricted
to head-scale area, so read its magnitude cautiously; the slope trend is the clean
result.)

## The lesson

Three independent tricks (ridge curvature; cross-divide area/order contrast;
catchment multiplicity) all failed the same control for the same reason: at a
fixed drainage area, a channel head looks like any other point — because **area
(discharge) is the head-forming variable, and divide geometry doesn't add to it.**
This is consistent, not a dead end: it tells us divides are for robust *structure*
and *coarse-resolution survival*, while the head's *position* is a discharge–slope
process threshold that chi-z reads. Pulling a fourth purely-geometric trick is
unlikely to beat area; the productive next step is to lean on area/slope + chi-z
*within* divide-defined valleys (v2), and test that at 12 m.

## Figures

`figures/` — `apex_duality_*` (L_H sanity check), `apex_detect_*` (wrap, null),
`cross_divide_*` (contrast). Note the apex_detect/cross_divide head overlays were
from the buggy snap; the quantitative conclusions come from the corrected
`feature_screen.py` / `resolution_screen.py`.
