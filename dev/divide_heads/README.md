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

## Process-domain workflow (current synthesis, 2026-07)

The head hunt reframed: don't force one break onto Clubb's channel head — **map the
three process domains** (hillslope → colluvial → fluvial; Montgomery & Dietrich) and
their **two transitions + lateral input**, reporting offsets instead of errors. Each
transition wants its *own* cue, and the cues' reliability is set by where the
hillslope sits on the **threshold continuum** (Roering et al. 1999, 2007). The
ordering below is deliberately **non-circular**: the valley head is found *first*
from local curvature (needs no channel, no length), and everything else follows.

```
 z                                                       TRANSITIONS
 ^  crest                                                -----------
 |  /`--.__            HILLSLOPE                          (1) hillslope->colluvial
 | /        `--.__     convex; S<S_c (low E*)                 = VALLEY HEAD  (V)
 |/     (V)       `--._  or planar at S_c (high E*)           curvature: convex->
 |      v             `--.__  COLLUVIAL                       convergent  [local:
 | ......................   `-._ debris-flow; low S-A         channel-free, length-free]
 |                    (C)      `-.__ concavity (theta~0)  (2) colluvial->fluvial
 |                     v          `-.___                      = CHANNEL HEAD (C)
 |                                     `-.__ FLUVIAL          incision onset; located
 |                             (T)        `-.__ S~A^-theta    LAST; report as offset
 |                              v            `-.___       (3) tributary = AREA STEP (T)
 +--------------------------------------------------->        = lateral fluvial input
      |<------- L_H ------->|          downstream dist
      hillslope length = crest -> V  (NOT crest -> C: V is found first, by curvature)


                             DEM
                              |
                     fill + flow routing
                 (receivers, area A, flow distance)
                              |
        +---------------------+---------------------+
        |                     |                     |
    CURVATURE              SLOPE S(A)              AREA A
     (local)                  |                 (monotonic)
        |            [used later, weighted           |
   ridge cells         by regime E*]             area STEPS
   -> C_HT = E/K            (cue 2)          = tributary junctions (T)
        |                                    = LATERAL INPUT; bound
   convex->concave                            the clean reach   (cue 3)
   = VALLEY HEAD (V)  <== ANCHOR: channel-free, length-free   (cue 1)
        |
   L_H = crest -> V
        |
   S_c : assume ~0.8 (tan 39 deg)   OR   calibrate via R*-E* collapse
        |                                (R*=R/(S_c L_H); uses crest->V geometry,
        |                                 same anchor -- NOT crest->channel)
        v
   E* = |C_HT| . L_H / S_c  =  L_H / (S_c/|C_HT|)          [ REGIME DIAL ]
        |
        |   low  E*  -> hillslopes sub-threshold : slope-area theta is informative
        |   high E*  -> slope saturates at S_c   : use curvature + area only
        v
   weight the cues -> locate CHANNEL HEAD (C)   [last; offset where no universal break]
        |
        v
   OUTPUTS:  valley heads (V) | channel heads (C)+offset | colluvial extent | catchments
```

**Why it isn't circular** ($L_H$ = crest→channel would be, since the channel is the
unknown): the transport hillslope ends at the **valley head**, not the channel head,
and V comes from local curvature convergence — so $L_H$=crest→V, then $E^*$, then C
*last*, with C never feeding back. $S_c$ is the one non-topographic input (assume, or
calibrate — but calibration via $R^*$ also uses crest→V geometry, not the channel).
In steep terrain slope *saturates* and carries no domain info, which is exactly why
**curvature (valley head) + area-steps (tributaries) are universal** and slope is
weighted down. Status: catchment decomposition (`catchment_decompose.py`,
`plot_catchments.py`) and the continuum test (`hillslope_continuum.py`) validated at
both MBR and Feather; the curvature valley-head detector is the next build.

## The payoff: resolution robustness vs the FIELD heads (`robustness.py`)

`r.fluvial.channelheads method=divides` (divide-defined valleys + chi-z) vs Clubb's
field heads, in-hull (density-controlled), valley scale held physically constant,
`min_segment_length` scaled to physical length, **with the coarse-resolution chi-z
profile fix** (below):

| res (m) | cells/valley | **divides** recall@50 | curvature recall@50 |
|--:|--:|--:|--:|
| 1–20 | 100→5 | **0.57–0.58** (flat) | 0.42 @1 m, dead by 3 m |
| **30** | 3.3 | **0.32** (was 0.04 without the fix) | 0.00 |
| 40 | 2.5 | 0.06 | 0.00 |
| 50–90 | ≤2 | 0.00 | 0.00 |

**The divide method is resolution-INVARIANT to 20 m (recall@50 ~0.57), holds
partially to 30 m (~0.32) with the fix, and is gone by ~40–50 m** &mdash;
curvature-DrEICH is dead by 3 m and divides beat it even at 1 m. The 1 m numbers
reproduce the earlier independent computations (943 heads @T=10000, unchanged by
the fix &mdash; no regression).

**Failure mechanism & the fix (`robustness.py`, `srtm30.py`).** The 30 m collapse
was **chi-z profile starvation, not loss of valley structure**: at 30 m a
10000 m&sup2; valley is ~3 cells across, so the hilltop&rarr;source profile is
~3&ndash;5 nodes &mdash; too few for the split. **Nyquist-like criterion: ~5 cells
across a first-order valley** for full recall. The fix extends the chi-z profile
downstream &mdash; but it works only at the **exact split minimum**
(`2*min_segment_length+1` nodes): diagnosed by head drainage area, one node *more*
drags the head **into the downstream trunk** across a ksn break (head area jumps
6&rarr;33 at 30 m, 12&rarr;141 at 20 m). So `min_profile_nodes` auto-defaults to
that minimum &mdash; a no-op at fine resolution, recovering 30 m to recall 0.32.

**Honest caveat on the "valleys found" count.** In `figures/robustness.png` the
valley count *rises* to ~229 at 70 m &mdash; that is an **artifact**, not persistence:
the valley-scale threshold `T = T_m&sup2;/res&sup2;` floors at 2 cells past ~50 m, so
`get_sources` counts grid-scale noise, not physical first-order valleys. The clean
gauge is **cells/valley = &radic;T_m&sup2; / res**: heads need ~5 (res &le; 20 m full,
~3 at 30 m partial); the physical first-order valley (~100 m) is itself resolved only
to ~2 cells (res &le; ~50 m). So the divides outlive the chi-z **modestly** (heads
~30–40 m, structure ~50 m), **not** to 90 m.

(Caveat: *aggregated* 1 m MBR, an intrinsic-resolution test; real 12/30 m data adds
acquisition effects and is the separate, north-star validation.)

This is the constructive capstone to the negative exploration below: divide
*geometry* doesn't locate heads, but divide *structure* (the valleys) makes the
chi-z head finder resolution-robust to ~30 m against ground truth.

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
