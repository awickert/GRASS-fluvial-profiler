# dev/ — diagnostic & prototyping scripts (r.fluvial.channelheads R&D)

Throwaway-quality R&D scripts kept around to **diagnose and reproduce** the
hollow / colluvial-to-fluvial-transition development. They are *not* part of the
addon and are *not* portable or run in CI.

**They hardcode local inputs:**
- GRASS locations `TrempealeauHills/chanheads_proto` and `MidBaileyRun`
  (under `~/grassdata/`),
- `/tmp/clubb_channel_heads.xlsx` (Fiona Clubb's 53 Mid Bailey Run heads),

and assume the system-python GRASS setup (see the project memory
`grass-gunittest-recipe` and `../r.fluvial.channelheads/DESIGN.md`). Typical run:

```
PROJ_DATA=/usr/share/proj PROJ_LIB=/usr/share/proj PYTHONPATH=. \
  GRASS_PYTHON=/usr/bin/python3 \
  grass ~/grassdata/TrempealeauHills/chanheads_proto --exec \
  /usr/bin/python3 dev/<script>.py
```

Some need the `r.stream.basins` addon (`g.extension extension=r.stream.basins`).

## Scripts

| script | what it does |
| --- | --- |
| `accumulate_prune_pass.py` | The accumulate-and-prune first pass: grow each branch's all-cell contributing-area slope–area until a rollover fires, mark the transition, prune the subwatershed's tributaries. Reports hollows, % pruned, local A\*. |
| `fluvial_over_hillshade.py` | Build the fluvial network pruned at hollows (keep channel cells with A ≥ local A\*) and plot it over a hillshade → PNG. |
| `false_negative_map.py` | Recall check: flag unresolved branches whose local (A, S) matches a *nearby* hollow's signature, then concentrate to branches and **re-fit** to separate real near-misses from coincidental matches. |
| `nearmiss_characterize.py` | Compare near-miss branches (R≥0.3 but failing the θ gate) to confirmed hollows on θ, A\*, transition slope, drainage area, elevation. |
| `floodplain_check.py` | Quantify floodplain flat/weird-slope contamination and the data-driven elevation-cut + slope-floor fix. |
| `subwatershed_scale_test.py` | All-cell subwatershed detection pass-rate vs drainage-area scale (`r.water.outlet` basins). |
| `mbr_slope_area.py` | Mid Bailey Run slope–area + broken-stick break, with Clubb's heads overlaid. |
| `mbr_vs_clubb.py` | Relationship between A\* and Clubb's mapped head areas (the ~0.5× / colluvial-zone result). |

### Per-flowline rollover (current direction — replaces the subwatershed-pooling pass)

These trace **down from every channel head** along the drainage-direction
raster (GRASS routing — *not* the old endpoint-hash topology) and look for the
slope–area rollover on each flowline. This sidesteps the pooling/topology
fragility diagnosed below.

| script | what it does |
| --- | --- |
| `rollover_peak_down_only.py` | First cut: walk down from each head to the slope–area peak (running-max + sustained-decline). Over-detects (≈32 k, some peaks land in the trunk). Superseded. |
| `rollover_peak_detector.py` | Robust version: per-head **up+down tent fit** (rising colluvial limb *and* falling fluvial limb, variance-reduction gate), then dedup to unique rollover cells. Trempealeau → 2855 transitions, median area ≈ 820 m², θ≈0.47, both bowls populated. **Open issue:** the tent can lock onto the *first* rise-then-fall rather than the persistent hillslope→channel regime change — under revision. |
| `fluvial_anchor_from_mouth.py` | Earlier failed attempt: fit `S∝A^-θ` to the developed reach and walk *up* to the departure. Misfires because a head's flowline runs through a confluence into the trunk, so the fit anchors on the trunk. Kept as the informative negative result. |

### Diagnostics for the "only some subcatchments have hollows" bug

The bug was **not** terrain (west/east are nearly identical) — it was the
subwatershed-pooling plumbing (hand-rolled topology + off-map-contaminated
accumulation). These scripts establish that, and motivated the per-flowline
rewrite above.

| script | what it does |
| --- | --- |
| `subwatershed_pooling_bug.py` | Per-candidate accumulate-and-prune outcome map (passed/failed/pruned/data-starved) split west vs east; forced fits of the largest west subwatersheds. |
| `ancestor_coverage_check.py` | Smoking gun: compares each segment's accumulation (`maxA`) to the cells its ancestors' basins actually cover — shows the pools don't match true drainage (off-map inflow + topology gaps). |
| `west_east_terrain_compare.py` | Matched-elevation-band slopes, west vs east — shows the terrain is the same (disproves the terrain explanation). |
| `west_east_slope_area_curves.py` | Binned slope–area curves for representative west vs east subwatersheds. |
| `slope_floor_sweep.py` | West/east pass counts vs the slope floor — disproves the slope-floor hypothesis. |

These produced the numbers and figures recorded in `../r.fluvial.channelheads/DESIGN.md`.
They duplicate boilerplate (region setup, the `fit`/broken-stick helper) on
purpose — each is self-contained for quick edit-and-rerun diagnosis.
