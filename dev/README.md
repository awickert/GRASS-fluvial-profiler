# dev/ — diagnostic & prototyping scripts (r.fluvial.hollow R&D)

Throwaway-quality R&D scripts kept around to **diagnose and reproduce** the
hollow / colluvial-to-fluvial-transition development. They are *not* part of the
addon and are *not* portable or run in CI.

**They hardcode local inputs:**
- GRASS locations `TrempealeauHills/chanheads_proto` and `MidBaileyRun`
  (under `~/grassdata/`),
- `/tmp/clubb_channel_heads.xlsx` (Fiona Clubb's 53 Mid Bailey Run heads),

and assume the system-python GRASS setup (see the project memory
`grass-gunittest-recipe` and `../r.fluvial.hollow/DESIGN.md`). Typical run:

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

These produced the numbers and figures recorded in `../r.fluvial.hollow/DESIGN.md`.
They duplicate boilerplate (region setup, the `fit`/broken-stick helper) on
purpose — each is self-contained for quick edit-and-rerun diagnosis.
