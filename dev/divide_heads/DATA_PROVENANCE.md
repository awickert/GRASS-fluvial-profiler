# Data provenance — channel-head validation (dev/divide_heads)

## Field-mapped channel heads (ground truth)
Clubb et al. (2014), *Objective extraction of channel heads from high-resolution
topographic data*, Water Resources Research, doi:10.1002/2013WR015167.
3-site head coordinates: `Channel_head_coords.xlsx` (104 heads: Mid Bailey Run OH
53, Indian Creek OH 36, Feather River CA 15).

## DEMs

### Feather River, CA (Feather South + Feather North)
Both boxes are subsets of the same OpenTopography lidar dataset (the Clubb 2014
Feather River survey), gridded to 1 m bare-earth via TIN (ground returns; grid
resolution 1 m, max triangle 50 m; GeoTIFF).

> **Dataset Citation:** Oroville, CA: Middle Fork of the Feather River, Sierra
> Nevada Foothills. Distributed by OpenTopography. https://doi.org/10.5069/G9FT8HZ1
> . Accessed 2026-07-01.
> **Use License:** CC BY 4.0

- Feather **South** (10 heads): UTM 10N E 644195–646406, N 4388463–4391321 (2211×2858, relief 299–1139 m).
- Feather **North** (5 heads):  UTM 10N E 648000–652000, N 4395000–4398000 (4000×3000, relief 377–1326 m).
  DEMs: NAD83/UTM 10N (EPSG:26910), 1 m, 0% NoData; field heads nominally UTM 10N
  (≤~1–2 m possible NAD83/WGS84 datum offset, well below the 30/50 m match tolerance).
- The 15 spreadsheet heads split S/N at N ≈ 4,394,000 (a ~6.2 km gap): 10 South, 5 North.

### Mid Bailey Run, OH and Indian Creek, OH
1 m bare-earth lidar DEMs from the Clubb et al. (2014) DrEICH test data
(`/tmp/dreich_algorithm/bailey_run_dem.flt`, `indian_creek_dem.flt`; UTM 17N).
*(Original OpenTopography DOIs to be added if we publish — TODO.)*
