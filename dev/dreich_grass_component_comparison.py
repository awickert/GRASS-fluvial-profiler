"""Compare the faithful LSDTopoTools DrEICH pipeline (method=lsdtt) against
GRASS-native drop-in replacements for its components, on the Bailey window.

Question (Andy): if we strip pieces out of the LSDTT port and use what GRASS
already provides, how similar are the channel heads -- and why or why not?

FINDINGS (2026-06-28):

1. TANGENTIAL CURVATURE -- GRASS r.param.scale does NOT reproduce LSDTT's.
   On the LSDTT filled DEM (7107x6266), correlation with the C++ tangential
   curvature (bailey_run_dem_tan_curv):
       my polyfit port : corr  1.00000   (faithful)
       r.param.scale crosc : corr -0.016  std ~3e6   (uncorrelated, blown-up scale)
       r.param.scale planc : corr  0.000  std ~1e3   (uncorrelated)
       r.param.scale longc : corr -0.510  std ~0.017 (moderate, comparable scale)
   r.param.scale uses a different surface fit (window weighting, normalisation,
   and curvature definitions) than LSDTT's circular-masked 2nd-order polyfit, and
   the crosc/planc denominators (~1/slope) explode on near-planar windows.
   Swapping the BEST analog (-longc, threshold re-tuned to the fitted scale) into
   the pipeline -- changing ONLY curvature, keeping the faithful routing -- gives:
       LSDTT-faithful curvature : 635 heads, 630/634 overlap with C++  (99.4%)
       GRASS -longc curvature   : 694 heads, 278/634 overlap with C++  (43.8%)
   i.e. ~56% of channel heads move. The faithful polyfit curvature is necessary.

2. FLOW ROUTING -- GRASS routing diverges from LSDTT's steepest-descent D8.
   Established in prior diagnostics (dev/dreich_routing_*; dreich-port-status):
   port-vs-reference head distance stays ~33-40 m across r.watershed (least-cost
   D8, 40.0 m), r.fill.dir (steepest-descent D8, 41.6 m), and even LSDTT's exact
   filled DEM + steepest-descent D8 (39.8 m). The faithful LSDFlowInfo port
   reproduces sources EXACTLY (112757) where GRASS routing does not.

CONCLUSION: there is no promising GRASS-component drop-in to promote to an
alternative method= option; faithful reproduction of DrEICH requires LSDTT's own
curvature and routing, which method=lsdtt ports. A *separate*, intentionally
GRASS-native "DrEICH-like" head finder (its own calibration, not claiming to
reproduce LSDTT) remains possible future work, but it is a different product.

Reproduce the curvature comparison (needs the /tmp reference rasters + GRASS):
    # 1. GRASS: import the filled DEM, run r.param.scale, export crosc/planc/longc
    #    (see the g.region / r.in.bin / r.param.scale calls in the session notes)
    # 2. python this script  (expects /tmp/bailey_{crosc,planc,longc}_grass.flt)
"""
import numpy as np
from rivernetworkx import dreich as D

NR, NC = 7107, 6266
ND = -9999.0
NDF = np.float32(ND)
ALG = '/tmp/dreich_algorithm'


def main():
    cpp = np.fromfile(f'{ALG}/bailey_run_dem_tan_curv.flt', dtype='<f4').reshape(NR, NC)
    mine = np.load('/tmp/dreich_tancurv_mine.npy')
    print('--- curvature field correlation vs C++ tangential curvature ---')
    for tag, path, arr in [('my polyfit port', None, mine),
                           ('r.param.scale crosc', '/tmp/bailey_crosc_grass.flt', None),
                           ('r.param.scale planc', '/tmp/bailey_planc_grass.flt', None),
                           ('r.param.scale longc', '/tmp/bailey_longc_grass.flt', None)]:
        g = arr if arr is not None else np.fromfile(path, dtype='<f4').reshape(NR, NC)
        both = (cpp != NDF) & (g != NDF)
        c = cpp[both].astype(np.float64); gg = g[both].astype(np.float64)
        a = np.polyfit(c, gg, 1)[0]
        print('  %-22s corr=%+.4f  scale=%.4g  std=%.4g'
              % (tag, np.corrcoef(c, gg)[0, 1], a, gg.std()))

    print('--- channel heads: swap ONLY curvature, faithful routing ---')
    filled = np.fromfile(f'{ALG}/bailey_run_dem_fill.flt', dtype='<f4').reshape(NR, NC)
    longc = np.fromfile('/tmp/bailey_longc_grass.flt', dtype='<f4').reshape(NR, NC)
    gcurv = np.where(longc != NDF, -longc, NDF).astype(np.float32)
    ref = np.fromfile(f'{ALG}/bailey_run_dem_CH.flt', dtype='<f4').reshape(NR, NC) != NDF
    for tag, curv, thr in [('LSDTT-faithful curvature', None, 0.1),
                           ('GRASS -longc curvature', gcurv, 0.059)]:
        heads = D.extract_channel_heads(filled, nodata=ND, cellsize=1.0,
                                        fill_dem=False, curvature=curv,
                                        tan_curv_threshold=thr)
        m = np.zeros((NR, NC), bool)
        for r, c in heads:
            m[r, c] = True
        print('  %-26s heads=%4d  overlap_with_C++=%3d / 634 (%.1f%%)'
              % (tag, len(heads), int((m & ref).sum()), 100 * (m & ref).sum() / 634))


if __name__ == '__main__':
    main()
