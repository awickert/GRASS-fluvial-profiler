"""Shared 1 m Mid Bailey Run setup for the divides->channel-heads prototype.

Builds the DrEICH FlowInfo ONCE on the 1 m MBR DEM (reusing the already-filled
DEM and tangential-curvature rasters on disk, so we skip the two slowest steps),
runs the standard curvature valley-gate to get the reference heads, and pickles
the FlowInfo + curvature + reference to /tmp so the A/B gate prototypes can
iterate in seconds instead of re-running the ~minutes-long flow setup.

It also INSTRUMENTS the curvature path so we build the divide-based gates against
verified topology semantics (what node does find_valleys tag? what span does the
chi-z profile cover?) rather than a guess.

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_setup.py
"""
import pickle
import time
import numpy as np
from rivernetworkx import dreich as D

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
RES = 1.0
THRESHOLD = 100           # DrEICH source area threshold (cells)
CACHE = '/tmp/divides_proto_fi.pkl'


def log(msg):
    print('[%6.1fs] %s' % (time.time() - T0, msg), flush=True)


T0 = time.time()
log('loading filled DEM + curvature (.flt)')
filled = np.fromfile('%s/bailey_run_dem_fill.flt' % ALG, '<f4').reshape(NR, NC)
tcurv = np.fromfile('%s/bailey_run_dem_tan_curv.flt' % ALG, '<f4').reshape(NR, NC)

log('build_flowinfo (steepest-descent D8 on filled surface)')
fi = D.build_flowinfo(filled, ND, RES)
log('contributing_area')
D.contributing_area(fi)
log('get_sources')
sources = D.get_sources(fi, THRESHOLD)
log('junction_network')
D.junction_network(fi, sources)
log('build_svector + distance_from_outlet')
D.build_svector(fi)
D.distance_from_outlet(fi)
log('%d sources, %d junctions' % (len(sources), len(fi['JunctionVector'])))

log('find_valleys (curvature gate, the step we are replacing)')
valley = D.find_valleys(fi, tcurv, sources, n_connecting_nodes=10,
                        tan_curv_threshold=0.1)
nz = np.where(valley != -9999)[0]
vjuncs = np.unique(valley[nz])
log('curvature gate: %d tagged nodes -> %d distinct valley junctions'
    % (len(nz), len(vjuncs)))

log('channel_heads_from_valleys (chi-z stage; UNCHANGED across A/B)')
_, final = D.channel_heads_from_valleys(fi, filled, valley,
                                        min_segment_length=10, A_0=1000.0,
                                        m_over_n=0.525)
log('curvature-DrEICH reference: %d heads' % len(final))

# ---- instrument the topology semantics the chi-z stage actually consumes ----
row_of, col_of, recv = fi['row_of'], fi['col_of'], fi['recv']
JV, ncontrib = fi['JunctionVector'], fi['ncontrib']
print('\n--- what find_valleys hands the chi-z stage (5 sample valley juncs) ---')
for jn in vjuncs[:5]:
    dn = int(JV[jn])
    hilltop = D._find_farthest_upslope(fi, dn)
    nodeseq, chi, elev = D._build_channel_chi(fi, filled, hilltop, dn,
                                              1000.0, 0.525)
    head = D._calculate_channel_head(nodeseq, chi, elev, 10)
    print('jn=%d  dn area=%d cells  hilltop area=%d  profile_len=%d  '
          'head_idx=%d  head area=%d'
          % (jn, ncontrib[dn], ncontrib[hilltop], len(nodeseq), head,
             ncontrib[nodeseq[head]] if 0 <= head < len(nodeseq) else -1))

# profile-length + junction-area distributions across ALL curvature valleys
plens, dn_areas = [], []
for jn in vjuncs:
    dn = int(JV[jn])
    hilltop = D._find_farthest_upslope(fi, dn)
    nodeseq, _, _ = D._build_channel_chi(fi, filled, hilltop, dn, 1000.0, 0.525)
    plens.append(len(nodeseq))
    dn_areas.append(int(ncontrib[dn]))
plens = np.array(plens); dn_areas = np.array(dn_areas)
print('\nprofile length (nodes): med=%d  p10=%d  p90=%d  min=%d  max=%d'
      % (np.median(plens), np.percentile(plens, 10), np.percentile(plens, 90),
         plens.min(), plens.max()))
print('valley-junction area (cells): med=%d  p10=%d  p90=%d'
      % (np.median(dn_areas), np.percentile(dn_areas, 10),
         np.percentile(dn_areas, 90)))

log('pickling FlowInfo + curvature + reference heads -> %s' % CACHE)
with open(CACHE, 'wb') as f:
    pickle.dump(dict(fi=fi, tcurv=tcurv, filled=filled,
                     sources=np.asarray(sources),
                     ref_heads=np.asarray(final, dtype=np.int64),
                     curv_valley=valley), f, protocol=4)
log('done')
