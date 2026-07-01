"""Valley -> nested sub-catchment decomposition keyed on tributary AREA-STEPS
(the universal signal: crisp at both MBR and Feather). For each first-order valley
(divide-bounded, above a T-source) walk the trunk DOWNSTREAM and record every
tributary junction as an area step:  at junction node j with trunk-upstream donor
u,  tributary cells = ncon[j] - ncon[u] - 1  (ncon counts self, source=1).

The FIRST significant tributary junction below the source bounds the clean
single-thread valley reach (downstream of it, lateral input contaminates a single
profile -- this is the profile bound for rollover/incision detection).

Step 1 here: just find + print the junctions for a few MBR sources and confirm they
match the area-distance steps we saw (~1e4->1e5 near 200 m). Run:
  PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 catchment_decompose.py <mbr|feather_s> [T]
"""
import sys
import numpy as np
from rivernetworkx import dreich as D
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
T = float(sys.argv[2]) if len(sys.argv) > 2 else 5000.0
TRIB_FRAC = 0.20                                          # tributary counts if its area >= this * trunk area at the junction
cfg = SITES[key]
z, west, north, res, nd = read_flt(cfg['flt']); z = np.where(z == nd, np.nan, z)
field = load_field_heads(cfg['key'], cfg.get('nfilter'))
nr, nc = z.shape
cr = ((north - field[:, 1]) / res).astype(int); cc = ((field[:, 0] - west) / res).astype(int)
ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < nc)
R0, R1 = max(0, cr[ins].min() - 400), min(nr, cr[ins].max() + 400)
C0, C1 = max(0, cc[ins].min() - 400), min(nc, cc[ins].max() + 400)
zc = z[R0:R1, C0:C1]
dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32)
print('routing %s (T=%.0f)...' % (cfg['name'], T), flush=True)
filled = D.fill(dem, nd, cellsize=res); fi = D.build_flowinfo(filled, nd, res)
D.contributing_area(fi); D.distance_from_outlet(fi)
ncon, recv, fd = fi['ncontrib'], fi['recv'], fi['fd']
pxa = res * res
Tc = T / pxa                                              # T in CELLS


def trunk_junctions(src):
    """Walk trunk downstream from src; yield (node, trunk_area_m2, trib_area_m2,
    dist_from_src_m) at every tributary junction (trib >= TRIB_FRAC of trunk)."""
    cur = src; d0 = fd[src]; out = []
    while True:
        nx = int(recv[cur])
        if nx == cur:
            break
        trib_cells = ncon[nx] - ncon[cur] - 1            # cells entering at nx that are NOT the trunk or nx itself
        if trib_cells > 0 and trib_cells >= TRIB_FRAC * ncon[cur]:
            out.append((nx, ncon[nx] * pxa, trib_cells * pxa, d0 - fd[nx]))
        cur = nx
    return out


srcs = D.get_sources(fi, Tc)
# order sources by area so the report is stable; take a handful
srcs = sorted(srcs.tolist(), key=lambda s: ncon[s])[:8]
for s in srcs:
    js = trunk_junctions(s)
    print('\nsource ncon=%.0f m2  (%d trunk junctions >= %.0f%% of trunk):'
          % (ncon[s] * pxa, len(js), 100 * TRIB_FRAC))
    for node, ta, tra, dd in js[:6]:
        print('   @%.0f m downstream: trunk=%.0f m2  +tributary=%.0f m2  (%.0f%% jump)'
              % (dd, ta, tra, 100 * tra / (ta - tra)))
    if js:
        node, ta, tra, dd = js[0]
        print('   -> FIRST tributary junction bounds the clean reach at %.0f m down, trunk=%.0f m2' % (dd, ta))
