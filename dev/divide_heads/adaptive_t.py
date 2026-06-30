"""Adaptive T-finder (Andy's idea): sweep the drainage-area threshold T and look
for the scale-invariant regime -- where the landscape is organized into stable
first-order catchments -- BEFORE the hillslope breakdown, with NO channel heads.

Hypothesis: as T decreases, the network first resolves the first-order
catchments, then 'settles' (a scale-invariant / power-law regime where lowering T
just subdivides self-similarly), then breaks down into tiny hillslope bits. The
right T sits at the transition out of the fluvial regime into hillslope.

Diagnostics vs T (all from the flow routing alone -- no heads):
  * n_src(T)   -- number of first-order sources (channel-initiation points)
  * n_chan(T)  -- cells with contributing area >= T (drainage-density / total
                  channel-length proxy)
  * log-log local slopes d ln n / d ln T -- a constant slope = scale-invariant
    'settled' regime; a break/knee = transition (basin scale at large T,
    hillslope breakdown at small T).
The field-validated valley scale (T~5000 m2 at MBR, from the field-head density)
is marked: does a statistical feature land there, with no field data used?

Run: PYTHONPATH=<repo>:<repo>/dev:<repo>/dev/divide_heads /usr/bin/python3 adaptive_t.py <mbr|indian>
"""
import sys
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from rivernetworkx import dreich as D
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
cfg = SITES[key]
buf = 400
T_FIELD = 5000.0                       # field-validated valley scale (m^2)

z, west, north, res, nd = read_flt(cfg['flt'])
z = np.where(z == nd, np.nan, z)
field = load_field_heads(cfg['key'])
nr, nc = z.shape
cr = ((north - field[:, 1]) / res).astype(int); cc = ((field[:, 0] - west) / res).astype(int)
ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < nc)
R0, R1 = max(0, cr[ins].min() - buf), min(nr, cr[ins].max() + buf)
C0, C1 = max(0, cc[ins].min() - buf), min(nc, cc[ins].max() + buf)
zc = z[R0:R1, C0:C1]
dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32)
area_km2 = np.isfinite(zc).sum() * res * res / 1e6
print('%s: clip %s = %.2f km2 valid; building flow routing...' % (cfg['name'], zc.shape, area_km2), flush=True)

filled = D.fill(dem, nd, cellsize=res)
fi = D.build_flowinfo(filled, nd, res)
D.contributing_area(fi)
ncontrib = fi['ncontrib']
print('  routed %d nodes; sweeping T...' % fi['N'], flush=True)

# log-spaced T (cells): from a few cells (hillslope/pixel) to large basins
T_list = np.unique(np.round(np.logspace(np.log10(4), np.log10(200000), 45)).astype(int))
rows = []
for T in T_list:
    nsrc = len(D.get_sources(fi, int(T)))
    nchan = int((ncontrib >= T).sum())
    rows.append((T, nsrc, nchan))
A = np.array(rows, float)
T_, nsrc, nchan = A[:, 0], A[:, 1], A[:, 2]

# local log-log slopes (central difference)
def slope(x, y):
    lx, ly = np.log(x), np.log(np.maximum(y, 1))
    s = np.full_like(lx, np.nan)
    s[1:-1] = (ly[2:] - ly[:-2]) / (lx[2:] - lx[:-2])
    return s
s_src = slope(T_, nsrc); s_chan = slope(T_, nchan)

print('\n  T_cell  n_src  n_chan  | dlnNsrc/dlnT  dlnNchan/dlnT')
for i in range(len(T_)):
    mark = '  <== field T~%g' % T_FIELD if T_[i - 1] < T_FIELD <= T_[i] else ''
    print('  %6d %6d %7d  |   %6.2f       %6.2f%s'
          % (T_[i], nsrc[i], nchan[i], s_src[i], s_chan[i], mark))

np.savez('/tmp/adaptive_t_%s.npz' % key, T=T_, nsrc=nsrc, nchan=nchan,
         s_src=s_src, s_chan=s_chan, area_km2=area_km2, T_field=T_FIELD, site=cfg['name'])
print('\n  field T~%g m2 -> source density %.0f /km2 (field heads ~63/km2)'
      % (T_FIELD, np.interp(T_FIELD, T_, nsrc) / area_km2))
print('saved /tmp/adaptive_t_%s.npz' % key)
