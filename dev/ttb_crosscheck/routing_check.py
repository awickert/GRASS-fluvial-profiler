"""Routing-identity check: does TopoToolbox's flow accumulation on FLOWobj(M)
match dreich's ncontrib cell-for-cell? If yes, routing is provably shared and the
divide residual is construction (per-link trim vs per-order getdivide), not
routing -- the LSDTopoTools lesson, ruled out by evidence. Run after
ttb_routing_check.m has written /tmp/ttb_flowacc.f64.

Run:  PYTHONPATH=/path/to/r.fluvial /usr/bin/python3 routing_check.py
"""
import numpy as np
import rasterio
from rivernetworkx import dreich as D

TIF = '/tmp/topotoolbox/DEMdata/srtm_bigtujunga30m_utm11.tif'
ND = -9999.0
R0, R1, C0, C1 = 150, 450, 300, 700        # same window as build_M.py

with rasterio.open(TIF) as src:
    dem = src.read(1).astype(np.float64)[R0:R1, C0:C1]
    cs = float(src.transform.a)
    if src.nodata is not None:
        dem[dem == src.nodata] = np.nan
nr, nc = dem.shape
filled = D.fill(dem.astype(np.float32), ND, 0.0001, cs)
fi = D.build_flowinfo(filled, ND, cs)
D.contributing_area(fi)
mine = np.zeros((nr, nc))
mine[fi['row_of'], fi['col_of']] = fi['ncontrib']

ttb = np.fromfile('/tmp/ttb_flowacc.f64', '<f8').reshape((nr, nc), order='F')

print('grid %dx%d' % (nr, nc))
print('mine  ncontrib: min=%g max=%g' % (mine.min(), mine.max()))
print('ttb   flowacc : min=%g max=%g' % (ttb.min(), ttb.max()))
d = ttb - mine
print('difference (ttb - mine): exact-equal cells %.4f%%, max|diff|=%g, mean|diff|=%g'
      % (100 * np.mean(d == 0), np.abs(d).max(), np.abs(d).mean()))
# common off-by-one? TopoToolbox flowacc counts self; check ttb == mine and ttb == mine (+/-)
for off in (0, 1, -1):
    print('  frac cells with ttb == mine%+d : %.4f' % (off, np.mean(ttb == mine + off)))
# where do they differ, and by how much?
nz = d != 0
if nz.any():
    print('differing cells: %d (%.3f%%); diff distribution:' % (nz.sum(), 100*nz.mean()))
    vals, cnts = np.unique(d[nz], return_counts=True)
    order = np.argsort(-cnts)
    for k in order[:8]:
        print('    diff=%+g : %d cells' % (vals[k], cnts[k]))
