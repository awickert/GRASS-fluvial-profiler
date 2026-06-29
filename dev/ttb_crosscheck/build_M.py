"""Cross-check harness, Python side. Read the TopoToolbox example DEM (Big
Tujunga 30 m), clip a clean no-nodata window, fill + route with dreich, and build
the sparse giver->receiver flow matrix M (MATLAB column-major linear indices) for
FLOWobj(M,...) in Octave -- so TopoToolbox DIVIDEobj routes on the IDENTICAL flow
directions as my divide pipeline. Also compute my divide network on the same
window and save it. Comparison is done in ttb_compare.py after the Octave run.

Run:  PYTHONPATH=. /usr/bin/python3 dev/ttb_compare_build.py
"""
import numpy as np
import rasterio
from scipy import sparse
from scipy.io import savemat
from rivernetworkx import dreich as D
from rivernetworkx import divides as DV

TIF = '/tmp/topotoolbox/DEMdata/srtm_bigtujunga30m_utm11.tif'
ND = -9999.0
THRESH = 200          # channel area threshold (cells), same in both tools
# clip a clean interior window (row0, row1, col0, col1)
R0, R1, C0, C1 = 150, 450, 300, 700

with rasterio.open(TIF) as src:
    full = src.read(1).astype(np.float64)
    cs = float(src.transform.a)
    src_nd = src.nodata
print('full DEM %s cellsize=%g nodata=%s' % (full.shape, cs, src_nd))

dem = full[R0:R1, C0:C1].copy()
if src_nd is not None:
    dem[dem == src_nd] = np.nan
print('window %s  nan cells: %d' % (dem.shape, np.isnan(dem).sum()))
assert not np.isnan(dem).any(), 'window has nodata; pick another'
nr, nc = dem.shape

filled = D.fill(dem.astype(np.float32), ND, 0.0001, cs)
fi = D.build_flowinfo(filled, ND, cs)
D.contributing_area(fi)
recv = fi['recv']; row = fi['row_of']; col = fi['col_of']
N = nr * nc

# MATLAB column-major linear index (0-based) for cell (r,c): c*nr + r
def mlin(r, c):
    return c * nr + r

g = mlin(row, col)
r_ = mlin(row[recv], col[recv])
moving = recv != np.arange(fi['N'])          # exclude baselevel self-loops
M = sparse.csc_matrix((np.ones(moving.sum()), (g[moving], r_[moving])), shape=(N, N))
print('M: %d nonzeros (givers), shape %dx%d' % (M.nnz, N, N))
savemat('/tmp/ttb_M.mat', {'M': M, 'nr': nr, 'nc': nc, 'cs': cs, 'thresh': THRESH},
        do_compression=True)

# ---- my divide network on the same window/routing ----
_, lab = D.drainage_divides(fi, threshold=THRESH)
stream = np.zeros((nr, nc), bool)
stream[row, col] = fi['ncontrib'] >= THRESH
edges, meta = DV.extract_divide_edges(lab, stream=stream)
fx = DV.diagonal_flow_corners(recv, row, col, nr, nc)
edges = DV.split_diagonal_crossings(edges, fx, meta['cw'])
gph = DV.build_divide_graph(edges)
order, dd = DV.order_divides(gph, scheme='topo')
# save my divide edges as (i,j) corner pairs + per-edge order
mine = []
for k, p in enumerate(gph['segments']):
    ii, jj = DV.decode_corner(np.asarray(p) // 2, meta['cw'])
    for a in range(len(p) - 1):
        mine.append((ii[a], jj[a], ii[a + 1], jj[a + 1], order[k]))
mine = np.array(mine, dtype=np.int64)
np.save('/tmp/ttb_mine_edges.npy', mine)
print('my divides: %d segments, %d edges, max order %d, %d unordered'
      % (len(gph['segments']), len(mine), order.max(), int((order == 0).sum())))
print('SAVED /tmp/ttb_M.mat and /tmp/ttb_mine_edges.npy')
