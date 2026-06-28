"""How far do DrEICH channel heads move if the D8 routing comes from
r.watershed -s (single-flow-direction, steepest descent) instead of the
faithful LSDFlowInfo steepest-descent-on-filled port?

This is the UNTESTED variant: prior diagnostics used r.watershed's DEFAULT
(MFD / least-cost AT search, ~40 m). The -s flag is the closest GRASS analog to
LSDTopoTools' routing, so it is the fair routing-swap to measure.

Method: build a FlowInfo from r.watershed -s drainage directions, then run the
SAME faithful head-finder (rivernetworkx.dreich), holding curvature and
elevation fixed at the LSDTT-faithful fields so ONLY routing changes. Compare
head positions to the 634-head C++ reference (the same metric as the other
dev/dreich_routing_*.py diagnostics).

Run (needs the /tmp reference rasters + system-python GRASS, see
grass-gunittest-recipe / dreich-port-status memory):
    PROJ_DATA=/usr/share/proj grass -c XY /tmp/rwtest/PERMANENT --exec \
        /usr/bin/python3 dev/dreich_routing_rwatershed_test.py [raw|filled]
"""
import sys
import numpy as np
from grass.script import run_command

from rivernetworkx import dreich as D

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0
RES = 1.0
NORTH = 4369656.1523925
WEST = 398767.32685636
# DrEICH parameters (LSDTT defaults), curvature threshold as in the faithful run.
THRESHOLD = 100
A0 = 1000.0
MN = 0.525
NCONNECT = 10
MINSEG = 10
CURV_THRESH = 0.1

# r.watershed drainage direction: 1..8 CCW from NE (verified empirically:
# east-flowing cells = 8). |dir| gives the neighbour; dir<0 marks a cell whose
# receiver leaves the region. abs(dir) -> (drow, dcol):
#   1 NE  2 N  3 NW  4 W  5 SW  6 S  7 SE  8 E
_DR = np.array([0, -1, -1, -1, 0, 1, 1, 1, 0])      # index by |dir| (0..8)
_DC = np.array([0, 1, 0, -1, -1, -1, 0, 1, 1])
_CARDINAL = {2, 4, 6, 8}


def fi_from_directions(dirgrid, dem, ndf=ND, res=RES):
    """Assemble a build_flowinfo-compatible FlowInfo from an external D8
    drainage-direction raster (r.watershed encoding). Only recv / flc come from
    the directions; everything else (node indexing, donor stack) is built exactly
    as rivernetworkx.dreich.build_flowinfo does."""
    z = dem.astype(np.float32)
    nr, nc = z.shape
    valid = z != np.float32(ndf)
    NodeIndex = np.full((nr, nc), -1, dtype=np.int64)
    rr, cc = np.where(valid)
    N = rr.size
    NodeIndex[rr, cc] = np.arange(N)

    a = np.abs(dirgrid[rr, cc]).astype(np.int64)
    a = np.where((a >= 1) & (a <= 8), a, 0)             # 0 / NULL dir -> outlet
    tr = rr + _DR[a]
    tc = cc + _DC[a]
    ingrid = (a > 0) & (tr >= 0) & (tr < nr) & (tc >= 0) & (tc < nc)
    trc = np.clip(tr, 0, nr - 1)
    tcc = np.clip(tc, 0, nc - 1)
    routes = ingrid & valid[trc, tcc]                  # receiver is a real cell
    recv = np.where(routes, NodeIndex[trc, tcc], np.arange(N))   # else self (outlet)
    flc = np.where(routes, np.where(np.isin(a, list(_CARDINAL)), 1, 2), 0).astype(np.int8)

    ndon = np.zeros(N, dtype=np.int64)
    np.add.at(ndon, recv, 1)
    delta = np.zeros(N + 1, dtype=np.int64)
    delta[1:] = np.cumsum(ndon)
    donorstack = np.argsort(recv, kind='stable')
    return dict(nr=nr, nc=nc, ndf=np.float32(ndf), res=float(res), N=N,
                NodeIndex=NodeIndex, row_of=rr.astype(np.int64),
                col_of=cc.astype(np.int64), recv=recv, flc=flc, ndon=ndon,
                delta=delta, donorstack=donorstack,
                baselevel=np.where(recv == np.arange(N))[0])


def heads_from_fi(fi, filled, tcurv):
    """Run the faithful DrEICH head-finder on a given FlowInfo."""
    D.contributing_area(fi)
    sources = D.get_sources(fi, THRESHOLD)
    D.junction_network(fi, sources)
    valley = D.find_valleys(fi, tcurv, sources, NCONNECT, CURV_THRESH)
    D.build_svector(fi)
    D.distance_from_outlet(fi)
    _, final = D.channel_heads_from_valleys(fi, filled, valley, MINSEG, A0, MN)
    return [(int(fi['row_of'][n]), int(fi['col_of'][n])) for n in final]


def main():
    surface = sys.argv[1] if len(sys.argv) > 1 else 'raw'
    routing = sys.argv[2] if len(sys.argv) > 2 else 'sfd'   # sfd (-s) or mfd (default)
    run_command('g.region', n=NORTH, s=4362549.1523925, w=WEST,
                e=405033.32685636, nsres=RES, ewres=RES, quiet=True)

    src_flt = ('bailey_run_dem.flt' if surface == 'raw'
               else 'bailey_run_dem_fill.flt')
    print('r.watershed (%s) on the %s DEM (%s)' % (routing, surface, src_flt))
    run_command('r.in.bin', flags='f', input='%s/%s' % (ALG, src_flt),
                output='bailey_dem', bytes=4, anull=ND,
                north=NORTH, south=4362549.1523925, east=405033.32685636,
                west=WEST, rows=NR, cols=NC, overwrite=True, quiet=True)
    rw_flags = 's' if routing == 'sfd' else ''
    run_command('r.watershed', flags=rw_flags, elevation='bailey_dem',
                drainage='rw_dir', overwrite=True, quiet=True)
    run_command('r.out.bin', input='rw_dir', output='/tmp/rw_dir.bin',
                bytes=4, null=0, quiet=True)
    dirgrid = np.fromfile('/tmp/rw_dir.bin', dtype='<i4').reshape(NR, NC)
    print('  directions: %d routed, %d off-edge(neg), %d null/zero'
          % (int((dirgrid > 0).sum()), int((dirgrid < 0).sum()),
             int((dirgrid == 0).sum())))

    filled = np.fromfile('%s/bailey_run_dem_fill.flt' % ALG,
                         dtype='<f4').reshape(NR, NC)
    tcurv = np.fromfile('%s/bailey_run_dem_tan_curv.flt' % ALG,
                        dtype='<f4').reshape(NR, NC)

    fi = fi_from_directions(dirgrid, filled)
    print('  FlowInfo: %d nodes, %d baselevel' % (fi['N'], fi['baselevel'].size))
    heads = heads_from_fi(fi, filled, tcurv)
    print('r.watershed (%s) heads: %d  (C++ reference 634)' % (routing, len(heads)))

    E = np.array([WEST + (c + 0.5) * RES for (r, c) in heads])
    Nn = np.array([NORTH - (r + 0.5) * RES for (r, c) in heads])
    ref = np.load('/tmp/dreich_ref_heads.npy')
    d = np.array([np.min(np.hypot(ref[:, 0] - e, ref[:, 1] - n))
                  for e, n in zip(E, Nn)])
    print('r.watershed(%s)-vs-reference (Euclidean): median %.1f m | '
          '<=5m %.0f%% | <=10m %.0f%% | <=20m %.0f%%'
          % (routing, np.median(d), 100 * np.mean(d <= 5),
             100 * np.mean(d <= 10), 100 * np.mean(d <= 20)))
    np.save('/tmp/dreich_rwatershed_%s_%s_heads.npy' % (surface, routing),
            np.column_stack([E, Nn]))


if __name__ == '__main__':
    main()
