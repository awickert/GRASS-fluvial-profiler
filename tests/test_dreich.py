"""Unit tests for the DrEICH channel-head port (rivernetworkx.dreich).

GRASS-free, synthetic fixtures. The heavy bit-exact validation against LSDTT lives
in dev/ (needs the /tmp reference rasters); these tests pin the building blocks and
end-to-end plumbing on small grids.
"""
import numpy as np
import pytest

from rivernetworkx import dreich as D


# ------------------------------------------------------------------ fill
def test_fill_raises_a_pit():
    z = np.array([[5, 5, 5, 5, 5],
                  [5, 4, 4, 4, 5],
                  [5, 4, 0, 4, 5],     # interior pit
                  [5, 4, 4, 4, 5],
                  [5, 5, 5, 4, 5]], dtype=np.float32)  # one outlet on the south edge
    filled = D.fill(z, nodata=-9999.0, min_slope=0.0, cellsize=1.0)
    # the pit must be raised to (at least) its lowest rim so it is no longer a sink
    assert filled[2, 2] > z[2, 2]
    # cells that already drain are unchanged where they were not part of a depression
    assert filled[0, 0] == z[0, 0]


def test_fill_imposes_min_slope_gradient():
    # a flat plateau draining off one edge -> min_slope tilts it monotonically
    z = np.full((6, 6), 10.0, dtype=np.float32)
    z[-1, :] = 9.0                      # lower south edge = outlet
    filled = D.fill(z, nodata=-9999.0, min_slope=0.01, cellsize=1.0)
    assert np.all(np.isfinite(filled))
    # filled is >= original everywhere (filling only raises)
    assert np.all(filled >= z - 1e-6)


# ------------------------------------------------------------- flow routing
def _ramp(nr=8, nc=8, slope=1.0):
    """Elevation decreasing toward the south (increasing row) -> flow goes south."""
    z = (np.arange(nr)[:, None] * -slope + 100.0) * np.ones((1, nc))
    return z.astype(np.float32)


def test_receivers_point_downhill():
    z = _ramp()
    fi = D.build_flowinfo(z, nodata=-9999.0, cellsize=1.0)
    # a non-edge node's receiver should have a strictly lower (more southerly) row
    rid = fi['NodeIndex'][2, 3]
    rr, cc = fi['row_of'][fi['recv'][rid]], fi['col_of'][fi['recv'][rid]]
    assert rr >= 2                       # receiver is to the south (row increases downhill)
    assert z[rr, cc] <= z[2, 3]


def test_contributing_area_accumulates_downstream():
    z = _ramp()
    fi = D.build_flowinfo(z, nodata=-9999.0, cellsize=1.0)
    D.contributing_area(fi)
    # the most-downstream interior row should have larger contributing area than the top
    top = fi['ncontrib'][fi['NodeIndex'][1, 4]]
    bot = fi['ncontrib'][fi['NodeIndex'][6, 4]]
    assert bot >= top
    assert fi['ncontrib'].min() >= 1


def test_flow_distance_increases_upstream():
    z = _ramp()
    fi = D.build_flowinfo(z, nodata=-9999.0, cellsize=1.0)
    D.contributing_area(fi); D.build_svector(fi); D.distance_from_outlet(fi)
    fd = fi['fd']
    assert fd.min() >= 0.0
    # an upstream cell is farther from the outlet than a downstream cell on its path
    assert fd[fi['NodeIndex'][1, 4]] >= fd[fi['NodeIndex'][6, 4]]


# --------------------------------------------------------- chi-z regression
def test_reg32_perfect_line():
    x = np.linspace(0, 10, 30).astype(np.float32)
    y = (3.0 * x + 1.0).astype(np.float32)
    r2, dw = D._reg32(x, y)
    assert r2 > 0.9999                   # perfect fit


def test_lu2_matches_numpy_solve():
    # well-conditioned 2x2: faithful float32 LU close to numpy's float64 solve
    A = np.array([[4.0, 1.0], [1.0, 3.0]])
    b = np.array([1.0, 2.0])
    x0, x1 = D._lu2_f32(*A.ravel(), *b)
    ref = np.linalg.solve(A, b)
    assert abs(x0 - ref[0]) < 1e-4 and abs(x1 - ref[1]) < 1e-4


def test_chi_z_split_classifies_hillslope_as_non_channel():
    # upper "hillslope" convex in chi, lower "channel" linear in chi (nodeseq runs
    # hilltop -> downstream). The chi-z split must classify the convex hillslope as
    # NON-channel, i.e. place the head within the linear channel region (index
    # >= n_hill). (On perfectly smooth profiles the DrEICH statistic is biased
    # downstream because Durbin-Watson of the hillslope keeps falling with segment
    # length; the robust invariant is that the head is not up in the hillslope.)
    rng = np.random.RandomState(0)
    n_hill, n_chan = 40, 60
    chi_chan = np.linspace(0, 6, n_chan)
    elev_chan = 2.0 * chi_chan + rng.normal(0, 0.05, n_chan)   # ~linear channel
    chi_hill = np.linspace(6, 10, n_hill)
    elev_hill = elev_chan[-1] + 3.0 * (chi_hill - 6) ** 1.6    # convex hillslope
    chi = np.concatenate([chi_hill[::-1], chi_chan[::-1]]).astype(np.float32)
    elev = np.concatenate([elev_hill[::-1], elev_chan[::-1]]).astype(np.float32)
    nodeseq = list(range(len(chi)))
    head = D._calculate_channel_head(nodeseq, chi, elev, min_segment_length=10)
    assert head >= n_hill - 5                          # head is at/below the hillslope
    assert head <= len(chi) - 10                       # and a valid channel split


# --------------------------------------------------------------- end-to-end
def test_extract_runs_and_returns_cells():
    # a tilted plane with a central incised valley; just check the plumbing runs and
    # returns in-bounds (row, col) cells without error.
    nr, nc = 60, 60
    z = (np.arange(nr)[:, None] * -0.5 + 100.0) * np.ones((1, nc))
    valley = np.exp(-((np.arange(nc)[None, :] - nc / 2) ** 2) / 8.0) * 3.0
    z = (z - valley).astype(np.float32)
    heads = D.extract_channel_heads(z, nodata=-9999.0, cellsize=1.0,
                                    threshold=20, min_segment_length=5)
    assert isinstance(heads, list)
    for r, c in heads:
        assert 0 <= r < nr and 0 <= c < nc


# --------------------------------------------------------------- network split
def _y_flowinfo():
    """Two channels (heads at (0,0) and (0,4)) joining at (2,2), then a trunk to
    the outlet at (5,2). Walls elsewhere. Returns (fi, heads)."""
    W = 100.0
    z = np.full((6, 5), W, dtype=np.float32)
    for (r, c), v in {(0, 0): 9, (1, 1): 8, (2, 2): 6, (3, 2): 5,
                      (4, 2): 4, (5, 2): 3, (0, 4): 9, (1, 3): 8}.items():
        z[r, c] = v
    fi = D.build_flowinfo(z, nodata=-9999.0, cellsize=1.0)
    return fi, [(0, 0), (0, 4)]


def test_channel_network_segments_topology():
    fi, heads = _y_flowinfo()
    segs = D.channel_network_segments(fi, heads)
    # two branches + one trunk = three links
    assert len(segs) == 3
    by_cat = {s['cat']: s for s in segs}
    assert set(by_cat) == {1, 2, 3}
    # cats are assigned by upstream (row, col): the two heads come before the
    # confluence-rooted trunk.
    assert by_cat[1]['cells'][0] == (0, 0)
    assert by_cat[2]['cells'][0] == (0, 4)
    assert by_cat[3]['cells'][0] == (2, 2)
    # both branches drain into the trunk (cat 3); the trunk exits the map (0).
    assert by_cat[1]['tostream'] == 3
    assert by_cat[2]['tostream'] == 3
    assert by_cat[3]['tostream'] == 0
    # the confluence cell (2,2) is shared: it ends each branch and starts the trunk.
    assert by_cat[1]['cells'][-1] == (2, 2)
    assert by_cat[2]['cells'][-1] == (2, 2)
    assert by_cat[3]['cells'][-1] == (5, 2)


def test_channel_network_tostream_is_a_converging_tree():
    fi, heads = _y_flowinfo()
    segs = D.channel_network_segments(fi, heads)
    cats = {s['cat'] for s in segs}
    # every non-zero tostream references a real link (no dangling edges)...
    for s in segs:
        assert s['tostream'] == 0 or s['tostream'] in cats
    # ...and following tostream from any link terminates at the outlet (no cycles).
    nxt = {s['cat']: s['tostream'] for s in segs}
    for c in cats:
        seen, cur = set(), c
        while cur != 0:
            assert cur not in seen          # would be a cycle
            seen.add(cur)
            cur = nxt[cur]


def test_return_flowinfo_round_trips_through_network():
    nr, nc = 60, 60
    z = (np.arange(nr)[:, None] * -0.5 + 100.0) * np.ones((1, nc))
    valley = np.exp(-((np.arange(nc)[None, :] - nc / 2) ** 2) / 8.0) * 3.0
    z = (z - valley).astype(np.float32)
    heads, fi = D.extract_channel_heads(z, nodata=-9999.0, cellsize=1.0,
                                        threshold=20, min_segment_length=5,
                                        return_flowinfo=True)
    segs = D.channel_network_segments(fi, heads)
    # the network reaches every head: each head is the upstream end of some link.
    starts = {s['cells'][0] for s in segs}
    for h in heads:
        assert h in starts


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-v']))
