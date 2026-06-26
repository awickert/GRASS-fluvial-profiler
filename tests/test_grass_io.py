"""
Tests for the pure (GRASS-free) pieces of rivernetworkx.grass_io:
sample_raster (coord -> index math) and assemble_records, plus an end-to-end
assemble -> build_graph with sampled values.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import rivernetworkx as rnx  # noqa: E402


def test_sample_raster_nearest():
    # 3 rows x 4 cols, 10 m cells; region west=0, north=30 -> covers x[0,40], y[0,30]
    arr = np.arange(12, dtype=float).reshape(3, 4)  # arr[0,0]=0 (NW), arr[2,3]=11 (SE)
    bounds = dict(west=0.0, north=30.0, nsres=10.0, ewres=10.0)
    # (5,25)->row0,col0=0 ; (35,5)->row2,col3=11 ; (15,15)->row1,col1=5
    vals = rnx.sample_raster(arr, [5.0, 35.0, 15.0], [25.0, 5.0, 15.0], **bounds)
    assert np.allclose(vals, [0.0, 11.0, 5.0])


def test_sample_raster_out_of_range_is_nan():
    arr = np.arange(12, dtype=float).reshape(3, 4)
    bounds = dict(west=0.0, north=30.0, nsres=10.0, ewres=10.0)
    vals = rnx.sample_raster(arr, [-5.0, 100.0], [25.0, 25.0], **bounds)
    assert np.all(np.isnan(vals))


def test_sample_raster_on_east_south_edge():
    # Points exactly on the outer east/south edges belong to the last cell,
    # not "outside" (regression: floor() lands one index past the end).
    arr = np.arange(12, dtype=float).reshape(3, 4)  # east=40, south=0
    bounds = dict(west=0.0, north=30.0, nsres=10.0, ewres=10.0)
    # (40,15) -> col clipped 3, row 1 -> arr[1,3]=7 ; (20,0) -> row clipped 2 -> arr[2,2]=10
    vals = rnx.sample_raster(arr, [40.0, 20.0], [15.0, 0.0], **bounds)
    assert np.allclose(vals, [7.0, 10.0])
    # just past the edge is still NaN
    assert np.all(np.isnan(rnx.sample_raster(arr, [40.001], [15.0], **bounds)))


def test_offmap_inflow_cats():
    # Vertices are ordered upstream -> downstream. Negative accumulation only at
    # the DOWNSTREAM end (the outlet leaving the map) is normal and allowed;
    # negative reaching the channel HEAD flags an incomplete catchment.
    recs = [
        {'cat': 1, 'A': np.array([5.0, 10.0, -3.0])},    # neg only at outlet: ok
        {'cat': 2, 'A': np.array([-1.0, 3.0, 8.0])},     # neg at head: off-map inflow
        {'cat': 3, 'A': np.array([4.0, 9.0])},           # all positive: ok
        {'cat': 4},                                       # no A sampled: ok
    ]
    assert rnx.offmap_inflow_cats(recs) == [2]


def test_assemble_records_skips_none():
    cats = [1, 3]
    tostream = {1: 3, 3: 0}
    geometry = {1: (np.array([-20.0, 0.0]), np.array([50.0, 30.0])),
                3: (np.array([0.0, 0.0]), np.array([30.0, 0.0]))}
    zmap = {1: np.array([9.0, 3.0]), 3: np.array([3.0, 0.0])}
    recs = rnx.assemble_records(cats, tostream, geometry, z=zmap, A=None)
    assert recs[0]['cat'] == 1 and recs[0]['tostream'] == 3
    assert 'z' in recs[0] and 'A' not in recs[0]


def test_assemble_then_build():
    # sampled records flow straight into the pure builder
    cats = [3, 1, 2]
    tostream = {3: 0, 1: 3, 2: 3}
    geometry = {
        3: (np.array([0.0, 0.0]), np.array([30.0, 0.0])),
        1: (np.array([-20.0, 0.0]), np.array([50.0, 30.0])),
        2: (np.array([20.0, 0.0]), np.array([50.0, 30.0])),
    }
    zmap = {3: np.array([3.0, 0.0]), 1: np.array([9.0, 3.0]), 2: np.array([9.0, 3.0])}
    recs = rnx.assemble_records(cats, tostream, geometry, z=zmap)
    G = rnx.build_graph(recs)
    assert set(G.edges) == {(1, 3), (2, 3), (3, 0)}
    assert np.isclose(G.nodes[3]['s'][0], 30.0)
    assert np.isclose(G.nodes[1]['s'][0], 30.0 + np.hypot(20.0, 20.0))


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("PASS", name)
    print("All grass_io (pure) tests passed.")
