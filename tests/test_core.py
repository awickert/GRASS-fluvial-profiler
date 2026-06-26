"""
Tests for rivernetworkx.core — pure NetworkX, no GRASS required.

Network used (segments = edges, junctions = nodes; cat -> tostream):

       1            2          headwaters (cats 1, 2) at y=50
        \\          /
         \\        /
          3 (junction node), at (0, 30)
          |
          0 (outlet / off-map), at (0, 0)

Segment lengths: seg 3 = 30 (vertical); segs 1 and 2 = hypot(20,20) = 28.284...
Cumulative distance upstream of the outlet (node 's'):
    node 0 = 0 ; node 3 = 30 ; nodes 1,2 = 30 + 28.284 = 58.284
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import rivernetworkx as rnx  # noqa: E402

DIAG = np.hypot(20.0, 20.0)  # 28.2842...

RECORDS = [
    # seg 3: junction (0,30) -> outlet (0,0); ordered upstream -> downstream
    {"cat": 3, "tostream": 0, "x": [0.0, 0.0], "y": [30.0, 0.0], "z": [3.0, 0.0]},
    # seg 1: headwater (-20,50) -> junction (0,30)
    {"cat": 1, "tostream": 3, "x": [-20.0, 0.0], "y": [50.0, 30.0], "z": [9.0, 3.0]},
    # seg 2: headwater (20,50) -> junction (0,30)
    {"cat": 2, "tostream": 3, "x": [20.0, 0.0], "y": [50.0, 30.0], "z": [9.0, 3.0]},
]


def test_segment_distances():
    s_down, s_up = rnx.segment_distances([0.0, 0.0], [30.0, 0.0])
    assert np.allclose(s_down, [0.0, 30.0])
    assert np.allclose(s_up, [30.0, 0.0])


def test_topology():
    G = rnx.build_graph(RECORDS)
    assert set(G.nodes) == {0, 1, 2, 3}
    assert set(G.edges) == {(1, 3), (2, 3), (3, 0)}
    # headwaters have nothing upstream; node 3 is a junction
    assert [n for n in G if G.in_degree(n) == 0] == [1, 2] or \
           sorted(n for n in G if G.in_degree(n) == 0) == [1, 2]
    assert G.in_degree(3) == 2


def test_cumulative_distance():
    G = rnx.build_graph(RECORDS)
    assert np.isclose(G.nodes[0]["s"][0], 0.0)
    assert np.isclose(G.nodes[3]["s"][0], 30.0)
    assert np.isclose(G.nodes[1]["s"][0], 30.0 + DIAG)
    assert np.isclose(G.nodes[2]["s"][0], 30.0 + DIAG)


def test_traversal():
    G = rnx.build_graph(RECORDS)
    assert rnx.upstream_subnetwork(G, 3).number_of_nodes() == 3  # nodes 1, 2, 3
    assert rnx.downstream_path(G, 1, outlet=0) == [1, 3, 0]


def test_assemble_downstream_profile():
    # Records indexed by cat; downstream path from headwater 1 -> junction 3.
    by_cat = {r["cat"]: r for r in RECORDS}
    G = rnx.build_graph(RECORDS)
    path = rnx.downstream_path(G, 1, outlet=0)[:-1]  # drop the off-map outlet
    assert path == [1, 3]
    prof = rnx.assemble_downstream_profile(by_cat, path, attrs=("z",))
    # seg 1 has 2 vertices, seg 3 has 2; one shared junction dropped -> 3 points
    assert len(prof["s"]) == 3 and len(prof["z"]) == 3
    # distance from mouth: headwater at 30 + DIAG, mouth at 0
    assert np.isclose(prof["s"].max(), 30.0 + DIAG)
    assert np.isclose(prof["s"].min(), 0.0)
    # elevation falls downstream: 9 (headwater) -> 0 (mouth)
    assert np.isclose(prof["z"][0], 9.0) and np.isclose(prof["z"][-1], 0.0)


def test_densify():
    s = np.array([0.0, 10.0])
    new_s, arr = rnx.densify(s, {"z": np.array([0.0, 100.0])}, dx_target=2.5)
    assert len(new_s) == 5  # ceil(10/2.5)+1
    assert np.allclose(new_s, [0.0, 2.5, 5.0, 7.5, 10.0])
    assert np.allclose(arr["z"], [0.0, 25.0, 50.0, 75.0, 100.0])  # linear interp


def test_moving_average():
    s = np.arange(5.0)            # 0,1,2,3,4
    y = np.array([0.0, 10.0, 0.0, 10.0, 0.0])
    sm = rnx.moving_average(s, y, window=2.0)  # +/- 1 each side
    assert np.isclose(sm[2], np.mean([10.0, 0.0, 10.0]))  # neighbors at 1,2,3
    # NaNs are ignored, not propagated
    y2 = np.array([0.0, np.nan, 4.0])
    sm2 = rnx.moving_average(np.arange(3.0), y2, window=10.0)
    assert np.isclose(sm2[0], 2.0)


def test_json_roundtrip(tmp_path=None):
    G = rnx.build_graph(RECORDS)
    path = os.path.join(tmp_path or "/tmp", "rnx_roundtrip.json")
    rnx.export_json(G, path)
    H = rnx.load_json(path)
    assert set(H.nodes) == set(G.nodes)
    assert set(H.edges) == set(G.edges)
    assert np.isclose(H.nodes[1]["s"][0], 30.0 + DIAG)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("PASS", name)
    print("All core tests passed.")
