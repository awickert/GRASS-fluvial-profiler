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

import json
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
    # each point is labelled with its source segment (junction point de-duped)
    assert list(prof["cat"]) == [1, 1, 3]


def test_densify_collapses_duplicate_stations():
    # a coincident vertex gives a repeated distance; densify must not let
    # np.interp see a non-increasing x (it would return an ill-defined value)
    new_s, arr = rnx.densify([0.0, 5.0, 5.0, 10.0],
                             {"z": [0.0, 5.0, 99.0, 10.0]}, dx_target=2.5)
    assert np.allclose(new_s, [0.0, 2.5, 5.0, 7.5, 10.0])
    assert np.allclose(arr["z"], [0.0, 2.5, 5.0, 7.5, 10.0])  # 99 spike collapsed


def test_channel_slope_no_inf_on_duplicate():
    S = rnx.channel_slope([0.0, 5.0, 5.0, 10.0], [0.0, 5.0, 99.0, 10.0])
    assert not np.any(np.isinf(S))  # zero spacing -> NaN, never +/-inf


def test_export_json_is_standard_and_roundtrips_nan(tmp_path=None):
    G = rnx.build_graph(RECORDS)  # off-map node 0 carries z=[nan]
    path = os.path.join(tmp_path or "/tmp", "rnx_nan.json")
    rnx.export_json(G, path)
    raw = open(path).read()
    assert "NaN" not in raw          # standards-compliant: no bare NaN token
    json.loads(raw)                  # a strict parser accepts it
    H = rnx.load_json(path)
    assert np.isnan(H.nodes[0]["z"][0])   # null restored to NaN on load


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


def test_channel_slope():
    # z falls 6 -> 0 over s 0 -> 20: uniform descending slope of 0.3
    S = rnx.channel_slope([0.0, 10.0, 20.0], [6.0, 3.0, 0.0])
    assert np.allclose(S, 0.3)
    # too-short segment -> NaN
    assert np.all(np.isnan(rnx.channel_slope([0.0], [5.0])))


def test_smooth_segment():
    rec = {"x": [0.0, 0.0, 0.0], "y": [20.0, 10.0, 0.0],
           "z": [0.0, 10.0, 0.0]}  # a spike at the middle vertex
    # window/2 = 25 >= the 20-long segment, so every vertex averages all three
    out = rnx.smooth_segment(rec, ("z",), window=50.0)
    assert np.allclose(out["z"], 10.0 / 3.0)


def test_slope_area():
    # one segment, z 2->0 over s 0->20 (S=0.1), area 10->20
    recs = [{"cat": 1, "tostream": 0, "x": [0.0, 0.0], "y": [20.0, 0.0],
             "z": [2.0, 0.0], "A": [10.0, 20.0]}]
    A, S = rnx.slope_area(recs)
    assert np.allclose(S, 0.1)
    assert np.allclose(A, [10.0, 20.0])
    logA, logS = rnx.slope_area(recs, log=True)
    assert np.allclose(logA, np.log10([10.0, 20.0]))
    assert np.allclose(logS, np.log10([0.1, 0.1]))
    # nonpositive area/slope dropped in log space
    recs2 = [{"cat": 1, "tostream": 0, "x": [0.0, 0.0], "y": [20.0, 0.0],
              "z": [0.0, 0.0], "A": [-1.0, 5.0]}]  # S=0 and one A<0
    logA2, logS2 = rnx.slope_area(recs2, log=True)
    assert len(logA2) == 0


def test_moving_average_matches_bruteforce():
    # the vectorized O(n log n) implementation must match the naive per-point
    # window mean, including NaN handling and unsorted input
    s = np.array([0.0, 1.0, 1.5, 3.0, 7.0, 7.2])
    y = np.array([1.0, np.nan, 3.0, 5.0, 2.0, 8.0])
    w = 2.0
    half = w / 2.0
    ref = []
    for si in s:
        sel = (s >= si - half) & (s <= si + half)
        vals = y[sel][~np.isnan(y[sel])]
        ref.append(np.mean(vals) if len(vals) else np.nan)
    ref = np.array(ref)
    assert np.allclose(rnx.moving_average(s, y, w), ref, equal_nan=True)
    idx = np.array([3, 0, 5, 1, 4, 2])  # unsorted input -> unsorted-matching out
    assert np.allclose(rnx.moving_average(s[idx], y[idx], w), ref[idx],
                       equal_nan=True)


def test_dense_network_degradation():
    # A 2-vertex reach shorter than dx_target / window: must not crash, must
    # keep >= 2 points, and smoothing just averages what little it has.
    s = np.array([0.0, 3.0])
    new_s, arr = rnx.densify(s, {"z": np.array([1.0, 4.0])}, dx_target=100.0)
    assert len(new_s) == 2 and np.allclose(arr["z"], [1.0, 4.0])  # endpoints kept
    rec = {"x": [0.0, 0.0], "y": [3.0, 0.0], "z": [1.0, 4.0]}
    out = rnx.smooth_segment(rec, ("z",), window=100.0)
    assert np.allclose(out["z"], 2.5)  # both averaged; no reach beyond the segment


def test_json_roundtrip(tmp_path=None):
    G = rnx.build_graph(RECORDS)
    path = os.path.join(tmp_path or "/tmp", "rnx_roundtrip.json")
    rnx.export_json(G, path)
    H = rnx.load_json(path)
    assert set(H.nodes) == set(G.nodes)
    assert set(H.edges) == set(G.edges)
    assert np.isclose(H.nodes[1]["s"][0], 30.0 + DIAG)


def _two_limb_cloud(logA_star=3.0, theta=0.5, h=-1.0, seed=0):
    """Synthetic slope-area cloud: flat hillslope below the knot, power-law above."""
    rng = np.random.default_rng(seed)
    logA = np.linspace(1.5, 5.0, 600)
    logS = np.where(logA <= logA_star, h, h - theta * (logA - logA_star))
    logS = logS + rng.normal(0.0, 0.05, logA.shape)  # mild scatter
    return logA, logS


def test_fit_sa_break_recovers_knot():
    logA, logS = _two_limb_cloud(logA_star=3.0, theta=0.5, h=-1.0)
    fit = rnx.fit_sa_break(logA, logS)
    assert fit is not None
    assert abs(fit["logA_star"] - 3.0) < 0.1     # knot within a grid step or two
    assert abs(fit["theta"] - 0.5) < 0.05        # fluvial concavity
    assert abs(fit["hillslope_logS"] - (-1.0)) < 0.05
    assert np.isclose(fit["A_star"], 10.0 ** fit["logA_star"])


def test_fit_sa_break_too_few_points():
    assert rnx.fit_sa_break([1.0, 2.0, 3.0], [0.0, -0.1, -0.2]) is None


def test_colluvial_fluvial_transition():
    recs = [
        # crosses A*=1200 between vertices at A=1000 (x=10) and A=2000 (x=20)
        {"cat": 5, "x": [0.0, 10.0, 20.0], "y": [0.0, 0.0, 0.0],
         "A": [100.0, 1000.0, 2000.0]},
        {"cat": 6, "x": [0.0, 5.0], "y": [0.0, 0.0], "A": [100.0, 800.0]},   # all below
        {"cat": 7, "x": [0.0, 5.0], "y": [0.0, 0.0], "A": [5000.0, 9000.0]},  # all above
        # reversed order (downstream first) still found and oriented
        {"cat": 8, "x": [20.0, 10.0], "y": [9.0, 9.0], "A": [2000.0, 1000.0]},
    ]
    pts = rnx.colluvial_fluvial_transition(recs, A_star=1200.0)
    by_cat = {c: (x, y) for x, y, c in pts}
    assert set(by_cat) == {5, 8}                 # only the crossing segments
    # interp: f = (1200-1000)/(2000-1000) = 0.2 -> x = 10 + 0.2*10 = 12
    assert np.isclose(by_cat[5][0], 12.0)
    assert np.isclose(by_cat[8][0], 12.0) and np.isclose(by_cat[8][1], 9.0)


def test_chi_constant_area_is_distance():
    # A == ref_area everywhere -> integrand == 1 -> chi = flow distance from the
    # downstream end, growing upstream.
    dist = np.array([0.0, 1.0, 3.0, 6.0])        # upstream -> downstream cumulative
    area = np.array([10.0, 10.0, 10.0, 10.0])    # downstream end is index -1
    c = rnx.chi(area, dist, theta=0.5, ref_area=10.0, base_chi=0.0)
    # downstream node (max-distance end here is index -1) anchors at 0;
    # chi at each node = distance back to that downstream end.
    assert np.isclose(c[-1], 0.0)
    assert np.allclose(c, [6.0, 5.0, 3.0, 0.0])


def test_chi_theta_zero_ignores_area():
    dist = np.array([0.0, 2.0, 5.0])
    area = np.array([1.0, 50.0, 900.0])          # increases downstream
    c = rnx.chi(area, dist, theta=0.0, ref_area=1.0)
    assert np.allclose(c, [5.0, 3.0, 0.0])       # = distance from downstream end


def test_chi_order_invariant():
    dist = np.array([0.0, 1.0, 2.5, 4.0])
    area = np.array([5.0, 20.0, 100.0, 500.0])   # upstream(small) -> downstream(big)
    c = rnx.chi(area, dist)
    # reversed inputs (downstream-first) must give the reversed same chi values
    c_rev = rnx.chi(area[::-1], dist[::-1])
    assert np.allclose(c, c_rev[::-1])
    # chi strictly increases upstream (toward the small-area end)
    assert np.all(np.diff(c) < 0)                # index 0 is upstream = largest chi


def test_linfit_r2_dw():
    from rivernetworkx.core import _linfit_r2_dw
    x = np.arange(10.0)
    r2, _ = _linfit_r2_dw(x, 2.0 * x + 1.0)      # perfect line
    assert r2 > 0.9999
    r2c, dwc = _linfit_r2_dw(x, x ** 2)          # convex -> line fit autocorrelated
    assert r2c < 1.0 and dwc < 2.0               # DW < 2 = positive autocorrelation


def test_channel_head_chi_split_matches_argmax_of_score():
    # Faithful implementation check: the returned index is the argmax of
    # test = R2_channel - (DW_hillslope - 2)/2 over all valid splits.
    from rivernetworkx.core import _linfit_r2_dw
    rng = np.random.RandomState(0)
    chi = np.linspace(8.0, 0.0, 50)
    z = 1.5 * chi + np.where(chi > 4.0, 0.4 * (chi - 4.0) ** 2, 0.0) \
        + rng.normal(0, 0.05, 50)
    m = 6
    best = max(
        ((_linfit_r2_dw(chi[h:], z[h:])[0]
          - (_linfit_r2_dw(chi[:h], z[:h])[1] - 2.0) / 2.0), h)
        for h in range(m, len(chi) - m + 1)
    )[1]
    assert rnx.channel_head_chi_split(chi, z, min_segment_length=m) == best


def test_channel_head_chi_split_too_short():
    assert rnx.channel_head_chi_split([3, 2, 1], [3, 2, 1], min_segment_length=10) == -1


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("PASS", name)
    print("All core tests passed.")
