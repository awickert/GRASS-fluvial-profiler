"""Verify the productionised divide-based head method
(extract_channel_heads(valleys='divides')) reproduces the validated v2 on 1 m
Mid Bailey Run: ~943 heads at the valley scale T=10000, with the same Clubb
field-head recall. If this matches, the curvature gate has been cleanly replaced
by divide-defined valleys, reusing DrEICH's chi-z stage.

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 test_method_divides.py
"""
import pickle
import time
import numpy as np
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
import divides_lib as L

WEST, NORTH = L.WEST, L.NORTH
T0 = time.time()
print('loading cache...', flush=True)
d = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))
filled = d['filled']
clubb = L.load_clubb()

for T in (5000, 10000):
    t = time.time()
    heads = D.extract_channel_heads(filled, nodata=-9999.0, cellsize=1.0,
                                    fill_dem=False, filled=filled,
                                    threshold=T, valleys='divides')
    hx = np.array([WEST + c + 0.5 for (r, c) in heads])
    hy = np.array([NORTH - r - 0.5 for (r, c) in heads])
    tt = cKDTree(np.column_stack([hx, hy]))
    dist, _ = tt.query(clubb)
    rec = {tol: float(np.mean(dist <= tol)) for tol in (10, 30, 50, 100)}
    print('T=%-6d  %5d heads  Clubb recall @30=%.2f @50=%.2f @100=%.2f  median=%.0fm  (%.0fs)'
          % (T, len(heads), rec[30], rec[50], rec[100], np.median(dist), time.time() - t),
          flush=True)
print('total %.0fs' % (time.time() - T0))
print('(v2 reference at T=10000: 943 heads, Clubb recall@50 ~0.57, median ~46 m)')
