"""Decisive 1 m check: does the divides->chi-z method find Clubb's 53 field-mapped
channel heads as well as curvature-DrEICH does? If yes, the modest divides-vs-
curvature agreement (~0.5) is moot -- both equally approximate ground truth, and
divides add resolution robustness.

Recall = fraction of Clubb's 53 heads with a method head within tol. (Precision
vs Clubb is not meaningful -- Clubb mapped only part of the basin, so most method
heads lie in unmapped area.)

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_vs_clubb.py
"""
import numpy as np
import openpyxl
from scipy.spatial import cKDTree
from rivernetworkx import dreich as D
import divides_lib as L

WEST, NORTH = 398767.32685636, 4369656.1523925

# Clubb's Mid Bailey Run field heads (projected E/N, same CRS as the DEM)
wb = openpyxl.load_workbook('/tmp/clubb_channel_heads.xlsx', data_only=True, read_only=True)
ws = wb['Sheet1']
hE, hN = [], []
for row in ws.iter_rows(min_row=3, values_only=True):
    if row[0] and 'Bailey' in str(row[0]):
        hN.append(float(row[2])); hE.append(float(row[3]))
clubb = np.column_stack([hE, hN])
print('Clubb field heads: %d' % len(clubb))

C = L.load_cache()
fi, filled = C['fi'], C['filled']


def xy(nodes):
    nodes = np.asarray(nodes, dtype=np.int64)
    return np.column_stack([WEST + (fi['col_of'][nodes] + 0.5),
                            NORTH - (fi['row_of'][nodes] + 0.5)])


def divide_heads_1m(T):
    heads = []
    for s in D.get_sources(fi, T):
        h, _s, _n = L.head_and_score(fi, filled, int(s), min_segment_length=10)
        heads.append(h)
    final = L.dedup_furthest_upstream(fi, filled, heads)
    return [h for h in final if h != 0]


methods = {'curvature(634)': np.asarray(C['ref_heads'])}
for T in (5000, 10000, 20000):
    methods['divides T=%d' % T] = np.asarray(divide_heads_1m(T))

ck = cKDTree(clubb)
print('\nmethod              heads   recall vs Clubb @  10m  30m  50m 100m   med_dist')
for name, nodes in methods.items():
    mxy = xy(nodes)
    mt = cKDTree(mxy)
    d_clubb, _ = mt.query(clubb)          # nearest method head to each Clubb head
    rec = [np.mean(d_clubb <= t) for t in (10, 30, 50, 100)]
    print('%-18s  %5d        %.2f %.2f %.2f %.2f      %5.0fm'
          % (name, len(nodes), rec[0], rec[1], rec[2], rec[3], np.median(d_clubb)))
