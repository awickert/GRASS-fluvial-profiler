"""Design B -- chi-z self-gating, validated against Clubb's field heads.

Design A (divides_valleys_v2): valley set = T-sources; a single global threshold
T sets which valleys exist; keep all. Design B: OVER-generate valleys at a fine T,
score each by the chi-z split quality, and keep only score >= tau. The question:
does the chi-z linearity score select the field heads better than a global T?

Candidates (head, score, area) are computed once at a fine T and cached, so the
tau sweep is free. Compared to Clubb's 53 field heads (recall + median distance),
head-to-head with Design A at matched count.

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_designB.py
"""
import os
import time
import numpy as np
from rivernetworkx import dreich as D
import divides_lib as L

CAND_T = 2000               # fine candidate scale (over-complete valley set)
CANDCACHE = '/tmp/divides_designB_cand.npz'

T0 = time.time()
def log(m): print('[%6.1fs] %s' % (time.time() - T0, m), flush=True)

C = L.load_cache()
fi, filled = C['fi'], C['filled']
clubb = L.load_clubb()
log('loaded; %d Clubb field heads' % len(clubb))

if os.path.exists(CANDCACHE):
    z = np.load(CANDCACHE)
    cand_head, cand_score = z['head'], z['score']
    log('loaded %d cached candidates (T=%d)' % (len(cand_head), CAND_T))
else:
    tsrc = D.get_sources(fi, CAND_T)
    log('scoring %d candidate valleys (T=%d)...' % (len(tsrc), CAND_T))
    heads = np.empty(len(tsrc), dtype=np.int64)
    scores = np.empty(len(tsrc), dtype=np.float64)
    for i, s in enumerate(tsrc):
        h, sc, _n = L.head_and_score(fi, filled, int(s), min_segment_length=10)
        heads[i] = h; scores[i] = sc
    cand_head, cand_score = heads, scores
    np.savez(CANDCACHE, head=heads, score=scores)
    log('scored + cached %d candidates' % len(heads))

valid = cand_head != 0
cand_head, cand_score = cand_head[valid], cand_score[valid]
finite = np.isfinite(cand_score)
print('score distribution (valid cands): min=%.2f p25=%.2f med=%.2f p75=%.2f max=%.2f'
      % (cand_score[finite].min(), np.percentile(cand_score[finite], 25),
         np.median(cand_score[finite]), np.percentile(cand_score[finite], 75),
         cand_score[finite].max()))

# -------- Design B: keep score >= tau, dedup, compare to Clubb
print('\nDESIGN B (chi-z score gate, candidate T=%d):' % CAND_T)
print('  tau    heads   r@30  r@50  r@100   med_dist')
for tau in (0.0, 0.5, 0.8, 0.9, 0.95, 0.99):
    keep = cand_head[cand_score >= tau]
    final = L.dedup_furthest_upstream(fi, filled, keep)
    final = [h for h in final if h != 0]
    rec, med = L.recall_vs_points(L.nodes_xy(fi, final), clubb)
    print('  %.2f  %6d   %.2f  %.2f  %.2f   %5.0fm'
          % (tau, len(final), rec[30], rec[50], rec[100], med), flush=True)

# -------- Design A (v2): single global T, keep all -- for the head-to-head
print('\nDESIGN A (global T threshold, keep all):')
print('     T   heads   r@30  r@50  r@100   med_dist')
for T in (5000, 8000, 10000, 15000):
    heads = []
    for s in D.get_sources(fi, T):
        h, _s, _n = L.head_and_score(fi, filled, int(s), min_segment_length=10)
        heads.append(h)
    final = [h for h in L.dedup_furthest_upstream(fi, filled, heads) if h != 0]
    rec, med = L.recall_vs_points(L.nodes_xy(fi, final), clubb)
    print('  %6d  %5d   %.2f  %.2f  %.2f   %5.0fm'
          % (T, len(final), rec[30], rec[50], rec[100], med), flush=True)
log('done')
