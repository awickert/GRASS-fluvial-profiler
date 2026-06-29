"""Divides->heads, correct realization: the divides define first-order VALLEYS
directly. The channel network at a valley-scale area threshold T partitions the
terrain into divide-bounded first-order basins; each basin's T-source is the
valley, and one chi-z head is found on the ridgetop->T-source profile.

T sets the valley scale (robust to coarsening -- ridge survival 0.84+); the chi-z
fit localizes the head within each valley. No curvature anywhere.

Compare to the 634 curvature-DrEICH heads (1 m reference). Sweep T.

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_valleys_v2.py
"""
import time
import numpy as np
from rivernetworkx import dreich as D
import divides_lib as L

T0 = time.time()
def log(m): print('[%6.1fs] %s' % (time.time() - T0, m), flush=True)

C = L.load_cache()
fi, filled, ref_heads = C['fi'], C['filled'], C['ref_heads']
ref_rc = L.heads_rc(fi, ref_heads)
log('loaded; reference = %d heads, med area = %d cells'
    % (len(ref_heads), np.median(fi['ncontrib'][ref_heads])))

print('\n     T   valleys  heads  hd_med_area   R@10  R@30  R@50  P@10  P@30  P@50')
for T in (5000, 10000, 20000, 50000, 100000):
    tsrc = D.get_sources(fi, T)               # first-order valley heads at scale T
    heads = []
    for s in tsrc:
        h, _score, _n = L.head_and_score(fi, filled, int(s),
                                         A_0=1000.0, m_over_n=0.525,
                                         min_segment_length=10)
        heads.append(h)
    final = L.dedup_furthest_upstream(fi, filled, heads)
    final = np.asarray([h for h in final if h != 0], dtype=np.int64)
    test_rc = L.heads_rc(fi, final)
    cmp = L.compare_heads(test_rc, ref_rc, res=1.0, tols=(10, 30, 50))
    hd_area = int(np.median(fi['ncontrib'][final])) if len(final) else -1
    print('%6d  %7d  %5d  %11d   %.2f  %.2f  %.2f  %.2f  %.2f  %.2f'
          % (T, len(tsrc), len(final), hd_area,
             cmp[10][0], cmp[30][0], cmp[50][0],
             cmp[10][1], cmp[30][1], cmp[50][1]), flush=True)
    log('T=%d done' % T)
log('done')
