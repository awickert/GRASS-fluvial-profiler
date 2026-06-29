"""Design A -- divide-shape valley gate (1 m Mid Bailey Run).

Replace the DrEICH cross-valley-curvature gate (find_valleys) with a divide-flank
gate: a cell is valley-like if flanked by drainage divides on both sides within W.
Walk each source down the divide-flanked valley bottom exactly as find_valleys
does, tag the valley junction, and run the UNCHANGED chi-z head stage. Compare the
resulting heads to the 634 curvature-DrEICH heads (the 1 m reference Andy chose).

Run:  PYTHONPATH=. /usr/bin/python3 dev/divides_designA.py
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
log('loaded cache: %d sources, reference = %d curvature heads'
    % (len(C['sources']), len(ref_heads)))

sources = C['sources']
# DrEICH ridge scale: basins of ~first-order-valley size (the ridge-survival test
# used 50000-cell basins -> 1.8% divide density). Sweep around it.
print('\n  Dthr     W   div%%  flank%%  valleys  heads   R@5  R@10  R@30  P@5  P@10  P@30')
for DIV_THRESHOLD in (20000, 50000, 100000):
    divide, _ = D.drainage_divides(fi, threshold=DIV_THRESHOLD)
    divpct = 100.0 * divide.sum() / divide.size
    log('drainage_divides(%d): %.2f%% divide cells' % (DIV_THRESHOLD, divpct))
    for W in (80, 120, 200):
        grid = L.valley_flank_grid(divide, W)
        gate = L.gate_from_grid(fi, grid)
        flankpct = 100.0 * grid.sum() / grid.size
        valley = L.find_valleys_gated(fi, gate, sources, n_connecting_nodes=10)
        nv = len(np.unique(valley[valley != L.ND]))
        _, final = D.channel_heads_from_valleys(fi, filled, valley,
                                                min_segment_length=10,
                                                A_0=1000.0, m_over_n=0.525)
        test_rc = L.heads_rc(fi, final)
        cmp = L.compare_heads(test_rc, ref_rc, res=1.0, tols=(5, 10, 30))
        print('%6d  %4d  %5.2f  %6.2f  %7d  %5d  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f'
              % (DIV_THRESHOLD, W, divpct, flankpct, nv, len(final),
                 cmp[5][0], cmp[10][0], cmp[30][0],
                 cmp[5][1], cmp[10][1], cmp[30][1]), flush=True)
log('done')
