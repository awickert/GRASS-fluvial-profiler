"""Diagnose the order-0 'cycle' segments in the MBR divide network: are they
tiny enclosed-basin noise loops (the first-order signal to examine), or an
artifact of the -1 drains-nowhere label regions / segment tracing?

Caches `lab` to /tmp so reruns skip the slow drainage_divides.
"""
import os
import numpy as np
from rivernetworkx import dreich as D
from rivernetworkx import divides as DV
import divides_lib as L

STREAM_T = 10000
LABCACHE = '/tmp/divides_lab_T%d.npy' % STREAM_T
C = L.load_cache()
fi = C['fi']
if os.path.exists(LABCACHE):
    lab = np.load(LABCACHE)
else:
    _, lab = D.drainage_divides(fi, threshold=STREAM_T)
    np.save(LABCACHE, lab)

valid = lab >= 0
print('cells: total=%d  label>=0=%d (%.1f%%)  label==-1=%d (%.1f%%)'
      % (lab.size, valid.sum(), 100 * valid.mean(),
         (lab == -1).sum(), 100 * (lab == -1).mean()))
print('distinct basin labels: %d' % len(np.unique(lab[valid])))

edges, meta = DV.extract_divide_edges(lab)
g = DV.build_divide_graph(edges)
segs = g['segments']
seglen = np.array([len(p) - 1 for p in segs])
# classify: cycle if first==last terminus OR a terminus not a break node
breaks = set(g['endpoints']) | set(g['junctions'])
is_cycle = np.array([(p[0] == p[-1]) or (p[0] not in breaks) or (p[-1] not in breaks)
                     for p in segs])
print('\nsegments: %d total, %d cycles (%.0f%%)'
      % (len(segs), is_cycle.sum(), 100 * is_cycle.mean()))
print('cycle seg length: med=%.0f p90=%.0f max=%d'
      % (np.median(seglen[is_cycle]), np.percentile(seglen[is_cycle], 90),
         seglen[is_cycle].max()))
print('non-cycle seg length: med=%.0f p90=%.0f'
      % (np.median(seglen[~is_cycle]), np.percentile(seglen[~is_cycle], 90)))

# how many cycle segments border a -1 region? (decode an edge of each, check cells)
cw = meta['cw']
def borders_nodata(p):
    for k in range(min(len(p) - 1, 5)):
        i0, j0 = p[k] // cw, p[k] % cw
        for di in (-1, 0):
            for dj in (-1, 0):
                ii, jj = i0 + di, j0 + dj
                if 0 <= ii < meta['nr'] and 0 <= jj < meta['nc'] and lab[ii, jj] == -1:
                    return True
    return False
cyc_idx = np.where(is_cycle)[0]
samp = cyc_idx[:: max(1, len(cyc_idx) // 500)]
nd_frac = np.mean([borders_nodata(segs[i]) for i in samp])
print('\nfraction of (sampled) cycle segs bordering a -1 region: %.2f' % nd_frac)
