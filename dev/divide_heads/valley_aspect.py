"""Length, width, and aspect ratio of first-order valleys (1 m Mid Bailey Run),
defined WITHOUT using catchment area (per Andy's objection that width != area/L).

For each divide-defined channel head:
  * rim       = boundary of the head's contributing area (head-bounding arc)
  * axis      = head -> hilltop (farthest-upslope rim point): the along-valley
                axis, i.e. "valley head down to channel head" (LENGTH direction)
  * LENGTH    = extent of the rim projected ON the axis
  * WIDTH     = extent of the rim projected PERPENDICULAR to the axis
  * aspect    = length / width
PCA on the rim gives an independent aspect-ratio cross-check (scale-free).

Question answered: is the aspect ratio ~constant across valleys (-> one length
scale, everything scales linearly, "no damage") or variable (-> length and width
are genuinely independent and both matter)?

Run: PYTHONPATH=<repo>:<repo>/dev /usr/bin/python3 valley_aspect.py
"""
import pickle
import numpy as np
from scipy import ndimage
from rivernetworkx import dreich as D

RES = 1.0
filled = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))['filled']
print('extracting heads + flowinfo (T=10000)...', flush=True)
heads, fi = D.extract_channel_heads(filled, nodata=-9999.0, cellsize=RES,
                                    fill_dem=False, filled=filled,
                                    threshold=10000, valleys='divides',
                                    return_flowinfo=True)
assert 900 <= len(heads) <= 980, 'anchor failed: %d heads' % len(heads)
print('  anchor OK: %d heads\n' % len(heads), flush=True)

NI, SV, SVI, ncon = fi['NodeIndex'], fi['SVector'], fi['SVectorIndex'], fi['ncontrib']
ro, co, fd = fi['row_of'], fi['col_of'], fi['fd']
STR8 = np.ones((3, 3), bool)


def rim_nodes(h):
    s = int(SVI[h]); n = int(ncon[h])
    U = SV[s:s + n]; rows = ro[U]; cols = co[U]
    r0, c0 = rows.min(), cols.min()
    H = rows.max() - r0 + 3; Wd = cols.max() - c0 + 3
    lr = rows - r0 + 1; lc = cols - c0 + 1
    mask = np.zeros((H, Wd), bool); mask[lr, lc] = True
    locn = np.full((H, Wd), -1, np.int64); locn[lr, lc] = U
    rim = mask & ~ndimage.binary_erosion(mask, structure=STR8, border_value=0)
    return locn[rim], n


rows = []
for (hr, hc) in heads:
    h = NI[hr, hc]
    rn, n = rim_nodes(int(h))
    if rn.size < 3:
        continue
    rr = ro[rn].astype(float); cc = co[rn].astype(float)
    flow = np.clip(fd[rn] - fd[h], 0, None)
    k = int(np.argmax(flow))                       # hilltop = farthest-upslope rim cell
    axis = np.array([cc[k] - hc, rr[k] - hr], float)
    nrm = np.hypot(*axis)
    if nrm < 1e-9:
        continue
    u = axis / nrm; up = np.array([-u[1], u[0]])
    P = np.column_stack([cc - hc, rr - hr])        # rim relative to head (cells)
    along = P @ u; perp = P @ up
    length = (along.max() - along.min()) * RES
    width = (perp.max() - perp.min()) * RES
    # independent aspect via PCA on the rim point cloud
    Q = np.column_stack([cc, rr]); Q = Q - Q.mean(0)
    ev = np.linalg.eigvalsh(Q.T @ Q / Q.shape[0])
    pca_aspect = float(np.sqrt(max(ev[1], 1e-9) / max(ev[0], 1e-9)))
    rows.append((length, width, length / width, pca_aspect,
                 float(flow[k]) * RES, n * RES * RES))

A = np.array(rows)
length, width, aspect, pca_aspect, flowlen, area = A.T

def q(a, p): return np.percentile(a, p)
print('=== first-order valley length & width (1 m MBR, %d heads) ===' % len(A))
print('LENGTH (axis, head->hilltop)   median %.0f m   IQR %.0f-%.0f' % (np.median(length), q(length,25), q(length,75)))
print('WIDTH  (perp to axis)          median %.0f m   IQR %.0f-%.0f' % (np.median(width), q(width,25), q(width,75)))
print('ASPECT length/width            median %.2f   IQR %.2f-%.2f   p10 %.2f p90 %.2f'
      % (np.median(aspect), q(aspect,25), q(aspect,75), q(aspect,10), q(aspect,90)))
print('ASPECT (PCA, independent)      median %.2f   IQR %.2f-%.2f   <- cross-check'
      % (np.median(pca_aspect), q(pca_aspect,25), q(pca_aspect,75)))
print()
cv = lambda a: np.std(a) / np.mean(a)
print('--- is there one length scale, or two? ---')
print('CV(aspect ratio)               %.2f   (small -> ~constant aspect -> one scale; large -> two)' % cv(aspect))
r = np.corrcoef(np.log(length), np.log(width))[0, 1]
print('corr(log length, log width)    %.2f   (->1 = they co-vary = one scale)' % r)
# does width predicted from a single scale (constant aspect) hold?
const_aspect = np.median(aspect)
pred_w = length / const_aspect
rel_err = np.abs(pred_w - width) / width
print('width from length*const-aspect: median rel. error %.0f%%' % (100 * np.median(rel_err)))
print()
print('compare to the BAD area/length width: %.0f m  vs proper projection width %.0f m'
      % (np.median(area / flowlen), np.median(width)))

np.savez('/tmp/valley_aspect_1m.npz',
         length=length, width=width, aspect=aspect, pca_aspect=pca_aspect,
         flowlen=flowlen, area=area)
print('\nper-head length/width/aspect -> /tmp/valley_aspect_1m.npz')
