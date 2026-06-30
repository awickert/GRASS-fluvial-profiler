"""Figure for the valley head length scale (run valley_length_scale.py first).
Panels: (1) per-head linear half-width + flowpath length distributions with the
independent 1/(2Dd) anchor; (2) heads mapped, coloured by linear half-width;
(3) one example head with its contributing-area rim, showing the head sits at the
downstream TIP of its bowl (why rim distances run 0 -> headwall)."""
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage
from rivernetworkx import dreich as D
import divides_lib as L

RES = 1.0
d = np.load('/tmp/valley_length_scale_1m.npz', allow_pickle=True)
cols = list(d['columns']); T = d['table']
def col(n): return T[:, cols.index(n)]
Lh_Dd = float(d['Lh_from_Dd'])
ml, mf = col('lin_median'), col('flow_median')

fig = plt.figure(figsize=(15, 5))

# --- panel 1: distributions ---------------------------------------------------
ax = fig.add_subplot(1, 3, 1)
bins = np.linspace(0, 180, 46)
ax.hist(ml, bins=bins, alpha=0.6, label='linear half-width', color='#3b6')
ax.hist(mf, bins=bins, alpha=0.5, label='flowpath length', color='#36b')
ax.axvline(np.median(ml), color='#3b6', ls='--', lw=1.5)
ax.axvline(np.median(mf), color='#36b', ls='--', lw=1.5)
ax.axvline(Lh_Dd, color='k', ls=':', lw=2, label='1/(2Dd) anchor = %.0f m' % Lh_Dd)
ax.set_xlabel('per-head distance (m)'); ax.set_ylabel('number of heads')
ax.set_title('valley head length scale (1 m MBR, %d heads)\nlinear med %.0f m, flowpath med %.0f m'
             % (len(ml), np.median(ml), np.median(mf)))
ax.legend(fontsize=8)

# --- panel 2: map -------------------------------------------------------------
ax = fig.add_subplot(1, 3, 2)
sc = ax.scatter(col('x'), col('y'), c=ml, s=14, cmap='viridis', vmin=40, vmax=120)
ax.set_aspect('equal'); ax.set_title('heads coloured by linear half-width (m)')
ax.set_xlabel('Easting (m)'); ax.set_ylabel('Northing (m)')
plt.colorbar(sc, ax=ax, shrink=0.8)

# --- panel 3: one example bowl ------------------------------------------------
ax = fig.add_subplot(1, 3, 3)
filled = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))['filled']
heads, fi = D.extract_channel_heads(filled, nodata=-9999.0, cellsize=RES,
                                    fill_dem=False, filled=filled,
                                    threshold=10000, valleys='divides',
                                    return_flowinfo=True)
NI, SV, SVI, nc = fi['NodeIndex'], fi['SVector'], fi['SVectorIndex'], fi['ncontrib']
ro, co = fi['row_of'], fi['col_of']
# pick a head near the median area
areas = np.array([nc[NI[r, c]] for (r, c) in heads])
hi = int(np.argsort(np.abs(areas - np.median(areas)))[0])
hr, hc = heads[hi]; h = NI[hr, hc]
U = SV[int(SVI[h]):int(SVI[h]) + int(nc[h])]
rows, cols2 = ro[U], co[U]
r0, c0 = rows.min() - 1, cols2.min() - 1
H = rows.max() - r0 + 2; W = cols2.max() - c0 + 2
mask = np.zeros((H, W), bool); mask[rows - r0, cols2 - c0] = True
er = ndimage.binary_erosion(mask, structure=np.ones((3, 3), bool), border_value=0)
rim = mask & ~er
ax.imshow(np.where(mask, 1.0, np.nan), cmap='Greys', alpha=0.35,
          extent=[c0, c0 + W, r0 + H, r0], interpolation='none')
rr, cc = np.where(rim)
ax.scatter(cc + c0 + 0.5, rr + r0 + 0.5, s=4, c='#c33', label='head-bounding rim')
ax.scatter([hc + 0.5], [hr + 0.5], s=80, marker='*', c='gold',
           edgecolor='k', zorder=5, label='channel head')
ax.set_aspect('equal'); ax.invert_yaxis()
ax.set_title('example head + contributing-area rim\n(head at downstream tip; area %d m$^2$)' % nc[h])
ax.set_xlabel('col'); ax.set_ylabel('row'); ax.legend(fontsize=8, loc='lower right')

plt.tight_layout()
plt.savefig('dev/divide_heads/figures/valley_length_scale.png', dpi=130)
print('wrote dev/divide_heads/figures/valley_length_scale.png')
