"""Figure: first-order valley aspect ratio (run valley_aspect.py first).
(1) aspect-ratio distribution; (2) length-vs-width with the constant-aspect line;
(3) one example valley with its head->hilltop axis and perpendicular width."""
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import ndimage
from rivernetworkx import dreich as D

RES = 1.0
d = np.load('/tmp/valley_aspect_1m.npz')
length, width, aspect = d['length'], d['width'], d['aspect']

fig = plt.figure(figsize=(15, 5))

ax = fig.add_subplot(1, 3, 1)
ax.hist(aspect, bins=np.linspace(1, 6, 41), color='#759', alpha=0.8)
ax.axvline(np.median(aspect), color='k', ls='--', label='median %.2f' % np.median(aspect))
ax.set_xlabel('aspect ratio  length / width'); ax.set_ylabel('number of valleys')
ax.set_title('valley aspect ratio (1 m MBR, %d valleys)\nIQR %.2f-%.2f, long tail = elongated hollows'
             % (len(aspect), np.percentile(aspect, 25), np.percentile(aspect, 75)))
ax.legend()

ax = fig.add_subplot(1, 3, 2)
ax.scatter(width, length, s=8, alpha=0.4, color='#367')
xx = np.linspace(width.min(), width.max(), 50)
ax.plot(xx, np.median(aspect) * xx, 'k--', label='constant aspect %.2f' % np.median(aspect))
ax.set_xlabel('width (m)'); ax.set_ylabel('length (m)')
ax.set_title('length vs width\ncorr(log,log) %.2f  -> co-vary but scatter is real'
             % np.corrcoef(np.log(length), np.log(width))[0, 1])
ax.legend()

# example valley with axis + width
ax = fig.add_subplot(1, 3, 3)
filled = pickle.load(open('/tmp/divides_proto_fi.pkl', 'rb'))['filled']
heads, fi = D.extract_channel_heads(filled, nodata=-9999.0, cellsize=RES, fill_dem=False,
                                    filled=filled, threshold=10000, valleys='divides',
                                    return_flowinfo=True)
NI, SV, SVI, nc = fi['NodeIndex'], fi['SVector'], fi['SVectorIndex'], fi['ncontrib']
ro, co, fd = fi['row_of'], fi['col_of'], fi['fd']
areas = np.array([nc[NI[r, c]] for (r, c) in heads])
hi = int(np.argsort(np.abs(areas - np.median(areas)))[0])
hr, hc = heads[hi]; h = NI[hr, hc]
U = SV[int(SVI[h]):int(SVI[h]) + int(nc[h])]; rows, cols = ro[U], co[U]
r0, c0 = rows.min() - 1, cols.min() - 1
H = rows.max() - r0 + 2; W = cols.max() - c0 + 2
mask = np.zeros((H, W), bool); mask[rows - r0, cols - c0] = True
rim = mask & ~ndimage.binary_erosion(mask, structure=np.ones((3, 3), bool), border_value=0)
ax.imshow(np.where(mask, 1.0, np.nan), cmap='Greys', alpha=0.3,
          extent=[c0, c0 + W, r0 + H, r0], interpolation='none')
rr, cc = np.where(rim); rimc = cc + c0 + 0.5; rimr = rr + r0 + 0.5
flow = np.clip(fd[NI[rr + r0, cc + c0]] - fd[h], 0, None)
k = int(np.argmax(flow))
ax.scatter(rimc, rimr, s=4, c='#c33')
ax.scatter([hc + 0.5], [hr + 0.5], s=120, marker='*', c='gold', edgecolor='k', zorder=6, label='channel head')
ax.scatter([rimc[k]], [rimr[k]], s=70, marker='^', c='#3a3', edgecolor='k', zorder=6, label='hilltop (valley head)')
ax.annotate('', xy=(rimc[k], rimr[k]), xytext=(hc + 0.5, hr + 0.5),
            arrowprops=dict(arrowstyle='->', color='k', lw=2))
ax.text(0.05, 0.95, 'length axis\n(head->hilltop)', transform=ax.transAxes, va='top', fontsize=8)
ax.set_aspect('equal'); ax.invert_yaxis()
ax.set_title('example valley: length axis + rim\n(width = rim extent perpendicular to axis)')
ax.set_xlabel('col'); ax.set_ylabel('row'); ax.legend(fontsize=8, loc='lower right')

plt.tight_layout()
plt.savefig('dev/divide_heads/figures/valley_aspect.png', dpi=130)
print('wrote dev/divide_heads/figures/valley_aspect.png')
