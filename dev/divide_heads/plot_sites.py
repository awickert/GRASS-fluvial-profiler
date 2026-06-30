"""Cross-site Arm-A figure: density-controlled field-head recall/precision vs
detected head density, MBR vs Indian Creek, with each site's field density and
the interpolated density-matched operating point marked.
Run site_experiment.py first (writes /tmp/site_<key>.npz)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SITES = [('mbr', 'Mid Bailey Run', 'C0'), ('indian', 'Indian Creek', 'C1')]
# rows cols: 0 T_m2, 1 Tcell, 2 nHd, 3 density, 4 rec30, 5 rec50, 6 prec50, 7 median
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
for key, name, c in SITES:
    d = np.load('/tmp/site_%s.npz' % key, allow_pickle=True)
    A = d['rows']; dens = A[:, 3]; o = np.argsort(dens)
    fd = float(d['field_dens']); m = d['matched']     # (r30, r50, prec, median)
    ax.plot(dens[o], A[o, 5], 'o-', color=c, label='%s recall@50' % name)
    ax.plot(dens[o], A[o, 6], 's--', color=c, alpha=0.55, label='%s precision@50' % name)
    ax.axvline(fd, color=c, ls=':', lw=1.2)
    ax.plot([fd], [m[1]], '*', color=c, ms=16, mec='k', zorder=6)
    ax.plot([fd], [m[2]], '*', color=c, ms=12, mec='k', alpha=0.6, zorder=6)
    ax2.plot(dens[o], A[o, 7], 'o-', color=c, label='%s median dist' % name)
    ax2.axvline(fd, color=c, ls=':', lw=1.2)
    ax2.plot([fd], [m[3]], '*', color=c, ms=16, mec='k', zorder=6)
ax.set_xlabel('detected head density (heads / km$^2$)'); ax.set_ylabel('recall / precision @ 50 m')
ax.set_ylim(0, 1.02); ax.set_title('Arm A: field-head recovery vs density\n(stars = interpolated field-density operating point)')
ax.legend(fontsize=8)
ax2.set_xlabel('detected head density (heads / km$^2$)'); ax2.set_ylabel('median nearest distance (m)')
ax2.set_title('median field-head distance vs density'); ax2.legend(fontsize=8)
plt.tight_layout()
plt.savefig('dev/divide_heads/figures/sites_armA.png', dpi=130)
print('wrote dev/divide_heads/figures/sites_armA.png')
