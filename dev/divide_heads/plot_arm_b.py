"""Arm B figures (run arm_b.py <site> first -> /tmp/arm_b_<site>.npz).
(1) recall vs absolute resolution (divides, curvature, valley structure) + 2R mark;
(2) divides recall@50 vs dimensionless res/(2R_site);
(3) per-valley breakdown: detection fraction vs res/(2R_i) (per-head local radius)
    vs the site-aggregate -- does it sharpen to a cliff at 1 (size variation) or
    stay gradual (graceful method)?"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
d = np.load('/tmp/arm_b_%s.npz' % key, allow_pickle=True)
rows = d['rows']; det = d['det']; fld_R = d['fld_R']; R_site = float(d['R_site'])
res = rows[:, 0]; site = str(d['site'])
twoR = 2 * R_site

fig, ax = plt.subplots(1, 3, figsize=(16, 5))

# (1) absolute resolution
ax[0].plot(res, rows[:, 3], 'o-', color='C0', label='divides recall@50')
ax[0].plot(res, rows[:, 2], 'o--', color='C0', alpha=0.5, label='divides recall@30')
ax[0].plot(res, rows[:, 5], 's-', color='C3', label='curvature recall@50')
ax[0].plot(res, rows[:, 6], '^-', color='C2', label='valley-structure recall')
ax[0].axvline(twoR, color='k', ls=':', label='2R = %.0f m (predicted breakdown)' % twoR)
ax[0].axhline(0.5, color='grey', lw=0.5)
ax[0].set_xlabel('DEM resolution (m)'); ax[0].set_ylabel('recall'); ax[0].set_ylim(0, 1.05)
ax[0].set_title('%s: Arm B vs resolution' % site); ax[0].legend(fontsize=7)

# (2) dimensionless res/(2R_site)
x = res / twoR
ax[1].plot(x, rows[:, 3], 'o-', color='C0', label='divides recall@50')
ax[1].plot(x, rows[:, 6], '^-', color='C2', label='valley-structure recall')
ax[1].axvline(1.0, color='k', ls=':', label='res = 2R (one pixel = valley-head width)')
ax[1].axhline(0.5, color='grey', lw=0.5)
ax[1].set_xlabel('res / 2R (pixels per valley-head width)'); ax[1].set_ylabel('recall'); ax[1].set_ylim(0, 1.05)
ax[1].set_xscale('log'); ax[1].set_title('dimensionless (site median R)'); ax[1].legend(fontsize=8)

# (3) per-valley normalized breakdown vs aggregate
ok = np.isfinite(fld_R) & (fld_R > 0)
xi = (res[:, None] / (2 * fld_R[None, :]))[:, ok].ravel()      # res/(2R_i) per (res, head)
di = det[:, ok].astype(float).ravel()
bins = np.logspace(np.log10(max(xi.min(), 1e-2)), np.log10(xi.max()), 16)
idx = np.digitize(xi, bins)
bx = [0.5 * (bins[i - 1] + bins[i]) for i in range(1, len(bins))]
by = [di[idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, len(bins))]
ax[2].plot(bx, by, 'o-', color='C4', label='per-valley (res/2R$_i$, own radius)')
ax[2].plot(x, rows[:, 3], 's-', color='C0', alpha=0.6, label='aggregate (res/2R$_{site}$)')
ax[2].axvline(1.0, color='k', ls=':', label='res = 2R$_i$')
ax[2].axhline(0.5, color='grey', lw=0.5)
ax[2].set_xlabel('res / 2R$_i$'); ax[2].set_ylabel('field-head detection fraction'); ax[2].set_ylim(0, 1.05)
ax[2].set_xscale('log'); ax[2].set_title('graceful method vs valley-size distribution'); ax[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig('dev/divide_heads/figures/arm_b_%s.png' % key, dpi=130)
print('wrote dev/divide_heads/figures/arm_b_%s.png' % key)
