"""Adaptive T-finder figure (run adaptive_t.py <site> first).
(1) n_src and n_chan vs T (log-log); (2) their log-log slopes vs T -- a constant
slope = scale-invariant 'settled' regime, a break/knee = transition. Field-
validated T marked: does a feature land there with no head data used?"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

keys = sys.argv[1:] or ['mbr']
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
for key in keys:
    d = np.load('/tmp/adaptive_t_%s.npz' % key, allow_pickle=True)
    T = d['T']; tf = float(d['T_field']); site = str(d['site'])
    ax[0].loglog(T, d['nsrc'], 'o-', ms=3, label='%s n_src' % site)
    ax[0].loglog(T, d['nchan'], 's--', ms=3, alpha=0.5, label='%s n_chan' % site)
    ax[1].semilogx(T, d['s_src'], 'o-', ms=3, label='%s d lnN_src/d lnT' % site)
for a in ax:
    a.axvline(tf, color='k', ls=':', label='field T~%g' % tf)
ax[0].set_xlabel('T (cells)'); ax[0].set_ylabel('count'); ax[0].legend(fontsize=7)
ax[0].set_title('network counts vs threshold')
ax[1].set_xlabel('T (cells)'); ax[1].set_ylabel('d ln(count) / d ln T')
ax[1].set_title('log-log slope (knee/plateau = transition)'); ax[1].legend(fontsize=7)
ax[1].axhline(0, color='grey', lw=0.5)
plt.tight_layout()
plt.savefig('dev/divide_heads/figures/adaptive_t.png', dpi=130)
print('wrote dev/divide_heads/figures/adaptive_t.png')
