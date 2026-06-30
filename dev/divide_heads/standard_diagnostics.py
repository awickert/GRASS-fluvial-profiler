"""Standard field-free channel/valley-initiation diagnostics (Montgomery &
Foufoula-Georgiou 1993): slope-area rollover and the cumulative area distribution
break -- from the flow routing alone, ONE network pass, no T-sweep. Compares the
break to the field-validated valley scale T~5000 m2.

  * slope-area: median local slope in log-area bins. Hillslopes: slope ~flat or
    rising with area; valleys/channels: slope falls (S ~ A^-theta). The ROLLOVER
    (peak median slope) marks the hillslope->valley transition = a field-free
    valley-initiation area.
  * cumulative area distribution P(A>=a): a break in its log-log slope marks the
    same transition (channelized self-similar tail vs hillslopes).

Run: PYTHONPATH=<repo>:<repo>/dev:<repo>/dev/divide_heads /usr/bin/python3 standard_diagnostics.py <mbr|indian>
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rivernetworkx import dreich as D
from site_experiment import read_flt, load_field_heads, SITES

key = sys.argv[1] if len(sys.argv) > 1 else 'mbr'
cfg = SITES[key]; buf = 400; T_FIELD = 5000.0

z, west, north, res, nd = read_flt(cfg['flt']); z = np.where(z == nd, np.nan, z)
field = load_field_heads(cfg['key']); nr, nc = z.shape
cr = ((north - field[:, 1]) / res).astype(int); cc = ((field[:, 0] - west) / res).astype(int)
ins = (cr >= 0) & (cr < nr) & (cc >= 0) & (cc < nc)
R0, R1 = max(0, cr[ins].min() - buf), min(nr, cr[ins].max() + buf)
C0, C1 = max(0, cc[ins].min() - buf), min(nc, cc[ins].max() + buf)
zc = z[R0:R1, C0:C1]; dem = np.where(np.isfinite(zc), zc, nd).astype(np.float32)
print('%s: clip %s; routing...' % (cfg['name'], zc.shape), flush=True)
filled = D.fill(dem, nd, cellsize=res); fi = D.build_flowinfo(filled, nd, res); D.contributing_area(fi)
recv, flc, row, col, ncon = fi['recv'], fi['flc'], fi['row_of'], fi['col_of'], fi['ncontrib']

zf = filled[row, col].astype(np.float64)
step = np.where(flc == 2, res * np.sqrt(2.0), res).astype(np.float64)
slope = np.clip((zf - zf[recv]) / np.where(step > 0, step, 1), 1e-5, None)
area = ncon.astype(np.float64) * res * res
internal = recv != np.arange(fi['N'])                 # drop baselevels
slope, area = slope[internal], area[internal]

# slope-area: median slope in log-area bins (need enough cells per bin)
la = np.log10(area); bins = np.linspace(la.min(), la.max(), 45)
idx = np.digitize(la, bins); ba, bs, bn = [], [], []
for i in range(1, len(bins)):
    m = idx == i
    if m.sum() >= 30:
        ba.append(10 ** (0.5 * (bins[i - 1] + bins[i]))); bs.append(np.median(slope[m])); bn.append(m.sum())
ba, bs = np.array(ba), np.array(bs)
A_roll = ba[np.argmax(bs)]                              # plateau PEAK (not the transition)
# valley/channel initiation is the KNEE: where, past the hillslope plateau, the
# local slope-area exponent theta = -d ln S / d ln A rises into the fluvial range.
lA, lS = np.log(ba), np.log(bs)
theta = np.full_like(lA, np.nan); theta[1:-1] = -(lS[2:] - lS[:-2]) / (lA[2:] - lA[:-2])
ipk = int(np.argmax(bs)); THR = 0.2
knee_idx = next((i for i in range(ipk, len(ba)) if np.isfinite(theta[i]) and theta[i] >= THR), ipk)
A_knee = ba[knee_idx]

# cumulative area distribution P(A>=a) and its local log-log slope
agrid = np.logspace(np.log10(area.min()), np.log10(area.max()), 60)
Pge = np.array([(area >= a).mean() for a in agrid])
lp = np.log(np.maximum(Pge, 1e-9)); lg = np.log(agrid)
beta = np.full_like(lg, np.nan); beta[1:-1] = -(lp[2:] - lp[:-2]) / (lg[2:] - lg[:-2])

print('  slope-area plateau peak at A = %.0f m2 (NOT the transition)' % A_roll)
print('  slope-area fluvial-onset KNEE (theta>=%.1f) at A = %.0f m2  (field valley scale T~%g m2)'
      % (THR, A_knee, T_FIELD))
print('  cum-area exponent beta near A=5000 -> %.2f (fluvial ~0.4-0.45)'
      % np.interp(5000, agrid, beta))

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
ax[0].loglog(ba, bs, 'o-', color='C0')
ax[0].axvline(A_knee, color='C1', ls='-', lw=2, label='fluvial-onset knee %.0f m$^2$' % A_knee)
ax[0].axvline(A_roll, color='C1', ls='--', alpha=0.5, label='plateau peak %.0f m$^2$' % A_roll)
ax[0].axvline(T_FIELD, color='k', ls=':', label='field T~%g m$^2$' % T_FIELD)
ax[0].set_xlabel('drainage area (m$^2$)'); ax[0].set_ylabel('median slope')
ax[0].set_title('%s: slope-area (rollover = valley initiation)' % cfg['name']); ax[0].legend(fontsize=8)
ax[1].semilogx(agrid, beta, 'o-', color='C2', ms=3)
ax[1].axhline(0.43, color='grey', ls='--', label='fluvial beta~0.43')
ax[1].axvline(T_FIELD, color='k', ls=':', label='field T~%g' % T_FIELD)
ax[1].set_xlabel('drainage area a (m$^2$)'); ax[1].set_ylabel('cum-dist exponent d lnP(A>=a)/d ln a')
ax[1].set_title('cumulative area distribution exponent'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('dev/divide_heads/figures/standard_diag_%s.png' % key, dpi=130)
np.savez('/tmp/standard_diag_%s.npz' % key, ba=ba, bs=bs, A_roll=A_roll, agrid=agrid,
         Pge=Pge, beta=beta, T_field=T_FIELD, site=cfg['name'])
print('  wrote dev/divide_heads/figures/standard_diag_%s.png' % key)
