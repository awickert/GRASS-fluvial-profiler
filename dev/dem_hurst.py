"""Independent estimate of the Hurst exponent H of the Mid Bailey Run 1 m DEM,
to test the self-affine basis of the DrEICH curvature-threshold scaling.

The threshold calibration (window = 7*res) gave tan_curv ~ res^(-1.12), which
under kappa ~ L^(H-2) (L = measurement scale) implies H = 0.88. If H measured
directly from the terrain (no channel-head info) agrees, the scaling is
physically grounded, not curve-fitting.

Two independent estimators over the 1-30 m band that matters here:
  - 2nd-order structure function:  gamma(L) = <(z(x+L)-z(x))^2> ~ L^(2H)
  - radially-averaged power spectrum: P(k) ~ k^(-(2H+2))
"""
import numpy as np

ALG = '/tmp/dreich_algorithm'
NR, NC = 7107, 6266
ND = -9999.0


def main():
    z = np.fromfile('%s/bailey_run_dem.flt' % ALG, '<f4').reshape(NR, NC).astype(np.float64)
    z = np.where(z == ND, np.nan, z)

    # --- structure function (cardinal x and y), lags in cells = metres ---
    lags = np.array([1, 2, 3, 4, 5, 8, 12, 20, 30, 50, 80, 120])
    gx = np.array([np.nanmean((z[:, L:] - z[:, :-L]) ** 2) for L in lags])
    gy = np.array([np.nanmean((z[L:, :] - z[:-L, :]) ** 2) for L in lags])
    fit = (lags >= 2) & (lags <= 30)                       # the resolution band
    for tag, g in [('x', gx), ('y', gy)]:
        s = np.polyfit(np.log(lags[fit]), np.log(g[fit]), 1)[0]
        print('structure fn (%s): slope=2H=%.3f -> H=%.3f' % (tag, s, s / 2))
    s_iso = np.polyfit(np.log(lags[fit]), np.log(0.5 * (gx + gy)[fit]), 1)[0]
    print('structure fn (iso): 2H=%.3f -> H=%.3f' % (s_iso, s_iso / 2))

    # --- power spectrum on the largest NaN-free square we can find ---
    S = 2048
    block = None
    for r0 in (1500, 2500, 3500):
        for c0 in (3200, 4000, 2400):
            b = z[r0:r0 + S, c0:c0 + S]
            if b.shape == (S, S) and not np.isnan(b).any():
                block = b
                break
        if block is not None:
            break
    if block is None:
        print('power spectrum: no NaN-free %dx%d block found, skipping' % (S, S))
        return
    b = block - block.mean()
    win = np.hanning(S)[:, None] * np.hanning(S)[None, :]   # reduce spectral leakage
    P = np.abs(np.fft.fftshift(np.fft.fft2(b * win))) ** 2
    ky = np.fft.fftshift(np.fft.fftfreq(S))[:, None] * np.ones((1, S))
    kx = np.ones((S, 1)) * np.fft.fftshift(np.fft.fftfreq(S))[None, :]
    kr = np.sqrt(kx ** 2 + ky ** 2)
    kbin = np.round(kr * S).astype(int)
    nb = np.bincount(kbin.ravel(), P.ravel())
    cnt = np.bincount(kbin.ravel())
    Pk = nb[1:] / np.maximum(cnt[1:], 1)
    k = np.arange(1, len(Pk) + 1) / S                      # cycles per metre
    band = (1.0 / 30 <= k) & (k <= 1.0 / 2)                # 2-30 m wavelengths
    alpha = -np.polyfit(np.log(k[band]), np.log(Pk[band]), 1)[0]
    print('power spectrum: P~k^-%.2f -> H=(alpha-2)/2=%.3f' % (alpha, (alpha - 2) / 2))
    print('(threshold-scaling fit implied H = 0.88)')


if __name__ == '__main__':
    main()
