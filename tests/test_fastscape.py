"""Unit tests for the FastScape solver (rivernetworkx.fastscape).

The decisive physical check: a uniformly-uplifted domain run to steady state must
satisfy the stream-power balance U = K*A^m*S^n, i.e. the channel slope-area relation
S ~ A^(-m/n). We also check implicit stability (huge dt does not blow up) and that
the solution converges (dz/dt -> 0).
"""
import numpy as np
import pytest

from rivernetworkx import fastscape as F
from rivernetworkx.dreich import build_flowinfo, contributing_area


def _run_to_steady(nr=25, nc=25, K=1e-4, m=0.5, n=1.0, U=1e-3, dt=2000.0, nsteps=400):
    rng = np.random.RandomState(1)
    z = (rng.rand(nr, nc) * 1.0).astype(np.float32)        # small initial roughness
    z[0, :] = z[-1, :] = z[:, 0] = z[:, -1] = 0.0          # base-level edges
    zf = F.evolve(z, nodata=-9999.0, cellsize=10.0, K=K, m=m, n=n,
                  uplift=U, dt=dt, nsteps=nsteps, fixed_boundary='edges')
    return zf


def test_evolve_runs_and_is_bounded():
    zf = _run_to_steady(nsteps=100)
    assert np.all(np.isfinite(zf))
    assert zf.max() > 0.0                                  # relief built up
    # implicit scheme is stable: no runaway with large dt
    assert zf.max() < 1e4


def test_converges_to_steady_state():
    # run, then run a bit more; steady state means the surface barely changes
    z1 = _run_to_steady(nsteps=400)
    z2 = F.evolve(z1, nodata=-9999.0, cellsize=10.0, K=1e-4, m=0.5, n=1.0,
                  uplift=1e-3, dt=2000.0, nsteps=50, fixed_boundary='edges')
    interior = np.ones_like(z1, bool)
    interior[0, :] = interior[-1, :] = interior[:, 0] = interior[:, -1] = False
    change = np.abs(z2[interior] - z1[interior]).max() / max(z1[interior].max(), 1e-9)
    assert change < 0.05                                   # < 5% relief change -> near steady


def test_steady_state_slope_area_exponent():
    m, n = 0.5, 1.0
    cellsize = 10.0
    zf = _run_to_steady(nr=31, nc=31, K=1e-4, m=m, n=n, U=1e-3, dt=2000.0, nsteps=700)
    # recompute slope + area on the final surface
    fi = build_flowinfo(zf, nodata=-9999.0, cellsize=cellsize)
    contributing_area(fi)
    recv = fi['recv']; flc = fi['flc']; res = cellsize
    diag = np.sqrt(2.0) * res
    znode = zf[fi['row_of'], fi['col_of']].astype(float)
    area = fi['ncontrib'].astype(float) * cellsize ** 2
    slope = np.zeros(fi['N'])
    for i in range(fi['N']):
        r = recv[i]
        if r == i:
            continue
        dx = diag if flc[i] == 2 else res
        slope[i] = max((znode[i] - znode[r]) / dx, 0.0)
    # channel cells: drop tiny-area hillslope cells and flats
    chan = (area > 20 * cellsize ** 2) & (slope > 1e-6)
    logA = np.log10(area[chan]); logS = np.log10(slope[chan])
    coef = np.polyfit(logA, logS, 1)
    # theoretical exponent is -m/n = -0.5; allow tolerance for a small, finite domain
    assert -0.5 - 0.2 < coef[0] < -0.5 + 0.2


def test_nonlinear_n_runs():
    # n != 1 exercises the Newton path; just check it stays finite/bounded
    zf = _run_to_steady(nr=20, nc=20, K=1e-4, m=0.45, n=2.0, U=1e-3, dt=1000.0, nsteps=150)
    assert np.all(np.isfinite(zf)) and zf.max() < 1e4


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-v']))
