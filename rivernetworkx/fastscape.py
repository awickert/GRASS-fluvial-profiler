"""
FastScape: implicit stream-power landscape-evolution solver (Braun & Willett, 2013,
Geomorphology, doi:10.1016/j.geomorph.2012.10.008).

Solves the detachment-limited stream-power incision model

    dz/dt = U - K * A^m * S^n

on a D8 drainage network. The implicit, stack-ordered update is unconditionally
stable for arbitrarily large timesteps -- its defining advantage. This reuses the
Braun & Willett FlowInfo machinery already built for the DrEICH port
(``rivernetworkx.dreich``): D8 receivers, the ordered node stack, and contributing
area.

The GRASS module ``r.fluvial.fastscape`` is a thin wrapper around :func:`evolve`.

Method per timestep:
  1. Add uplift to the interior (non-fixed) nodes:  z += U*dt.
  2. Route flow: depression-fill for receivers (so flow crosses pits), build the
     D8 FlowInfo, contributing area, and the ordered stack.
  3. Implicit incision in stack order (downstream node updated before its donors,
     so the receiver elevation z_r^{t+1} is already known):
        n == 1 (linear, closed form):
            z_i = (z_i* + C * z_r) / (1 + C),   C = K*dt*A^m / dx
        n != 1 (Newton-Raphson, Braun & Willett Appendix):
            solve  z_i - z_i* + K*dt*A^m*((z_i - z_r)/dx)^n = 0,   z_i >= z_r
     where z_i* is the post-uplift elevation and dx the flow length to the receiver.
  Fixed boundary nodes (default: the four domain edges) are held constant -- they
  are the base level the interior erodes toward.

A is drainage area in map units^2 (contributing cells * cellsize^2); dx is cellsize
for a cardinal receiver, cellsize*sqrt(2) for a diagonal one.
"""
import numpy as np

from .dreich import build_flowinfo, contributing_area, build_svector, fill


def _erode_stack(zn, fi, area, K, m, n, dt, fixed):
    """In-place implicit stream-power incision on node-indexed elevations ``zn``,
    processed in ordered-stack (downstream-first) order."""
    SVector = fi['SVector']; recv = fi['recv']; flc = fi['flc']; res = fi['res']
    diag = np.sqrt(2.0) * res
    Km = float(K) * float(dt)
    for node in SVector:
        r = recv[node]
        if r == node or fixed[node]:
            continue                                   # base level / fixed: held
        dx = diag if flc[node] == 2 else res
        C = Km * (area[node] ** m) / (dx ** n)
        zr = zn[r]
        zstar = zn[node]                               # already uplifted
        if zstar <= zr:
            zn[node] = zr                              # pit / lake: deposit to receiver level
            continue
        if n == 1.0:
            zn[node] = (zstar + C * zr) / (1.0 + C)
        else:
            # Newton-Raphson; start from the no-erosion elevation
            z = zstar
            for _ in range(40):
                s = (z - zr) / dx
                if s < 0:
                    z = zr
                    break
                f = z - zstar + Km * (area[node] ** m) * (s ** n)
                fp = 1.0 + Km * (area[node] ** m) * n * (s ** (n - 1)) / dx
                dz = f / fp
                z -= dz
                if z < zr:
                    z = zr
                if abs(dz) < 1e-6:
                    break
            zn[node] = z


def evolve(z, nodata=-9999.0, cellsize=1.0, *, K=1e-5, m=0.5, n=1.0,
           uplift=1e-3, dt=1000.0, nsteps=100, fixed_boundary='edges',
           min_slope=1e-4, record_every=0):
    """Evolve a DEM under uniform-or-spatial uplift and stream-power incision.

    Parameters
    ----------
    z : 2-D ndarray
        Initial elevation (row 0 = north); ``nodata`` marks no-data cells.
    nodata, cellsize : float
        No-data sentinel and (square) cell size in map units.
    K : float or 2-D ndarray
        Erodibility (scalar or per-cell).
    m, n : float
        Stream-power area and slope exponents (e.g. 0.5, 1.0).
    uplift : float or 2-D ndarray
        Rock-uplift rate (length/time), applied to interior nodes.
    dt : float
        Timestep (time units consistent with K and uplift).
    nsteps : int
        Number of timesteps.
    fixed_boundary : {'edges', None} or 2-D bool ndarray
        Nodes held at constant elevation (base level). 'edges' fixes the domain
        border; pass a boolean mask for custom outlets; None fixes only natural
        base-level pits.
    min_slope : float
        Depression-fill minimum gradient used for flow routing each step.
    record_every : int
        If > 0, also return a list of (step, z.copy()) snapshots at that interval.

    Returns
    -------
    z_final : 2-D ndarray
        Evolved elevation. If ``record_every`` > 0, returns ``(z_final, history)``.
    """
    ndf = np.float32(nodata)
    z = z.astype(np.float64).copy()
    nr, nc = z.shape
    data = z != float(ndf)

    if fixed_boundary is None:
        fixed_grid = np.zeros((nr, nc), bool)
    elif isinstance(fixed_boundary, str) and fixed_boundary == 'edges':
        fixed_grid = np.zeros((nr, nc), bool)
        fixed_grid[0, :] = fixed_grid[-1, :] = True
        fixed_grid[:, 0] = fixed_grid[:, -1] = True
        fixed_grid &= data
    else:
        fixed_grid = np.asarray(fixed_boundary, bool)
    interior = data & ~fixed_grid

    U = uplift if np.isscalar(uplift) else np.asarray(uplift, float)
    cell_area = cellsize * cellsize
    history = []

    for step in range(nsteps):
        # 1. uplift the interior
        if np.isscalar(U):
            z[interior] += U * dt
        else:
            z[interior] += U[interior] * dt

        # 2. route: fill for receivers, build FlowInfo + stack + area
        zf = z.astype(np.float32).copy()
        zf[~data] = ndf
        filled = fill(zf, nodata, min_slope, cellsize)
        fi = build_flowinfo(filled, nodata, cellsize)
        contributing_area(fi)
        build_svector(fi)
        area = fi['ncontrib'].astype(np.float64) * cell_area

        # 3. implicit incision on the actual (un-filled) elevations
        zn = z[fi['row_of'], fi['col_of']].astype(np.float64)
        Kn = (K if np.isscalar(K) else np.asarray(K, float)[fi['row_of'], fi['col_of']])
        fixed_n = fixed_grid[fi['row_of'], fi['col_of']]
        if np.isscalar(Kn):
            _erode_stack(zn, fi, area, float(Kn), m, n, dt, fixed_n)
        else:
            _erode_stack_varK(zn, fi, area, Kn, m, n, dt, fixed_n)
        z[fi['row_of'], fi['col_of']] = zn

        if record_every and (step % record_every == 0 or step == nsteps - 1):
            history.append((step, z.copy()))

    out = z.astype(np.float32)
    out[~data] = ndf
    return (out, history) if record_every else out


def _erode_stack_varK(zn, fi, area, Kn, m, n, dt, fixed):
    """Spatially-variable-K variant of :func:`_erode_stack`."""
    SVector = fi['SVector']; recv = fi['recv']; flc = fi['flc']; res = fi['res']
    diag = np.sqrt(2.0) * res
    for node in SVector:
        r = recv[node]
        if r == node or fixed[node]:
            continue
        dx = diag if flc[node] == 2 else res
        C = float(Kn[node]) * dt * (area[node] ** m) / (dx ** n)
        zr = zn[r]; zstar = zn[node]
        if zstar <= zr:
            zn[node] = zr
            continue
        if n == 1.0:
            zn[node] = (zstar + C * zr) / (1.0 + C)
        else:
            z = zstar
            Km = float(Kn[node]) * dt
            for _ in range(40):
                s = (z - zr) / dx
                if s < 0:
                    z = zr; break
                f = z - zstar + Km * (area[node] ** m) * (s ** n)
                fp = 1.0 + Km * (area[node] ** m) * n * (s ** (n - 1)) / dx
                dz = f / fp
                z -= dz
                if z < zr:
                    z = zr
                if abs(dz) < 1e-6:
                    break
            zn[node] = z


def steady_state_relief(K, m, n, uplift, area, slope_exponent_area=None):
    """Analytic steady-state channel slope S = (U/K)^(1/n) * A^(-m/n) at a point
    (dz/dt = 0 in the stream-power model). Useful for sanity checks / initial
    guesses. ``area`` is drainage area in map units^2."""
    return (uplift / K) ** (1.0 / n) * area ** (-m / n)
