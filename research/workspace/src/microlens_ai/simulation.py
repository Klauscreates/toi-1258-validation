from __future__ import annotations

import math
from collections.abc import Iterable


def paczynski_magnification(time: float, t0: float, t_e: float, u0: float) -> float:
    """Point-source point-lens microlensing magnification."""

    if t_e <= 0:
        raise ValueError("Einstein crossing time t_e must be positive")
    if u0 < 0:
        raise ValueError("impact parameter u0 must be non-negative")

    tau = (time - t0) / t_e
    u = math.sqrt(u0 * u0 + tau * tau)
    return (u * u + 2.0) / (u * math.sqrt(u * u + 4.0))


def paczynski_flux(
    time: float,
    t0: float,
    t_e: float,
    u0: float,
    source_flux: float,
    blend_flux: float,
) -> float:
    """Observed blended flux for a simple microlensing event."""

    return source_flux * paczynski_magnification(time, t0, t_e, u0) + blend_flux


def paczynski_series(
    times: Iterable[float],
    t0: float,
    t_e: float,
    u0: float,
    source_flux: float,
    blend_flux: float,
) -> list[float]:
    """Generate a deterministic microlensing light curve."""

    return [
        paczynski_flux(time, t0, t_e, u0, source_flux, blend_flux)
        for time in times
    ]
