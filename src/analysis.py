\
"""Analysis metrics used to compare the two sides of the project."""

from __future__ import annotations

import numpy as np


def normalize_rows(
    values: np.ndarray,
    eps: float = 1e-6,
    *,
    invalid_value: float = np.nan,
    return_mask: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Normalize each time slice so site values sum to one.

    Why is this needed?
    -------------------
    Biological single-excitation populations naturally sum to 1.

    QoolQit Rydberg occupations do NOT have to sum to 1 because the driven
    many-body system may contain zero, one, or multiple Rydberg excitations.

    Therefore we normalize QoolQit occupations only when comparing their
    *spatial distribution* with the biological population distribution.

    This is a visualization/comparison transformation, not a claim of
    physical equivalence.
    """
    x = np.asarray(values, dtype=float)
    sums = x.sum(axis=1, keepdims=True)

    if x.ndim != 2:
        raise ValueError("values must be a 2D time-by-site array.")
    if eps < 0:
        raise ValueError("eps must be non-negative.")

    normalized = np.full_like(x, invalid_value)
    valid = sums[:, 0] > eps
    normalized[valid] = x[valid] / sums[valid]
    if return_mask:
        return normalized, valid
    return normalized


def resample_on_relative_time(
    values: np.ndarray,
    n_output: int,
) -> np.ndarray:
    """Linearly resample time-dependent site values to a common grid."""
    values = np.asarray(values, dtype=float)

    old_t = np.linspace(0.0, 1.0, len(values))
    new_t = np.linspace(0.0, 1.0, n_output)

    result = np.zeros((n_output, values.shape[1]), dtype=float)

    for site in range(values.shape[1]):
        result[:, site] = np.interp(new_t, old_t, values[:, site])

    return result


def spatial_profile_rmse(
    biological_populations: np.ndarray,
    rydberg_occupations: np.ndarray,
    *,
    occupation_threshold: float = 1e-6,
) -> float:
    """RMSE between normalized spatial profiles on a common relative-time grid.

    This metric deliberately compares *shape/pattern*, not absolute physical
    time and not total excitation number.
    """
    bio = normalize_rows(biological_populations)
    ryd, valid = normalize_rows(
        rydberg_occupations,
        eps=occupation_threshold,
        return_mask=True,
    )

    n = max(len(bio), len(ryd))
    bio_r = resample_on_relative_time(bio, n)
    ryd_r = resample_on_relative_time(ryd, n)

    valid_r = resample_on_relative_time(valid[:, None].astype(float), n)[:, 0] > 0.999
    finite = valid_r & np.all(np.isfinite(ryd_r), axis=1)
    if not np.any(finite):
        raise ValueError("No Rydberg profiles exceed the occupation threshold.")
    return float(np.sqrt(np.mean((bio_r[finite] - ryd_r[finite]) ** 2)))


def split_profile_rmse(
    biological_populations: np.ndarray,
    rydberg_occupations: np.ndarray,
    *,
    occupation_threshold: float = 1e-6,
) -> dict[str, float | int]:
    """Return all-point, tuning-point, and held-out profile RMSE values.

    The Rydberg data define the common time grid. Even valid time indices are
    tuning points and odd valid indices are held out from optimization.
    """
    bio = resample_on_relative_time(
        np.asarray(normalize_rows(biological_populations), dtype=float),
        len(rydberg_occupations),
    )
    ryd, valid = normalize_rows(
        rydberg_occupations,
        eps=occupation_threshold,
        return_mask=True,
    )
    valid = valid & np.all(np.isfinite(ryd), axis=1)
    indices = np.arange(len(ryd))

    def _rmse(mask: np.ndarray) -> float:
        selected = valid & mask
        if not np.any(selected):
            return float("nan")
        return float(np.sqrt(np.mean((bio[selected] - ryd[selected]) ** 2)))

    return {
        "rmse_all": _rmse(np.ones(len(ryd), dtype=bool)),
        "rmse_tuning": _rmse(indices % 2 == 0),
        "rmse_held_out": _rmse(indices % 2 == 1),
        "valid_time_points": int(np.sum(valid)),
        "excluded_time_points": int(np.sum(~valid)),
    }


def target_site_peak(
    values: np.ndarray,
    target_site: int,
) -> tuple[float, int]:
    """Return maximum target-site value and its time-index."""
    values = np.asarray(values, dtype=float)

    if not 0 <= target_site < values.shape[1]:
        raise ValueError("target_site is outside the site range.")

    column = values[:, target_site]
    index = int(np.argmax(column))

    return float(column[index]), index


def inverse_participation_ratio(populations: np.ndarray) -> np.ndarray:
    """Compute IPR(t) = sum_i P_i(t)^2.

    Interpretation for normalized populations:
        IPR near 1   -> localized on one site
        smaller IPR  -> spread across several sites
    """
    p = normalize_rows(populations)
    return np.sum(p**2, axis=1)
