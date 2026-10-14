"""Bounded, reproducible QoolQit parameter optimization."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Callable

import numpy as np
from scipy.optimize import differential_evolution

from .analysis import split_profile_rmse
from .neutral_atom_mapping import triangle_coordinates, validate_coordinates
from .qoolqit_simulation import QoolQitRun, run_occupation_dynamics


@dataclass(frozen=True)
class Controls:
    duration: float = 10.0
    amplitude: float = 0.8
    detuning: float = 0.0


@dataclass
class OptimizationResult:
    parameters: np.ndarray
    objective: float
    evaluations: int
    elapsed_seconds: float
    history: list[float]
    seed: int


def evaluate(
    biological_populations: np.ndarray,
    coordinates: list[tuple[float, float]],
    controls: Controls,
    *,
    n_times: int = 25,
    occupation_threshold: float = 1e-6,
) -> tuple[QoolQitRun, dict[str, float | int]]:
    """Execute QoolQit and compute masked tuning/held-out metrics."""
    run = run_occupation_dynamics(
        coordinates,
        duration=controls.duration,
        amplitude=controls.amplitude,
        detuning=controls.detuning,
        n_times=n_times,
        initial_site=0,
        occupation_threshold=occupation_threshold,
    )
    metrics = split_profile_rmse(
        biological_populations,
        run.occupations,
        occupation_threshold=occupation_threshold,
    )
    return run, metrics


def _optimize(
    objective: Callable[[np.ndarray], float],
    bounds: list[tuple[float, float]],
    *,
    seed: int,
    maxiter: int,
    popsize: int,
) -> OptimizationResult:
    history: list[float] = []

    def callback(intermediate_result) -> bool:
        history.append(float(intermediate_result.fun))
        return False

    start = perf_counter()
    result = differential_evolution(
        objective,
        bounds=bounds,
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        polish=False,
        updating="immediate",
        workers=1,
        callback=callback,
        tol=1e-4,
    )
    return OptimizationResult(
        parameters=np.asarray(result.x, dtype=float),
        objective=float(result.fun),
        evaluations=int(result.nfev),
        elapsed_seconds=perf_counter() - start,
        history=history,
        seed=seed,
    )


def optimize_geometry(
    biological_populations: np.ndarray,
    controls: Controls,
    *,
    seed: int = 11,
    maxiter: int = 4,
    popsize: int = 4,
    n_times: int = 25,
) -> OptimizationResult:
    """Optimize the third vertex of a normalized three-atom triangle."""

    def objective(params: np.ndarray) -> float:
        try:
            coords = triangle_coordinates(float(params[0]), float(params[1]))
            validate_coordinates(coords, maximum_extent=4.0)
            _, metrics = evaluate(
                biological_populations, coords, controls, n_times=n_times
            )
            value = float(metrics["rmse_tuning"])
            return value if np.isfinite(value) else 10.0
        except (ValueError, RuntimeError):
            return 10.0

    return _optimize(
        objective,
        bounds=[(-1.0, 2.0), (0.35, 2.5)],
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
    )


def optimize_controls(
    biological_populations: np.ndarray,
    coordinates: list[tuple[float, float]],
    *,
    seed: int = 23,
    maxiter: int = 4,
    popsize: int = 4,
    n_times: int = 25,
) -> OptimizationResult:
    """Optimize constant amplitude, detuning, and dimensionless duration."""

    def objective(params: np.ndarray) -> float:
        controls = Controls(
            amplitude=float(params[0]),
            detuning=float(params[1]),
            duration=float(params[2]),
        )
        try:
            _, metrics = evaluate(
                biological_populations, coordinates, controls, n_times=n_times
            )
            value = float(metrics["rmse_tuning"])
            return value if np.isfinite(value) else 10.0
        except (ValueError, RuntimeError):
            return 10.0

    return _optimize(
        objective,
        bounds=[(0.2, 1.8), (-1.5, 1.5), (4.0, 16.0)],
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
    )
