"""Run the reproducible BioRydberg 2.0 contest study.

This experiment performs real QoolQit local emulation, bounded geometry and
control optimization, held-out evaluation, and a limited seven-site benchmark.
"""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.analysis import normalize_rows, target_site_peak
from src.biological_model import load_hamiltonian_csv, simulate_unitary_transport, three_site_example, two_site_example
from src.neutral_atom_mapping import couplings_to_qoolqit_coordinates, rydberg_interaction_matrix, triangle_coordinates
from src.optimization import Controls, evaluate, optimize_controls, optimize_geometry
from src.qoolqit_simulation import run_occupation_dynamics
from src.visualization import (
    apply_publication_style,
    plot_baseline_optimized,
    plot_interaction_matrix,
    plot_network,
    plot_raw_and_total,
    plot_site_populations,
    save_figure_pair,
)


def _serializable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)}")


def run_study(*, maxiter: int = 3, popsize: int = 4, n_times: int = 25) -> dict:
    apply_publication_style()
    figure_dir = ROOT / "results" / "figures"
    table_dir = ROOT / "results" / "tables"
    config_dir = ROOT / "results" / "configs"
    for directory in (figure_dir, table_dir, config_dir):
        directory.mkdir(parents=True, exist_ok=True)

    study_start = perf_counter()
    two_h, two_t, two_p = two_site_example()
    three_h, bio_t, bio_p = three_site_example()
    conservation = {
        "two_site": float(np.max(np.abs(two_p.sum(axis=1) - 1.0))),
        "three_site": float(np.max(np.abs(bio_p.sum(axis=1) - 1.0))),
    }
    peak, peak_index = target_site_peak(bio_p, 2)

    fig = plot_site_populations(bio_t, bio_p, "Three-site biological excitation dynamics", "Excitation probability")
    save_figure_pair(fig, figure_dir / "01_three_site_biological")
    plt.close(fig)

    couplings = three_h.copy()
    np.fill_diagonal(couplings, 0.0)
    baseline_coords = couplings_to_qoolqit_coordinates(couplings, seed=42)
    baseline_controls = Controls()
    baseline_run, baseline_metrics = evaluate(bio_p, baseline_coords, baseline_controls, n_times=n_times)

    fig = plot_network(baseline_coords, baseline_run.interaction_matrix, "Coupling-informed three-atom register")
    save_figure_pair(fig, figure_dir / "02_baseline_geometry")
    plt.close(fig)
    fig = plot_interaction_matrix(baseline_run.interaction_matrix, "Baseline Rydberg interaction matrix")
    save_figure_pair(fig, figure_dir / "03_interaction_matrix")
    plt.close(fig)
    fig = plot_raw_and_total(baseline_run.relative_times, baseline_run.occupations)
    save_figure_pair(fig, figure_dir / "04_raw_rydberg_occupations")
    plt.close(fig)

    geometry_trials = [
        optimize_geometry(bio_p, baseline_controls, seed=seed, maxiter=maxiter, popsize=popsize, n_times=n_times)
        for seed in (11, 29)
    ]
    best_geometry = min(geometry_trials, key=lambda item: item.objective)
    optimized_coords = triangle_coordinates(*best_geometry.parameters)
    geometry_run, geometry_metrics = evaluate(bio_p, optimized_coords, baseline_controls, n_times=n_times)

    control_trials = [
        optimize_controls(bio_p, optimized_coords, seed=seed, maxiter=maxiter, popsize=popsize, n_times=n_times)
        for seed in (23, 41)
    ]
    best_control = min(control_trials, key=lambda item: item.objective)
    optimized_controls = Controls(
        amplitude=float(best_control.parameters[0]),
        detuning=float(best_control.parameters[1]),
        duration=float(best_control.parameters[2]),
    )
    optimized_run, optimized_metrics = evaluate(bio_p, optimized_coords, optimized_controls, n_times=n_times)

    baseline_spatial = normalize_rows(baseline_run.occupations)
    optimized_spatial = normalize_rows(optimized_run.occupations)
    fig = plot_baseline_optimized(
        bio_p,
        baseline_spatial,
        optimized_spatial,
        float(baseline_metrics["rmse_held_out"]),
        float(optimized_metrics["rmse_held_out"]),
    )
    save_figure_pair(fig, figure_dir / "05_baseline_vs_optimized")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.6))
    labels = []
    for name, trials in (("Geometry", geometry_trials), ("Controls", control_trials)):
        for trial in trials:
            labels.append(f"{name} seed {trial.seed}")
            history = trial.history or [trial.objective]
            ax.plot(np.arange(1, len(history) + 1), history, marker="o", label=labels[-1])
    ax.axhline(float(baseline_metrics["rmse_tuning"]), color="#111827", linestyle="--", label="Baseline")
    ax.set_xlabel("Differential-evolution generation")
    ax.set_ylabel("Tuning-point spatial RMSE")
    ax.set_title("Reproducible optimization histories")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    save_figure_pair(fig, figure_dir / "06_optimization_history")
    plt.close(fig)

    # Limited seven-site benchmark: one compiled local-emulator run, no search.
    seven_h = load_hamiltonian_csv(ROOT / "data" / "fmo_inspired_7site.csv")
    seven_c = seven_h.copy()
    np.fill_diagonal(seven_c, 0.0)
    seven_t = np.linspace(0.0, 12.0, 240)
    seven_bio = simulate_unitary_transport(seven_h, seven_t, initial_site=0)
    seven_coords = couplings_to_qoolqit_coordinates(seven_c, seed=42)
    seven_run = run_occupation_dynamics(
        seven_coords,
        duration=optimized_controls.duration,
        amplitude=optimized_controls.amplitude,
        detuning=optimized_controls.detuning,
        n_times=17,
        initial_site=0,
    )
    from src.analysis import split_profile_rmse
    seven_metrics = split_profile_rmse(seven_bio, seven_run.occupations)
    seven_peak, seven_peak_idx = target_site_peak(seven_bio, 6)
    fig = plot_raw_and_total(seven_run.relative_times, seven_run.occupations)
    fig.axes[0].set_title("Seven-site illustrative benchmark: raw QoolQit occupations")
    save_figure_pair(fig, figure_dir / "07_seven_site_benchmark")
    plt.close(fig)

    rows = [
        {"model": "baseline", **baseline_metrics, "amplitude": baseline_controls.amplitude, "detuning": baseline_controls.detuning, "duration": baseline_controls.duration},
        {"model": "geometry_optimized", **geometry_metrics, "amplitude": baseline_controls.amplitude, "detuning": baseline_controls.detuning, "duration": baseline_controls.duration},
        {"model": "geometry_and_controls_optimized", **optimized_metrics, "amplitude": optimized_controls.amplitude, "detuning": optimized_controls.detuning, "duration": optimized_controls.duration},
        {"model": "seven_site_limited", **seven_metrics, "amplitude": optimized_controls.amplitude, "detuning": optimized_controls.detuning, "duration": optimized_controls.duration},
    ]
    pd.DataFrame(rows).to_csv(table_dir / "study_metrics.csv", index=False)

    summary = {
        "title": "BioRydberg 2.0: Geometry-Aware Quantum Transport with Rydberg Atom Arrays",
        "qoolqit_version": version("qoolqit"),
        "python_version": sys.version.split()[0],
        "seeds": {"layout": 42, "geometry": [11, 29], "controls": [23, 41]},
        "n_times": n_times,
        "optimizer": {"name": "scipy differential_evolution", "maxiter": maxiter, "popsize": popsize},
        "normalization_threshold": 1e-6,
        "split_policy": "even valid time indices tune; odd valid time indices held out",
        "probability_conservation_max_error": conservation,
        "biological_target_site_3_peak": {"population": peak, "time": float(bio_t[peak_index])},
        "baseline": {"coordinates": baseline_coords, "controls": baseline_controls.__dict__, "metrics": baseline_metrics, "emulator_seconds": baseline_run.elapsed_seconds},
        "geometry_optimized": {"coordinates": optimized_coords, "metrics": geometry_metrics, "best_seed": best_geometry.seed, "evaluations": best_geometry.evaluations},
        "fully_optimized": {"controls": optimized_controls.__dict__, "metrics": optimized_metrics, "best_seed": best_control.seed, "evaluations": best_control.evaluations, "emulator_seconds": optimized_run.elapsed_seconds},
        "seven_site": {"description": "illustrative dataset; one limited benchmark without optimization", "metrics": seven_metrics, "emulator_seconds": seven_run.elapsed_seconds, "target_site_7_bio_peak": seven_peak, "target_site_7_bio_peak_time": float(seven_t[seven_peak_idx])},
        "total_study_seconds": perf_counter() - study_start,
        "scientific_scope": "Normalized occupations compare spatial patterns only; they are not biological excitation probabilities or quantum-state fidelities.",
    }
    with (ROOT / "results" / "study_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, default=_serializable)
    with (config_dir / "contest_study.json").open("w", encoding="utf-8") as handle:
        json.dump({key: summary[key] for key in ("qoolqit_version", "python_version", "seeds", "n_times", "optimizer", "normalization_threshold", "split_policy")}, handle, indent=2)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=3)
    parser.add_argument("--popsize", type=int, default=4)
    parser.add_argument("--n-times", type=int, default=25)
    args = parser.parse_args()
    summary = run_study(maxiter=args.maxiter, popsize=args.popsize, n_times=args.n_times)
    print(json.dumps(summary, indent=2, default=_serializable))


if __name__ == "__main__":
    main()
