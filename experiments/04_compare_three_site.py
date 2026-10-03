\
"""Experiment 04 — Compare biological and QoolQit spatial population patterns.

Run:
    python experiments\04_compare_three_site.py

This is the main competition-story experiment.

We compare:
- biological single-excitation population distribution
- normalized QoolQit Rydberg occupation distribution

We compare on RELATIVE time rather than claiming the two models have the same
physical timescale.

The resulting RMSE is a pattern-similarity metric, NOT proof that the
Hamiltonians are equivalent.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.biological_model import three_site_example
from src.neutral_atom_mapping import couplings_to_qoolqit_coordinates
from src.qoolqit_simulation import run_occupation_dynamics
from src.analysis import (
    normalize_rows,
    spatial_profile_rmse,
)
from src.visualization import plot_comparison, save_figure


def main():
    H, bio_times, bio_populations = three_site_example()

    couplings = H.copy()
    np.fill_diagonal(couplings, 0.0)

    coordinates = couplings_to_qoolqit_coordinates(couplings)

    qrun = run_occupation_dynamics(
        coordinates=coordinates,
        duration=10.0,
        amplitude=0.8,
        detuning=0.0,
        n_times=31,
    )

    # QoolQit can have total Rydberg occupation != 1.
    # Normalize only to compare WHERE the occupation is distributed.
    qoolqit_spatial = normalize_rows(qrun.occupations)

    rmse = spatial_profile_rmse(
        bio_populations,
        qrun.occupations,
    )

    print(f"Spatial-profile RMSE: {rmse:.6f}")
    print(
        "Interpretation: lower means the normalized site-distribution curves "
        "are more similar. It does not mean the Hamiltonians are identical."
    )

    fig = plot_comparison(
        biological=bio_populations,
        rydberg_normalized=qoolqit_spatial,
    )

    figure_path = (
        ROOT / "results" / "figures" / "04_biological_vs_qoolqit.png"
    )
    save_figure(fig, figure_path)

    # Save a very small summary table that can later be used in the slides.
    table = pd.DataFrame(
        {
            "metric": ["spatial_profile_rmse"],
            "value": [rmse],
            "note": [
                "Comparison of normalized spatial distributions on relative time"
            ],
        }
    )

    table_path = ROOT / "results" / "tables" / "comparison_metrics.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(table_path, index=False)

    print(f"Saved figure to: {figure_path}")
    print(f"Saved metric table to: {table_path}")

    plt.show()


if __name__ == "__main__":
    main()
