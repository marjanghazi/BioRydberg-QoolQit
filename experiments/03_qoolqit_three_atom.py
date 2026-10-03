\
"""Experiment 03 — Build the corresponding three-atom QoolQit program.

Run:
    python experiments\03_qoolqit_three_atom.py

Learning goals:
1. Convert biological coupling strengths into a geometry heuristic.
2. Build a QoolQit Register from coordinates.
3. Apply a global Drive.
4. Build and compile a QuantumProgram.
5. Execute with LocalEmulator.
6. Measure Rydberg occupation on every atom.

Scientific warning:
The Rydberg occupations are NOT the same quantity as a conserved biological
single-excitation population. This experiment is an analogue/inspired mapping.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np

from src.biological_model import three_site_example
from src.neutral_atom_mapping import (
    couplings_to_qoolqit_coordinates,
    rydberg_interaction_matrix,
)
from src.qoolqit_simulation import run_occupation_dynamics
from src.visualization import (
    plot_network,
    plot_site_populations,
    save_figure,
)


def main():
    # Recreate the coupling matrix from our biological Hamiltonian.
    H, _, _ = three_site_example()

    # The off-diagonal elements are biological hopping couplings.
    couplings = H.copy()
    np.fill_diagonal(couplings, 0.0)

    coordinates = couplings_to_qoolqit_coordinates(couplings)

    print("QoolQit dimensionless coordinates:")
    for i, xy in enumerate(coordinates):
        print(f"Atom {i + 1}: {xy}")

    print("\nOur calculated 1/r^6 matrix:")
    print(rydberg_interaction_matrix(coordinates))

    run = run_occupation_dynamics(
        coordinates=coordinates,
        duration=10.0,
        amplitude=0.8,
        detuning=0.0,
        n_times=31,
    )

    print("\nQoolQit Register interaction matrix:")
    print(run.interaction_matrix)

    # Plot the geometry.
    fig_network = plot_network(
        coordinates,
        coupling_matrix=couplings,
        title="Biological Couplings Mapped to a QoolQit Atom Geometry",
    )

    save_figure(
        fig_network,
        ROOT / "results" / "figures" / "03_qoolqit_geometry.png",
    )

    # Plot raw Rydberg occupations.
    fig_occ = plot_site_populations(
        run.relative_times,
        run.occupations,
        title="Experiment 03 — QoolQit Rydberg Occupation Dynamics",
        ylabel="Rydberg occupation <n_i>",
    )

    save_figure(
        fig_occ,
        ROOT / "results" / "figures" / "03_qoolqit_occupations.png",
    )

    plt.show()


if __name__ == "__main__":
    main()
