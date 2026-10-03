\
"""Experiment 01 — Learn excitation transfer with the smallest possible model.

Run:
    python experiments\01_two_site_reference.py

What to observe:
1. At t=0, population is entirely on Site 1.
2. Coupling causes the quantum amplitude to move between Site 1 and Site 2.
3. The total probability remains approximately 1 at every time.
"""

from pathlib import Path
import sys

# Allow this script to import src/ when executed from the repository root.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np

from src.biological_model import two_site_example
from src.visualization import plot_site_populations, save_figure


def main():
    H, times, populations = two_site_example()

    print("Two-site Hamiltonian:")
    print(H)
    print()

    # A useful correctness check:
    # a single-excitation unitary model should conserve total probability.
    total_probability = populations.sum(axis=1)
    max_error = np.max(np.abs(total_probability - 1.0))
    print(f"Maximum probability-conservation error: {max_error:.3e}")

    fig = plot_site_populations(
        times,
        populations,
        title="Experiment 01 — Two-Site Excitation Transport",
    )

    output = ROOT / "results" / "figures" / "01_two_site_transport.png"
    save_figure(fig, output)

    print(f"Saved figure to: {output}")
    plt.show()


if __name__ == "__main__":
    main()
