\
"""Experiment 02 — A small biological energy-transport network.

Run:
    python experiments\02_three_site_reference.py

This experiment introduces:
- different site energies
- unequal couplings
- a target site
- interference between multiple possible pathways
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np

from src.biological_model import three_site_example
from src.analysis import inverse_participation_ratio, target_site_peak
from src.visualization import plot_site_populations, save_figure


def main():
    H, times, populations = three_site_example()

    print("Three-site biological Hamiltonian:")
    print(H)
    print()

    peak, index = target_site_peak(populations, target_site=2)
    print(
        "Target Site 3 maximum population: "
        f"{peak:.4f} at t = {times[index]:.4f}"
    )

    ipr = inverse_participation_ratio(populations)
    print(f"Minimum IPR (most spread-out point): {ipr.min():.4f}")

    fig = plot_site_populations(
        times,
        populations,
        title="Experiment 02 — Three-Site Biological Reference Model",
    )

    output = ROOT / "results" / "figures" / "02_three_site_reference.png"
    save_figure(fig, output)

    print(f"Saved figure to: {output}")
    plt.show()


if __name__ == "__main__":
    main()
