\
"""Plotting utilities for BioRydberg.

All functions return matplotlib Figure objects so they work both in normal
Python scripts and in Jupyter notebooks.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx


def plot_site_populations(
    times: np.ndarray,
    populations: np.ndarray,
    title: str,
    ylabel: str = "Population",
):
    """Plot one curve per site."""
    fig, ax = plt.subplots(figsize=(9, 5))

    for site in range(populations.shape[1]):
        ax.plot(times, populations[:, site], label=f"Site {site + 1}")

    ax.set_xlabel("Time")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    return fig


def plot_network(
    coordinates: list[tuple[float, float]],
    coupling_matrix: np.ndarray | None = None,
    title: str = "Network",
):
    """Draw a site/atom network.

    If a coupling matrix is supplied, edge width roughly follows coupling
    magnitude. Coordinates are used exactly as supplied.
    """
    coords = np.asarray(coordinates, dtype=float)
    n = len(coords)

    graph = nx.Graph()
    graph.add_nodes_from(range(n))

    widths = []

    if coupling_matrix is not None:
        coupling_matrix = np.asarray(coupling_matrix, dtype=float)
        max_weight = np.max(np.abs(coupling_matrix))

        for i in range(n):
            for j in range(i + 1, n):
                w = abs(coupling_matrix[i, j])
                if w > 0:
                    graph.add_edge(i, j, weight=w)
                    widths.append(
                        0.8 + 4.0 * (w / max_weight if max_weight > 0 else 0)
                    )
    else:
        # If no coupling matrix is given, just connect all atom pairs lightly.
        for i in range(n):
            for j in range(i + 1, n):
                graph.add_edge(i, j)
                widths.append(1.0)

    pos = {i: tuple(coords[i]) for i in range(n)}

    fig, ax = plt.subplots(figsize=(6, 6))

    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=900,
        ax=ax,
    )

    nx.draw_networkx_labels(
        graph,
        pos,
        labels={i: str(i + 1) for i in range(n)},
        ax=ax,
    )

    nx.draw_networkx_edges(
        graph,
        pos,
        width=widths,
        alpha=0.65,
        ax=ax,
    )

    ax.set_title(title)
    ax.set_axis_off()
    fig.tight_layout()

    return fig


def plot_comparison(
    biological: np.ndarray,
    rydberg_normalized: np.ndarray,
):
    """Compare site patterns using relative time from 0 to 1."""
    bio_t = np.linspace(0.0, 1.0, len(biological))
    ryd_t = np.linspace(0.0, 1.0, len(rydberg_normalized))

    fig, ax = plt.subplots(figsize=(10, 5))

    for site in range(biological.shape[1]):
        ax.plot(
            bio_t,
            biological[:, site],
            label=f"Bio site {site + 1}",
        )

        ax.plot(
            ryd_t,
            rydberg_normalized[:, site],
            linestyle="--",
            label=f"QoolQit atom {site + 1}",
        )

    ax.set_xlabel("Relative evolution time")
    ax.set_ylabel("Normalized spatial population")
    ax.set_title(
        "Biological reference vs QoolQit spatial occupation pattern"
    )
    ax.grid(alpha=0.25)
    ax.legend(ncol=2)
    fig.tight_layout()

    return fig


def save_figure(fig, path: str | Path) -> None:
    """Create parent directory and save a figure."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
