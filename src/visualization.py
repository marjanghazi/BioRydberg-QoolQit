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


SITE_COLORS = ("#00B8D9", "#7C3AED", "#F59E0B", "#10B981", "#EF4444", "#3B82F6", "#EC4899")


def apply_publication_style() -> None:
    """Apply a restrained, readable style shared by all submission plots."""
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "figure.dpi": 120,
            "savefig.dpi": 240,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def plot_site_populations(
    times: np.ndarray,
    populations: np.ndarray,
    title: str,
    ylabel: str = "Population",
):
    """Plot one curve per site."""
    fig, ax = plt.subplots(figsize=(9, 5))

    for site in range(populations.shape[1]):
        ax.plot(
            times,
            populations[:, site],
            label=f"Site {site + 1}",
            color=SITE_COLORS[site % len(SITE_COLORS)],
            linewidth=2.0,
        )

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


def save_figure_pair(fig, stem: str | Path) -> tuple[Path, Path]:
    """Save a figure as high-resolution PNG and scalable SVG."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    svg = stem.with_suffix(".svg")
    fig.savefig(png, dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    return png, svg


def plot_interaction_matrix(matrix: np.ndarray, title: str):
    """Plot the dimensionless 1/r^6 interaction matrix."""
    matrix = np.asarray(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(5.8, 4.8))
    image = ax.imshow(matrix, cmap="magma", vmin=0.0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{matrix[i, j]:.3f}", ha="center", va="center", color="white" if matrix[i, j] > 0.45 else "black", fontsize=8)
    ax.set_xticks(range(matrix.shape[0]), [f"A{i + 1}" for i in range(matrix.shape[0])])
    ax.set_yticks(range(matrix.shape[0]), [f"A{i + 1}" for i in range(matrix.shape[0])])
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label=r"Dimensionless interaction $1/r^6$")
    fig.tight_layout()
    return fig


def plot_raw_and_total(relative_times: np.ndarray, occupations: np.ndarray):
    """Show site-resolved occupations and their non-conserved total."""
    fig, (ax, total_ax) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    for site in range(occupations.shape[1]):
        ax.plot(relative_times, occupations[:, site], color=SITE_COLORS[site % len(SITE_COLORS)], linewidth=2, label=f"Atom {site + 1}")
    ax.set_ylabel(r"Rydberg occupation $\langle n_i\rangle$")
    ax.set_title("QoolQit site occupations from a localized initial excitation")
    ax.legend(ncol=min(occupations.shape[1], 4))
    ax.grid(alpha=0.2)
    total_ax.plot(relative_times, occupations.sum(axis=1), color="#111827", linewidth=2)
    total_ax.axhline(1.0, color="#9CA3AF", linestyle="--", linewidth=1)
    total_ax.set_xlabel("Relative program time")
    total_ax.set_ylabel("Total")
    total_ax.grid(alpha=0.2)
    fig.tight_layout()
    return fig


def plot_baseline_optimized(
    biological: np.ndarray,
    baseline_spatial: np.ndarray,
    optimized_spatial: np.ndarray,
    baseline_rmse: float,
    optimized_rmse: float,
):
    """Compare biological profiles with baseline and optimized QoolQit runs."""
    n = len(baseline_spatial)
    bio = np.zeros((n, biological.shape[1]))
    old = np.linspace(0.0, 1.0, len(biological))
    new = np.linspace(0.0, 1.0, n)
    for site in range(biological.shape[1]):
        bio[:, site] = np.interp(new, old, biological[:, site])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True, sharey=True)
    for ax, data, name, score in (
        (axes[0], baseline_spatial, "Coupling-informed baseline", baseline_rmse),
        (axes[1], optimized_spatial, "Optimized geometry and controls", optimized_rmse),
    ):
        for site in range(bio.shape[1]):
            color = SITE_COLORS[site % len(SITE_COLORS)]
            ax.plot(new, bio[:, site], color=color, linewidth=2, label=f"Bio {site + 1}")
            ax.plot(new, data[:, site], color=color, linewidth=1.8, linestyle="--", label=f"Atom {site + 1}")
        ax.set_title(f"{name}\nheld-out RMSE = {score:.3f}")
        ax.set_xlabel("Relative evolution time")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Normalized spatial profile")
    axes[1].legend(ncol=2, fontsize=8, loc="upper right")
    fig.suptitle("Biological reference and QoolQit spatial profiles", fontsize=15)
    fig.tight_layout()
    return fig
