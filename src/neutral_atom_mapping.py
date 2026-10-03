\
"""Map a biological coupling graph to a neutral-atom register geometry.

Physical idea
-------------
In the biological reference Hamiltonian, J_ij is a *hopping coupling*.

In QoolQit's Rydberg model, the interaction coefficient is approximately

    J_ij^(Rydberg) = 1 / r_ij^6

where r_ij is the dimensionless atom separation.

These are not the same physical interaction. We therefore use the biological
coupling magnitudes only as a *geometry-design heuristic*:

    stronger biological connection -> place atoms relatively closer
    weaker biological connection   -> place atoms relatively farther

The function below uses NetworkX's spring layout to create a 2D embedding,
then rescales the geometry so that the closest atom pair has distance 1,
which is the convenient QoolQit dimensionless convention.
"""

from __future__ import annotations

import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist


def coupling_matrix_to_graph(couplings: np.ndarray) -> nx.Graph:
    """Convert a coupling matrix into a weighted NetworkX graph."""
    couplings = np.asarray(couplings, dtype=float)
    n = couplings.shape[0]

    if couplings.shape != (n, n):
        raise ValueError("couplings must be a square matrix.")

    graph = nx.Graph()
    graph.add_nodes_from(range(n))

    for i in range(n):
        for j in range(i + 1, n):
            weight = abs(couplings[i, j])

            # Skip exactly-zero biological couplings.
            if weight > 0:
                graph.add_edge(i, j, weight=weight)

    return graph


def graph_to_qoolqit_coordinates(
    graph: nx.Graph,
    seed: int = 42,
) -> list[tuple[float, float]]:
    """Generate 2D atom coordinates and normalize minimum spacing to 1.

    NetworkX's spring layout uses larger edge weights as stronger attractions,
    so strongly coupled biological sites tend to be placed closer together.
    """
    if len(graph.nodes) < 2:
        raise ValueError("At least two sites are required.")

    pos = nx.spring_layout(graph, weight="weight", seed=seed, dim=2)

    # Preserve deterministic node order 0,1,2,...
    coords = np.array([pos[i] for i in sorted(graph.nodes)], dtype=float)

    # Move the geometric center to the origin. This is not physically required,
    # but it makes plots easier to read.
    coords -= coords.mean(axis=0, keepdims=True)

    distances = pdist(coords)
    min_distance = distances.min()

    if min_distance <= 0:
        raise ValueError("Generated layout contains overlapping atoms.")

    # QoolQit recommends using unit spacing for the closest atom pair in the
    # dimensionless model. After this scaling, min r_ij = 1.
    coords /= min_distance

    return [tuple(map(float, xy)) for xy in coords]


def couplings_to_qoolqit_coordinates(
    couplings: np.ndarray,
    seed: int = 42,
) -> list[tuple[float, float]]:
    """One-call helper: coupling matrix -> graph -> QoolQit coordinates."""
    graph = coupling_matrix_to_graph(couplings)
    return graph_to_qoolqit_coordinates(graph, seed=seed)


def rydberg_interaction_matrix(
    coordinates: list[tuple[float, float]],
) -> np.ndarray:
    """Compute the ideal dimensionless 1/r^6 interaction matrix ourselves.

    This is useful for learning and for comparing our geometry with the
    interaction matrix reported by QoolQit's Register class.
    """
    coords = np.asarray(coordinates, dtype=float)
    n = len(coords)
    matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            r = np.linalg.norm(coords[i] - coords[j])
            interaction = 1.0 / (r**6)
            matrix[i, j] = interaction
            matrix[j, i] = interaction

    return matrix
