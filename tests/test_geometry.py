import numpy as np
import pytest
from scipy.spatial.distance import pdist

from src.neutral_atom_mapping import rydberg_interaction_matrix, triangle_coordinates, validate_coordinates


def test_triangle_has_unit_minimum_separation_and_symmetric_interactions():
    coordinates = triangle_coordinates(0.4, 1.2)
    assert np.isclose(np.min(pdist(np.asarray(coordinates))), 1.0)
    validate_coordinates(coordinates)
    matrix = rydberg_interaction_matrix(coordinates)
    assert np.allclose(matrix, matrix.T)
    assert np.allclose(np.diag(matrix), 0.0)


def test_overlapping_atoms_are_rejected():
    with pytest.raises(ValueError):
        validate_coordinates([(0.0, 0.0), (0.0, 0.0)])
