import numpy as np
import pytest

from src.biological_model import build_hamiltonian, localized_initial_state, three_site_example, two_site_example


@pytest.mark.parametrize("factory", [two_site_example, three_site_example])
def test_reference_models_are_hermitian_and_conserve_probability(factory):
    hamiltonian, _, populations = factory()
    assert np.allclose(hamiltonian, hamiltonian.T.conj())
    assert np.max(np.abs(populations.sum(axis=1) - 1.0)) < 1e-12


def test_hamiltonian_rejects_asymmetric_couplings():
    with pytest.raises(ValueError, match="Hermitian"):
        build_hamiltonian(np.zeros(2), np.array([[0.0, 1.0], [0.0, 0.0]]))


def test_localized_state_validation():
    assert np.array_equal(localized_initial_state(3, 1), np.array([0.0, 1.0, 0.0]))
    with pytest.raises(ValueError):
        localized_initial_state(3, 3)
