\
"""Reference models for biological excitation-energy transport.

This module implements the *biological* side of the project.

We use the single-excitation basis:

    |0> = excitation on pigment/site 0
    |1> = excitation on pigment/site 1
    ...

In this basis, an excitonic Hamiltonian is represented by an N x N matrix:

                 [ E0   J01  J02 ... ]
    H_bio   =     [ J01  E1   J12 ... ]
                 [ J02  J12  E2  ... ]
                 [ ...              ]

where:
    E_i  = site energy of pigment i
    J_ij = excitation hopping coupling between pigments i and j

The evolution is unitary:

    |psi(t)> = exp(-i H t) |psi(0)>

and the probability of finding the excitation at site i is:

    P_i(t) = |psi_i(t)|^2

We use dimensionless/arbitrary units here. That is deliberate: first we want
to understand the transport physics before introducing experimental units.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import expm


def build_hamiltonian(
    site_energies: np.ndarray,
    couplings: np.ndarray,
) -> np.ndarray:
    """Build an excitonic Hamiltonian.

    Parameters
    ----------
    site_energies:
        1D array with one energy E_i per site.

    couplings:
        Symmetric N x N matrix. couplings[i, j] is J_ij.
        The diagonal should normally be zero.

    Returns
    -------
    numpy.ndarray
        Hermitian N x N Hamiltonian.

    Why do we validate the matrix?
    ------------------------------
    Quantum Hamiltonians must be Hermitian. With real-valued site energies
    and couplings, that means the matrix must be symmetric.
    """
    site_energies = np.asarray(site_energies, dtype=float)
    couplings = np.asarray(couplings, dtype=float)

    n = len(site_energies)

    if couplings.shape != (n, n):
        raise ValueError(
            f"Coupling matrix must have shape {(n, n)}, "
            f"but got {couplings.shape}."
        )

    if not np.allclose(couplings, couplings.T):
        raise ValueError("Coupling matrix must be symmetric.")

    # Put site energies on the diagonal.
    H = np.diag(site_energies)

    # Add off-diagonal hopping couplings.
    H = H + couplings

    if not np.allclose(H, H.T.conj()):
        raise ValueError("Hamiltonian is not Hermitian.")

    return H


def localized_initial_state(n_sites: int, initial_site: int = 0) -> np.ndarray:
    """Create |psi(0)> with one excitation localized at one site."""
    if not 0 <= initial_site < n_sites:
        raise ValueError("initial_site is outside the valid site range.")

    psi0 = np.zeros(n_sites, dtype=complex)
    psi0[initial_site] = 1.0 + 0.0j
    return psi0


def simulate_unitary_transport(
    hamiltonian: np.ndarray,
    times: np.ndarray,
    initial_site: int = 0,
) -> np.ndarray:
    """Simulate excitation populations P_i(t).

    Parameters
    ----------
    hamiltonian:
        N x N excitonic Hamiltonian.

    times:
        1D array of simulation times.

    initial_site:
        Site where the excitation begins.

    Returns
    -------
    numpy.ndarray
        Shape = (number_of_times, number_of_sites).
        Entry [k, i] is the population of site i at times[k].
    """
    H = np.asarray(hamiltonian, dtype=complex)
    times = np.asarray(times, dtype=float)

    n_sites = H.shape[0]
    psi0 = localized_initial_state(n_sites, initial_site)

    populations = np.zeros((len(times), n_sites), dtype=float)

    for k, t in enumerate(times):
        # U(t) = exp(-i H t), because we choose hbar = 1.
        U = expm(-1j * H * t)

        # Evolve the wavefunction.
        psi_t = U @ psi0

        # Born rule: probability = |amplitude|^2.
        populations[k] = np.abs(psi_t) ** 2

    return populations


def two_site_example() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a simple two-site model used for the first lesson.

    We choose equal site energies and one coupling J=1.
    This produces clear oscillatory transfer between the two sites.
    """
    site_energies = np.array([0.0, 0.0])

    couplings = np.array(
        [
            [0.0, 1.0],
            [1.0, 0.0],
        ]
    )

    H = build_hamiltonian(site_energies, couplings)
    times = np.linspace(0.0, 10.0, 300)
    populations = simulate_unitary_transport(H, times, initial_site=0)

    return H, times, populations


def three_site_example() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a small biological-network-inspired three-site model.

    Site 0 is our initial donor-like site.
    Site 2 acts as a target/acceptor-like site.

    The couplings are deliberately unequal, so the network is not trivial.
    """
    site_energies = np.array([0.00, 0.20, -0.10])

    couplings = np.array(
        [
            [0.00, 1.00, 0.25],
            [1.00, 0.00, 0.70],
            [0.25, 0.70, 0.00],
        ]
    )

    H = build_hamiltonian(site_energies, couplings)
    times = np.linspace(0.0, 12.0, 300)
    populations = simulate_unitary_transport(H, times, initial_site=0)

    return H, times, populations
