"""QoolQit simulation helpers for the neutral-atom side of BioRydberg.

This module contains the main QoolQit part of the project.

The overall workflow is:

    Atom Coordinates
          |
          v
       Register
          +
        Drive
          |
          v
    QuantumProgram
          |
          v
    compile_to(AnalogDevice)
          |
          v
     LocalEmulator
          |
          v
       Occupation
          |
          v
    Rydberg dynamics


IMPORTANT PHYSICS NOTE
======================

The biological excitation-transport Hamiltonian and QoolQit's standard
Rydberg Hamiltonian are NOT the same Hamiltonian.

The biological model contains excitation-hopping terms such as

    J_ij |i><j|

which transfer a single excitation between biological sites.

The QoolQit Rydberg model contains terms involving

    n_i n_j

together with a laser drive and detuning.

Therefore, the QoolQit simulation in this project should be interpreted as
a neutral-atom transport-like / excitation-spreading analogue.

We are NOT claiming that this is an exact simulation of the FMO Hamiltonian.

Later in the project we compare normalized spatial population patterns
between the biological reference model and the neutral-atom model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# ---------------------------------------------------------------------------
# QoolQit imports
# ---------------------------------------------------------------------------
#
# AnalogDevice
#     Represents the neutral-atom device model to which our abstract program
#     will be compiled.
#
# ConstantWaveform
#     Creates a control signal whose value stays constant during the sequence.
#
# Drive
#     Describes the laser control applied to the atoms.
#
# QuantumProgram
#     Combines:
#
#         Register + Drive
#
#     into one analog quantum program.
#
# Register
#     Stores the locations of our neutral atoms.
# ---------------------------------------------------------------------------

from qoolqit import (
    AnalogDevice,
    ConstantWaveform,
    Drive,
    QuantumProgram,
    Register,
)

# ---------------------------------------------------------------------------
# Execution-related QoolQit classes
# ---------------------------------------------------------------------------
#
# EmulationConfig
#     Tells the emulator what quantities we want to measure.
#
# LocalEmulator
#     Simulates the neutral-atom program locally on our computer.
#
# Occupation
#     Measures the Rydberg-state occupation <n_i> of every atom.
# ---------------------------------------------------------------------------

from qoolqit.execution import (
    EmulationConfig,
    LocalEmulator,
    Occupation,
)


# ===========================================================================
# DATA CONTAINER
# ===========================================================================


@dataclass
class QoolQitRun:
    """Store the important outputs of one QoolQit simulation.

    Attributes
    ----------
    relative_times:
        Times at which QoolQit actually stored the occupation observable.

        These are relative times:

            0.0 = beginning of sequence
            0.5 = halfway through
            1.0 = end of sequence

    occupations:
        Rydberg occupation of every atom at every stored time.

        Shape:

            (number_of_times, number_of_atoms)

        For example:

            occupations[5, 2]

        means:

            occupation of atom 3 at the sixth stored time.

    coordinates:
        Dimensionless 2D neutral-atom positions.

    interaction_matrix:
        QoolQit's Rydberg interaction matrix.

        Approximately:

            J_ij = 1 / r_ij^6

    program:
        The compiled QoolQit QuantumProgram.

        We keep it because later we may want to inspect:
        - register
        - sequence
        - compilation
        - device information
    """

    relative_times: np.ndarray
    occupations: np.ndarray
    coordinates: list[tuple[float, float]]
    interaction_matrix: np.ndarray
    program: QuantumProgram


# ===========================================================================
# BUILD THE QOOLQIT PROGRAM
# ===========================================================================


def build_program(
    coordinates: list[tuple[float, float]],
    duration: float = 10.0,
    amplitude: float = 0.8,
    detuning: float = 0.0,
) -> QuantumProgram:
    """Build and compile a QoolQit Rydberg analog program.

    Parameters
    ----------
    coordinates:
        Dimensionless 2D positions of our neutral atoms.

        Example:

            [
                (0.0, 0.0),
                (1.0, 0.0),
                (0.5, 1.0),
            ]

        The distances between these atoms determine their Rydberg
        interaction strengths.

    duration:
        Duration of the dimensionless analog program.

        Default:

            10.0

        This is NOT yet directly a physical time such as microseconds.
        QoolQit's compiler maps the dimensionless description to a device.

    amplitude:
        Global Rabi-drive amplitude Omega.

        The Rabi amplitude controls how strongly the laser couples:

            |g> <-> |r>

        where:

            |g> = atomic ground state
            |r> = Rydberg excited state

        Larger Omega generally means stronger/faster laser-driven dynamics.

    detuning:
        Global laser detuning delta.

        Detuning tells us how far the laser frequency is from the exact
        atomic transition frequency.

        Roughly:

            delta = 0

        means resonant driving.

    Returns
    -------
    QuantumProgram
        A compiled QoolQit QuantumProgram.

    Notes
    -----
    All parameters are written in QoolQit's dimensionless representation
    before compilation.
    """

    # -----------------------------------------------------------------------
    # STEP 1
    # Create the neutral-atom register.
    # -----------------------------------------------------------------------
    #
    # The Register describes WHERE our neutral atoms are located.
    #
    # This is crucial because the Rydberg interaction depends strongly
    # on distance:
    #
    #       J_ij ~ 1 / r_ij^6
    #
    # Therefore geometry is effectively part of the Hamiltonian.
    # -----------------------------------------------------------------------

    register = Register.from_coordinates(coordinates)

    # -----------------------------------------------------------------------
    # STEP 2
    # Create the Rabi-amplitude waveform.
    # -----------------------------------------------------------------------
    #
    # We start with the simplest possible control:
    #
    #       Omega(t) = constant
    #
    # during the entire experiment.
    #
    # Later we can experiment with:
    #
    #       ramps
    #       pulses
    #       optimized waveforms
    #
    # but a constant waveform is much easier to understand first.
    # -----------------------------------------------------------------------

    amplitude_waveform = ConstantWaveform(
        duration,
        amplitude,
    )

    # -----------------------------------------------------------------------
    # STEP 3
    # Create the detuning waveform.
    # -----------------------------------------------------------------------
    #
    # Again, we currently use:
    #
    #       delta(t) = constant
    #
    # throughout the simulation.
    # -----------------------------------------------------------------------

    detuning_waveform = ConstantWaveform(
        duration,
        detuning,
    )

    # -----------------------------------------------------------------------
    # STEP 4
    # Combine amplitude and detuning into one Drive.
    # -----------------------------------------------------------------------
    #
    # Think of Drive as:
    #
    #       "How do I control the neutral atoms with my laser?"
    #
    # -----------------------------------------------------------------------

    drive = Drive(
        amplitude=amplitude_waveform,
        detuning=detuning_waveform,
    )

    # -----------------------------------------------------------------------
    # STEP 5
    # Build the QuantumProgram.
    # -----------------------------------------------------------------------
    #
    # Conceptually:
    #
    #       Register
    #          +
    #        Drive
    #          =
    #     QuantumProgram
    #
    # The register describes the system.
    # The drive describes its time-dependent control.
    # -----------------------------------------------------------------------

    program = QuantumProgram(
        register=register,
        drive=drive,
    )

    # -----------------------------------------------------------------------
    # STEP 6
    # Compile the abstract program to the neutral-atom AnalogDevice.
    # -----------------------------------------------------------------------
    #
    # Before compilation our model is dimensionless.
    #
    # Compilation maps the abstract description onto device-compatible
    # parameters.
    #
    # "max_energy" asks QoolQit to scale the dimensionless model so it
    # effectively uses the available device energy scale.
    # -----------------------------------------------------------------------

    program.compile_to(
        AnalogDevice(),
        profile="max_energy",
    )

    return program


# ===========================================================================
# RUN THE QOOLQIT SIMULATION
# ===========================================================================


def run_occupation_dynamics(
    coordinates: list[tuple[float, float]],
    duration: float = 10.0,
    amplitude: float = 0.8,
    detuning: float = 0.0,
    n_times: int = 31,
) -> QoolQitRun:
    """Run QoolQit and obtain Rydberg occupation on every neutral atom.

    Parameters
    ----------
    coordinates:
        Neutral-atom positions.

    duration:
        Duration of the dimensionless quantum program.

    amplitude:
        Constant global Rabi-drive amplitude Omega.

    detuning:
        Constant global detuning delta.

    n_times:
        Number of observation times requested between:

            0.0 and 1.0

        These are relative sequence times.

    Returns
    -------
    QoolQitRun
        Object containing:

        - actual stored relative times
        - occupation dynamics
        - neutral-atom coordinates
        - Rydberg interaction matrix
        - compiled QoolQit program
    """

    # -----------------------------------------------------------------------
    # Basic validation
    # -----------------------------------------------------------------------

    if n_times < 2:
        raise ValueError(
            "n_times must be at least 2 so that the beginning and end "
            "of the sequence can both be observed."
        )

    if duration <= 0:
        raise ValueError("duration must be greater than zero.")

    if len(coordinates) < 2:
        raise ValueError(
            "At least two neutral atoms are required for this experiment."
        )

    # -----------------------------------------------------------------------
    # STEP 1
    # Build and compile the QoolQit program.
    # -----------------------------------------------------------------------

    program = build_program(
        coordinates=coordinates,
        duration=duration,
        amplitude=amplitude,
        detuning=detuning,
    )

    # -----------------------------------------------------------------------
    # STEP 2
    # Choose the times at which we WANT to measure occupation.
    # -----------------------------------------------------------------------
    #
    # np.linspace generates evenly-spaced numbers:
    #
    # Example with n_times = 5:
    #
    #       [0.00, 0.25, 0.50, 0.75, 1.00]
    #
    # These are RELATIVE times.
    #
    #       0.0 = beginning
    #       1.0 = end
    #
    # -----------------------------------------------------------------------

    requested_times = np.linspace(
        0.0,
        1.0,
        n_times,
    )

    # -----------------------------------------------------------------------
    # STEP 3
    # Define our observable.
    # -----------------------------------------------------------------------
    #
    # We want the Rydberg occupation:
    #
    #       <n_i>
    #
    # for every atom i.
    #
    # Very roughly:
    #
    #       <n_i> near 0
    #           -> atom is mostly not Rydberg excited
    #
    #       <n_i> near 1
    #           -> atom is mostly Rydberg excited
    #
    # -----------------------------------------------------------------------

    occupation_observable = Occupation(
        evaluation_times=requested_times.tolist(),
    )

    # -----------------------------------------------------------------------
    # STEP 4
    # Configure the emulator.
    # -----------------------------------------------------------------------
    #
    # We tell QoolQit:
    #
    #       "During the simulation, calculate the occupation observable."
    #
    # with_modulation=False means we initially use the ideal programmed
    # waveform.
    #
    # Later we can enable hardware modulation to investigate more realistic
    # device behavior.
    # -----------------------------------------------------------------------

    emulation_config = EmulationConfig(
        observables=(occupation_observable,),
        with_modulation=False,
    )

    # -----------------------------------------------------------------------
    # STEP 5
    # Create the local quantum emulator.
    # -----------------------------------------------------------------------

    emulator = LocalEmulator(
        emulation_config=emulation_config,
    )

    # -----------------------------------------------------------------------
    # STEP 6
    # Run the simulation.
    # -----------------------------------------------------------------------
    #
    # This is where the actual neutral-atom quantum dynamics are simulated.
    # -----------------------------------------------------------------------

    job = emulator.run(program)

    # -----------------------------------------------------------------------
    # STEP 7
    # Obtain the Results object.
    # -----------------------------------------------------------------------

    results = job.results()

    # -----------------------------------------------------------------------
    # STEP 8
    # Verify that the occupation observable exists.
    # -----------------------------------------------------------------------
    #
    # Usually the result tags should contain:
    #
    #       ["occupation"]
    #
    # This defensive check gives us a useful error message if the QoolQit
    # execution API changes or if the observable was not generated.
    # -----------------------------------------------------------------------

    result_tags = results.get_result_tags()

    if "occupation" not in result_tags:
        raise RuntimeError(
            "QoolQit did not return an 'occupation' observable.\n"
            f"Available result tags are: {result_tags}"
        )

    # -----------------------------------------------------------------------
    # STEP 9
    # Retrieve the ACTUAL times stored inside Pulser/QoolQit.
    # -----------------------------------------------------------------------
    #
    # VERY IMPORTANT:
    #
    # Previously we did this:
    #
    #       for t in requested_times:
    #           results.get_result("occupation", time=float(t))
    #
    # That caused the error:
    #
    #       'occupation' is not available at time
    #       0.43333333333333335
    #
    # Why?
    #
    # Floating-point numbers cannot always be stored exactly.
    #
    # For example, one calculation might produce:
    #
    #       0.43333333333333335
    #
    # while another part of the software stores:
    #
    #       0.4333333333333333
    #
    # Those represent essentially the same physical point, but Python's exact
    # list lookup treats them as different.
    #
    # Pulser's Results.get_result() performs an exact time lookup.
    #
    # Therefore we ask QoolQit/Pulser:
    #
    #       "Which times did YOU actually store?"
    #
    # -----------------------------------------------------------------------

    stored_times = np.asarray(
        results.get_result_times("occupation"),
        dtype=float,
    )

    # -----------------------------------------------------------------------
    # STEP 10
    # Retrieve all occupation values directly.
    # -----------------------------------------------------------------------
    #
    # get_tagged_results() returns something like:
    #
    # {
    #     "occupation": [
    #
    #         [atom1, atom2, atom3],   # first time
    #
    #         [atom1, atom2, atom3],   # second time
    #
    #         ...
    #     ]
    # }
    #
    # Therefore there is NO reason for us to repeatedly search using
    # floating-point times.
    #
    # This completely avoids the error you encountered.
    # -----------------------------------------------------------------------

    tagged_results = results.get_tagged_results()

    if "occupation" not in tagged_results:
        raise RuntimeError(
            "The result tags included 'occupation', but no occupation "
            "values were found in get_tagged_results()."
        )

    occupations = np.asarray(
        tagged_results["occupation"],
        dtype=float,
    )

    # -----------------------------------------------------------------------
    # STEP 11
    # Check that the returned data has the shape we expect.
    # -----------------------------------------------------------------------
    #
    # We expect:
    #
    #       rows    = stored times
    #       columns = atoms
    #
    # Example for 31 times and 3 atoms:
    #
    #       occupations.shape == (31, 3)
    #
    # -----------------------------------------------------------------------

    expected_atoms = len(coordinates)

    if occupations.ndim != 2:
        raise RuntimeError(
            "Unexpected occupation result shape.\n"
            f"Expected a 2D array, but received shape "
            f"{occupations.shape}."
        )

    if occupations.shape[0] != len(stored_times):
        raise RuntimeError(
            "Number of occupation time slices does not match the number "
            "of stored evaluation times.\n"
            f"Times: {len(stored_times)}\n"
            f"Occupation rows: {occupations.shape[0]}"
        )

    if occupations.shape[1] != expected_atoms:
        raise RuntimeError(
            "Number of occupation values per time step does not match "
            "the number of neutral atoms.\n"
            f"Atoms: {expected_atoms}\n"
            f"Occupation columns: {occupations.shape[1]}"
        )

    # -----------------------------------------------------------------------
    # STEP 12
    # Obtain QoolQit's neutral-atom interaction matrix.
    # -----------------------------------------------------------------------
    #
    # For our dimensionless Rydberg register:
    #
    #       J_ij ~ 1 / r_ij^6
    #
    # So if:
    #
    #       r_12 = 1
    #
    # then:
    #
    #       J_12 = 1
    #
    # while if:
    #
    #       r_12 = 2
    #
    # then:
    #
    #       J_12 = 1 / 64
    #
    # This extreme distance dependence is one of the central features
    # of Rydberg neutral-atom physics.
    # -----------------------------------------------------------------------

    interaction_matrix = np.asarray(
        program.register.interaction_matrix(),
        dtype=float,
    )

    # -----------------------------------------------------------------------
    # STEP 13
    # Package all simulation information into QoolQitRun.
    # -----------------------------------------------------------------------

    return QoolQitRun(
        relative_times=stored_times,
        occupations=occupations,
        coordinates=coordinates,
        interaction_matrix=interaction_matrix,
        program=program,
    )