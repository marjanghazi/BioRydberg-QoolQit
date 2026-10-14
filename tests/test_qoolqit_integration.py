import numpy as np

from src.qoolqit_simulation import run_occupation_dynamics


def test_local_emulator_executes_compiled_qoolqit_program():
    run = run_occupation_dynamics(
        [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)],
        duration=2.0,
        amplitude=0.5,
        detuning=0.0,
        n_times=5,
        initial_site=0,
    )
    assert run.program.is_compiled
    assert run.occupations.shape == (5, 3)
    assert np.allclose(run.occupations[0], [1.0, 0.0, 0.0], atol=1e-8)
    assert np.allclose(run.interaction_matrix, run.interaction_matrix.T)
    assert run.qoolqit_version == "1.4.0"
