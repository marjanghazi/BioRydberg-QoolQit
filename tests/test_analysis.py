import numpy as np
import pytest

from src.analysis import normalize_rows, split_profile_rmse, spatial_profile_rmse


def test_normalization_marks_low_occupation_as_undefined():
    values = np.array([[0.0, 0.0], [0.2, 0.3]])
    normalized, valid = normalize_rows(values, eps=1e-3, return_mask=True)
    assert not valid[0]
    assert np.isnan(normalized[0]).all()
    assert valid[1]
    assert np.allclose(normalized[1], [0.4, 0.6])


def test_profile_metric_excludes_invalid_time_points():
    bio = np.array([[1.0, 0.0], [0.4, 0.6], [0.2, 0.8]])
    ryd = np.array([[0.0, 0.0], [0.2, 0.3], [0.1, 0.4]])
    metrics = split_profile_rmse(bio, ryd, occupation_threshold=1e-3)
    assert metrics["excluded_time_points"] == 1
    assert metrics["valid_time_points"] == 2
    assert np.isfinite(metrics["rmse_all"])
    assert np.isfinite(spatial_profile_rmse(bio, ryd, occupation_threshold=1e-3))


def test_profile_metric_rejects_all_invalid():
    with pytest.raises(ValueError):
        spatial_profile_rmse(np.array([[1.0, 0.0]]), np.zeros((1, 2)))
