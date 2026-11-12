"""Exercise the generated STC API through the Python and C bindings."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tensor_omics import build_kd_index, tox_stc
from tensor_omics.error_handling import ToxInputError


def test_radius_variants():
    vectors = np.array([[0.0, 0.1, 0.2, 50.0, 55.0, 60.0]])
    dimension_order = np.array([1], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    labels = tox_stc.calculate_labels_as_density(vectors, 0.15, dimension_order, kd_indices, n_tiles=2)
    np.testing.assert_array_equal(labels, [2, 3, 2, 1, 1, 1])

    # Hold the density labels and the seeds fixed to isolate the growth radius.
    density_labels = np.array([8., 8., 8., 4., 4., 4.])
    seed_indices = np.array([1, 4], dtype=np.int32)
    inputs = (vectors, dimension_order, kd_indices, density_labels, seed_indices)
    options = dict(t_observables=10, n_tiles=2)
    global_result = tox_stc.obtain_ensembles_global_radius(*inputs, r=0.15, **options)
    local_result = tox_stc.obtain_ensembles_seed_radius(*inputs, seed_radii=np.array([0.15, 5.5]), **options)
    assert "n_ensembles" not in global_result
    assert "n_ensembles" not in local_result
    np.testing.assert_array_equal(global_result["ensemble_matrix"].sum(axis=0), [3, 1])
    np.testing.assert_array_equal(local_result["ensemble_matrix"].sum(axis=0), [3, 3])
    np.testing.assert_array_equal(global_result["stop_reasons"], [1, 4])
    np.testing.assert_array_equal(local_result["stop_reasons"], [1, 1])
    assert global_result["mad_ambient"] == local_result["mad_ambient"] == 2.0

    # The generator trims the merged output to the number of surviving columns.
    merged = tox_stc.merge_ensembles(global_result["ensemble_matrix"])
    np.testing.assert_array_equal(merged, global_result["ensemble_matrix"][:, :1])


def test_invalid_observable_history_rows():
    mask = np.array([True, False])
    labels = np.array([2., 1.])
    history = np.zeros((4, 2), dtype=np.float64, order="F")
    with np.testing.assert_raises(ToxInputError):
        tox_stc.compute_ensemble_observable(mask, labels, 0, 1, history)
    with np.testing.assert_raises(ToxInputError):
        tox_stc.accept_ensemble(history, 2)


if __name__ == "__main__":
    test_radius_variants()
    test_invalid_observable_history_rows()
