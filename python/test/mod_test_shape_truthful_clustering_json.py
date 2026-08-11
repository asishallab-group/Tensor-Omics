"""
Python test suite for tox_stc_json (serialize_stc_results_as_json,
write_stc_interactive_html_report), mirroring test/mod_test_shape_truthful_clustering_json.F90
"""

import numpy as np
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import serialize_stc_results_as_json, write_stc_interactive_html_report
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_EMPTY_INPUT

STOP_REASON_FIXED_POINT = 4


def _fixture():
    """D=2, N=4, 2 seeds/ensembles: {1,2,3} (seed=1, d=1) and {2,3,4} (seed=4, d=0), which
    overlap on {2,3} -- Overlap Coefficient 2/3 -- and are merged into one super-ensemble."""
    vectors = np.asfortranarray([[0.0, 1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 0.0]], dtype=np.float64)
    dim_names = ['x', 'y']
    seed_selection_mask = np.array([True, False, False, True], dtype=np.bool_)

    ensemble_masks = np.zeros((4, 2), dtype=np.bool_, order='F')
    ensemble_masks[0:3, 0] = True
    ensemble_masks[1:4, 1] = True

    ensemble_stop_reason = np.array([STOP_REASON_FIXED_POINT, STOP_REASON_FIXED_POINT], dtype=np.int32)
    ensemble_growth_radii = np.array([1.0, 1.0], dtype=np.float64)

    ensemble_k_history = np.asfortranarray([[2, 2], [3, 3]], dtype=np.int32)
    ensemble_d_history = np.asfortranarray([[0, 0], [1, 0]], dtype=np.int32)
    ensemble_G_history = np.asfortranarray([[2.0, 2.0], [1.5, 1.5]], dtype=np.float64)

    ensemble_mu_history = np.zeros((2, 2, 2), dtype=np.float64, order='F')
    ensemble_mu_history[:, 0, 0] = [0.5, 0.0]
    ensemble_mu_history[:, 1, 0] = [1.0, 0.0]
    ensemble_mu_history[:, 0, 1] = [2.5, 0.0]
    ensemble_mu_history[:, 1, 1] = [3.0, 0.0]

    ensemble_S_history = np.zeros((2, 2, 2), dtype=np.float64, order='F')
    ensemble_S_history[0, 1, 0] = 0.5

    ensemble_U_history = np.zeros((2, 2, 2, 2), dtype=np.float64, order='F')
    ensemble_U_history[:, 0, 1, 0] = [1.0, 0.0]

    ensemble_low_confidence_masks = np.zeros((4, 2), dtype=np.bool_, order='F')
    ensemble_low_confidence_masks[0, 0] = True

    super_ensembles = np.asfortranarray([[1, 0], [2, 0]], dtype=np.int32)

    return (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
            ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
            super_ensembles)


def _common_kwargs():
    return dict(
        k_min=3, k_density=4, chordal_dist_max_as_prcnt_of_range=0.1, d_max=1,
        G_max=2.0, RMSE_change_max=0.5, f_max=0.8, a=3,
        exclusion_radius_percentile=50.0, bandwidth_percentile=68.0,
        reconciliation_mode='merge_overlap_coefficient', min_overlap_coefficient=0.5,
    )


def test_json_two_ensembles_with_overlap():
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
     super_ensembles) = _fixture()

    filename = "test_stc_two_ensembles_py.json"
    serialize_stc_results_as_json(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
        super_ensembles, **_common_kwargs())

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"dim_names":["x","y"]' in content
    assert '"n_vectors":4' in content
    assert '"n_dimensions":2' in content
    assert '"n_ensembles":2' in content
    assert '"k_min":3' in content
    assert '"reconciliation_mode":"merge_overlap_coefficient"' in content

    assert ('"id":1,"coords":[0.0000000000000000E+000,0.0000000000000000E+000],'
            '"n_ensembles":1,"n_low_confidence_ensembles":1,"ensembles":[1],'
            '"low_confidence_ensembles":[1],"seed_of":[1]') in content
    assert ('"id":2,"coords":[1.0000000000000000E+000,0.0000000000000000E+000],'
            '"n_ensembles":2,"n_low_confidence_ensembles":0,"ensembles":[1,2],'
            '"low_confidence_ensembles":[],"seed_of":[]') in content

    assert ('{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'
            '"size":3,"d":1,"G":1.5000000000000000E+000,"mu":[1.0000000000000000E+000,0.0000000000000000E+000],'
            '"u1":[1.0000000000000000E+000,0.0000000000000000E+000],"s1":5.0000000000000000E-001,'
            '"super_ensemble_id":1}') in content
    assert ('{"id":2,"seed_point_id":4,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'
            '"size":3,"d":0,"G":1.5000000000000000E+000,"mu":[3.0000000000000000E+000,0.0000000000000000E+000],'
            '"super_ensemble_id":1}') in content
    assert '"u2"' not in content

    assert '"super_ensembles":[{"group_id":1,"ensemble_ids":[1,2]}]' in content
    assert ('"overlap_coefficient_matrix":[{"a":1,"b":2,"overlap_coefficient":'
            '6.6666666666666663E-001}]') in content


def test_json_estimated_params_included():
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
     super_ensembles) = _fixture()

    filename = "test_stc_estimated_params_py.json"
    serialize_stc_results_as_json(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
        super_ensembles,
        estimated_k_min=5, estimated_k_density=6, estimated_density_quantile=0.75,
        estimated_chordal_dist_max_as_prcnt_of_range=0.2, estimated_G_max=3.0, estimated_d_max=2,
        **_common_kwargs())

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"estimated_k_min":5' in content
    assert '"estimated_k_density":6' in content
    assert '"estimated_d_max":2' in content


def test_json_zero_ensembles():
    """n_selected_seed=0, derived from an all-False seed_selection_mask, is a valid,
    well-defined 'no ensembles' input, not an error."""
    vectors = np.asfortranarray([[0.0, 1.0], [0.0, 0.0]], dtype=np.float64)
    dim_names = ['x', 'y']
    seed_selection_mask = np.array([False, False], dtype=np.bool_)
    ensemble_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    ensemble_stop_reason = np.zeros((0,), dtype=np.int32)
    ensemble_growth_radii = np.zeros((0,), dtype=np.float64)
    ensemble_U_history = np.zeros((2, 2, 1, 0), dtype=np.float64, order='F')
    ensemble_S_history = np.zeros((2, 1, 0), dtype=np.float64, order='F')
    ensemble_d_history = np.zeros((1, 0), dtype=np.int32, order='F')
    ensemble_G_history = np.zeros((1, 0), dtype=np.float64, order='F')
    ensemble_mu_history = np.zeros((2, 1, 0), dtype=np.float64, order='F')
    ensemble_k_history = np.zeros((1, 0), dtype=np.int32, order='F')
    ensemble_low_confidence_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    super_ensembles = np.zeros((2, 0), dtype=np.int32, order='F')

    filename = "test_stc_zero_ensembles_py.json"
    serialize_stc_results_as_json(
        filename, 0, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
        super_ensembles, **_common_kwargs())

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"ensembles":[]' in content
    assert '"super_ensembles":[]' in content
    assert '"overlap_coefficient_matrix":[]' in content


def test_html_report_wraps_json_in_template_and_d3():
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
     super_ensembles) = _fixture()

    filename = "test_stc_report_py.html"
    write_stc_interactive_html_report(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history, ensemble_low_confidence_masks,
        super_ensembles, **_common_kwargs())

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert content.startswith('<!DOCTYPE html>')
    assert 'd3js.org' in content
    assert 'const DATA = {"dim_names":["x","y"]' in content
    assert content.rstrip().endswith('</html>')


def test_json_invalid_n_dimensions():
    vectors = np.zeros((0, 2), dtype=np.float64, order='F')
    dim_names = []
    seed_selection_mask = np.array([False, False], dtype=np.bool_)
    ensemble_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    ensemble_stop_reason = np.zeros((0,), dtype=np.int32)
    ensemble_growth_radii = np.zeros((0,), dtype=np.float64)
    ensemble_U_history = np.zeros((0, 0, 1, 0), dtype=np.float64, order='F')
    ensemble_S_history = np.zeros((0, 1, 0), dtype=np.float64, order='F')
    ensemble_d_history = np.zeros((1, 0), dtype=np.int32, order='F')
    ensemble_G_history = np.zeros((1, 0), dtype=np.float64, order='F')
    ensemble_mu_history = np.zeros((0, 1, 0), dtype=np.float64, order='F')
    ensemble_k_history = np.zeros((1, 0), dtype=np.int32, order='F')
    ensemble_low_confidence_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    super_ensembles = np.zeros((2, 0), dtype=np.int32, order='F')

    assert_error(
        lambda: serialize_stc_results_as_json(
            'test_stc_invalid_py.json', 0, vectors, dim_names, seed_selection_mask, ensemble_masks,
            ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history,
            ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history,
            ensemble_low_confidence_masks, super_ensembles, **_common_kwargs()),
        "n_dimensions=0 must be rejected", ERR_EMPTY_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
