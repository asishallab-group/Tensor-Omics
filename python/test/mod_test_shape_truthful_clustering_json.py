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
STOP_REASON_REJECTED_AFTER_STABLE = 2


def _fixture():
    """D=2, N=4, 2 seeds/ensembles: {1,2,3} (seed=1, d=1) and {2,3,4} (seed=4, d=0), which
    overlap on {2,3} -- Overlap Coefficient 2/3 -- and are merged into one super-ensemble.
    Both ensembles are fully accepted at both retained history columns -- the
    rejected-trailing-column case gets its own dedicated fixture below. Every S_history entry
    that feeds an observable_history rmse is chosen so normal_error is a perfect square (4.0 or
    9.0), avoiding any rounding-tie risk from the +epsilon guard around a normal_error of
    exactly 0."""
    vectors = np.asfortranarray([[0.0, 1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 0.0]], dtype=np.float64)
    dim_names = ['x', 'y']
    seed_selection_mask = np.array([True, False, False, True], dtype=np.bool_)

    ensemble_masks = np.zeros((4, 2), dtype=np.bool_, order='F')
    ensemble_masks[0:3, 0] = True
    ensemble_masks[1:4, 1] = True

    ensemble_stop_reason = np.array([STOP_REASON_FIXED_POINT, STOP_REASON_FIXED_POINT], dtype=np.int32)
    ensemble_growth_radii = np.array([1.0, 1.0], dtype=np.float64)

    # k-1=1 at every retained column of both ensembles, to keep the eigenvalue arithmetic
    # (S**2/(k-1)) below simple; k_history's exact value is otherwise unasserted.
    ensemble_k_history = np.asfortranarray([[2, 2], [2, 2]], dtype=np.int32)
    ensemble_d_history = np.asfortranarray([[0, 0], [1, 0]], dtype=np.int32)
    ensemble_G_history = np.asfortranarray([[2.0, 2.0], [1.5, 1.5]], dtype=np.float64)

    ensemble_mu_history = np.zeros((2, 2, 2), dtype=np.float64, order='F')
    ensemble_mu_history[:, 0, 0] = [0.5, 0.0]
    ensemble_mu_history[:, 1, 0] = [1.0, 0.0]
    ensemble_mu_history[:, 0, 1] = [2.5, 0.0]
    ensemble_mu_history[:, 1, 1] = [3.0, 0.0]

    ensemble_S_history = np.zeros((2, 2, 2), dtype=np.float64, order='F')
    ensemble_S_history[0, 0, 0] = 2.0  # ensemble 1, iter 1 (d=0): normal_error=4.0, rmse=2.0
    ensemble_S_history[0, 1, 0] = 0.5  # ensemble 1, iter 2 (d=1): s1, excluded from normal_error
    ensemble_S_history[1, 1, 0] = 3.0  # ensemble 1, iter 2 (d=1): normal_error=9.0, rmse=3.0
    ensemble_S_history[0, 0, 1] = 2.0  # ensemble 2, iter 1 (d=0): normal_error=4.0, rmse=2.0
    ensemble_S_history[0, 1, 1] = 3.0  # ensemble 2, iter 2 (d=0): normal_error=9.0, rmse=3.0

    ensemble_U_history = np.zeros((2, 2, 2, 2), dtype=np.float64, order='F')
    ensemble_U_history[:, 0, 1, 0] = [1.0, 0.0]  # ensemble 1's u1 (tangent)
    ensemble_U_history[:, 1, 1, 0] = [0.0, 1.0]  # ensemble 1's normal direction (residual_length)
    ensemble_U_history[:, 0, 1, 1] = [1.0, 0.0]  # ensemble 2 has d=0: only ever read by residual_length
    ensemble_U_history[:, 1, 1, 1] = [0.0, 1.0]

    # Both ensembles accepted at every retained column -- no rejection in this fixture.
    ensemble_accepted_history = np.ones((2, 2), dtype=np.bool_, order='F')

    # Ensemble 1: seed = point 1; point 2 joins at iteration 1, point 3 at iteration 2; point 4
    # never joins. T (max) = 2, matching idx=2.
    # Ensemble 2: seed = point 4; point 3 joins at iteration 1, point 2 at iteration 2; point 1
    # never joins. T (max) = 2.
    ensemble_member_added_at_step = np.asfortranarray([[0, -1], [1, 2], [2, 1], [-1, 0]], dtype=np.int32)

    # Bootstrap (iteration 1) basis, duplicating history column 1 exactly -- both d=0 here, so
    # U_first is never actually read (stc_chordal_distance requires d>0 to be "applicable").
    ensemble_U_first = np.zeros((2, 2, 2), dtype=np.float64, order='F')
    ensemble_U_first[:, :, 0] = ensemble_U_history[:, :, 0, 0]
    ensemble_U_first[:, :, 1] = ensemble_U_history[:, :, 0, 1]
    ensemble_d_first = np.array([ensemble_d_history[0, 0], ensemble_d_history[0, 1]], dtype=np.int32)

    ensemble_low_confidence_masks = np.zeros((4, 2), dtype=np.bool_, order='F')
    ensemble_low_confidence_masks[0, 0] = True

    super_ensembles = np.asfortranarray([[1, 0], [2, 0]], dtype=np.int32)

    return (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
            ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
            ensemble_G_history, ensemble_mu_history, ensemble_k_history,
            ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
            ensemble_U_first, ensemble_d_first, super_ensembles)


def _common_kwargs():
    return dict(
        k_min=3, k_density=4, chordal_dist_max_as_prcnt_of_range=0.1, d_max=1,
        G_max=2.0, RMSE_change_max=0.5, f_max=0.8, a=3,
        exclusion_radius_percentile=50.0, bandwidth_percentile=68.0,
        reconciliation_mode='merge_overlap_coefficient', min_overlap_coefficient=0.5,
    )


def _all_eligible_kwargs(n_selected_seed):
    """The new required per-ensemble reconciliation-eligibility arrays, all-eligible unless a
    test overrides them -- mirrors test/mod_test_shape_truthful_clustering_json.F90's own
    per-test `ensemble_eligible = .true.` convention."""
    return dict(
        ensemble_eligible=np.ones(n_selected_seed, dtype=np.bool_),
        ensemble_eligible_by_stop_condition=np.ones(n_selected_seed, dtype=np.bool_),
        ensemble_eligible_by_dimension=np.ones(n_selected_seed, dtype=np.bool_),
        ensemble_eligible_by_var_explained=np.ones(n_selected_seed, dtype=np.bool_),
    )


def test_json_two_ensembles_with_overlap():
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history,
     ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
     ensemble_U_first, ensemble_d_first, super_ensembles) = _fixture()

    filename = "test_stc_two_ensembles_py.json"
    serialize_stc_results_as_json(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(2))

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"dim_names":["x","y"]' in content
    assert '"n_vectors":4' in content
    assert '"n_dimensions":2' in content
    assert '"n_ensembles":2' in content
    assert '"k_min":3' in content
    assert '"reconciliation_mode":"merge_overlap_coefficient"' in content

    # residual_length is the point's distance off each ensemble's *normal* subspace at its
    # final retained basis (identity here, see _fixture); point 1 sits exactly at ensemble 1's
    # mean (residual 0); point 2 sits at ensemble 1's mean too (residual 0) but 2.0 off
    # ensemble 2's mean along its own (identity) normal directions.
    assert ('"id":1,"coords":[0.0000000000000000E+000,0.0000000000000000E+000],'
            '"n_ensembles":1,"n_low_confidence_ensembles":1,'
            '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000}],'
            '"low_confidence_ensembles":[1],"seed_of":[1]') in content
    assert ('"id":2,"coords":[1.0000000000000000E+000,0.0000000000000000E+000],'
            '"n_ensembles":2,"n_low_confidence_ensembles":0,'
            '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000},'
            '{"id":2,"residual_length":2.0000000000000000E+000}],'
            '"low_confidence_ensembles":[],"seed_of":[]') in content

    # ensembles: #1 has d=1 (u1/s1/line_start/line_end present), #2 has d=0 (no tangent
    # direction or line at all). Neither ensemble's drift is ever computable here (ensemble 1's
    # iter1-vs-iter2 d mismatches, 0 vs 1; ensemble 2's two iterations are both d=0, and
    # stc_chordal_distance requires d>0 to be "applicable" at all). Both ensembles are fully
    # eligible (see _all_eligible_kwargs), so reconciliation_eligible/excluded_by close out
    # each ensemble object.
    assert ('{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'
            '"size":3,"t_final":2,"d":1,"G":1.5000000000000000E+000,'
            '"mu":[1.0000000000000000E+000,0.0000000000000000E+000],'
            '"u1":[1.0000000000000000E+000,0.0000000000000000E+000],"s1":5.0000000000000000E-001,'
            '"line_start":[0.0000000000000000E+000,0.0000000000000000E+000],'
            '"line_end":[2.0000000000000000E+000,0.0000000000000000E+000],'
            '"observable_history":[{"iteration":1,"g":2.0000000000000000E+000,'
            '"rmse":2.0000000000000000E+000,"drift":null},{"iteration":2,"g":1.5000000000000000E+000,'
            '"rmse":3.0000000000000000E+000,"drift":null}],'
            '"super_ensemble_id":1,"final_chordal_distance":null,'
            '"reconciliation_eligible":true,"excluded_by":[]}') in content
    assert ('{"id":2,"seed_point_id":4,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'
            '"size":3,"t_final":2,"d":0,"G":1.5000000000000000E+000,'
            '"mu":[3.0000000000000000E+000,0.0000000000000000E+000],'
            '"observable_history":[{"iteration":1,"g":2.0000000000000000E+000,'
            '"rmse":2.0000000000000000E+000,"drift":null},{"iteration":2,"g":1.5000000000000000E+000,'
            '"rmse":3.0000000000000000E+000,"drift":null}],'
            '"super_ensemble_id":1,"final_chordal_distance":null,'
            '"reconciliation_eligible":true,"excluded_by":[]}') in content
    assert '"u2"' not in content

    assert '"super_ensembles":[{"group_id":1,"ensemble_ids":[1,2]}]' in content
    assert ('"overlap_coefficient_matrix":[{"a":1,"b":2,"overlap_coefficient":'
            '6.6666666666666663E-001}]') in content


def test_json_estimated_params_included():
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history,
     ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
     ensemble_U_first, ensemble_d_first, super_ensembles) = _fixture()

    filename = "test_stc_estimated_params_py.json"
    serialize_stc_results_as_json(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles,
        estimated_k_min=5, estimated_k_density=6, estimated_density_quantile=0.75,
        estimated_chordal_dist_max_as_prcnt_of_range=0.2, estimated_G_max=3.0, estimated_d_max=2,
        **_common_kwargs(), **_all_eligible_kwargs(2))

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
    ensemble_accepted_history = np.zeros((1, 0), dtype=np.bool_, order='F')
    ensemble_member_added_at_step = np.zeros((2, 0), dtype=np.int32, order='F')
    ensemble_low_confidence_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    ensemble_U_first = np.zeros((2, 2, 0), dtype=np.float64, order='F')
    ensemble_d_first = np.zeros((0,), dtype=np.int32)
    super_ensembles = np.zeros((2, 0), dtype=np.int32, order='F')

    filename = "test_stc_zero_ensembles_py.json"
    serialize_stc_results_as_json(
        filename, 0, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(0))

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"ensembles":[]' in content
    assert '"super_ensembles":[]' in content
    assert '"overlap_coefficient_matrix":[]' in content


def test_json_no_history_ensemble_omits_observable_keys():
    """An ensemble whose k_history is entirely zero never produced an observable at all (only
    possible for STOP_REASON_MAX_SIZE firing at the bootstrap step itself). Its
    t_final/d/G/mu/tangent/observable_history keys must be omitted, not merely null --
    final_chordal_distance is the one exception, always present (null here). The one point
    that is a member of this ensemble still gets its usual `ensembles` entry, with a
    residual_length of 0.0 (the documented fallback when no retained basis exists at all)."""
    vectors = np.asfortranarray([[0.0, 1.0], [0.0, 0.0]], dtype=np.float64)
    dim_names = ['x', 'y']
    seed_selection_mask = np.array([True, False], dtype=np.bool_)
    ensemble_masks = np.asfortranarray([[True], [False]], dtype=np.bool_)
    ensemble_stop_reason = np.array([STOP_REASON_FIXED_POINT], dtype=np.int32)
    ensemble_growth_radii = np.array([1.0], dtype=np.float64)
    ensemble_k_history = np.asfortranarray([[0]], dtype=np.int32)
    ensemble_d_history = np.zeros((1, 1), dtype=np.int32, order='F')
    ensemble_G_history = np.zeros((1, 1), dtype=np.float64, order='F')
    ensemble_mu_history = np.zeros((2, 1, 1), dtype=np.float64, order='F')
    ensemble_S_history = np.zeros((2, 1, 1), dtype=np.float64, order='F')
    ensemble_U_history = np.zeros((2, 2, 1, 1), dtype=np.float64, order='F')
    ensemble_accepted_history = np.ones((1, 1), dtype=np.bool_, order='F')
    ensemble_member_added_at_step = np.asfortranarray([[0], [-1]], dtype=np.int32)
    ensemble_U_first = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    ensemble_d_first = np.array([0], dtype=np.int32)
    ensemble_low_confidence_masks = np.zeros((2, 1), dtype=np.bool_, order='F')
    super_ensembles = np.zeros((2, 0), dtype=np.int32, order='F')

    filename = "test_stc_no_history_py.json"
    serialize_stc_results_as_json(
        filename, 0, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(1))

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert ('{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'
            '"size":1,"super_ensemble_id":null,"final_chordal_distance":null,'
            '"reconciliation_eligible":true,"excluded_by":[]}') in content
    assert '"t_final"' not in content
    assert '"d":' not in content
    assert '"G":' not in content
    assert '"mu":' not in content
    assert '"line_start"' not in content
    assert '"observable_history"' not in content
    assert '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000}]' in content


def test_html_report_wraps_json_in_template_and_d3():
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history,
     ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
     ensemble_U_first, ensemble_d_first, super_ensembles) = _fixture()

    filename = "test_stc_report_py.html"
    write_stc_interactive_html_report(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(2))

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
    ensemble_accepted_history = np.zeros((1, 0), dtype=np.bool_, order='F')
    ensemble_member_added_at_step = np.zeros((2, 0), dtype=np.int32, order='F')
    ensemble_low_confidence_masks = np.zeros((2, 0), dtype=np.bool_, order='F')
    ensemble_U_first = np.zeros((0, 0, 0), dtype=np.float64, order='F')
    ensemble_d_first = np.zeros((0,), dtype=np.int32)
    super_ensembles = np.zeros((2, 0), dtype=np.int32, order='F')

    assert_error(
        lambda: serialize_stc_results_as_json(
            'test_stc_invalid_py.json', 0, vectors, dim_names, seed_selection_mask, ensemble_masks,
            ensemble_stop_reason, ensemble_growth_radii, ensemble_U_history, ensemble_S_history,
            ensemble_d_history, ensemble_G_history, ensemble_mu_history, ensemble_k_history,
            ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
            ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(0)),
        "n_dimensions=0 must be rejected", ERR_EMPTY_INPUT)


def test_json_rejected_trailing_column_uses_last_accepted():
    """Regression test for a real, pre-existing bug: stc_push_ensemble_history also pushes a
    *rejected* final candidate into the trailing history window right before growth halts on
    STOP_REASON_REJECTED_IMMEDIATELY/STOP_REASON_REJECTED_AFTER_STABLE -- so the last
    *populated* history column is not always the ensemble's real last *accepted* state. This
    fixture's column 1 is the true accepted state (d=1, G=1.0, mu=[0.5,0.0], u1=[1,0]); its
    column 2 is a rejected candidate with deliberately different, wrong-if-leaked values
    (d=0, G=99.0, mu=[9.0,9.0], accepted_history=False). Every reported "current state" field
    must reflect column 1, and observable_history must contain exactly that one entry."""
    vectors = np.asfortranarray([[0.0, 1.0, 9.0], [0.0, 0.0, 9.0]], dtype=np.float64)
    dim_names = ['x', 'y']
    seed_selection_mask = np.array([True, False, False], dtype=np.bool_)
    ensemble_masks = np.asfortranarray([[True], [True], [False]], dtype=np.bool_)
    ensemble_stop_reason = np.array([STOP_REASON_REJECTED_AFTER_STABLE], dtype=np.int32)
    ensemble_growth_radii = np.array([1.0], dtype=np.float64)

    ensemble_k_history = np.asfortranarray([[2], [3]], dtype=np.int32)
    ensemble_d_history = np.asfortranarray([[1], [0]], dtype=np.int32)
    ensemble_G_history = np.asfortranarray([[1.0], [99.0]], dtype=np.float64)
    ensemble_mu_history = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    ensemble_mu_history[:, 0, 0] = [0.5, 0.0]
    ensemble_mu_history[:, 1, 0] = [9.0, 9.0]
    ensemble_S_history = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    ensemble_S_history[:, 0, 0] = [0.5, 2.0]  # eigen(2)=4.0/(2-1)=4.0 -> rmse=2.0
    ensemble_S_history[:, 1, 0] = [7.0, 7.0]
    ensemble_U_history = np.zeros((2, 2, 2, 1), dtype=np.float64, order='F')
    ensemble_U_history[:, 0, 0, 0] = [1.0, 0.0]
    ensemble_U_history[:, 1, 0, 0] = [0.0, 1.0]
    ensemble_accepted_history = np.asfortranarray([[True], [False]], dtype=np.bool_)

    # T=1: only the bootstrap was ever actually accepted.
    ensemble_member_added_at_step = np.asfortranarray([[0], [1], [-1]], dtype=np.int32)
    ensemble_U_first = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    ensemble_U_first[:, :, 0] = ensemble_U_history[:, :, 0, 0]
    ensemble_d_first = np.array([1], dtype=np.int32)
    ensemble_low_confidence_masks = np.zeros((3, 1), dtype=np.bool_, order='F')
    super_ensembles = np.zeros((1, 0), dtype=np.int32, order='F')

    filename = "test_stc_rejected_trailing_py.json"
    serialize_stc_results_as_json(
        filename, 0, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(1))

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"stop_reason":"rejected_after_stable"' in content
    assert '"size":2' in content
    assert '"t_final":1' in content  # the rejected step never counted
    assert '"d":1' in content  # from the accepted column, not the rejected column's 0
    assert '"G":1.0000000000000000E+000' in content  # not the rejected column's 99.0
    assert '"mu":[5.0000000000000000E-001,0.0000000000000000E+000]' in content  # not [9,9]
    assert '"u1":[1.0000000000000000E+000,0.0000000000000000E+000]' in content  # rejected column had d=0
    assert ('"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,'
            '"rmse":2.0000000000000000E+000,"drift":null}]') in content
    assert '9.9000000000000000E+001' not in content  # rejected G=99.0 never appears anywhere
    # Point 3's own coords legitimately contain [9,9], so check the ensemble's mu specifically.
    assert '"mu":[9.0000000000000000E+000,9.0000000000000000E+000]' not in content


def test_json_drift_and_final_chordal_distance():
    """A 3-retained-iteration, d=1-throughout ensemble whose tangent direction rotates 90
    degrees, then back: u1 = [1,0] (iter 1 = bootstrap = U_first), [0,1] (iter 2), [1,0]
    (iter 3). Every consecutive pair is orthogonal, so the chordal distance
    sqrt(1 - cos(theta)**2) is exactly 1.0 for both drift entries. final_chordal_distance
    compares iteration 3 against {U_first, iter1, iter2}: 0.0, 0.0, 1.0 -- max = 1.0."""
    vectors = np.asfortranarray([[0.0, 5.0], [0.0, 5.0]], dtype=np.float64)
    dim_names = ['x', 'y']
    seed_selection_mask = np.array([True, False], dtype=np.bool_)
    ensemble_masks = np.asfortranarray([[True], [True]], dtype=np.bool_)
    ensemble_stop_reason = np.array([STOP_REASON_FIXED_POINT], dtype=np.int32)
    ensemble_growth_radii = np.array([1.0], dtype=np.float64)

    ensemble_k_history = np.asfortranarray([[2], [2], [2]], dtype=np.int32)
    ensemble_d_history = np.asfortranarray([[1], [1], [1]], dtype=np.int32)
    ensemble_G_history = np.asfortranarray([[1.0], [2.0], [3.0]], dtype=np.float64)
    ensemble_mu_history = np.zeros((2, 3, 1), dtype=np.float64, order='F')

    # s1=1.0 throughout (unasserted); the second singular value drives rmse via eigen(2) (d=1
    # excludes eigen(1)): 2**2/1=4 -> rmse=2, 3**2/1=9 -> rmse=3, 4**2/1=16 -> rmse=4.
    ensemble_S_history = np.zeros((2, 3, 1), dtype=np.float64, order='F')
    ensemble_S_history[:, 0, 0] = [1.0, 2.0]
    ensemble_S_history[:, 1, 0] = [1.0, 3.0]
    ensemble_S_history[:, 2, 0] = [1.0, 4.0]

    ensemble_U_history = np.zeros((2, 2, 3, 1), dtype=np.float64, order='F')
    ensemble_U_history[:, 0, 0, 0] = [1.0, 0.0]
    ensemble_U_history[:, 0, 1, 0] = [0.0, 1.0]
    ensemble_U_history[:, 0, 2, 0] = [1.0, 0.0]

    ensemble_accepted_history = np.ones((3, 1), dtype=np.bool_, order='F')
    ensemble_member_added_at_step = np.asfortranarray([[0], [3]], dtype=np.int32)  # T=3

    ensemble_U_first = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    ensemble_U_first[:, 0, 0] = [1.0, 0.0]  # matches iteration 1 exactly
    ensemble_d_first = np.array([1], dtype=np.int32)
    ensemble_low_confidence_masks = np.zeros((2, 1), dtype=np.bool_, order='F')
    super_ensembles = np.zeros((1, 0), dtype=np.int32, order='F')

    filename = "test_stc_drift_py.json"
    serialize_stc_results_as_json(
        filename, 0, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(1))

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"t_final":3' in content
    assert ('"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,'
            '"rmse":2.0000000000000000E+000,"drift":null},'
            '{"iteration":2,"g":2.0000000000000000E+000,"rmse":3.0000000000000000E+000,'
            '"drift":1.0000000000000000E+000},'
            '{"iteration":3,"g":3.0000000000000000E+000,"rmse":4.0000000000000000E+000,'
            '"drift":1.0000000000000000E+000}]') in content
    assert '"final_chordal_distance":1.0000000000000000E+000' in content


def test_json_observable_history_rmse_null_when_size_le_one():
    """observable_history's rmse key is present but null when that column's ensemble size was
    <=1 (k_history entry of 1 -- normal_error/(k-1) would divide by zero)."""
    vectors = np.asfortranarray([[0.0], [0.0]], dtype=np.float64)
    dim_names = ['x', 'y']
    seed_selection_mask = np.array([True], dtype=np.bool_)
    ensemble_masks = np.asfortranarray([[True]], dtype=np.bool_)
    ensemble_stop_reason = np.array([STOP_REASON_FIXED_POINT], dtype=np.int32)
    ensemble_growth_radii = np.array([1.0], dtype=np.float64)
    ensemble_k_history = np.asfortranarray([[1]], dtype=np.int32)  # degenerate: k-1=0
    ensemble_d_history = np.zeros((1, 1), dtype=np.int32, order='F')
    ensemble_G_history = np.asfortranarray([[1.0]], dtype=np.float64)
    ensemble_mu_history = np.zeros((2, 1, 1), dtype=np.float64, order='F')
    ensemble_S_history = np.ones((2, 1, 1), dtype=np.float64, order='F')
    ensemble_U_history = np.zeros((2, 2, 1, 1), dtype=np.float64, order='F')
    ensemble_accepted_history = np.ones((1, 1), dtype=np.bool_, order='F')
    # T=1, so idx=1 maps to iteration 1, not the seed sentinel 0.
    ensemble_member_added_at_step = np.asfortranarray([[1]], dtype=np.int32)
    ensemble_U_first = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    ensemble_d_first = np.array([0], dtype=np.int32)
    ensemble_low_confidence_masks = np.zeros((1, 1), dtype=np.bool_, order='F')
    super_ensembles = np.zeros((1, 0), dtype=np.int32, order='F')

    filename = "test_stc_rmse_null_py.json"
    serialize_stc_results_as_json(
        filename, 0, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(1))

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert ('"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,'
            '"rmse":null,"drift":null}]') in content


def test_json_stop_reason_filter_excludes_pair():
    """Reuses _fixture()'s two intersecting ensembles (OC = 2/3), but overrides ensemble 2's
    stop reason to STOP_REASON_REJECTED_AFTER_STABLE. Without allowed_stop_reasons, both
    ensembles are eligible, the (1,2) pair appears in overlap_coefficient_matrix as usual, and
    params.excluded_stop_reasons is empty. With allowed_stop_reasons excluding
    STOP_REASON_REJECTED_AFTER_STABLE, ensemble 2 is eligible=False (this module honors the
    ensemble_eligible* arguments directly, not derived internally from allowed_stop_reasons --
    that argument is reported for transparency only, in params.excluded_stop_reasons), so the
    (1,2) pair must vanish from overlap_coefficient_matrix entirely."""
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history,
     ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
     ensemble_U_first, ensemble_d_first, super_ensembles) = _fixture()
    ensemble_stop_reason[1] = STOP_REASON_REJECTED_AFTER_STABLE

    # -- baseline: no filter, both ensembles eligible ------------------------------------
    filename = "test_stc_filter_baseline_py.json"
    serialize_stc_results_as_json(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles, **_common_kwargs(), **_all_eligible_kwargs(2))

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"excluded_stop_reasons":[]' in content
    assert ('"overlap_coefficient_matrix":[{"a":1,"b":2,"overlap_coefficient":'
            '6.6666666666666663E-001}]') in content
    assert '"reconciliation_eligible":true,"excluded_by":[]' in content

    # -- filtered: ensemble 2 ineligible by stop condition -------------------------------
    allowed = np.array([True, False, True, True], dtype=np.bool_)
    eligible = np.array([True, False], dtype=np.bool_)
    eligible_by_stop_condition = np.array([True, False], dtype=np.bool_)
    eligible_by_dimension = np.ones(2, dtype=np.bool_)
    eligible_by_var_explained = np.ones(2, dtype=np.bool_)

    filename = "test_stc_filter_excluded_py.json"
    serialize_stc_results_as_json(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles,
        allowed_stop_reasons=allowed, ensemble_eligible=eligible,
        ensemble_eligible_by_stop_condition=eligible_by_stop_condition,
        ensemble_eligible_by_dimension=eligible_by_dimension,
        ensemble_eligible_by_var_explained=eligible_by_var_explained,
        **_common_kwargs())

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert '"excluded_stop_reasons":["rejected_after_stable"]' in content
    assert '"overlap_coefficient_matrix":[]' in content
    # The excluded ensemble itself must still be fully reported -- the filter only ever
    # suppresses pairing, never the ensemble's own existence in the output.
    assert '"id":2,"seed_point_id":4,"stop_reason":"rejected_after_stable"' in content


def test_json_reconciliation_eligible_and_excluded_by_exact_strings():
    """Exact-string-match regression test: ensemble 1 is fully eligible (excluded_by an empty
    array); ensemble 2 is excluded by two criteria at once (stop condition and dimension, not
    variance explained). Also proves the JSON key is genuinely `excluded_by` (11 chars), never
    the 27-char `reconciliation_excluded_by` that silently truncated in `ensemble_keys`'s
    character(len=24) buffer during this module's own development."""
    (vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
     ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
     ensemble_G_history, ensemble_mu_history, ensemble_k_history,
     ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
     ensemble_U_first, ensemble_d_first, super_ensembles) = _fixture()

    eligible = np.array([True, False], dtype=np.bool_)
    eligible_by_stop_condition = np.array([True, False], dtype=np.bool_)
    eligible_by_dimension = np.array([True, False], dtype=np.bool_)
    eligible_by_var_explained = np.ones(2, dtype=np.bool_)

    filename = "test_stc_eligible_excluded_by_py.json"
    serialize_stc_results_as_json(
        filename, 1, vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason,
        ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history,
        ensemble_G_history, ensemble_mu_history, ensemble_k_history,
        ensemble_accepted_history, ensemble_member_added_at_step, ensemble_low_confidence_masks,
        ensemble_U_first, ensemble_d_first, super_ensembles,
        ensemble_eligible=eligible, ensemble_eligible_by_stop_condition=eligible_by_stop_condition,
        ensemble_eligible_by_dimension=eligible_by_dimension,
        ensemble_eligible_by_var_explained=eligible_by_var_explained,
        **_common_kwargs())

    with open(filename, "r") as f:
        content = f.read()
    os.remove(filename)

    assert ('"super_ensemble_id":1,"final_chordal_distance":null,'
            '"reconciliation_eligible":true,"excluded_by":[]}') in content
    assert ('"reconciliation_eligible":false,"excluded_by":["stop_condition","dimension"]') in content
    assert 'reconciliation_excluded_by' not in content


if __name__ == "__main__":
    run_all_tests(globals().values())
