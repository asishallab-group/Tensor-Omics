!> Unit test suite for tox_stc_json (serialize_stc_results_as_json,
!| write_stc_interactive_html_report). Fixtures are hand-crafted raw result arrays -- the same
!| shape ensemble_identification_merged/ensemble_reconciliation produce -- rather than run
!| through the full pipeline, since these two entry points only ever consume that shape and
!| never call into seeding/growth/observable/accept themselves.
module mod_test_shape_truthful_clustering_json
    use tox_stc_json, only: serialize_stc_results_as_json, write_stc_interactive_html_report
    use tox_shape_truthful_clustering_kernel, only: STOP_REASON_FIXED_POINT, &
        STOP_REASON_REJECTED_AFTER_STABLE
    use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_MERGE_OVERLAP_COEFFICIENT
    use tox_errors, only: is_ok, is_err
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_json() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(11))

        all_tests(1) = test_case("test_json_two_ensembles_with_overlap", test_json_two_ensembles_with_overlap)
        all_tests(2) = test_case("test_json_estimated_params_included", test_json_estimated_params_included)
        all_tests(3) = test_case("test_json_zero_ensembles", test_json_zero_ensembles)
        all_tests(4) = test_case("test_json_no_history_ensemble_omits_observable_keys", &
                                 test_json_no_history_ensemble_omits_observable_keys)
        all_tests(5) = test_case("test_html_report_wraps_json_in_template_and_d3", test_html_report_wraps_json_in_template_and_d3)
        all_tests(6) = test_case("test_json_invalid_n_dimensions", test_json_invalid_n_dimensions)
        all_tests(7) = test_case("test_json_rejected_trailing_column_uses_last_accepted", &
                                 test_json_rejected_trailing_column_uses_last_accepted)
        all_tests(8) = test_case("test_json_drift_and_final_chordal_distance", &
                                 test_json_drift_and_final_chordal_distance)
        all_tests(9) = test_case("test_json_observable_history_rmse_null_when_size_le_one", &
                                 test_json_observable_history_rmse_null_when_size_le_one)
        all_tests(10) = test_case("test_json_stop_reason_filter_excludes_pair", &
                                  test_json_stop_reason_filter_excludes_pair)
        all_tests(11) = test_case("test_json_reconciliation_eligible_and_excluded_by_exact_strings", &
                                  test_json_reconciliation_eligible_and_excluded_by_exact_strings)
    end function get_all_tests_shape_truthful_clustering_json

    !> D=2, N=4, 2 seeds/ensembles: {1,2,3} (seed=1, d=1) and {2,3,4} (seed=4, d=0), which
    !| overlap on {2,3} -- Overlap Coefficient 2/3 -- and are merged into one super-ensemble.
    !| Exercises: points' membership/seed_of/coords/residual_length, an ensemble with a tangent
    !| direction (u1/s1/line_start/line_end present, u2/s2 absent since d=1), an ensemble with
    !| d=0 (mu present, no tangent direction or line at all), super_ensembles, and the
    !| overlap_coefficient_matrix. Both ensembles are fully accepted at both retained history
    !| columns (ensemble_accepted_history all .true.) -- the rejected-trailing-column case gets
    !| its own dedicated fixture below. Every S_history entry that feeds an observable_history
    !| rmse is chosen so normal_error is a perfect square (4.0 or 9.0): sqrt(x + epsilon(1.0d0))
    !| for x this large rounds cleanly back to sqrt(x) (no rounding-tie risk the way x=0 would
    !| have, see the epsilon-guard note in accept_ensemble's own RMSE-drift criterion).
    subroutine build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                             vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                             ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                             ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                             ensemble_accepted_history, ensemble_member_added_at_step, &
                             ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles)
        integer(int32), intent(out) :: n_dimensions, n_vectors, n_selected_seed, o, max_group_size
        integer(int32), intent(out) :: n_super_ensembles
        real(real64), intent(out) :: vectors(2, 4)
        character(len=1), intent(out) :: dim_names(2)
        logical, intent(out) :: seed_selection_mask(4)
        logical, intent(out) :: ensemble_masks(4, 2)
        integer(int32), intent(out) :: ensemble_stop_reason(2)
        real(real64), intent(out) :: ensemble_growth_radii(2)
        real(real64), intent(out) :: ensemble_U_history(2, 2, 2, 2)
        real(real64), intent(out) :: ensemble_S_history(2, 2, 2)
        integer(int32), intent(out) :: ensemble_d_history(2, 2)
        real(real64), intent(out) :: ensemble_G_history(2, 2)
        real(real64), intent(out) :: ensemble_mu_history(2, 2, 2)
        integer(int32), intent(out) :: ensemble_k_history(2, 2)
        logical, intent(out) :: ensemble_accepted_history(2, 2)
        integer(int32), intent(out) :: ensemble_member_added_at_step(4, 2)
        logical, intent(out) :: ensemble_low_confidence_masks(4, 2)
        real(real64), intent(out) :: ensemble_U_first(2, 2, 2)
        integer(int32), intent(out) :: ensemble_d_first(2)
        integer(int32), intent(out) :: super_ensembles(2, 2)

        n_dimensions = 2
        n_vectors = 4
        n_selected_seed = 2
        o = 2
        max_group_size = 2
        n_super_ensembles = 1

        vectors(:, 1) = [0.0d0, 0.0d0]
        vectors(:, 2) = [1.0d0, 0.0d0]
        vectors(:, 3) = [2.0d0, 0.0d0]
        vectors(:, 4) = [3.0d0, 0.0d0]
        dim_names = ['x', 'y']

        seed_selection_mask = [.true., .false., .false., .true.]

        ensemble_masks(:, 1) = [.true., .true., .true., .false.]
        ensemble_masks(:, 2) = [.false., .true., .true., .true.]

        ensemble_stop_reason = [STOP_REASON_FIXED_POINT, STOP_REASON_FIXED_POINT]
        ensemble_growth_radii = [1.0d0, 1.0d0]

        ! k-1=1 at every retained column of both ensembles, to keep the eigenvalue arithmetic
        ! (S**2/(k-1)) below simple; k_history's exact value is otherwise unasserted.
        ensemble_k_history(:, 1) = [2, 2]
        ensemble_k_history(:, 2) = [2, 2]

        ensemble_d_history(:, 1) = [0, 1]
        ensemble_d_history(:, 2) = [0, 0]

        ensemble_G_history(:, 1) = [2.0d0, 1.5d0]
        ensemble_G_history(:, 2) = [2.0d0, 1.5d0]

        ensemble_mu_history(:, 1, 1) = [0.5d0, 0.0d0]
        ensemble_mu_history(:, 2, 1) = [1.0d0, 0.0d0]
        ensemble_mu_history(:, 1, 2) = [2.5d0, 0.0d0]
        ensemble_mu_history(:, 2, 2) = [3.0d0, 0.0d0]

        ensemble_S_history = 0.0d0
        ensemble_S_history(1, 1, 1) = 2.0d0 ! ensemble 1, iter 1 (d=0): normal_error = 4.0, rmse = 2.0
        ensemble_S_history(1, 2, 1) = 0.5d0 ! ensemble 1, iter 2 (d=1): s1, excluded from normal_error
        ensemble_S_history(2, 2, 1) = 3.0d0 ! ensemble 1, iter 2 (d=1): normal_error = 9.0, rmse = 3.0
        ensemble_S_history(1, 1, 2) = 2.0d0 ! ensemble 2, iter 1 (d=0): normal_error = 4.0, rmse = 2.0
        ensemble_S_history(1, 2, 2) = 3.0d0 ! ensemble 2, iter 2 (d=0): normal_error = 9.0, rmse = 3.0

        ensemble_U_history = 0.0d0
        ensemble_U_history(:, 1, 2, 1) = [1.0d0, 0.0d0] ! ensemble 1's u1 (tangent)
        ensemble_U_history(:, 2, 2, 1) = [0.0d0, 1.0d0] ! ensemble 1's normal direction (residual_length)
        ensemble_U_history(:, 1, 2, 2) = [1.0d0, 0.0d0] ! ensemble 2 has d=0, so this is only ever read by
        ensemble_U_history(:, 2, 2, 2) = [0.0d0, 1.0d0] ! residual_length (full identity basis, not "u1")

        ! Both ensembles accepted at every retained column -- no rejection in this fixture.
        ensemble_accepted_history = .true.

        ! Ensemble 1: seed = point 1 (vector index 1); point 2 joins at growth iteration 1,
        ! point 3 at iteration 2; point 4 never joins. T (max) = 2, matching idx=2 below.
        ensemble_member_added_at_step(:, 1) = [0, 1, 2, -1]
        ! Ensemble 2: seed = point 4; point 3 joins at iteration 1, point 2 at iteration 2;
        ! point 1 never joins. T (max) = 2.
        ensemble_member_added_at_step(:, 2) = [-1, 2, 1, 0]

        ! Bootstrap (iteration 1) basis, duplicating history column 1 exactly (as the real
        ! pipeline's own invariant guarantees, see "First growth step" in misc/mod_STC.md) --
        ! both d=0 here, so U_first is never actually read (see stc_chordal_distance's own
        ! d_a>0 applicability gate), only its shape must be consistent.
        ensemble_U_first(:, :, 1) = ensemble_U_history(:, :, 1, 1)
        ensemble_U_first(:, :, 2) = ensemble_U_history(:, :, 1, 2)
        ensemble_d_first(1) = ensemble_d_history(1, 1)
        ensemble_d_first(2) = ensemble_d_history(1, 2)

        ensemble_low_confidence_masks = .false.
        ensemble_low_confidence_masks(1, 1) = .true.

        super_ensembles(:, 1) = [1, 2]
        super_ensembles(:, 2) = [0, 0]
    end subroutine build_fixture

    subroutine test_json_two_ensembles_with_overlap()
        integer(int32) :: n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles
        real(real64) :: vectors(2, 4)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        logical :: ensemble_accepted_history(2, 2)
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        integer(int32) :: ensemble_member_added_at_step(4, 2), ensemble_d_first(2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2), ensemble_U_first(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        logical :: ensemble_eligible(2), ensemble_eligible_by_stop_condition(2), ensemble_eligible_by_dimension(2)
        logical :: ensemble_eligible_by_var_explained(2)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_accepted_history, ensemble_member_added_at_step, &
                           ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles)
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_two_ensembles.json'
        call serialize_stc_results_as_json(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "two ensembles: serialize_stc_results_as_json must not fail")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, '"dim_names":["x","y"]', "two ensembles: dim_names")
        call assert_string_contains(content, '"n_vectors":4', "two ensembles: n_vectors")
        call assert_string_contains(content, '"n_dimensions":2', "two ensembles: n_dimensions")
        call assert_string_contains(content, '"n_ensembles":2', "two ensembles: n_ensembles")
        call assert_string_contains(content, '"k_min":3', "two ensembles: k_min")
        call assert_string_contains(content, '"reconciliation_mode":"merge_overlap_coefficient"', &
                                    "two ensembles: reconciliation_mode name")

        ! points: residual_length is the point's distance off each ensemble's *normal* subspace
        ! at its final retained basis (identity here, see build_fixture) -- point 1 sits exactly
        ! at ensemble 1's mean (residual 0); point 2 sits at ensemble 1's mean too (residual 0)
        ! but 2.0 off ensemble 2's mean along its own (identity) normal directions.
        call assert_string_contains(content, &
            '"id":1,"coords":[0.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"n_ensembles":1,"n_low_confidence_ensembles":1,'//&
            '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000}],'//&
            '"low_confidence_ensembles":[1],"seed_of":[1]', "two ensembles: point 1")
        call assert_string_contains(content, &
            '"id":2,"coords":[1.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"n_ensembles":2,"n_low_confidence_ensembles":0,'//&
            '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000},'//&
            '{"id":2,"residual_length":2.0000000000000000E+000}],'//&
            '"low_confidence_ensembles":[],"seed_of":[]', "two ensembles: point 2")

        ! ensembles: #1 has d=1 (u1/s1/line_start/line_end present), #2 has d=0 (no tangent
        ! direction or line at all). Both ensembles' observable_history spans iterations 1-2
        ! (t_final=2, idx=2, so 1:idx=1:2 -- no gap). Neither ensemble's drift is ever
        ! computable here (ensemble 1's iter1-vs-iter2 d mismatches, 0 vs 1; ensemble 2's two
        ! iterations are both d=0, and stc_chordal_distance requires d>0 to be "applicable" at
        ! all) -- both final_chordal_distance values are null for the same reason (every
        ! reference comparison is d=0 vs d=1, or d=0 vs d=0).
        call assert_string_contains(content, &
            '{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'//&
            '"size":3,"t_final":2,"d":1,"G":1.5000000000000000E+000,'//&
            '"mu":[1.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"u1":[1.0000000000000000E+000,0.0000000000000000E+000],"s1":5.0000000000000000E-001,'//&
            '"line_start":[0.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"line_end":[2.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"observable_history":[{"iteration":1,"g":2.0000000000000000E+000,'//&
            '"rmse":2.0000000000000000E+000,"drift":null},{"iteration":2,"g":1.5000000000000000E+000,'//&
            '"rmse":3.0000000000000000E+000,"drift":null}],'//&
            '"super_ensemble_id":1,"final_chordal_distance":null,"reconciliation_eligible":true,"excluded_by":[]}', &
            "two ensembles: ensemble 1 (d=1)")
        call assert_string_contains(content, &
            '{"id":2,"seed_point_id":4,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'//&
            '"size":3,"t_final":2,"d":0,"G":1.5000000000000000E+000,'//&
            '"mu":[3.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"observable_history":[{"iteration":1,"g":2.0000000000000000E+000,'//&
            '"rmse":2.0000000000000000E+000,"drift":null},{"iteration":2,"g":1.5000000000000000E+000,'//&
            '"rmse":3.0000000000000000E+000,"drift":null}],'//&
            '"super_ensemble_id":1,"final_chordal_distance":null,"reconciliation_eligible":true,"excluded_by":[]}', &
            "two ensembles: ensemble 2 (d=0)")
        call assert_true(index(content, '"u2"') == 0, "two ensembles: no ensemble has d>=2, so no u2 key anywhere")

        call assert_string_contains(content, '"super_ensembles":[{"group_id":1,"ensemble_ids":[1,2]}]', &
                                    "two ensembles: super_ensembles")
        call assert_string_contains(content, &
            '"overlap_coefficient_matrix":[{"a":1,"b":2,"overlap_coefficient":6.6666666666666663E-001}]', &
            "two ensembles: overlap_coefficient_matrix")
    end subroutine test_json_two_ensembles_with_overlap

    subroutine test_json_estimated_params_included()
        integer(int32) :: n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles
        real(real64) :: vectors(2, 4)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        logical :: ensemble_accepted_history(2, 2)
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        integer(int32) :: ensemble_member_added_at_step(4, 2), ensemble_d_first(2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2), ensemble_U_first(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        logical :: ensemble_eligible(2), ensemble_eligible_by_stop_condition(2), ensemble_eligible_by_dimension(2)
        logical :: ensemble_eligible_by_var_explained(2)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_accepted_history, ensemble_member_added_at_step, &
                           ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles)
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_estimated_params.json'
        call serialize_stc_results_as_json(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            estimated_k_min=5_int32, estimated_k_density=6_int32, estimated_density_quantile=0.75d0, &
            estimated_chordal_dist_max_as_prcnt_of_range=0.2d0, estimated_G_max=3.0d0, estimated_d_max=2_int32, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "estimated params: serialize_stc_results_as_json must not fail")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, '"estimated_k_min":5', "estimated params: estimated_k_min present")
        call assert_string_contains(content, '"estimated_k_density":6', "estimated params: estimated_k_density present")
        call assert_string_contains(content, '"estimated_d_max":2', "estimated params: estimated_d_max present")
    end subroutine test_json_estimated_params_included

    !> `n_selected_seed = 0` is documented as a valid, well-defined "no ensembles" input, not
    !| an error -- every ensemble-shaped array collapses to a genuinely empty JSON array.
    subroutine test_json_zero_ensembles()
        real(real64) :: vectors(2, 2)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(2)
        logical :: ensemble_masks(2, 0), ensemble_low_confidence_masks(2, 0), ensemble_accepted_history(1, 0)
        integer(int32) :: ensemble_stop_reason(0), ensemble_d_history(1, 0), ensemble_k_history(1, 0)
        integer(int32) :: ensemble_member_added_at_step(2, 0), ensemble_d_first(0)
        real(real64) :: ensemble_growth_radii(0), ensemble_U_history(2, 2, 1, 0), ensemble_S_history(2, 1, 0)
        real(real64) :: ensemble_G_history(1, 0), ensemble_mu_history(2, 1, 0), ensemble_U_first(2, 2, 0)
        integer(int32) :: super_ensembles(2, 0)
        logical :: ensemble_eligible(0), ensemble_eligible_by_stop_condition(0), ensemble_eligible_by_dimension(0)
        logical :: ensemble_eligible_by_var_explained(0)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        vectors(:, 1) = [0.0d0, 0.0d0]
        vectors(:, 2) = [1.0d0, 0.0d0]
        dim_names = ['x', 'y']
        seed_selection_mask = .false.

        filename = 'test_stc_zero_ensembles.json'
        call serialize_stc_results_as_json(filename, &
            2_int32, 2_int32, 0_int32, 1_int32, 2_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "zero ensembles: must not be an error")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, '"ensembles":[]', "zero ensembles: ensembles array is empty")
        call assert_string_contains(content, '"super_ensembles":[]', "zero ensembles: super_ensembles array is empty")
        call assert_string_contains(content, '"overlap_coefficient_matrix":[]', &
                                    "zero ensembles: overlap_coefficient_matrix array is empty")
        call assert_string_contains(content, '"n_ensembles":0,"n_low_confidence_ensembles":0,'//&
                                    '"ensembles":[],"low_confidence_ensembles":[],"seed_of":[]', &
                                    "zero ensembles: point has no memberships at all")
    end subroutine test_json_zero_ensembles

    !> An ensemble whose k_history is entirely zero never produced an observable at all (only
    !| possible for STOP_REASON_MAX_SIZE firing at the bootstrap step itself). Its
    !| `t_final`/`d`/`G`/`mu`/tangent/`observable_history` keys must be omitted (not emitted as
    !| JSON null), not merely have null values -- `final_chordal_distance` is the one exception,
    !| always present (null here), since it is computed outside the "idx > 0" guard. The one
    !| point that is a member of this ensemble still gets its usual `ensembles` entry, with a
    !| residual_length of 0.0 (the documented fallback when no retained basis exists to project
    !| against at all).
    subroutine test_json_no_history_ensemble_omits_observable_keys()
        real(real64) :: vectors(2, 2)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(2)
        logical :: ensemble_masks(2, 1), ensemble_low_confidence_masks(2, 1), ensemble_accepted_history(1, 1)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(1, 1), ensemble_k_history(1, 1)
        integer(int32) :: ensemble_member_added_at_step(2, 1), ensemble_d_first(1)
        real(real64) :: ensemble_growth_radii(1), ensemble_U_history(2, 2, 1, 1), ensemble_S_history(2, 1, 1)
        real(real64) :: ensemble_G_history(1, 1), ensemble_mu_history(2, 1, 1), ensemble_U_first(2, 2, 1)
        integer(int32) :: super_ensembles(2, 0)
        logical :: ensemble_eligible(1), ensemble_eligible_by_stop_condition(1), ensemble_eligible_by_dimension(1)
        logical :: ensemble_eligible_by_var_explained(1)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        vectors(:, 1) = [0.0d0, 0.0d0]
        vectors(:, 2) = [1.0d0, 0.0d0]
        dim_names = ['x', 'y']
        seed_selection_mask = [.true., .false.]
        ensemble_masks(:, 1) = [.true., .false.]
        ensemble_stop_reason(1) = STOP_REASON_FIXED_POINT
        ensemble_growth_radii(1) = 1.0d0
        ensemble_k_history(:, 1) = [0]
        ensemble_d_history = 0
        ensemble_G_history = 0.0d0
        ensemble_mu_history = 0.0d0
        ensemble_S_history = 0.0d0
        ensemble_U_history = 0.0d0
        ensemble_accepted_history = .true.
        ensemble_member_added_at_step(:, 1) = [0, -1]
        ensemble_U_first = 0.0d0
        ensemble_d_first(1) = 0
        ensemble_low_confidence_masks = .false.
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_no_history.json'
        call serialize_stc_results_as_json(filename, &
            2_int32, 2_int32, 1_int32, 1_int32, 2_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "no history: must not be an error")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, &
            '{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'//&
            '"size":1,"super_ensemble_id":null,"final_chordal_distance":null,'//&
            '"reconciliation_eligible":true,"excluded_by":[]}', &
            "no history: d/G/mu/tangent/observable_history keys all omitted")
        call assert_true(index(content, '"t_final"') == 0, "no history: no 't_final' key anywhere")
        call assert_true(index(content, '"d":') == 0, "no history: no 'd' key anywhere")
        call assert_true(index(content, '"G":') == 0, "no history: no 'G' key anywhere")
        call assert_true(index(content, '"mu":') == 0, "no history: no 'mu' key anywhere")
        call assert_true(index(content, '"line_start"') == 0, "no history: no 'line_start' key anywhere")
        call assert_true(index(content, '"observable_history"') == 0, "no history: no 'observable_history' key anywhere")
        call assert_string_contains(content, &
            '"ensembles":[{"id":1,"residual_length":0.0000000000000000E+000}]', &
            "no history: member point still reports a fallback zero residual_length")
    end subroutine test_json_no_history_ensemble_omits_observable_keys

    subroutine test_html_report_wraps_json_in_template_and_d3()
        integer(int32) :: n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles
        real(real64) :: vectors(2, 4)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        logical :: ensemble_accepted_history(2, 2)
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        integer(int32) :: ensemble_member_added_at_step(4, 2), ensemble_d_first(2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2), ensemble_U_first(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        logical :: ensemble_eligible(2), ensemble_eligible_by_stop_condition(2), ensemble_eligible_by_dimension(2)
        logical :: ensemble_eligible_by_var_explained(2)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_accepted_history, ensemble_member_added_at_step, &
                           ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles)
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_report.html'
        call write_stc_interactive_html_report(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "html report: write_stc_interactive_html_report must not fail")

        ! Unformatted stream access, unlike formatted stream access, has no record structure
        ! at all -- reading this multi-line file with a formatted '(A)' read would stop at the
        ! first embedded newline instead of returning the whole content.
        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='unformatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, pos=1) content
        close (unit, status='delete')

        call assert_true(content(1:15) == '<!DOCTYPE html>', "html report: starts with the template head")
        call assert_string_contains(content, 'd3js.org', "html report: embeds the vendored D3 bundle")
        call assert_string_contains(content, 'const DATA = {"dim_names":["x","y"]', &
                                    "html report: embeds this run's JSON payload")
        call assert_string_contains(content, '</html>', "html report: ends with the template tail")
    end subroutine test_html_report_wraps_json_in_template_and_d3

    subroutine test_json_invalid_n_dimensions()
        real(real64) :: vectors(0, 2)
        character(len=1) :: dim_names(0)
        logical :: seed_selection_mask(2)
        logical :: ensemble_masks(2, 0), ensemble_low_confidence_masks(2, 0), ensemble_accepted_history(1, 0)
        integer(int32) :: ensemble_stop_reason(0), ensemble_d_history(1, 0), ensemble_k_history(1, 0)
        integer(int32) :: ensemble_member_added_at_step(2, 0), ensemble_d_first(0)
        real(real64) :: ensemble_growth_radii(0), ensemble_U_history(0, 0, 1, 0), ensemble_S_history(0, 1, 0)
        real(real64) :: ensemble_G_history(1, 0), ensemble_mu_history(0, 1, 0), ensemble_U_first(0, 0, 0)
        integer(int32) :: super_ensembles(2, 0)
        logical :: ensemble_eligible(0), ensemble_eligible_by_stop_condition(0), ensemble_eligible_by_dimension(0)
        logical :: ensemble_eligible_by_var_explained(0)
        integer(int32) :: ierr

        seed_selection_mask = .false.

        call serialize_stc_results_as_json('test_stc_invalid.json', &
            0_int32, 2_int32, 0_int32, 1_int32, 2_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_err(ierr), "invalid n_dimensions: n_dimensions=0 must be rejected")
    end subroutine test_json_invalid_n_dimensions

    !> Regression test for a real, pre-existing bug: `stc_push_ensemble_history` also pushes a
    !| *rejected* final candidate into the trailing history window right before growth halts on
    !| `STOP_REASON_REJECTED_IMMEDIATELY`/`STOP_REASON_REJECTED_AFTER_STABLE` -- so the last
    !| *populated* history column is not always the ensemble's real last *accepted* state. This
    !| fixture's column 1 is the true accepted state (d=1, G=1.0, mu=[0.5,0.0], u1=[1,0]); its
    !| column 2 is a rejected candidate with deliberately different, wrong-if-leaked values
    !| (d=0, G=99.0, mu=[9.0,9.0], accepted_history=.false.). Every reported "current state"
    !| field must reflect column 1, and observable_history must contain exactly that one
    !| (accepted) iteration, not two.
    subroutine test_json_rejected_trailing_column_uses_last_accepted()
        real(real64) :: vectors(2, 3)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(3)
        logical :: ensemble_masks(3, 1), ensemble_low_confidence_masks(3, 1), ensemble_accepted_history(2, 1)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(2, 1), ensemble_k_history(2, 1)
        integer(int32) :: ensemble_member_added_at_step(3, 1), ensemble_d_first(1)
        real(real64) :: ensemble_growth_radii(1), ensemble_U_history(2, 2, 2, 1), ensemble_S_history(2, 2, 1)
        real(real64) :: ensemble_G_history(2, 1), ensemble_mu_history(2, 2, 1), ensemble_U_first(2, 2, 1)
        integer(int32) :: super_ensembles(1, 0)
        logical :: ensemble_eligible(1), ensemble_eligible_by_stop_condition(1), ensemble_eligible_by_dimension(1)
        logical :: ensemble_eligible_by_var_explained(1)
        character(len=40) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        vectors(:, 1) = [0.0d0, 0.0d0]
        vectors(:, 2) = [1.0d0, 0.0d0]
        vectors(:, 3) = [9.0d0, 9.0d0] ! never a member -- the rejected candidate's would-be addition
        dim_names = ['x', 'y']
        seed_selection_mask = [.true., .false., .false.]
        ensemble_masks(:, 1) = [.true., .true., .false.]
        ensemble_stop_reason(1) = STOP_REASON_REJECTED_AFTER_STABLE
        ensemble_growth_radii(1) = 1.0d0

        ! Column 1: the real, accepted iteration 1 (bootstrap).
        ensemble_k_history(1, 1) = 2
        ensemble_d_history(1, 1) = 1
        ensemble_G_history(1, 1) = 1.0d0
        ensemble_mu_history(:, 1, 1) = [0.5d0, 0.0d0]
        ensemble_S_history(:, 1, 1) = [0.5d0, 2.0d0] ! eigen(2) = 4.0/(2-1) = 4.0 -> rmse = 2.0
        ensemble_U_history(:, 1, 1, 1) = [1.0d0, 0.0d0]
        ensemble_U_history(:, 2, 1, 1) = [0.0d0, 1.0d0]
        ensemble_accepted_history(1, 1) = .true.

        ! Column 2: the rejected candidate for iteration 2 -- deliberately different in every
        ! field, so any leak of this data into the JSON is immediately obvious.
        ensemble_k_history(2, 1) = 3
        ensemble_d_history(2, 1) = 0
        ensemble_G_history(2, 1) = 99.0d0
        ensemble_mu_history(:, 2, 1) = [9.0d0, 9.0d0]
        ensemble_S_history(:, 2, 1) = [7.0d0, 7.0d0]
        ensemble_U_history(:, :, 2, 1) = 0.0d0
        ensemble_accepted_history(2, 1) = .false.

        ! T = 1: only the bootstrap was ever actually accepted.
        ensemble_member_added_at_step(:, 1) = [0, 1, -1]

        ensemble_U_first(:, :, 1) = ensemble_U_history(:, :, 1, 1)
        ensemble_d_first(1) = ensemble_d_history(1, 1)
        ensemble_low_confidence_masks = .false.
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_rejected_trailing.json'
        call serialize_stc_results_as_json(filename, &
            2_int32, 3_int32, 1_int32, 2_int32, 1_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "rejected trailing column: must not fail")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, '"stop_reason":"rejected_after_stable"', &
                                    "rejected trailing column: stop_reason name")
        call assert_string_contains(content, '"size":2', "rejected trailing column: size from final_ensemble_mask")
        call assert_string_contains(content, '"t_final":1', &
                                    "rejected trailing column: T is 1, not 2 -- the rejected step never counted")
        call assert_string_contains(content, '"d":1', "rejected trailing column: d from the accepted column, not 0")
        call assert_string_contains(content, '"G":1.0000000000000000E+000', &
                                    "rejected trailing column: G from the accepted column, not 99.0")
        call assert_string_contains(content, '"mu":[5.0000000000000000E-001,0.0000000000000000E+000]', &
                                    "rejected trailing column: mu from the accepted column, not [9,9]")
        call assert_string_contains(content, '"u1":[1.0000000000000000E+000,0.0000000000000000E+000]', &
                                    "rejected trailing column: u1 present at all (the rejected column had d=0)")
        call assert_string_contains(content, &
            '"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,'//&
            '"rmse":2.0000000000000000E+000,"drift":null}]', &
            "rejected trailing column: observable_history has exactly the one accepted entry")
        call assert_true(index(content, '9.9000000000000000E+001') == 0, &
                         "rejected trailing column: the rejected G=99.0 never appears anywhere")
        call assert_true(index(content, '"mu":[9.0000000000000000E+000,9.0000000000000000E+000]') == 0, &
                         "rejected trailing column: the rejected mu=[9,9] never appears as an ensemble's mu " // &
                         "(point 3's own coords legitimately contain [9,9] too, so a bare substring check " // &
                         "would false-positive)")
    end subroutine test_json_rejected_trailing_column_uses_last_accepted

    !> A 3-retained-iteration, d=1-throughout ensemble whose tangent direction rotates 90
    !| degrees, then back: u1 = [1,0] (iter 1 = bootstrap = U_first), [0,1] (iter 2), [1,0]
    !| (iter 3). Every consecutive pair is orthogonal, so `stc_chordal_distance`'s
    !| sqrt(1 - cos(theta)**2) is exactly 1.0 for both drift entries -- no irrational
    !| intermediate values to hand-verify. `final_chordal_distance` compares iteration 3
    !| against {U_first, iter1, iter2}: 0.0 against the (identical) first two, 1.0 against
    !| iter2 -- max = 1.0.
    subroutine test_json_drift_and_final_chordal_distance()
        real(real64) :: vectors(2, 2)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(2)
        logical :: ensemble_masks(2, 1), ensemble_low_confidence_masks(2, 1), ensemble_accepted_history(3, 1)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(3, 1), ensemble_k_history(3, 1)
        integer(int32) :: ensemble_member_added_at_step(2, 1), ensemble_d_first(1)
        real(real64) :: ensemble_growth_radii(1), ensemble_U_history(2, 2, 3, 1), ensemble_S_history(2, 3, 1)
        real(real64) :: ensemble_G_history(3, 1), ensemble_mu_history(2, 3, 1), ensemble_U_first(2, 2, 1)
        integer(int32) :: super_ensembles(1, 0)
        logical :: ensemble_eligible(1), ensemble_eligible_by_stop_condition(1), ensemble_eligible_by_dimension(1)
        logical :: ensemble_eligible_by_var_explained(1)
        character(len=40) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        vectors(:, 1) = [0.0d0, 0.0d0]
        vectors(:, 2) = [5.0d0, 5.0d0]
        dim_names = ['x', 'y']
        seed_selection_mask = [.true., .false.]
        ensemble_masks(:, 1) = [.true., .true.]
        ensemble_stop_reason(1) = STOP_REASON_FIXED_POINT
        ensemble_growth_radii(1) = 1.0d0

        ensemble_k_history(:, 1) = [2, 2, 2]
        ensemble_d_history(:, 1) = [1, 1, 1]
        ensemble_G_history(:, 1) = [1.0d0, 2.0d0, 3.0d0]
        ensemble_mu_history(:, 1, 1) = [0.0d0, 0.0d0]
        ensemble_mu_history(:, 2, 1) = [0.0d0, 0.0d0]
        ensemble_mu_history(:, 3, 1) = [0.0d0, 0.0d0]

        ! s1 = 1.0 throughout (unasserted); the second singular value drives rmse via eigen(2)
        ! (d=1 excludes eigen(1) from normal_error): 2.0**2/1=4.0 -> rmse=2.0,
        ! 3.0**2/1=9.0 -> rmse=3.0, 4.0**2/1=16.0 -> rmse=4.0.
        ensemble_S_history(:, 1, 1) = [1.0d0, 2.0d0]
        ensemble_S_history(:, 2, 1) = [1.0d0, 3.0d0]
        ensemble_S_history(:, 3, 1) = [1.0d0, 4.0d0]

        ensemble_U_history(:, 1, 1, 1) = [1.0d0, 0.0d0]
        ensemble_U_history(:, 1, 2, 1) = [0.0d0, 1.0d0]
        ensemble_U_history(:, 1, 3, 1) = [1.0d0, 0.0d0]

        ensemble_accepted_history(:, 1) = [.true., .true., .true.]
        ensemble_member_added_at_step(:, 1) = [0, 3] ! T=3; which real member joined "when" is irrelevant here

        ensemble_U_first(:, :, 1) = 0.0d0
        ensemble_U_first(:, 1, 1) = [1.0d0, 0.0d0] ! matches iteration 1 exactly, as the real pipeline guarantees
        ensemble_d_first(1) = 1

        ensemble_low_confidence_masks = .false.
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_drift.json'
        call serialize_stc_results_as_json(filename, &
            2_int32, 2_int32, 1_int32, 3_int32, 1_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "drift/final chordal distance: must not fail")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, '"t_final":3', "drift: t_final")
        call assert_string_contains(content, &
            '"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,'//&
            '"rmse":2.0000000000000000E+000,"drift":null},'//&
            '{"iteration":2,"g":2.0000000000000000E+000,"rmse":3.0000000000000000E+000,'//&
            '"drift":1.0000000000000000E+000},'//&
            '{"iteration":3,"g":3.0000000000000000E+000,"rmse":4.0000000000000000E+000,'//&
            '"drift":1.0000000000000000E+000}]', &
            "drift: iteration 1 has no drift (nothing precedes it); iterations 2/3 are both exactly 1.0 " // &
            "(consecutive u1 vectors are orthogonal, sqrt(1-0**2))")
        call assert_string_contains(content, '"final_chordal_distance":1.0000000000000000E+000', &
                                    "drift: final_chordal_distance is the max over {U_first, iter1, iter2} " // &
                                    "vs iter3 -- 0.0, 0.0, 1.0 -- i.e. 1.0")
    end subroutine test_json_drift_and_final_chordal_distance

    !> `observable_history`'s `rmse` key is present but null when that column's ensemble size
    !| was <=1 (k_history entry of 1 -- normal_error/(k-1) would divide by zero). This can only
    !| ever occur at a genuinely degenerate retained column; this module does not itself
    !| validate that k_history is consistent with any real growth trajectory, so the fixture
    !| constructs it directly.
    subroutine test_json_observable_history_rmse_null_when_size_le_one()
        real(real64) :: vectors(2, 1)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(1)
        logical :: ensemble_masks(1, 1), ensemble_low_confidence_masks(1, 1), ensemble_accepted_history(1, 1)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(1, 1), ensemble_k_history(1, 1)
        integer(int32) :: ensemble_member_added_at_step(1, 1), ensemble_d_first(1)
        real(real64) :: ensemble_growth_radii(1), ensemble_U_history(2, 2, 1, 1), ensemble_S_history(2, 1, 1)
        real(real64) :: ensemble_G_history(1, 1), ensemble_mu_history(2, 1, 1), ensemble_U_first(2, 2, 1)
        integer(int32) :: super_ensembles(1, 0)
        logical :: ensemble_eligible(1), ensemble_eligible_by_stop_condition(1), ensemble_eligible_by_dimension(1)
        logical :: ensemble_eligible_by_var_explained(1)
        character(len=40) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        vectors(:, 1) = [0.0d0, 0.0d0]
        dim_names = ['x', 'y']
        seed_selection_mask = [.true.]
        ensemble_masks(:, 1) = [.true.]
        ensemble_stop_reason(1) = STOP_REASON_FIXED_POINT
        ensemble_growth_radii(1) = 1.0d0
        ensemble_k_history(1, 1) = 1 ! degenerate: k-1 = 0
        ensemble_d_history(1, 1) = 0
        ensemble_G_history(1, 1) = 1.0d0
        ensemble_mu_history(:, 1, 1) = [0.0d0, 0.0d0]
        ensemble_S_history(:, 1, 1) = [1.0d0, 1.0d0]
        ensemble_U_history(:, :, 1, 1) = 0.0d0
        ensemble_accepted_history(1, 1) = .true.
        ensemble_member_added_at_step(1, 1) = 1 ! T=1, so idx=1 maps to iteration 1, not the seed sentinel 0
        ensemble_U_first(:, :, 1) = 0.0d0
        ensemble_d_first(1) = 0
        ensemble_low_confidence_masks = .false.
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_rmse_null.json'
        call serialize_stc_results_as_json(filename, &
            2_int32, 1_int32, 1_int32, 1_int32, 1_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "rmse null: must not fail (no division-by-zero crash)")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, &
            '"observable_history":[{"iteration":1,"g":1.0000000000000000E+000,"rmse":null,"drift":null}]', &
            "rmse null: rmse key present but null when k_history entry is 1")
    end subroutine test_json_observable_history_rmse_null_when_size_le_one

    !> Reuses `build_fixture`'s two intersecting ensembles (OC = 2/3, see
    !| test_json_two_ensembles_with_overlap), but overrides ensemble 2's stop reason to
    !| STOP_REASON_REJECTED_AFTER_STABLE. Without `allowed_stop_reasons`, the (1,2) pair
    !| appears in `overlap_coefficient_matrix` as usual and `params.excluded_stop_reasons` is
    !| empty. With `allowed_stop_reasons` excluding STOP_REASON_REJECTED_AFTER_STABLE, ensemble
    !| 2 becomes ineligible, so the (1,2) pair must vanish from `overlap_coefficient_matrix`
    !| entirely (not merely zeroed) and `params.excluded_stop_reasons` must report exactly
    !| `["rejected_after_stable"]`. Mirrors this module's own filter -- see
    !| `tox_shape_truthful_clustering_reconciliation_kernel`'s `allowed_stop_reasons`, which
    !| this JSON-layer filter is deliberately kept consistent with (misc/mod_STC.md, "Filtering
    !| ensembles out of reconciliation by Stop Condition").
    subroutine test_json_stop_reason_filter_excludes_pair()
        integer(int32) :: n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles
        real(real64) :: vectors(2, 4)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        logical :: ensemble_accepted_history(2, 2)
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        integer(int32) :: ensemble_member_added_at_step(4, 2), ensemble_d_first(2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2), ensemble_U_first(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        logical :: allowed(4)
        logical :: ensemble_eligible(2), ensemble_eligible_by_stop_condition(2), ensemble_eligible_by_dimension(2)
        logical :: ensemble_eligible_by_var_explained(2)
        character(len=40) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_accepted_history, ensemble_member_added_at_step, &
                           ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles)
        ensemble_stop_reason(2) = STOP_REASON_REJECTED_AFTER_STABLE

        ! -- baseline: no filter, both ensembles eligible ------------------------------------
        ensemble_eligible = .true.
        ensemble_eligible_by_stop_condition = .true.
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_filter_baseline.json'
        call serialize_stc_results_as_json(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "stop reason filter: baseline must not fail")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, '"excluded_stop_reasons":[]', &
                                    "stop reason filter: baseline reports no exclusions")
        call assert_string_contains(content, &
            '"overlap_coefficient_matrix":[{"a":1,"b":2,"overlap_coefficient":6.6666666666666663E-001}]', &
            "stop reason filter: baseline still has the (1,2) pair")
        call assert_string_contains(content, '"reconciliation_eligible":true,"excluded_by":[]', &
                                    "stop reason filter: baseline -- reconciliation_eligible/excluded_by present " // &
                                    "at least once, both ensembles fully eligible")
        deallocate (content)

        ! -- filtered: ensemble 2 ineligible by stop condition (allowed_stop_reasons here is
        ! reported for transparency only, in params.excluded_stop_reasons -- the actual
        ! eligibility this module honors comes from the ensemble_eligible* arguments below,
        ! exactly as `ensemble_reconciliation` itself would have computed them) --------------
        allowed = .true.
        allowed(STOP_REASON_REJECTED_AFTER_STABLE) = .false.
        ensemble_eligible = [.true., .false.]
        ensemble_eligible_by_stop_condition = [.true., .false.]
        ensemble_eligible_by_dimension = .true.
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_filter_excluded.json'
        call serialize_stc_results_as_json(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            allowed_stop_reasons=allowed, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "stop reason filter: filtered call must not fail")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, '"excluded_stop_reasons":["rejected_after_stable"]', &
                                    "stop reason filter: filtered call reports the exclusion")
        call assert_string_contains(content, '"overlap_coefficient_matrix":[]', &
                                    "stop reason filter: (1,2) pair is gone once ensemble 2 is ineligible")
        ! Ensemble 2 itself must still be fully reported -- the filter only ever suppresses
        ! pairing, never the ensemble's own existence in the output (see the doc note above).
        call assert_string_contains(content, '"id":2,"seed_point_id":4,"stop_reason":"rejected_after_stable"', &
                                    "stop reason filter: excluded ensemble 2 still fully reported")
    end subroutine test_json_stop_reason_filter_excludes_pair

    !> Exact-string-match regression test for `reconciliation_eligible`/`excluded_by`: ensemble 1
    !| is fully eligible (`excluded_by` an empty array); ensemble 2 is excluded by *two*
    !| criteria at once (stop condition and dimension, not variance explained) -- proving
    !| `excluded_by` lists every failing criterion, not just the first, and that the JSON key is
    !| genuinely `excluded_by` (11 chars), never the 27-char `reconciliation_excluded_by` that
    !| silently truncated in `ensemble_keys`'s `character(len=24)` buffer during this module's
    !| own development (see `tox_stc_json.F90`'s comment at the `reconciliation_eligible` key
    !| assignment).
    subroutine test_json_reconciliation_eligible_and_excluded_by_exact_strings()
        integer(int32) :: n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles
        real(real64) :: vectors(2, 4)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        logical :: ensemble_accepted_history(2, 2)
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        integer(int32) :: ensemble_member_added_at_step(4, 2), ensemble_d_first(2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2), ensemble_U_first(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        logical :: ensemble_eligible(2), ensemble_eligible_by_stop_condition(2), ensemble_eligible_by_dimension(2)
        logical :: ensemble_eligible_by_var_explained(2)
        character(len=48) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_accepted_history, ensemble_member_added_at_step, &
                           ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles)

        ensemble_eligible = [.true., .false.]
        ensemble_eligible_by_stop_condition = [.true., .false.]
        ensemble_eligible_by_dimension = [.true., .false.]
        ensemble_eligible_by_var_explained = .true.

        filename = 'test_stc_eligible_excluded_by.json'
        call serialize_stc_results_as_json(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_accepted_history, ensemble_member_added_at_step, &
            ensemble_low_confidence_masks, ensemble_U_first, ensemble_d_first, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ensemble_eligible=ensemble_eligible, ensemble_eligible_by_stop_condition=ensemble_eligible_by_stop_condition, &
            ensemble_eligible_by_dimension=ensemble_eligible_by_dimension, &
            ensemble_eligible_by_var_explained=ensemble_eligible_by_var_explained, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "eligible/excluded_by: must not fail")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, &
            '"super_ensemble_id":1,"final_chordal_distance":null,"reconciliation_eligible":true,"excluded_by":[]}', &
            "eligible/excluded_by: ensemble 1 fully eligible, excluded_by is an empty array")
        call assert_string_contains(content, &
            '"reconciliation_eligible":false,"excluded_by":["stop_condition","dimension"]', &
            "eligible/excluded_by: ensemble 2 lists both failing criteria, in stop/dimension/var order")
        call assert_true(index(content, 'reconciliation_excluded_by') == 0, &
                         "eligible/excluded_by: the truncated 27-char key name never appears")
    end subroutine test_json_reconciliation_eligible_and_excluded_by_exact_strings

end module mod_test_shape_truthful_clustering_json
