!> Unit test suite for tox_stc_json (serialize_stc_results_as_json,
!| write_stc_interactive_html_report). Fixtures are hand-crafted raw result arrays -- the same
!| shape ensemble_identification_merged/ensemble_reconciliation produce -- rather than run
!| through the full pipeline, since these two entry points only ever consume that shape and
!| never call into seeding/growth/observable/accept themselves.
module mod_test_shape_truthful_clustering_json
    use tox_stc_json, only: serialize_stc_results_as_json, write_stc_interactive_html_report
    use tox_shape_truthful_clustering_kernel, only: STOP_REASON_FIXED_POINT
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
        allocate (all_tests(6))

        all_tests(1) = test_case("test_json_two_ensembles_with_overlap", test_json_two_ensembles_with_overlap)
        all_tests(2) = test_case("test_json_estimated_params_included", test_json_estimated_params_included)
        all_tests(3) = test_case("test_json_zero_ensembles", test_json_zero_ensembles)
        all_tests(4) = test_case("test_json_no_history_ensemble_omits_observable_keys", &
                                 test_json_no_history_ensemble_omits_observable_keys)
        all_tests(5) = test_case("test_html_report_wraps_json_in_template_and_d3", test_html_report_wraps_json_in_template_and_d3)
        all_tests(6) = test_case("test_json_invalid_n_dimensions", test_json_invalid_n_dimensions)
    end function get_all_tests_shape_truthful_clustering_json

    !> D=2, N=4, 2 seeds/ensembles: {1,2,3} (seed=1, d=1) and {2,3,4} (seed=4, d=0), which
    !| overlap on {2,3} -- Overlap Coefficient 2/3 -- and are merged into one super-ensemble.
    !| Exercises: points' membership/seed_of/coords, an ensemble with a tangent direction
    !| (u1/s1 present, u2/s2 absent since d=1), an ensemble with d=0 (mu present, no tangent
    !| direction at all), super_ensembles, and the overlap_coefficient_matrix.
    subroutine build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                             vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                             ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                             ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                             ensemble_low_confidence_masks, super_ensembles)
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
        logical, intent(out) :: ensemble_low_confidence_masks(4, 2)
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

        ensemble_k_history(:, 1) = [2, 3]
        ensemble_k_history(:, 2) = [2, 3]

        ensemble_d_history(:, 1) = [0, 1]
        ensemble_d_history(:, 2) = [0, 0]

        ensemble_G_history(:, 1) = [2.0d0, 1.5d0]
        ensemble_G_history(:, 2) = [2.0d0, 1.5d0]

        ensemble_mu_history(:, 1, 1) = [0.5d0, 0.0d0]
        ensemble_mu_history(:, 2, 1) = [1.0d0, 0.0d0]
        ensemble_mu_history(:, 1, 2) = [2.5d0, 0.0d0]
        ensemble_mu_history(:, 2, 2) = [3.0d0, 0.0d0]

        ensemble_S_history = 0.0d0
        ensemble_S_history(1, 2, 1) = 0.5d0

        ensemble_U_history = 0.0d0
        ensemble_U_history(:, 1, 2, 1) = [1.0d0, 0.0d0]

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
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_low_confidence_masks, super_ensembles)

        filename = 'test_stc_two_ensembles.json'
        call serialize_stc_results_as_json(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_low_confidence_masks, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
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

        ! points
        call assert_string_contains(content, '"id":1,"coords":[0.0000000000000000E+000,0.0000000000000000E+000],'//&
                                    '"n_ensembles":1,"n_low_confidence_ensembles":1,"ensembles":[1],'//&
                                    '"low_confidence_ensembles":[1],"seed_of":[1]', "two ensembles: point 1")
        call assert_string_contains(content, '"id":2,"coords":[1.0000000000000000E+000,0.0000000000000000E+000],'//&
                                    '"n_ensembles":2,"n_low_confidence_ensembles":0,"ensembles":[1,2],'//&
                                    '"low_confidence_ensembles":[],"seed_of":[]', "two ensembles: point 2")

        ! ensembles: #1 has d=1 (u1/s1 present), #2 has d=0 (no tangent direction at all)
        call assert_string_contains(content, &
            '{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'//&
            '"size":3,"d":1,"G":1.5000000000000000E+000,"mu":[1.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"u1":[1.0000000000000000E+000,0.0000000000000000E+000],"s1":5.0000000000000000E-001,'//&
            '"super_ensemble_id":1}', "two ensembles: ensemble 1 (d=1)")
        call assert_string_contains(content, &
            '{"id":2,"seed_point_id":4,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'//&
            '"size":3,"d":0,"G":1.5000000000000000E+000,"mu":[3.0000000000000000E+000,0.0000000000000000E+000],'//&
            '"super_ensemble_id":1}', "two ensembles: ensemble 2 (d=0, no tangent direction)")
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
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_low_confidence_masks, super_ensembles)

        filename = 'test_stc_estimated_params.json'
        call serialize_stc_results_as_json(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_low_confidence_masks, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
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
        logical :: ensemble_masks(2, 0), ensemble_low_confidence_masks(2, 0)
        integer(int32) :: ensemble_stop_reason(0), ensemble_d_history(1, 0), ensemble_k_history(1, 0)
        real(real64) :: ensemble_growth_radii(0), ensemble_U_history(2, 2, 1, 0), ensemble_S_history(2, 1, 0)
        real(real64) :: ensemble_G_history(1, 0), ensemble_mu_history(2, 1, 0)
        integer(int32) :: super_ensembles(2, 0)
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
            ensemble_low_confidence_masks, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
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
    !| possible for STOP_REASON_MAX_SIZE firing at the bootstrap step itself). Its `d`/`G`/`mu`
    !| keys must be omitted (not emitted as JSON null), not merely have null values.
    subroutine test_json_no_history_ensemble_omits_observable_keys()
        real(real64) :: vectors(2, 2)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(2)
        logical :: ensemble_masks(2, 1), ensemble_low_confidence_masks(2, 1)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(1, 1), ensemble_k_history(1, 1)
        real(real64) :: ensemble_growth_radii(1), ensemble_U_history(2, 2, 1, 1), ensemble_S_history(2, 1, 1)
        real(real64) :: ensemble_G_history(1, 1), ensemble_mu_history(2, 1, 1)
        integer(int32) :: super_ensembles(2, 0)
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
        ensemble_low_confidence_masks = .false.

        filename = 'test_stc_no_history.json'
        call serialize_stc_results_as_json(filename, &
            2_int32, 2_int32, 1_int32, 1_int32, 2_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_low_confidence_masks, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ierr=ierr)
        call assert_true(is_ok(ierr), "no history: must not be an error")

        inquire (file=trim(filename), size=filesize)
        open (newunit=unit, file=trim(filename), access='stream', form='formatted', status='old', action='read')
        allocate (character(len=filesize) :: content)
        read (unit, '(A)') content
        close (unit, status='delete')

        call assert_string_contains(content, &
            '{"id":1,"seed_point_id":1,"stop_reason":"fixed_point","growth_radius":1.0000000000000000E+000,'//&
            '"size":1,"super_ensemble_id":null}', "no history: d/G/mu/tangent keys all omitted")
        call assert_true(index(content, '"d":') == 0, "no history: no 'd' key anywhere")
        call assert_true(index(content, '"G":') == 0, "no history: no 'G' key anywhere")
        call assert_true(index(content, '"mu":') == 0, "no history: no 'mu' key anywhere")
    end subroutine test_json_no_history_ensemble_omits_observable_keys

    subroutine test_html_report_wraps_json_in_template_and_d3()
        integer(int32) :: n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles
        real(real64) :: vectors(2, 4)
        character(len=1) :: dim_names(2)
        logical :: seed_selection_mask(4), ensemble_masks(4, 2), ensemble_low_confidence_masks(4, 2)
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(2, 2), ensemble_k_history(2, 2)
        real(real64) :: ensemble_growth_radii(2), ensemble_U_history(2, 2, 2, 2), ensemble_S_history(2, 2, 2)
        real(real64) :: ensemble_G_history(2, 2), ensemble_mu_history(2, 2, 2)
        integer(int32) :: super_ensembles(2, 2)
        character(len=32) :: filename
        character(len=:), allocatable :: content
        integer(int32) :: ierr, unit, filesize

        call build_fixture(n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
                           vectors, dim_names, seed_selection_mask, ensemble_masks, ensemble_stop_reason, &
                           ensemble_growth_radii, ensemble_U_history, ensemble_S_history, ensemble_d_history, &
                           ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
                           ensemble_low_confidence_masks, super_ensembles)

        filename = 'test_stc_report.html'
        call write_stc_interactive_html_report(filename, &
            n_dimensions, n_vectors, n_selected_seed, o, max_group_size, n_super_ensembles, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_low_confidence_masks, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
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
        logical :: ensemble_masks(2, 0), ensemble_low_confidence_masks(2, 0)
        integer(int32) :: ensemble_stop_reason(0), ensemble_d_history(1, 0), ensemble_k_history(1, 0)
        real(real64) :: ensemble_growth_radii(0), ensemble_U_history(0, 0, 1, 0), ensemble_S_history(0, 1, 0)
        real(real64) :: ensemble_G_history(1, 0), ensemble_mu_history(0, 1, 0)
        integer(int32) :: super_ensembles(2, 0)
        integer(int32) :: ierr

        seed_selection_mask = .false.

        call serialize_stc_results_as_json('test_stc_invalid.json', &
            0_int32, 2_int32, 0_int32, 1_int32, 2_int32, 0_int32, &
            vectors, dim_names, seed_selection_mask, &
            ensemble_masks, ensemble_stop_reason, ensemble_growth_radii, &
            ensemble_U_history, ensemble_S_history, ensemble_d_history, &
            ensemble_G_history, ensemble_mu_history, ensemble_k_history, &
            ensemble_low_confidence_masks, super_ensembles, &
            k_min=3_int32, k_density=4_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=1_int32, &
            G_max=2.0d0, RMSE_change_max=0.5d0, f_max=0.8d0, a=3_int32, &
            exclusion_radius_percentile=50.0d0, bandwidth_percentile=68.0d0, &
            reconciliation_mode=MODE_MERGE_OVERLAP_COEFFICIENT, min_overlap_coefficient=0.5d0, &
            ierr=ierr)
        call assert_true(is_err(ierr), "invalid n_dimensions: n_dimensions=0 must be rejected")
    end subroutine test_json_invalid_n_dimensions

end module mod_test_shape_truthful_clustering_json
