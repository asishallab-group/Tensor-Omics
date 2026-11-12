!> Unit test suite for shatter clustering data calculations.
module mod_test_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use asserts
    use test_suite, only: test_case
    use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, get_err_arg_pos
    use f42_kd_tree, only: build_kd_index_expert
    use f42_kd_tree_impl, only: KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH
    use tox_stc, only: calculate_density_radius, calculate_density_radius_expert, &
                       calculate_labels_as_density, calculate_labels_as_density_expert, &
                       identify_ensemble_seeds, identify_ensemble_seeds_expert, &
                       grow_ensemble, grow_ensemble_expert, compute_ensemble_observable, accept_ensemble, &
                       obtain_ensembles_global_radius, obtain_ensembles_seed_radius, &
                       obtain_ensembles_seed_radius_expert, merge_ensembles, merge_ensembles_expert
    use tox_stc_impl, only: grow_ensemble_impl, obtain_ensembles_impl, MODE_GLOBAL_RADIUS, grow_single_seed_helper, compute_ensemble_observable_helper, accept_ensemble_helper, &
                            compute_ambient_density_stats_helper => compute_ambient_density_stats_impl, &
                            STC_STOP_FIXED_POINT, STC_STOP_REJECT_AFTER_ACCEPT, &
                            STC_STOP_NEVER_ACCEPTED, STC_STOP_NO_CANDIDATES

    implicit none
    public

contains

    function get_all_tests_shatter_cluster_data() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(70))
        all_tests(1) = test_case("test_density_radius_basic", test_density_radius_basic)
        all_tests(2) = test_case("test_density_labels_basic", test_density_labels_basic)
        all_tests(3) = test_case("test_density_radius_invalid_quantile", test_density_radius_invalid_quantile)
        all_tests(4) = test_case("test_density_labels_invalid_r", test_density_labels_invalid_r)
        all_tests(5) = test_case("test_density_labels_invalid_kd_indices", test_density_labels_invalid_kd_indices)
        all_tests(6) = test_case("test_density_radius_default_quantile", test_density_radius_default_quantile)
        all_tests(7) = test_case("test_density_labels_small_radius", test_density_labels_small_radius)
        all_tests(8) = test_case("test_density_labels_single_vector", test_density_labels_single_vector)
        all_tests(9) = test_case("test_density_labels_identical_vectors", test_density_labels_identical_vectors)
        all_tests(10) = test_case("test_density_radius_plain", test_density_radius_plain)
        all_tests(11) = test_case("test_density_labels_plain", test_density_labels_plain)
        all_tests(12) = test_case("test_identify_ensemble_seeds_basic", test_identify_ensemble_seeds_basic)
        all_tests(13) = test_case("test_identify_ensemble_seeds_k_seeding_effect", test_identify_ensemble_seeds_k_seeding_effect)
        all_tests(14) = test_case("test_identify_ensemble_seeds_invalid_k_seeding", test_identify_ensemble_seeds_invalid_k_seeding)
        all_tests(15) = test_case("test_identify_ensemble_seeds_plain", test_identify_ensemble_seeds_plain)
        all_tests(16) = test_case("test_identify_ensemble_seeds_single_vector_invalid", test_identify_ensemble_seeds_single_vector_invalid)
        all_tests(17) = test_case("test_identify_ensemble_seeds_invalid_inputs", test_identify_ensemble_seeds_invalid_inputs)
        all_tests(18) = test_case("test_identify_ensemble_seeds_identical_vectors", test_identify_ensemble_seeds_identical_vectors)
        all_tests(19) = test_case("test_grow_ensemble_basic", test_grow_ensemble_basic)
        all_tests(20) = test_case("test_grow_ensemble_no_growth", test_grow_ensemble_no_growth)
        all_tests(21) = test_case("test_grow_ensemble_invalid_inputs", test_grow_ensemble_invalid_inputs)
        all_tests(22) = test_case("test_grow_ensemble_plain", test_grow_ensemble_plain)
        all_tests(23) = test_case("test_grow_ensemble_multistep", test_grow_ensemble_multistep)
        all_tests(24) = test_case("test_grow_ensemble_alpha_sensitivity", test_grow_ensemble_alpha_sensitivity)
        all_tests(25) = test_case("test_grow_ensemble_zero_mad", test_grow_ensemble_zero_mad)
        all_tests(26) = test_case("test_grow_ensemble_mask_states", test_grow_ensemble_mask_states)
        all_tests(27) = test_case("test_compute_ensemble_observable_basic", test_compute_ensemble_observable_basic)
        all_tests(28) = test_case("test_compute_ensemble_observable_sliding_window", test_compute_ensemble_observable_sliding_window)
        all_tests(29) = test_case("test_compute_ensemble_observable_history_width", test_compute_ensemble_observable_history_width)
        all_tests(30) = test_case("test_compute_ensemble_observable_empty_and_zero", test_compute_ensemble_observable_empty_and_zero)
        all_tests(31) = test_case("test_compute_ensemble_observable_plain", test_compute_ensemble_observable_plain)
        all_tests(32) = test_case("test_compute_ensemble_observable_invalid_inputs", test_compute_ensemble_observable_invalid_inputs)
        all_tests(33) = test_case("test_accept_ensemble_default", test_accept_ensemble_default)
        all_tests(34) = test_case("test_accept_ensemble_alpha_sensitivity", test_accept_ensemble_alpha_sensitivity)
        all_tests(35) = test_case("test_accept_ensemble_boundary_and_window", test_accept_ensemble_boundary_and_window)
        all_tests(36) = test_case("test_accept_ensemble_invalid_inputs", test_accept_ensemble_invalid_inputs)
        all_tests(37) = test_case("test_accept_ensemble_helper_basic", test_accept_ensemble_helper_basic)
        all_tests(38) = test_case("test_obtain_ensembles_basic", test_obtain_ensembles_basic)
        all_tests(39) = test_case("test_obtain_ensembles_unmerged", test_obtain_ensembles_unmerged)
        all_tests(40) = test_case("test_obtain_ensembles_zero_seeds", test_obtain_ensembles_zero_seeds)
        all_tests(41) = test_case("test_obtain_ensembles_invalid_inputs", test_obtain_ensembles_invalid_inputs)
        all_tests(42) = test_case("test_obtain_ensembles_500_vectors", test_obtain_ensembles_500_vectors)
        all_tests(43) = test_case("test_merge_ensembles_basic", test_merge_ensembles_basic)
        all_tests(44) = test_case("test_merge_ensembles_zero_seeds", test_merge_ensembles_zero_seeds)
        all_tests(45) = test_case("test_merge_ensembles_no_overlap", test_merge_ensembles_no_overlap)
        all_tests(46) = test_case("test_merge_ensembles_invalid_inputs", test_merge_ensembles_invalid_inputs)
        all_tests(47) = test_case("test_merge_ensembles_filters_final_singletons", &
                                  test_merge_ensembles_filters_final_singletons)
        all_tests(48) = test_case("test_obtain_ensembles_tile_count_invariance", &
                                  test_obtain_ensembles_tile_count_invariance)
        all_tests(49) = test_case("test_obtain_ensembles_invalid_n_tiles", &
                                  test_obtain_ensembles_invalid_n_tiles)
        all_tests(50) = test_case("test_compute_ambient_density_stats", &
                                  test_compute_ambient_density_stats)
        all_tests(51) = test_case("test_obtain_ensembles_layer_n_tiles_agreement", &
                                  test_obtain_ensembles_layer_n_tiles_agreement)
        all_tests(52) = test_case("test_grow_ensemble_rejects_negative_density_labels", &
                                  test_grow_ensemble_rejects_negative_density_labels)
        all_tests(53) = test_case("test_observable_history_matrix", &
                                  test_observable_history_matrix)
        all_tests(54) = test_case("test_obtain_ensembles_rejects_invalid_history", &
                                  test_obtain_ensembles_rejects_invalid_history)
        all_tests(55) = test_case("test_obtain_ensembles_stop_reasons", &
                                  test_obtain_ensembles_stop_reasons)
        all_tests(56) = test_case("test_obtain_ensembles_observable_window_match", &
                                  test_obtain_ensembles_observable_window_match)
        all_tests(57) = test_case("test_obtain_ensembles_reports_degenerate_mad", &
                                  test_obtain_ensembles_reports_degenerate_mad)
        all_tests(58) = test_case("test_merge_ensembles_overlap_coefficient_regimes", &
                                  test_merge_ensembles_overlap_coefficient_regimes)
        all_tests(59) = test_case("test_merge_ensembles_clears_unused_columns", &
                                  test_merge_ensembles_clears_unused_columns)
        all_tests(60) = test_case("test_density_labels_tiles", test_density_labels_tiles)
        all_tests(61) = test_case("test_density_labels_invalid_tiles", test_density_labels_invalid_tiles)
        all_tests(62) = test_case("test_growth_commit_reference", test_growth_commit_reference)
        all_tests(63) = test_case("test_growth_commit_rejected_history", test_growth_commit_rejected_history)
        all_tests(64) = test_case("test_growth_incremental_surface", test_growth_incremental_surface)
        all_tests(65) = test_case("test_identify_ensemble_seeds_radii", test_identify_ensemble_seeds_radii)
        all_tests(66) = test_case("test_identify_ensemble_seeds_radii_percentile", &
                                  test_identify_ensemble_seeds_radii_percentile)
        all_tests(67) = test_case("test_obtain_ensembles_per_seed_radii", test_obtain_ensembles_per_seed_radii)
        all_tests(68) = test_case("test_observable_window_rejects_after_shift", test_observable_window_rejects_after_shift)
        all_tests(69) = test_case("test_history_requires_two_columns", test_history_requires_two_columns)
        all_tests(70) = test_case("test_zero_mad_error_positions", test_zero_mad_error_positions)
    end function get_all_tests_shatter_cluster_data

    subroutine test_density_radius_basic()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), tmp_mean_vec(n_dims), tmp_distances(n_vecs), radius
        integer(int32) :: tmp_perm(n_vecs), ierr, i_setup

        do i_setup = 1, n_vecs
            vectors(:, i_setup) = [real(i_setup, real64), real(i_setup, real64)]
        end do

        call calculate_density_radius_expert(vectors, n_dims, n_vecs, tmp_mean_vec, tmp_distances, tmp_perm, &
            0.50_real64, radius, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_density_radius_basic: ierr success control check")
        call assert_equal_real(radius, 3.5355339059327378_real64, 1.0e-12_real64, &
                               "test_density_radius_basic: exact analytical radius verification")
    end subroutine test_density_radius_basic

    subroutine test_density_labels_basic()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, 1)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, i_chk
        character(len=128) :: assert_msg

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.1_real64]
        vectors(:, 3) = [0.1_real64, 0.0_real64]
        vectors(:, 4) = [0.1_real64, 0.1_real64]
        vectors(:, 5) = [0.0_real64, 0.2_real64]
        vectors(:, 6) = [0.2_real64, 0.0_real64]
        vectors(:, 7) = [10.0_real64, 10.0_real64]
        vectors(:, 8) = [20.0_real64, 20.0_real64]
        vectors(:, 9) = [30.0_real64, 30.0_real64]
        vectors(:, 10) = [40.0_real64, 40.0_real64]

        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_basic: tree construction check")

        call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, 0.5_real64, dimension_order, kd_indices, &
            size(tmp_stack, 3, kind=int32), tmp_stack, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_basic: ierr execution check")

        do i_chk = 1, 6
            write (assert_msg, '(A,I0)') "Dense cluster tracking verification at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), 6.0_real64, 0.0_real64, assert_msg)
        end do

        do i_chk = 7, 10
            write (assert_msg, '(A,I0)') "Isolated point density verification at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), 1.0_real64, 0.0_real64, assert_msg)
        end do
    end subroutine test_density_labels_basic

    subroutine test_density_radius_invalid_quantile()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 5_int32

        real(real64) :: vectors(n_dims, n_vecs), tmp_mean_vec(n_dims), tmp_distances(n_vecs), radius
        integer(int32) :: tmp_perm(n_vecs), ierr, i_setup

        do i_setup = 1, n_vecs
            vectors(:, i_setup) = [real(i_setup, real64), real(i_setup, real64)]
        end do

        call calculate_density_radius_expert(vectors, n_dims, n_vecs, tmp_mean_vec, tmp_distances, tmp_perm, &
            -0.1_real64, radius, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_radius_invalid_quantile: Negative quantile must fail")

        call calculate_density_radius_expert(vectors, n_dims, n_vecs, tmp_mean_vec, tmp_distances, tmp_perm, &
            1.1_real64, radius, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_radius_invalid_quantile: Quantile > 1.0 must fail")
    end subroutine test_density_radius_invalid_quantile

    subroutine test_density_labels_invalid_r()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, 1), ierr

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 1.0_real64]
        vectors(:, 3) = [2.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, -1.0_real64, dimension_order, kd_indices, &
            size(tmp_stack, 3, kind=int32), tmp_stack, label_densities, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_labels_invalid_r: Negative search radius must fail")
    end subroutine test_density_labels_invalid_r

    subroutine test_density_labels_invalid_kd_indices()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, 1), ierr

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 1.0_real64]
        vectors(:, 3) = [2.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        kd_indices = [1_int32, 0_int32, 3_int32]
        call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, 0.5_real64, dimension_order, kd_indices, &
            size(tmp_stack, 3, kind=int32), tmp_stack, label_densities, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_labels_invalid_kd_indices: Zero index must fail")

        kd_indices = [1_int32, 2_int32, 4_int32]
        call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, 0.5_real64, dimension_order, kd_indices, &
            size(tmp_stack, 3, kind=int32), tmp_stack, label_densities, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_labels_invalid_kd_indices: Out of bounds index must fail")
    end subroutine test_density_labels_invalid_kd_indices

    subroutine test_density_radius_default_quantile()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), tmp_mean_vec(n_dims), tmp_distances(n_vecs)
        real(real64) :: radius_default, radius_explicit
        integer(int32) :: tmp_perm(n_vecs), ierr, i_setup

        do i_setup = 1, n_vecs
            vectors(:, i_setup) = [real(i_setup, real64), real(i_setup, real64)]
        end do

        call calculate_density_radius_expert(vectors, n_dims, n_vecs, tmp_mean_vec, tmp_distances, tmp_perm, &
            radius=radius_default, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_radius_default_quantile: default call check")

        call calculate_density_radius_expert(vectors, n_dims, n_vecs, tmp_mean_vec, tmp_distances, tmp_perm, &
            0.15_real64, radius_explicit, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_radius_default_quantile: explicit call check")

        call assert_equal_real(radius_default, radius_explicit, 1.0e-12_real64, &
                               "test_density_radius_default_quantile: default matches 0.15 option")
    end subroutine test_density_radius_default_quantile

    subroutine test_density_labels_small_radius()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, 1)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, i_chk
        character(len=128) :: assert_msg

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.1_real64]
        vectors(:, 3) = [0.1_real64, 0.0_real64]
        vectors(:, 4) = [0.1_real64, 0.1_real64]
        vectors(:, 5) = [0.0_real64, 0.2_real64]
        vectors(:, 6) = [0.2_real64, 0.0_real64]
        vectors(:, 7) = [10.0_real64, 10.0_real64]
        vectors(:, 8) = [20.0_real64, 20.0_real64]
        vectors(:, 9) = [30.0_real64, 30.0_real64]
        vectors(:, 10) = [40.0_real64, 40.0_real64]

        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_small_radius: tree construction check")

        call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, 0.05_real64, dimension_order, kd_indices, &
            size(tmp_stack, 3, kind=int32), tmp_stack, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_small_radius: ierr execution check")

        do i_chk = 1, 10
            write (assert_msg, '(A,I0)') "test_density_labels_small_radius: pruning verification check at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), 1.0_real64, 0.0_real64, assert_msg)
        end do
    end subroutine test_density_labels_small_radius

    subroutine test_density_labels_single_vector()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 1_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, 1)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr

        vectors(:, 1) = [1.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_single_vector: tree build check")

        call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, 0.5_real64, dimension_order, kd_indices, &
            size(tmp_stack, 3, kind=int32), tmp_stack, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_single_vector: execution check")

        call assert_equal_real(label_densities(1), 1.0_real64, 0.0_real64, &
                               "test_density_labels_single_vector: single point density tracking")
    end subroutine test_density_labels_single_vector

    subroutine test_density_labels_identical_vectors()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 5_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, 1)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, i_chk
        character(len=128) :: assert_msg

        do i_chk = 1, n_vecs
            vectors(:, i_chk) = [4.2_real64, 4.2_real64]
        end do
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_identical_vectors: tree build check")

        call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, 0.1_real64, dimension_order, kd_indices, &
            size(tmp_stack, 3, kind=int32), tmp_stack, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_identical_vectors: execution check")

        do i_chk = 1, n_vecs
            write (assert_msg, '(A,I0)') "test_density_labels_identical_vectors: identical coordinates check at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), real(n_vecs, real64), 0.0_real64, assert_msg)
        end do
    end subroutine test_density_labels_identical_vectors

    subroutine test_density_radius_plain()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), radius
        integer(int32) :: ierr, i_setup

        do i_setup = 1, n_vecs
            vectors(:, i_setup) = [real(i_setup, real64), real(i_setup, real64)]
        end do

        call calculate_density_radius(vectors, n_dims, n_vecs, 0.50_real64, radius, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_density_radius_plain: ierr success check")
        call assert_equal_real(radius, 3.5355339059327378_real64, 1.0e-12_real64, &
                               "test_density_radius_plain: analytical radius check")
    end subroutine test_density_radius_plain

    subroutine test_density_labels_plain()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, i_chk
        character(len=128) :: assert_msg

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.1_real64]
        vectors(:, 3) = [0.1_real64, 0.0_real64]
        vectors(:, 4) = [0.1_real64, 0.1_real64]
        vectors(:, 5) = [0.0_real64, 0.2_real64]
        vectors(:, 6) = [0.2_real64, 0.0_real64]
        vectors(:, 7) = [10.0_real64, 10.0_real64]
        vectors(:, 8) = [20.0_real64, 20.0_real64]
        vectors(:, 9) = [30.0_real64, 30.0_real64]
        vectors(:, 10) = [40.0_real64, 40.0_real64]

        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_plain: tree build check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, dimension_order, kd_indices, &
            min(8_int32, n_vecs), label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_plain: execution check")

        do i_chk = 1, 6
            write (assert_msg, '(A,I0)') "test_density_labels_plain: dense cluster check at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), 6.0_real64, 0.0_real64, assert_msg)
        end do

        do i_chk = 7, 10
            write (assert_msg, '(A,I0)') "test_density_labels_plain: isolated points check at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), 1.0_real64, 0.0_real64, assert_msg)
        end do
    end subroutine test_density_labels_plain

    subroutine test_density_labels_tiles()
        integer(int32), parameter :: n_dims = 2_int32, n_vecs = 11_int32
        integer(int32), parameter :: tile_counts(5) = [1_int32, 2_int32, 8_int32, n_vecs, n_vecs + 3_int32]
        real(real64), parameter :: radii(3) = [0.0_real64, 1.0_real64, 100.0_real64]
        real(real64) :: vectors(n_dims, n_vecs), labels(n_vecs), expected(n_vecs), tmp_values(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_workspace(n_vecs)
        integer(int32) :: tmp_perm(n_vecs), tmp_rec_stack(3, n_vecs)
        integer(int32), allocatable :: tmp_stack(:, :, :)
        integer(int32) :: ierr, i_vec, i_neighbor, i_radius, i_case, n_tiles

        vectors(1, :) = [4.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 2.0_real64, &
                         8.0_real64, 3.0_real64, 1.0_real64, 20.0_real64, 4.0_real64, 30.0_real64]
        vectors(2, :) = 0.0_real64
        dimension_order = [1_int32, 2_int32]
        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_values, tmp_perm, tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "density tiles: tree construction")

        do i_radius = 1, size(radii)
            expected = 0.0_real64
            do i_vec = 1, n_vecs
                do i_neighbor = 1, n_vecs
                    if (sum((vectors(:, i_vec) - vectors(:, i_neighbor))**2) <= radii(i_radius)**2) then
                        expected(i_vec) = expected(i_vec) + 1.0_real64
                    end if
                end do
            end do

            call calculate_labels_as_density(vectors, n_dims, n_vecs, radii(i_radius), dimension_order, kd_indices, &
                min(8_int32, n_vecs), labels, ierr)
            call assert_equal_int(ierr, ERR_OK, "density tiles: default allocation")
            call assert_true(all(labels == expected), "density tiles: default matches all-pairs counts")

            do i_case = 1, size(tile_counts)
                n_tiles = tile_counts(i_case)
                allocate (tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, n_tiles))
                tmp_stack = -99_int32
                labels = -1.0_real64
                call calculate_labels_as_density_expert(vectors, n_dims, n_vecs, radii(i_radius), dimension_order, &
                    kd_indices, n_tiles, tmp_stack, labels, ierr)
                call assert_equal_int(ierr, ERR_OK, "density tiles: explicit workspace")
                call assert_true(all(labels == expected), "density tiles: explicit matches all-pairs counts")
                if (n_tiles > n_vecs) then
                    call assert_true(all(tmp_stack(:, :, n_vecs + 1:) == -99_int32), &
                                     "density tiles: idle stacks remain untouched")
                end if
                deallocate (tmp_stack)

                labels = -1.0_real64
                call calculate_labels_as_density(vectors, n_dims, n_vecs, radii(i_radius), dimension_order, kd_indices, &
                    n_tiles, labels, ierr)
                call assert_equal_int(ierr, ERR_OK, "density tiles: requested allocation")
                call assert_true(all(labels == expected), "density tiles: allocation matches all-pairs counts")
            end do
        end do

        call calculate_labels_as_density(vectors(:, 1:1), n_dims, 1_int32, 0.0_real64, dimension_order, [1_int32], &
            huge(1_int32), labels(1:1), ierr)
        call assert_equal_int(ierr, ERR_OK, "density tiles: allocation clamps excess tiles")
        call assert_equal_real(labels(1), 1.0_real64, 0.0_real64, "density tiles: singleton allocation")
    end subroutine test_density_labels_tiles

    subroutine test_density_labels_invalid_tiles()
        real(real64) :: vectors(1, 1), labels(1)
        integer(int32) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, 1)
        integer(int32) :: ierr, n_tiles

        vectors = 0.0_real64
        do n_tiles = -1_int32, 0_int32
            call calculate_labels_as_density(vectors, 1_int32, 1_int32, 0.0_real64, [1_int32], [1_int32], n_tiles, &
                labels, ierr)
            call assert_true(ierr /= ERR_OK, "density tiles: allocating entry rejects nonpositive tiles")
            call calculate_labels_as_density_expert(vectors, 1_int32, 1_int32, 0.0_real64, [1_int32], [1_int32], &
                n_tiles, tmp_stack, labels, ierr)
            call assert_true(ierr /= ERR_OK, "density tiles: validated entry rejects nonpositive tiles")
        end do
    end subroutine test_density_labels_invalid_tiles

    subroutine test_growth_commit_reference()
        integer(int32), parameter :: n_dims = 2_int32, n_vecs = 9_int32
        integer(int32), parameter :: windows(3) = [2_int32, 3_int32, 0_int32]
        integer(int32), parameter :: seeds(3) = [1_int32, 4_int32, 7_int32]
        real(real64), parameter :: thresholds(4) = [0.0_real64, 0.1_real64, 0.5_real64, 10.0_real64]
        real(real64), parameter :: mad_factors(3) = [0.0_real64, 1.0_real64, 100.0_real64]
        real(real64) :: vectors(n_dims, n_vecs), densities(n_vecs), tmp_values(n_vecs)
        real(real64) :: tmp_abs_diff(n_vecs), median_ambient, mad_ambient
        real(real64) :: member_densities(n_vecs), center_density, swap_density
        real(real64), allocatable :: history(:, :), reference_history(:, :)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_workspace(n_vecs)
        integer(int32) :: tmp_perm(n_vecs), tmp_rec_stack(3, n_vecs)
        integer(int32) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        integer(int32) :: ierr, i_seed, i_window, i_threshold, i_factor, n_cols
        integer(int32) :: i_vec, i_member, i_sort, n_members, n_candidates, i_iter, stop_reason, reference_stop
        logical(c_bool) :: tmp_vicinity(n_vecs), tmp_surface(n_vecs), current_mask(n_vecs), result_mask(n_vecs)
        logical(c_bool) :: reference_mask(n_vecs), previous_mask(n_vecs), accepted
        logical :: saw_fixed, saw_isolated, saw_first_reject, saw_later_reject

        vectors(1, :) = [0.0_real64, 0.4_real64, 0.8_real64, 1.2_real64, 1.6_real64, &
                         2.0_real64, 10.0_real64, 20.0_real64, 30.0_real64]
        vectors(2, :) = 0.0_real64
        densities = [4.0_real64, 4.0_real64, 4.0_real64, 16.0_real64, 16.0_real64, &
                     20.0_real64, 2.0_real64, 1.0_real64, 3.0_real64]
        dimension_order = [1_int32, 2_int32]
        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_values, tmp_perm, tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "commit reference: build tree")
        call compute_ambient_density_stats_helper(densities, n_vecs, tmp_perm, tmp_abs_diff, &
                                                  median_ambient, mad_ambient)
        saw_fixed = .false.
        saw_isolated = .false.
        saw_first_reject = .false.
        saw_later_reject = .false.

        do i_window = 1, size(windows)
            n_cols = windows(i_window)
            if (n_cols == 0_int32) n_cols = n_vecs
            allocate (history(5, n_cols), reference_history(5, n_cols))
            do i_factor = 1, size(mad_factors)
                do i_threshold = 1, size(thresholds)
                    do i_seed = 1, size(seeds)

                        call grow_single_seed_helper(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                            densities, seeds(i_seed), 0.5_real64, mad_factors(i_factor), &
                            mad_ambient, thresholds(i_threshold), &
                            tmp_stack, tmp_vicinity, tmp_surface, tmp_perm, tmp_abs_diff, &
                            history, current_mask, result_mask, stop_reason)
                        reference_mask = .false.
                        reference_mask(seeds(i_seed)) = .true.
                        reference_history = 0.0_real64
                        i_iter = 1_int32
                        call compute_ensemble_observable_helper(reference_mask, densities, n_vecs, &
                            0_int32, i_iter, reference_history)
                        do
                            previous_mask = reference_mask
                            n_members = 0_int32
                            do i_vec = 1, n_vecs
                                if (.not. previous_mask(i_vec)) cycle
                                n_members = n_members + 1_int32
                                member_densities(n_members) = densities(i_vec)
                            end do

                            do i_vec = 2, n_members
                                i_sort = i_vec
                                do while (i_sort > 1_int32)
                                    if (member_densities(i_sort - 1) <= member_densities(i_sort)) exit
                                    swap_density = member_densities(i_sort - 1)
                                    member_densities(i_sort - 1) = member_densities(i_sort)
                                    member_densities(i_sort) = swap_density
                                    i_sort = i_sort - 1_int32
                                end do
                            end do
                            center_density = (member_densities((n_members + 1)/2) + &
                                              member_densities(n_members/2 + 1))/2.0_real64
                            do i_vec = 1, n_vecs
                                if (previous_mask(i_vec)) cycle
                                if (abs(densities(i_vec) - center_density) > mad_factors(i_factor)*mad_ambient) cycle
                                do i_member = 1, n_vecs
                                    if (.not. previous_mask(i_member)) cycle
                                    if (sum((vectors(:, i_vec) - vectors(:, i_member))**2) <= 0.25_real64) then
                                        reference_mask(i_vec) = .true.
                                        exit
                                    end if
                                end do
                            end do
                            n_candidates = count(reference_mask) - n_members
                            if (n_candidates == 0_int32) then
                                reference_stop = STC_STOP_FIXED_POINT
                                if (i_iter == 1_int32) reference_stop = STC_STOP_NO_CANDIDATES
                                exit
                            end if
                            i_iter = i_iter + 1_int32
                            call compute_ensemble_observable_helper(reference_mask, densities, n_vecs, &
                                n_candidates, i_iter, reference_history)
                            call accept_ensemble_helper(reference_history, i_iter, thresholds(i_threshold), accepted)
                            if (.not. accepted) then
                                reference_mask = previous_mask
                                reference_stop = STC_STOP_REJECT_AFTER_ACCEPT
                                if (i_iter == 2_int32) reference_stop = STC_STOP_NEVER_ACCEPTED
                                exit
                            end if
                        end do
                        call assert_true(logical(all(result_mask .eqv. reference_mask)), "commit reference: final membership")
                        call assert_true(logical(all(current_mask .eqv. reference_mask)), "commit reference: accepted state")
                        call assert_equal_int(stop_reason, reference_stop, "commit reference: stop reason")
                        call assert_true(all(abs(history - reference_history) < 1.0e-12_real64), &
                                         "commit reference: all five observables, including rejected trial and window shifts")
                        saw_fixed = saw_fixed .or. stop_reason == STC_STOP_FIXED_POINT
                        saw_isolated = saw_isolated .or. stop_reason == STC_STOP_NO_CANDIDATES
                        saw_first_reject = saw_first_reject .or. stop_reason == STC_STOP_NEVER_ACCEPTED
                        saw_later_reject = saw_later_reject .or. stop_reason == STC_STOP_REJECT_AFTER_ACCEPT
                    end do
                end do
            end do
            deallocate (history, reference_history)
        end do
        call assert_true(saw_fixed .and. saw_isolated .and. saw_first_reject .and. saw_later_reject, &
                         "commit reference: fixture exercises all four stop reasons")
    end subroutine test_growth_commit_reference

    subroutine test_growth_commit_rejected_history()
        integer(int32), parameter :: n_dims = 1_int32, n_vecs = 6_int32
        real(real64) :: vectors(n_dims, n_vecs), densities(n_vecs), tmp_values(n_vecs), tmp_abs_diff(n_vecs)
        real(real64) :: history(5, 3), expected_trial(5)
        integer(int32) :: kd_indices(n_vecs), tmp_workspace(n_vecs), tmp_perm(n_vecs), tmp_rec_stack(3, n_vecs)
        integer(int32) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH), ierr, stop_reason
        logical(c_bool) :: tmp_vicinity(n_vecs), tmp_surface(n_vecs), current_mask(n_vecs), result_mask(n_vecs)

        vectors(1, :) = [0.0_real64, 0.4_real64, 0.8_real64, 1.2_real64, 10.0_real64, 20.0_real64]
        densities = [4.0_real64, 4.0_real64, 4.0_real64, 16.0_real64, 1.0_real64, 2.0_real64]
        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, [1_int32], &
                                   tmp_workspace, tmp_values, tmp_perm, tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "commit history: build tree")
        call grow_single_seed_helper(vectors, n_dims, n_vecs, [1_int32], kd_indices, densities, &
                                     1_int32, 0.5_real64, 100.0_real64, 1.0_real64, 0.1_real64, &
                                     tmp_stack, tmp_vicinity, tmp_surface, tmp_perm, tmp_abs_diff, &
                                     history, current_mask, result_mask, stop_reason)
        call assert_equal_int(stop_reason, STC_STOP_REJECT_AFTER_ACCEPT, "commit history: late rejection")
        call assert_true(logical(all(result_mask(1:3)) .and. .not. any(result_mask(4:))), &
                         "commit history: exactly the previously accepted chain remains")
        call assert_true(logical(all(current_mask .eqv. result_mask)), "commit history: accepted workspace unchanged by rejection")
        call assert_true(logical(all(tmp_surface(1:4)) .and. .not. any(tmp_surface(5:))), &
                         "commit history: spatial union includes the rejected candidate")
        expected_trial = [7.0_real64, 64.0_real64/13.0_real64, 91.0_real64/64.0_real64, 4.0_real64, 1.0_real64]
        call assert_true(all(abs(history(:, 3) - expected_trial) < 1.0e-12_real64), &
                         "commit history: rejected trial observables are preserved after window shift")
        call assert_true(all(history(4, :) == [2.0_real64, 3.0_real64, 4.0_real64]), &
                         "commit history: rolling membership counts include accepted and rejected steps")
    end subroutine test_growth_commit_rejected_history

    subroutine test_growth_incremental_surface()
        integer(int32), parameter :: n_dims = 1_int32, n_vecs = 6_int32
        real(real64) :: vectors(n_dims, n_vecs), densities(n_vecs), tmp_values(n_vecs), tmp_abs_diff(n_vecs)
        real(real64) :: history(5, n_vecs)
        integer(int32) :: kd_indices(n_vecs), tmp_workspace(n_vecs), tmp_perm(n_vecs), tmp_rec_stack(3, n_vecs)
        integer(int32) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH), ierr, stop_reason
        logical(c_bool) :: tmp_vicinity(n_vecs), tmp_surface(n_vecs), current_mask(n_vecs), result_mask(n_vecs)

        vectors(1, :) = [0.0_real64, 0.9_real64, -0.9_real64, -1.8_real64, 0.1_real64, 10.0_real64]
        densities = [10.0_real64, 14.0_real64, 16.0_real64, 16.0_real64, 1.0_real64, 2.0_real64]
        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, [1_int32], &
                                   tmp_workspace, tmp_values, tmp_perm, tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "incremental surface: build tree")
        tmp_surface = .true.
        current_mask = .true.
        call grow_single_seed_helper(vectors, n_dims, n_vecs, [1_int32], kd_indices, densities, &
                                     1_int32, 1.0_real64, 4.0_real64, 1.0_real64, 10.0_real64, &
                                     tmp_stack, tmp_vicinity, tmp_surface, tmp_perm, tmp_abs_diff, &
                                     history, current_mask, result_mask, stop_reason)
        call assert_equal_int(stop_reason, STC_STOP_FIXED_POINT, "incremental surface: reaches fixed point")
        call assert_true(logical(all(result_mask(1:4)) .and. .not. any(result_mask(5:))), &
                         "incremental surface: reevaluates cached neighbors while retaining accepted members")
        call assert_true(logical(all(tmp_surface(1:5)) .and. .not. tmp_surface(6)), &
                         "incremental surface: cache retains incompatible spatial neighbors and excludes distant points")
        call assert_true(all(history(4, 1:4) == [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]), &
                         "incremental surface: newly accepted members expand geometry only in the following round")
        call assert_true(all(history(5, 1:4) == [0.0_real64, 1.0_real64, 1.0_real64, 1.0_real64]), &
                         "incremental surface: exactly one new candidate per round")
        call assert_true(all(history(:, 5:) == 0.0_real64), "incremental surface: no spurious extra growth")
    end subroutine test_growth_incremental_surface

    subroutine test_identify_ensemble_seeds_basic()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: k_seeding = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_distances(n_vecs)
        real(real64) :: tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds, ierr
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs)
        logical(c_bool) :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)
        real(real64) :: seed_radii(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [10.0_real64, 0.0_real64]
        vectors(:, 5) = [11.0_real64, 0.0_real64]
        vectors(:, 6) = [12.0_real64, 0.0_real64]

        density_labels = [60.0_real64, 50.0_real64, 40.0_real64, &
                          30.0_real64, 20.0_real64, 10.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_basic: tree build")

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_basic: execution")
        call assert_true(all(sorted_perm == [1_int32, 2_int32, 3_int32, &
                                             4_int32, 5_int32, 6_int32]), &
                         "test_identify_ensemble_seeds_basic: density ranking")
        call assert_equal_int(n_seeds, 2_int32, &
                              "test_identify_ensemble_seeds_basic: one seed per separated region")
        call assert_true(logical(seed_mask(1)), &
                         "test_identify_ensemble_seeds_basic: densest point in first region is seed")
        call assert_false(logical(seed_mask(2)), &
                          "test_identify_ensemble_seeds_basic: covered dense point is skipped")
        call assert_false(logical(seed_mask(3)), &
                          "test_identify_ensemble_seeds_basic: covered first-region point is skipped")
        call assert_true(logical(seed_mask(4)), &
                         "test_identify_ensemble_seeds_basic: densest uncovered second region becomes seed")
        call assert_false(logical(seed_mask(5)), &
                          "test_identify_ensemble_seeds_basic: covered second-region point is skipped")
        call assert_false(logical(seed_mask(6)), &
                          "test_identify_ensemble_seeds_basic: covered second-region point is skipped")
        call assert_true(logical(all(tmp_visited_mask)), &
                         "test_identify_ensemble_seeds_basic: all vectors are visited by seeding")
        call assert_equal_int(count(seed_mask, kind=int32), n_seeds, &
                              "test_identify_ensemble_seeds_basic: seed count matches seed mask")
    end subroutine test_identify_ensemble_seeds_basic

    subroutine test_identify_ensemble_seeds_k_seeding_effect()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_distances(n_vecs)
        real(real64) :: tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds_small, n_seeds_large, ierr
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs)
        logical(c_bool) :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs)
        logical(c_bool) :: seed_mask_small(n_vecs), seed_mask_large(n_vecs)
        real(real64) :: seed_radii(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [10.0_real64, 0.0_real64]
        vectors(:, 5) = [11.0_real64, 0.0_real64]
        vectors(:, 6) = [12.0_real64, 0.0_real64]

        density_labels = [60.0_real64, 50.0_real64, 40.0_real64, &
                          30.0_real64, 20.0_real64, 10.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_k_seeding_effect: tree build")

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            1_int32, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds_small, seed_mask_small, seed_radii, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=1 execution")
        call assert_equal_int(n_seeds_small, 4_int32, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=1 seed count")
        call assert_true(logical(seed_mask_small(1)), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 first seed")
        call assert_true(logical(seed_mask_small(3)), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 uncovered first-region tail")
        call assert_true(logical(seed_mask_small(4)), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 second-region seed")
        call assert_true(logical(seed_mask_small(6)), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 uncovered second-region tail")

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            3_int32, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds_large, seed_mask_large, seed_radii, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=3 execution")
        call assert_equal_int(n_seeds_large, 2_int32, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=3 seed count")
        call assert_true(logical(seed_mask_large(1)), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=3 first-region seed")
        call assert_true(logical(seed_mask_large(4)), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=3 second-region seed")
        call assert_true(n_seeds_large < n_seeds_small, &
                         "test_identify_ensemble_seeds_k_seeding_effect: larger k gives broader coverage")
    end subroutine test_identify_ensemble_seeds_k_seeding_effect

    subroutine test_identify_ensemble_seeds_invalid_k_seeding()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_distances(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds, ierr
        logical(c_bool) :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)
        real(real64) :: seed_radii(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [3.0_real64, 0.0_real64]

        density_labels = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32, 4_int32]

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            0_int32, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_k_seeding: k=0 must fail")

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            n_vecs, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_k_seeding: k>=n_vectors must fail")
    end subroutine test_identify_ensemble_seeds_invalid_k_seeding

    subroutine test_identify_ensemble_seeds_plain()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: k_seeding = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), sorted_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), n_seeds, ierr
        logical(c_bool) :: seed_mask(n_vecs)
        real(real64) :: seed_radii(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [10.0_real64, 0.0_real64]
        vectors(:, 5) = [11.0_real64, 0.0_real64]
        vectors(:, 6) = [12.0_real64, 0.0_real64]

        density_labels = [60.0_real64, 50.0_real64, 40.0_real64, &
                          30.0_real64, 20.0_real64, 10.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_plain: tree build")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, sorted_perm, n_seeds, seed_mask, seed_radii, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_plain: execution")
        call assert_equal_int(n_seeds, 2_int32, &
                              "test_identify_ensemble_seeds_plain: expected two seeds")
        call assert_true(logical(seed_mask(1)), &
                         "test_identify_ensemble_seeds_plain: first-region density centre is seed")
        call assert_true(logical(seed_mask(4)), &
                         "test_identify_ensemble_seeds_plain: second-region density centre is seed")
        call assert_false(logical(seed_mask(2)), &
                          "test_identify_ensemble_seeds_plain: covered first-region point is skipped")
        call assert_false(logical(seed_mask(3)), &
                          "test_identify_ensemble_seeds_plain: covered first-region point is skipped")
        call assert_false(logical(seed_mask(5)), &
                          "test_identify_ensemble_seeds_plain: covered second-region point is skipped")
        call assert_false(logical(seed_mask(6)), &
                          "test_identify_ensemble_seeds_plain: covered second-region point is skipped")
    end subroutine test_identify_ensemble_seeds_plain

    subroutine test_identify_ensemble_seeds_single_vector_invalid()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 1_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_distances(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds, ierr
        logical(c_bool) :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)
        real(real64) :: seed_radii(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        density_labels = [1.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32]

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            1_int32, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)

        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_single_vector_invalid: nearest neighbors require n_vectors>=2")
    end subroutine test_identify_ensemble_seeds_single_vector_invalid

    subroutine test_identify_ensemble_seeds_invalid_inputs()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: k_seeding = 2_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_distances(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds, ierr
        logical(c_bool) :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)
        real(real64) :: seed_radii(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [3.0_real64, 0.0_real64]

        density_labels = [4.0_real64, 3.0_real64, 2.0_real64, -1.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32, 4_int32]

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_inputs: negative density must fail")

        density_labels = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        kd_indices = [0_int32, 2_int32, 3_int32, 4_int32]

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_inputs: invalid kd index must fail")

        kd_indices = [1_int32, 2_int32, 3_int32, 4_int32]
        dimension_order = [0_int32, 2_int32]

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_inputs: invalid dimension order must fail")
    end subroutine test_identify_ensemble_seeds_invalid_inputs

    subroutine test_identify_ensemble_seeds_identical_vectors()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: k_seeding = 2_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_distances(n_vecs)
        real(real64) :: tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds, ierr
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs)
        logical(c_bool) :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)
        real(real64) :: seed_radii(n_vecs)

        vectors = 5.0_real64
        density_labels = [1.0_real64, 4.0_real64, 3.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_identical_vectors: tree build")

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_identical_vectors: execution")
        call assert_equal_int(sorted_perm(1), 2_int32, &
                              "test_identify_ensemble_seeds_identical_vectors: densest vector ranks first")
        call assert_equal_int(n_seeds, 1_int32, &
                              "test_identify_ensemble_seeds_identical_vectors: zero-radius coverage gives one seed")
        call assert_true(logical(seed_mask(2)), &
                         "test_identify_ensemble_seeds_identical_vectors: densest vector is the seed")
        call assert_false(logical(seed_mask(1)), &
                          "test_identify_ensemble_seeds_identical_vectors: coincident point is covered")
        call assert_false(logical(seed_mask(3)), &
                          "test_identify_ensemble_seeds_identical_vectors: coincident point is covered")
        call assert_false(logical(seed_mask(4)), &
                          "test_identify_ensemble_seeds_identical_vectors: coincident point is covered")
        call assert_true(logical(all(tmp_visited_mask)), &
                         "test_identify_ensemble_seeds_identical_vectors: all coincident points visited")
    end subroutine test_identify_ensemble_seeds_identical_vectors

    subroutine test_grow_ensemble_basic()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs), tmp_surface_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.2_real64]
        vectors(:, 4) = [0.0_real64, 0.4_real64]

        density_labels = [5.0_real64, 5.0_real64, 8.0_real64, 1.0_real64]
        ensemble_mask = [.true., .false., .false., .false.]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_basic: tree construction check")

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_basic: execution ierr check")
        call assert_true(logical(ensemble_mask(1)), "test_grow_ensemble_basic: vector 1 must remain a member")
        call assert_true(logical(ensemble_mask(2)), "test_grow_ensemble_basic: vector 2 must be added")
        call assert_false(logical(ensemble_mask(3)), "test_grow_ensemble_basic: vector 3 too far away")
        call assert_false(logical(ensemble_mask(4)), "test_grow_ensemble_basic: vector 4 density too low")
    end subroutine test_grow_ensemble_basic

    subroutine test_grow_ensemble_no_growth()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs), tmp_surface_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [100.0_real64, 1.0_real64, 2.0_real64]
        ensemble_mask = [.true., .false., .false.]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_no_growth: execute ierr success")
        call assert_true(logical(ensemble_mask(1)), "test_grow_ensemble_no_growth: vector 1 stays member")
        call assert_false(logical(ensemble_mask(2)), "test_grow_ensemble_no_growth: vector 2 rejected due to density gap")
        call assert_false(logical(ensemble_mask(3)), "test_grow_ensemble_no_growth: vector 3 rejected")
    end subroutine test_grow_ensemble_no_growth

    subroutine test_grow_ensemble_invalid_inputs()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64), tmp_perm(n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs), tmp_surface_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false.]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, -0.1_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(ierr /= ERR_OK, "test_grow_ensemble_invalid_inputs: negative r must fail")

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.5_real64, dimension_order, kd_indices, &
            density_labels, -1.0_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(ierr /= ERR_OK, "test_grow_ensemble_invalid_inputs: negative alpha_mad must fail")
    end subroutine test_grow_ensemble_invalid_inputs

    subroutine test_grow_ensemble_plain()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)
        real(real64) :: density_radius

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]
        vectors(:, 4) = [100.0_real64, 100.0_real64]
        vectors(:, 5) = [200.0_real64, 200.0_real64]
        vectors(:, 6) = [300.0_real64, 300.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64, &
                          1.0_real64, 2.0_real64, 3.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_plain: build tree check")

        call calculate_density_radius(vectors, n_dims, n_vecs, radius=density_radius, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "density radius calculation")

        ensemble_mask = .false.
        ensemble_mask(1) = .true.
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_plain: explicit r and alpha_mad ierr check")
        call assert_true(logical(ensemble_mask(1)), "Case 1: vector 1 is active")
        call assert_true(logical(ensemble_mask(2)), "Case 1: vector 2 is successfully added")
        call assert_false(logical(ensemble_mask(3)), "Case 1: vector 3 is too far away")

        ensemble_mask = .false.
        ensemble_mask(1) = .true.
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, density_radius, dimension_order, kd_indices, &
            density_labels, 0.5_real64, grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_plain: computed radius ierr check")

        ensemble_mask = .false.
        ensemble_mask(1) = .true.
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, grown_mask=grown_mask, ierr=ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_plain: omitted alpha_mad ierr check")

        ensemble_mask = .false.
        ensemble_mask(1) = .true.
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, density_radius, dimension_order, kd_indices, &
            density_labels, grown_mask=grown_mask, ierr=ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_plain: computed radius and default MAD ierr check")
    end subroutine test_grow_ensemble_plain

    subroutine test_grow_ensemble_multistep()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]
        vectors(:, 4) = [100.0_real64, 100.0_real64]
        vectors(:, 5) = [200.0_real64, 200.0_real64]
        vectors(:, 6) = [300.0_real64, 300.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64, &
                          1.0_real64, 2.0_real64, 3.0_real64]
        ensemble_mask = .false.
        ensemble_mask(1) = .true.
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(logical(ensemble_mask(2)), "step 1: vector 2 absorbed")
        call assert_false(logical(ensemble_mask(3)), "step 1: vector 3 not yet reachable")

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(logical(ensemble_mask(3)), "step 2: vector 3 pulled in by vector 2 surface")
    end subroutine test_grow_ensemble_multistep

    subroutine test_grow_ensemble_alpha_sensitivity()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs), tmp_surface_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [10.0_real64, 8.5_real64, 12.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.25_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, "alpha 0.25 test ierr control")
        call assert_false(logical(ensemble_mask(2)), "alpha 0.25 should reject vector 2")

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 1.00_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, "alpha 1.00 test ierr control")
        call assert_true(logical(ensemble_mask(2)), "alpha 1.00 should accept vector 2")
    end subroutine test_grow_ensemble_alpha_sensitivity

    subroutine test_grow_ensemble_zero_mad()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs), tmp_surface_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false.]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask

        call assert_true(ierr /= ERR_OK, &
                         "test_grow_ensemble_zero_mad: zero ambient MAD must be rejected")

        ensemble_mask = [.true., .false., .false.]

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(ierr /= ERR_OK, &
                         "test_grow_ensemble_zero_mad: allocating wrapper rejects a zero ambient MAD")

        density_labels = [10.0_real64, 6.0_real64, 2.0_real64]
        ensemble_mask = [.true., .false., .false.]

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_zero_mad: spread densities accepted")
        call assert_true(logical(ensemble_mask(1)), "test_grow_ensemble_zero_mad: vector 1 active")
    end subroutine test_grow_ensemble_zero_mad

    subroutine test_grow_ensemble_mask_states()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs), tmp_surface_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 6.0_real64, 9.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        ensemble_mask = [.false., .false., .false.]
        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, "empty mask check ierr control")
        call assert_false(logical(any(ensemble_mask)), "empty mask should stay empty")

        ensemble_mask = [.true., .true., .true.]
        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, "full mask check ierr control")
        call assert_true(logical(all(ensemble_mask)), "full mask should stay full")
    end subroutine test_grow_ensemble_mask_states

    subroutine test_compute_ensemble_observable_basic()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 10)
        logical(c_bool) :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 5.0_real64, 20.0_real64]
        ensemble_mask = [.true., .true., .false.]
        observables = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 2_int32, 1_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)

        call assert_equal_int(ierr, ERR_OK, "basic observable ierr control")
        call assert_equal_real(observables(1, 1), 7.5_real64, 1.0e-10_real64, "rho_arith metric check")
        call assert_equal_real(observables(2, 1), 6.666666666666667_real64, 1.0e-10_real64, "rho_harm metric check")
        call assert_equal_real(observables(3, 1), 1.125_real64, 1.0e-10_real64, "H_rho ratio check")
        call assert_equal_real(observables(4, 1), 2.0_real64, 1.0e-10_real64, "active size check")
        call assert_equal_real(observables(5, 1), 2.0_real64, 1.0e-10_real64, "candidate count recorded in row 5")
    end subroutine test_compute_ensemble_observable_basic

    subroutine test_compute_ensemble_observable_sliding_window()
        integer(int32), parameter :: n_vecs = 4_int32
        real(real64) :: density_labels(n_vecs), observables(5, 2)
        logical(c_bool) :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 5.0_real64, 20.0_real64, 2.0_real64]
        observables = 0.0_real64

        ensemble_mask = [.true., .true., .false., .false.]
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 1_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)

        ensemble_mask = [.true., .true., .true., .false.]
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 1_int32, 2_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)

        call assert_equal_real(observables(1, 1), 7.5_real64, 1.0e-10_real64, "col 1 has iter 1")
        call assert_equal_real(observables(1, 2), 11.666666666666666_real64, 1.0e-10_real64, "col 2 has iter 2")

        ensemble_mask = [.true., .true., .true., .true.]
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 1_int32, 3_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)

        call assert_equal_int(ierr, ERR_OK, "sliding window ierr control")
        call assert_equal_real(observables(1, 1), 11.666666666666666_real64, 1.0e-10_real64, "col 1 shifted to iter 2")
        call assert_equal_real(observables(1, 2), 9.25_real64, 1.0e-10_real64, "col 2 updated to iter 3")
    end subroutine test_compute_ensemble_observable_sliding_window

    subroutine test_compute_ensemble_observable_history_width()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 3)
        logical(c_bool) :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64]
        ensemble_mask = [.true., .false., .false.]
        observables = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 1_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)
        call assert_equal_int(ierr, ERR_OK, "history width: iter 1 ierr check")
        call assert_equal_real(observables(4, 1), 1.0_real64, 1.0e-10_real64, "iter 1 stored in col 1")

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 2_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)
        call assert_equal_int(ierr, ERR_OK, "history width: iter 2 ierr check")
        call assert_equal_real(observables(4, 2), 1.0_real64, 1.0e-10_real64, "iter 2 stored in col 2")
    end subroutine test_compute_ensemble_observable_history_width

    !> Regression: the observable writer and the acceptance reader share one window, the history width.
    !| Before, the writer wrapped at `t_observables` (default 10) while the reader wrapped at the width,
    !| so a history wider than 10 compared never-written zero columns and accepted every later step.
    subroutine test_observable_window_rejects_after_shift()
        integer(int32), parameter :: n_vecs = 4_int32
        real(real64) :: density_labels(n_vecs), narrow(5, 3), wide(5, 12)
        logical(c_bool) :: ensemble_mask(n_vecs), is_accepted
        integer(int32) :: ierr, i_iter

        ! H_rho stays 1 for the first three members and jumps to 7/3 when the fourth joins.
        density_labels = [1.0_real64, 1.0_real64, 1.0_real64, 9.0_real64]
        narrow = 0.0_real64
        wide = 0.0_real64

        do i_iter = 1, 4
            ensemble_mask = .false.
            ensemble_mask(1:i_iter) = .true.
            call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 1_int32, i_iter, narrow, &
                size(narrow, 2, kind=int32), size(narrow, 1, kind=int32), ierr)
            call assert_equal_int(ierr, ERR_OK, "window regression: narrow update")
            call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 1_int32, i_iter, wide, &
                size(wide, 2, kind=int32), size(wide, 1, kind=int32), ierr)
            call assert_equal_int(ierr, ERR_OK, "window regression: wide update")
        end do

        ! Iteration 4 is past the 3-column window, so the history has shifted.
        call accept_ensemble(narrow, size(narrow, 2, kind=int32), size(narrow, 1, kind=int32), &
            4_int32, 0.5_real64, is_accepted, ierr)
        call assert_equal_int(ierr, ERR_OK, "window regression: narrow accept ierr")
        call assert_false(logical(is_accepted), "shifted window still compares the last two steps")

        call accept_ensemble(wide, size(wide, 2, kind=int32), size(wide, 1, kind=int32), &
            4_int32, 0.5_real64, is_accepted, ierr)
        call assert_equal_int(ierr, ERR_OK, "window regression: wide accept ierr")
        call assert_false(logical(is_accepted), "unshifted window compares the same two steps")

        ! A history wider than the old default window of 10 keeps writing where acceptance reads.
        wide = 0.0_real64
        density_labels = [1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64]
        do i_iter = 1, 10
            ensemble_mask = .true.
            call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 1_int32, i_iter, wide, &
                size(wide, 2, kind=int32), size(wide, 1, kind=int32), ierr)
        end do
        density_labels(4) = 9.0_real64
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 1_int32, 11_int32, wide, &
            size(wide, 2, kind=int32), size(wide, 1, kind=int32), ierr)
        call accept_ensemble(wide, size(wide, 2, kind=int32), size(wide, 1, kind=int32), &
            11_int32, 0.5_real64, is_accepted, ierr)
        call assert_equal_int(ierr, ERR_OK, "window regression: iteration 11 ierr")
        call assert_false(logical(is_accepted), "iteration 11 is compared against iteration 10, not an empty column")
    end subroutine test_observable_window_rejects_after_shift

    !> A zero ambient MAD is reported against `density_labels`: at its own position by a direct call
    !| to the implementation, and with the position cleared by the generated entry points.
    subroutine test_zero_mad_error_positions()
        integer(int32), parameter :: n_dims = 1_int32, n_vecs = 3_int32, n_seeds = 1_int32
        real(real64) :: vectors(n_dims, n_vecs), labels(n_vecs), abs_diff(n_vecs, 1), obs(5, 2, 1), mad
        integer(int32) :: dims(n_dims), kd(n_vecs), seeds(n_seeds), stack(3, 64, 1), perm(n_vecs, 1)
        integer(int32) :: stops(n_seeds), ierr
        logical(c_bool) :: mask(n_vecs), grown(n_vecs), vicinity(n_vecs, 1), surface(n_vecs, 1)
        logical(c_bool) :: current(n_vecs, 1), ensembles(n_vecs, n_seeds)

        vectors(1, :) = [0.0_real64, 1.0_real64, 2.0_real64]
        labels = 5.0_real64
        dims = 1_int32
        kd = [1_int32, 2_int32, 3_int32]
        seeds = 1_int32
        mask = [.true., .false., .false.]

        call grow_ensemble_impl(vectors, n_dims, n_vecs, mask, 1.0_real64, dims, kd, labels, &
            0.5_real64, stack(:, :, 1), vicinity(:, 1), surface(:, 1), perm(:, 1), abs_diff(:, 1), grown, ierr)
        call assert_equal_int(get_err_arg_pos(ierr), 8_int32, "direct growth names density_labels argument")

        call obtain_ensembles_impl(vectors, n_dims, n_vecs, dims, kd, labels, seeds, n_seeds, &
            MODE_GLOBAL_RADIUS, r=1.0_real64, t_observables=2_int32, n_tiles=1_int32, tmp_stack=stack, &
            tmp_vicinity_mask=vicinity, tmp_surface_mask=surface, tmp_perm=perm, tmp_abs_diff=abs_diff, &
            tmp_observables=obs, tmp_current_mask=current, ensemble_matrix=ensembles, stop_reasons=stops, &
            mad_ambient=mad, ierr=ierr)
        call assert_equal_int(get_err_arg_pos(ierr), 6_int32, "direct ensemble growth names density_labels argument")

        call grow_ensemble(vectors, n_dims, n_vecs, mask, 1.0_real64, dims, kd, labels, &
            0.5_real64, grown, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "generated growth clears implementation argument positions")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dims, kd, labels, seeds, n_seeds, &
            1.0_real64, t_observables=2_int32, n_tiles=1_int32, ensemble_matrix=ensembles, &
            stop_reasons=stops, mad_ambient=mad, ierr=ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "generated ensembles clear implementation argument positions")
    end subroutine test_zero_mad_error_positions

    !> Acceptance compares two consecutive steps, so every entry point requires a history of two columns.
    subroutine test_history_requires_two_columns()
        integer(int32), parameter :: n_dims = 2_int32, n_vecs = 4_int32, n_seeds = 2_int32
        real(real64) :: density_labels(n_vecs), single(5, 1), vectors(n_dims, n_vecs), tmp_val_buf(n_vecs)
        real(real64) :: mad_ambient
        logical(c_bool) :: ensemble_mask(n_vecs), is_accepted, ensemble_matrix(n_vecs, n_seeds)
        integer(int32) :: ierr, n_ensembles, seed_indices(n_seeds), stop_reasons(n_seeds)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs)

        density_labels = [6.0_real64, 6.0_real64, 4.0_real64, 4.0_real64]
        ensemble_mask = [.true., .false., .false., .false.]
        single = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 1_int32, single, &
            size(single, 2, kind=int32), size(single, 1, kind=int32), ierr)
        call assert_true(ierr /= ERR_OK, "a one-column history is rejected by the observable update")

        call accept_ensemble(single, size(single, 2, kind=int32), size(single, 1, kind=int32), &
            2_int32, 0.5_real64, is_accepted, ierr)
        call assert_true(ierr /= ERR_OK, "a one-column history is rejected by acceptance")

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]
        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 3_int32]
        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "two columns: tree construction check")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=1_int32, n_tiles=1_int32, &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_true(ierr /= ERR_OK, "a one-iteration history window is rejected by growth")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=2_int32, n_tiles=1_int32, &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, "a two-iteration history window is accepted")
    end subroutine test_history_requires_two_columns

    subroutine test_compute_ensemble_observable_empty_and_zero()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 5)
        logical(c_bool) :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        observables = 0.0_real64

        density_labels = [10.0_real64, 5.0_real64, 2.0_real64]
        ensemble_mask = [.false., .false., .false.]

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 1_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)
        call assert_equal_int(ierr, ERR_OK, "empty ensemble ierr check")
        call assert_equal_real(observables(1, 1), 0.0_real64, 1.0e-10_real64, "empty ensemble outputs zeros")

        density_labels = [0.0_real64, 0.0_real64, 10.0_real64]
        ensemble_mask = [.true., .true., .false.]

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 2_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)
        call assert_true(ierr /= ERR_OK, "zero density labels must be rejected")

        density_labels = [1.0_real64, 1.0_real64, 10.0_real64]
        ensemble_mask = [.true., .true., .false.]

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 2_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)
        call assert_equal_int(ierr, ERR_OK, "unit density labels are accepted")
        call assert_equal_real(observables(1, 2), 1.0_real64, 1.0e-10_real64, "rho_arith is 1.0")
        call assert_equal_real(observables(2, 2), 1.0_real64, 1.0e-10_real64, "rho_harm is 1.0")
        call assert_equal_real(observables(3, 2), 1.0_real64, 1.0e-10_real64, &
                               "H_rho is 1.0 for a genuinely homogeneous ensemble")
    end subroutine test_compute_ensemble_observable_empty_and_zero

    subroutine test_compute_ensemble_observable_plain()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs)
        real(real64) :: observables(5, 5)
        logical(c_bool) :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64]
        ensemble_mask = [.true., .true., .false.]

        observables = 0.0_real64
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 1_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)

        call assert_equal_int(ierr, ERR_OK, "observable history update succeeds")
        call assert_equal_real(observables(1, 1), 10.0_real64, 1.0e-10_real64, "arithmetic mean check")
    end subroutine test_compute_ensemble_observable_plain

    subroutine test_compute_ensemble_observable_invalid_inputs()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 5)
        logical(c_bool) :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64]
        ensemble_mask = [.true., .false., .false.]
        observables = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, -1_int32, 1_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)
        call assert_true(ierr /= ERR_OK, "negative candidate_count must fail")

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 0_int32, observables, &
            size(observables, 2, kind=int32), size(observables, 1, kind=int32), ierr)
        call assert_true(ierr /= ERR_OK, "current_iter <= 0 must fail")
    end subroutine test_compute_ensemble_observable_invalid_inputs

    subroutine test_accept_ensemble_default()
        real(real64) :: observables(5, 5)
        logical(c_bool) :: is_accepted
        integer(int32) :: ierr

        observables = 0.0_real64
        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.2_real64

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            2_int32, is_accepted=is_accepted, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "accept_ensemble default ierr control")
        call assert_true(logical(is_accepted), "log2 fold change 0.263 < default 0.5 should be accepted")

        observables(3, 3) = 1.8_real64
        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            3_int32, is_accepted=is_accepted, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "accept_ensemble rejection ierr control")
        call assert_false(logical(is_accepted), "log2 fold change 0.585 >= default 0.5 should be rejected")
    end subroutine test_accept_ensemble_default

    subroutine test_accept_ensemble_alpha_sensitivity()
        real(real64) :: observables(5, 5)
        logical(c_bool) :: is_accepted
        integer(int32) :: ierr

        observables = 0.0_real64
        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.5_real64

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            2_int32, 0.2_real64, is_accepted, ierr)
        call assert_false(logical(is_accepted), "alpha 0.2 rejects 0.585 change")

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            2_int32, 0.4_real64, is_accepted, ierr)
        call assert_false(logical(is_accepted), "alpha 0.4 rejects 0.585 change")

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            2_int32, 0.8_real64, is_accepted, ierr)
        call assert_true(logical(is_accepted), "alpha 0.8 accepts 0.585 change")

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            2_int32, 1.0_real64, is_accepted, ierr)
        call assert_true(logical(is_accepted), "alpha 1.0 accepts 0.585 change")

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            2_int32, 1.5_real64, is_accepted, ierr)
        call assert_true(logical(is_accepted), "alpha 1.5 accepts 0.585 change")
    end subroutine test_accept_ensemble_alpha_sensitivity

    subroutine test_accept_ensemble_boundary_and_window()
        real(real64) :: observables(5, 2)
        logical(c_bool) :: is_accepted
        integer(int32) :: ierr

        observables = 0.0_real64

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            1_int32, 0.5_real64, is_accepted, ierr)
        call assert_true(logical(is_accepted), "iter 1 must always be accepted")

        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.1_real64
        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            3_int32, 0.5_real64, is_accepted, ierr)
        call assert_true(logical(is_accepted), "sliding window check accepts smooth step")
    end subroutine test_accept_ensemble_boundary_and_window

    subroutine test_accept_ensemble_invalid_inputs()
        real(real64) :: observables(5, 5)
        logical(c_bool) :: is_accepted
        integer(int32) :: ierr

        observables = 1.0_real64

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            0_int32, is_accepted=is_accepted, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "current_iter <= 0 must fail")

        call accept_ensemble(observables, size(observables, 2, kind=int32), size(observables, 1, kind=int32), &
            2_int32, -0.1_real64, is_accepted, ierr)
        call assert_true(ierr /= ERR_OK, "negative alpha_accept must fail")
    end subroutine test_accept_ensemble_invalid_inputs

    subroutine test_accept_ensemble_helper_basic()

        real(real64) :: observables(5, 4)
        logical(c_bool) :: is_accepted

        observables = 0.0_real64

        call accept_ensemble_helper(observables, 1_int32, alpha_accept=0.5_real64, &
                                    is_accepted=is_accepted)
        call assert_true(logical(is_accepted), "helper iter 1 must always accept")

        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.2_real64

        call accept_ensemble_helper(observables, 2_int32, alpha_accept=0.5_real64, &
                                    is_accepted=is_accepted)
        call assert_true(logical(is_accepted), "helper accepts smooth fold change")

        observables(3, 3) = 2.0_real64

        call accept_ensemble_helper(observables, 3_int32, alpha_accept=0.5_real64, &
                                    is_accepted=is_accepted)
        call assert_false(logical(is_accepted), "helper rejects sharp heterogeneity jump")

    end subroutine test_accept_ensemble_helper_basic

    subroutine test_obtain_ensembles_basic()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: n_seeds = 2_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_raw, n_merged
        integer(int32) :: seed_indices(n_seeds)
        logical(c_bool), allocatable :: raw_matrix(:, :)
        logical(c_bool) :: merged_matrix(n_vecs, n_seeds)
        real(real64) :: mad_ambient
        integer(int32) :: stop_reasons(n_seeds)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]
        vectors(:, 5) = [100.0_real64, 100.0_real64]
        vectors(:, 6) = [200.0_real64, 200.0_real64]

        density_labels = [10.0_real64, 10.0_real64, 8.0_real64, 8.0_real64, &
                          1.0_real64, 1.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        seed_indices = [1_int32, 3_int32]

        allocate (raw_matrix(n_vecs, n_seeds))

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=raw_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, &
            ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(raw_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "obtain_ensembles execution check")

        call merge_ensembles(raw_matrix, n_vecs, n_raw, 0.0_real64, merged_matrix, n_merged, ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles execution check")
        call assert_equal_int(n_merged, 2_int32, "expected 2 distinct merged clusters")

        call assert_false(logical(any(merged_matrix(5, :))), "vector 5 is background noise")
        call assert_false(logical(any(merged_matrix(6, :))), "vector 6 is background noise")

    end subroutine test_obtain_ensembles_basic

    subroutine test_obtain_ensembles_unmerged()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: n_seeds = 4_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs)
        integer(int32) :: seed_indices(n_seeds), ierr, n_ensembles
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds)
        real(real64) :: mad_ambient
        integer(int32) :: stop_reasons(n_seeds)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.1_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.1_real64]

        density_labels = [6.0_real64, 6.0_real64, 4.0_real64, 4.0_real64]
        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 2_int32, 3_int32, 4_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "unmerged obtain_ensembles execution check")
        call assert_equal_int(n_ensembles, n_seeds, &
                              "unmerged run keeps all grown seed ensembles")
        call assert_equal_int(size(ensemble_matrix, 1), n_vecs, &
                              "row dimension matches n_vectors")
        call assert_equal_int(size(ensemble_matrix, 2), n_seeds, &
                              "column dimension matches n_seeds")

    end subroutine test_obtain_ensembles_unmerged

    subroutine test_obtain_ensembles_zero_seeds()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32
        integer(int32), parameter :: n_seeds = 0_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: seed_indices(1), ierr, n_ensembles
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds)
        real(real64) :: mad_ambient
        integer(int32) :: stop_reasons(n_seeds)

        vectors = 0.0_real64
        density_labels = 1.0_real64
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]
        seed_indices = [1_int32]

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "zero seeds execution check")
        call assert_equal_int(n_ensembles, 0_int32, "n_ensembles must equal 0")
        call assert_equal_int(size(ensemble_matrix, 1), n_vecs, &
                              "row count must equal n_vectors")
        call assert_equal_int(size(ensemble_matrix, 2), 0_int32, &
                              "column count must be 0")

    end subroutine test_obtain_ensembles_zero_seeds

    subroutine test_obtain_ensembles_invalid_inputs()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32
        integer(int32), parameter :: n_seeds = 1_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: seed_indices(n_seeds), ierr, n_ensembles
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds)
        real(real64) :: mad_ambient
        integer(int32) :: stop_reasons(n_seeds)

        vectors = 0.0_real64
        density_labels = 1.0_real64
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]
        seed_indices = [1_int32]

        seed_indices(1) = 4_int32

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_true(ierr /= ERR_OK, "seed_indices out of bounds must fail")

        seed_indices(1) = 1_int32

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, -0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_true(ierr /= ERR_OK, "negative r must fail")

    end subroutine test_obtain_ensembles_invalid_inputs

    subroutine test_obtain_ensembles_500_vectors()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 500_int32
        integer(int32), parameter :: n_seeds = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        real(real64) :: angle, radius
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs)
        integer(int32) :: ierr, n_raw, n_merged, i_vec
        integer(int32) :: seed_indices(n_seeds)
        logical(c_bool), allocatable :: raw_matrix(:, :)
        logical(c_bool) :: merged_matrix(n_vecs, n_seeds)
        character(len=128) :: assert_msg
        real(real64) :: mad_ambient
        integer(int32) :: stop_reasons(n_seeds)

        do i_vec = 1, 150
            angle = real(i_vec, real64)*0.15_real64
            radius = real(mod(i_vec, 15_int32) + 1_int32, real64)*0.1_real64
            vectors(1, i_vec) = radius*cos(angle)
            vectors(2, i_vec) = radius*sin(angle)
        end do

        do i_vec = 151, 300
            angle = real(i_vec, real64)*0.15_real64
            radius = real(mod(i_vec, 15_int32) + 1_int32, real64)*0.1_real64
            vectors(1, i_vec) = 50.0_real64 + radius*cos(angle)
            vectors(2, i_vec) = 50.0_real64 + radius*sin(angle)
        end do

        do i_vec = 301, 450
            angle = real(i_vec, real64)*0.15_real64
            radius = real(mod(i_vec, 15_int32) + 1_int32, real64)*0.1_real64
            vectors(1, i_vec) = -50.0_real64 + radius*cos(angle)
            vectors(2, i_vec) = 50.0_real64 + radius*sin(angle)
        end do

        do i_vec = 451, 500
            vectors(1, i_vec) = 200.0_real64 + &
                                real(i_vec - 450_int32, real64)*10.0_real64
            vectors(2, i_vec) = 200.0_real64 + &
                                real(i_vec - 450_int32, real64)*10.0_real64
        end do

        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, dimension_order, kd_indices, &
            min(8_int32, n_vecs), density_labels, ierr)

        seed_indices = [15_int32, 165_int32, 315_int32]

        allocate (raw_matrix(n_vecs, n_seeds))

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=raw_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, &
            ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(raw_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "obtain_ensembles 500 execution check")

        call merge_ensembles(raw_matrix, n_vecs, n_raw, 0.0_real64, merged_matrix, n_merged, ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles 500 execution check")
        call assert_equal_int(n_merged, 3_int32, &
                              "expected 3 distinct merged clusters")

        do i_vec = 451, 500
            write (assert_msg, '(A,I0)') &
                "test_obtain_ensembles_500_vectors: background noise verification at index: ", &
                i_vec
            call assert_false(logical(any(merged_matrix(i_vec, :))), assert_msg)
        end do

    end subroutine test_obtain_ensembles_500_vectors

    subroutine test_merge_ensembles_filters_final_singletons()
        integer(int32), parameter :: n_vectors = 5_int32
        integer(int32), parameter :: n_seeds = 4_int32

        logical(c_bool) :: raw_masks(n_vectors, n_seeds)
        logical(c_bool) :: merged_matrix(n_vectors, n_seeds)
        integer(int32) :: n_ensembles
        integer(int32) :: ierr

        raw_masks = .false.

        raw_masks(1, 1) = .true.

        raw_masks(1, 2) = .true.
        raw_masks(2, 2) = .true.

        raw_masks(3, 3) = .true.

        raw_masks(4, 4) = .true.
        raw_masks(5, 4) = .true.

        call merge_ensembles(raw_masks, n_vectors, n_seeds, 0.0_real64, merged_matrix, n_ensembles, ierr)

        call assert_equal_int(ierr, ERR_OK, &
                              "singleton filtering: merge execution")
        call assert_equal_int(n_ensembles, 2_int32, &
                              "singleton filtering: isolated singleton removed")
        call assert_true(logical(any(merged_matrix(1, :))), &
                         "singleton filtering: mergeable singleton member retained")
        call assert_true(logical(any(merged_matrix(2, :))), &
                         "singleton filtering: merged partner retained")
        call assert_false(logical(any(merged_matrix(3, :))), &
                          "singleton filtering: isolated singleton omitted")
        call assert_true(logical(any(merged_matrix(4, :))), &
                         "singleton filtering: non-singleton cluster retained")
        call assert_true(logical(any(merged_matrix(5, :))), &
                         "singleton filtering: second non-singleton member retained")
    end subroutine test_merge_ensembles_filters_final_singletons

    subroutine test_merge_ensembles_basic()
        integer(int32), parameter :: n_vecs = 5_int32
        integer(int32), parameter :: n_seeds = 3_int32

        logical(c_bool) :: raw_masks(n_vecs, n_seeds)
        logical(c_bool) :: merged_matrix(n_vecs, n_seeds)
        integer(int32) :: n_ensembles, ierr

        raw_masks(:, 1) = [.true., .true., .false., .false., .false.]
        raw_masks(:, 2) = [.false., .true., .true., .false., .false.]
        raw_masks(:, 3) = [.false., .false., .false., .true., .true.]

        call merge_ensembles(raw_masks, n_vecs, n_seeds, 0.0_real64, merged_matrix, n_ensembles, ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles_basic ierr check")
        call assert_equal_int(n_ensembles, 2_int32, "masks 1 & 2 merge into 1, leaving 2 unique clusters")

        call assert_true(logical(merged_matrix(1, 1)), "merged 1 contains elem 1")
        call assert_true(logical(merged_matrix(2, 1)), "merged 1 contains elem 2")
        call assert_true(logical(merged_matrix(3, 1)), "merged 1 contains elem 3")
        call assert_false(logical(merged_matrix(4, 1)), "merged 1 excludes elem 4")
    end subroutine test_merge_ensembles_basic

    subroutine test_merge_ensembles_zero_seeds()
        integer(int32), parameter :: n_vecs = 4_int32

        logical(c_bool) :: raw_masks(n_vecs, 0)
        logical(c_bool) :: merged_matrix(n_vecs, 0)
        integer(int32) :: n_ensembles, ierr

        call merge_ensembles(raw_masks, n_vecs, 0_int32, 0.0_real64, merged_matrix, n_ensembles, ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles_zero_seeds ierr check")
        call assert_equal_int(n_ensembles, 0_int32, "n_ensembles is 0")
        call assert_equal_int(size(merged_matrix, 2), 0_int32, "columns count is 0")
    end subroutine test_merge_ensembles_zero_seeds

    subroutine test_merge_ensembles_no_overlap()
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: n_seeds = 2_int32

        logical(c_bool) :: raw_masks(n_vecs, n_seeds)
        logical(c_bool) :: merged_matrix(n_vecs, n_seeds)
        integer(int32) :: n_ensembles, ierr

        raw_masks(:, 1) = [.true., .true., .false., .false.]
        raw_masks(:, 2) = [.false., .false., .true., .true.]

        call merge_ensembles(raw_masks, n_vecs, n_seeds, 0.0_real64, merged_matrix, n_ensembles, ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles_no_overlap ierr check")
        call assert_equal_int(n_ensembles, 2_int32, "disjoint masks do not merge")
    end subroutine test_merge_ensembles_no_overlap

    subroutine test_merge_ensembles_invalid_inputs()
        integer(int32), parameter :: n_vecs = 3_int32
        integer(int32), parameter :: n_seeds = 2_int32

        logical(c_bool) :: raw_masks(n_vecs, n_seeds)
        logical(c_bool) :: merged_matrix(n_vecs, n_seeds)
        integer(int32) :: n_ensembles, ierr

        raw_masks = .true.

        call merge_ensembles(raw_masks, n_vecs, n_seeds, -0.1_real64, merged_matrix, n_ensembles, ierr)
        call assert_true(ierr /= ERR_OK, "min_overlap_coefficient < 0 must fail")

        call merge_ensembles(raw_masks, n_vecs, n_seeds, 1.1_real64, merged_matrix, n_ensembles, ierr)
        call assert_true(ierr /= ERR_OK, "min_overlap_coefficient > 1 must fail")
    end subroutine test_merge_ensembles_invalid_inputs

    subroutine test_merge_ensembles_overlap_coefficient_regimes()
        integer(int32), parameter :: n_vecs = 5_int32
        integer(int32), parameter :: n_seeds = 3_int32
        integer(int32), parameter :: n_chain = 4_int32

        logical(c_bool) :: raw_masks(n_vecs, n_seeds)
        logical(c_bool) :: chain_masks(n_vecs, n_chain)
        logical(c_bool) :: merged_matrix(n_vecs, n_chain)
        integer(int32) :: n_ensembles, ierr

        raw_masks = .false.
        raw_masks(1, 1) = .true.; raw_masks(2, 1) = .true.; raw_masks(3, 1) = .true.
        raw_masks(1, 2) = .true.; raw_masks(2, 2) = .true.; raw_masks(4, 2) = .true.
        raw_masks(3, 3) = .true.; raw_masks(4, 3) = .true.; raw_masks(5, 3) = .true.

        call merge_ensembles(raw_masks, n_vecs, n_seeds, 0.0_real64, merged_matrix, n_ensembles, ierr)
        call assert_equal_int(ierr, ERR_OK, "transitivity regimes: threshold 1 execution check")
        call assert_equal_int(n_ensembles, 1_int32, "threshold 1 merges all three ensembles")
        call assert_equal_int(count(merged_matrix(:, 1)), n_vecs, &
                              "threshold 1 union covers every vector")

        ! OC(E1,E2) = 2/3, OC(E1,E3) = OC(E2,E3) = 1/3. A threshold between the two keeps E3 out, and --
        ! unlike an absolute count on growing masks -- the E1/E2 union cannot pull it back in afterwards.
        call merge_ensembles(raw_masks, n_vecs, n_seeds, 0.6_real64, merged_matrix, n_ensembles, ierr)
        call assert_equal_int(ierr, ERR_OK, "overlap coefficient regimes: threshold 0.6 execution check")
        call assert_equal_int(n_ensembles, 2_int32, &
                              "threshold 0.6 merges E1 and E2 but leaves E3 separate")
        call assert_equal_int(count(merged_matrix(:, 1)), 4_int32, &
                              "threshold 0.6 merges E1 and E2 into four vectors")
        call assert_equal_int(count(merged_matrix(:, 2)), 3_int32, &
                              "threshold 0.6 leaves E3 at its own three vectors")

        ! An Overlap Coefficient of 1 demands full containment of the smaller ensemble in the larger.
        call merge_ensembles(raw_masks, n_vecs, n_seeds, 1.0_real64, merged_matrix, n_ensembles, ierr)
        call assert_equal_int(ierr, ERR_OK, "overlap coefficient regimes: threshold 1.0 execution check")
        call assert_equal_int(n_ensembles, 3_int32, &
                              "threshold 1.0 merges nothing here, no ensemble contains another")

        chain_masks = .false.
        chain_masks(1, 1) = .true.; chain_masks(2, 1) = .true.
        chain_masks(2, 2) = .true.; chain_masks(3, 2) = .true.
        chain_masks(3, 3) = .true.; chain_masks(4, 3) = .true.
        chain_masks(4, 4) = .true.; chain_masks(5, 4) = .true.

        call merge_ensembles(chain_masks, n_vecs, n_chain, 0.0_real64, merged_matrix, n_ensembles, ierr)
        call assert_equal_int(ierr, ERR_OK, "transitivity regimes: chain execution check")
        call assert_equal_int(n_ensembles, 1_int32, "a chain of overlaps collapses to one ensemble")
        call assert_equal_int(count(merged_matrix(:, 1)), n_vecs, &
                              "the chain union covers every vector")
    end subroutine test_merge_ensembles_overlap_coefficient_regimes

    subroutine test_merge_ensembles_clears_unused_columns()
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: n_seeds = 4_int32

        logical(c_bool) :: raw_masks(n_vecs, n_seeds)
        logical(c_bool) :: merged_masks(n_vecs, n_seeds)
        logical(c_bool) :: tmp_active_flag(n_seeds)
        integer(int32) :: tmp_parent(n_seeds)
        integer(int32) :: n_ensembles, ierr, i_col

        raw_masks = .false.
        raw_masks(1, 1) = .true.; raw_masks(2, 1) = .true.
        raw_masks(2, 2) = .true.; raw_masks(3, 2) = .true.
        raw_masks(4, 3) = .true.; raw_masks(5, 3) = .true.
        raw_masks(5, 4) = .true.; raw_masks(6, 4) = .true.

        merged_masks = .false.
        tmp_active_flag = .false.
        tmp_parent = 0_int32

        call merge_ensembles_expert(raw_masks, n_vecs, n_seeds, 0.0_real64, merged_masks, tmp_active_flag, &
            tmp_parent, n_ensembles, ierr)

        call assert_equal_int(ierr, ERR_OK, "clear unused columns: execution check")
        call assert_equal_int(n_ensembles, 2_int32, "two ensembles survive the merge")

        call assert_equal_int(count(merged_masks(:, 1)), 3_int32, "first ensemble holds 3 vectors")
        call assert_equal_int(count(merged_masks(:, 2)), 3_int32, "second ensemble holds 3 vectors")

        do i_col = n_ensembles + 1_int32, n_seeds
            call assert_equal_int(count(merged_masks(:, i_col)), 0_int32, &
                                  "columns past n_ensembles are cleared")
        end do
    end subroutine test_merge_ensembles_clears_unused_columns

    subroutine test_obtain_ensembles_tile_count_invariance()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 12_int32
        integer(int32), parameter :: n_seeds = 5_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_ensembles
        integer(int32) :: seed_indices(n_seeds), i_vec, i_tiles
        logical(c_bool) :: reference_matrix(n_vecs, n_seeds), tiled_matrix(n_vecs, n_seeds)
        real(real64) :: mad_ambient
        integer(int32) :: stop_reasons(n_seeds)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [0.0_real64, 0.4_real64]
        vectors(:, 4) = [10.0_real64, 10.0_real64]
        vectors(:, 5) = [10.0_real64, 10.2_real64]
        vectors(:, 6) = [10.0_real64, 10.4_real64]
        vectors(:, 7) = [20.0_real64, 20.0_real64]
        vectors(:, 8) = [20.0_real64, 20.2_real64]
        vectors(:, 9) = [30.0_real64, 30.0_real64]
        vectors(:, 10) = [30.0_real64, 30.2_real64]
        vectors(:, 11) = [100.0_real64, 100.0_real64]
        vectors(:, 12) = [200.0_real64, 200.0_real64]

        density_labels = [8.0_real64, 8.0_real64, 8.0_real64, 8.0_real64, &
                          8.0_real64, 8.0_real64, 7.0_real64, 7.0_real64, &
                          7.0_real64, 7.0_real64, 1.0_real64, 1.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "tile invariance: tree construction check")

        seed_indices = [1_int32, 4_int32, 7_int32, 9_int32, 11_int32]

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=1_int32, &
            ensemble_matrix=reference_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(reference_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "tile invariance: serial reference run check")
        call assert_equal_int(n_ensembles, n_seeds, "tile invariance: reference ensemble count")

        do i_tiles = 2, 7
            tiled_matrix = .false.

            call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
                seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=i_tiles, &
                ensemble_matrix=tiled_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
            if (ierr == ERR_OK) n_ensembles = size(tiled_matrix, dim=2, kind=int32)

            call assert_equal_int(ierr, ERR_OK, "tile invariance: tiled run execution check")
            call assert_equal_int(n_ensembles, n_seeds, "tile invariance: tiled ensemble count")

            do i_vec = 1, n_vecs
                call assert_true(logical(all(tiled_matrix(i_vec, :) .eqv. reference_matrix(i_vec, :))), &
                                 "tile invariance: tiled result must match the serial reference")
            end do
        end do

        call assert_true(all(count(reference_matrix, dim=1) >= 1_int32), &
                         "tile invariance: every seed produces a non-empty ensemble")

    end subroutine test_obtain_ensembles_tile_count_invariance

    subroutine test_obtain_ensembles_invalid_n_tiles()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: n_seeds = 2_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_ensembles
        integer(int32) :: seed_indices(n_seeds)
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds)
        real(real64) :: mad_ambient
        integer(int32) :: stop_reasons(n_seeds)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]

        density_labels = [6.0_real64, 6.0_real64, 4.0_real64, 4.0_real64]
        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 3_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=0_int32, &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_true(ierr /= ERR_OK, "n_tiles = 0 must fail")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=-3_int32, &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_true(ierr /= ERR_OK, "negative n_tiles must fail")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=64_int32, &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, "n_tiles above n_seeds must be accepted, not rejected")
        call assert_equal_int(n_ensembles, n_seeds, "surplus-tile run still grows every seed")

    end subroutine test_obtain_ensembles_invalid_n_tiles

    subroutine test_obtain_ensembles_layer_n_tiles_agreement()

        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: n_seeds = 2_int32
        integer(int32), parameter :: n_tiles = 5_int32
        integer(int32), parameter :: t_obs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_ensembles, n_ensembles_alloc
        integer(int32) :: seed_indices(n_seeds), i_vec

        integer(int32) :: tmp_stack(3, 64, n_tiles), tmp_perm(n_vecs, n_tiles)
        logical(c_bool) :: tmp_vicinity_mask(n_vecs, n_tiles), tmp_surface_mask(n_vecs, n_tiles)
        real(real64) :: tmp_abs_diff(n_vecs, n_tiles)
        real(real64) :: tmp_observables(5, min(t_obs, n_vecs), n_tiles)
        logical(c_bool) :: tmp_current_mask(n_vecs, n_tiles)
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds), alloc_matrix(n_vecs, n_seeds)
        integer(int32) :: stop_reasons(n_seeds)
        real(real64) :: mad_ambient

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]
        vectors(:, 5) = [100.0_real64, 100.0_real64]
        vectors(:, 6) = [200.0_real64, 200.0_real64]

        density_labels = [10.0_real64, 10.0_real64, 8.0_real64, 8.0_real64, &
                          1.0_real64, 1.0_real64]
        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 3_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "layer agreement: tree construction check")

        call obtain_ensembles_seed_radius_expert(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
            density_labels, seed_indices, n_seeds, spread(0.5_real64, 1_int32, n_seeds), 0.5_real64, 0.5_real64, &
            t_obs, n_tiles, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, tmp_observables, &
            tmp_current_mask, ensemble_matrix, stop_reasons, mad_ambient, ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, &
                              "layer agreement: preallocated layer must accept n_tiles > n_seeds")
        call assert_equal_int(n_ensembles, n_seeds, "layer agreement: every seed grown")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, 0.5_real64, 0.5_real64, t_obs, n_tiles, alloc_matrix, stop_reasons, &
            mad_ambient, ierr)
        if (ierr == ERR_OK) n_ensembles_alloc = size(alloc_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, &
                              "layer agreement: allocating layer must accept n_tiles > n_seeds")
        call assert_equal_int(n_ensembles_alloc, n_ensembles, &
                              "layer agreement: both layers report the same ensemble count")

        do i_vec = 1, n_vecs
            call assert_true(logical(all(alloc_matrix(i_vec, :) .eqv. ensemble_matrix(i_vec, :))), &
                             "layer agreement: both layers produce identical ensembles")
        end do

    end subroutine test_obtain_ensembles_layer_n_tiles_agreement

    subroutine test_grow_ensemble_rejects_negative_density_labels()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical(c_bool) :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs), tmp_surface_mask(n_vecs)
        logical(c_bool) :: grown_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.2_real64]
        vectors(:, 4) = [0.0_real64, 0.4_real64]

        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "negative density labels: tree construction check")

        density_labels = [5.0_real64, 5.0_real64, -1.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false., .false.]

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(ierr /= ERR_OK, &
                         "grow_ensemble must reject negative density labels")

        ensemble_mask = [.true., .false., .false., .false.]

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(ierr /= ERR_OK, &
                         "grow_ensemble_alloc must reject negative density labels")

        density_labels = [5.0_real64, 5.0_real64, 0.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false., .false.]

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(ierr /= ERR_OK, &
                         "grow_ensemble must reject a zero density label")

        ensemble_mask = [.true., .false., .false., .false.]

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_true(ierr /= ERR_OK, &
                         "grow_ensemble_alloc must reject a zero density label")

        density_labels = [5.0_real64, 8.0_real64, 1.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false., .false.]

        call grow_ensemble_expert(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, dimension_order, kd_indices, &
            density_labels, 0.5_real64, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            grown_mask, ierr)
        if (ierr == ERR_OK) ensemble_mask = grown_mask
        call assert_equal_int(ierr, ERR_OK, &
                              "grow_ensemble must accept a unit density label")
    end subroutine test_grow_ensemble_rejects_negative_density_labels

    subroutine test_observable_history_matrix()
        integer(int32), parameter :: n_vecs = 3_int32

        logical(c_bool) :: ensemble_mask(n_vecs)
        real(real64) :: density_labels(n_vecs)
        real(real64) :: too_few_rows(4, 6), correct_rows(5, 6)
        logical(c_bool) :: is_accepted
        integer(int32) :: ierr

        ensemble_mask = .true._c_bool
        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        too_few_rows = 0.0_real64
        correct_rows = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 1_int32, &
                                         too_few_rows, 6_int32, 4_int32, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "observable update rejects a four-row matrix")

        call accept_ensemble(too_few_rows, 6_int32, 4_int32, 2_int32, is_accepted=is_accepted, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "ensemble acceptance rejects a four-row matrix")

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, 0_int32, 1_int32, correct_rows, &
            size(correct_rows, 2, kind=int32), size(correct_rows, 1, kind=int32), ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "compute_ensemble_observable must accept a 5-row matrix")

        call accept_ensemble(correct_rows, size(correct_rows, 2, kind=int32), size(correct_rows, 1, kind=int32), &
            2_int32, 0.5_real64, is_accepted, ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "accept_ensemble must accept a 5-row matrix")
    end subroutine test_observable_history_matrix

    subroutine test_obtain_ensembles_rejects_invalid_history()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: n_seeds = 2_int32
        integer(int32), parameter :: n_tiles = 2_int32
        integer(int32), parameter :: t_obs = 6_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_ensembles
        integer(int32) :: seed_indices(n_seeds)

        integer(int32) :: tmp_stack(3, 64, n_tiles), tmp_perm(n_vecs, n_tiles)
        logical(c_bool) :: tmp_vicinity_mask(n_vecs, n_tiles), tmp_surface_mask(n_vecs, n_tiles)
        real(real64) :: tmp_abs_diff(n_vecs, n_tiles)
        real(real64) :: good_observables(5, min(t_obs, n_vecs), n_tiles)
        logical(c_bool) :: tmp_current_mask(n_vecs, n_tiles)
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds)
        integer(int32) :: stop_reasons(n_seeds)
        real(real64) :: mad_ambient

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]

        density_labels = [6.0_real64, 6.0_real64, 4.0_real64, 4.0_real64]
        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 3_int32]
        good_observables = 0.0_real64

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "observable rows: tree construction check")

        call obtain_ensembles_seed_radius_expert(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
            density_labels, seed_indices, n_seeds, spread(0.5_real64, 1_int32, n_seeds), 0.5_real64, 0.5_real64, &
            0_int32, n_tiles, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            good_observables, tmp_current_mask, ensemble_matrix, stop_reasons, mad_ambient, ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_true(ierr /= ERR_OK, &
                         "obtain_ensembles must reject a nonpositive history window")

        call obtain_ensembles_seed_radius_expert(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
            density_labels, seed_indices, n_seeds, spread(0.5_real64, 1_int32, n_seeds), 0.5_real64, 0.5_real64, &
            t_obs, n_tiles, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            good_observables, tmp_current_mask, ensemble_matrix, stop_reasons, mad_ambient, ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, &
                              "obtain_ensembles must accept a correctly sized observable workspace")
        call assert_equal_int(n_ensembles, n_seeds, "every seed still grown")
    end subroutine test_obtain_ensembles_rejects_invalid_history

    subroutine test_obtain_ensembles_stop_reasons()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 8_int32
        integer(int32), parameter :: n_seeds = 2_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_raw
        integer(int32) :: seed_indices(n_seeds), stop_reasons(n_seeds)
        logical(c_bool) :: raw_matrix(n_vecs, n_seeds)
        real(real64) :: mad_ambient

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.4_real64]
        vectors(:, 3) = [0.0_real64, 0.8_real64]
        vectors(:, 4) = [0.0_real64, 1.2_real64]
        vectors(:, 5) = [50.0_real64, 50.0_real64]
        vectors(:, 6) = [60.0_real64, 60.0_real64]
        vectors(:, 7) = [70.0_real64, 70.0_real64]
        vectors(:, 8) = [80.0_real64, 80.0_real64]

        density_labels = [4.0_real64, 4.0_real64, 4.0_real64, 4.0_real64, &
                          1.0_real64, 2.0_real64, 3.0_real64, 5.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "stop reasons: tree construction check")

        seed_indices = [1_int32, 5_int32]
        stop_reasons = 0_int32

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=raw_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, &
            ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(raw_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "stop reasons: execution check")

        call assert_equal_int(stop_reasons(1), STC_STOP_FIXED_POINT, &
                              "chain seed reaches a fixed point after growing")
        call assert_true(count(raw_matrix(:, 1)) > 1_int32, &
                         "fixed-point ensemble grew beyond its seed")

        call assert_equal_int(stop_reasons(2), STC_STOP_NO_CANDIDATES, &
                              "isolated seed reports no candidates ever found")
        call assert_equal_int(count(raw_matrix(:, 2)), 1_int32, &
                              "isolated ensemble is still the bare seed")

        call assert_true(STC_STOP_NO_CANDIDATES /= STC_STOP_NEVER_ACCEPTED, &
                         "isolated seed and first-batch rejection are distinct codes")

        stop_reasons = 0_int32

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, alpha_accept=0.0_real64, t_observables=10_int32, &
            n_tiles=max(1_int32, min(8_int32, n_seeds)), ensemble_matrix=raw_matrix, stop_reasons=stop_reasons, &
            mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(raw_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "stop reasons: zero alpha_accept execution check")
        call assert_equal_int(stop_reasons(1), STC_STOP_NEVER_ACCEPTED, &
                              "zero alpha_accept rejects the first batch")
        call assert_equal_int(count(raw_matrix(:, 1)), 1_int32, &
                              "rejected trial leaves the accepted ensemble at the bare seed")

        call assert_equal_int(stop_reasons(2), STC_STOP_NO_CANDIDATES, &
                              "isolated seed still reports no candidates under zero alpha_accept")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=raw_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, &
            ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(raw_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, "stop reasons: omitting the argument still works")
    end subroutine test_obtain_ensembles_stop_reasons

    subroutine test_obtain_ensembles_observable_window_match()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: n_seeds = 2_int32
        integer(int32), parameter :: n_tiles = 2_int32
        integer(int32), parameter :: t_obs = 6_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_ensembles
        integer(int32) :: seed_indices(n_seeds), stop_reasons(n_seeds)

        integer(int32) :: tmp_stack(3, 64, n_tiles), tmp_perm(n_vecs, n_tiles)
        logical(c_bool) :: tmp_vicinity_mask(n_vecs, n_tiles), tmp_surface_mask(n_vecs, n_tiles)
        real(real64) :: tmp_abs_diff(n_vecs, n_tiles)
        real(real64) :: mad_ambient
        real(real64) :: exact_observables(5, min(t_obs, n_vecs), n_tiles)
        real(real64) :: narrow_observables(5, 2, n_tiles)
        logical(c_bool) :: tmp_current_mask(n_vecs, n_tiles)
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]

        density_labels = [6.0_real64, 6.0_real64, 4.0_real64, 4.0_real64]
        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 3_int32]
        exact_observables = 0.0_real64
        narrow_observables = 0.0_real64

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "observable window: tree construction check")

        call obtain_ensembles_seed_radius_expert(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
            density_labels, seed_indices, n_seeds, spread(0.5_real64, 1_int32, n_seeds), 0.5_real64, 0.5_real64, &
            t_obs, n_tiles, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            exact_observables, tmp_current_mask, ensemble_matrix, stop_reasons, mad_ambient, ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, &
                              "workspace matches the history capped by the input size")

        call obtain_ensembles_seed_radius_expert(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
            density_labels, seed_indices, n_seeds, spread(0.5_real64, 1_int32, n_seeds), 0.5_real64, 0.5_real64, &
            2_int32, n_tiles, tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, tmp_abs_diff, &
            narrow_observables, tmp_current_mask, ensemble_matrix, stop_reasons, mad_ambient, ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, &
                              "two-column observable history is accepted")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, &
                              "allocating wrapper still works when n_vectors < t_observables")
    end subroutine test_obtain_ensembles_observable_window_match

    subroutine test_obtain_ensembles_reports_degenerate_mad()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: n_seeds = 2_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, n_raw
        integer(int32) :: seed_indices(n_seeds)
        real(real64) :: mad_ambient, mad_quarter, mad_full
        logical(c_bool) :: raw_matrix(n_vecs, n_seeds), quarter_matrix(n_vecs, n_seeds)
        logical(c_bool) :: full_matrix(n_vecs, n_seeds)
        integer(int32) :: i_vec
        integer(int32) :: stop_reasons(n_seeds)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]
        vectors(:, 5) = [100.0_real64, 100.0_real64]
        vectors(:, 6) = [200.0_real64, 200.0_real64]

        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 3_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                                   tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "degenerate mad: tree construction check")

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64, 10.0_real64, &
                          1.0_real64, 1.0_real64]

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, 0.25_real64, t_observables=10_int32, &
            n_tiles=max(1_int32, min(8_int32, n_seeds)), ensemble_matrix=quarter_matrix, stop_reasons=stop_reasons, &
            mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(quarter_matrix, dim=2, kind=int32)
        call assert_true(ierr /= ERR_OK, &
                         "a degenerate ambient MAD must be rejected, not clustered under")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, 1.0_real64, t_observables=10_int32, &
            n_tiles=max(1_int32, min(8_int32, n_seeds)), ensemble_matrix=full_matrix, stop_reasons=stop_reasons, &
            mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(full_matrix, dim=2, kind=int32)
        call assert_true(ierr /= ERR_OK, &
                         "a degenerate ambient MAD is rejected for every alpha_mad")

        density_labels = [10.0_real64, 10.0_real64, 8.0_real64, 8.0_real64, &
                          1.0_real64, 1.0_real64]
        mad_ambient = -1.0_real64

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=raw_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, &
            ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(raw_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, "spread densities execution check")
        call assert_true(mad_ambient > 0.0_real64, &
                         "spread density labels give a non-zero ambient MAD")
        call assert_true(count(raw_matrix(:, 1)) > 1_int32, &
                         "growth proceeds once the threshold can separate")

        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.5_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=raw_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, &
            ierr=ierr)
        if (ierr == ERR_OK) n_raw = size(raw_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, "growth succeeds with the same explicit MAD output")
    end subroutine test_obtain_ensembles_reports_degenerate_mad

    subroutine test_compute_ambient_density_stats()

        integer(int32), parameter :: n_vecs = 5_int32

        real(real64) :: density_labels(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: tmp_perm(n_vecs)
        real(real64) :: median_ambient, mad_ambient

        density_labels = [8.0_real64, 2.0_real64, 4.0_real64, 10.0_real64, 1.0_real64]

        call compute_ambient_density_stats_helper(density_labels, n_vecs, tmp_perm, &
                                                  tmp_abs_diff, median_ambient, mad_ambient)

        call assert_equal_real(median_ambient, 4.0_real64, 1.0e-12_real64, &
                               "ambient median of odd sample")
        call assert_equal_real(mad_ambient, 3.0_real64, 1.0e-12_real64, &
                               "ambient MAD of odd sample")

        call assert_equal_real(density_labels(1), 8.0_real64, 1.0e-12_real64, &
                               "density labels stay unmodified")
        call assert_equal_real(density_labels(5), 1.0_real64, 1.0e-12_real64, &
                               "density labels stay unmodified")

    end subroutine test_compute_ambient_density_stats

    !> Every selected seed carries a positive locally fitted radius; non-seeds carry none.
    subroutine test_identify_ensemble_seeds_radii()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: k_seeding = 1_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        real(real64) :: tmp_distances(n_vecs), seed_radii(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_rec_stack(3, n_vecs)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds, ierr, i_vec
        integer(int32) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        logical(c_bool) :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)

        ! A dense triple around the origin and a sparse triple far away and widely spaced.
        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.1_real64, 0.0_real64]
        vectors(:, 3) = [0.2_real64, 0.0_real64]
        vectors(:, 4) = [50.0_real64, 0.0_real64]
        vectors(:, 5) = [55.0_real64, 0.0_real64]
        vectors(:, 6) = [60.0_real64, 0.0_real64]

        density_labels = [9.0_real64, 8.0_real64, 7.0_real64, &
                          6.0_real64, 5.0_real64, 4.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_rec_stack, ierr)

        call identify_ensemble_seeds_expert(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, tmp_perm, tmp_distances, tmp_stack, tmp_visited_mask, tmp_newly_covered_mask, sorted_perm, &
            n_seeds, seed_mask, seed_radii, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "seed radii: execution")

        do i_vec = 1, n_vecs
            if (seed_mask(i_vec)) then
                call assert_true(seed_radii(i_vec) > 0.0_real64, &
                                 "seed radii: every seed has a positive growth radius")
            else
                call assert_equal_real(seed_radii(i_vec), 0.0_real64, 0.0_real64, &
                                       "seed radii: non-seeds carry no radius")
            end if
        end do

        ! The dense region's seed must be fitted tighter than the sparse region's seed.
        call assert_true(logical(seed_mask(1)) .and. logical(seed_mask(4)), &
                         "seed radii: one seed per region")
        call assert_true(seed_radii(1) < seed_radii(4), &
                         "seed radii: dense region grows tighter than sparse region")

    end subroutine test_identify_ensemble_seeds_radii

    !> The radius is the requested percentile of the k nearest-neighbor distances.
    subroutine test_identify_ensemble_seeds_radii_percentile()
        integer(int32), parameter :: n_dims = 1_int32
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: k_seeding = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        real(real64) :: seed_radii(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_rec_stack(3, n_vecs)
        integer(int32) :: sorted_perm(n_vecs), n_seeds, ierr
        logical(c_bool) :: seed_mask(n_vecs)

        ! Neighbor distances of vector 1 are 1, 2 and 4.
        vectors(1, :) = [0.0_real64, 1.0_real64, 2.0_real64, 4.0_real64]
        density_labels = [9.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        dimension_order = [1_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_rec_stack, ierr)

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, sorted_perm, n_seeds, seed_mask, seed_radii, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "seed radii percentile: median execution")
        call assert_equal_real(seed_radii(1), 2.0_real64, 1.0e-12_real64, &
                               "seed radii percentile: default is the median neighbor distance")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, sorted_perm, n_seeds, seed_mask, seed_radii, 1.0_real64, ierr)
        call assert_equal_int(ierr, ERR_OK, "seed radii percentile: max execution")
        call assert_equal_real(seed_radii(1), 4.0_real64, 1.0e-12_real64, &
                               "seed radii percentile: p=1 is the farthest of the k neighbors")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, dimension_order, kd_indices, &
            k_seeding, sorted_perm, n_seeds, seed_mask, seed_radii, 1.5_real64, ierr)
        call assert_true(ierr /= ERR_OK, &
                         "seed radii percentile: a quantile outside [0,1] must fail")

    end subroutine test_identify_ensemble_seeds_radii_percentile

    !> Per-seed radii let one seed grow where a uniform radius large enough to do so
    !! would have merged the dense region into a single blob.
    subroutine test_obtain_ensembles_per_seed_radii()
        integer(int32), parameter :: n_dims = 1_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: n_seeds = 2_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        real(real64) :: per_seed_radii(n_seeds)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_rec_stack(3, n_vecs)
        integer(int32) :: seed_indices(n_seeds), ierr, n_ensembles
        logical(c_bool) :: ensemble_matrix(n_vecs, n_seeds), generated_matrix(n_vecs, n_seeds)
        integer(int32) :: stop_reasons(n_seeds), generated_stops(n_seeds), generated_count
        real(real64) :: mad_ambient, generated_mad

        ! Dense cluster 1..3 spaced 0.1 apart, sparse cluster 4..6 spaced 5 apart.
        vectors(1, :) = [0.0_real64, 0.1_real64, 0.2_real64, &
                         50.0_real64, 55.0_real64, 60.0_real64]
        density_labels = [8.0_real64, 8.0_real64, 8.0_real64, &
                          4.0_real64, 4.0_real64, 4.0_real64]
        dimension_order = [1_int32]
        seed_indices = [1_int32, 4_int32]

        call build_kd_index_expert(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_rec_stack, ierr)

        ! A uniform radius tight enough for the dense cluster leaves the sparse seed isolated.
        call obtain_ensembles_global_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, 0.15_real64, t_observables=10_int32, n_tiles=max(1_int32, min(8_int32, n_seeds)), &
            ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "per-seed radii: uniform radius execution")
        call assert_equal_int(count(ensemble_matrix(:, 2)), 1_int32, &
                              "per-seed radii: uniform tight radius strands the sparse seed")


        call obtain_ensembles_seed_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                                            density_labels, seed_indices, n_seeds, seed_radii=spread(0.15_real64, 1, n_seeds), &
                                            t_observables=10_int32, n_tiles=2_int32, ensemble_matrix=generated_matrix, &
                                            stop_reasons=generated_stops, mad_ambient=generated_mad, ierr=ierr)
        if (ierr == ERR_OK) generated_count = size(generated_matrix, dim=2, kind=int32)
        call assert_equal_int(ierr, ERR_OK, "equal per-seed radius execution")
        call assert_equal_int(generated_count, n_ensembles, "equal per-seed radius ensemble count")
        call assert_true(all(generated_matrix .eqv. ensemble_matrix), "equal per-seed radius ensemble membership")
        call assert_true(all(generated_stops == stop_reasons), "equal per-seed radius stop reasons")
        call assert_true(generated_mad == mad_ambient, "equal per-seed radius ambient MAD")

        ! The same run with a locally fitted radius per seed grows both ensembles.
        per_seed_radii = [0.15_real64, 5.5_real64]

        call obtain_ensembles_seed_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, per_seed_radii, t_observables=10_int32, &
            n_tiles=max(1_int32, min(8_int32, n_seeds)), ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, &
            mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_equal_int(ierr, ERR_OK, "per-seed radii: per-seed radius execution")
        call assert_equal_int(count(ensemble_matrix(:, 1)), 3_int32, &
                              "per-seed radii: dense ensemble still grows tightly")
        call assert_true(count(ensemble_matrix(:, 2)) > 1_int32, &
                         "per-seed radii: sparse ensemble grows at its own scale")


        ! A negative radius is rejected rather than silently clamped.
        per_seed_radii = [0.15_real64, -1.0_real64]

        call obtain_ensembles_seed_radius(vectors, n_dims, n_vecs, dimension_order, kd_indices, density_labels, &
            seed_indices, n_seeds, per_seed_radii, t_observables=10_int32, &
            n_tiles=max(1_int32, min(8_int32, n_seeds)), ensemble_matrix=ensemble_matrix, stop_reasons=stop_reasons, &
            mad_ambient=mad_ambient, ierr=ierr)
        if (ierr == ERR_OK) n_ensembles = size(ensemble_matrix, dim=2, kind=int32)

        call assert_true(ierr /= ERR_OK, "per-seed radii: negative radius must fail")

    end subroutine test_obtain_ensembles_per_seed_radii

end module mod_test_shatter_cluster_data
