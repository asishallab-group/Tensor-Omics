!> Unit test suite for shatter clustering data calculations.
module mod_test_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use asserts
    use test_suite, only: test_case
    use tox_errors, only: ERR_OK
    use f42_kd_tree, only: build_kd_index
    use tox_shatter_cluster_data, only: calculate_density_radius, &
                                        calculate_labels_as_density, &
                                        calculate_density_radius_alloc, &
                                        calculate_labels_as_density_alloc, &
                                        identify_ensemble_seeds, &
                                        identify_ensemble_seeds_alloc, &
                                        grow_ensemble, &
                                        grow_ensemble_alloc, &
                                        compute_ensemble_observable, &
                                        compute_ensemble_observable_alloc, &
                                        accept_ensemble, &
                                        accept_ensemble_helper, &
                                        obtain_ensembles_alloc, &
                                        obtain_ensembles, &
                                        merge_ensembles_alloc, &
                                        merge_ensembles

    implicit none
    public

contains

    function get_all_tests_shatter_cluster_data() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(46))
        all_tests(1) = test_case("test_density_radius_basic", test_density_radius_basic)
        all_tests(2) = test_case("test_density_labels_basic", test_density_labels_basic)
        all_tests(3) = test_case("test_density_radius_invalid_quantile", test_density_radius_invalid_quantile)
        all_tests(4) = test_case("test_density_labels_invalid_r", test_density_labels_invalid_r)
        all_tests(5) = test_case("test_density_labels_invalid_kd_indices", test_density_labels_invalid_kd_indices)
        all_tests(6) = test_case("test_density_radius_default_quantile", test_density_radius_default_quantile)
        all_tests(7) = test_case("test_density_labels_small_radius", test_density_labels_small_radius)
        all_tests(8) = test_case("test_density_labels_single_vector", test_density_labels_single_vector)
        all_tests(9) = test_case("test_density_labels_identical_vectors", test_density_labels_identical_vectors)
        all_tests(10) = test_case("test_density_radius_alloc", test_density_radius_alloc)
        all_tests(11) = test_case("test_density_labels_alloc", test_density_labels_alloc)
        all_tests(12) = test_case("test_identify_ensemble_seeds_basic", test_identify_ensemble_seeds_basic)
        all_tests(13) = test_case("test_identify_ensemble_seeds_k_seeding_effect", test_identify_ensemble_seeds_k_seeding_effect)
        all_tests(14) = test_case("test_identify_ensemble_seeds_invalid_k_seeding", test_identify_ensemble_seeds_invalid_k_seeding)
        all_tests(15) = test_case("test_identify_ensemble_seeds_alloc", test_identify_ensemble_seeds_alloc)
        all_tests(16) = test_case("test_identify_ensemble_seeds_single_vector_invalid", test_identify_ensemble_seeds_single_vector_invalid)
        all_tests(17) = test_case("test_identify_ensemble_seeds_invalid_inputs", test_identify_ensemble_seeds_invalid_inputs)
        all_tests(18) = test_case("test_identify_ensemble_seeds_identical_vectors", test_identify_ensemble_seeds_identical_vectors)
        all_tests(19) = test_case("test_grow_ensemble_basic", test_grow_ensemble_basic)
        all_tests(20) = test_case("test_grow_ensemble_no_growth", test_grow_ensemble_no_growth)
        all_tests(21) = test_case("test_grow_ensemble_invalid_inputs", test_grow_ensemble_invalid_inputs)
        all_tests(22) = test_case("test_grow_ensemble_alloc", test_grow_ensemble_alloc)
        all_tests(23) = test_case("test_grow_ensemble_multistep", test_grow_ensemble_multistep)
        all_tests(24) = test_case("test_grow_ensemble_alpha_sensitivity", test_grow_ensemble_alpha_sensitivity)
        all_tests(25) = test_case("test_grow_ensemble_zero_mad", test_grow_ensemble_zero_mad)
        all_tests(26) = test_case("test_grow_ensemble_mask_states", test_grow_ensemble_mask_states)
        all_tests(27) = test_case("test_compute_ensemble_observable_basic", test_compute_ensemble_observable_basic)
        all_tests(28) = test_case("test_compute_ensemble_observable_sliding_window", test_compute_ensemble_observable_sliding_window)
        all_tests(29) = test_case("test_compute_ensemble_observable_infinite_history", test_compute_ensemble_observable_infinite_history)
        all_tests(30) = test_case("test_compute_ensemble_observable_empty_and_zero", test_compute_ensemble_observable_empty_and_zero)
        all_tests(31) = test_case("test_compute_ensemble_observable_alloc", test_compute_ensemble_observable_alloc)
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
    end function get_all_tests_shatter_cluster_data

    subroutine test_density_radius_basic()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), tmp_mean_vec(n_dims), tmp_distances(n_vecs), radius
        integer(int32) :: tmp_perm(n_vecs), ierr, i_setup

        do i_setup = 1, n_vecs
            vectors(:, i_setup) = [real(i_setup, real64), real(i_setup, real64)]
        end do

        call calculate_density_radius(vectors, n_dims, n_vecs, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius, 0.50_real64, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_density_radius_basic: ierr success control check")
        call assert_equal_real(radius, 3.5355339059327378_real64, 1.0e-12_real64, &
                               "test_density_radius_basic: exact analytical radius verification")
    end subroutine test_density_radius_basic

    subroutine test_density_labels_basic()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs)
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

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_basic: tree construction check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)
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

        call calculate_density_radius(vectors, n_dims, n_vecs, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius, -0.1_real64, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_radius_invalid_quantile: Negative quantile must fail")

        call calculate_density_radius(vectors, n_dims, n_vecs, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius, 1.1_real64, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_radius_invalid_quantile: Quantile > 1.0 must fail")
    end subroutine test_density_radius_invalid_quantile

    subroutine test_density_labels_invalid_r()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), ierr

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 1.0_real64]
        vectors(:, 3) = [2.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call calculate_labels_as_density(vectors, n_dims, n_vecs, -1.0_real64, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_labels_invalid_r: Negative search radius must fail")
    end subroutine test_density_labels_invalid_r

    subroutine test_density_labels_invalid_kd_indices()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), ierr

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 1.0_real64]
        vectors(:, 3) = [2.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        kd_indices = [1_int32, 0_int32, 3_int32]
        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_labels_invalid_kd_indices: Zero index must fail")

        kd_indices = [1_int32, 2_int32, 4_int32]
        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)
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

        call calculate_density_radius(vectors, n_dims, n_vecs, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius_default, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_radius_default_quantile: default call check")

        call calculate_density_radius(vectors, n_dims, n_vecs, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius_explicit, 0.15_real64, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_radius_default_quantile: explicit call check")

        call assert_equal_real(radius_default, radius_explicit, 1.0e-12_real64, &
                               "test_density_radius_default_quantile: default matches 0.15 option")
    end subroutine test_density_radius_default_quantile

    subroutine test_density_labels_small_radius()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs)
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

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_small_radius: tree construction check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.05_real64, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)
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
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr

        vectors(:, 1) = [1.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_single_vector: tree build check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_single_vector: execution check")

        call assert_equal_real(label_densities(1), 1.0_real64, 0.0_real64, &
                               "test_density_labels_single_vector: single point density tracking")
    end subroutine test_density_labels_single_vector

    subroutine test_density_labels_identical_vectors()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 5_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr, i_chk
        character(len=128) :: assert_msg

        do i_chk = 1, n_vecs
            vectors(:, i_chk) = [4.2_real64, 4.2_real64]
        end do
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_identical_vectors: tree build check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.1_real64, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_identical_vectors: execution check")

        do i_chk = 1, n_vecs
            write (assert_msg, '(A,I0)') "test_density_labels_identical_vectors: identical coordinates check at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), real(n_vecs, real64), 0.0_real64, assert_msg)
        end do
    end subroutine test_density_labels_identical_vectors

    subroutine test_density_radius_alloc()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs), radius
        integer(int32) :: ierr, i_setup

        do i_setup = 1, n_vecs
            vectors(:, i_setup) = [real(i_setup, real64), real(i_setup, real64)]
        end do

        call calculate_density_radius_alloc(vectors, n_dims, n_vecs, &
                                            radius, 0.50_real64, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_density_radius_alloc: ierr success check")
        call assert_equal_real(radius, 3.5355339059327378_real64, 1.0e-12_real64, &
                               "test_density_radius_alloc: analytical radius check")
    end subroutine test_density_radius_alloc

    subroutine test_density_labels_alloc()
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

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_alloc: tree build check")

        call calculate_labels_as_density_alloc(vectors, n_dims, n_vecs, 0.5_real64, &
                                               dimension_order, kd_indices, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_alloc: execution check")

        do i_chk = 1, 6
            write (assert_msg, '(A,I0)') "test_density_labels_alloc: dense cluster check at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), 6.0_real64, 0.0_real64, assert_msg)
        end do

        do i_chk = 7, 10
            write (assert_msg, '(A,I0)') "test_density_labels_alloc: isolated points check at index: ", i_chk
            call assert_equal_real(label_densities(i_chk), 1.0_real64, 0.0_real64, assert_msg)
        end do
    end subroutine test_density_labels_alloc

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
        logical :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [10.0_real64, 0.0_real64]
        vectors(:, 5) = [11.0_real64, 0.0_real64]
        vectors(:, 6) = [12.0_real64, 0.0_real64]

        density_labels = [60.0_real64, 50.0_real64, 40.0_real64, &
                          30.0_real64, 20.0_real64, 10.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_basic: tree build")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, k_seeding, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_basic: execution")
        call assert_true(all(sorted_perm == [1_int32, 2_int32, 3_int32, &
                                             4_int32, 5_int32, 6_int32]), &
                         "test_identify_ensemble_seeds_basic: density ranking")
        call assert_equal_int(n_seeds, 2_int32, &
                              "test_identify_ensemble_seeds_basic: one seed per separated region")
        call assert_true(seed_mask(1), &
                         "test_identify_ensemble_seeds_basic: densest point in first region is seed")
        call assert_false(seed_mask(2), &
                          "test_identify_ensemble_seeds_basic: covered dense point is skipped")
        call assert_false(seed_mask(3), &
                          "test_identify_ensemble_seeds_basic: covered first-region point is skipped")
        call assert_true(seed_mask(4), &
                         "test_identify_ensemble_seeds_basic: densest uncovered second region becomes seed")
        call assert_false(seed_mask(5), &
                          "test_identify_ensemble_seeds_basic: covered second-region point is skipped")
        call assert_false(seed_mask(6), &
                          "test_identify_ensemble_seeds_basic: covered second-region point is skipped")
        call assert_true(all(tmp_visited_mask), &
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
        logical :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs)
        logical :: seed_mask_small(n_vecs), seed_mask_large(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [10.0_real64, 0.0_real64]
        vectors(:, 5) = [11.0_real64, 0.0_real64]
        vectors(:, 6) = [12.0_real64, 0.0_real64]

        density_labels = [60.0_real64, 50.0_real64, 40.0_real64, &
                          30.0_real64, 20.0_real64, 10.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_k_seeding_effect: tree build")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, 1_int32, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds_small, seed_mask_small, ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=1 execution")
        call assert_equal_int(n_seeds_small, 4_int32, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=1 seed count")
        call assert_true(seed_mask_small(1), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 first seed")
        call assert_true(seed_mask_small(3), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 uncovered first-region tail")
        call assert_true(seed_mask_small(4), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 second-region seed")
        call assert_true(seed_mask_small(6), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=1 uncovered second-region tail")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, 3_int32, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds_large, seed_mask_large, ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=3 execution")
        call assert_equal_int(n_seeds_large, 2_int32, &
                              "test_identify_ensemble_seeds_k_seeding_effect: k=3 seed count")
        call assert_true(seed_mask_large(1), &
                         "test_identify_ensemble_seeds_k_seeding_effect: k=3 first-region seed")
        call assert_true(seed_mask_large(4), &
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
        logical :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [3.0_real64, 0.0_real64]

        density_labels = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32, 4_int32]

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, 0_int32, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_k_seeding: k=0 must fail")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, n_vecs, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_k_seeding: k>=n_vectors must fail")
    end subroutine test_identify_ensemble_seeds_invalid_k_seeding

    subroutine test_identify_ensemble_seeds_alloc()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 6_int32
        integer(int32), parameter :: k_seeding = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), sorted_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), n_seeds, ierr
        logical :: seed_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [10.0_real64, 0.0_real64]
        vectors(:, 5) = [11.0_real64, 0.0_real64]
        vectors(:, 6) = [12.0_real64, 0.0_real64]

        density_labels = [60.0_real64, 50.0_real64, 40.0_real64, &
                          30.0_real64, 20.0_real64, 10.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_alloc: tree build")

        call identify_ensemble_seeds_alloc(vectors, n_dims, n_vecs, density_labels, &
                                           dimension_order, kd_indices, k_seeding, &
                                           sorted_perm, n_seeds, seed_mask, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_identify_ensemble_seeds_alloc: execution")
        call assert_equal_int(n_seeds, 2_int32, &
                              "test_identify_ensemble_seeds_alloc: expected two seeds")
        call assert_true(seed_mask(1), &
                         "test_identify_ensemble_seeds_alloc: first-region density centre is seed")
        call assert_true(seed_mask(4), &
                         "test_identify_ensemble_seeds_alloc: second-region density centre is seed")
        call assert_false(seed_mask(2), &
                          "test_identify_ensemble_seeds_alloc: covered first-region point is skipped")
        call assert_false(seed_mask(3), &
                          "test_identify_ensemble_seeds_alloc: covered first-region point is skipped")
        call assert_false(seed_mask(5), &
                          "test_identify_ensemble_seeds_alloc: covered second-region point is skipped")
        call assert_false(seed_mask(6), &
                          "test_identify_ensemble_seeds_alloc: covered second-region point is skipped")
    end subroutine test_identify_ensemble_seeds_alloc

    subroutine test_identify_ensemble_seeds_single_vector_invalid()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 1_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_distances(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64)
        integer(int32) :: tmp_perm(n_vecs), sorted_perm(n_vecs), n_seeds, ierr
        logical :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        density_labels = [1.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32]

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, 1_int32, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)

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
        logical :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 0.0_real64]
        vectors(:, 3) = [2.0_real64, 0.0_real64]
        vectors(:, 4) = [3.0_real64, 0.0_real64]

        density_labels = [4.0_real64, 3.0_real64, 2.0_real64, -1.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32, 4_int32]

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, k_seeding, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_inputs: negative density must fail")

        density_labels = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        kd_indices = [0_int32, 2_int32, 3_int32, 4_int32]

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, k_seeding, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)
        call assert_true(ierr /= ERR_OK, &
                         "test_identify_ensemble_seeds_invalid_inputs: invalid kd index must fail")

        kd_indices = [1_int32, 2_int32, 3_int32, 4_int32]
        dimension_order = [0_int32, 2_int32]

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, k_seeding, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)
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
        logical :: tmp_visited_mask(n_vecs), tmp_newly_covered_mask(n_vecs), seed_mask(n_vecs)

        vectors = 5.0_real64
        density_labels = [1.0_real64, 4.0_real64, 3.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_identical_vectors: tree build")

        call identify_ensemble_seeds(vectors, n_dims, n_vecs, density_labels, &
                                     dimension_order, kd_indices, k_seeding, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)

        call assert_equal_int(ierr, ERR_OK, &
                              "test_identify_ensemble_seeds_identical_vectors: execution")
        call assert_equal_int(sorted_perm(1), 2_int32, &
                              "test_identify_ensemble_seeds_identical_vectors: densest vector ranks first")
        call assert_equal_int(n_seeds, 1_int32, &
                              "test_identify_ensemble_seeds_identical_vectors: zero-radius coverage gives one seed")
        call assert_true(seed_mask(2), &
                         "test_identify_ensemble_seeds_identical_vectors: densest vector is the seed")
        call assert_false(seed_mask(1), &
                          "test_identify_ensemble_seeds_identical_vectors: coincident point is covered")
        call assert_false(seed_mask(3), &
                          "test_identify_ensemble_seeds_identical_vectors: coincident point is covered")
        call assert_false(seed_mask(4), &
                          "test_identify_ensemble_seeds_identical_vectors: coincident point is covered")
        call assert_true(all(tmp_visited_mask), &
                         "test_identify_ensemble_seeds_identical_vectors: all coincident points visited")
    end subroutine test_identify_ensemble_seeds_identical_vectors

    subroutine test_grow_ensemble_basic()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 4_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.2_real64]
        vectors(:, 4) = [0.0_real64, 0.4_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64, 1.0_real64]
        ensemble_mask = [.true., .false., .false., .false.]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_basic: tree construction check")

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, &
                           dimension_order, kd_indices, density_labels, 0.5_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_basic: execution ierr check")
        call assert_true(ensemble_mask(1), "test_grow_ensemble_basic: vector 1 must remain a member")
        call assert_true(ensemble_mask(2), "test_grow_ensemble_basic: vector 2 must be added")
        call assert_false(ensemble_mask(3), "test_grow_ensemble_basic: vector 3 too far away")
        call assert_false(ensemble_mask(4), "test_grow_ensemble_basic: vector 4 density too low")
    end subroutine test_grow_ensemble_basic

    subroutine test_grow_ensemble_no_growth()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [100.0_real64, 1.0_real64, 1.0_real64]
        ensemble_mask = [.true., .false., .false.]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, &
                           dimension_order, kd_indices, density_labels, 0.5_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_no_growth: execute ierr success")
        call assert_true(ensemble_mask(1), "test_grow_ensemble_no_growth: vector 1 stays member")
        call assert_false(ensemble_mask(2), "test_grow_ensemble_no_growth: vector 2 rejected due to density gap")
        call assert_false(ensemble_mask(3), "test_grow_ensemble_no_growth: vector 3 rejected")
    end subroutine test_grow_ensemble_no_growth

    subroutine test_grow_ensemble_invalid_inputs()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), tmp_perm(n_vecs), ierr
        logical :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false.]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, -0.1_real64, &
                           dimension_order, kd_indices, density_labels, 0.5_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)
        call assert_true(ierr /= ERR_OK, "test_grow_ensemble_invalid_inputs: negative r must fail")

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.5_real64, &
                           dimension_order, kd_indices, density_labels, -1.0_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)
        call assert_true(ierr /= ERR_OK, "test_grow_ensemble_invalid_inputs: negative alpha_mad must fail")
    end subroutine test_grow_ensemble_invalid_inputs

    subroutine test_grow_ensemble_alloc()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical :: ensemble_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_alloc: build tree check")

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble_alloc(vectors, n_dims, n_vecs, ensemble_mask, &
                                 dimension_order, kd_indices, density_labels, &
                                 r=0.6_real64, alpha_mad=0.5_real64, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_alloc: explicit r and alpha_mad ierr check")
        call assert_true(ensemble_mask(1), "Case 1: vector 1 is active")
        call assert_true(ensemble_mask(2), "Case 1: vector 2 is successfully added")
        call assert_false(ensemble_mask(3), "Case 1: vector 3 is too far away")

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble_alloc(vectors, n_dims, n_vecs, ensemble_mask, &
                                 dimension_order, kd_indices, density_labels, &
                                 alpha_mad=0.5_real64, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_alloc: omitted r ierr check")

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble_alloc(vectors, n_dims, n_vecs, ensemble_mask, &
                                 dimension_order, kd_indices, density_labels, &
                                 r=0.6_real64, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_alloc: omitted alpha_mad ierr check")

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble_alloc(vectors, n_dims, n_vecs, ensemble_mask, &
                                 dimension_order, kd_indices, density_labels, &
                                 ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_alloc: fully default call ierr check")
    end subroutine test_grow_ensemble_alloc

    subroutine test_grow_ensemble_multistep()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical :: ensemble_mask(n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false.]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        call grow_ensemble_alloc(vectors, n_dims, n_vecs, ensemble_mask, &
                                 dimension_order, kd_indices, density_labels, &
                                 r=0.6_real64, alpha_mad=0.5_real64, ierr=ierr)
        call assert_true(ensemble_mask(2), "step 1: vector 2 absorbed")
        call assert_false(ensemble_mask(3), "step 1: vector 3 not yet reachable")

        call grow_ensemble_alloc(vectors, n_dims, n_vecs, ensemble_mask, &
                                 dimension_order, kd_indices, density_labels, &
                                 r=0.6_real64, alpha_mad=0.5_real64, ierr=ierr)
        call assert_true(ensemble_mask(3), "step 2: vector 3 pulled in by vector 2 surface")
    end subroutine test_grow_ensemble_multistep

    subroutine test_grow_ensemble_alpha_sensitivity()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [10.0_real64, 8.5_real64, 12.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, &
                           dimension_order, kd_indices, density_labels, 0.25_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)
        call assert_equal_int(ierr, ERR_OK, "alpha 0.25 test ierr control")
        call assert_false(ensemble_mask(2), "alpha 0.25 should reject vector 2")

        ensemble_mask = [.true., .false., .false.]
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, &
                           dimension_order, kd_indices, density_labels, 1.00_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)
        call assert_equal_int(ierr, ERR_OK, "alpha 1.00 test ierr control")
        call assert_true(ensemble_mask(2), "alpha 1.00 should accept vector 2")
    end subroutine test_grow_ensemble_alpha_sensitivity

    subroutine test_grow_ensemble_zero_mad()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        ensemble_mask = [.true., .false., .false.]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, &
                           dimension_order, kd_indices, density_labels, 0.5_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_grow_ensemble_zero_mad: execution check")
        call assert_true(ensemble_mask(1), "test_grow_ensemble_zero_mad: vector 1 active")
        call assert_true(ensemble_mask(2), "test_grow_ensemble_zero_mad: vector 2 absorbed under zero MAD")
    end subroutine test_grow_ensemble_zero_mad

    subroutine test_grow_ensemble_mask_states()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), density_labels(n_vecs), tmp_val_buf(n_vecs), tmp_abs_diff(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), tmp_perm(n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs)
        integer(int32) :: tmp_rec_stack(3, n_vecs), ierr
        logical :: ensemble_mask(n_vecs), tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.5_real64]
        vectors(:, 3) = [0.0_real64, 1.0_real64]

        density_labels = [5.0_real64, 5.0_real64, 5.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        ensemble_mask = [.false., .false., .false.]
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, &
                           dimension_order, kd_indices, density_labels, 0.5_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)
        call assert_equal_int(ierr, ERR_OK, "empty mask check ierr control")
        call assert_false(any(ensemble_mask), "empty mask should stay empty")

        ensemble_mask = [.true., .true., .true.]
        call grow_ensemble(vectors, n_dims, n_vecs, ensemble_mask, 0.6_real64, &
                           dimension_order, kd_indices, density_labels, 0.5_real64, &
                           tmp_stack, tmp_vicinity_mask, tmp_perm, tmp_abs_diff, ierr)
        call assert_equal_int(ierr, ERR_OK, "full mask check ierr control")
        call assert_true(all(ensemble_mask), "full mask should stay full")
    end subroutine test_grow_ensemble_mask_states

    subroutine test_compute_ensemble_observable_basic()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 10)
        logical :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 5.0_real64, 20.0_real64]
        ensemble_mask = [.true., .true., .false.]
        observables = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         candidate_count=2_int32, current_iter=1_int32, &
                                         observables=observables, t_observables=10_int32, ierr=ierr)

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
        logical :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 5.0_real64, 20.0_real64, 2.0_real64]
        observables = 0.0_real64

        ensemble_mask = [.true., .true., .false., .false.]
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         0_int32, 1_int32, observables, &
                                         t_observables=2_int32, ierr=ierr)

        ensemble_mask = [.true., .true., .true., .false.]
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         1_int32, 2_int32, observables, &
                                         t_observables=2_int32, ierr=ierr)

        call assert_equal_real(observables(1, 1), 7.5_real64, 1.0e-10_real64, "col 1 has iter 1")
        call assert_equal_real(observables(1, 2), 11.666666666666666_real64, 1.0e-10_real64, "col 2 has iter 2")

        ensemble_mask = [.true., .true., .true., .true.]
        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         1_int32, 3_int32, observables, &
                                         t_observables=2_int32, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "sliding window ierr control")
        call assert_equal_real(observables(1, 1), 11.666666666666666_real64, 1.0e-10_real64, "col 1 shifted to iter 2")
        call assert_equal_real(observables(1, 2), 9.25_real64, 1.0e-10_real64, "col 2 updated to iter 3")
    end subroutine test_compute_ensemble_observable_sliding_window

    subroutine test_compute_ensemble_observable_infinite_history()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 3)
        logical :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64]
        ensemble_mask = [.true., .false., .false.]
        observables = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         0_int32, 1_int32, observables, &
                                         t_observables=0_int32, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "infinite history (0) ierr check")
        call assert_equal_real(observables(4, 1), 1.0_real64, 1.0e-10_real64, "iter 1 stored in col 1")

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         0_int32, 2_int32, observables, &
                                         t_observables=-1_int32, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "infinite history (-1) ierr check")
        call assert_equal_real(observables(4, 2), 1.0_real64, 1.0e-10_real64, "iter 2 stored in col 2")
    end subroutine test_compute_ensemble_observable_infinite_history

    subroutine test_compute_ensemble_observable_empty_and_zero()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 5)
        logical :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        observables = 0.0_real64

        density_labels = [10.0_real64, 5.0_real64, 2.0_real64]
        ensemble_mask = [.false., .false., .false.]

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         0_int32, 1_int32, observables, &
                                         t_observables=5_int32, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "empty ensemble ierr check")
        call assert_equal_real(observables(1, 1), 0.0_real64, 1.0e-10_real64, "empty ensemble outputs zeros")

        density_labels = [0.0_real64, 0.0_real64, 10.0_real64]
        ensemble_mask = [.true., .true., .false.]

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         0_int32, 2_int32, observables, &
                                         t_observables=5_int32, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "zero density ierr check")
        call assert_equal_real(observables(1, 2), 0.0_real64, 1.0e-10_real64, "rho_arith is 0.0")
        call assert_equal_real(observables(2, 2), 0.0_real64, 1.0e-10_real64, "rho_harm guarded to 0.0")
        call assert_equal_real(observables(3, 2), 1.0_real64, 1.0e-10_real64, "H_rho defaults safely to 1.0")
    end subroutine test_compute_ensemble_observable_empty_and_zero

    subroutine test_compute_ensemble_observable_alloc()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs)
        real(real64), allocatable :: observables(:, :)
        logical :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64]
        ensemble_mask = [.true., .true., .false.]

        call compute_ensemble_observable_alloc(ensemble_mask, density_labels, n_vecs, &
                                               0_int32, 1_int32, observables, &
                                               t_observables=5_int32, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "alloc wrapper ierr check")
        call assert_true(allocated(observables), "observables matrix auto-allocated")
        call assert_equal_int(size(observables, 1), 5_int32, "observables rows = 5")
        call assert_equal_int(size(observables, 2), 5_int32, "observables cols = 5")
        call assert_equal_real(observables(1, 1), 10.0_real64, 1.0e-10_real64, "arithmetic mean check")
    end subroutine test_compute_ensemble_observable_alloc

    subroutine test_compute_ensemble_observable_invalid_inputs()
        integer(int32), parameter :: n_vecs = 3_int32
        real(real64) :: density_labels(n_vecs), observables(5, 5)
        logical :: ensemble_mask(n_vecs)
        integer(int32) :: ierr

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64]
        ensemble_mask = [.true., .false., .false.]
        observables = 0.0_real64

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         candidate_count=-1_int32, current_iter=1_int32, &
                                         observables=observables, t_observables=5_int32, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "negative candidate_count must fail")

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vecs, &
                                         candidate_count=0_int32, current_iter=0_int32, &
                                         observables=observables, t_observables=5_int32, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "current_iter <= 0 must fail")
    end subroutine test_compute_ensemble_observable_invalid_inputs

    subroutine test_accept_ensemble_default()
        real(real64) :: observables(5, 5)
        logical :: is_accepted
        integer(int32) :: ierr

        observables = 0.0_real64
        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.2_real64

        call accept_ensemble(observables, 2_int32, is_accepted=is_accepted, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "accept_ensemble default ierr control")
        call assert_true(is_accepted, "log2 fold change 0.263 < default 0.5 should be accepted")

        observables(3, 3) = 1.8_real64
        call accept_ensemble(observables, 3_int32, is_accepted=is_accepted, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "accept_ensemble rejection ierr control")
        call assert_false(is_accepted, "log2 fold change 0.585 >= default 0.5 should be rejected")
    end subroutine test_accept_ensemble_default

    subroutine test_accept_ensemble_alpha_sensitivity()
        real(real64) :: observables(5, 5)
        logical :: is_accepted
        integer(int32) :: ierr

        observables = 0.0_real64
        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.5_real64

        call accept_ensemble(observables, 2_int32, alpha_accept=0.2_real64, is_accepted=is_accepted, ierr=ierr)
        call assert_false(is_accepted, "alpha 0.2 rejects 0.585 change")

        call accept_ensemble(observables, 2_int32, alpha_accept=0.4_real64, is_accepted=is_accepted, ierr=ierr)
        call assert_false(is_accepted, "alpha 0.4 rejects 0.585 change")

        call accept_ensemble(observables, 2_int32, alpha_accept=0.8_real64, is_accepted=is_accepted, ierr=ierr)
        call assert_true(is_accepted, "alpha 0.8 accepts 0.585 change")

        call accept_ensemble(observables, 2_int32, alpha_accept=1.0_real64, is_accepted=is_accepted, ierr=ierr)
        call assert_true(is_accepted, "alpha 1.0 accepts 0.585 change")

        call accept_ensemble(observables, 2_int32, alpha_accept=1.5_real64, is_accepted=is_accepted, ierr=ierr)
        call assert_true(is_accepted, "alpha 1.5 accepts 0.585 change")
    end subroutine test_accept_ensemble_alpha_sensitivity

    subroutine test_accept_ensemble_boundary_and_window()
        real(real64) :: observables(5, 2)
        logical :: is_accepted
        integer(int32) :: ierr

        observables = 0.0_real64

        call accept_ensemble(observables, 1_int32, alpha_accept=0.5_real64, is_accepted=is_accepted, ierr=ierr)
        call assert_true(is_accepted, "iter 1 must always be accepted")

        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.1_real64
        call accept_ensemble(observables, 3_int32, alpha_accept=0.5_real64, is_accepted=is_accepted, ierr=ierr)
        call assert_true(is_accepted, "sliding window check accepts smooth step")
    end subroutine test_accept_ensemble_boundary_and_window

    subroutine test_accept_ensemble_invalid_inputs()
        real(real64) :: observables(5, 5)
        logical :: is_accepted
        integer(int32) :: ierr

        observables = 1.0_real64

        call accept_ensemble(observables, current_iter=0_int32, is_accepted=is_accepted, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "current_iter <= 0 must fail")

        call accept_ensemble(observables, current_iter=2_int32, alpha_accept=-0.1_real64, &
                             is_accepted=is_accepted, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "negative alpha_accept must fail")
    end subroutine test_accept_ensemble_invalid_inputs

    subroutine test_accept_ensemble_helper_basic()

        real(real64) :: observables(5, 4)
        logical :: is_accepted

        observables = 0.0_real64

        call accept_ensemble_helper(observables, 1_int32, alpha_accept=0.5_real64, &
                                    is_accepted=is_accepted)
        call assert_true(is_accepted, "helper iter 1 must always accept")

        observables(3, 1) = 1.0_real64
        observables(3, 2) = 1.2_real64

        call accept_ensemble_helper(observables, 2_int32, alpha_accept=0.5_real64, &
                                    is_accepted=is_accepted)
        call assert_true(is_accepted, "helper accepts smooth fold change")

        observables(3, 3) = 2.0_real64

        call accept_ensemble_helper(observables, 3_int32, alpha_accept=0.5_real64, &
                                    is_accepted=is_accepted)
        call assert_false(is_accepted, "helper rejects sharp heterogeneity jump")

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
        logical, allocatable :: raw_matrix(:, :), merged_matrix(:, :)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.2_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.2_real64]
        vectors(:, 5) = [100.0_real64, 100.0_real64]
        vectors(:, 6) = [200.0_real64, 200.0_real64]

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64, 10.0_real64, &
                          1.0_real64, 1.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        seed_indices = [1_int32, 3_int32]

        allocate (raw_matrix(n_vecs, n_seeds))

        call obtain_ensembles_alloc(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                                    density_labels, seed_indices, n_seeds, &
                                    r=0.5_real64, ensemble_matrix=raw_matrix, &
                                    n_ensembles=n_raw, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "obtain_ensembles execution check")

        call merge_ensembles_alloc(raw_matrix, n_vecs, n_raw, min_intersection=1_int32, &
                                   merged_matrix=merged_matrix, n_ensembles=n_merged, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles execution check")
        call assert_equal_int(n_merged, 2_int32, "expected 2 distinct merged clusters")

        call assert_false(any(merged_matrix(5, :)), "vector 5 is background noise")
        call assert_false(any(merged_matrix(6, :)), "vector 6 is background noise")

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
        logical :: ensemble_matrix(n_vecs, n_seeds)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.1_real64]
        vectors(:, 3) = [10.0_real64, 10.0_real64]
        vectors(:, 4) = [10.0_real64, 10.1_real64]

        density_labels = [10.0_real64, 10.0_real64, 10.0_real64, 10.0_real64]
        dimension_order = [1_int32, 2_int32]
        seed_indices = [1_int32, 2_int32, 3_int32, 4_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        call obtain_ensembles_alloc(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                                    density_labels, seed_indices, n_seeds, &
                                    r=0.5_real64, ensemble_matrix=ensemble_matrix, &
                                    n_ensembles=n_ensembles, ierr=ierr)

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
        logical :: ensemble_matrix(n_vecs, n_seeds)

        vectors = 0.0_real64
        density_labels = 1.0_real64
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]
        seed_indices = [1_int32]

        call obtain_ensembles_alloc(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                                    density_labels, seed_indices, n_seeds, &
                                    r=0.5_real64, ensemble_matrix=ensemble_matrix, &
                                    n_ensembles=n_ensembles, ierr=ierr)

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
        logical :: ensemble_matrix(n_vecs, n_seeds)

        vectors = 0.0_real64
        density_labels = 1.0_real64
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]
        seed_indices = [1_int32]

        seed_indices(1) = 4_int32

        call obtain_ensembles_alloc(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                                    density_labels, seed_indices, n_seeds, &
                                    r=0.5_real64, ensemble_matrix=ensemble_matrix, &
                                    n_ensembles=n_ensembles, ierr=ierr)

        call assert_true(ierr /= ERR_OK, "seed_indices out of bounds must fail")

        seed_indices(1) = 1_int32

        call obtain_ensembles_alloc(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                                    density_labels, seed_indices, n_seeds, &
                                    r=-0.5_real64, ensemble_matrix=ensemble_matrix, &
                                    n_ensembles=n_ensembles, ierr=ierr)

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
        logical, allocatable :: raw_matrix(:, :), merged_matrix(:, :)
        character(len=128) :: assert_msg

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

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, &
                            tmp_rec_stack, ierr)

        call calculate_labels_as_density_alloc(vectors, n_dims, n_vecs, 0.5_real64, &
                                               dimension_order, kd_indices, &
                                               density_labels, ierr)

        seed_indices = [15_int32, 165_int32, 315_int32]

        allocate (raw_matrix(n_vecs, n_seeds))

        call obtain_ensembles_alloc(vectors, n_dims, n_vecs, dimension_order, kd_indices, &
                                    density_labels, seed_indices, n_seeds, &
                                    r=0.5_real64, ensemble_matrix=raw_matrix, &
                                    n_ensembles=n_raw, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "obtain_ensembles 500 execution check")

        call merge_ensembles_alloc(raw_matrix, n_vecs, n_raw, &
                                   min_intersection=1_int32, &
                                   merged_matrix=merged_matrix, &
                                   n_ensembles=n_merged, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles 500 execution check")
        call assert_equal_int(n_merged, 3_int32, &
                              "expected 3 distinct merged clusters")

        do i_vec = 451, 500
            write (assert_msg, '(A,I0)') &
                "test_obtain_ensembles_500_vectors: background noise verification at index: ", &
                i_vec
            call assert_false(any(merged_matrix(i_vec, :)), assert_msg)
        end do

    end subroutine test_obtain_ensembles_500_vectors

    subroutine test_merge_ensembles_basic()
        integer(int32), parameter :: n_vecs = 5_int32
        integer(int32), parameter :: n_seeds = 3_int32

        logical :: raw_masks(n_vecs, n_seeds)
        logical, allocatable :: merged_matrix(:, :)
        integer(int32) :: n_ensembles, ierr

        raw_masks(:, 1) = [.true., .true., .false., .false., .false.]
        raw_masks(:, 2) = [.false., .true., .true., .false., .false.]
        raw_masks(:, 3) = [.false., .false., .false., .true., .true.]

        call merge_ensembles_alloc(raw_masks, n_vecs, n_seeds, min_intersection=1_int32, &
                                   merged_matrix=merged_matrix, n_ensembles=n_ensembles, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles_basic ierr check")
        call assert_equal_int(n_ensembles, 2_int32, "masks 1 & 2 merge into 1, leaving 2 unique clusters")
        call assert_equal_int(size(merged_matrix, 2), 2_int32, "output matrix has 2 columns")

        call assert_true(merged_matrix(1, 1), "merged 1 contains elem 1")
        call assert_true(merged_matrix(2, 1), "merged 1 contains elem 2")
        call assert_true(merged_matrix(3, 1), "merged 1 contains elem 3")
        call assert_false(merged_matrix(4, 1), "merged 1 excludes elem 4")
    end subroutine test_merge_ensembles_basic

    subroutine test_merge_ensembles_zero_seeds()
        integer(int32), parameter :: n_vecs = 4_int32

        logical :: raw_masks(n_vecs, 0)
        logical, allocatable :: merged_matrix(:, :)
        integer(int32) :: n_ensembles, ierr

        call merge_ensembles_alloc(raw_masks, n_vecs, 0_int32, min_intersection=1_int32, &
                                   merged_matrix=merged_matrix, n_ensembles=n_ensembles, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles_zero_seeds ierr check")
        call assert_true(allocated(merged_matrix), "merged_matrix allocated for 0 seeds")
        call assert_equal_int(n_ensembles, 0_int32, "n_ensembles is 0")
        call assert_equal_int(size(merged_matrix, 2), 0_int32, "columns count is 0")
    end subroutine test_merge_ensembles_zero_seeds

    subroutine test_merge_ensembles_no_overlap()
        integer(int32), parameter :: n_vecs = 4_int32
        integer(int32), parameter :: n_seeds = 2_int32

        logical :: raw_masks(n_vecs, n_seeds)
        logical, allocatable :: merged_matrix(:, :)
        integer(int32) :: n_ensembles, ierr

        raw_masks(:, 1) = [.true., .true., .false., .false.]
        raw_masks(:, 2) = [.false., .false., .true., .true.]

        call merge_ensembles_alloc(raw_masks, n_vecs, n_seeds, min_intersection=1_int32, &
                                   merged_matrix=merged_matrix, n_ensembles=n_ensembles, ierr=ierr)

        call assert_equal_int(ierr, ERR_OK, "merge_ensembles_no_overlap ierr check")
        call assert_equal_int(n_ensembles, 2_int32, "disjoint masks do not merge")
    end subroutine test_merge_ensembles_no_overlap

    subroutine test_merge_ensembles_invalid_inputs()
        integer(int32), parameter :: n_vecs = 3_int32
        integer(int32), parameter :: n_seeds = 2_int32

        logical :: raw_masks(n_vecs, n_seeds)
        logical, allocatable :: merged_matrix(:, :)
        integer(int32) :: n_ensembles, ierr

        raw_masks = .true.

        call merge_ensembles_alloc(raw_masks, n_vecs, n_seeds, min_intersection=0_int32, &
                                   merged_matrix=merged_matrix, n_ensembles=n_ensembles, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "min_intersection = 0 must fail")

        call merge_ensembles_alloc(raw_masks, n_vecs, n_seeds, min_intersection=n_vecs + 1_int32, &
                                   merged_matrix=merged_matrix, n_ensembles=n_ensembles, ierr=ierr)
        call assert_true(ierr /= ERR_OK, "min_intersection > n_vectors must fail")
    end subroutine test_merge_ensembles_invalid_inputs

end module mod_test_shatter_cluster_data
