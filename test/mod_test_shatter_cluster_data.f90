! filepath: test/mod_test_shatter_cluster_data.f90
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
                                        calculate_labels_as_density_alloc

    implicit none
    public

contains

    function get_all_tests_shatter_cluster_data() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(11))
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
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_l_stack(n_vecs)
        integer(int32) :: tmp_r_stack(n_vecs), tmp_rec_stack(3, n_vecs), ierr, i_chk
        logical :: tmp_vicinity_mask(n_vecs, n_vecs)
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
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_l_stack, tmp_r_stack, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_basic: tree construction check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask, label_densities, ierr)
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
        logical :: tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 1.0_real64]
        vectors(:, 3) = [2.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call calculate_labels_as_density(vectors, n_dims, n_vecs, -1.0_real64, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask, label_densities, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_labels_invalid_r: Negative search radius must fail")
    end subroutine test_density_labels_invalid_r

    subroutine test_density_labels_invalid_kd_indices()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 3_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs), ierr
        logical :: tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [1.0_real64, 1.0_real64]
        vectors(:, 3) = [2.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        kd_indices = [1_int32, 0_int32, 3_int32]
        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask, label_densities, ierr)
        call assert_true(ierr /= ERR_OK, "test_density_labels_invalid_kd_indices: Zero index must fail")

        kd_indices = [1_int32, 2_int32, 4_int32]
        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask, label_densities, ierr)
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
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_l_stack(n_vecs)
        integer(int32) :: tmp_r_stack(n_vecs), tmp_rec_stack(3, n_vecs), ierr, i_chk
        logical :: tmp_vicinity_mask(n_vecs, n_vecs)
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
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_l_stack, tmp_r_stack, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_small_radius: tree construction check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.05_real64, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask, label_densities, ierr)
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
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_l_stack(n_vecs)
        integer(int32) :: tmp_r_stack(n_vecs), tmp_rec_stack(3, n_vecs), ierr
        logical :: tmp_vicinity_mask(n_vecs, n_vecs)

        vectors(:, 1) = [1.0_real64, 2.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_l_stack, tmp_r_stack, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_single_vector: tree build check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask, label_densities, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_single_vector: execution check")

        call assert_equal_real(label_densities(1), 1.0_real64, 0.0_real64, &
                               "test_density_labels_single_vector: single point density tracking")
    end subroutine test_density_labels_single_vector

    subroutine test_density_labels_identical_vectors()
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 5_int32

        real(real64) :: vectors(n_dims, n_vecs), label_densities(n_vecs), tmp_val_buf(n_vecs)
        integer(int32) :: dimension_order(n_dims), kd_indices(n_vecs), tmp_stack(3, 64, n_vecs)
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_l_stack(n_vecs)
        integer(int32) :: tmp_r_stack(n_vecs), tmp_rec_stack(3, n_vecs), ierr, i_chk
        logical :: tmp_vicinity_mask(n_vecs, n_vecs)
        character(len=128) :: assert_msg

        do i_chk = 1, n_vecs
            vectors(:, i_chk) = [4.2_real64, 4.2_real64]
        end do
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(vectors, n_dims, n_vecs, kd_indices, dimension_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_l_stack, tmp_r_stack, &
                            tmp_rec_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_identical_vectors: tree build check")

        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.1_real64, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask, label_densities, ierr)
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
        integer(int32) :: tmp_workspace(n_vecs), tmp_perm_kd(n_vecs), tmp_l_stack(n_vecs)
        integer(int32) :: tmp_r_stack(n_vecs), tmp_rec_stack(3, n_vecs), ierr, i_chk
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
                            tmp_workspace, tmp_val_buf, tmp_perm_kd, tmp_l_stack, tmp_r_stack, &
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

end module mod_test_shatter_cluster_data
