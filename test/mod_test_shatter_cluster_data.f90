! filepath: test/mod_test_shatter_cluster_data.f90
!> Unit test suite for shatter clustering data calculations.
module mod_test_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use asserts
    use test_suite, only: test_case
    use tox_errors, only: is_ok, ERR_OK
    use tox_shatter_cluster_data, only: calculate_density_radius, &
                                        calculate_labels_as_density

    implicit none
    public

contains

    !> Get array registry of all available suite tests.
    function get_all_tests_shatter_cluster_data() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(2))

        all_tests(1) = test_case("test_density_radius_basic", test_density_radius_basic)
        all_tests(2) = test_case("test_density_labels_basic", test_density_labels_basic)
    end function get_all_tests_shatter_cluster_data

    !> Validates baseline density radius configuration with 10 vectors.
    subroutine test_density_radius_basic()
        ! 1. Arrange
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs)
        real(real64) :: tmp_mean_vec(n_dims)
        real(real64) :: tmp_distances(n_vecs)
        integer(int32) :: tmp_perm(n_vecs)
        real(real64) :: radius
        integer(int32) :: ierr
        integer(int32) :: i_setup

        ! Establish a uniform sequential spacing track from (1,1) up to (10,10)
        do i_setup = 1, n_vecs
            vectors(:, i_setup) = [real(i_setup, real64), real(i_setup, real64)]
        end do

        ! 2. Act
        call calculate_density_radius(vectors, n_dims, n_vecs, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius, 0.50_real64, ierr)

        ! 3. Assert
        call assert_equal_int(ierr, ERR_OK, "test_density_radius_basic: ierr success control check")
        call assert_in_range_real(radius, 0.0_real64, 10.0_real64, "test_density_radius_basic: radius range verification")
    end subroutine test_density_radius_basic

    !> Validates cluster neighborhood density counting logic with a 10-vector landscape.
    subroutine test_density_labels_basic()
        ! 1. Arrange
        integer(int32), parameter :: n_dims = 2_int32
        integer(int32), parameter :: n_vecs = 10_int32

        real(real64) :: vectors(n_dims, n_vecs)
        integer(int32) :: kd_indices(n_vecs)
        real(real64) :: label_densities(n_vecs)
        integer(int32) :: ierr
        integer(int32) :: i_chk

        ! Member 1 to 6: A tightly packed cluster near the origin
        vectors(:, 1) = [0.0_real64, 0.0_real64]
        vectors(:, 2) = [0.0_real64, 0.1_real64]
        vectors(:, 3) = [0.1_real64, 0.0_real64]
        vectors(:, 4) = [0.1_real64, 0.1_real64]
        vectors(:, 5) = [0.0_real64, 0.2_real64]
        vectors(:, 6) = [0.2_real64, 0.0_real64]

        ! Member 7 to 10: Completely isolated sparse points placed far away
        vectors(:, 7) = [10.0_real64, 10.0_real64]
        vectors(:, 8) = [20.0_real64, 20.0_real64]
        vectors(:, 9) = [30.0_real64, 30.0_real64]
        vectors(:, 10) = [40.0_real64, 40.0_real64]

        ! Map linear KD indexing track
        do i_chk = 1, n_vecs
            kd_indices(i_chk) = i_chk
        end do

        ! 2. Act
        call calculate_labels_as_density(vectors, n_dims, n_vecs, 0.5_real64, &
                                         kd_indices, label_densities, ierr)

        ! 3. Assert
        call assert_equal_int(ierr, ERR_OK, "test_density_labels_basic: ierr execution check")

        ! Members 1 to 6 should find exactly 6 neighbors (including itself)
        ! Enforces absolute zero-tolerance on exact integer density counts
        do i_chk = 1, 6
            call assert_equal_real(label_densities(i_chk), 6.0_real64, 0.0_real64, "Dense cluster tracking verification")
        end do

        ! Members 7 to 10 should find exactly 1 neighbor (only itself)
        ! Enforces absolute zero-tolerance on exact integer density counts
        do i_chk = 7, 10
            call assert_equal_real(label_densities(i_chk), 1.0_real64, 0.0_real64, "Isolated point density verification")
        end do
    end subroutine test_density_labels_basic

end module mod_test_shatter_cluster_data
