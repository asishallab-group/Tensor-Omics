!> The `get_family_centroid` cases: family k's centroid is column k of `family_centroids`; a
!| family index outside 1..n_families is refused with ERR_INVALID_INPUT, an `out_centroid` of
!| another length than the samples with ERR_SIZE_MISMATCH, in that order; a refusal leaves
!| `out_centroid` alone.
module mod_test_get_family_centroid
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_data_accessors, only: get_family_centroid
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_data_accessors_fixtures, only: N_SAMPLES, N_FAMILIES, UNWRITTEN, family_centroids
    implicit none
    private
    public :: get_all_tests_get_family_centroid

contains

    !> Every case of `get_family_centroid`.
    function get_all_tests_get_family_centroid() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(4))
        all_tests(1) = test_case("test_get_family_centroid_values", test_get_family_centroid_values)
        all_tests(2) = test_case("test_get_family_centroid_family_bound", test_get_family_centroid_family_bound)
        all_tests(3) = test_case("test_get_family_centroid_length", test_get_family_centroid_length)
        all_tests(4) = test_case("test_get_family_centroid_no_sample", test_get_family_centroid_no_sample)
    end function get_all_tests_get_family_centroid

    !> Kept from the old `test_data_accessors` (family 2's centroid of the shipped data, compared
    !| with the column it was read from): each family's column, [100*k + 1, 100*k + 2]: [101, 102]
    !| and [201, 202].
    subroutine test_get_family_centroid_values()
        real(real64) :: centroids(N_SAMPLES, N_FAMILIES), centroid(N_SAMPLES), expected(N_SAMPLES, N_FAMILIES)
        integer(int32) :: i_family, ierr

        centroids = family_centroids()
        expected = reshape([101.0_real64, 102.0_real64, 201.0_real64, 202.0_real64], [N_SAMPLES, N_FAMILIES])
        do i_family = 1, N_FAMILIES
            centroid = UNWRITTEN
            call get_family_centroid(i_family, centroids, centroid, ierr)
            call assert_err(ierr, ERR_OK, "test_get_family_centroid_values: ierr of family "//i_family)
            call assert_equal_array_real(centroid, expected(:, i_family), N_SAMPLES, 0.0_real64, &
                                         "test_get_family_centroid_values: family "//i_family)
        end do
    end subroutine test_get_family_centroid_values

    !> 0 and n_families + 1 are refused with ERR_INVALID_INPUT and leave `out_centroid` alone;
    !| so is the family index when the centroid's length is wrong too, which is checked second.
    subroutine test_get_family_centroid_family_bound()
        real(real64) :: centroids(N_SAMPLES, N_FAMILIES), centroid(N_SAMPLES), expected(N_SAMPLES), too_long(N_SAMPLES + 1)
        integer(int32) :: ierr

        centroids = family_centroids()
        expected = UNWRITTEN
        centroid = UNWRITTEN
        call get_family_centroid(0, centroids, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_family_centroid_family_bound: 0")
        call assert_equal_array_real(centroid, expected, N_SAMPLES, 0.0_real64, "test_get_family_centroid_family_bound: 0 writes")
        call get_family_centroid(N_FAMILIES + 1, centroids, centroid, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_family_centroid_family_bound: n_families + 1")
        call assert_equal_array_real(centroid, expected, N_SAMPLES, 0.0_real64, &
                                     "test_get_family_centroid_family_bound: n_families + 1 writes")
        call get_family_centroid(0, centroids, too_long, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_family_centroid_family_bound: 0, and a centroid too long")
    end subroutine test_get_family_centroid_family_bound

    !> An `out_centroid` one longer or one shorter than the samples is refused with
    !| ERR_SIZE_MISMATCH.
    subroutine test_get_family_centroid_length()
        real(real64) :: centroids(N_SAMPLES, N_FAMILIES), too_long(N_SAMPLES + 1), too_short(N_SAMPLES - 1)
        integer(int32) :: ierr

        centroids = family_centroids()
        call get_family_centroid(1, centroids, too_long, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_family_centroid_length: one too long")
        call get_family_centroid(1, centroids, too_short, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_get_family_centroid_length: one too short")
    end subroutine test_get_family_centroid_length

    !> Families with no sample have a centroid of no value: an empty `out_centroid` is the right
    !| length.
    subroutine test_get_family_centroid_no_sample()
        real(real64) :: centroids(0, N_FAMILIES), centroid(0)
        integer(int32) :: ierr

        call get_family_centroid(1, centroids, centroid, ierr)
        call assert_err(ierr, ERR_OK, "test_get_family_centroid_no_sample")
    end subroutine test_get_family_centroid_no_sample
end module mod_test_get_family_centroid
