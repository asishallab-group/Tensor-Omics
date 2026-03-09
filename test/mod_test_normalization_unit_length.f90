! filepath: test/mod_test_normalization_unit_length.f90
!> Unit test suite for normalization_unit_length routine.
module mod_test_normalization_unit_length
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_normalization, only: normalize_unit_length
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use tox_errors
    use mod_test_suite, only: test_case
    implicit none

   
    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests_normalization_unit_length() result(all_tests)
        type(test_case),allocatable :: all_tests(:)

        allocate(all_tests(1))
        all_tests(1) = test_case("test_normalization_unit_length", test_normalize_unit_length)
    end function get_all_tests_normalization_unit_length

    
    

    !> Test the normalize_unit_length function with various cases.
    subroutine test_normalize_unit_length()
        integer(int32) :: ierr
        integer(int32), parameter :: n_dims = 3
        real(real64) :: vector(n_dims)
        real(real64) :: expected(n_dims)

        ! -------------------------------
        ! Case 1: Normal vector
        ! -------------------------------
        vector = [3.0_real64, 4.0_real64, -123.0_real64]
        expected = vector / sqrt(sum(vector**2))

        call normalize_unit_length(vector, n_dims, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_unit_length: normal vector ierr")
        call assert_equal_array_real(vector, expected, n_dims, TOL, "test_normalize_unit_length: normal vector result")

        ! -------------------------------
        ! Case 2: Zero vector
        ! -------------------------------
        vector = [0.0_real64, 0.0_real64, 0.0_real64]

        call normalize_unit_length(vector, n_dims, ierr)
        call assert_equal_int(ierr, ERR_DIVISION_BY_ZERO, "test_normalize_unit_length: zero vector should trigger ERR_DIVISION_BY_ZERO")

        ! -------------------------------
        ! Case 3: Already normalized
        ! -------------------------------
        vector = [0.61647708797651190566181131895072_real64, -0.42_real64, 0.666_real64]
        expected = vector

        call normalize_unit_length(vector, n_dims, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_unit_length: already normalized vector ierr")
        call assert_equal_array_real(vector, expected, n_dims, 0.0_real64, "test_normalize_unit_length: already normalized vector result")

        ! -------------------------------
        ! Case 4: NaN
        ! -------------------------------
        vector = [0.61647708797651190566181131895072_real64, -0.42_real64, ieee_value(1.0_real64, ieee_quiet_nan)]
        expected = vector

        call normalize_unit_length(vector, n_dims, ierr)
        call assert_equal_int(ierr, ERR_NAN_INF, "test_normalize_unit_length: vector with NaN should trigger ERR_NAN_INF")

        ! -------------------------------
        ! Case 5: Infinity
        ! -------------------------------
        vector = [0.61647708797651190566181131895072_real64, -0.42_real64, ieee_value(1.0_real64, ieee_positive_inf)]
        expected = vector

        call normalize_unit_length(vector, n_dims, ierr)
        call assert_equal_int(ierr, ERR_NAN_INF, "test_normalize_unit_length: vector with Infinity should trigger ERR_NAN_INF")
    end subroutine test_normalize_unit_length

    
end module mod_test_normalization_unit_length
