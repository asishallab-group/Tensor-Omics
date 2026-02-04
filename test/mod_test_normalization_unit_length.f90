!> Unit test suite for normalization_unit_length routine.
module mod_test_normalization_unit_length
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_normalization, only: normalize_unit_length
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use tox_errors
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(1)

        all_tests(1) = test_case("test_normalization_unit_length", test_normalize_unit_length)
    end function get_all_tests

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

    !> Run all normalization_unit_length tests.
    subroutine run_all_tests_normalization_unit_length
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All normalization_unit_length tests passed successfully."
    end subroutine run_all_tests_normalization_unit_length

    !> Run specific normalization_unit_length tests by name.
    subroutine run_named_tests_normalization_unit_length(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_normalization_unit_length
end module mod_test_normalization_unit_length
