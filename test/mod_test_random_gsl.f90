!> Unit test suite for random_gsl routine.
module mod_test_random_gsl
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_associated, c_ptr
    use f42_random_gsl
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

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(1)
        all_tests(1) = test_case("test", test)
    end function get_all_tests

    !> Run all random_gsl tests.
    subroutine run_all_tests_random_gsl
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All random_gsl tests passed successfully."
    end subroutine run_all_tests_random_gsl

    !> Run specific random_gsl tests by name.
    subroutine run_named_tests_random_gsl(test_names)
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
    end subroutine run_named_tests_random_gsl

    subroutine test
        type(c_ptr) :: rng
        rng = create_rng()
        print *, random_uniform(rng)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_hypergeom(rng, 5, 6, 2)
        print *, random_binomial(rng, 0.6_real64, 5_int32)
        print *, random_binomial(rng, 0.6_real64, 5_int32)
        print *, random_binomial(rng, 0.6_real64, 5_int32)
        print *, random_binomial(rng, 0.6_real64, 5_int32)
        print *, rand_range(0.6_real64, 1.5_real64, rng)
        print *, rand_range(0.6_real64, 1.5_real64, rng)
        print *, rand_range(0.6_real64, 1.5_real64, rng)
        print *, rand_range(0.6_real64, 1.5_real64, rng)
        print *, rand_range(0.6_real64, 1.5_real64, rng)
        call free_rng(rng)
    end subroutine test
end module mod_test_random_gsl
