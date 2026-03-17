#include "macros.h"

!> Unit test suite for random_gsl routine.
module mod_test_random_gsl
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_ptr
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
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

    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(2)

        all_tests(1) = test_case("test_multiv_hypergeom", test_multiv_hypergeom)
        all_tests(2) = test_case("test_multinomial", test_multinomial)
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

    subroutine test_multiv_hypergeom()
        integer(int32), parameter :: n_pop = 5
        integer(int32) :: pop(n_pop), drawn(n_pop), drawn2(n_pop), expected_pop(n_pop)
        integer(int32) :: total_pop, n_to_draw, i, iter
        type(rng_t) :: rng, rng1, rng2
        real(real64) :: rand

        ! ============================================================
        ! BASIC TEST
        ! ============================================================
        pop = [10, 5, 8, 7, 6]
        expected_pop = pop
        total_pop = sum(pop)
        n_to_draw = 15

        rng = create_rng()

        rand = random_uniform(rng)
        do i = 1, 1000
            call reset_rng(rng)
            call assert_equal_real(random_uniform(rng), rand, 0.0_real64, "test_multiv_hypergeom: resetting rng doesn't result in same values")
        end do
        call random_multiv_hypergeom(rng, pop, n_pop, total_pop, n_to_draw, drawn)

        call assert_equal_int(sum(drawn), n_to_draw, "test_multiv_hypergeom: basic sum(drawn) == n_to_draw")

        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multiv_hypergeom: basic no negative draws")
            ! call assert_true(drawn(i) <= pop(i), "test_multiv_hypergeom: basic cannot exceed population")
        end do

        ! ============================================================
        ! EDGE CASE: ZERO DRAWS
        ! ============================================================
        pop = [5, 7, 9, 3, 1]
        expected_pop = pop
        rng = create_rng()
        call random_multiv_hypergeom(rng, pop, n_pop, sum(pop), 0, drawn)

        call assert_equal_int(sum(drawn), 0, "test_multiv_hypergeom: zero draws sum=0")

        ! ============================================================
        ! EDGE CASE: FULL DRAW
        ! ============================================================
        pop = [3, 4, 5, 6, 7]
        expected_pop = pop
        rng = create_rng()
        call random_multiv_hypergeom(rng, pop, n_pop, sum(pop), sum(pop), drawn)

        call assert_equal_array_int(drawn, expected_pop, n_pop, "test_multiv_hypergeom: full draw returns population")

        ! ============================================================
        ! EDGE CASE: SINGLE POPULATION
        ! ============================================================
        pop = 12
        rng = create_rng()
        call random_multiv_hypergeom(rng, pop, 1, 12, 7, drawn)

        call assert_equal_int(drawn(1), 7, "test_multiv_hypergeom: single population")

        ! ============================================================
        ! REPRODUCIBILITY
        ! ============================================================
        pop = [6, 4, 10, 3, 2]
        expected_pop = pop
        rng1 = create_rng()
        rng2 = create_rng()

        call random_multiv_hypergeom(rng1, pop, n_pop, sum(pop), 8, drawn)
        pop = expected_pop
        call random_multiv_hypergeom(rng2, pop, n_pop, sum(pop), 8, drawn2)

        call assert_equal_array_int(drawn, drawn2, n_pop, "test_multiv_hypergeom: reproducibility")

        ! ============================================================
        ! LARGE POPULATION
        ! ============================================================
        pop = [1000000, 2000000, 1500000, 500000, 250000]
        expected_pop = pop
        total_pop = sum(pop)
        n_to_draw = 1000000

        rng = create_rng()
        call random_multiv_hypergeom(rng, pop, n_pop, total_pop, n_to_draw, drawn)

        call assert_equal_int(sum(drawn), n_to_draw, "test_multiv_hypergeom: large sum correct")

        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multiv_hypergeom: large no negative draws")
            call assert_true(drawn(i) <= pop(i), "test_multiv_hypergeom: large cannot exceed population")
        end do
    end subroutine test_multiv_hypergeom

    subroutine test_multinomial()
        integer(int32), parameter :: n_pop = 4
        integer(int32) :: pop(n_pop), drawn(n_pop), drawn2(n_pop)
        integer(int32) :: total_pop, n_to_draw, i, iter, count_big
        type(rng_t) :: rng, rng1, rng2

        ! ============================================================
        ! BASIC TEST
        ! ============================================================
        pop = [3, 5, 2, 4]
        total_pop = sum(pop)
        n_to_draw = 12

        rng = create_rng()
        call random_multinomial(rng, pop, n_pop, total_pop, n_to_draw, drawn)

        call assert_equal_int(sum(drawn), n_to_draw, "test_multinomial: basic sum(drawn) == n_to_draw")

        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multinomial: basic no negative draws")
        end do

        ! ============================================================
        ! EDGE CASE: ZERO DRAWS
        ! ============================================================
        pop = [1,2,3,4]
        rng = create_rng()
        call random_multinomial(rng, pop, n_pop, sum(pop), 0, drawn)

        call assert_equal_int(sum(drawn), 0, "test_multinomial: zero draws sum=0")

        ! ============================================================
        ! EDGE CASE: SINGLE POPULATION
        ! ============================================================
        pop = 10
        rng = create_rng()
        call random_multinomial(rng, pop, 1, 10, 7, drawn)

        call assert_equal_int(drawn(1), 7, "test_multinomial: single population")

        ! ============================================================
        ! REPRODUCIBILITY
        ! ============================================================
        pop = [4, 3, 3, 2]
        rng1 = create_rng()
        rng2 = create_rng()

        call random_multinomial(rng1, pop, n_pop, sum(pop), 12, drawn)
        call random_multinomial(rng2, pop, n_pop, sum(pop), 12, drawn2)

        call assert_equal_array_int(drawn, drawn2, n_pop, "test_multinomial: reproducibility")

        ! ============================================================
        ! STATISTICAL SANITY CHECK
        ! Category 2 has largest population → should dominate often
        ! ============================================================
        pop = [3, 10, 2, 5]
        count_big = 0
        rng = create_rng()

        do iter = 1, 200
            call random_multinomial(rng, pop, n_pop, sum(pop), 20, drawn)
            if (drawn(2) >= max(drawn(1), drawn(3), drawn(4))) count_big = count_big + 1
        end do

        call assert_true(count_big > 80, "test_multinomial: statistical sanity check")

        ! ============================================================
        ! LARGE POPULATION
        ! ============================================================
        pop = [1000000, 2000000, 3000000, 4000000]
        total_pop = sum(pop)
        n_to_draw = 5000000

        rng = create_rng()
        call random_multinomial(rng, pop, n_pop, total_pop, n_to_draw, drawn)

        call assert_equal_int(sum(drawn), n_to_draw, "test_multinomial: large sum correct")

        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multinomial: large no negative draws")
        end do

    end subroutine test_multinomial

end module mod_test_random_gsl
