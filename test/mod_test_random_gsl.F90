#include <src/macros.h>

!> Unit test suite for f42_random_gsl.
module mod_test_random_gsl
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use f42_random_gsl
    use tox_errors
    use test_suite, only: test_case
    implicit none

contains

    !> Get array of all available tests.
    function get_all_tests_random_gsl() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(19))

        all_tests(1) = test_case("test_create_rng_ok_and_reproducible", test_create_rng_ok_and_reproducible)
        all_tests(2) = test_case("test_reset_rng_restores_stream", test_reset_rng_restores_stream)
        all_tests(3) = test_case("test_uninitialized_rng_reports_pointer_null", test_uninitialized_rng_reports_pointer_null)
        all_tests(4) = test_case("test_destroyed_rng_reports_pointer_null", test_destroyed_rng_reports_pointer_null)
        all_tests(5) = test_case("test_random_uniform_range", test_random_uniform_range)
        all_tests(6) = test_case("test_random_hypergeom_rejects_negative_arguments", &
                                 test_random_hypergeom_rejects_negative_arguments)
        all_tests(7) = test_case("test_random_hypergeom_bounds", test_random_hypergeom_bounds)
        all_tests(8) = test_case("test_random_binomial_rejects_invalid_input", test_random_binomial_rejects_invalid_input)
        all_tests(9) = test_case("test_random_binomial_bounds", test_random_binomial_bounds)
        all_tests(10) = test_case("test_rand_range_bounds_and_rejects_nan", test_rand_range_bounds_and_rejects_nan)
        all_tests(11) = test_case("test_multinomial", test_multinomial)
        all_tests(12) = test_case("test_multiv_hypergeom", test_multiv_hypergeom)
        all_tests(13) = test_case("test_multinomial_invalid_n_populations_arg_pos", &
                                  test_multinomial_invalid_n_populations_arg_pos)
        all_tests(14) = test_case("test_multinomial_invalid_population_sizes_arg_pos", &
                                  test_multinomial_invalid_population_sizes_arg_pos)
        all_tests(15) = test_case("test_multinomial_total_population_mismatch_arg_pos", &
                                  test_multinomial_total_population_mismatch_arg_pos)
        all_tests(16) = test_case("test_multiv_hypergeom_invalid_n_populations_arg_pos", &
                                  test_multiv_hypergeom_invalid_n_populations_arg_pos)
        all_tests(17) = test_case("test_multiv_hypergeom_invalid_population_sizes_arg_pos", &
                                  test_multiv_hypergeom_invalid_population_sizes_arg_pos)
        all_tests(18) = test_case("test_multiv_hypergeom_total_population_mismatch_arg_pos", &
                                  test_multiv_hypergeom_total_population_mismatch_arg_pos)
        all_tests(19) = test_case("test_multiv_hypergeom_pop_sizes_reduction_known_limitation", &
                                  test_multiv_hypergeom_pop_sizes_reduction_known_limitation)
    end function get_all_tests_random_gsl

    !> A generator created with an explicit seed reproduces the exact same first draw as one
    !| created with the same seed independently, and `create_rng` reports success.
    subroutine test_create_rng_ok_and_reproducible()
        type(rng_t) :: rng1, rng2
        integer(int32) :: ierr

        rng1 = create_rng(seed=11_int32, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_create_rng_ok_and_reproducible: default create_rng succeeds")

        rng2 = create_rng(seed=11_int32, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_create_rng_ok_and_reproducible: explicit-seed create_rng succeeds")

        call assert_equal_real(random_uniform(rng1, ierr=ierr), random_uniform(rng2, ierr=ierr), 0.0_real64, &
            "test_create_rng_ok_and_reproducible: two generators seeded alike agree on their first draw")
    end subroutine test_create_rng_ok_and_reproducible

    !> Resetting a generator to the seed it was created with reproduces its first draw, every time.
    subroutine test_reset_rng_restores_stream()
        type(rng_t) :: rng
        integer(int32) :: ierr, i
        real(real64) :: first_draw

        rng = create_rng(seed=7_int32, ierr=ierr)
        first_draw = random_uniform(rng, ierr=ierr)

        do i = 1, 5
            call reset_rng(rng, seed=7_int32, ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_reset_rng_restores_stream: reset_rng succeeds")
            call assert_equal_real(random_uniform(rng, ierr=ierr), first_draw, 0.0_real64, &
                "test_reset_rng_restores_stream: resetting to the same seed reproduces the first draw")
        end do
    end subroutine test_reset_rng_restores_stream

    !> A freshly-declared `rng_t` default-initializes to a null handle, and every procedure that
    !| takes one other than `create_rng` rejects it with `ERR_POINTER_NULL` at argument 1 rather
    !| than dereferencing it.
    subroutine test_uninitialized_rng_reports_pointer_null()
        type(rng_t) :: rng
        integer(int32) :: ierr, n_drawn, pop(2), drawn(2)
        real(real64) :: r

        r = random_uniform(rng, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_uninitialized_rng_reports_pointer_null: random_uniform", &
                        arg_pos=1_int32)
        call assert_equal_real(r, 0.0_real64, 0.0_real64, &
            "test_uninitialized_rng_reports_pointer_null: random_uniform result is zeroed")

        call reset_rng(rng, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_uninitialized_rng_reports_pointer_null: reset_rng", &
                        arg_pos=1_int32)

        n_drawn = random_hypergeom(rng, 5_int32, 5_int32, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_uninitialized_rng_reports_pointer_null: random_hypergeom", &
                        arg_pos=1_int32)

        n_drawn = random_binomial(rng, 0.5_real64, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_uninitialized_rng_reports_pointer_null: random_binomial", &
                        arg_pos=1_int32)

        r = rand_range(rng, 0.0_real64, 1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_uninitialized_rng_reports_pointer_null: rand_range", &
                        arg_pos=1_int32)

        pop = [3, 4]
        call random_multinomial(rng, 2_int32, pop, 7_int32, 5_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_uninitialized_rng_reports_pointer_null: random_multinomial", &
                        arg_pos=1_int32)

        pop = [3, 4]
        call random_multiv_hypergeom(rng, 2_int32, pop, 7_int32, 5_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_uninitialized_rng_reports_pointer_null: random_multiv_hypergeom", &
                        arg_pos=1_int32)
    end subroutine test_uninitialized_rng_reports_pointer_null

    !> A destroyed handle reverts to the same rejected state as one that was never created, and
    !| destroying it a second time must not crash (no double-free of the GSL handle).
    subroutine test_destroyed_rng_reports_pointer_null()
        type(rng_t) :: rng
        integer(int32) :: ierr
        real(real64) :: r

        rng = create_rng(ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_destroyed_rng_reports_pointer_null: create_rng succeeds")

        call destroy_rng(rng)
        r = random_uniform(rng, ierr=ierr)
        call assert_err(ierr, ERR_POINTER_NULL, "test_destroyed_rng_reports_pointer_null: use after destroy", &
                        arg_pos=1_int32)

        call destroy_rng(rng)  ! must be a no-op, not a double-free
    end subroutine test_destroyed_rng_reports_pointer_null

    subroutine test_random_uniform_range()
        type(rng_t) :: rng
        integer(int32) :: ierr, i
        real(real64) :: r

        rng = create_rng(ierr=ierr)
        do i = 1, 200
            r = random_uniform(rng, ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_random_uniform_range: succeeds")
            call assert_true(r >= 0.0_real64 .and. r < 1.0_real64, "test_random_uniform_range: value in [0, 1)")
        end do
    end subroutine test_random_uniform_range

    subroutine test_random_hypergeom_rejects_negative_arguments()
        type(rng_t) :: rng
        integer(int32) :: ierr, n_drawn

        rng = create_rng(ierr=ierr)

        n_drawn = random_hypergeom(rng, -1_int32, 5_int32, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_random_hypergeom_rejects_negative_arguments: n_population1 < 0", &
                        arg_pos=2_int32)
        call assert_equal_int(n_drawn, 0_int32, "test_random_hypergeom_rejects_negative_arguments: result is zeroed")

        n_drawn = random_hypergeom(rng, 5_int32, -1_int32, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_random_hypergeom_rejects_negative_arguments: n_population2 < 0", &
                        arg_pos=3_int32)

        n_drawn = random_hypergeom(rng, 5_int32, 5_int32, -1_int32, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_random_hypergeom_rejects_negative_arguments: n_samples < 0", &
                        arg_pos=4_int32)

        n_drawn = random_hypergeom(rng, 5_int32, 5_int32, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_random_hypergeom_rejects_negative_arguments: valid input succeeds")
    end subroutine test_random_hypergeom_rejects_negative_arguments

    subroutine test_random_hypergeom_bounds()
        type(rng_t) :: rng
        integer(int32) :: ierr, n_drawn, i

        rng = create_rng(ierr=ierr)
        do i = 1, 100
            n_drawn = random_hypergeom(rng, 10_int32, 5_int32, 8_int32, ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_random_hypergeom_bounds: succeeds")
            call assert_true(n_drawn >= max(0_int32, 8_int32 - 5_int32), &
                "test_random_hypergeom_bounds: at least max(0, n_samples - n_population2)")
            call assert_true(n_drawn <= min(8_int32, 10_int32), &
                "test_random_hypergeom_bounds: at most min(n_samples, n_population1)")
        end do
    end subroutine test_random_hypergeom_bounds

    subroutine test_random_binomial_rejects_invalid_input()
        type(rng_t) :: rng
        integer(int32) :: ierr, n_drawn

        rng = create_rng(ierr=ierr)

        n_drawn = random_binomial(rng, -0.1_real64, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_random_binomial_rejects_invalid_input: p < 0", arg_pos=2_int32)
        call assert_equal_int(n_drawn, 0_int32, "test_random_binomial_rejects_invalid_input: result is zeroed")

        n_drawn = random_binomial(rng, 1.1_real64, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_random_binomial_rejects_invalid_input: p > 1", arg_pos=2_int32)

        n_drawn = random_binomial(rng, 0.5_real64, -1_int32, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_random_binomial_rejects_invalid_input: n_samples < 0", &
                        arg_pos=3_int32)

        n_drawn = random_binomial(rng, 0.5_real64, 3_int32, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_random_binomial_rejects_invalid_input: valid input succeeds")
    end subroutine test_random_binomial_rejects_invalid_input

    subroutine test_random_binomial_bounds()
        type(rng_t) :: rng
        integer(int32) :: ierr, n_drawn, i

        rng = create_rng(ierr=ierr)
        do i = 1, 100
            n_drawn = random_binomial(rng, 0.5_real64, 20_int32, ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_random_binomial_bounds: succeeds")
            call assert_true(n_drawn >= 0_int32 .and. n_drawn <= 20_int32, &
                "test_random_binomial_bounds: value in [0, n_samples]")
        end do
    end subroutine test_random_binomial_bounds

    subroutine test_rand_range_bounds_and_rejects_nan()
        type(rng_t) :: rng
        integer(int32) :: ierr, i
        real(real64) :: r, nan

        nan = ieee_value(1.0_real64, ieee_quiet_nan)
        rng = create_rng(ierr=ierr)

        do i = 1, 100
            r = rand_range(rng, 2.0_real64, 5.0_real64, ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_rand_range_bounds_and_rejects_nan: succeeds")
            call assert_true(r >= 2.0_real64 .and. r < 5.0_real64, "test_rand_range_bounds_and_rejects_nan: in [min, max)")
        end do

        r = rand_range(rng, 3.0_real64, 3.0_real64, ierr=ierr)
        call assert_equal_real(r, 3.0_real64, 0.0_real64, "test_rand_range_bounds_and_rejects_nan: min == max returns min")

        r = rand_range(rng, nan, 5.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "test_rand_range_bounds_and_rejects_nan: NaN min is rejected", arg_pos=2_int32)

        r = rand_range(rng, 2.0_real64, nan, ierr=ierr)
        call assert_err(ierr, ERR_NAN_INF, "test_rand_range_bounds_and_rejects_nan: NaN max is rejected", arg_pos=3_int32)
    end subroutine test_rand_range_bounds_and_rejects_nan

    subroutine test_multinomial()
        integer(int32), parameter :: n_pop = 4
        integer(int32) :: pop(n_pop), drawn(n_pop), drawn2(n_pop)
        integer(int32) :: total_pop, n_to_draw, i, iter, count_big, ierr
        type(rng_t) :: rng, rng1, rng2

        ! ============================================================
        ! BASIC TEST
        ! ============================================================
        pop = [3, 5, 2, 4]
        total_pop = sum(pop)
        n_to_draw = 12

        rng = create_rng(ierr=ierr)
        call random_multinomial(rng, n_pop, pop, total_pop, n_to_draw, drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multinomial: basic call succeeds")
        call assert_equal_int(sum(drawn), n_to_draw, "test_multinomial: basic sum(drawn) == n_to_draw")
        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multinomial: basic no negative draws")
        end do

        ! ============================================================
        ! EDGE CASE: ZERO DRAWS
        ! ============================================================
        pop = [1, 2, 3, 4]
        call random_multinomial(rng, n_pop, pop, sum(pop), 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multinomial: zero draws succeeds")
        call assert_equal_int(sum(drawn), 0, "test_multinomial: zero draws sum == 0")

        ! ============================================================
        ! EDGE CASE: FULL DRAW
        ! ============================================================
        ! Unlike random_multiv_hypergeom (sampling *without* replacement, so drawing every
        ! remaining element forces an exact match), random_multinomial draws *with* replacement
        ! via a sequence of independent binomial splits -- n_to_draw == total_population does
        ! not make the outcome deterministic, only sum(drawn) == n_to_draw and non-negativity.
        pop = [4, 5, 6, 7]
        call random_multinomial(rng, n_pop, pop, sum(pop), sum(pop), drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multinomial: full draw succeeds")
        call assert_equal_int(sum(drawn), sum(pop), "test_multinomial: full draw: sum(drawn) == n_to_draw")
        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multinomial: full draw: no negative draws")
        end do

        ! ============================================================
        ! EDGE CASE: SINGLE POPULATION
        ! ============================================================
        call random_multinomial(rng, 1_int32, [10_int32], 10_int32, 7_int32, drawn(1:1), ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multinomial: single population succeeds")
        call assert_equal_int(drawn(1), 7, "test_multinomial: single population gets everything")

        ! ============================================================
        ! REPRODUCIBILITY
        ! ============================================================
        pop = [4, 3, 3, 2]
        rng1 = create_rng(seed=17_int32, ierr=ierr)
        rng2 = create_rng(seed=17_int32, ierr=ierr)

        call random_multinomial(rng1, n_pop, pop, sum(pop), 12_int32, drawn, ierr=ierr)
        call random_multinomial(rng2, n_pop, pop, sum(pop), 12_int32, drawn2, ierr=ierr)
        call assert_equal_array_int(drawn, drawn2, n_pop, "test_multinomial: reproducibility with the same seed")

        ! ============================================================
        ! STATISTICAL SANITY CHECK: the largest subpopulation should dominate often
        ! ============================================================
        pop = [3, 10, 2, 5]
        count_big = 0
        do iter = 1, 200
            call random_multinomial(rng, n_pop, pop, sum(pop), 20_int32, drawn, ierr=ierr)
            if (drawn(2) >= max(drawn(1), drawn(3), drawn(4))) count_big = count_big + 1
        end do
        call assert_true(count_big > 80, "test_multinomial: statistical sanity check")

        ! ============================================================
        ! LARGE POPULATION
        ! ============================================================
        pop = [1000000, 2000000, 3000000, 4000000]
        total_pop = sum(pop)
        n_to_draw = 5000000
        call random_multinomial(rng, n_pop, pop, total_pop, n_to_draw, drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multinomial: large population succeeds")
        call assert_equal_int(sum(drawn), n_to_draw, "test_multinomial: large sum correct")
        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multinomial: large no negative draws")
        end do
    end subroutine test_multinomial

    subroutine test_multiv_hypergeom()
        integer(int32), parameter :: n_pop = 5
        integer(int32) :: pop(n_pop), drawn(n_pop), drawn2(n_pop), original_pop(n_pop)
        integer(int32) :: total_pop, n_to_draw, i, ierr
        type(rng_t) :: rng, rng1, rng2

        ! ============================================================
        ! BASIC TEST
        ! ============================================================
        pop = [10, 5, 8, 7, 6]
        total_pop = sum(pop)
        n_to_draw = 15

        rng = create_rng(ierr=ierr)
        call random_multiv_hypergeom(rng, n_pop, pop, total_pop, n_to_draw, drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multiv_hypergeom: basic call succeeds")
        call assert_equal_int(sum(drawn), n_to_draw, "test_multiv_hypergeom: basic sum(drawn) == n_to_draw")
        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multiv_hypergeom: basic no negative draws")
        end do

        ! ============================================================
        ! EDGE CASE: ZERO DRAWS
        ! ============================================================
        pop = [5, 7, 9, 3, 1]
        call random_multiv_hypergeom(rng, n_pop, pop, sum(pop), 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multiv_hypergeom: zero draws succeeds")
        call assert_equal_int(sum(drawn), 0, "test_multiv_hypergeom: zero draws sum == 0")

        ! ============================================================
        ! EDGE CASE: FULL DRAW
        ! ============================================================
        pop = [3, 4, 5, 6, 7]
        original_pop = pop
        call random_multiv_hypergeom(rng, n_pop, pop, sum(pop), sum(pop), drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multiv_hypergeom: full draw succeeds")
        call assert_equal_array_int(drawn, original_pop, n_pop, "test_multiv_hypergeom: full draw returns the population")

        ! ============================================================
        ! EDGE CASE: SINGLE POPULATION
        ! ============================================================
        pop(1:1) = [12]
        call random_multiv_hypergeom(rng, 1_int32, pop(1:1), 12_int32, 7_int32, drawn(1:1), ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multiv_hypergeom: single population succeeds")
        call assert_equal_int(drawn(1), 7, "test_multiv_hypergeom: single population gets everything")

        ! ============================================================
        ! REPRODUCIBILITY
        ! ============================================================
        pop = [6, 4, 10, 3, 2]
        original_pop = pop
        rng1 = create_rng(seed=99_int32, ierr=ierr)
        rng2 = create_rng(seed=99_int32, ierr=ierr)

        call random_multiv_hypergeom(rng1, n_pop, pop, sum(pop), 8_int32, drawn, ierr=ierr)
        pop = original_pop
        call random_multiv_hypergeom(rng2, n_pop, pop, sum(pop), 8_int32, drawn2, ierr=ierr)
        call assert_equal_array_int(drawn, drawn2, n_pop, "test_multiv_hypergeom: reproducibility with the same seed")

        ! ============================================================
        ! LARGE POPULATION
        ! ============================================================
        pop = [1000000, 2000000, 1500000, 500000, 250000]
        original_pop = pop
        total_pop = sum(pop)
        n_to_draw = 1000000
        call random_multiv_hypergeom(rng, n_pop, pop, total_pop, n_to_draw, drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multiv_hypergeom: large population succeeds")
        call assert_equal_int(sum(drawn), n_to_draw, "test_multiv_hypergeom: large sum correct")
        do i = 1, n_pop
            call assert_true(drawn(i) >= 0, "test_multiv_hypergeom: large no negative draws")
            call assert_true(drawn(i) <= original_pop(i), "test_multiv_hypergeom: large cannot exceed subpopulation")
        end do
    end subroutine test_multiv_hypergeom

    !> Regression test for a bug found (not introduced) during this port, faithfully preserved
    !| from the ported algorithm rather than fixed, per this stage's "port near-verbatim" rule:
    !| `population_sizes` is documented as the remaining pool after the draw, and elements
    !| `1..n_populations-1` really are reduced by `drawn(i)` -- but the internal loop only ever
    !| runs over `1..n_populations-1` and mutates `population_sizes` through an alias scoped to
    !| that range, so the *last* element, `population_sizes(n_populations)`, is never written
    !| back and comes out unchanged, however much of it was drawn.
    !| Known limitation: see [[f42_random_gsl(module):random_multiv_hypergeom(subroutine)]]'s own
    !| doc comment.
    subroutine test_multiv_hypergeom_pop_sizes_reduction_known_limitation()
        integer(int32), parameter :: n_pop = 4
        integer(int32) :: pop(n_pop), original_pop(n_pop), drawn(n_pop), ierr, i
        type(rng_t) :: rng

        pop = [10, 20, 15, 5]
        original_pop = pop
        rng = create_rng(ierr=ierr)

        call random_multiv_hypergeom(rng, n_pop, pop, sum(original_pop), 25_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_multiv_hypergeom_pop_sizes_reduction_known_limitation: succeeds")

        do i = 1, n_pop - 1
            call assert_equal_int(pop(i), original_pop(i) - drawn(i), &
                "test_multiv_hypergeom_pop_sizes_reduction_known_limitation: "// &
                "population_sizes(i) is reduced by drawn(i) for every element but the last")
        end do

        ! Guard against a vacuous pass: the limitation is only visible if the last
        ! subpopulation actually had something drawn from it.
        call assert_true(drawn(n_pop) > 0_int32, &
            "test_multiv_hypergeom_pop_sizes_reduction_known_limitation: "// &
            "fixture draws a nonzero amount from the last subpopulation")

        ! The known limitation itself: the last element is returned unchanged, not reduced.
        call assert_equal_int(pop(n_pop), original_pop(n_pop), &
            "test_multiv_hypergeom_pop_sizes_reduction_known_limitation: "// &
            "population_sizes(n_populations) is NOT reduced by drawn(n_populations) -- known limitation")
    end subroutine test_multiv_hypergeom_pop_sizes_reduction_known_limitation

    !> Regression test for a bug caught during this port: `random_multinomial`'s dummy list had
    !| to be reordered (`n_populations` before `population_sizes`) so its declaration precedes
    !| the array it sizes, and the very first attempt at that reorder left the `arg_pos=`
    !| literals pointing at the *old* positions -- `n_populations`'s own bad-value check blamed
    !| argument 3 (`population_sizes`'s new position) instead of its own, argument 2. Pins the
    !| correct position so a future reordering cannot silently swap them back unnoticed.
    subroutine test_multinomial_invalid_n_populations_arg_pos()
        integer(int32) :: drawn(0), pop(0), ierr
        type(rng_t) :: rng

        rng = create_rng(ierr=ierr)

        call random_multinomial(rng, 0_int32, pop, 0_int32, 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, &
            "test_multinomial_invalid_n_populations_arg_pos: n_populations == 0 blames argument 2", arg_pos=2_int32)

        call random_multinomial(rng, -1_int32, pop, 0_int32, 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
            "test_multinomial_invalid_n_populations_arg_pos: n_populations < 0 blames argument 2", arg_pos=2_int32)
    end subroutine test_multinomial_invalid_n_populations_arg_pos

    !> Companion regression test: a negative element of `population_sizes` (`n_populations`
    !| itself valid) must blame argument 3, `population_sizes`'s own, current position.
    subroutine test_multinomial_invalid_population_sizes_arg_pos()
        integer(int32), parameter :: n_pop = 3
        integer(int32) :: pop(n_pop), drawn(n_pop), ierr
        type(rng_t) :: rng

        rng = create_rng(ierr=ierr)
        pop = [1, -2, 3]

        call random_multinomial(rng, n_pop, pop, sum(pop), 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
            "test_multinomial_invalid_population_sizes_arg_pos: negative element blames argument 3", arg_pos=3_int32)
        call assert_equal_int(sum(drawn), 0_int32, "test_multinomial_invalid_population_sizes_arg_pos: drawn is zeroed")
    end subroutine test_multinomial_invalid_population_sizes_arg_pos

    !> `total_population` that does not equal `sum(population_sizes)` is a dimension mismatch,
    !| not silently trusted -- the draw itself would otherwise read past the pool with no
    !| diagnostic. Blames `total_population`'s own position (4), the derived cross-check value,
    !| in the spirit of codegen_guide.md §5.3's `n_selected_<arg>` convention.
    subroutine test_multinomial_total_population_mismatch_arg_pos()
        integer(int32), parameter :: n_pop = 3
        integer(int32) :: pop(n_pop), drawn(n_pop), ierr
        type(rng_t) :: rng

        rng = create_rng(ierr=ierr)
        pop = [1, 2, 3]

        call random_multinomial(rng, n_pop, pop, sum(pop) + 1_int32, 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
            "test_multinomial_total_population_mismatch_arg_pos: total_population /= sum(population_sizes)", &
            arg_pos=4_int32)
    end subroutine test_multinomial_total_population_mismatch_arg_pos

    !> Same regression as `test_multinomial_invalid_n_populations_arg_pos`, for the multivariate
    !| hypergeometric sibling.
    subroutine test_multiv_hypergeom_invalid_n_populations_arg_pos()
        integer(int32) :: drawn(0), pop(0), ierr
        type(rng_t) :: rng

        rng = create_rng(ierr=ierr)

        call random_multiv_hypergeom(rng, 0_int32, pop, 0_int32, 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, &
            "test_multiv_hypergeom_invalid_n_populations_arg_pos: n_populations == 0 blames argument 2", arg_pos=2_int32)
    end subroutine test_multiv_hypergeom_invalid_n_populations_arg_pos

    !> Same regression as `test_multinomial_invalid_population_sizes_arg_pos`, for the
    !| multivariate hypergeometric sibling.
    subroutine test_multiv_hypergeom_invalid_population_sizes_arg_pos()
        integer(int32), parameter :: n_pop = 3
        integer(int32) :: pop(n_pop), drawn(n_pop), ierr
        type(rng_t) :: rng

        rng = create_rng(ierr=ierr)
        pop = [1, -2, 3]

        call random_multiv_hypergeom(rng, n_pop, pop, sum(pop), 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
            "test_multiv_hypergeom_invalid_population_sizes_arg_pos: negative element blames argument 3", &
            arg_pos=3_int32)
    end subroutine test_multiv_hypergeom_invalid_population_sizes_arg_pos

    !> Same regression as `test_multinomial_total_population_mismatch_arg_pos`, for the
    !| multivariate hypergeometric sibling.
    subroutine test_multiv_hypergeom_total_population_mismatch_arg_pos()
        integer(int32), parameter :: n_pop = 3
        integer(int32) :: pop(n_pop), drawn(n_pop), ierr
        type(rng_t) :: rng

        rng = create_rng(ierr=ierr)
        pop = [1, 2, 3]

        call random_multiv_hypergeom(rng, n_pop, pop, sum(pop) + 1_int32, 0_int32, drawn, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
            "test_multiv_hypergeom_total_population_mismatch_arg_pos: total_population /= sum(population_sizes)", &
            arg_pos=4_int32)
    end subroutine test_multiv_hypergeom_total_population_mismatch_arg_pos

end module mod_test_random_gsl
