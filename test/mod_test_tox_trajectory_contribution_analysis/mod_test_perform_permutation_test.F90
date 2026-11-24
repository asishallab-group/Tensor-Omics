!> The `perform_permutation_test` cases: each permutation pairs the factor of `sample_idx` with
!| the dependent of a randomly drawn *other* sample, so the cases check which sample was drawn,
!| never which random number came out.
!|
!| Most cases make the draw visible in the result: the factor is [1, 0] in every sample and the
!| dependent of sample s is [s, 0], so with MODE_BASELINE_RAW a permutation's contributions are
!| [s, 0] and its total is s, the sample it drew. Every case that draws seeds the generator, so
!| it is deterministic; see `test_perform_permutation_test_random_seed` for why that is still a
!| test of the randomness.
module mod_test_perform_permutation_test
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: perform_permutation_test, perform_permutation_test_expert
    use tox_trajectory_contribution_analysis_impl, only: MODE_BASELINE_RAW, MODE_BASELINE_MEAN
    use tox_errors
    use f42_random_impl, only: init_random
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_perform_permutation_test

    ! The layout of `fill_sample_marking`: factor 1 is the factor, factor 2 the dependent.
    integer(int32), parameter :: MARK_N_FACTORS = 2, MARK_N_TIMEPOINTS = 2
    integer(int32), parameter :: MARK_FACTOR_IDX = 1, MARK_DEPENDENT_IDX = 2

contains

    !> Every case of `perform_permutation_test`.
    function get_all_tests_perform_permutation_test() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(10))
        all_tests(1) = test_case("test_perform_permutation_test_values", test_perform_permutation_test_values)
        all_tests(2) = test_case("test_perform_permutation_test_expert", test_perform_permutation_test_expert)
        all_tests(3) = test_case("test_perform_permutation_test_two_samples", &
                                 test_perform_permutation_test_two_samples)
        all_tests(4) = test_case("test_perform_permutation_test_draws_every_other_sample", &
                                 test_perform_permutation_test_draws_every_other_sample)
        all_tests(5) = test_case("test_perform_permutation_test_random_seed", &
                                 test_perform_permutation_test_random_seed)
        all_tests(6) = test_case("test_perform_permutation_test_dimensions", test_perform_permutation_test_dimensions)
        all_tests(7) = test_case("test_perform_permutation_test_index_bounds", &
                                 test_perform_permutation_test_index_bounds)
        all_tests(8) = test_case("test_perform_permutation_test_invalid_mode", &
                                 test_perform_permutation_test_invalid_mode)
        all_tests(9) = test_case("test_perform_permutation_test_rejects_nan_and_inf", &
                                 test_perform_permutation_test_rejects_nan_and_inf)
        all_tests(10) = test_case("test_perform_permutation_test_single_sample", &
                                  test_perform_permutation_test_single_sample)
    end function get_all_tests_perform_permutation_test

    !> The sample-marking data: in each of `n_samples` samples the factor is [1, 0] and the
    !| dependent of sample s is [s, 0].
    subroutine fill_sample_marking(n_samples, trajectories)
        integer(int32), intent(in) :: n_samples
        real(real64), intent(out) :: trajectories(MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS)
        integer(int32) :: i_sample

        do i_sample = 1, n_samples
            trajectories(MARK_FACTOR_IDX, i_sample, :) = [1.0_real64, 0.0_real64]
            trajectories(MARK_DEPENDENT_IDX, i_sample, 1) = real(i_sample, real64)
            trajectories(MARK_DEPENDENT_IDX, i_sample, 2) = 0.0_real64
        end do
    end subroutine fill_sample_marking

    !> Runs the plain tier, raw, on the sample-marking data of `n_samples` samples, and returns
    !| the sample each permutation drew. Checks `ierr` OK and that every permutation's
    !| contributions are [s, 0] with total s for some sample s, so the draw is read correctly.
    subroutine draw_samples(n_samples, sample_idx, n_permutations, random_seed, drawn, msg)
        integer(int32), intent(in) :: n_samples, sample_idx, n_permutations
        integer(int32), intent(in), optional :: random_seed
        integer(int32), intent(out) :: drawn(n_permutations)
        character(*), intent(in) :: msg
        real(real64) :: trajectories(MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS)
        real(real64) :: local_contributions(MARK_N_TIMEPOINTS, n_permutations), total_contributions(n_permutations)
        integer(int32) :: ierr, i_perm

        call fill_sample_marking(n_samples, trajectories)
        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN
        call perform_permutation_test(trajectories, MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS, MARK_FACTOR_IDX, &
                                      MARK_DEPENDENT_IDX, sample_idx, MODE_BASELINE_RAW, n_permutations, &
                                      local_contributions, total_contributions, random_seed, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")

        do i_perm = 1, n_permutations
            drawn(i_perm) = nint(total_contributions(i_perm), int32)
            call assert_equal_real(total_contributions(i_perm), real(drawn(i_perm), real64), 0.0_real64, &
                                   msg//": a total is a sample number")
            call assert_equal_real(local_contributions(1, i_perm), real(drawn(i_perm), real64), 0.0_real64, &
                                   msg//": the first contribution is the total")
            call assert_equal_real(local_contributions(2, i_perm), 0.0_real64, 0.0_real64, &
                                   msg//": the second contribution is 0")
        end do
    end subroutine draw_samples

    !> The values of a permutation, kept from the old suite and made to tell the samples apart:
    !| four samples, the factor (factor 1) of sample 2 is [1, 2, 3, 6] and the dependent
    !| (factor 2) is [8, 2, 4, 10] in every sample but sample 2. With MODE_BASELINE_MEAN the
    !| deviations are [-2, -1, 0, 3] (mean 3) and [2, -4, -2, 4] (mean 6), so every permutation,
    !| whichever other sample it draws, is [-4, 4, 0, 12] with total 12. Sample 2's own dependent
    !| [0, 0, 0, 100] would give [-50, -25, 0, 75], total 0, and the factor of any other sample,
    !| [5, 5, 5, 5], would give zeros: drawing the own sample or the wrong factor shows.
    subroutine fill_values(trajectories, expected_local)
        real(real64), intent(out) :: trajectories(2, 4, 4), expected_local(4)
        integer(int32) :: i_sample

        do i_sample = 1, 4
            trajectories(1, i_sample, :) = 5.0_real64
            trajectories(2, i_sample, :) = [8.0_real64, 2.0_real64, 4.0_real64, 10.0_real64]
        end do
        trajectories(1, 2, :) = [1.0_real64, 2.0_real64, 3.0_real64, 6.0_real64]
        trajectories(2, 2, :) = [0.0_real64, 0.0_real64, 0.0_real64, 100.0_real64]
        expected_local = [-4.0_real64, 4.0_real64, 0.0_real64, 12.0_real64]
    end subroutine fill_values

    !> Five permutations of the data of `fill_values`, without a seed: whatever the generator
    !| draws, each is [-4, 4, 0, 12] with total 12.
    subroutine test_perform_permutation_test_values()
        integer(int32), parameter :: n_permutations = 5
        real(real64) :: trajectories(2, 4, 4), expected_local(4), expected_locals(4, n_permutations)
        real(real64) :: expected_totals(n_permutations)
        real(real64) :: local_contributions(4, n_permutations), total_contributions(n_permutations)
        integer(int32) :: ierr, i_perm

        call fill_values(trajectories, expected_local)
        do i_perm = 1, n_permutations
            expected_locals(:, i_perm) = expected_local
        end do
        expected_totals = 12.0_real64
        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN

        call perform_permutation_test(trajectories, 2, 4, 4, 1, 2, 2, MODE_BASELINE_MEAN, n_permutations, &
                                      local_contributions, total_contributions, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_perform_permutation_test_values: ierr")
        call assert_equal_array_real(local_contributions, expected_locals, 4*n_permutations, 0.0_real64, &
                                     "test_perform_permutation_test_values: local_contributions", 4)
        call assert_equal_array_real(total_contributions, expected_totals, n_permutations, 0.0_real64, &
                                     "test_perform_permutation_test_values: total_contributions")
    end subroutine test_perform_permutation_test_values

    !> The `_expert` tier on the data of `fill_values`, seeded, with work arrays that start out
    !| full of a stale value: the same [-4, 4, 0, 12], total 12, for every permutation.
    subroutine test_perform_permutation_test_expert()
        integer(int32), parameter :: n_permutations = 5
        real(real64) :: trajectories(2, 4, 4), expected_local(4), expected_locals(4, n_permutations)
        real(real64) :: expected_totals(n_permutations), tmp_factor(4), tmp_dependent(4)
        real(real64) :: local_contributions(4, n_permutations), total_contributions(n_permutations)
        integer(int32) :: ierr, i_perm

        call fill_values(trajectories, expected_local)
        do i_perm = 1, n_permutations
            expected_locals(:, i_perm) = expected_local
        end do
        expected_totals = 12.0_real64
        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN
        tmp_factor = 1.0e300_real64
        tmp_dependent = -1.0e300_real64

        call perform_permutation_test_expert(trajectories, 2, 4, 4, 1, 2, 2, MODE_BASELINE_MEAN, n_permutations, &
                                             local_contributions, total_contributions, tmp_factor, tmp_dependent, &
                                             random_seed=42_int32, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_perform_permutation_test_expert: ierr")
        call assert_equal_array_real(local_contributions, expected_locals, 4*n_permutations, 0.0_real64, &
                                     "test_perform_permutation_test_expert: local_contributions", 4)
        call assert_equal_array_real(total_contributions, expected_totals, n_permutations, 0.0_real64, &
                                     "test_perform_permutation_test_expert: total_contributions")
    end subroutine test_perform_permutation_test_expert

    !> n_samples = 2, the smallest that has another sample: the draw is forced. With sample 1
    !| every permutation draws sample 2, with sample 2 every one draws sample 1 (kept from the old
    !| `select_random_sample` case).
    subroutine test_perform_permutation_test_two_samples()
        integer(int32), parameter :: n_permutations = 8
        integer(int32) :: drawn(n_permutations), expected(n_permutations)

        call draw_samples(2, 1, n_permutations, 42_int32, drawn, "test_perform_permutation_test_two_samples: sample 1")
        expected = 2
        call assert_equal_array_int(drawn, expected, n_permutations, &
                                    "test_perform_permutation_test_two_samples: sample 1 draws sample 2")

        call draw_samples(2, 2, n_permutations, 42_int32, drawn, "test_perform_permutation_test_two_samples: sample 2")
        expected = 1
        call assert_equal_array_int(drawn, expected, n_permutations, &
                                    "test_perform_permutation_test_two_samples: sample 2 draws sample 1")
    end subroutine test_perform_permutation_test_two_samples

    !> Five samples, 200 seeded permutations each for sample 1, 3 and 5 (both ends of the range,
    !| and one in the middle, where the draw has to step over the own sample): no permutation
    !| draws the own sample, and every other sample is drawn (kept from the old
    !| `select_random_sample` case). The seed fixes the draws; that a uniform draw misses one of
    !| four samples 200 times has the probability 4 * (3/4)**200, about 4e-25, so no reasonable
    !| generator fails this by chance.
    subroutine test_perform_permutation_test_draws_every_other_sample()
        integer(int32), parameter :: n_samples = 5, n_permutations = 200
        integer(int32) :: drawn(n_permutations), sample_indices(3), sample_idx, i_case, i_sample
        character(len=2) :: label

        sample_indices = [1, 3, 5]
        do i_case = 1, size(sample_indices)
            sample_idx = sample_indices(i_case)
            write (label, "(I0)") sample_idx
            call draw_samples(n_samples, sample_idx, n_permutations, 42_int32, drawn, &
                              "test_perform_permutation_test_draws_every_other_sample: sample "//trim(label))
            do i_sample = 1, n_samples
                if (i_sample == sample_idx) then
                    call assert_equal_int(count(drawn == i_sample, kind=int32), 0_int32, &
                                          "test_perform_permutation_test_draws_every_other_sample: sample " &
                                          //trim(label)//" never draws itself")
                else
                    call assert_true(any(drawn == i_sample), &
                                     "test_perform_permutation_test_draws_every_other_sample: sample " &
                                     //trim(label)//" draws every other sample")
                end if
            end do
        end do
    end subroutine test_perform_permutation_test_draws_every_other_sample

    !> What `random_seed` promises, each part deterministic. 64 permutations of sample 3 out of
    !| five samples are a sequence of draws from the four others; with the seed 42 as reference:
    !| - the same seed gives the same draws, even after the caller's generator has moved on: the
    !|   procedure reseeds it;
    !| - without `random_seed`, the draws continue the caller's generator, so after
    !|   `init_random(42)` they are the seeded ones;
    !| - another seed (43) gives other draws: a seed that changed nothing would be no seed. Two
    !|   uniform sequences of 64 draws from four samples coincide with the probability 4**-64,
    !|   about 3e-39.
    !| This replaces the old "Case 2 should not be reproducible", which compared two unseeded runs
    !| of two permutations out of nineteen samples and failed whenever both drew the same, with
    !| the probability (1/19)**2, about 0.28%.
    subroutine test_perform_permutation_test_random_seed()
        integer(int32), parameter :: n_samples = 5, sample_idx = 3, n_permutations = 64
        integer(int32) :: reference(n_permutations), drawn(n_permutations)
        real(real64) :: unrelated_draws(7)

        call draw_samples(n_samples, sample_idx, n_permutations, 42_int32, reference, &
                          "test_perform_permutation_test_random_seed: seed 42")

        call init_random(7_int32)
        call random_number(unrelated_draws)
        call draw_samples(n_samples, sample_idx, n_permutations, 42_int32, drawn, &
                          "test_perform_permutation_test_random_seed: seed 42 again")
        call assert_equal_array_int(drawn, reference, n_permutations, &
                                    "test_perform_permutation_test_random_seed: the same seed draws the same")

        call init_random(42_int32)
        call draw_samples(n_samples, sample_idx, n_permutations, drawn=drawn, &
                          msg="test_perform_permutation_test_random_seed: no seed")
        call assert_equal_array_int(drawn, reference, n_permutations, &
                                    "test_perform_permutation_test_random_seed: no seed continues init_random(42)")

        call draw_samples(n_samples, sample_idx, n_permutations, 43_int32, drawn, &
                          "test_perform_permutation_test_random_seed: seed 43")
        call assert_true(any(drawn /= reference), &
                         "test_perform_permutation_test_random_seed: another seed draws otherwise")
    end subroutine test_perform_permutation_test_random_seed

    !> Each dimension on its own, the others valid: `n_factors` (argument 2), `n_timepoints` (4)
    !| and `n_permutations` (9) are ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1. `n_samples`
    !| (3) has the minimum 2, so 0 and -1 are both ERR_INVALID_INPUT; 1 is in
    !| `test_perform_permutation_test_single_sample`.
    subroutine test_perform_permutation_test_dimensions()
        integer(int32), parameter :: n_samples = 3, n_permutations = 2
        integer(int32), parameter :: SIZE_VALUES(2) = [0, -1]
        integer(int32), parameter :: EXPECTED_CODES(2) = [ERR_EMPTY_INPUT, ERR_INVALID_INPUT]
        character(*), parameter :: DIM_NAMES(4) = ["n_factors     ", "n_samples     ", "n_timepoints  ", &
                                                   "n_permutations"]
        real(real64) :: trajectories(MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS)
        real(real64) :: local_contributions(MARK_N_TIMEPOINTS, n_permutations), total_contributions(n_permutations)
        integer(int32) :: dims(4), positions(4), i_dim, i_value, expected_code, ierr

        call fill_sample_marking(n_samples, trajectories)
        positions = [2, 3, 4, 9]
        do i_dim = 1, 4
            do i_value = 1, 2
                dims = [MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS, n_permutations]
                dims(i_dim) = SIZE_VALUES(i_value)
                call perform_permutation_test(trajectories, dims(1), dims(2), dims(3), 1, 1, 1, MODE_BASELINE_RAW, &
                                              dims(4), local_contributions, total_contributions, 42_int32, ierr)
                expected_code = merge(ERR_INVALID_INPUT, EXPECTED_CODES(i_value), i_dim == 2)
                call assert_err(ierr, expected_code, "test_perform_permutation_test_dimensions: " &
                                //trim(DIM_NAMES(i_dim))//" = "//merge(" 0", "-1", i_value == 1), positions(i_dim))
            end do
        end do
    end subroutine test_perform_permutation_test_dimensions

    !> `factor_idx` (argument 5) and `dependent_idx` (6) lie in [1, n_factors], `sample_idx` (7)
    !| in [1, n_samples]: 0 and the upper bound + 1, just outside, are ERR_INVALID_INPUT. The
    !| bounds themselves are accepted in `test_perform_permutation_test_draws_every_other_sample`
    !| (samples 1 and 5 of five) and wherever the factor (1) and the dependent (2 of two) are
    !| taken.
    subroutine test_perform_permutation_test_index_bounds()
        integer(int32), parameter :: n_samples = 3, n_permutations = 2
        real(real64) :: trajectories(MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS)
        real(real64) :: local_contributions(MARK_N_TIMEPOINTS, n_permutations), total_contributions(n_permutations)
        integer(int32) :: indices(3), upper(3), positions(3), i_index, i_side, ierr
        character(*), parameter :: INDEX_NAMES(3) = ["factor_idx   ", "dependent_idx", "sample_idx   "]

        call fill_sample_marking(n_samples, trajectories)
        upper = [MARK_N_FACTORS, MARK_N_FACTORS, n_samples]
        positions = [5, 6, 7]
        do i_index = 1, 3
            do i_side = 1, 2
                indices = [MARK_FACTOR_IDX, MARK_DEPENDENT_IDX, 1]
                indices(i_index) = merge(0_int32, upper(i_index) + 1_int32, i_side == 1)
                call perform_permutation_test(trajectories, MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS, &
                                              indices(1), indices(2), indices(3), MODE_BASELINE_RAW, n_permutations, &
                                              local_contributions, total_contributions, 42_int32, ierr)
                call assert_err(ierr, ERR_INVALID_INPUT, "test_perform_permutation_test_index_bounds: " &
                                //trim(INDEX_NAMES(i_index))//merge(" = 0        ", " = upper + 1", i_side == 1), &
                                positions(i_index))
            end do
        end do
    end subroutine test_perform_permutation_test_index_bounds

    !> `baseline_mode` (argument 8) is one of the three MODE_BASELINE_* values (1 to 3); 0 and 4,
    !| just outside, are ERR_INVALID_INPUT.
    subroutine test_perform_permutation_test_invalid_mode()
        integer(int32), parameter :: n_samples = 3, n_permutations = 2
        real(real64) :: trajectories(MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS)
        real(real64) :: local_contributions(MARK_N_TIMEPOINTS, n_permutations), total_contributions(n_permutations)
        integer(int32) :: ierr

        call fill_sample_marking(n_samples, trajectories)
        call perform_permutation_test(trajectories, MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS, MARK_FACTOR_IDX, &
                                      MARK_DEPENDENT_IDX, 1, 0_int32, n_permutations, local_contributions, &
                                      total_contributions, 42_int32, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_perform_permutation_test_invalid_mode: 0", 8)
        call perform_permutation_test(trajectories, MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS, MARK_FACTOR_IDX, &
                                      MARK_DEPENDENT_IDX, 1, 4_int32, n_permutations, local_contributions, &
                                      total_contributions, 42_int32, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_perform_permutation_test_invalid_mode: 4", 8)
    end subroutine test_perform_permutation_test_invalid_mode

    !> NaN, +Inf and -Inf in `trajectories` (argument 1) are ERR_NAN_INF, even where no
    !| permutation reads: the value sits in the factor of the last sample, while the factor is
    !| taken from sample 1.
    subroutine test_perform_permutation_test_rejects_nan_and_inf()
        integer(int32), parameter :: n_samples = 3, n_permutations = 2
        real(real64) :: trajectories(MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS)
        real(real64) :: local_contributions(MARK_N_TIMEPOINTS, n_permutations), total_contributions(n_permutations)
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            call fill_sample_marking(n_samples, trajectories)
            trajectories(MARK_FACTOR_IDX, n_samples, MARK_N_TIMEPOINTS) = non_finite_value(i_value)
            call perform_permutation_test(trajectories, MARK_N_FACTORS, n_samples, MARK_N_TIMEPOINTS, &
                                          MARK_FACTOR_IDX, MARK_DEPENDENT_IDX, 1, MODE_BASELINE_RAW, n_permutations, &
                                          local_contributions, total_contributions, 42_int32, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_perform_permutation_test_rejects_nan_and_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_perform_permutation_test_rejects_nan_and_inf

    !> n_samples = 1 leaves no other sample to draw (issue #99 asks for a sample "!= sample"), so
    !| it is ERR_INVALID_INPUT at `n_samples` (argument 3), checked before any draw: n_samples has
    !| the minimum 2 (FES). Before, the draw picked sample 2 and read `trajectories` past its end.
    subroutine test_perform_permutation_test_single_sample()
        real(real64) :: trajectories(MARK_N_FACTORS, 1, MARK_N_TIMEPOINTS)
        real(real64) :: local_contributions(MARK_N_TIMEPOINTS, 2), total_contributions(2)
        integer(int32) :: ierr

        call fill_sample_marking(1, trajectories)
        call perform_permutation_test(trajectories, MARK_N_FACTORS, 1, MARK_N_TIMEPOINTS, MARK_FACTOR_IDX, &
                                      MARK_DEPENDENT_IDX, 1, MODE_BASELINE_RAW, 2, local_contributions, &
                                      total_contributions, 42_int32, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_perform_permutation_test_single_sample", 3)
    end subroutine test_perform_permutation_test_single_sample

end module mod_test_perform_permutation_test
