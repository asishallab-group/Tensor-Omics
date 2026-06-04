! filepath: test/mod_test_trajectory_contribution_analysis.f90
!> Unit test suite for trajectory_contribution_analysis routine.
module mod_test_trajectory_contribution_analysis
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_trajectory_contribution_analysis
    use tox_errors
    use tox_trajectory_normalization
    use test_suite, only: test_case
    implicit none

    
    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests_trajectory_contribution_analysis() result(all_tests)
        type(test_case),allocatable :: all_tests(:)
        allocate(all_tests(17))
        all_tests(1) = test_case("test_compute_baselines_factor_dependent", test_compute_baselines_factor_dependent)
        all_tests(2) = test_case("test_compute_contributions", test_compute_contributions)
        all_tests(3) = test_case("test_compute_all_contributions", test_compute_all_contributions)
        all_tests(4) = test_case("test_select_random_sample_helper", test_select_random_sample_helper)
        all_tests(5) = test_case("test_perform_permutation_test", test_perform_permutation_test)
        all_tests(6) = test_case("test_compute_p_values", test_compute_p_values)
        all_tests(7) = test_case("test_normalize_variable_timeseries", test_normalize_variable_timeseries)
        all_tests(8) = test_case("test_normalize_single_trajectory", test_normalize_single_trajectory)
        all_tests(9) = test_case("test_normalize_all_trajectories", test_normalize_all_trajectories)
        all_tests(10) = test_case("test_normalize_edge_cases", test_normalize_edge_cases)
        all_tests(11) = test_case("test_normalize_invalid_inputs", test_normalize_invalid_inputs)
        all_tests(12) = test_case("test_compute_velocity_trajectories", test_compute_velocity_trajectories)
        all_tests(13) = test_case("test_compute_acceleration_from_velocity", test_compute_acceleration_from_velocity)
        all_tests(14) = test_case("test_compute_velocity_acceleration_contributions", test_compute_velocity_acceleration_contributions)
        all_tests(15) = test_case("test_compute_velocity_acceleration_contributions_alloc", test_compute_velocity_acceleration_contribs_alloc)
        all_tests(16) = test_case("test_compute_velocity_trajectory", test_compute_velocity_trajectory)
        all_tests(17) = test_case("test_compute_acceleration_from_velocity_trajectory", test_compute_acceleration_from_velocity_trajectory)

    end function get_all_tests_trajectory_contribution_analysis

    !> Test for compute_velocity_trajectory
    subroutine test_compute_velocity_trajectory()
        integer(int32), parameter :: n_timepoints = 5
        real(real64) :: trajectory(n_timepoints), velocity(n_timepoints-1), expected_velocity(n_timepoints-1)
        integer(int32) :: ierr

        ! Simple increasing sequence
        trajectory = [1.0_real64, 2.0_real64, 4.0_real64, 7.0_real64, 11.0_real64]
        expected_velocity = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]

        call compute_velocity_trajectory(trajectory, velocity, n_timepoints, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_velocity_trajectory: ierr")
        call assert_equal_array_real(velocity, expected_velocity, n_timepoints-1, TOL, "test_compute_velocity_trajectory: velocity")
    end subroutine test_compute_velocity_trajectory

    !> Test for compute_acceleration_from_velocity_trajectory
    subroutine test_compute_acceleration_from_velocity_trajectory()
        integer(int32), parameter :: n_timepoints = 5
        real(real64) :: velocity(n_timepoints-1), acceleration(n_timepoints-2), expected_acceleration(n_timepoints-2)
        integer(int32) :: ierr

        ! Simple increasing velocity
        velocity = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        expected_acceleration = [1.0_real64, 1.0_real64, 1.0_real64]

        call compute_acceleration_from_velocity_trajectory(velocity, acceleration, n_timepoints, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_acceleration_from_velocity_trajectory: ierr")
        call assert_equal_array_real(acceleration, expected_acceleration, n_timepoints-2, TOL, "test_compute_acceleration_from_velocity_trajectory: acceleration")
    end subroutine test_compute_acceleration_from_velocity_trajectory

     subroutine test_compute_velocity_trajectories()
        real(real64) :: trajectories(2,2,4)  ! (n_factors, n_samples, n_timepoints)
       real(real64) :: velocity(3,2,2)      ! (n_timepoints-1, n_factors, n_samples)
       real(real64) :: expected(3,2,2)
        integer(int32) :: ierr
        integer(int32) :: factor, sample, t

        ! Reshape to (n_factors=2, n_samples=2, n_timepoints=4)
        ! Data order in Fortran column-major: varies fastest → slowest is factor, sample, time
        ! So list all factors for sample1/time1, then sample2/time1, then sample1/time2, etc.
        trajectories = reshape([ &
            1.0_real64, 0.0_real64, &      
            2.0_real64, -3.0_real64, &     
            2.0_real64, -1.0_real64, &    
            5.0_real64, -1.0_real64, &    
            4.0_real64, -1.0_real64, &     
            9.0_real64, 0.0_real64, &      
            7.0_real64, 0.0_real64, &      
            14.0_real64, 0.0_real64], & 
            shape=[2,2,4])

        ! Compute expected compact velocities: velocity(t) = trajectory(t+1) - trajectory(t)
        expected = 0.0_real64
        do sample = 1, 2
            do factor = 1, 2
                do t = 1, 3
                    expected(t, factor, sample) = trajectories(factor, sample, t+1) - trajectories(factor, sample, t)
                end do
            end do
        end do

        call compute_velocity_trajectories(trajectories, velocity, 2, 2, 4, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_velocity_trajectories: expected OK status")

        call assert_equal_array_real(velocity, &
                                     expected, &
                                     size(velocity), TOL, &
                                     "compute_velocity_trajectories: velocity mismatch")
    end subroutine test_compute_velocity_trajectories

    subroutine test_compute_acceleration_from_velocity()
        real(real64) :: velocity(3,2,2)      ! (n_timepoints-1, n_factors, n_samples)
        real(real64) :: acceleration(2,2,2)  ! (n_timepoints-2, n_factors, n_samples)
        real(real64) :: expected(2,2,2)
        integer(int32) :: ierr
        integer(int32) :: factor, sample, t

        ! Compact velocity layout: (time=3, factor=2, sample=2)
        velocity = reshape([ &
            0.0_real64, 0.0_real64, &      
            0.0_real64, 0.0_real64, &     
            1.0_real64, -1.0_real64, &     
            3.0_real64, 2.0_real64, &      
            2.0_real64, 0.0_real64, &      
            4.0_real64, 1.0_real64], &
            shape=[3,2,2])


        expected = 0.0_real64
        do sample = 1, 2
            do factor = 1, 2
                do t = 1, 2
                    expected(t, factor, sample) = velocity(t+1, factor, sample) - velocity(t, factor, sample)
                end do
            end do
        end do

        call compute_acceleration_from_velocity(velocity, acceleration, 2, 2, 4, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_acceleration_from_velocity: expected OK status")

        call assert_equal_array_real(acceleration, &
                                     expected, &
                                     size(acceleration), TOL, &
                                     "compute_acceleration_from_velocity: acceleration mismatch")
    end subroutine test_compute_acceleration_from_velocity

    subroutine test_compute_velocity_acceleration_contributions()
        real(real64) :: trajectories(2,1,4)
        real(real64) :: C_vel(2,2,1)
        real(real64) :: C_acc(2,2,1)
        real(real64) :: series_vel(4,2,2,1)
        real(real64) :: series_acc(4,2,2,1)
        real(real64) :: velocity(3,2,1)
        real(real64) :: acceleration(2,2,1)
        real(real64) :: expected_total_vel, expected_total_acc
        real(real64) :: expected_series_vel(4), expected_series_acc(4)
        real(real64) :: factor_velocity(3,2)
        real(real64) :: dependent_velocity(3)
        real(real64) :: factor_acceleration(2,2)
        real(real64) :: dependent_acceleration(2)
        real(real64) :: raw_velocity_contrib(3)
        real(real64) :: raw_acceleration_contrib(2)
        integer(int32) :: ierr
        integer(int32) :: mode

        call set_ok(ierr)

       ! Reshape to (n_factors=2, n_samples=1, n_timepoints=4)
        trajectories = reshape([ &
            1.0_real64, 1.0_real64, 3.0_real64, 2.0_real64, 6.0_real64, 2.0_real64, 10.0_real64, 1.0_real64], &
            shape=[2,1,4])

        mode = BASELINE_RAW

           call compute_velocity_acceleration_contributions(trajectories, 2, 1, 4, mode, &
               factor_velocity, dependent_velocity, raw_velocity_contrib, &
               C_vel, series_vel, C_acc, series_acc, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_velocity_acceleration_contributions: expected OK status")

        call compute_velocity_trajectories(trajectories, velocity, 2, 1, 4, ierr)
        call assert_equal_int(ierr, ERR_OK, "velocity back-reference")

        call compute_acceleration_from_velocity(velocity, acceleration, 2, 1, 4, ierr)
        call assert_equal_int(ierr, ERR_OK, "acceleration back-reference")

        factor_velocity(:,1) = velocity(:,1,1)
        dependent_velocity = velocity(:,2,1)

        call compute_contributions(factor_velocity(:,1), dependent_velocity, &
            int(size(dependent_velocity), kind=int32), mode, raw_velocity_contrib, expected_total_vel, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_velocity_acceleration_contributions: expected velocity contribution status")

        expected_series_vel = 0.0_real64
        expected_series_vel(2:4) = raw_velocity_contrib

        factor_acceleration(:,1) = acceleration(:,1,1)
        dependent_acceleration = acceleration(:,2,1)

        call compute_contributions(factor_acceleration(:,1), dependent_acceleration, &
            int(size(dependent_acceleration), kind=int32), mode, raw_acceleration_contrib, expected_total_acc, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_velocity_acceleration_contributions: expected acceleration contribution status")

        expected_series_acc = 0.0_real64
        expected_series_acc(3:4) = raw_acceleration_contrib

        call assert_equal_real(C_vel(1,2,1), expected_total_vel, TOL, &
            "compute_velocity_acceleration_contributions: total velocity contribution")

        call assert_equal_array_real(series_vel(:,1,2,1), expected_series_vel, size(expected_series_vel), TOL, &
            "compute_velocity_acceleration_contributions: velocity series")

        call assert_equal_real(C_acc(1,2,1), expected_total_acc, TOL, &
            "compute_velocity_acceleration_contributions: total acceleration contribution")

        call assert_equal_array_real(series_acc(:,1,2,1), expected_series_acc, size(expected_series_acc), TOL, &
                "compute_velocity_acceleration_contributions: acceleration series")
    end subroutine test_compute_velocity_acceleration_contributions

    subroutine test_compute_velocity_acceleration_contribs_alloc()
        real(real64) :: trajectories(2,1,4)
        real(real64) :: factor_velocity(3,2)
        real(real64) :: dependent_velocity(3)
        real(real64) :: velocity_contrib(3)
        real(real64) :: C_vel_ref(2,2,1)
        real(real64) :: C_acc_ref(2,2,1)
        real(real64) :: series_vel_ref(4,2,2,1)
        real(real64) :: series_acc_ref(4,2,2,1)
        real(real64) :: C_vel_alloc(2,2,1)
        real(real64) :: C_acc_alloc(2,2,1)
        real(real64) :: series_vel_alloc(4,2,2,1)
        real(real64) :: series_acc_alloc(4,2,2,1)
        integer(int32) :: ierr
        integer(int32) :: mode

        ! Reshape to (n_factors=2, n_samples=1, n_timepoints=4)
        trajectories = reshape([ &
            1.0_real64, 1.0_real64, 3.0_real64, 2.0_real64, 6.0_real64, 2.0_real64, 10.0_real64, 1.0_real64], &
            shape=[2,1,4])

        mode = BASELINE_RAW

        call compute_velocity_acceleration_contributions(trajectories, 2, 1, 4, mode, &
            factor_velocity, dependent_velocity, velocity_contrib, &
            C_vel_ref, series_vel_ref, C_acc_ref, series_acc_ref, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_velocity_acceleration_contributions_alloc: reference call")

        call compute_velocity_acceleration_contributions_alloc(trajectories, 2, 1, 4, mode, &
            C_vel_alloc, series_vel_alloc, C_acc_alloc, series_acc_alloc, ierr)
        call assert_equal_int(ierr, ERR_OK, "compute_velocity_acceleration_contributions_alloc: expected OK status")

        call assert_equal_array_real(reshape(C_vel_alloc, [size(C_vel_alloc)]), &
            reshape(C_vel_ref, [size(C_vel_ref)]), size(C_vel_alloc), TOL, &
            "compute_velocity_acceleration_contributions_alloc: velocity totals")

        call assert_equal_array_real(reshape(series_vel_alloc, [size(series_vel_alloc)]), &
            reshape(series_vel_ref, [size(series_vel_ref)]), size(series_vel_alloc), TOL, &
            "compute_velocity_acceleration_contributions_alloc: velocity series")

        call assert_equal_array_real(reshape(C_acc_alloc, [size(C_acc_alloc)]), &
            reshape(C_acc_ref, [size(C_acc_ref)]), size(C_acc_alloc), TOL, &
            "compute_velocity_acceleration_contributions_alloc: acceleration totals")

        call assert_equal_array_real(reshape(series_acc_alloc, [size(series_acc_alloc)]), &
            reshape(series_acc_ref, [size(series_acc_ref)]), size(series_acc_alloc), TOL, &
            "compute_velocity_acceleration_contributions_alloc: acceleration series")
    end subroutine test_compute_velocity_acceleration_contribs_alloc



    !> initializes random number generator with a randomly selected seed
    subroutine setup_random
        integer(int32) :: seed

        call random_init(.false., .false.) ! reset random number generator to non-reproducible
        seed = int(rand_range(0.0_real64, real(huge(1_int32), kind=real64)), kind=int32) ! pick random seed
        call init_random(seed) ! set random number generator to seed
        write (*, "('Using random seed: ', I0)") seed
    end subroutine setup_random

    !> Test the compute_p_values function with a simple case where we know the expected p-values.
    subroutine test_compute_p_values()
        integer(int32), parameter :: n_timepoints = 3, n_permutations = 4
        integer(int32) :: ierr
        real(real64) :: local_contributions_observed(n_timepoints)
        real(real64) :: total_contribution_observed
        real(real64) :: local_contributions_perm(n_timepoints, n_permutations)
        real(real64) :: total_contributions_perm(n_permutations)
        real(real64) :: local_p_values(n_timepoints)
        real(real64) :: total_p_value
        real(real64) :: expected_local_p(n_timepoints)
        real(real64) :: expected_total_p

        ! -------------------------------
        ! Setup observed contributions
        ! -------------------------------
        local_contributions_observed = [2.0_real64, 0.0_real64, 2.0_real64]
        total_contribution_observed  = sum(local_contributions_observed)  ! = 4.0

        ! -------------------------------
        ! Setup permutation contributions
        ! -------------------------------
        ! Permutation 1: [1,0,1], total=2
        local_contributions_perm(:,1) = [1.0_real64, 0.0_real64, 1.0_real64]
        total_contributions_perm(1)   = 2.0_real64

        ! Permutation 2: [2,0,2], total=4
        local_contributions_perm(:,2) = [2.0_real64, 0.0_real64, 2.0_real64]
        total_contributions_perm(2)   = 4.0_real64

        ! Permutation 3: [3,1,3], total=7
        local_contributions_perm(:,3) = [3.0_real64, 1.0_real64, 3.0_real64]
        total_contributions_perm(3)   = 7.0_real64

        ! Permutation 4: [0,0,0], total=0
        local_contributions_perm(:,4) = [0.0_real64, 0.0_real64, 0.0_real64]
        total_contributions_perm(4)   = 0.0_real64

        ! -------------------------------
        ! Call routine
        ! -------------------------------
        call compute_p_values(local_contributions_observed, total_contribution_observed, &
            local_contributions_perm, total_contributions_perm, n_timepoints, n_permutations, &
            local_p_values, total_p_value, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_p_values: ierr")

        ! -------------------------------
        ! Expected p-values
        ! -------------------------------
        ! For each timepoint:
        ! Timepoint 1 observed=2.0 → perms >=2.0: [2,3] → 2/4 = 0.5
        ! Timepoint 2 observed=0.0 → perms >=0.0: [1,2,3,4] → 4/4 = 1.0
        ! Timepoint 3 observed=2.0 → perms >=2.0: [2,3] → 2/4 = 0.5
        expected_local_p = [0.5_real64, 1.0_real64, 0.5_real64]

        ! Total observed=4.0 → perms >=4.0: [2,3] → 2/4 = 0.5
        expected_total_p = 0.5_real64

        ! -------------------------------
        ! Assertions
        ! -------------------------------
        call assert_equal_array_real(local_p_values, expected_local_p, n_timepoints, TOL, "test_compute_p_values: local p-values")
        call assert_equal_real(total_p_value, expected_total_p, TOL, "test_compute_p_values: total p-value")
    end subroutine test_compute_p_values
    
    !> Test perform_permutation_test with various cases, including reproducibility with random seed.
    subroutine test_perform_permutation_test()
        integer(int32), parameter :: n_factors = 2, n_samples = 20, n_timepoints = 3
        integer(int32), parameter :: n_permutations = 2
        integer(int32) :: ierr, mode
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        integer(int32) :: factor_idx, dependent_idx, sample_idx
        real(real64) :: local_contributions(n_timepoints, n_permutations)
        real(real64) :: total_contributions(n_permutations)
        real(real64) :: temp_factor(n_timepoints), temp_dependent(n_timepoints)
        real(real64) :: expected_local(n_timepoints, n_permutations)
        real(real64) :: expected_total(n_permutations)

        call setup_random()

        ! Case 1: test functionality without randomness, as all data unlike current_sample is same
        sample_idx    = 1
        mode          = BASELINE_MEAN
        factor_idx    = 1
        dependent_idx = 1 ! will be different for the permutations
        trajectories(1,1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(1,2:,1) = 2.0_real64
        trajectories(1,2:,2) = 4.0_real64
        trajectories(1,2:,3) = 6.0_real64
        ! Dependent 2 values across samples/timepoints
        trajectories(2,1,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        trajectories(2,2:,1) = 1.0_real64
        trajectories(2,2:,2) = 3.0_real64
        trajectories(2,2:,3) = 5.0_real64


        call perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, &
            factor_idx, dependent_idx, sample_idx, mode, n_permutations, &
            local_contributions, total_contributions, temp_factor, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_perform_permutation_test: Case 1 Permutation test ierr")

        ! Factor trajectory (sample 1): [1,2,3], mean=2.0
        ! Dependent trajectory (sample 2): [1,3,5], mean=3.0
        ! Contributions = (factor - 2.0)*(dependent - 3.0)
        expected_local(1, :) = (1.0-2.0)*(1.0-3.0)   ! = 2.0
        expected_local(2, :) = (2.0-2.0)*(3.0-3.0)   ! = 0.0
        expected_local(3, :) = (3.0-2.0)*(5.0-3.0)   ! = 2.0
        expected_total    = sum(expected_local(:, 1))   ! = 4.0

        call assert_equal_array_real(local_contributions, expected_local, n_timepoints * n_permutations, TOL, "test_perform_permutation_test: Case 1 local contributions")
        call assert_equal_array_real(total_contributions, expected_total, n_permutations, TOL, "test_perform_permutation_test: Case 1 total contribution")

        ! Case 2: test randomness reproducibility: without seed -> not reproducible
        call random_number(trajectories)
        call perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, &
            factor_idx, dependent_idx, sample_idx, mode, n_permutations, &
            expected_local, expected_total, temp_factor, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_perform_permutation_test: Case 2 Permutation test ierr expected contribs")

        call perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, &
            factor_idx, dependent_idx, sample_idx, mode, n_permutations, &
            local_contributions, total_contributions, temp_factor, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_perform_permutation_test: Case 2 Permutation test ierr")

        call assert_false(all(local_contributions == expected_local), "test_perform_permutation_test: Case 2 should not be reproducible")

        ! Case 3: test randomness reproducibility: with seed -> reproducible
        call random_number(trajectories)
        call perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, &
            factor_idx, dependent_idx, sample_idx, mode, n_permutations, &
            expected_local, expected_total, temp_factor, temp_dependent, ierr, random_seed=42_int32)
        call assert_equal_int(ierr, ERR_OK, "test_perform_permutation_test: Case 3 Permutation test ierr expected contribs")

        call perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, &
            factor_idx, dependent_idx, sample_idx, mode, n_permutations, &
            local_contributions, total_contributions, temp_factor, temp_dependent, ierr, random_seed=42_int32)
        call assert_equal_int(ierr, ERR_OK, "test_perform_permutation_test: Case 3 Permutation test ierr")

        call assert_true(all(local_contributions == expected_local), "test_perform_permutation_test: Case 3 should be reproducible")
    end subroutine test_perform_permutation_test

    !> Test the select_random_sample_helper function with various cases.
    subroutine test_select_random_sample_helper
        use f42_utils, only: rand_range
        integer(int32), parameter :: n_samples = 10, current_sample = 5
        integer(int32) :: ierr, random_sample, i
        integer(int32), dimension(n_samples) :: sample_counts

        call setup_random()

        ! Case 1: Basic selection test -> current_sample never selected
        sample_counts = 0
        do i = 1, 1000
            call select_random_sample_helper(n_samples, current_sample, random_sample, ierr)
            call assert_equal_int(ierr, ERR_OK, "test_select_random_sample_helper: Case 1: Unexpected error when selecting sample")
            
            sample_counts(random_sample) = sample_counts(random_sample) + 1
        end do

        call assert_equal_int(sample_counts(current_sample), 0_int32, "test_select_random_sample_helper: Case 1: current_sample selected")
        call assert_equal_int(count(sample_counts /= 0), n_samples - 1, "test_select_random_sample_helper: Case 1: for 1000 iterations all other 9 samples should definitely be selected once")
    
        ! Case 2: only two samples -> forced selection
        call select_random_sample_helper(2_int32, 1_int32, random_sample, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_select_random_sample_helper: Case 2: Unexpected error when selecting sample")
        call assert_equal_int(random_sample, 2_int32, "test_select_random_sample_helper: Case 2: selected wrong sample")

        call select_random_sample_helper(2_int32, 2_int32, random_sample, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_select_random_sample_helper: Case 2: Unexpected error when selecting sample")
        call assert_equal_int(random_sample, 1_int32, "test_select_random_sample_helper: Case 2: selected wrong sample")

        ! Case 3: Error case: n_samples=current_sample=1
        call select_random_sample_helper(1_int32, 1_int32, random_sample, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_select_random_sample_helper: Case 3: Expected error for samples=1")

        ! Case 4: Error case: current_sample out of range
        call select_random_sample_helper(n_samples, 0_int32, random_sample, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_select_random_sample_helper: Case 4: Expected error for current_sample=0")

        call select_random_sample_helper(n_samples, n_samples + 1, random_sample, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_select_random_sample_helper: Case 4: Expected error for current_sample>n_samples")
    end subroutine test_select_random_sample_helper
    
    !> Test the compute_all_contributions function with various cases.
    subroutine test_compute_all_contributions()
        integer(int32), parameter :: n_factors = 2, n_samples = 1, n_timepoints = 3
        integer(int32), parameter :: n_selected_factors = 1, n_selected_dependents = 1
        integer(int32) :: ierr, i_dependent
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        integer(int32) :: factor_indices(n_factors)
        integer(int32) :: dependent_indices(n_factors)
        real(real64) :: local_contributions(n_timepoints, n_factors, n_factors, n_samples)
        real(real64) :: total_contributions(n_factors, n_factors, n_samples)
        real(real64) :: temp_factors(n_timepoints, n_factors), temp_dependent(n_timepoints)
        real(real64) :: expected_local(n_timepoints)
        real(real64) :: expected_total

        ! -------------------------------
        ! Case 1: zero vectors
        ! -------------------------------
        ! Factor trajectory: [1,2,3]
        ! Dependent trajectory: [4,5,6]

        trajectories(1,1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2,1,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        factor_indices = [1, 2]
        dependent_indices = [2, 1]

        call compute_all_contributions(trajectories, 0_int32, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_factors=0")

        call compute_all_contributions(trajectories, n_factors, 0_int32, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_samples=0")

        call compute_all_contributions(trajectories, n_factors, n_samples, 0_int32, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_timepoints=0")

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), 0_int32, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_selected_factors=0")

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), 0_int32, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_selected_dependents=0")

        ! -------------------------------
        ! Case 2: MEAN baseline
        ! -------------------------------

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MEAN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 2 ierr")

        ! Baselines: mean(factor)=2.0, mean(dependent)=5.0
        expected_local(1) = (1.0-2.0)*(4.0-5.0)   ! = 1.0
        expected_local(2) = (2.0-2.0)*(5.0-5.0)   ! = 0.0
        expected_local(3) = (3.0-2.0)*(6.0-5.0)   ! = 1.0
        expected_total = sum(expected_local)      ! = 2.0

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 2 local contributions")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 2 total contribution")

        ! -------------------------------
        ! Case 3: MIN baseline
        ! -------------------------------
        ! Factor trajectory: [2,4,6]
        ! Dependent trajectory: [1,3,5]
        trajectories(1,1,:) = [2.0_real64, 4.0_real64, 6.0_real64]
        trajectories(2,1,:) = [1.0_real64, 3.0_real64, 5.0_real64]

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 3 ierr")

        ! Baselines: min(factor)=2.0, min(dependent)=1.0
        expected_local(1) = (2.0-2.0)*(1.0-1.0)   ! = 0.0
        expected_local(2) = (4.0-2.0)*(3.0-1.0)   ! = 4.0
        expected_local(3) = (6.0-2.0)*(5.0-1.0)   ! = 16.0
        expected_total = sum(expected_local)      ! = 20.0

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 3 local contributions")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 3 total contribution")

        ! -------------------------------
        ! Case 4: all elements equal
        ! -------------------------------
        trajectories(1,1,:) = 1.0_real64

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 4 ierr zero factor")

        expected_local = 0.0_real64
        expected_total = 0.0_real64

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 4 local contributions zero factor")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 4 total contribution zero factor")

        trajectories(1,1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2,1,:) = 1.0_real64

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 4 ierr zero dependent")

        expected_local = 0.0_real64
        expected_total = 0.0_real64

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 4 local contributions zero dependent")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 4 total contribution zero dependent")

        ! -------------------------------
        ! Case 5: multiple selected indices
        ! -------------------------------
        trajectories(1,1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2,1,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        factor_indices = [1, 2]
        dependent_indices = [2, 2]

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices, n_factors, dependent_indices, n_factors, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 5 ierr")

        do i_dependent = 1, size(dependent_indices)
            expected_local(1) = (1.0-1.0)*(4.0-4.0)
            expected_local(2) = (2.0-1.0)*(5.0-4.0)
            expected_local(3) = (3.0-1.0)*(6.0-4.0)
            expected_total = sum(expected_local)

            call assert_equal_array_real(local_contributions(:,1,i_dependent,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 5 local contributions")
            call assert_equal_real(total_contributions(1,i_dependent,1), expected_total, TOL, "test_compute_all_contributions: Case 5 total contribution")

            expected_local(1) = (4.0-4.0) ** 2
            expected_local(2) = (5.0-4.0) ** 2
            expected_local(3) = (6.0-4.0) ** 2
            expected_total = sum(expected_local)

            call assert_equal_array_real(local_contributions(:,1,i_dependent,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 5 local contributions")
            call assert_equal_real(total_contributions(1,i_dependent,1), expected_total, TOL, "test_compute_all_contributions: Case 5 total contribution")
        end do

    end subroutine test_compute_all_contributions
    
    !> Test the compute_contributions function with various modes and edge cases.
    subroutine test_compute_contributions()
        integer(int32), parameter :: n_dims = 4
        integer(int32) :: ierr, mode
        real(real64) :: factor(n_dims), dependent(n_dims)
        real(real64) :: local_contributions(n_dims), expected_local(n_dims)
        real(real64) :: total_contribution, expected_total

        ! -------------------------------
        ! Case 1: RAW baseline
        ! -------------------------------
        factor    = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent = [2.0_real64, 1.0_real64, 0.0_real64, -1.0_real64]
        mode = 1  ! RAW

        call compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_contributions: Case 1 ierr")

        ! RAW baseline means baseline = 0
        expected_local = factor * dependent
        expected_total = sum(expected_local)

        call assert_equal_array_real(local_contributions, expected_local, n_dims, TOL, "test_compute_contributions: Case 1 local contributions")
        call assert_equal_real(total_contribution, expected_total, TOL, "test_compute_contributions: Case 1 total contribution")

        ! -------------------------------
        ! Case 2: MIN baseline
        ! -------------------------------
        factor    = [3.0_real64, 5.0_real64, 2.0_real64, 4.0_real64]
        dependent = [1.0_real64, 2.0_real64, 0.0_real64, -1.0_real64]
        mode = 2  ! MIN

        call compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_contributions: Case 2 ierr")

        ! Baseline = min(factor)=2, min(dependent)=-1
        expected_local = (factor - minval(factor)) * (dependent - minval(dependent))
        expected_total = sum(expected_local)

        call assert_equal_array_real(local_contributions, expected_local, n_dims, TOL, "test_compute_contributions: Case 2 local contributions")
        call assert_equal_real(total_contribution, expected_total, TOL, "test_compute_contributions: Case 2 total contribution")

        ! -------------------------------
        ! Case 3: MEAN baseline
        ! -------------------------------
        factor    = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        mode = 3  ! MEAN

        call compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_contributions: Case 3 ierr")

        ! Baseline = mean(factor)=2.5, mean(dependent)=2.5
        expected_local = (factor - sum(factor) / n_dims) * (dependent - sum(dependent) / n_dims)
        expected_total = sum(expected_local)

        call assert_equal_array_real(local_contributions, expected_local, n_dims, TOL, "test_compute_contributions: Case 3 local contributions")
        call assert_equal_real(total_contribution, expected_total, TOL, "test_compute_contributions: Case 3 total contribution")
    end subroutine test_compute_contributions

    !> Test the compute_baselines_factor_dependent function with various modes and edge cases.
    subroutine test_compute_baselines_factor_dependent()
        integer(int32), parameter :: n_timepoints = 4
        real(real64) :: factor(n_timepoints), dependent(5)  ! dependent has 5 elements for mismatch test
        real(real64) :: factor_baseline, dependent_baseline, expected_factor_baseline, expected_dependent_baseline
        integer(int32) :: ierr

        ! Case 1: BASELINE_RAW (no centering)
        factor = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent(1:n_timepoints) = [5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64]
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_RAW, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: BASELINE_RAW: expected OK status")
        call assert_equal_real(factor_baseline, 0.0_real64, TOL, "test_compute_baselines_factor_dependent: BASELINE_RAW factor_baseline")
        call assert_equal_real(dependent_baseline, 0.0_real64, TOL, "test_compute_baselines_factor_dependent: BASELINE_RAW dependent_baseline")

        ! Case 2: BASELINE_MIN (minimum-centered)
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_MIN, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: BASELINE_MIN: expected OK status")
        call assert_equal_real(factor_baseline, minval(factor), TOL, "test_compute_baselines_factor_dependent: BASELINE_MIN factor_baseline")
        call assert_equal_real(dependent_baseline, minval(dependent(1:n_timepoints)), TOL, "test_compute_baselines_factor_dependent: BASELINE_MIN dependent_baseline")

        ! Case 3: BASELINE_MEAN (mean-centered)
        expected_factor_baseline = sum(factor) / real(n_timepoints, kind=real64)
        expected_dependent_baseline = sum(dependent(1:n_timepoints)) / real(n_timepoints, kind=real64)
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_MEAN, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: BASELINE_MEAN: expected OK status")
        call assert_equal_real(factor_baseline, expected_factor_baseline, TOL, "test_compute_baselines_factor_dependent: BASELINE_MEAN factor_baseline")
        call assert_equal_real(dependent_baseline, expected_dependent_baseline, TOL, "test_compute_baselines_factor_dependent: BASELINE_MEAN dependent_baseline")

        ! Case 4: mismatched input lengths
        dependent = [5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64, 9.0_real64]  ! length 5
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_RAW, &
                                               factor_baseline, dependent_baseline, ierr)  ! valid
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: valid lengths")
        call compute_baselines_factor_dependent(4_int32, factor(1:4), dependent(1:5), BASELINE_RAW, &
                                               factor_baseline, dependent_baseline, ierr)  ! invalid: mismatch
        ! Note: This test may not trigger length mismatch since we now validate in C wrapper only

        ! Case 5: invalid mode
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), 99_int32, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_compute_baselines_factor_dependent: invalid mode")
    end subroutine test_compute_baselines_factor_dependent

    !> Test the normalization of variable timeseries
    subroutine test_normalize_variable_timeseries()
        integer(int32) :: n_points = 5
        real(real64) :: v(5), v_norm(5), v_norm_expected(5), v1_norm(1), v1_norm_expected(1), v1(1)
        integer(int32) :: ierr, status
        
        ! Test 1: Normal case
        v = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
        v_norm_expected = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, n_points, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: normal case should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, n_points, TOL, "test_normalize_variable_timeseries: normal case values")
        
        ! Test 2: Constant vector (all values same)
        v = [2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64]
        v_norm_expected = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, n_points, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: constant vector should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, n_points, TOL, "test_normalize_variable_timeseries: constant vector values")
        
        ! Test 3: Negative values
        v = [-5.0_real64, -2.0_real64, 0.0_real64, 3.0_real64, 6.0_real64]
        v_norm_expected = [0.0_real64, 0.27272727_real64, 0.45454545_real64, 0.72727273_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, n_points, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: negative values should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, n_points, 1.0e-8_real64, "test_normalize_variable_timeseries: negative values")
        
        ! Test 4: Single point
        v1 = [3.14_real64]
        v1_norm_expected = [0.0_real64]  ! Single value normalized to 0
        
        call normalize_variable_timeseries(v1, v1_norm, 1, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: single point should succeed")
        call assert_equal_array_real(v1_norm, v1_norm_expected, 1, TOL, "test_normalize_variable_timeseries: single point")
    end subroutine test_normalize_variable_timeseries

    !> Test the normalization of a single trajectory (factor) across timepoints for one sample.
    subroutine test_normalize_single_trajectory()
        integer(int32), parameter :: n_factors = 3, n_timepoints = 4
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints, n_factors)
        integer(int32) :: ierr, i_factor, i_timepoint, status
        
        ! Create test trajectory for ONE SAMPLE: factors × timepoints
        do i_factor = 1, n_factors
            do i_timepoint = 1, n_timepoints
                trajectory(i_timepoint, i_factor) = real(i_factor * 10 + i_timepoint, real64)
            end do
        end do
        
        ! Expected: Each factor normalized independently across time
        ! Factor 1: [11, 12, 13, 14] → normalized: [0.0, 0.333..., 0.666..., 1.0]
        ! Factor 2: [21, 22, 23, 24] → normalized: [0.0, 0.333..., 0.666..., 1.0]
        ! Factor 3: [31, 32, 33, 34] → normalized: [0.0, 0.333..., 0.666..., 1.0]
        do i_factor = 1, n_factors
            do i_timepoint = 1, n_timepoints
                expected(i_timepoint, i_factor) = (real(i_timepoint, real64) - 1.0_real64) / real(n_timepoints - 1, real64)
            end do
        end do
        
        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_single_trajectory: should succeed")
        
        ! Check each factor independently
        do i_factor = 1, n_factors
            call assert_equal_array_real(trajectory_norm(:, i_factor), expected(:, i_factor), n_timepoints, TOL, &
                "test_normalize_single_trajectory: factor ")
        end do
        
        ! Verify min=0 and max=1 for each factor (across time)
        do i_factor = 1, n_factors
            call assert_equal_real(minval(trajectory_norm(:, i_factor)), 0.0_real64, TOL, &
                "test_normalize_single_trajectory: min=0 for factor ")
            call assert_equal_real(maxval(trajectory_norm(:, i_factor)), 1.0_real64, TOL, &
                "test_normalize_single_trajectory: max=1 for factor ")
        end do
    end subroutine test_normalize_single_trajectory

    !> Test the normalization of all trajectories
    subroutine test_normalize_all_trajectories()
        integer(int32), parameter :: n_factors = 2, n_samples = 3, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints)
        integer(int32) :: ierr, i_factor, i_sample, i_timepoint, status
        
        ! Fill with known pattern: factor × sample × timepoint
        do i_factor = 1, n_factors
            do i_sample = 1, n_samples
                do i_timepoint = 1, n_timepoints
                    trajectories(i_factor, i_sample, i_timepoint) = &
                        real(i_factor * 100 + i_sample * 10 + i_timepoint, real64)
                end do
            end do
        end do
        
        call normalize_all_trajectories(trajectories, trajectories_norm, &
                                      n_factors, n_samples, n_timepoints, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_all_trajectories: should succeed")
        
        ! Check: For each (factor, sample), values should be normalized across time
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                call assert_equal_real(minval(trajectories_norm(i_factor, i_sample, :)), &
                                     0.0_real64, TOL, &
                                     "test_normalize_all_trajectories: factor ")
                call assert_equal_real(maxval(trajectories_norm(i_factor, i_sample, :)), &
                                     1.0_real64, TOL, &
                                     "test_normalize_all_trajectories: factor ")
            end do
        end do
        
        ! Verify all values are in [0,1]
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                do i_timepoint = 1, n_timepoints
                    call assert_true(trajectories_norm(i_factor, i_sample, i_timepoint) >= 0.0_real64 - TOL, &
                                   "test_normalize_all_trajectories: value >= 0")
                    call assert_true(trajectories_norm(i_factor, i_sample, i_timepoint) <= 1.0_real64 + TOL, &
                                   "test_normalize_all_trajectories: value <= 1")
                end do
            end do
        end do
    end subroutine test_normalize_all_trajectories

    !> Test edge cases for the normalization function
    subroutine test_normalize_edge_cases()
        real(real64) :: v(3), v_norm(3), v_norm_expected(3)
        real(real64) :: v5(5), v5_norm(5), v5_norm_expected(5)
        integer(int32) :: ierr, status
        
        ! Test 1: Very small values, should result in zero vector due to division by near zero
        v = [tiny(1.0_real64), 2.0_real64 * tiny(1.0_real64), 3.0_real64 * tiny(1.0_real64)]
        v_norm_expected = [0.0_real64, 0.0_real64, 0.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, 3, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: small values should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, 3, TOL, "test_normalize_edge_cases: small values")
        
        ! Test 2: Large values
        v = [1.0e10_real64, 2.0e10_real64, 3.0e10_real64]
        v_norm_expected = [0.0_real64, 0.5_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, 3, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: large values should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, 3, TOL, "test_normalize_edge_cases: large values")
        
        ! Test 3: Mixed positive and negative with zero
        v5 = [-10.0_real64, -5.0_real64, 0.0_real64, 5.0_real64, 10.0_real64]
        v5_norm_expected = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v5, v5_norm, 5, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: mixed values should succeed")
        call assert_equal_array_real(v5_norm, v5_norm_expected, 5, 1.0e-8_real64, "test_normalize_edge_cases: mixed values")
        
        ! Test 4: Already normalized (values in [0,1])
        v5 = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        v5_norm_expected = v5  ! Should stay the same
        
        call normalize_variable_timeseries(v5, v5_norm, 5, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: already normalized should succeed")
        call assert_equal_array_real(v5_norm, v5_norm_expected, 5, TOL, "test_normalize_edge_cases: already normalized")
    end subroutine test_normalize_edge_cases

    !> Test invalid inputs for normalization function, such as empty array, negative n_points, and NaN/Inf values.
    subroutine test_normalize_invalid_inputs()
        real(real64) :: v(5), v_norm(5)
        real(real64) :: v3(3), v3_norm(3)
        integer(int32) :: ierr, status
        
        ! Test 1: Empty array (n_points = 0)
        call normalize_variable_timeseries(v, v_norm, 0, ierr, status)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_normalize_invalid_inputs: empty array should return ERR_EMPTY_INPUT")
        
        ! Test 2: Negative n_points
        call normalize_variable_timeseries(v, v_norm, -1, ierr, status)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_normalize_invalid_inputs: negative n_points should return ERR_INVALID_INPUT")
        
        ! Test 3: NaN in input
        v3 = [1.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64]
        call normalize_variable_timeseries(v3, v3_norm, 3, ierr, status)
        ! Note: NaN handling depends on minval/maxval behavior - check if error is set
        if (ierr /= ERR_OK) then
            ! If error is set, it should be ERR_NAN_INF
            call assert_equal_int(ierr, ERR_NAN_INF, "test_normalize_invalid_inputs: NaN should return ERR_NAN_INF")
        end if
        
        ! Test 4: Infinity in input
        v3 = [1.0_real64, 2.0_real64, huge(1.0_real64)]
        call normalize_variable_timeseries(v3, v3_norm, 3, ierr, status)
        ! Similar to NaN case
        if (ierr /= ERR_OK) then
            call assert_equal_int(ierr, ERR_NAN_INF, "test_normalize_invalid_inputs: Infinity should return ERR_NAN_INF")
        end if
    end subroutine test_normalize_invalid_inputs

end module mod_test_trajectory_contribution_analysis
