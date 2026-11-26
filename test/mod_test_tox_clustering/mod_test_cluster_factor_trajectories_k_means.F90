!> The `cluster_factor_trajectories_k_means` cases, its expert entry included: two groups of
!| samples whose trajectories evolve alike, the flattening order the documentation states, one
!| iteration, the smallest sizes, the bounds of `n_clusters`, the dimension checks and the NaN/Inf
!| checks.
!|
!| One sample is one point: its trajectory flattened as `factor1_t1, factor2_t1, factor1_t2, ...`.
!| The k-means loop itself is `mod_test_k_means_clustering`'s subject; these cases make sure the
!| right numbers reach it. So the samples have 2 factors and 3 timepoints (two different extents,
!| so a swap shows), and every value differs from its neighbours, so a centroid read in any other
!| order than the documented one differs from the expected one. The four samples, flattened:
!| - s1 = [1, 10, 3, 30, 5, 50] and s3 = [3, 12, 5, 32, 7, 52], s3 = s1 + 2 everywhere;
!| - s2 = [20, 2, 10, 0, 0, 1] and s4 = [21, 3, 11, 3, 4, 3].
!| s1 and s3 are at 24 from each other and s2 and s4 at 32, and every pair across is above 3000.
module mod_test_cluster_factor_trajectories_k_means
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_clustering, only: cluster_factor_trajectories_k_means, cluster_factor_trajectories_k_means_expert
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_clustering_fixtures, only: UNWRITTEN_INT, UNWRITTEN_REAL, N_NON_FINITE, NON_FINITE_NAMES, &
                                                non_finite_value
    implicit none
    private
    public :: get_all_tests_cluster_factor_trajectories_k_means

    !> The fixture's extents: factors, samples, timepoints, and the length of a flattened trajectory.
    integer(int32), parameter :: FACTORS = 2, SAMPLES = 4, TIMEPOINTS = 3, FLAT_LENGTH = FACTORS*TIMEPOINTS
    !> The four samples flattened, one per column, as the module documentation above lists them.
    real(real64), parameter :: FLAT(FLAT_LENGTH, SAMPLES) = reshape([ &
                               1.0_real64, 10.0_real64, 3.0_real64, 30.0_real64, 5.0_real64, 50.0_real64, &
                               20.0_real64, 2.0_real64, 10.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, &
                               3.0_real64, 12.0_real64, 5.0_real64, 32.0_real64, 7.0_real64, 52.0_real64, &
                               21.0_real64, 3.0_real64, 11.0_real64, 3.0_real64, 4.0_real64, 3.0_real64], &
                               [FLAT_LENGTH, SAMPLES])

contains

    !> Every case of `cluster_factor_trajectories_k_means` and its expert entry.
    function get_all_tests_cluster_factor_trajectories_k_means() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(9))
        all_tests(1) = test_case("test_cluster_factor_trajectories_k_means_values", &
                                 test_cluster_factor_trajectories_k_means_values)
        all_tests(2) = test_case("test_cluster_factor_trajectories_k_means_expert_values", &
                                 test_cluster_factor_trajectories_k_means_expert_values)
        all_tests(3) = test_case("test_cluster_factor_trajectories_k_means_one_iteration", &
                                 test_cluster_factor_trajectories_k_means_one_iteration)
        all_tests(4) = test_case("test_cluster_factor_trajectories_k_means_smallest_sizes", &
                                 test_cluster_factor_trajectories_k_means_smallest_sizes)
        all_tests(5) = test_case("test_cluster_factor_trajectories_k_means_n_clusters_bounds", &
                                 test_cluster_factor_trajectories_k_means_n_clusters_bounds)
        all_tests(6) = test_case("test_cluster_factor_trajectories_k_means_dimensions", &
                                 test_cluster_factor_trajectories_k_means_dimensions)
        all_tests(7) = test_case("test_cluster_factor_trajectories_k_means_rejects_nan_and_inf", &
                                 test_cluster_factor_trajectories_k_means_rejects_nan_and_inf)
        all_tests(8) = test_case("test_cluster_factor_trajectories_k_means_expert_refuses_alike", &
                                 test_cluster_factor_trajectories_k_means_expert_refuses_alike)
        all_tests(9) = test_case("test_cluster_factor_trajectories_k_means_iterations_below_one", &
                                 test_cluster_factor_trajectories_k_means_iterations_below_one)
    end function get_all_tests_cluster_factor_trajectories_k_means

    !> The four samples as `trajectories(factor, sample, timepoint)`, written out value by value
    !| from the list in the module documentation.
    subroutine fill_trajectories(trajectories)
        real(real64), intent(out) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS)

        ! s1: t1 (1, 10), t2 (3, 30), t3 (5, 50)
        trajectories(:, 1, 1) = [1.0_real64, 10.0_real64]
        trajectories(:, 1, 2) = [3.0_real64, 30.0_real64]
        trajectories(:, 1, 3) = [5.0_real64, 50.0_real64]
        ! s2: t1 (20, 2), t2 (10, 0), t3 (0, 1)
        trajectories(:, 2, 1) = [20.0_real64, 2.0_real64]
        trajectories(:, 2, 2) = [10.0_real64, 0.0_real64]
        trajectories(:, 2, 3) = [0.0_real64, 1.0_real64]
        ! s3: t1 (3, 12), t2 (5, 32), t3 (7, 52)
        trajectories(:, 3, 1) = [3.0_real64, 12.0_real64]
        trajectories(:, 3, 2) = [5.0_real64, 32.0_real64]
        trajectories(:, 3, 3) = [7.0_real64, 52.0_real64]
        ! s4: t1 (21, 3), t2 (11, 3), t3 (4, 3)
        trajectories(:, 4, 1) = [21.0_real64, 3.0_real64]
        trajectories(:, 4, 2) = [11.0_real64, 3.0_real64]
        trajectories(:, 4, 3) = [4.0_real64, 3.0_real64]
    end subroutine fill_trajectories

    !> Checks `ierr` OK, the labels, the counts and the final centroids exactly; `labels` and
    !| `label_counts` have one element more than the call is told, which must stay UNWRITTEN_INT.
    subroutine check_result(ierr, labels, label_counts, centroids, expected_labels, expected_counts, expected_centroids, msg)
        integer(int32), intent(in) :: ierr, labels(SAMPLES + 1), label_counts(3)
        real(real64), intent(in) :: centroids(FLAT_LENGTH, 2), expected_centroids(FLAT_LENGTH, 2)
        integer(int32), intent(in) :: expected_labels(SAMPLES), expected_counts(2)
        character(*), intent(in) :: msg

        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_int(labels, expected_labels, SAMPLES, msg//": labels")
        call assert_equal_array_int(label_counts, expected_counts, 2_int32, msg//": label_counts")
        call assert_equal_array_real(centroids, expected_centroids, FLAT_LENGTH*2_int32, 0.0_real64, msg//": centroids", &
                                     n_rows=FLAT_LENGTH)
        call assert_equal_int(labels(SAMPLES + 1), UNWRITTEN_INT, msg//": no label past the last sample")
        call assert_equal_int(label_counts(3), UNWRITTEN_INT, msg//": no count past the last cluster")
    end subroutine check_result

    !> Kept and given values from the old R suite's `cluster_factor_trajectories_k_means` smoke test
    !| (random starting centroids, only shapes checked). From s1 and s2: s3 is at 24 from s1 and
    !| s4 at 32 from s2, so the labels are [1, 2, 1, 2], and the centroids the means
    !| (s1 + s3)/2 = [2, 11, 4, 31, 6, 51] and (s2 + s4)/2 = [20.5, 2.5, 10.5, 1.5, 2, 2], flattened
    !| as documented; the second iteration changes nothing.
    subroutine test_cluster_factor_trajectories_k_means_values()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH, 2), expected(FLAT_LENGTH, 2)
        integer(int32) :: labels(SAMPLES + 1), label_counts(3), ierr

        call fill_trajectories(trajectories)
        centroids(:, 1) = FLAT(:, 1)
        centroids(:, 2) = FLAT(:, 2)
        expected(:, 1) = [2.0_real64, 11.0_real64, 4.0_real64, 31.0_real64, 6.0_real64, 51.0_real64]
        expected(:, 2) = [20.5_real64, 2.5_real64, 10.5_real64, 1.5_real64, 2.0_real64, 2.0_real64]
        labels = UNWRITTEN_INT
        label_counts = UNWRITTEN_INT

        call cluster_factor_trajectories_k_means(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, &
                                                 labels, label_counts, 300_int32, ierr)
        call check_result(ierr, labels, label_counts, centroids, [1, 2, 1, 2], [2, 2], expected, &
                          "test_cluster_factor_trajectories_k_means_values")
    end subroutine test_cluster_factor_trajectories_k_means_values

    !> The expert entry on the same input gives the same result, and leaves every sample's
    !| flattened trajectory in its column of the work array, as the work array's documentation
    !| says ("each sample's trajectory as one column").
    subroutine test_cluster_factor_trajectories_k_means_expert_values()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH, 2), expected(FLAT_LENGTH, 2)
        real(real64) :: work(FLAT_LENGTH, SAMPLES), factors_work(FLAT_LENGTH, 2)
        integer(int32) :: labels(SAMPLES + 1), label_counts(3), ierr

        call fill_trajectories(trajectories)
        centroids(:, 1) = FLAT(:, 1)
        centroids(:, 2) = FLAT(:, 2)
        expected(:, 1) = [2.0_real64, 11.0_real64, 4.0_real64, 31.0_real64, 6.0_real64, 51.0_real64]
        expected(:, 2) = [20.5_real64, 2.5_real64, 10.5_real64, 1.5_real64, 2.0_real64, 2.0_real64]
        labels = UNWRITTEN_INT
        label_counts = UNWRITTEN_INT
        work = UNWRITTEN_REAL

        call cluster_factor_trajectories_k_means_expert(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, &
                                                        centroids, labels, label_counts, work, factors_work, 300_int32, ierr)
        call check_result(ierr, labels, label_counts, centroids, [1, 2, 1, 2], [2, 2], expected, &
                          "test_cluster_factor_trajectories_k_means_expert_values")
        call assert_equal_array_real(work, FLAT, FLAT_LENGTH*SAMPLES, 0.0_real64, &
                                     "test_cluster_factor_trajectories_k_means_expert_values: flattened trajectories", &
                                     n_rows=FLAT_LENGTH)
    end subroutine test_cluster_factor_trajectories_k_means_expert_values

    !> `max_iterations` reaches the loop: from s1 and s3, which belong together, one iteration
    !| assigns s2 and s4 to s1 as well (s2 at 3800 from s1 and 4088 from s3, s4 at 3452 and
    !| 3692), so the labels are [1, 1, 2, 1], and the centroids (s1 + s2 + s4)/3 =
    !| [14, 5, 8, 11, 3, 18] and s3. With more iterations the run would go on to [1, 2, 1, 2].
    subroutine test_cluster_factor_trajectories_k_means_one_iteration()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH, 2), expected(FLAT_LENGTH, 2)
        integer(int32) :: labels(SAMPLES + 1), label_counts(3), ierr

        call fill_trajectories(trajectories)
        centroids(:, 1) = FLAT(:, 1)
        centroids(:, 2) = FLAT(:, 3)
        expected(:, 1) = [14.0_real64, 5.0_real64, 8.0_real64, 11.0_real64, 3.0_real64, 18.0_real64]
        expected(:, 2) = FLAT(:, 3)
        labels = UNWRITTEN_INT
        label_counts = UNWRITTEN_INT

        call cluster_factor_trajectories_k_means(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, &
                                                 labels, label_counts, 1_int32, ierr)
        call check_result(ierr, labels, label_counts, centroids, [1, 1, 2, 1], [3, 1], expected, &
                          "test_cluster_factor_trajectories_k_means_one_iteration")
    end subroutine test_cluster_factor_trajectories_k_means_one_iteration

    !> One factor, one sample, one timepoint, one cluster: the trajectory is the single value 7,
    !| the label 1, the count 1, and the centroid moves onto it from -7.
    subroutine test_cluster_factor_trajectories_k_means_smallest_sizes()
        real(real64) :: trajectories(1, 1, 1), centroids(1, 1)
        integer(int32) :: labels(1), label_counts(1), ierr

        trajectories = 7.0_real64
        centroids = -7.0_real64
        call cluster_factor_trajectories_k_means(1_int32, trajectories, 1_int32, 1_int32, 1_int32, centroids, labels, &
                                                 label_counts, 300_int32, ierr)
        call assert_err(ierr, ERR_OK, "test_cluster_factor_trajectories_k_means_smallest_sizes: ierr")
        call assert_equal_int(labels(1), 1_int32, "test_cluster_factor_trajectories_k_means_smallest_sizes: label")
        call assert_equal_int(label_counts(1), 1_int32, "test_cluster_factor_trajectories_k_means_smallest_sizes: count")
        call assert_equal_real(centroids(1, 1), 7.0_real64, 0.0_real64, &
                               "test_cluster_factor_trajectories_k_means_smallest_sizes: centroid")
    end subroutine test_cluster_factor_trajectories_k_means_smallest_sizes

    !> Calls the plain entry on input it must refuse; checks the error and that the centroids,
    !| which it may overwrite, come back as they went in, NaN included.
    subroutine check_refused(n_clusters, trajectories, n_factors, n_samples, n_timepoints, initial_centroids, &
                             n_centroid_elements, expected_code, expected_arg_pos, msg)
        integer(int32), intent(in) :: n_clusters, n_factors, n_samples, n_timepoints, n_centroid_elements
        real(real64), intent(in) :: trajectories(*), initial_centroids(*)
        integer(int32), intent(in) :: expected_code, expected_arg_pos
        character(*), intent(in) :: msg
        real(real64) :: centroids(FLAT_LENGTH*(SAMPLES + 1))
        integer(int32) :: labels(SAMPLES + 1), label_counts(SAMPLES + 1), ierr

        centroids = 0.0_real64
        centroids(1:n_centroid_elements) = initial_centroids(1:n_centroid_elements)
        call cluster_factor_trajectories_k_means(n_clusters, trajectories, n_factors, n_samples, n_timepoints, centroids, &
                                                 labels, label_counts, 300_int32, ierr)
        call assert_err(ierr, expected_code, msg, arg_pos=expected_arg_pos)
        call assert_identical_array_real(centroids, initial_centroids, n_centroid_elements, msg//": centroids unchanged")
    end subroutine check_refused

    !> `n_clusters` lies in [1, n_samples], so 0, -1 and n_samples + 1 are refused, naming it.
    subroutine test_cluster_factor_trajectories_k_means_n_clusters_bounds()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH, SAMPLES + 1)

        call fill_trajectories(trajectories)
        centroids = 1.0_real64
        call check_refused(0_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, 0_int32, &
                           ERR_INVALID_INPUT, 1_int32, "test_cluster_factor_trajectories_k_means_n_clusters_bounds: 0")
        call check_refused(-1_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, 0_int32, &
                           ERR_INVALID_INPUT, 1_int32, "test_cluster_factor_trajectories_k_means_n_clusters_bounds: -1")
        call check_refused(SAMPLES + 1_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, &
                           FLAT_LENGTH*(SAMPLES + 1_int32), ERR_INVALID_INPUT, 1_int32, &
                           "test_cluster_factor_trajectories_k_means_n_clusters_bounds: n_samples + 1")
    end subroutine test_cluster_factor_trajectories_k_means_n_clusters_bounds

    !> No factors or no timepoints is empty input, a negative count invalid, each naming itself.
    !| No samples is refused naming `n_clusters`, whose bound [1, n_samples] is checked first and
    !| cannot hold for n_samples < 1.
    subroutine test_cluster_factor_trajectories_k_means_dimensions()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH)
        character(*), parameter :: msg = "test_cluster_factor_trajectories_k_means_dimensions: "

        call fill_trajectories(trajectories)
        centroids = 1.0_real64
        call check_refused(1_int32, trajectories, 0_int32, SAMPLES, TIMEPOINTS, centroids, 0_int32, &
                           ERR_EMPTY_INPUT, 3_int32, msg//"n_factors 0")
        call check_refused(1_int32, trajectories, -1_int32, SAMPLES, TIMEPOINTS, centroids, 0_int32, &
                           ERR_INVALID_INPUT, 3_int32, msg//"n_factors -1")
        call check_refused(1_int32, trajectories, FACTORS, SAMPLES, 0_int32, centroids, 0_int32, &
                           ERR_EMPTY_INPUT, 5_int32, msg//"n_timepoints 0")
        call check_refused(1_int32, trajectories, FACTORS, SAMPLES, -1_int32, centroids, 0_int32, &
                           ERR_INVALID_INPUT, 5_int32, msg//"n_timepoints -1")
        call check_refused(1_int32, trajectories, FACTORS, 0_int32, TIMEPOINTS, centroids, FLAT_LENGTH, &
                           ERR_INVALID_INPUT, 1_int32, msg//"n_samples 0")
        call check_refused(1_int32, trajectories, FACTORS, -1_int32, TIMEPOINTS, centroids, FLAT_LENGTH, &
                           ERR_INVALID_INPUT, 1_int32, msg//"n_samples -1")
    end subroutine test_cluster_factor_trajectories_k_means_dimensions

    !> NaN, +Inf and -Inf are refused in the trajectories (argument 2) and in the starting
    !| centroids (argument 6), at the first and the last element each.
    subroutine test_cluster_factor_trajectories_k_means_rejects_nan_and_inf()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH, 2)
        integer(int32) :: i_value
        character(*), parameter :: msg = "test_cluster_factor_trajectories_k_means_rejects_nan_and_inf: "

        do i_value = 1, N_NON_FINITE
            centroids(:, 1) = FLAT(:, 1)
            centroids(:, 2) = FLAT(:, 2)
            call fill_trajectories(trajectories)
            trajectories(1, 1, 1) = non_finite_value(i_value)
            call check_refused(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, FLAT_LENGTH*2_int32, &
                               ERR_NAN_INF, 2_int32, msg//trim(NON_FINITE_NAMES(i_value))//" first in trajectories")
            call fill_trajectories(trajectories)
            trajectories(FACTORS, SAMPLES, TIMEPOINTS) = non_finite_value(i_value)
            call check_refused(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, FLAT_LENGTH*2_int32, &
                               ERR_NAN_INF, 2_int32, msg//trim(NON_FINITE_NAMES(i_value))//" last in trajectories")

            call fill_trajectories(trajectories)
            centroids(1, 1) = non_finite_value(i_value)
            call check_refused(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, FLAT_LENGTH*2_int32, &
                               ERR_NAN_INF, 6_int32, msg//trim(NON_FINITE_NAMES(i_value))//" first in centroids")
            centroids(1, 1) = FLAT(1, 1)
            centroids(FLAT_LENGTH, 2) = non_finite_value(i_value)
            call check_refused(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, centroids, FLAT_LENGTH*2_int32, &
                               ERR_NAN_INF, 6_int32, msg//trim(NON_FINITE_NAMES(i_value))//" last in centroids")
        end do
    end subroutine test_cluster_factor_trajectories_k_means_rejects_nan_and_inf

    !> The expert entry validates as the plain one does: `n_clusters` out of range (argument 1),
    !| no factors (3), no timepoints (5), a NaN trajectory (2), an infinite centroid (6).
    subroutine test_cluster_factor_trajectories_k_means_expert_refuses_alike()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH, 2), work(FLAT_LENGTH, SAMPLES)
        real(real64) :: factors_work(FLAT_LENGTH, 2)
        integer(int32) :: labels(SAMPLES), label_counts(2), ierr
        character(*), parameter :: msg = "test_cluster_factor_trajectories_k_means_expert_refuses_alike: "

        call fill_trajectories(trajectories)
        centroids = 1.0_real64
        call cluster_factor_trajectories_k_means_expert(0_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, &
                                                        centroids, labels, label_counts, work, factors_work, 300_int32, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, msg//"n_clusters 0", arg_pos=1_int32)
        call cluster_factor_trajectories_k_means_expert(2_int32, trajectories, 0_int32, SAMPLES, TIMEPOINTS, &
                                                        centroids, labels, label_counts, work, factors_work, 300_int32, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, msg//"n_factors 0", arg_pos=3_int32)
        call cluster_factor_trajectories_k_means_expert(2_int32, trajectories, FACTORS, SAMPLES, 0_int32, &
                                                        centroids, labels, label_counts, work, factors_work, 300_int32, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, msg//"n_timepoints 0", arg_pos=5_int32)
        trajectories(2, 3, 2) = non_finite_value(1_int32)
        call cluster_factor_trajectories_k_means_expert(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, &
                                                        centroids, labels, label_counts, work, factors_work, 300_int32, ierr)
        call assert_err(ierr, ERR_NAN_INF, msg//"NaN trajectory", arg_pos=2_int32)
        call fill_trajectories(trajectories)
        centroids(3, 2) = non_finite_value(2_int32)
        call cluster_factor_trajectories_k_means_expert(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, &
                                                        centroids, labels, label_counts, work, factors_work, 300_int32, ierr)
        call assert_err(ierr, ERR_NAN_INF, msg//"+Inf centroid", arg_pos=6_int32)
    end subroutine test_cluster_factor_trajectories_k_means_expert_refuses_alike

    !> `max_iterations` is required here and lies in [1, ...) (DM_MIN(1)): 0 and -1 are refused,
    !| naming it -- argument 9 of the plain entry, 11 of the expert one, after its two work arrays --
    !| and the centroids come back unchanged. Below 1 the k-means loop never ran, and the run used to
    !| report success with every label 0 and `label_counts` never written.
    subroutine test_cluster_factor_trajectories_k_means_iterations_below_one()
        real(real64) :: trajectories(FACTORS, SAMPLES, TIMEPOINTS), centroids(FLAT_LENGTH, 2), work(FLAT_LENGTH, SAMPLES)
        real(real64) :: factors_work(FLAT_LENGTH, 2)
        real(real64) :: initial_centroids(FLAT_LENGTH, 2)
        integer(int32) :: labels(SAMPLES), label_counts(2), ierr, max_iterations, i_entry
        character(*), parameter :: ENTRY_NAMES(2) = ["plain ", "expert"]
        integer(int32), parameter :: ARG_POS(2) = [9_int32, 11_int32]
        character(*), parameter :: msg = "test_cluster_factor_trajectories_k_means_iterations_below_one: "

        call fill_trajectories(trajectories)
        initial_centroids(:, 1) = FLAT(:, 1)
        initial_centroids(:, 2) = FLAT(:, 2)
        do i_entry = 1, 2
            do max_iterations = 0, -1, -1
                centroids = initial_centroids
                if (i_entry == 1) then
                    call cluster_factor_trajectories_k_means(2_int32, trajectories, FACTORS, SAMPLES, TIMEPOINTS, &
                                                             centroids, labels, label_counts, max_iterations, ierr)
                else
                    call cluster_factor_trajectories_k_means_expert(2_int32, trajectories, FACTORS, SAMPLES, &
                                                                    TIMEPOINTS, centroids, labels, label_counts, work, &
                                                                    factors_work, max_iterations, ierr)
                end if
                call assert_err(ierr, ERR_INVALID_INPUT, msg//trim(ENTRY_NAMES(i_entry))//" entry, max_iterations "// &
                                max_iterations, arg_pos=ARG_POS(i_entry))
                call assert_equal_array_real(centroids, initial_centroids, FLAT_LENGTH*2_int32, 0.0_real64, &
                                             msg//trim(ENTRY_NAMES(i_entry))//" entry, max_iterations "// &
                                             max_iterations//": centroids unchanged", n_rows=FLAT_LENGTH)
            end do
        end do
    end subroutine test_cluster_factor_trajectories_k_means_iterations_below_one
end module mod_test_cluster_factor_trajectories_k_means
