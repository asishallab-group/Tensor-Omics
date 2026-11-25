!> The `k_means_clustering` cases: two well-separated groups from a known initialisation, the
!| iteration limit, one cluster, as many clusters as points, one point, equal points, ties, many
!| points, the bounds of `n_clusters`, the dimension checks and the NaN/Inf checks.
!|
!| The initialisation is the caller's: `centroids` goes in as the starting centroids and comes back
!| as the final ones, and nothing is random. So every expected label, count and centroid is derived
!| by hand from the documented loop -- assign every point to its closest centroid, recompute every
!| centroid as the mean of its points, stop when no label changes or after `max_iterations`. The
!| points are small integers whose group means are sums divided by a power of two (or exact
!| quotients), so every comparison has tolerance 0.
module mod_test_k_means_clustering
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_clustering, only: k_means_clustering
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_clustering_fixtures, only: UNWRITTEN_INT, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    implicit none
    private
    public :: get_all_tests_k_means_clustering

contains

    !> Every case of `k_means_clustering`.
    function get_all_tests_k_means_clustering() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(11))
        all_tests(1) = test_case("test_k_means_clustering_values", test_k_means_clustering_values)
        all_tests(2) = test_case("test_k_means_clustering_max_iterations", test_k_means_clustering_max_iterations)
        all_tests(3) = test_case("test_k_means_clustering_one_cluster", test_k_means_clustering_one_cluster)
        all_tests(4) = test_case("test_k_means_clustering_as_many_clusters_as_points", &
                                 test_k_means_clustering_as_many_clusters_as_points)
        all_tests(5) = test_case("test_k_means_clustering_one_point", test_k_means_clustering_one_point)
        all_tests(6) = test_case("test_k_means_clustering_equal_points", test_k_means_clustering_equal_points)
        all_tests(7) = test_case("test_k_means_clustering_ties_go_to_the_later_cluster", &
                                 test_k_means_clustering_ties_go_to_the_later_cluster)
        all_tests(8) = test_case("test_k_means_clustering_many_points", test_k_means_clustering_many_points)
        all_tests(9) = test_case("test_k_means_clustering_n_clusters_bounds", test_k_means_clustering_n_clusters_bounds)
        all_tests(10) = test_case("test_k_means_clustering_dimensions", test_k_means_clustering_dimensions)
        all_tests(11) = test_case("test_k_means_clustering_rejects_nan_and_inf", test_k_means_clustering_rejects_nan_and_inf)
    end function get_all_tests_k_means_clustering

    !> Runs k-means from `initial_centroids` and checks `ierr` OK, every label, every count and
    !| every final centroid exactly. `labels` and `label_counts` start at UNWRITTEN_INT and have one
    !| element more than the call is told, so an element the procedure skips and a write past the
    !| end both show. `max_iterations` is passed on only if present, so its default is used otherwise.
    subroutine check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, expected_labels, &
                             expected_counts, expected_centroids, msg, max_iterations)
        integer(int32), intent(in) :: n_clusters, n_points, n_dims
        real(real64), intent(in) :: data_points(n_dims, n_points), initial_centroids(n_dims, n_clusters)
        integer(int32), intent(in) :: expected_labels(n_points), expected_counts(n_clusters)
        real(real64), intent(in) :: expected_centroids(n_dims, n_clusters)
        character(*), intent(in) :: msg
        integer(int32), intent(in), optional :: max_iterations
        real(real64) :: centroids(n_dims, n_clusters)
        integer(int32) :: labels(n_points + 1), label_counts(n_clusters + 1), ierr

        centroids = initial_centroids
        labels = UNWRITTEN_INT
        label_counts = UNWRITTEN_INT
        call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, &
                                max_iterations=max_iterations, ierr=ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_int(labels, expected_labels, n_points, msg//": labels")
        call assert_equal_array_int(label_counts, expected_counts, n_clusters, msg//": label_counts")
        call assert_equal_array_real(centroids, expected_centroids, n_dims*n_clusters, 0.0_real64, msg//": centroids", &
                                     n_rows=n_dims)
        call assert_equal_int(labels(n_points + 1), UNWRITTEN_INT, msg//": no label past the last point")
        call assert_equal_int(label_counts(n_clusters + 1), UNWRITTEN_INT, msg//": no count past the last cluster")
    end subroutine check_k_means

    !> Runs k-means on input it must refuse and checks the error and that `centroids`, which the
    !| procedure may overwrite, comes back exactly as it went in, NaN included.
    subroutine check_refused(n_clusters, data_points, n_points, n_dims, initial_centroids, expected_code, &
                             expected_arg_pos, msg)
        integer(int32), intent(in) :: n_clusters, n_points, n_dims
        real(real64), intent(in) :: data_points(*), initial_centroids(*)
        integer(int32), intent(in) :: expected_code, expected_arg_pos
        character(*), intent(in) :: msg
        ! sized for every refused call in this suite: at most 2 x 3 centroids
        real(real64) :: centroids(6)
        integer(int32) :: labels(4), label_counts(4), ierr, n_centroid_elements

        n_centroid_elements = max(0_int32, min(6_int32, max(n_dims, 0_int32)*max(n_clusters, 0_int32)))
        centroids = 0.0_real64
        centroids(1:n_centroid_elements) = initial_centroids(1:n_centroid_elements)
        call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr=ierr)
        call assert_err(ierr, expected_code, msg, arg_pos=expected_arg_pos)
        call assert_identical_array_real(centroids, initial_centroids, n_centroid_elements, msg//": centroids unchanged")
    end subroutine check_refused

    !> Kept from the old suite (`test_k_means_clustering`, first run): (1, 1), (2, 2), (8, 8),
    !| (10, 10), started from the first two points (the old comment said "first and last").
    !| - iteration 1, centroids (1, 1), (2, 2): point 1 is at 0 from centroid 1; points 2, 3, 4 are
    !|   closer to (2, 2) (0, 72, 128 against 2, 98, 162), so labels [1, 2, 2, 2] and centroid 2
    !|   moves to (20/3, 20/3);
    !| - iteration 2: point 2 is at 2 from (1, 1) and at 2*(14/3)**2 from (20/3, 20/3), so it changes
    !|   to cluster 1; labels [1, 1, 2, 2], centroids (1.5, 1.5) and (9, 9);
    !| - iteration 3: no label changes, so the run stops with those centroids.
    subroutine test_k_means_clustering_values()
        integer(int32), parameter :: n_dims = 2, n_points = 4, n_clusters = 2
        real(real64) :: data_points(n_dims, n_points), initial_centroids(n_dims, n_clusters)
        real(real64) :: expected_centroids(n_dims, n_clusters)

        data_points = reshape([1.0_real64, 1.0_real64, 2.0_real64, 2.0_real64, 8.0_real64, 8.0_real64, &
                               10.0_real64, 10.0_real64], shape(data_points))
        initial_centroids = data_points(:, 1:2)
        expected_centroids = reshape([1.5_real64, 1.5_real64, 9.0_real64, 9.0_real64], shape(expected_centroids))

        call check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, [1, 1, 2, 2], [2, 2], &
                           expected_centroids, "test_k_means_clustering_values")
    end subroutine test_k_means_clustering_values

    !> Kept from the old suite (`test_k_means_clustering`, second run, max_iterations 1, whose
    !| centroid 20/3 was not exact): the last point is (11, 11) instead of (10, 10), so every mean
    !| is exact. From (1, 1) and (2, 2):
    !| - one iteration assigns [1, 2, 2, 2], recomputes (1, 1) and (7, 7) = (2 + 8 + 11)/3, and stops
    !|   with the labels of that assignment;
    !| - two iterations reassign point 2 to (1, 1) (2 against 50) and recompute (1.5, 1.5) and
    !|   (9.5, 9.5); the run stops at the limit, though a label changed;
    !| - the default (300) reaches the same end through a third iteration in which nothing changes.
    subroutine test_k_means_clustering_max_iterations()
        integer(int32), parameter :: n_dims = 2, n_points = 4, n_clusters = 2
        real(real64) :: data_points(n_dims, n_points), initial_centroids(n_dims, n_clusters)
        real(real64) :: after_one(n_dims, n_clusters), converged(n_dims, n_clusters)

        data_points = reshape([1.0_real64, 1.0_real64, 2.0_real64, 2.0_real64, 8.0_real64, 8.0_real64, &
                               11.0_real64, 11.0_real64], shape(data_points))
        initial_centroids = data_points(:, 1:2)
        after_one = reshape([1.0_real64, 1.0_real64, 7.0_real64, 7.0_real64], shape(after_one))
        converged = reshape([1.5_real64, 1.5_real64, 9.5_real64, 9.5_real64], shape(converged))

        call check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, [1, 2, 2, 2], [1, 3], &
                           after_one, "test_k_means_clustering_max_iterations: one iteration", max_iterations=1_int32)
        call check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, [1, 1, 2, 2], [2, 2], &
                           converged, "test_k_means_clustering_max_iterations: two iterations", max_iterations=2_int32)
        call check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, [1, 1, 2, 2], [2, 2], &
                           converged, "test_k_means_clustering_max_iterations: the default")
    end subroutine test_k_means_clustering_max_iterations

    !> One cluster takes every point, and its centroid is their mean wherever it starts:
    !| (1, 2), (3, 6), (5, 10), (7, -2) from (100, -100) give the mean (16/4, 16/4) = (4, 4).
    subroutine test_k_means_clustering_one_cluster()
        integer(int32), parameter :: n_dims = 2, n_points = 4, n_clusters = 1
        real(real64) :: data_points(n_dims, n_points), initial_centroids(n_dims, n_clusters)
        real(real64) :: expected_centroids(n_dims, n_clusters)

        data_points = reshape([1.0_real64, 2.0_real64, 3.0_real64, 6.0_real64, 5.0_real64, 10.0_real64, &
                               7.0_real64, -2.0_real64], shape(data_points))
        initial_centroids(:, 1) = [100.0_real64, -100.0_real64]
        expected_centroids(:, 1) = [4.0_real64, 4.0_real64]

        call check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, [1, 1, 1, 1], [4], &
                           expected_centroids, "test_k_means_clustering_one_cluster")
    end subroutine test_k_means_clustering_one_cluster

    !> The upper bound of `n_clusters`, `n_points`: started from the points in reverse order, every
    !| point is at 0 from its own centroid and farther from every other, so the labels are [4, 3, 2, 1],
    !| every count is 1, and the centroids stay where they started. 1-D points 0, 3, 7, 12.
    subroutine test_k_means_clustering_as_many_clusters_as_points()
        integer(int32), parameter :: n_dims = 1, n_points = 4, n_clusters = 4
        real(real64) :: data_points(n_dims, n_points), initial_centroids(n_dims, n_clusters)

        data_points(1, :) = [0.0_real64, 3.0_real64, 7.0_real64, 12.0_real64]
        initial_centroids(1, :) = [12.0_real64, 7.0_real64, 3.0_real64, 0.0_real64]

        call check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, [4, 3, 2, 1], [1, 1, 1, 1], &
                           initial_centroids, "test_k_means_clustering_as_many_clusters_as_points")
    end subroutine test_k_means_clustering_as_many_clusters_as_points

    !> The smallest sizes, one point in one dimension and one cluster: label 1, count 1, and the
    !| centroid moves from -3 onto the point 5.
    subroutine test_k_means_clustering_one_point()
        real(real64) :: data_points(1, 1), initial_centroids(1, 1), expected_centroids(1, 1)

        data_points = 5.0_real64
        initial_centroids = -3.0_real64
        expected_centroids = 5.0_real64

        call check_k_means(1_int32, data_points, 1_int32, 1_int32, initial_centroids, [1], [1], expected_centroids, &
                           "test_k_means_clustering_one_point")
    end subroutine test_k_means_clustering_one_point

    !> Points that are all equal, (5, 5) three times, from (5, 5) and (7, 7): every point is at 0
    !| from centroid 1, so labels [1, 1, 1], counts [3, 0], and centroid 1 stays (5, 5). Where the
    !| empty cluster's centroid ends up the documentation does not say, so it is not checked here.
    subroutine test_k_means_clustering_equal_points()
        integer(int32), parameter :: n_dims = 2, n_points = 3, n_clusters = 2
        real(real64) :: data_points(n_dims, n_points), centroids(n_dims, n_clusters)
        integer(int32) :: labels(n_points), label_counts(n_clusters), ierr

        data_points = 5.0_real64
        centroids(:, 1) = [5.0_real64, 5.0_real64]
        centroids(:, 2) = [7.0_real64, 7.0_real64]

        call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_k_means_clustering_equal_points: ierr")
        call assert_equal_array_int(labels, [1, 1, 1], n_points, "test_k_means_clustering_equal_points: labels")
        call assert_equal_array_int(label_counts, [3, 0], n_clusters, "test_k_means_clustering_equal_points: label_counts")
        call assert_equal_array_real(centroids(:, 1), [5.0_real64, 5.0_real64], n_dims, 0.0_real64, &
                                     "test_k_means_clustering_equal_points: centroid 1")
    end subroutine test_k_means_clustering_equal_points

    !> Kept from the old suite (`test_k_means_assign_cluster_helper`, point 2, which tested the
    !| internal helper directly): a point exactly as far from two centroids goes to the later one.
    !| The documentation does not say how a tie is broken; this pins the rule the code has (`<=`
    !| against the closest so far).
    !| - 1-D points 0, 5, 10 from 0 and 10, one iteration: 5 is at 25 from both, so labels [1, 2, 2]
    !|   and the centroids become 0 and 7.5;
    !| - three equal points from two equal centroids (the documentation says the centroids "should be
    !|   unique", and that this is not checked): every point ties, so all go to cluster 2, which
    !|   ends at their value 4.
    subroutine test_k_means_clustering_ties_go_to_the_later_cluster()
        real(real64) :: data_points(1, 3), initial_centroids(1, 2), expected_centroids(1, 2)
        real(real64) :: centroids(1, 2)
        integer(int32) :: labels(3), label_counts(2), ierr

        data_points(1, :) = [0.0_real64, 5.0_real64, 10.0_real64]
        initial_centroids(1, :) = [0.0_real64, 10.0_real64]
        expected_centroids(1, :) = [0.0_real64, 7.5_real64]
        call check_k_means(2_int32, data_points, 3_int32, 1_int32, initial_centroids, [1, 2, 2], [1, 2], &
                           expected_centroids, "test_k_means_clustering_ties_go_to_the_later_cluster: midpoint", &
                           max_iterations=1_int32)

        data_points = 4.0_real64
        centroids = 4.0_real64
        call k_means_clustering(2_int32, data_points, 3_int32, 1_int32, centroids, labels, label_counts, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_k_means_clustering_ties_go_to_the_later_cluster: equal centroids, ierr")
        call assert_equal_array_int(labels, [2, 2, 2], 3_int32, &
                                    "test_k_means_clustering_ties_go_to_the_later_cluster: equal centroids, labels")
        call assert_equal_array_int(label_counts, [0, 3], 2_int32, &
                                    "test_k_means_clustering_ties_go_to_the_later_cluster: equal centroids, label_counts")
        call assert_equal_real(centroids(1, 2), 4.0_real64, 0.0_real64, &
                               "test_k_means_clustering_ties_go_to_the_later_cluster: equal centroids, centroid 2")
    end subroutine test_k_means_clustering_ties_go_to_the_later_cluster

    !> 30 points in 3 dimensions, three groups of 10 around (0, 0, 0), (100, 0, 0) and (0, 100, 0),
    !| interleaved (point i is in group mod(i - 1, 3) + 1), so the labels are [1, 2, 3, 1, 2, 3, ...]
    !| only if each point is read as its own column. Each group is its centre plus the offsets
    !| +-e_x, +-e_y, +-e_z, +-(1, 1, 1) and 0 twice, which sum to 0, so the centroids end exactly at
    !| the centres. Started from the centre plus e_x of each group: the groups are 100 apart, so the
    !| first assignment is final.
    subroutine test_k_means_clustering_many_points()
        integer(int32), parameter :: n_dims = 3, n_points = 30, n_clusters = 3, n_per_group = 10
        real(real64) :: data_points(n_dims, n_points), initial_centroids(n_dims, n_clusters)
        real(real64) :: centres(n_dims, n_clusters), offsets(n_dims, n_per_group)
        integer(int32) :: expected_labels(n_points), i_point, i_group

        centres = reshape([0.0_real64, 0.0_real64, 0.0_real64, 100.0_real64, 0.0_real64, 0.0_real64, &
                           0.0_real64, 100.0_real64, 0.0_real64], shape(centres))
        offsets = reshape([1.0_real64, 0.0_real64, 0.0_real64, -1.0_real64, 0.0_real64, 0.0_real64, &
                           0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, -1.0_real64, 0.0_real64, &
                           0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, -1.0_real64, &
                           1.0_real64, 1.0_real64, 1.0_real64, -1.0_real64, -1.0_real64, -1.0_real64, &
                           0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64], shape(offsets))
        do i_point = 1, n_points
            i_group = mod(i_point - 1_int32, n_clusters) + 1_int32
            expected_labels(i_point) = i_group
            data_points(:, i_point) = centres(:, i_group) + offsets(:, (i_point - 1)/n_clusters + 1)
        end do
        initial_centroids = data_points(:, 1:n_clusters)

        call check_k_means(n_clusters, data_points, n_points, n_dims, initial_centroids, expected_labels, &
                           [n_per_group, n_per_group, n_per_group], centres, "test_k_means_clustering_many_points")
    end subroutine test_k_means_clustering_many_points

    !> `n_clusters` lies in [1, n_points] (DM_MIN(1), DM_MAX(n_points)), so 0, -1 and n_points + 1
    !| are refused, naming `n_clusters`; the bounds themselves are accepted
    !| (`test_k_means_clustering_one_cluster`, `test_k_means_clustering_as_many_clusters_as_points`).
    subroutine test_k_means_clustering_n_clusters_bounds()
        real(real64) :: data_points(2, 2), initial_centroids(2, 3)

        data_points = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64], shape(data_points))
        initial_centroids = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64], &
                                    shape(initial_centroids))

        call check_refused(0_int32, data_points, 2_int32, 2_int32, initial_centroids, ERR_INVALID_INPUT, 1_int32, &
                           "test_k_means_clustering_n_clusters_bounds: n_clusters 0")
        call check_refused(-1_int32, data_points, 2_int32, 2_int32, initial_centroids, ERR_INVALID_INPUT, 1_int32, &
                           "test_k_means_clustering_n_clusters_bounds: n_clusters -1")
        call check_refused(3_int32, data_points, 2_int32, 2_int32, initial_centroids, ERR_INVALID_INPUT, 1_int32, &
                           "test_k_means_clustering_n_clusters_bounds: n_clusters n_points + 1")
    end subroutine test_k_means_clustering_n_clusters_bounds

    !> No dimensions is empty input and a negative count invalid, naming `n_dims`. No points is
    !| refused too, but naming `n_clusters`: its bound [1, n_points] is checked first and cannot
    !| hold for n_points < 1, so `n_points`' own check is never what reports it.
    subroutine test_k_means_clustering_dimensions()
        real(real64) :: data_points(2, 2), initial_centroids(2, 1)

        data_points = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64], shape(data_points))
        initial_centroids(:, 1) = [1.0_real64, 2.0_real64]

        call check_refused(1_int32, data_points, 2_int32, 0_int32, initial_centroids, ERR_EMPTY_INPUT, 4_int32, &
                           "test_k_means_clustering_dimensions: n_dims 0")
        call check_refused(1_int32, data_points, 2_int32, -1_int32, initial_centroids, ERR_INVALID_INPUT, 4_int32, &
                           "test_k_means_clustering_dimensions: n_dims -1")
        call check_refused(1_int32, data_points, 0_int32, 2_int32, initial_centroids, ERR_INVALID_INPUT, 1_int32, &
                           "test_k_means_clustering_dimensions: n_points 0")
        call check_refused(1_int32, data_points, -1_int32, 2_int32, initial_centroids, ERR_INVALID_INPUT, 1_int32, &
                           "test_k_means_clustering_dimensions: n_points -1")
    end subroutine test_k_means_clustering_dimensions

    !> NaN, +Inf and -Inf are refused in the points (argument 2) and in the starting centroids
    !| (argument 5), first and last element each, and the centroids come back untouched.
    subroutine test_k_means_clustering_rejects_nan_and_inf()
        real(real64) :: data_points(2, 3), initial_centroids(2, 2)
        integer(int32) :: i_value, i_position, position

        do i_value = 1, N_NON_FINITE
            do i_position = 1, 2
                data_points = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64], &
                                      shape(data_points))
                initial_centroids = reshape([1.0_real64, 2.0_real64, 5.0_real64, 6.0_real64], shape(initial_centroids))
                position = merge(1_int32, 6_int32, i_position == 1)
                data_points(mod(position - 1, 2) + 1, (position - 1)/2 + 1) = non_finite_value(i_value)
                call check_refused(2_int32, data_points, 3_int32, 2_int32, initial_centroids, ERR_NAN_INF, 2_int32, &
                                   "test_k_means_clustering_rejects_nan_and_inf: "//trim(NON_FINITE_NAMES(i_value))// &
                                   " in data_points, element "//position)

                data_points = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64], &
                                      shape(data_points))
                position = merge(1_int32, 4_int32, i_position == 1)
                initial_centroids(mod(position - 1, 2) + 1, (position - 1)/2 + 1) = non_finite_value(i_value)
                call check_refused(2_int32, data_points, 3_int32, 2_int32, initial_centroids, ERR_NAN_INF, 5_int32, &
                                   "test_k_means_clustering_rejects_nan_and_inf: "//trim(NON_FINITE_NAMES(i_value))// &
                                   " in centroids, element "//position)
            end do
        end do
    end subroutine test_k_means_clustering_rejects_nan_and_inf
end module mod_test_k_means_clustering
