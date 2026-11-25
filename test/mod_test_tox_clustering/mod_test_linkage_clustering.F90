!> The `linkage_clustering` cases: for each of the three methods (average/UPGMA, weighted/WPGMA,
!| Ward), the merges, heights and cluster sizes of matrices whose merge order is derived by hand;
!| equal distances, two points, one point, duplicate points, all points equal, many points; the
!| range of `method`, the dimension checks, the distance-matrix checks and the NaN/Inf checks.
!| The matrix must come back unchanged every time, as the documentation promises.
!|
!| A merge is checked as the pair {smaller label, larger label}: positive labels are points,
!| `-k` the cluster made by merge `k`. Which of the two goes into `merge_i` and which into
!| `merge_j` is not documented, so it is not pinned. The documentation does not say how a tie
!| between two equally close pairs is broken either; where a case has ties, the merges pin the
!| rule the code has -- the pair with the lowest column, then the lowest row, of the lower
!| triangle, where a cluster sits at the lowest index of its points.
!|
!| UPGMA and WPGMA average distances that are integers or halves, so their heights are exact and
!| compared with tolerance 0. Ward's update takes square roots, so its heights are compared within
!| WARD_RTOL; the exact values are closed forms derived by hand.
module mod_test_linkage_clustering
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_clustering, only: linkage_clustering
    ! the method codes are the implementation's; the wrapper does not re-export them
    use tox_clustering_impl, only: METHOD_AVERAGE, METHOD_WEIGHTED, METHOD_WARD
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_clustering_fixtures, only: UNWRITTEN_INT, UNWRITTEN_REAL, N_NON_FINITE, NON_FINITE_NAMES, &
                                                non_finite_value
    implicit none
    private
    public :: get_all_tests_linkage_clustering

    !> Ward's heights are square roots of Lance-Williams updates whose order of operations is the
    !| implementation's, not the contract's, so each rounds differently from the closed form; a few
    !| ulps, relative.
    real(real64), parameter :: WARD_RTOL = 4*epsilon(1.0_real64)

    integer(int32), parameter :: N_METHODS = 3
    integer(int32), parameter :: METHODS(N_METHODS) = [METHOD_AVERAGE, METHOD_WEIGHTED, METHOD_WARD]
    character(*), parameter :: METHOD_NAMES(N_METHODS) = ["UPGMA", "WPGMA", "Ward "]

    !> The example of https://en.wikipedia.org/wiki/UPGMA (and /WPGMA), points a..e as 1..5.
    integer(int32), parameter :: WIKI_N = 5
    real(real64), parameter :: WIKI_DISTANCES(WIKI_N, WIKI_N) = reshape([ &
                               0.0_real64, 17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
                               17.0_real64, 0.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
                               21.0_real64, 30.0_real64, 0.0_real64, 28.0_real64, 39.0_real64, &
                               31.0_real64, 34.0_real64, 28.0_real64, 0.0_real64, 43.0_real64, &
                               23.0_real64, 21.0_real64, 39.0_real64, 43.0_real64, 0.0_real64], [WIKI_N, WIKI_N])

    !> Six points on a line, at 0, 10, 40, 11, 60, 62: the distances are the differences. The
    !| first merge, 2 with 4, has point 3 between them and points 5 and 6 after them, so it
    !| updates rows before, between and after the merged pair.
    integer(int32), parameter :: LINE_N = 6
    real(real64), parameter :: LINE_POSITIONS(LINE_N) = [0.0_real64, 10.0_real64, 40.0_real64, 11.0_real64, &
                                                         60.0_real64, 62.0_real64]

contains

    !> Every case of `linkage_clustering`.
    function get_all_tests_linkage_clustering() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(14))
        all_tests(1) = test_case("test_linkage_clustering_wikipedia_example", test_linkage_clustering_wikipedia_example)
        all_tests(2) = test_case("test_linkage_clustering_points_on_a_line", test_linkage_clustering_points_on_a_line)
        all_tests(3) = test_case("test_linkage_clustering_equal_distances", test_linkage_clustering_equal_distances)
        all_tests(4) = test_case("test_linkage_clustering_two_points", test_linkage_clustering_two_points)
        all_tests(5) = test_case("test_linkage_clustering_one_point", test_linkage_clustering_one_point)
        all_tests(6) = test_case("test_linkage_clustering_duplicate_points", test_linkage_clustering_duplicate_points)
        all_tests(7) = test_case("test_linkage_clustering_all_points_equal", test_linkage_clustering_all_points_equal)
        all_tests(8) = test_case("test_linkage_clustering_many_points", test_linkage_clustering_many_points)
        all_tests(9) = test_case("test_linkage_clustering_method_range", test_linkage_clustering_method_range)
        all_tests(10) = test_case("test_linkage_clustering_dimensions", test_linkage_clustering_dimensions)
        all_tests(11) = test_case("test_linkage_clustering_distance_matrix_checks", &
                                  test_linkage_clustering_distance_matrix_checks)
        all_tests(12) = test_case("test_linkage_clustering_rejects_nan_and_inf", test_linkage_clustering_rejects_nan_and_inf)
        all_tests(13) = test_case("test_linkage_clustering_check_order", test_linkage_clustering_check_order)
        all_tests(14) = test_case("test_linkage_clustering_ward_squares_below_tiny", &
                                  test_linkage_clustering_ward_squares_below_tiny)
    end function get_all_tests_linkage_clustering

    !> Clusters a copy of `distances` and checks `ierr` OK, that the copy comes back identical,
    !| every merge as the pair {lower, higher} (unless `check_merges` is false), every height
    !| (exactly when `rtol` is 0, else within `rtol`, relative) and every cluster size (with the
    !| merges). The outputs have one element more than the call is told, which must stay unwritten.
    subroutine check_linkage(n_points, distances, method, expected_lower, expected_higher, expected_heights, rtol, &
                             expected_sizes, check_merges, msg)
        integer(int32), intent(in) :: n_points, method
        real(real64), intent(in) :: distances(n_points, n_points)
        integer(int32), intent(in) :: expected_lower(n_points - 1), expected_higher(n_points - 1)
        integer(int32), intent(in) :: expected_sizes(n_points - 1)
        real(real64), intent(in) :: expected_heights(n_points - 1), rtol
        logical, intent(in) :: check_merges
        character(*), intent(in) :: msg
        real(real64) :: passed(n_points, n_points), heights(n_points)
        integer(int32) :: merge_i(n_points), merge_j(n_points), cluster_sizes(n_points), ierr, n_merges
        integer(int32) :: lower_labels(n_points), higher_labels(n_points)

        n_merges = n_points - 1
        passed = distances
        merge_i = UNWRITTEN_INT
        merge_j = UNWRITTEN_INT
        cluster_sizes = UNWRITTEN_INT
        heights = UNWRITTEN_REAL
        call linkage_clustering(passed, n_points, merge_i, merge_j, heights, cluster_sizes, method, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_identical_array_real(passed, distances, n_points*n_points, msg//": the matrix comes back unchanged")
        if (check_merges) then
            lower_labels = min(merge_i, merge_j)
            higher_labels = max(merge_i, merge_j)
            call assert_equal_array_int(lower_labels, expected_lower, n_merges, msg//": lower label of each merge")
            call assert_equal_array_int(higher_labels, expected_higher, n_merges, msg//": higher label of each merge")
            call assert_equal_array_int(cluster_sizes, expected_sizes, n_merges, msg//": cluster_sizes")
        end if
        if (rtol > 0.0_real64) then
            call assert_allclose_array_real(heights, expected_heights, n_merges, rtol, 0.0_real64, msg//": heights")
        else
            call assert_equal_array_real(heights, expected_heights, n_merges, 0.0_real64, msg//": heights")
        end if
        call assert_equal_int(merge_i(n_points), UNWRITTEN_INT, msg//": nothing written past the last merge_i")
        call assert_equal_int(merge_j(n_points), UNWRITTEN_INT, msg//": nothing written past the last merge_j")
        call assert_equal_int(cluster_sizes(n_points), UNWRITTEN_INT, msg//": nothing written past the last size")
        call assert_equal_real(heights(n_points), UNWRITTEN_REAL, 0.0_real64, msg//": nothing written past the last height")
    end subroutine check_linkage

    !> Calls `linkage_clustering` on input it must refuse; checks the error and that the matrix
    !| comes back identical, NaN included.
    subroutine check_refused(n_points, distances, n_elements, method, expected_code, expected_arg_pos, msg)
        integer(int32), intent(in) :: n_points, n_elements, method, expected_code, expected_arg_pos
        real(real64), intent(in) :: distances(n_elements)
        character(*), intent(in) :: msg
        real(real64) :: passed(n_elements), heights(4)
        integer(int32) :: merge_i(4), merge_j(4), cluster_sizes(4), ierr

        passed = distances
        call linkage_clustering(passed, n_points, merge_i, merge_j, heights, cluster_sizes, method, ierr)
        call assert_err(ierr, expected_code, msg, arg_pos=expected_arg_pos)
        call assert_identical_array_real(passed, distances, n_elements, msg//": the matrix comes back unchanged")
    end subroutine check_refused

    !> Fills `distances` with the matrix of `n_points` points at the same distance `distance` from
    !| each other.
    pure subroutine fill_equidistant(n_points, distance, distances)
        integer(int32), intent(in) :: n_points
        real(real64), intent(in) :: distance
        real(real64), intent(out) :: distances(n_points, n_points)
        integer(int32) :: i_point

        distances = distance
        do i_point = 1, n_points
            distances(i_point, i_point) = 0.0_real64
        end do
    end subroutine fill_equidistant

    !> Kept from the old suite (`test_linkage_methods` case 1, and the internal helpers'
    !| `test_linkage_helpers_xPGMA` and `test_linkage_helpers_ward`, whose intermediate matrices
    !| decide the heights below). Every method merges {1, 2} at 17, then {-1, 5}, then {3, 4} at
    !| 28, then {-3, -2}; sizes 2, 3, 2, 5.
    !| - UPGMA: d(12, c) = (21 + 30)/2 = 25.5, d(12, d) = 32.5, d(12, e) = (23 + 21)/2 = 22, the
    !|   second height; d(12e, c) = (2*25.5 + 39)/3 = 30, d(12e, d) = (2*32.5 + 43)/3 = 36; after
    !|   {3, 4} at 28, the last is (30 + 36)/2 = 33.
    !| - WPGMA: the same 22 and 28; d(12e, c) = (25.5 + 39)/2 = 32.25, d(12e, d) = (32.5 + 43)/2 =
    !|   37.75, the last (32.25 + 37.75)/2 = 35.
    !| - Ward, d(AB, C)**2 = ((nA + nC) d(A, C)**2 + (nB + nC) d(B, C)**2 - nC d(A, B)**2)/(nA + nB + nC):
    !|   d(12, c)**2 = (2*441 + 2*900 - 289)/3 = 2393/3, d(12, d)**2 = (2*961 + 2*1156 - 289)/3 = 1315,
    !|   d(12, e)**2 = (2*529 + 2*441 - 289)/3 = 1651/3, the second height (23.46);
    !|   d(12e, c)**2 = (3*2393/3 + 2*1521 - 1651/3)/4 = 7327/6 (34.9), d(12e, d)**2 = (3*1315 + 2*1849
    !|   - 1651/3)/4 = 10639/6 (42.1), so {3, 4} at 28 comes third; the last is
    !|   ((3 + 1)*7327/6 + (3 + 1)*10639/6 - 3*784)/5 = 28876/15 (43.88).
    subroutine test_linkage_clustering_wikipedia_example()
        integer(int32), parameter :: lower(WIKI_N - 1) = [1, -1, 3, -3], higher(WIKI_N - 1) = [2, 5, 4, -2]
        integer(int32), parameter :: sizes(WIKI_N - 1) = [2, 3, 2, 5]

        call check_linkage(WIKI_N, WIKI_DISTANCES, METHOD_AVERAGE, lower, higher, &
                           [17.0_real64, 22.0_real64, 28.0_real64, 33.0_real64], 0.0_real64, sizes, .true., &
                           "test_linkage_clustering_wikipedia_example: UPGMA")
        call check_linkage(WIKI_N, WIKI_DISTANCES, METHOD_WEIGHTED, lower, higher, &
                           [17.0_real64, 22.0_real64, 28.0_real64, 35.0_real64], 0.0_real64, sizes, .true., &
                           "test_linkage_clustering_wikipedia_example: WPGMA")
        call check_linkage(WIKI_N, WIKI_DISTANCES, METHOD_WARD, lower, higher, &
                           [17.0_real64, sqrt(1651.0_real64/3.0_real64), 28.0_real64, sqrt(28876.0_real64/15.0_real64)], &
                           WARD_RTOL, sizes, .true., "test_linkage_clustering_wikipedia_example: Ward")
    end subroutine test_linkage_clustering_wikipedia_example

    !> Six points on a line at 0, 10, 40, 11, 60, 62. Every method merges {2, 4} at 1 (the cluster
    !| then sits at index 2), {5, 6} at 2 (at 5), {-1, 1} (at 1), {-2, 3} (at 3), {-4, -3}; sizes
    !| 2, 2, 3, 3, 6.
    !| - UPGMA, the mean of all distances between the two clusters: d(24, 1) = 10.5,
    !|   d(56, 3) = (20 + 22)/2 = 21, and {1, 2, 4} against {3, 5, 6}: (40 + 60 + 62 + 30 + 50 + 52
    !|   + 29 + 49 + 51)/9 = 47.
    !| - WPGMA: 10.5 and 21 as well; d(124, 3) = (40 + 29.5)/2 = 34.75, d(124, 56) = (61 + 50.5)/2 =
    !|   55.75, the last (55.75 + 34.75)/2 = 45.25.
    !| - Ward on distances between points of a line is sqrt(2 nA nB/(nA + nB)) times the distance
    !|   of the centroids: {2, 4} (at 10.5) to {1}: sqrt(4/3)*10.5 = 7 sqrt(3); {5, 6} (at 61) to
    !|   {3}: sqrt(4/3)*21 = 14 sqrt(3), below {1, 2, 4} to {3}, sqrt(3/2)*33; and {1, 2, 4} (at 7)
    !|   to {3, 5, 6} (at 54): sqrt(3)*47.
    subroutine test_linkage_clustering_points_on_a_line()
        integer(int32), parameter :: lower(LINE_N - 1) = [2, 5, -1, -2, -4], higher(LINE_N - 1) = [4, 6, 1, 3, -3]
        integer(int32), parameter :: sizes(LINE_N - 1) = [2, 2, 3, 3, 6]
        real(real64) :: distances(LINE_N, LINE_N)
        integer(int32) :: i_row, i_col

        do i_col = 1, LINE_N
            do i_row = 1, LINE_N
                distances(i_row, i_col) = abs(LINE_POSITIONS(i_row) - LINE_POSITIONS(i_col))
            end do
        end do

        call check_linkage(LINE_N, distances, METHOD_AVERAGE, lower, higher, &
                           [1.0_real64, 2.0_real64, 10.5_real64, 21.0_real64, 47.0_real64], 0.0_real64, sizes, .true., &
                           "test_linkage_clustering_points_on_a_line: UPGMA")
        call check_linkage(LINE_N, distances, METHOD_WEIGHTED, lower, higher, &
                           [1.0_real64, 2.0_real64, 10.5_real64, 21.0_real64, 45.25_real64], 0.0_real64, sizes, .true., &
                           "test_linkage_clustering_points_on_a_line: WPGMA")
        call check_linkage(LINE_N, distances, METHOD_WARD, lower, higher, &
                           [1.0_real64, 2.0_real64, 7*sqrt(3.0_real64), 14*sqrt(3.0_real64), 47*sqrt(3.0_real64)], &
                           WARD_RTOL, sizes, .true., "test_linkage_clustering_points_on_a_line: Ward")
    end subroutine test_linkage_clustering_points_on_a_line

    !> Kept from the old suite (`test_linkage_methods` case 2): four points all at distance 1. Every
    !| merge is at 1 for every method (for Ward: any two clusters of a regular simplex with unit
    !| edges are at 1). Every merge is a tie, broken by the lowest column and row: {1, 2}, then
    !| the cluster at index 1 with 3, then with 4; sizes 2, 3, 4. Ward's updated distances are 1
    !| only to within rounding, which can break the ties either way, so for Ward only the heights
    !| are checked.
    subroutine test_linkage_clustering_equal_distances()
        integer(int32), parameter :: n_points = 4
        real(real64) :: distances(n_points, n_points)
        integer(int32) :: i_method

        call fill_equidistant(n_points, 1.0_real64, distances)
        do i_method = 1, N_METHODS
            call check_linkage(n_points, distances, METHODS(i_method), [1, -1, -2], [2, 3, 4], &
                               [1.0_real64, 1.0_real64, 1.0_real64], merge(WARD_RTOL, 0.0_real64, i_method == 3), &
                               [2, 3, 4], i_method /= 3, &
                               "test_linkage_clustering_equal_distances: "//trim(METHOD_NAMES(i_method)))
        end do
    end subroutine test_linkage_clustering_equal_distances

    !> Kept from the old suite (`test_linkage_methods` case 3): two points at 5 merge once, {1, 2}
    !| at 5, size 2, for every method.
    subroutine test_linkage_clustering_two_points()
        real(real64) :: distances(2, 2)
        integer(int32) :: i_method

        call fill_equidistant(2_int32, 5.0_real64, distances)
        do i_method = 1, N_METHODS
            call check_linkage(2_int32, distances, METHODS(i_method), [1], [2], [5.0_real64], &
                               0.0_real64, [2], .true., "test_linkage_clustering_two_points: "//trim(METHOD_NAMES(i_method)))
        end do
    end subroutine test_linkage_clustering_two_points

    !> Kept from the old suite (`test_linkage_methods` case 4, which checked only `ierr`): one point
    !| merges nothing; `ierr` is OK, the outputs (no elements) are not written at all, and the
    !| matrix [0] comes back.
    subroutine test_linkage_clustering_one_point()
        real(real64) :: distances(1, 1)
        integer(int32) :: i_method
        integer(int32) :: no_merges(0)
        real(real64) :: no_heights(0)

        distances = 0.0_real64
        do i_method = 1, N_METHODS
            call check_linkage(1_int32, distances, METHODS(i_method), no_merges, no_merges, &
                               no_heights, 0.0_real64, no_merges, .true., &
                               "test_linkage_clustering_one_point: "//trim(METHOD_NAMES(i_method)))
        end do
    end subroutine test_linkage_clustering_one_point

    !> Points 1 and 3 coincide (distance 0), and point 2 is at 4 from both. {1, 3} merges first,
    !| at 0, then {-1, 2}: UPGMA and WPGMA at (4 + 4)/2 = 4, Ward at sqrt((2*16 + 2*16 - 0)/3) =
    !| 8/sqrt(3); sizes 2, 3.
    subroutine test_linkage_clustering_duplicate_points()
        real(real64) :: distances(3, 3)

        distances = reshape([0.0_real64, 4.0_real64, 0.0_real64, 4.0_real64, 0.0_real64, 4.0_real64, &
                             0.0_real64, 4.0_real64, 0.0_real64], shape(distances))
        call check_linkage(3_int32, distances, METHOD_AVERAGE, [1, -1], [3, 2], [0.0_real64, 4.0_real64], 0.0_real64, &
                           [2, 3], .true., "test_linkage_clustering_duplicate_points: UPGMA")
        call check_linkage(3_int32, distances, METHOD_WEIGHTED, [1, -1], [3, 2], [0.0_real64, 4.0_real64], 0.0_real64, &
                           [2, 3], .true., "test_linkage_clustering_duplicate_points: WPGMA")
        call check_linkage(3_int32, distances, METHOD_WARD, [1, -1], [3, 2], [0.0_real64, 8.0_real64/sqrt(3.0_real64)], &
                           WARD_RTOL, [2, 3], .true., "test_linkage_clustering_duplicate_points: Ward")
    end subroutine test_linkage_clustering_duplicate_points

    !> Four points that all coincide: every distance is 0, so every height is exactly 0 for every
    !| method, Ward's square root included, and the ties go as in
    !| `test_linkage_clustering_equal_distances`: {1, 2}, {-1, 3}, {-2, 4}; sizes 2, 3, 4.
    subroutine test_linkage_clustering_all_points_equal()
        real(real64) :: distances(4, 4)
        integer(int32) :: i_method

        distances = 0.0_real64
        do i_method = 1, N_METHODS
            call check_linkage(4_int32, distances, METHODS(i_method), [1, -1, -2], [2, 3, 4], &
                               [0.0_real64, 0.0_real64, 0.0_real64], 0.0_real64, [2, 3, 4], .true., &
                               "test_linkage_clustering_all_points_equal: "//trim(METHOD_NAMES(i_method)))
        end do
    end subroutine test_linkage_clustering_all_points_equal

    !> 64 points all at distance 1, the equal-distance case at a size where the cluster labels and
    !| sizes run past a few: every height is 1 for every method; for UPGMA and WPGMA, whose averages
    !| of ones are exactly 1, the ties chain every point onto the first cluster -- merge 1 is
    !| {1, 2}, merge k > 1 is {-(k - 1), k + 1} -- and the size after merge k is k + 1.
    subroutine test_linkage_clustering_many_points()
        integer(int32), parameter :: n_points = 64
        integer(int32) :: lower(n_points - 1), higher(n_points - 1), sizes(n_points - 1), i_merge, i_method
        real(real64) :: heights(n_points - 1), distances(n_points, n_points)

        do i_merge = 1, n_points - 1
            lower(i_merge) = merge(1_int32, -(i_merge - 1_int32), i_merge == 1)
            higher(i_merge) = i_merge + 1_int32
            sizes(i_merge) = i_merge + 1_int32
        end do
        heights = 1.0_real64
        call fill_equidistant(n_points, 1.0_real64, distances)

        do i_method = 1, N_METHODS
            call check_linkage(n_points, distances, METHODS(i_method), lower, higher, heights, &
                               merge(WARD_RTOL, 0.0_real64, i_method == 3), sizes, i_method /= 3, &
                               "test_linkage_clustering_many_points: "//trim(METHOD_NAMES(i_method)))
        end do
    end subroutine test_linkage_clustering_many_points

    !> `method` is one of the three codes 0, 1, 2; -1 and 3 are refused, naming `method`.
    subroutine test_linkage_clustering_method_range()
        real(real64) :: distances(3, 3)

        call fill_equidistant(3_int32, 1.0_real64, distances)
        call check_refused(3_int32, distances, 9_int32, -1_int32, ERR_INVALID_INPUT, 7_int32, &
                           "test_linkage_clustering_method_range: -1")
        call check_refused(3_int32, distances, 9_int32, 3_int32, ERR_INVALID_INPUT, 7_int32, &
                           "test_linkage_clustering_method_range: 3")
    end subroutine test_linkage_clustering_method_range

    !> No points is empty input, a negative count invalid, both naming `n_points`. (For n_points
    !| = -2 the distance check reads (-2)**2 = 4 elements of a matrix that has none; the actual
    !| passed here has 9, so that read stays inside it.)
    subroutine test_linkage_clustering_dimensions()
        real(real64) :: distances(3, 3)

        call fill_equidistant(3_int32, 1.0_real64, distances)
        call check_refused(0_int32, distances, 9_int32, METHOD_AVERAGE, ERR_EMPTY_INPUT, 2_int32, &
                           "test_linkage_clustering_dimensions: n_points 0")
        call check_refused(-2_int32, distances, 9_int32, METHOD_AVERAGE, ERR_INVALID_INPUT, 2_int32, &
                           "test_linkage_clustering_dimensions: n_points -2")
    end subroutine test_linkage_clustering_dimensions

    !> The distance-matrix convention: a negative distance (kept from the old suite's
    !| `test_linkage_methods` case 6), a non-zero diagonal and an asymmetric matrix are each
    !| refused, naming `distances`, for every method.
    subroutine test_linkage_clustering_distance_matrix_checks()
        real(real64) :: distances(3, 3)
        integer(int32) :: i_method
        character(*), parameter :: prefix = "test_linkage_clustering_distance_matrix_checks: "

        do i_method = 1, N_METHODS
            call fill_equidistant(3_int32, 1.0_real64, distances)
            distances(1, 3) = -1.0_real64
            distances(3, 1) = -1.0_real64
            call check_refused(3_int32, distances, 9_int32, METHODS(i_method), ERR_INVALID_INPUT, 1_int32, &
                               prefix//trim(METHOD_NAMES(i_method))//": negative distance")
            call fill_equidistant(3_int32, 1.0_real64, distances)
            distances(2, 2) = 0.5_real64
            call check_refused(3_int32, distances, 9_int32, METHODS(i_method), ERR_INVALID_INPUT, 1_int32, &
                               prefix//trim(METHOD_NAMES(i_method))//": non-zero diagonal")
            call fill_equidistant(3_int32, 1.0_real64, distances)
            distances(3, 2) = 2.0_real64
            call check_refused(3_int32, distances, 9_int32, METHODS(i_method), ERR_INVALID_INPUT, 1_int32, &
                               prefix//trim(METHOD_NAMES(i_method))//": asymmetric")
        end do
    end subroutine test_linkage_clustering_distance_matrix_checks

    !> Kept from the old suite (`test_linkage_methods` case 5, NaN only): NaN, +Inf and -Inf, in
    !| one triangle only and in both, are refused naming `distances`, for every method.
    subroutine test_linkage_clustering_rejects_nan_and_inf()
        real(real64) :: distances(3, 3)
        integer(int32) :: i_method, i_value

        do i_method = 1, N_METHODS
            do i_value = 1, N_NON_FINITE
                call fill_equidistant(3_int32, 1.0_real64, distances)
                distances(3, 1) = non_finite_value(i_value)
                call check_refused(3_int32, distances, 9_int32, METHODS(i_method), ERR_NAN_INF, 1_int32, &
                                   "test_linkage_clustering_rejects_nan_and_inf: "//trim(METHOD_NAMES(i_method))//": "// &
                                   trim(NON_FINITE_NAMES(i_value))//" in the lower triangle")
                distances(1, 3) = non_finite_value(i_value)
                call check_refused(3_int32, distances, 9_int32, METHODS(i_method), ERR_NAN_INF, 1_int32, &
                                   "test_linkage_clustering_rejects_nan_and_inf: "//trim(METHOD_NAMES(i_method))//": "// &
                                   trim(NON_FINITE_NAMES(i_value))//" in both triangles")
            end do
        end do
    end subroutine test_linkage_clustering_rejects_nan_and_inf

    !> Which refusal wins when several apply: `n_points` before `method` before `distances`. A
    !| NaN matrix with method 3 names `method`; with no points as well, `n_points`.
    subroutine test_linkage_clustering_check_order()
        real(real64) :: distances(3, 3)

        call fill_equidistant(3_int32, 1.0_real64, distances)
        distances(2, 1) = non_finite_value(1_int32)
        call check_refused(3_int32, distances, 9_int32, 3_int32, ERR_INVALID_INPUT, 7_int32, &
                           "test_linkage_clustering_check_order: method 3 and NaN")
        call check_refused(0_int32, distances, 9_int32, 3_int32, ERR_EMPTY_INPUT, 2_int32, &
                           "test_linkage_clustering_check_order: no points, method 3 and NaN")
    end subroutine test_linkage_clustering_check_order

    !> BUG (#223, "Ward update"): Ward's update squares the distances, so three points at
    !| s = 2**-600 from each other, all normal numbers, get squares below the smallest subnormal.
    !| Every method merges {1, 2} at s, then {-1, 3} at s: for Ward, d(12, 3)**2 = (2 s**2 + 2 s**2
    !| - s**2)/3 = s**2. Today Ward's second height is 0.
    subroutine test_linkage_clustering_ward_squares_below_tiny()
        real(real64), parameter :: s = 2.0_real64**(-600)
        real(real64) :: distances(3, 3)

        call fill_equidistant(3_int32, s, distances)
        call check_linkage(3_int32, distances, METHOD_WARD, [1, -1], [2, 3], [s, s], WARD_RTOL, [2, 3], &
                           .true., "test_linkage_clustering_ward_squares_below_tiny")
    end subroutine test_linkage_clustering_ward_squares_below_tiny
end module mod_test_linkage_clustering
