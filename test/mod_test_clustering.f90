! filepath: test/mod_test_clustering.f90
!> Unit test suite for clustering routine.
module mod_test_clustering
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf  
    use tox_clustering  
    use tox_errors
    use mod_test_suite, only: test_case
    implicit none
    public
    

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests_clustering() result(all_tests)
        type(test_case),allocatable :: all_tests(:)

        allocate(all_tests(6))

        all_tests(1) = test_case("test_clustering_k_means_assign_cluster_helper", test_k_means_assign_cluster_helper)
        all_tests(2) = test_case("test_clustering_k_means_recompute_cluster_centroids_helper", test_k_means_recompute_cluster_centroids_helper)
        all_tests(3) = test_case("test_clustering_k_means_clustering", test_k_means_clustering)
        all_tests(4) = test_case("test_clustering_linkage_helpers_xPGMA", test_linkage_helpers_xPGMA)
        all_tests(5) = test_case("test_clustering_linkage_helpers_ward", test_linkage_helpers_ward)
        all_tests(6) = test_case("test_clustering_linkage_methods", test_linkage_methods)
    end function get_all_tests_clustering

    
    

    !> Test the Ward linkage helper functions with a known example.
    subroutine test_linkage_helpers_ward()
        integer(int32), parameter :: n_points = 5
        real(real64) :: dist(n_points, n_points)
        real(real64) :: expected_dist(n_points, n_points)
        integer(int32) :: cluster_sizes(n_points - 1)
        real(real64) :: INF, min_dist
        integer(int32) :: row_idx, col_idx, ierr

        call set_ok(ierr)
        INF = ieee_value(1.0_real64, ieee_positive_inf)

        ! -------------------------------
        ! Example from https://en.wikipedia.org/wiki/UPGMA 29.10.2025
        ! -------------------------------
        dist = reshape([ &
            0.0, 17.0, 21.0, 31.0, 23.0, &
            17.0, 0.0, 30.0, 34.0, 21.0, &
            21.0, 30.0, 0.0, 28.0, 39.0, &
            31.0, 34.0, 28.0, 0.0, 43.0, &
            23.0, 21.0, 39.0, 43.0, 0.0 &
        ], shape(dist), order=[2,1])
        cluster_sizes = [0, 0, 0, 0]

        expected_dist = reshape([ &
            1.0_real64,  17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
            INF,         -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
            28.24299323136035_real64, 30.0_real64,  0.0_real64, 28.0_real64, 39.0_real64, &
            36.26292872893749_real64, 34.0_real64, 28.0_real64,  0.0_real64, 43.0_real64, &
            23.459184413217212_real64, 21.0_real64, 39.0_real64, 43.0_real64,  0.0_real64 &
        ], shape(dist), order=[2,1])

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 1_int32, "test_linkage_helpers_ward: first merge, wrong min column index")
        call assert_equal_int(row_idx, 2_int32, "test_linkage_helpers_ward: first merge, wrong min row index")
        call assert_equal_real(min_dist, 17.0_real64, TOL, "test_linkage_helpers_ward: first merge, wrong min value")
        call merge_distances_ward_linkage_helper(dist, n_points, row_idx, col_idx, 1_int32, 1_int32, 1_int32, cluster_sizes, min_dist)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_ward: first merge")

        expected_dist = reshape([ &
            2.0_real64,  17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
            INF,         -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
            34.94519518713076_real64, 30.0_real64,  0.0_real64, 28.0_real64, 39.0_real64, &
            42.10898558106888_real64, 34.0_real64, 28.0_real64,  0.0_real64, 43.0_real64, &
                    INF, 21.0_real64,         INF,         INF, -1.0_real64 &
        ], shape(dist), order=[2,1])
        cluster_sizes = [2, 0, 0, 0]

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 1_int32, "test_linkage_helpers_ward: second merge, wrong min column index")
        call assert_equal_int(row_idx, 5_int32, "test_linkage_helpers_ward: second merge, wrong min row index")
        call assert_equal_real(min_dist, 23.459184413217212_real64, TOL, "test_linkage_helpers_ward: second merge, wrong min value")
        call merge_distances_ward_linkage_helper(dist, n_points, row_idx, col_idx, 1_int32, 2_int32, 2_int32, cluster_sizes, min_dist)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_ward: second merge")

        expected_dist = reshape([ &
             2.0_real64, 17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
                    INF, -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
            43.87558166755931_real64, 30.0_real64,  3.0_real64, 28.0_real64, 39.0_real64, &
                    INF, 34.0_real64,         INF, -1.0_real64, 43.0_real64, &
                    INF, 21.0_real64,         INF,         INF, -1.0_real64 &
        ], shape(dist), order=[2,1])
        cluster_sizes = [2, 3, 0, 0]

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 3_int32, "test_linkage_helpers_ward: third merge, wrong min column index")
        call assert_equal_int(row_idx, 4_int32, "test_linkage_helpers_ward: third merge, wrong min row index")
        call assert_equal_real(min_dist, 28.0_real64, TOL, "test_linkage_helpers_ward: third merge, wrong min value")
        call merge_distances_ward_linkage_helper(dist, n_points, row_idx, col_idx, 1_int32, 1_int32, 3_int32, cluster_sizes, min_dist)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_ward: third merge")

        expected_dist = reshape([ &
             4.0_real64, 17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
                    INF, -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
                    INF, 30.0_real64, -1.0_real64, 28.0_real64, 39.0_real64, &
                    INF, 34.0_real64,         INF, -1.0_real64, 43.0_real64, &
                    INF, 21.0_real64,         INF,         INF, -1.0_real64 &
        ], shape(dist), order=[2,1])
        cluster_sizes = [2, 3, 5, 0]

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 1_int32, "test_linkage_helpers_ward: last merge, wrong min column index")
        call assert_equal_int(row_idx, 3_int32, "test_linkage_helpers_ward: last merge, wrong min row index")
        call assert_equal_real(min_dist, 43.87558166755931_real64, TOL, "test_linkage_helpers_ward: last merge, wrong min value")
        call merge_distances_ward_linkage_helper(dist, n_points, row_idx, col_idx, 2_int32, 3_int32, 4_int32, cluster_sizes, min_dist)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_ward: last merge")
    end subroutine test_linkage_helpers_ward

    !> Test the xPGMA linkage helper functions with a known example.
    subroutine test_linkage_helpers_xPGMA()
        integer(int32), parameter :: n_points = 5
        real(real64) :: dist(n_points, n_points)
        real(real64) :: expected_dist(n_points, n_points)
        real(real64) :: INF, min_dist
        integer(int32) :: row_idx, col_idx, ierr

        call set_ok(ierr)
        INF = ieee_value(1.0_real64, ieee_positive_inf)

        ! -------------------------------
        ! Example from https://en.wikipedia.org/wiki/UPGMA 29.10.2025
        ! -------------------------------
        dist = reshape([ &
            0.0, 17.0, 21.0, 31.0, 23.0, &
            17.0, 0.0, 30.0, 34.0, 21.0, &
            21.0, 30.0, 0.0, 28.0, 39.0, &
            31.0, 34.0, 28.0, 0.0, 43.0, &
            23.0, 21.0, 39.0, 43.0, 0.0 &
        ], shape(dist), order=[2,1])

        expected_dist = reshape([ &
            1.0_real64,  17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
            INF,         -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
            25.5_real64, 30.0_real64,  0.0_real64, 28.0_real64, 39.0_real64, &
            32.5_real64, 34.0_real64, 28.0_real64,  0.0_real64, 43.0_real64, &
            22.0_real64, 21.0_real64, 39.0_real64, 43.0_real64,  0.0_real64 &
        ], shape(dist), order=[2,1])

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 1_int32, "test_linkage_helpers_xPGMA: first merge, wrong min column index")
        call assert_equal_int(row_idx, 2_int32, "test_linkage_helpers_xPGMA: first merge, wrong min row index")
        call assert_equal_real(min_dist, 17.0_real64, TOL, "test_linkage_helpers_xPGMA: first merge, wrong min value")
        call merge_distances_xPGMA_linkage_helper(dist, n_points, row_idx, col_idx, 1_int32, 1_int32, 1_int32)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_xPGMA: first merge")

        expected_dist = reshape([ &
            2.0_real64,  17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
            INF,         -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
            30.0_real64, 30.0_real64,  0.0_real64, 28.0_real64, 39.0_real64, &
            36.0_real64, 34.0_real64, 28.0_real64,  0.0_real64, 43.0_real64, &
                    INF, 21.0_real64,         INF,         INF, -1.0_real64 &
        ], shape(dist), order=[2,1])

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 1_int32, "test_linkage_helpers_xPGMA: second merge, wrong min column index")
        call assert_equal_int(row_idx, 5_int32, "test_linkage_helpers_xPGMA: second merge, wrong min row index")
        call assert_equal_real(min_dist, 22.0_real64, TOL, "test_linkage_helpers_xPGMA: second merge, wrong min value")
        call merge_distances_xPGMA_linkage_helper(dist, n_points, row_idx, col_idx, 1_int32, 2_int32, 2_int32)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_xPGMA: second merge")

        expected_dist = reshape([ &
             2.0_real64, 17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
                    INF, -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
            33.0_real64, 30.0_real64,  3.0_real64, 28.0_real64, 39.0_real64, &
                    INF, 34.0_real64,         INF, -1.0_real64, 43.0_real64, &
                    INF, 21.0_real64,         INF,         INF, -1.0_real64 &
        ], shape(dist), order=[2,1])

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 3_int32, "test_linkage_helpers_xPGMA: third merge, wrong min column index")
        call assert_equal_int(row_idx, 4_int32, "test_linkage_helpers_xPGMA: third merge, wrong min row index")
        call assert_equal_real(min_dist, 28.0_real64, TOL, "test_linkage_helpers_xPGMA: third merge, wrong min value")
        call merge_distances_xPGMA_linkage_helper(dist, n_points, row_idx, col_idx, 1_int32, 1_int32, 3_int32)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_xPGMA: third merge")

        expected_dist = reshape([ &
             4.0_real64, 17.0_real64, 21.0_real64, 31.0_real64, 23.0_real64, &
                    INF, -1.0_real64, 30.0_real64, 34.0_real64, 21.0_real64, &
                    INF, 30.0_real64, -1.0_real64, 28.0_real64, 39.0_real64, &
                    INF, 34.0_real64,         INF, -1.0_real64, 43.0_real64, &
                    INF, 21.0_real64,         INF,         INF, -1.0_real64 &
        ], shape(dist), order=[2,1])

        call get_min_distance_indices_helper(dist, n_points, row_idx, col_idx, min_dist)
        call assert_equal_int(col_idx, 1_int32, "test_linkage_helpers_xPGMA: last merge, wrong min column index")
        call assert_equal_int(row_idx, 3_int32, "test_linkage_helpers_xPGMA: last merge, wrong min row index")
        call assert_equal_real(min_dist, 33.0_real64, TOL, "test_linkage_helpers_xPGMA: last merge, wrong min value")
        call merge_distances_xPGMA_linkage_helper(dist, n_points, row_idx, col_idx, 2_int32, 3_int32, 4_int32)
        call assert_equal_array_real(dist, expected_dist, n_points ** 2, TOL, "test_linkage_helpers_xPGMA: last merge")
    end subroutine test_linkage_helpers_xPGMA

    !> Test the k-means assign_cluster_helper function with a simple example.
    subroutine test_linkage_methods
        integer(int32), parameter :: max_n_points = 5
        integer(int32) :: n_points
        real(real64), allocatable :: orig_dist(:, :), passed_dist(:, :)
        integer(int32) :: merge_i(max_n_points - 1), merge_j(max_n_points - 1)
        real(real64) :: heights(max_n_points - 1)
        integer(int32) :: cluster_sizes(max_n_points - 1)
        integer(int32), dimension(max_n_points - 1, 3) :: expected_merge_i
        integer(int32), dimension(max_n_points - 1, 3) :: expected_merge_j
        real(real64), dimension(max_n_points - 1, 3) :: expected_heights
        integer(int32), dimension(max_n_points - 1, 3) :: expected_cluster_sizes
        integer(int32), dimension(3), parameter :: methods = [METHOD_AVERAGE, METHOD_WEIGHTED, METHOD_WARD]
        character(len=5), dimension(3), parameter :: method_names = ["UPGMA", "WPGMA", "Ward "]
        character(len=:), allocatable :: method_name
        integer(int32) :: ierr, i, i_method

        call set_ok(ierr)

        expected_merge_i = reshape([ &
            1, -1,3, -3, & ! UPGMA
            1, -1,3, -3, & ! WPGMA
            1, -1,3, -3 & ! ward
        ], shape(expected_merge_i))
        expected_merge_j = reshape([ &
            2, 5, 4, -2, & ! UPGMA
            2, 5, 4, -2, & ! WPGMA
            2, 5, 4, -2 & ! ward
        ], shape(expected_merge_j))
        expected_heights = reshape([ &
            17.0_real64, 22.0_real64, 28.0_real64, 33.0_real64, & ! UPGMA
            17.0_real64, 22.0_real64, 28.0_real64, 35.0_real64, & ! WPGMA
            17.0_real64, 23.459184413217212_real64, 28.0_real64, 43.87558166755931_real64 & ! ward
        ], shape(expected_heights))
        expected_cluster_sizes = reshape([ &
            2, 3, 2, 5, & ! UPGMA
            2, 3, 2, 5, & ! WPGMA
            2, 3, 2, 5 & ! ward
        ], shape(expected_cluster_sizes))
        allocate(orig_dist(5,5), passed_dist(5,5))

        do i_method = 1, size(methods)
            method_name = trim(method_names(i_method))
            ! -------------------------------
            ! Case 1: Example from https://en.wikipedia.org/wiki/WPGMA 29.10.2025
            ! -------------------------------
            n_points = 5
            deallocate(orig_dist, passed_dist)
            allocate(orig_dist(5,5), passed_dist(5,5))
            orig_dist = reshape([ &
                0.0, 17.0, 21.0, 31.0, 23.0, &
                17.0, 0.0, 30.0, 34.0, 21.0, &
                21.0, 30.0, 0.0, 28.0, 39.0, &
                31.0, 34.0, 28.0, 0.0, 43.0, &
                23.0, 21.0, 39.0, 43.0, 0.0 &
            ], [5,5])
            passed_dist = orig_dist

            call linkage_clustering(passed_dist, n_points, merge_i, merge_j, heights, cluster_sizes, methods(i_method), ierr)
            call assert_equal_int(ierr, ERR_OK, "test_linkage_methods: "//method_name//": reference case ierr")
            call assert_equal_array_real(passed_dist, orig_dist, size(orig_dist, kind=int32), 0.0_real64, "test_linkage_methods: "//method_name//": reference output matrix doesn't match input matrix")
            do i = 1, n_points - 1
                call assert_equal_int(min(merge_i(i), merge_j(i)), expected_merge_i(i, i_method), "test_linkage_methods: "//method_name//": reference merge_i")
                call assert_equal_int(max(merge_i(i), merge_j(i)), expected_merge_j(i, i_method), "test_linkage_methods: "//method_name//": reference merge_j")
                call assert_not_equal_int(merge_i(i), merge_j(i), "test_linkage_methods: "//method_name//": reference merge_j")
            end do
            call assert_equal_array_real(heights, expected_heights(:, i_method), size(expected_heights, 1, kind=int32), TOL, "test_linkage_methods: "//method_name//": reference heights")
            call assert_equal_array_int(cluster_sizes, expected_cluster_sizes(:, i_method), size(expected_cluster_sizes, 1, kind=int32), "test_linkage_methods: "//method_name//": reference cluster_sizes")

            ! -------------------------------
            ! Case 2: Equal distances
            ! -------------------------------
            n_points = 4
            expected_heights(:n_points-1, i_method) = [1.0_real64, 1.0_real64, 1.0_real64]

            deallocate(orig_dist, passed_dist)
            allocate(orig_dist(4,4), passed_dist(4,4))
            orig_dist = reshape([ &
                0.0, 1.0, 1.0, 1.0, &
                1.0, 0.0, 1.0, 1.0, &
                1.0, 1.0, 0.0, 1.0, &
                1.0, 1.0, 1.0, 0.0 &
            ], [4,4])
            passed_dist = orig_dist

            call linkage_clustering(passed_dist, n_points, merge_i, merge_j, heights, cluster_sizes, methods(i_method), ierr)
            call assert_equal_int(ierr, ERR_OK, "test_linkage_methods: "//method_name//": equal-distance case ierr")
            call assert_equal_array_real(passed_dist, orig_dist, size(orig_dist, kind=int32), 0.0_real64, "test_linkage_methods: "//method_name//": equal-distance output matrix doesn't match input matrix")
            call assert_equal_array_real(heights, expected_heights(:, i_method), size(expected_heights, 1, kind=int32), TOL, "test_linkage_methods: "//method_name//": equal-distance heights")

            ! -------------------------------
            ! Case 3: Two points
            ! -------------------------------
            n_points = 2
            expected_heights(:n_points-1, i_method) = [5.0_real64]

            deallocate(orig_dist, passed_dist)
            allocate(orig_dist(2,2), passed_dist(2,2))
            orig_dist = reshape([ &
                0.0, 5.0, &
                5.0, 0.0 &
            ], [2,2])
            passed_dist = orig_dist

            call linkage_clustering(passed_dist, n_points, merge_i, merge_j, heights, cluster_sizes, methods(i_method), ierr)
            call assert_equal_int(ierr, ERR_OK, "test_linkage_methods: "//method_name//": two-point case ierr")
            call assert_equal_array_real(passed_dist, orig_dist, size(orig_dist, kind=int32), 0.0_real64, "test_linkage_methods: "//method_name//": two-point output matrix doesn't match input matrix")
            call assert_equal_array_real(heights, expected_heights(:, i_method), size(expected_heights, 1, kind=int32), TOL, "test_linkage_methods: "//method_name//": two-point height")

            ! -------------------------------
            ! Case 4: Single point
            ! -------------------------------
            n_points = 1

            deallocate(orig_dist, passed_dist)
            allocate(orig_dist(1,1), passed_dist(1,1))
            orig_dist = reshape([ &
                0.0 &
            ], [1,1])
            passed_dist = orig_dist

            call linkage_clustering(passed_dist, n_points, merge_i, merge_j, heights, cluster_sizes, methods(i_method), ierr)
            call assert_equal_int(ierr, ERR_OK, "test_linkage_methods: "//method_name//": single-point case ierr")

            ! -------------------------------
            ! Case 5: NaN in distance matrix
            ! -------------------------------
            n_points = 3

            deallocate(orig_dist, passed_dist)
            allocate(orig_dist(3,3), passed_dist(3,3))
            orig_dist = reshape([ &
                                            0.0_real64, 1.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), &
                                            1.0_real64, 0.0_real64, 1.0_real64, &
                ieee_value(1.0_real64, ieee_quiet_nan), 1.0_real64, 0.0_real64 &
            ], [3,3])
            passed_dist = orig_dist

            call linkage_clustering(passed_dist, n_points, merge_i, merge_j, heights, cluster_sizes, methods(i_method), ierr)
            call assert_equal_int(ierr, ERR_NAN_INF, "test_linkage_methods: "//method_name//": NaN case should trigger ERR_NAN_INF")
            call assert_equal_array_real(passed_dist, orig_dist, size(orig_dist, kind=int32), 0.0_real64, "test_linkage_methods: "//method_name//": NaN case should output matrix doesn't match input matrix")

            ! -------------------------------
            ! Case 6: Negative value in distance matrix
            ! -------------------------------
            n_points = 3

            deallocate(orig_dist, passed_dist)
            allocate(orig_dist(3,3), passed_dist(3,3))
            orig_dist = reshape([ &
                0.0_real64, 1.0_real64, -1.0_real64, &
                1.0_real64, 0.0_real64, 1.0_real64, &
                -1.0_real64, 1.0_real64, 0.0_real64 &
            ], [3,3])
            passed_dist = orig_dist

            call linkage_clustering(passed_dist, n_points, merge_i, merge_j, heights, cluster_sizes, methods(i_method), ierr)
            call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_linkage_methods: "//method_name//": Negative distance case should trigger ERR_INVALID_INPUT")
            call assert_equal_array_real(passed_dist, orig_dist, size(orig_dist, kind=int32), 0.0_real64, "test_linkage_methods: "//method_name//": Negative distance case should output matrix doesn't match input matrix")
        end do
    end subroutine test_linkage_methods

    !> Test the k-means clustering function with a simple example.
    subroutine test_k_means_clustering()
        integer(int32), parameter :: n_dims = 2, n_points = 4, n_clusters = 2
        real(real64) :: data_points(n_dims, n_points)
        real(real64) :: centroids(n_dims, n_clusters), expected_centroids(n_dims, n_clusters)
        integer(int32) :: labels(n_points), label_counts(n_clusters), ierr

        ! Define 2D points
        ! Cluster 1: (1,1), (2,2)
        ! Cluster 2: (8,8), (10,10)
        data_points(:,1) = [1.0_real64, 1.0_real64]
        data_points(:,2) = [2.0_real64, 2.0_real64]
        data_points(:,3) = [8.0_real64, 8.0_real64]
        data_points(:,4) = [10.0_real64, 10.0_real64]

        ! Initialize centroids to first and last points
        centroids(:,1) = data_points(:,1)
        centroids(:,2) = data_points(:,2)

        ! Expected final centroids
        expected_centroids(:,1) = [1.5_real64, 1.5_real64]
        expected_centroids(:,2) = [9.0_real64, 9.0_real64]

        ! Call clustering routine
        call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_k_means_clustering: expected OK status")

        ! Validate final centroids
        call assert_equal_array_real(centroids(:,1), expected_centroids(:,1), n_dims, TOL, "test_k_means_clustering: centroid(:,1) mismatch")
        call assert_equal_array_real(centroids(:,2), expected_centroids(:,2), n_dims, TOL, "test_k_means_clustering: centroid(:,2) mismatch")

        ! Validate label assignments
        call assert_equal_int(labels(1), 1_int32, "test_k_means_clustering: labels(1) should be 1")
        call assert_equal_int(labels(2), 1_int32, "test_k_means_clustering: labels(2) should be 1")
        call assert_equal_int(labels(3), 2_int32, "test_k_means_clustering: labels(3) should be 2")
        call assert_equal_int(labels(4), 2_int32, "test_k_means_clustering: labels(4) should be 2")

        ! Validate label counts
        call assert_equal_int(label_counts(1), 2_int32, "test_k_means_clustering: label_counts(1) mismatch")
        call assert_equal_int(label_counts(2), 2_int32, "test_k_means_clustering: label_counts(2) mismatch")

        ! Test max_iterations=1, should be different -> not clusters [1,2] [3,4], instead [1] [2,3,4]
        ! Initialize centroids to first and last points
        centroids(:,1) = data_points(:,1)
        centroids(:,2) = data_points(:,2)

        ! Expected final centroids
        expected_centroids(:,1) = [1.0_real64, 1.0_real64]
        expected_centroids(:,2) = [20.0_real64 / 3, 20.0_real64 / 3]

        ! Call clustering routine
        call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr, 1_int32)
        call assert_equal_int(ierr, ERR_OK, "test_k_means_clustering: expected OK status")

        ! Validate final centroids
        call assert_equal_array_real(centroids(:,1), expected_centroids(:,1), n_dims, TOL, "test_k_means_clustering: centroid(:,1) mismatch")
        call assert_equal_array_real(centroids(:,2), expected_centroids(:,2), n_dims, TOL, "test_k_means_clustering: centroid(:,2) mismatch")

        ! Validate label assignments
        call assert_equal_int(labels(1), 1, "test_k_means_clustering: labels(1) should be 1")
        call assert_equal_int(labels(2), 2, "test_k_means_clustering: labels(2) should be 1")
        call assert_equal_int(labels(3), 2, "test_k_means_clustering: labels(3) should be 2")
        call assert_equal_int(labels(4), 2, "test_k_means_clustering: labels(4) should be 2")

        ! Validate label counts
        call assert_equal_int(label_counts(1), 1, "test_k_means_clustering: label_counts(1) mismatch")
        call assert_equal_int(label_counts(2), 3, "test_k_means_clustering: label_counts(2) mismatch")
    end subroutine test_k_means_clustering

    !> Test the k-means recompute_cluster_centroids_helper function with a simple example.
    subroutine test_k_means_recompute_cluster_centroids_helper()
        integer(int32), parameter :: n_dims = 2, n_points = 4, n_clusters = 3
        real(real64) :: data_points(n_dims, n_points)
        real(real64) :: centroids(n_dims, n_clusters), expected_centroids(n_dims, n_clusters)
        integer(int32) :: labels(n_points), label_counts(n_clusters)

        ! Define 2D points
        ! Cluster 1: (1,1), (2,2)
        ! Cluster 2: (8,8), (10,10)
        ! Cluster 3: no points assigned
        data_points(:,1) = [1.0_real64, 1.0_real64]
        data_points(:,2) = [2.0_real64, 2.0_real64]
        data_points(:,3) = [8.0_real64, 8.0_real64]
        data_points(:,4) = [10.0_real64, 10.0_real64]

        ! Assign labels manually
        labels = [1, 1, 2, 2]
        label_counts = [2, 2, 0]

        ! Expected centroids:
        expected_centroids(:,1) = [1.5_real64, 1.5_real64]  ! mean of (1,1) and (2,2)
        expected_centroids(:,2) = [9.0_real64, 9.0_real64]  ! mean of (8,8) and (10,10)
        expected_centroids(:,3) = [0.0_real64, 0.0_real64]  ! no points assigned

        ! Initialize centroids to sentinel
        centroids = -999.0_real64

        ! Call routine
        call k_means_recompute_cluster_centroids_helper(data_points, n_points, n_dims, centroids, n_clusters, labels, label_counts)

        ! Validate centroids
        call assert_equal_array_real(centroids(:,1), expected_centroids(:,1), n_dims, TOL, "test_k_means_recompute_cluster_centroids_helper: centroid(:,1) mismatch")
        call assert_equal_array_real(centroids(:,2), expected_centroids(:,2), n_dims, TOL, "test_k_means_recompute_cluster_centroids_helper: centroid(:,2) mismatch")
        call assert_equal_array_real(centroids(:,3), expected_centroids(:,3), n_dims, TOL, "test_k_means_recompute_cluster_centroids_helper: centroid(:,3) should be zero")

        ! Validate label counts
        call assert_equal_int(label_counts(1), 2, "test_k_means_recompute_cluster_centroids_helper: label_counts(1) mismatch")
        call assert_equal_int(label_counts(2), 2, "test_k_means_recompute_cluster_centroids_helper: label_counts(2) mismatch")
        call assert_equal_int(label_counts(3), 0, "test_k_means_recompute_cluster_centroids_helper: label_counts(3) should be zero")

    end subroutine test_k_means_recompute_cluster_centroids_helper

    !> Test the k-means assign_cluster_helper function with a simple example.
    subroutine test_k_means_assign_cluster_helper()
        integer(int32), parameter :: n_dims = 2, n_points = 3, n_clusters = 2
        real(real64) :: data_points(n_dims, n_points)
        real(real64) :: centroids(n_dims, n_clusters)
        integer(int32) :: i_point, label, expected_label

        ! Define 2D points
        ! Point 1: (1.0, 1.0)
        ! Point 2: (5.0, 5.0)
        ! Point 3: (9.0, 9.0)
        data_points(:,1) = [1.0_real64, 1.0_real64]
        data_points(:,2) = [5.0_real64, 5.0_real64]
        data_points(:,3) = [9.0_real64, 9.0_real64]

        ! Define centroids
        ! Cluster 1: (0.0, 0.0)
        ! Cluster 2: (10.0, 10.0)
        centroids(:,1) = [0.0_real64, 0.0_real64]
        centroids(:,2) = [10.0_real64, 10.0_real64]

        ! Point 1 should be closer to Cluster 1
        i_point = 1
        expected_label = 1
        call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        call assert_equal_int(label, expected_label, "test_k_means_assign_cluster_helper: Point 1 should be assigned to Cluster 1")

        ! Point 2 has same distance to both -> should be assigned to last
        i_point = 2
        expected_label = 2
        call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        call assert_equal_int(label, expected_label, "test_k_means_assign_cluster_helper: Point 2 should be assigned to Cluster 2")

        ! Point 3 should be closer to Cluster 2
        i_point = 3
        expected_label = 2
        call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        call assert_equal_int(label, expected_label, "test_k_means_assign_cluster_helper: Point 3 should be assigned to Cluster 2")

    end subroutine test_k_means_assign_cluster_helper

    
end module mod_test_clustering
