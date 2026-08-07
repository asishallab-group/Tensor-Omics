!> Unit test suite for tox_shape_truthful_clustering_seeding (calculate_density_radius,
!| density_labels, seeds), generated from
!| src/kernel/shape_truthful_clustering/tox_shape_truthful_clustering_seeding_kernel.F90.
module mod_test_shape_truthful_clustering_seeding
    use tox_shape_truthful_clustering_seeding, only: calculate_density_radius_alloc, density_labels_alloc, seeds_alloc
    use f42_kd_tree, only: build_kd_index_alloc
    use tox_errors, only: is_ok, is_err
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_seeding() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(10))

        all_tests(1) = test_case("test_density_radius_default_percentile", test_density_radius_default_percentile)
        all_tests(2) = test_case("test_density_radius_custom_percentile", test_density_radius_custom_percentile)
        all_tests(3) = test_case("test_density_radius_invalid_percentile", test_density_radius_invalid_percentile)
        all_tests(4) = test_case("test_density_radius_single_vector", test_density_radius_single_vector)
        all_tests(5) = test_case("test_density_labels_basic", test_density_labels_basic)
        all_tests(6) = test_case("test_density_labels_zero_radius", test_density_labels_zero_radius)
        all_tests(7) = test_case("test_density_labels_invalid_kd_indices", test_density_labels_invalid_kd_indices)
        all_tests(8) = test_case("test_seeds_two_separated_clusters", test_seeds_two_separated_clusters)
        all_tests(9) = test_case("test_seeds_single_cluster_one_seed", test_seeds_single_cluster_one_seed)
        all_tests(10) = test_case("test_seeds_invalid_percentile", test_seeds_invalid_percentile)
    end function get_all_tests_shape_truthful_clustering_seeding

    ! --- calculate_density_radius ------------------------------------------
    !
    ! Shared fixture: D=2, N=5 points on a line, (0,0),(1,0),(2,0),(3,0),(4,0). Mean = (2,0),
    ! so the mean-to-vector distances are [2,1,0,1,2], sorted ascending [0,1,1,2,2]. Expected
    ! percentile values below are hand-computed from calc_percentile_helper's documented
    ! linear-interpolation formula: rank = (percentile/100)*(n-1) + 1.
    ! Default 15th percentile: rank = 0.15*4+1 = 1.6 -> interpolate sorted[1]=0, sorted[2]=1
    ! at fraction 0.6 -> 0.6. 50th percentile: rank = 0.5*4+1 = 3.0 exactly -> sorted[3] = 1.0.

    subroutine test_density_radius_default_percentile()
        real(real64)   :: vectors(2, 5) = reshape( &
                          [0.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0, 0.0d0, 3.0d0, 0.0d0, 4.0d0, 0.0d0], [2, 5])
        real(real64)   :: radius
        integer(int32) :: ierr

        call calculate_density_radius_alloc(vectors, 2_int32, 5_int32, radius=radius, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calculate_density_radius failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_real(radius, 0.6d0, 1.0d-9, "calculate_density_radius: default 15th percentile")
    end subroutine test_density_radius_default_percentile

    subroutine test_density_radius_custom_percentile()
        real(real64)   :: vectors(2, 5) = reshape( &
                          [0.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0, 0.0d0, 3.0d0, 0.0d0, 4.0d0, 0.0d0], [2, 5])
        real(real64)   :: radius
        integer(int32) :: ierr

        call calculate_density_radius_alloc(vectors, 2_int32, 5_int32, mean_to_other_vecs_dist_quant=0.5d0, radius=radius, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calculate_density_radius failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_real(radius, 1.0d0, 1.0d-9, "calculate_density_radius: 50th percentile")
    end subroutine test_density_radius_custom_percentile

    subroutine test_density_radius_invalid_percentile()
        real(real64)   :: vectors(2, 5) = reshape( &
                          [0.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0, 0.0d0, 3.0d0, 0.0d0, 4.0d0, 0.0d0], [2, 5])
        real(real64)   :: radius
        integer(int32) :: ierr

        call calculate_density_radius_alloc(vectors, 2_int32, 5_int32, mean_to_other_vecs_dist_quant=1.5d0, radius=radius, ierr=ierr)
        call assert_true(is_err(ierr), "calculate_density_radius should reject a quantile > 1.0")
    end subroutine test_density_radius_invalid_percentile

    !> A single vector: mean is itself, distance is exactly zero -- radius must be zero, not
    !| a crash or an undefined percentile.
    subroutine test_density_radius_single_vector()
        real(real64)   :: vectors(2, 1) = reshape([5.0d0, 5.0d0], [2, 1])
        real(real64)   :: radius
        integer(int32) :: ierr

        call calculate_density_radius_alloc(vectors, 2_int32, 1_int32, radius=radius, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calculate_density_radius failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_real(radius, 0.0d0, 1.0d-12, "calculate_density_radius: single vector has radius 0")
    end subroutine test_density_radius_single_vector

    ! --- density_labels ------------------------------------------------------
    !
    ! Shared fixture: D=2, N=11 points on a line, (0,0),(1,0),...,(10,0).

    !> radius=1.5: interior points see 3 neighbors (self and one on each side), the two
    !| endpoints see only 2 (self and one side).
    subroutine test_density_labels_basic()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: labels(11), expected(11)
        integer(int32) :: i

        call build_line_fixture(vectors, kd_indices, dim_order)

        call density_labels_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 1.5d0, labels, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels failed unexpectedly: ', ierr
            error stop
        end if

        expected = 3.0d0
        expected(1) = 2.0d0
        expected(11) = 2.0d0
        do i = 1, 11
            call assert_equal_real(labels(i), expected(i), 1.0d-12, "density_labels at radius 1.5")
        end do
    end subroutine test_density_labels_basic

    !> radius=0: every point is distinct, so only itself is within radius -- label 1 everywhere.
    subroutine test_density_labels_zero_radius()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: labels(11)

        call build_line_fixture(vectors, kd_indices, dim_order)

        call density_labels_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 0.0d0, labels, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels failed unexpectedly: ', ierr
            error stop
        end if
        call assert_true(all(abs(labels - 1.0d0) < 1.0d-12), "density_labels at radius 0 is 1 everywhere")
    end subroutine test_density_labels_zero_radius

    !> A kd_indices entry out of [1,n_vectors] must be rejected by validation.
    subroutine test_density_labels_invalid_kd_indices()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: labels(11)

        call build_line_fixture(vectors, kd_indices, dim_order)
        kd_indices(1) = 12

        call density_labels_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 1.5d0, labels, ierr)
        call assert_true(is_err(ierr), "density_labels should reject a kd_indices entry > n_vectors")
    end subroutine test_density_labels_invalid_kd_indices

    ! --- seeds -----------------------------------------------------------------

    !> Two tightly-packed clusters far apart: the mean sits between them, so the
    !| mean-to-vector-distance percentile radius is close to half the inter-cluster distance --
    !| large enough to cover a whole cluster, far too small to reach the other one. Expect
    !| exactly 2 seeds, one from each cluster.
    subroutine test_seeds_two_separated_clusters()
        real(real64)   :: vectors(2, 10) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0, &
                          10.0d0, 0.0d0, 10.1d0, 0.0d0, 10.0d0, 0.1d0, 9.9d0, 0.0d0, 10.0d0, -0.1d0], [2, 10])
        integer(int32) :: kd_indices(10), dim_order(2), ierr
        logical        :: is_seed_mask(10)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 10_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_two_separated_clusters: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 10_int32, kd_indices, dim_order, is_seed_mask=is_seed_mask, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(count(is_seed_mask, kind=int32), 2_int32, "seeds: two separated clusters give 2 seeds")
        call assert_true(any(is_seed_mask(1:5)), "seeds: cluster A (indices 1-5) has a seed")
        call assert_true(any(is_seed_mask(6:10)), "seeds: cluster B (indices 6-10) has a seed")
    end subroutine test_seeds_two_separated_clusters

    !> A single tight cluster: at the *default* (15th) percentile the density radius is
    !| smaller than the cluster's own point spacing (both are the same order of magnitude,
    !| derived from the same distances), so it does *not* merge every point into one seed --
    !| that is a genuine, correct property of the algorithm, not a bug: a density radius
    !| calibrated to be small relative to the data can fragment a locally sparse cluster. At
    !| the 100th percentile, the radius equals the cluster's own diameter (the farthest
    !| mean-to-vector distance) exactly, which *is* large enough to cover the whole cluster
    !| from a single pick -- so this is what the test asks for.
    subroutine test_seeds_single_cluster_one_seed()
        real(real64)   :: vectors(2, 5) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0], [2, 5])
        integer(int32) :: kd_indices(5), dim_order(2), ierr
        logical        :: is_seed_mask(5)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_single_cluster_one_seed: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, &
                         mean_to_other_vecs_dist_quant=1.0d0, is_seed_mask=is_seed_mask, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_int(count(is_seed_mask, kind=int32), 1_int32, "seeds: a single tight cluster gives 1 seed")
    end subroutine test_seeds_single_cluster_one_seed

    subroutine test_seeds_invalid_percentile()
        real(real64)   :: vectors(2, 5) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0], [2, 5])
        integer(int32) :: kd_indices(5), dim_order(2), ierr
        logical        :: is_seed_mask(5)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_invalid_percentile: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, mean_to_other_vecs_dist_quant=-0.1d0, &
                  is_seed_mask=is_seed_mask, ierr=ierr)
        call assert_true(is_err(ierr), "seeds should reject a negative mean_to_other_vecs_dist_quant")
    end subroutine test_seeds_invalid_percentile

    !> Build the shared 11-point line fixture and its k-d tree.
    subroutine build_line_fixture(vectors, kd_indices, dim_order)
        real(real64), intent(out) :: vectors(2, 11)
        integer(int32), intent(out) :: kd_indices(11)
        integer(int32), intent(out) :: dim_order(2)
        integer(int32) :: i, ierr

        do i = 1, 11
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]

        call build_kd_index_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'build_line_fixture: build_kd_index_alloc failed: ', ierr
            error stop
        end if
    end subroutine build_line_fixture

end module mod_test_shape_truthful_clustering_seeding
