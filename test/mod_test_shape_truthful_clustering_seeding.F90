!> Unit test suite for tox_shape_truthful_clustering_seeding (density_labels, seeds),
!| generated from
!| src/kernel/shape_truthful_clustering/tox_shape_truthful_clustering_seeding_kernel.F90.
module mod_test_shape_truthful_clustering_seeding
    use tox_shape_truthful_clustering_seeding, only: density_labels_alloc, seeds_alloc
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
        allocate (all_tests(15))

        all_tests(1) = test_case("test_density_labels_hand_computed", test_density_labels_hand_computed)
        all_tests(2) = test_case("test_density_labels_bandwidth_percentile_median", &
                                 test_density_labels_bandwidth_percentile_median)
        all_tests(3) = test_case("test_density_labels_invalid_bandwidth_percentile", &
                                 test_density_labels_invalid_bandwidth_percentile)
        all_tests(4) = test_case("test_density_labels_symmetric_neighborhood_does_not_underflow", &
                                 test_density_labels_symmetric_neighborhood_does_not_underflow)
        all_tests(5) = test_case("test_density_labels_uniform_interior_points_agree", &
                                 test_density_labels_uniform_interior_points_agree)
        all_tests(6) = test_case("test_density_labels_dense_vs_sparse", test_density_labels_dense_vs_sparse)
        all_tests(7) = test_case("test_density_labels_invalid_kd_indices", test_density_labels_invalid_kd_indices)
        all_tests(8) = test_case("test_density_labels_k_density_too_large", test_density_labels_k_density_too_large)
        all_tests(9) = test_case("test_density_labels_omitted_k_density_is_clamped", &
                                 test_density_labels_omitted_k_density_is_clamped)
        all_tests(10) = test_case("test_seeds_two_separated_clusters", test_seeds_two_separated_clusters)
        all_tests(11) = test_case("test_seeds_single_cluster_one_seed", test_seeds_single_cluster_one_seed)
        all_tests(12) = test_case("test_seeds_invalid_k_density", test_seeds_invalid_k_density)
        all_tests(13) = test_case("test_seeds_omitted_k_density_is_clamped", test_seeds_omitted_k_density_is_clamped)
        all_tests(14) = test_case("test_seeds_exclusion_radius_percentile_widens_coverage", &
                                  test_seeds_exclusion_radius_percentile_widens_coverage)
        all_tests(15) = test_case("test_seeds_invalid_exclusion_radius_percentile", &
                                  test_seeds_invalid_exclusion_radius_percentile)
    end function get_all_tests_shape_truthful_clustering_seeding

    ! --- density_labels ---------------------------------------------------------
    !
    ! 3 points on a line, (0,0),(1,0),(3,0), k_density=2 (every other point), default
    ! bandwidth_percentile=68.27. Hand-computed (and cross-checked against an independent
    ! Python re-implementation of the same formula, including calc_percentile_helper's own
    ! linear-interpolation rule: rank = (percentile/100)*(n-1) + 1) expected densities:
    !   point 1 (x=0): neighbor distances [1,3] -> rank = 0.6827*1+1 = 1.6827 -> interpolate
    !     1 and 3 at fraction 0.6827 -> bandwidth = 2.3654
    !     -> density = (exp(-1/(2*2.3654^2)) + exp(-9/(2*2.3654^2))) / 2.3654^2 ~= 0.2434133437
    !   point 2 (x=1): neighbor distances [1,2] -> rank = 1.6827 -> bandwidth = 1.6827
    !     -> density ~= 0.4702740525
    !   point 3 (x=3): neighbor distances [2,3] -> rank = 1.6827 -> bandwidth = 2.6827
    !     -> density ~= 0.1795903795

    subroutine test_density_labels_hand_computed()
        real(real64)   :: vectors(2, 3) = reshape([0.0d0, 0.0d0, 1.0d0, 0.0d0, 3.0d0, 0.0d0], [2, 3])
        integer(int32) :: kd_indices(3), dim_order(2), ierr
        real(real64)   :: labels(3)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 3_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_density_labels_hand_computed: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call density_labels_alloc(vectors, 2_int32, 3_int32, kd_indices, dim_order, k_density=2_int32, &
                                  labels=labels, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_real(labels(1), 0.2434133437d0, 1.0d-6, "density_labels: point at x=0")
        call assert_equal_real(labels(2), 0.4702740525d0, 1.0d-6, "density_labels: point at x=1")
        call assert_equal_real(labels(3), 0.1795903795d0, 1.0d-6, "density_labels: point at x=3")
    end subroutine test_density_labels_hand_computed

    !> Same fixture as test_density_labels_hand_computed, but bandwidth_percentile=50.0 (the
    !| median) instead of the default 68.27 -- confirms the parameter is actually wired
    !| through end to end (not just accepted and ignored), via the same hand-computation
    !| approach: rank = (50/100)*(2-1)+1 = 1.5 -> interpolate at fraction 0.5.
    !|   point 1 (x=0): neighbor distances [1,3] -> bandwidth = 2.0
    !|     -> density = (exp(-1/8) + exp(-9/8)) / 4 ~= 0.3017873425
    !|   point 2 (x=1): neighbor distances [1,2] -> bandwidth = 1.5 -> density ~= 0.5385998637
    !|   point 3 (x=3): neighbor distances [2,3] -> bandwidth = 2.5 -> density ~= 0.1940642069
    subroutine test_density_labels_bandwidth_percentile_median()
        real(real64)   :: vectors(2, 3) = reshape([0.0d0, 0.0d0, 1.0d0, 0.0d0, 3.0d0, 0.0d0], [2, 3])
        integer(int32) :: kd_indices(3), dim_order(2), ierr
        real(real64)   :: labels(3)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 3_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_density_labels_bandwidth_percentile_median: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call density_labels_alloc(vectors, 2_int32, 3_int32, kd_indices, dim_order, k_density=2_int32, &
                                  bandwidth_percentile=50.0d0, labels=labels, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_real(labels(1), 0.3017873425d0, 1.0d-6, "density_labels: median bandwidth, point at x=0")
        call assert_equal_real(labels(2), 0.5385998637d0, 1.0d-6, "density_labels: median bandwidth, point at x=1")
        call assert_equal_real(labels(3), 0.1940642069d0, 1.0d-6, "density_labels: median bandwidth, point at x=3")
    end subroutine test_density_labels_bandwidth_percentile_median

    subroutine test_density_labels_invalid_bandwidth_percentile()
        real(real64)   :: vectors(2, 3) = reshape([0.0d0, 0.0d0, 1.0d0, 0.0d0, 3.0d0, 0.0d0], [2, 3])
        integer(int32) :: kd_indices(3), dim_order(2), ierr
        real(real64)   :: labels(3)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 3_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_density_labels_invalid_bandwidth_percentile: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call density_labels_alloc(vectors, 2_int32, 3_int32, kd_indices, dim_order, k_density=2_int32, &
                                  bandwidth_percentile=101.0d0, labels=labels, ierr=ierr)
        call assert_true(is_err(ierr), "density_labels should reject bandwidth_percentile > 100")
    end subroutine test_density_labels_invalid_bandwidth_percentile

    !> The center of an evenly-spaced plus shape has all 4 of its k_density=4 neighbors at
    !| the identical distance 0.1 -- any percentile of a constant sample is that same
    !| constant, so the bandwidth here is simply 0.1, not a degenerate value needing a floor
    !| the way the earlier MAD-based formula did (MAD of an all-equal sample is exactly 0).
    !| Kept as a regression test even though the mechanism it originally guarded against
    !| (an underflow to a true 0.0) no longer applies to this formula at all: it still checks
    !| the basic sanity property that a symmetric neighborhood produces a genuine, strictly
    !| positive label.
    subroutine test_density_labels_symmetric_neighborhood_does_not_underflow()
        real(real64)   :: vectors(2, 5) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0], [2, 5])
        integer(int32) :: kd_indices(5), dim_order(2), ierr
        real(real64)   :: labels(5)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_density_labels_symmetric_neighborhood_does_not_underflow: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call density_labels_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, k_density=4_int32, &
                                  labels=labels, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(labels(1) > 0.0d0, "density_labels: a perfectly symmetric neighborhood is not exactly 0")
    end subroutine test_density_labels_symmetric_neighborhood_does_not_underflow

    !> On an evenly-spaced 11-point line, every interior point's k_density=4 nearest
    !| neighbors form the identical distance pattern [1,1,2,2] by translation symmetry, so
    !| all interior points (3..9) must get exactly the same density label -- an
    !| implementation-independent invariance check that does not rely on the exact formula.
    subroutine test_density_labels_uniform_interior_points_agree()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: labels(11)
        integer(int32) :: i

        call build_line_fixture(vectors, kd_indices, dim_order)

        call density_labels_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, k_density=4_int32, &
                                  labels=labels, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels failed unexpectedly: ', ierr
            error stop
        end if

        do i = 4, 9
            call assert_equal_real(labels(i), labels(3), 1.0d-9, "density_labels: interior points agree by symmetry")
        end do
    end subroutine test_density_labels_uniform_interior_points_agree

    !> A dense cluster (spacing 0.1) and a sparse cluster (spacing 2.0), far enough apart
    !| that k_density=2 never crosses between them: the dense cluster's adaptive bandwidth is
    !| far smaller, so its members must get a strictly higher density label.
    subroutine test_density_labels_dense_vs_sparse()
        real(real64)   :: vectors(2, 6) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.2d0, 0.0d0, &
                          100.0d0, 0.0d0, 102.0d0, 0.0d0, 104.0d0, 0.0d0], [2, 6])
        integer(int32) :: kd_indices(6), dim_order(2), ierr
        real(real64)   :: labels(6)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 6_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_density_labels_dense_vs_sparse: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call density_labels_alloc(vectors, 2_int32, 6_int32, kd_indices, dim_order, k_density=2_int32, &
                                  labels=labels, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(labels(2) > labels(5), "density_labels: the dense cluster's label exceeds the sparse cluster's")
    end subroutine test_density_labels_dense_vs_sparse

    !> A kd_indices entry out of [1,n_vectors] must be rejected by validation.
    subroutine test_density_labels_invalid_kd_indices()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: labels(11)

        call build_line_fixture(vectors, kd_indices, dim_order)
        kd_indices(1) = 12

        call density_labels_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, k_density=4_int32, &
                                  labels=labels, ierr=ierr)
        call assert_true(is_err(ierr), "density_labels should reject a kd_indices entry > n_vectors")
    end subroutine test_density_labels_invalid_kd_indices

    subroutine test_density_labels_k_density_too_large()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: labels(11)

        call build_line_fixture(vectors, kd_indices, dim_order)

        call density_labels_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, k_density=11_int32, &
                                  labels=labels, ierr=ierr)
        call assert_true(is_err(ierr), "density_labels should reject k_density > n_vectors - 1")
    end subroutine test_density_labels_k_density_too_large

    !> Regression test for a genuine crash: on branch smoothing-onward, calling
    !| `density_labels_alloc`/`seeds_alloc` with `k_density` *omitted* on a dataset smaller
    !| than the default (30) corrupted memory (an out-of-bounds k-NN query for 31 neighbors
    !| among 5 points) and crashed later, in unrelated code -- see
    !| `misc/code_gen_footgun.md`'s third entry for the generator-level root cause (an omitted
    !| optional's default is never validated against a runtime-dependent DM_MAX the way an
    !| explicit value is) and `density_labels_kernel`'s own `min(actual_k_density,
    !| n_vectors - 1)` clamp for the fix. This asserts the fix, not just its absence of a
    !| crash: omitting `k_density` here must resolve to exactly `n_vectors - 1` -- the same
    !| result an explicit `k_density = n_vectors - 1` call already produces.
    subroutine test_density_labels_omitted_k_density_is_clamped()
        real(real64)   :: vectors(2, 5) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0], [2, 5])
        integer(int32) :: kd_indices(5), dim_order(2), ierr
        real(real64)   :: labels_omitted(5), labels_explicit(5)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_density_labels_omitted_k_density_is_clamped: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call density_labels_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, labels=labels_omitted, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels_alloc (k_density omitted) failed unexpectedly: ', ierr
            error stop
        end if

        call density_labels_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, k_density=4_int32, &
                                  labels=labels_explicit, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'density_labels_alloc (k_density=4) failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_array_real(labels_omitted, labels_explicit, 5_int32, 1.0d-12, &
                                     "density_labels: an omitted k_density on N=5 clamps to exactly k_density=4")
    end subroutine test_density_labels_omitted_k_density_is_clamped

    ! --- seeds -------------------------------------------------------------------

    !> Two separated 2-point clusters, k_density=1: each point's single nearest neighbor is
    !| always its own cluster-mate (0.1 apart), never the other cluster (10 apart), so both
    !| the density ranking and the coverage radius (calc_ensemble_growth_radius on
    !| k_density=1, i.e. simply "distance to that one neighbor") stay entirely local -- and
    !| covering that one neighbor's exact distance is, by construction, enough to cover the
    !| whole (2-point) cluster from a single pick. Deliberately k_density=1 and 2-point
    !| clusters, not the larger, more "natural-looking" symmetric clusters an earlier version
    !| of this test used: with k_density>1, a cluster member's own median-of-k-neighbors
    !| coverage radius is generally *smaller* than the cluster's full diameter (a median
    !| undershoots a max), so a single seed does not reliably cover a larger cluster's
    !| farthest member -- see `misc/STC-experiments/README.md` for where this was first
    !| noticed on real data, and the discussion in the mod_STC.md seeding section for why
    !| that is an expected property of median-based coverage, not a bug to test around here.
    subroutine test_seeds_two_separated_clusters()
        real(real64)   :: vectors(2, 4) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 10.0d0, 0.0d0, 10.1d0, 0.0d0], [2, 4])
        integer(int32) :: kd_indices(4), dim_order(2), ierr
        logical        :: is_seed_mask(4)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 4_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_two_separated_clusters: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 4_int32, kd_indices, dim_order, k_density=1_int32, &
                         is_seed_mask=is_seed_mask, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(count(is_seed_mask, kind=int32), 2_int32, "seeds: two separated clusters give 2 seeds")
        call assert_true(any(is_seed_mask(1:2)), "seeds: cluster A (indices 1-2) has a seed")
        call assert_true(any(is_seed_mask(3:4)), "seeds: cluster B (indices 3-4) has a seed")
    end subroutine test_seeds_two_separated_clusters

    !> A single tight, 2-point cluster with k_density = n_vectors - 1 = 1: whichever point is
    !| picked first, its coverage radius (distance to its one neighbor) exactly covers that
    !| neighbor -- the whole cluster -- from a single pick. See
    !| test_seeds_two_separated_clusters above for why this stays a 2-point fixture.
    subroutine test_seeds_single_cluster_one_seed()
        real(real64)   :: vectors(2, 2) = reshape([0.0d0, 0.0d0, 0.1d0, 0.0d0], [2, 2])
        integer(int32) :: kd_indices(2), dim_order(2), ierr
        logical        :: is_seed_mask(2)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 2_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_single_cluster_one_seed: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 2_int32, kd_indices, dim_order, k_density=1_int32, &
                         is_seed_mask=is_seed_mask, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_int(count(is_seed_mask, kind=int32), 1_int32, "seeds: a single tight cluster gives 1 seed")
    end subroutine test_seeds_single_cluster_one_seed

    !> The shared 11-point line fixture, k_density=4. A wider exclusion radius (100th
    !| percentile of the k_density distances, i.e. the farthest neighbor -- 2.0, vs. the
    !| default 50th-percentile median of 1.5) suppresses more of the line per seed, so fewer
    !| seeds are needed to cover it. Both outcomes cross-checked against the actual Python
    !| binding's output on this exact fixture before being hardcoded here.
    subroutine test_seeds_exclusion_radius_percentile_widens_coverage()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        logical        :: mask_default(11), mask_wide(11), expected_default(11), expected_wide(11)

        call build_line_fixture(vectors, kd_indices, dim_order)

        call seeds_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, k_density=4_int32, &
                         is_seed_mask=mask_default, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds_alloc (default exclusion_radius_percentile) failed unexpectedly: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, k_density=4_int32, &
                         exclusion_radius_percentile=100.0d0, is_seed_mask=mask_wide, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds_alloc (exclusion_radius_percentile=100) failed unexpectedly: ', ierr
            error stop
        end if

        expected_default = .false.
        expected_default([2, 4, 6, 8, 10]) = .true.
        expected_wide = .false.
        expected_wide([1, 4, 8, 11]) = .true.

        call assert_equal_array_logical(mask_default, expected_default, 11_int32, &
                                        "seeds: default exclusion_radius_percentile (median) gives 5 seeds")
        call assert_equal_array_logical(mask_wide, expected_wide, 11_int32, &
                                        "seeds: exclusion_radius_percentile=100 (max) gives 4 seeds")
        call assert_true(count(mask_wide, kind=int32) < count(mask_default, kind=int32), &
                         "seeds: a wider exclusion radius needs fewer seeds to cover the same line")
    end subroutine test_seeds_exclusion_radius_percentile_widens_coverage

    subroutine test_seeds_invalid_exclusion_radius_percentile()
        real(real64)   :: vectors(2, 5) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0], [2, 5])
        integer(int32) :: kd_indices(5), dim_order(2), ierr
        logical        :: is_seed_mask(5)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_invalid_exclusion_radius_percentile: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, k_density=4_int32, &
                         exclusion_radius_percentile=101.0d0, is_seed_mask=is_seed_mask, ierr=ierr)
        call assert_true(is_err(ierr), "seeds should reject exclusion_radius_percentile > 100")
    end subroutine test_seeds_invalid_exclusion_radius_percentile

    subroutine test_seeds_invalid_k_density()
        real(real64)   :: vectors(2, 5) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0], [2, 5])
        integer(int32) :: kd_indices(5), dim_order(2), ierr
        logical        :: is_seed_mask(5)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_invalid_k_density: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, k_density=0_int32, &
                         is_seed_mask=is_seed_mask, ierr=ierr)
        call assert_true(is_err(ierr), "seeds should reject k_density < 1")
    end subroutine test_seeds_invalid_k_density

    !> Regression test for the same crash as
    !| test_density_labels_omitted_k_density_is_clamped above, exercised through `seeds`
    !| itself (which resolves `k_density` a second time, independently, for
    !| `calc_ensemble_growth_radius_kernel`'s own coverage-radius call): omitting `k_density`
    !| on N=5 must produce exactly the same seed selection as explicitly passing
    !| `k_density = n_vectors - 1 = 4`, not a crash.
    subroutine test_seeds_omitted_k_density_is_clamped()
        real(real64)   :: vectors(2, 5) = reshape([ &
                          0.0d0, 0.0d0, 0.1d0, 0.0d0, 0.0d0, 0.1d0, -0.1d0, 0.0d0, 0.0d0, -0.1d0], [2, 5])
        integer(int32) :: kd_indices(5), dim_order(2), ierr
        logical        :: mask_omitted(5), mask_explicit(5)

        dim_order = [1, 2]
        call build_kd_index_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_seeds_omitted_k_density_is_clamped: build_kd_index_alloc failed: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, is_seed_mask=mask_omitted, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds_alloc (k_density omitted) failed unexpectedly: ', ierr
            error stop
        end if

        call seeds_alloc(vectors, 2_int32, 5_int32, kd_indices, dim_order, k_density=4_int32, &
                         is_seed_mask=mask_explicit, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'seeds_alloc (k_density=4) failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_array_logical(mask_omitted, mask_explicit, 5_int32, &
                                        "seeds: an omitted k_density on N=5 clamps to exactly k_density=4")
    end subroutine test_seeds_omitted_k_density_is_clamped

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
