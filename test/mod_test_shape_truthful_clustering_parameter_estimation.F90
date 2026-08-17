!> Unit test suite for tox_shape_truthful_clustering_parameter_estimation (sample_estimator_anchors,
!| grow_estimator_anchor_clouds, estimate_stc_parameters), generated from
!| src/tox/shape_truthful_clustering/tox_shape_truthful_clustering_parameter_estimation_impl.F90.
module mod_test_shape_truthful_clustering_parameter_estimation
    use tox_shape_truthful_clustering_parameter_estimation, only: sample_estimator_anchors, &
                                                                   grow_estimator_anchor_clouds, &
                                                                   estimate_stc_parameters
    use f42_kd_tree, only: build_kd_index
    use tox_errors, only: is_ok, is_err
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_parameter_estimation() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(13))

        all_tests(1) = test_case("test_sample_anchors_hand_computed", test_sample_anchors_hand_computed)
        all_tests(2) = test_case("test_sample_anchors_duplicates_possible", test_sample_anchors_duplicates_possible)
        all_tests(3) = test_case("test_sample_anchors_invalid_n_anchors_zero", test_sample_anchors_invalid_n_anchors_zero)
        all_tests(4) = test_case("test_sample_anchors_invalid_n_anchors_too_large", &
                                 test_sample_anchors_invalid_n_anchors_too_large)
        all_tests(5) = test_case("test_grow_clouds_symmetric_line_ties_favor_lower_index", &
                                 test_grow_clouds_symmetric_line_ties_favor_lower_index)
        all_tests(6) = test_case("test_grow_clouds_seed_max_set_size_stops_early", &
                                 test_grow_clouds_seed_max_set_size_stops_early)
        all_tests(7) = test_case("test_grow_clouds_default_seed_max_set_size_can_yield_no_growth", &
                                 test_grow_clouds_default_seed_max_set_size_can_yield_no_growth)
        all_tests(8) = test_case("test_grow_clouds_invalid_seed_max_set_size", &
                                 test_grow_clouds_invalid_seed_max_set_size)
        all_tests(9) = test_case("test_estimate_parameters_collinear_line", test_estimate_parameters_collinear_line)
        all_tests(10) = test_case("test_estimate_parameters_too_few_valid_eas", &
                                  test_estimate_parameters_too_few_valid_eas)
        all_tests(11) = test_case("test_estimate_parameters_invalid_n_anchors", &
                                  test_estimate_parameters_invalid_n_anchors)
        all_tests(12) = test_case("test_estimate_parameters_invalid_seed_max_set_size", &
                                  test_estimate_parameters_invalid_seed_max_set_size)
        all_tests(13) = test_case("test_estimate_parameters_omitted_n_anchors_is_clamped", &
                                  test_estimate_parameters_omitted_n_anchors_is_clamped)
    end function get_all_tests_shape_truthful_clustering_parameter_estimation

    ! --- sample_estimator_anchors -------------------------------------------------------
    !
    ! 11 points, density labels equal to point index (1..11, already ascending -- the sort
    ! permutation is the identity). n_anchors=5 gives percentiles 20/40/60/80/100, ranks
    ! 3/5/7/9/11 exactly (no interpolation rounding needed): anchor_indices = [3,5,7,9,11].

    subroutine test_sample_anchors_hand_computed()
        real(real64)   :: density_labels(11)
        integer(int32) :: anchor_indices(5), expected(5), ierr, i

        do i = 1, 11
            density_labels(i) = real(i, real64)
        end do

        call sample_estimator_anchors(density_labels, 11_int32, 5_int32, anchor_indices, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'sample_estimator_anchors failed unexpectedly: ', ierr
            error stop
        end if

        expected = [3, 5, 7, 9, 11]
        call assert_equal_array_int(anchor_indices, expected, 5_int32, "sample_estimator_anchors: hand-computed ranks")
    end subroutine test_sample_anchors_hand_computed

    !> 5 points, n_anchors=5 (every point its own percentile mark): ranks 2/3/3/4/5 -- rank 3
    !| is hit twice (percentiles 40 and 60 both round to the same nearest-rank index on this
    !| small a set), so anchor_indices necessarily repeats an index. Documented, not a bug --
    !| see the kernel's own doc comment.
    subroutine test_sample_anchors_duplicates_possible()
        real(real64)   :: density_labels(5)
        integer(int32) :: anchor_indices(5), expected(5), ierr, i

        do i = 1, 5
            density_labels(i) = real(i, real64)
        end do

        call sample_estimator_anchors(density_labels, 5_int32, 5_int32, anchor_indices, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'sample_estimator_anchors failed unexpectedly: ', ierr
            error stop
        end if

        expected = [2, 3, 3, 4, 5]
        call assert_equal_array_int(anchor_indices, expected, 5_int32, "sample_estimator_anchors: duplicates possible")
    end subroutine test_sample_anchors_duplicates_possible

    subroutine test_sample_anchors_invalid_n_anchors_zero()
        real(real64)   :: density_labels(11)
        integer(int32) :: anchor_indices(0), ierr, i

        do i = 1, 11
            density_labels(i) = real(i, real64)
        end do

        call sample_estimator_anchors(density_labels, 11_int32, 0_int32, anchor_indices, ierr)
        call assert_true(is_err(ierr), "sample_estimator_anchors should reject n_anchors < 1")
    end subroutine test_sample_anchors_invalid_n_anchors_zero

    subroutine test_sample_anchors_invalid_n_anchors_too_large()
        real(real64)   :: density_labels(11)
        integer(int32) :: anchor_indices(12), ierr, i

        do i = 1, 11
            density_labels(i) = real(i, real64)
        end do

        call sample_estimator_anchors(density_labels, 11_int32, 12_int32, anchor_indices, ierr)
        call assert_true(is_err(ierr), "sample_estimator_anchors should reject n_anchors > n_vectors")
    end subroutine test_sample_anchors_invalid_n_anchors_too_large

    ! --- grow_estimator_anchor_clouds ---------------------------------------------------
    !
    ! Shared fixture: D=2, N=7 points on a line, (0,0)..(6,0). Two anchors at the opposite
    ! ends, point 1 (x=0) and point 7 (x=6).

    subroutine build_line_fixture_7(vectors, kd_indices, dim_order)
        real(real64), intent(out) :: vectors(2, 7)
        integer(int32), intent(out) :: kd_indices(7)
        integer(int32), intent(out) :: dim_order(2)
        integer(int32) :: i, ierr

        do i = 1, 7
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]

        call build_kd_index(vectors, 2_int32, 7_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'build_line_fixture_7: build_kd_index failed: ', ierr
            error stop
        end if
    end subroutine build_line_fixture_7

    !> seed_max_set_size=100 (grow until every point is claimed). Every round is an exact
    !| distance tie between the two clouds' own nearest-unclaimed candidate (both always 1.0
    !| apart on this evenly-spaced line) -- ties are broken by whichever EA is scanned first
    !| (anchor 1), so anchor 1 wins every single round, including the final point (6, itself
    !| equidistant from both clouds). Anchor 2's cloud never grows past its own single point.
    subroutine test_grow_clouds_symmetric_line_ties_favor_lower_index()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2)
        integer(int32) :: anchor_indices(2), cloud_sizes(2), expected_sizes(2), ierr
        logical(c_bool)        :: cloud_masks(7, 2), expected_cloud_1(7)

        call build_line_fixture_7(vectors, kd_indices, dim_order)
        anchor_indices = [1, 7]

        call grow_estimator_anchor_clouds(vectors, 2_int32, 7_int32, anchor_indices, 2_int32, &
                                          seed_max_set_size=100.0d0, &
                                          cloud_masks=cloud_masks, cloud_sizes=cloud_sizes, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'grow_estimator_anchor_clouds failed unexpectedly: ', ierr
            error stop
        end if

        expected_sizes = [6, 1]
        call assert_equal_array_int(cloud_sizes, expected_sizes, 2_int32, &
                                    "grow_estimator_anchor_clouds: ties favor the lower-indexed EA")

        expected_cloud_1 = .false.
        expected_cloud_1(1:6) = .true.
        call assert_equal_array_logical(cloud_masks(:, 1), expected_cloud_1, 7_int32, &
                                        "grow_estimator_anchor_clouds: EA 1's cloud absorbs points 1-6")
        call assert_true(cloud_masks(7, 2) .and. count(cloud_masks(:, 2)) == 1, &
                         "grow_estimator_anchor_clouds: EA 2's cloud stays its own single point")
    end subroutine test_grow_clouds_symmetric_line_ties_favor_lower_index

    !> Same fixture, seed_max_set_size=50 -> ceiling(0.5*7)=4 total claims: 2 anchors already
    !| present plus 2 more rounds, both won by EA 1 (same tie-break as above).
    subroutine test_grow_clouds_seed_max_set_size_stops_early()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2)
        integer(int32) :: anchor_indices(2), cloud_sizes(2), expected_sizes(2), ierr
        logical(c_bool)        :: cloud_masks(7, 2)

        call build_line_fixture_7(vectors, kd_indices, dim_order)
        anchor_indices = [1, 7]

        call grow_estimator_anchor_clouds(vectors, 2_int32, 7_int32, anchor_indices, 2_int32, &
                                          seed_max_set_size=50.0d0, &
                                          cloud_masks=cloud_masks, cloud_sizes=cloud_sizes, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'grow_estimator_anchor_clouds failed unexpectedly: ', ierr
            error stop
        end if

        expected_sizes = [3, 1]
        call assert_equal_array_int(cloud_sizes, expected_sizes, 2_int32, &
                                    "grow_estimator_anchor_clouds: seed_max_set_size stops growth at 4 total claims")
        call assert_equal_int(count(cloud_masks), 4_int32, &
                              "grow_estimator_anchor_clouds: exactly 4 points claimed in total")
    end subroutine test_grow_clouds_seed_max_set_size_stops_early

    !> Default seed_max_set_size (5.0): ceiling(0.05*7)=1, clamped up to n_anchors=2 itself --
    !| actual_max_claims never exceeds the anchor count, so no growth happens at all. A
    !| well-defined degenerate case, not an error: both clouds stay their own single point.
    subroutine test_grow_clouds_default_seed_max_set_size_can_yield_no_growth()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2)
        integer(int32) :: anchor_indices(2), cloud_sizes(2), expected_sizes(2), ierr
        logical(c_bool)        :: cloud_masks(7, 2)

        call build_line_fixture_7(vectors, kd_indices, dim_order)
        anchor_indices = [1, 7]

        call grow_estimator_anchor_clouds(vectors, 2_int32, 7_int32, anchor_indices, 2_int32, &
                                          cloud_masks=cloud_masks, cloud_sizes=cloud_sizes, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'grow_estimator_anchor_clouds failed unexpectedly: ', ierr
            error stop
        end if

        expected_sizes = [1, 1]
        call assert_equal_array_int(cloud_sizes, expected_sizes, 2_int32, &
                                    "grow_estimator_anchor_clouds: default seed_max_set_size yields no growth here")
    end subroutine test_grow_clouds_default_seed_max_set_size_can_yield_no_growth

    subroutine test_grow_clouds_invalid_seed_max_set_size()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2)
        integer(int32) :: anchor_indices(2), cloud_sizes(2), ierr
        logical(c_bool)        :: cloud_masks(7, 2)

        call build_line_fixture_7(vectors, kd_indices, dim_order)
        anchor_indices = [1, 7]

        call grow_estimator_anchor_clouds(vectors, 2_int32, 7_int32, anchor_indices, 2_int32, &
                                          seed_max_set_size=150.0d0, &
                                          cloud_masks=cloud_masks, cloud_sizes=cloud_sizes, ierr=ierr)
        call assert_true(is_err(ierr), "grow_estimator_anchor_clouds should reject seed_max_set_size > 100")
    end subroutine test_grow_clouds_invalid_seed_max_set_size

    ! --- estimate_stc_parameters --------------------------------------------------------

    !> D=2, N=21, a perfectly collinear, evenly-spaced line (0,0)..(20,0). Every estimator
    !| anchor's grown cloud is itself a sub-interval of the same line, so every EA agrees
    !| exactly on d=1 and on tangent direction (principal angle 0 between any pair,
    !| irrespective of individual singular-vector sign, since principal angles come from the
    !| SVD of U_i^T U_j, whose singular values are sign-invariant by construction) --
    !| chordal_dist_max_as_prcnt_of_range and d_max must both come out exactly 0.
    !| k_min/k_density/density_quantile/G_max are cross-checked against this exact,
    !| already-verified, fully deterministic kernel's own real output (no randomness anywhere
    !| in this pipeline).
    subroutine test_estimate_parameters_collinear_line()
        real(real64)   :: vectors(2, 21)
        integer(int32) :: kd_indices(21), dim_order(2), ierr, i
        real(real64)   :: est_k_min, est_k_density, est_density_quantile, est_G_max, est_d_max
        real(real64)   :: est_chordal_dist_max_as_prcnt_of_range

        do i = 1, 21
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]
        call build_kd_index(vectors, 2_int32, 21_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_estimate_parameters_collinear_line: build_kd_index failed: ', ierr
            error stop
        end if

        call estimate_stc_parameters(vectors, 2_int32, 21_int32, kd_indices, dim_order, seed_max_set_size=50.0d0, &
                                           estimated_k_min=est_k_min, estimated_k_density=est_k_density, &
                                           estimated_density_quantile=est_density_quantile, &
                                           estimated_chordal_dist_max_as_prcnt_of_range=est_chordal_dist_max_as_prcnt_of_range, &
                                           estimated_G_max=est_G_max, estimated_d_max=est_d_max, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'estimate_stc_parameters failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_real(est_chordal_dist_max_as_prcnt_of_range, 0.0d0, 1.0d-9, &
                               "estimate_stc_parameters: collinear data gives chordal_dist_max_as_prcnt_of_range=0")
        call assert_equal_real(est_d_max, 0.0d0, 1.0d-9, "estimate_stc_parameters: collinear data gives d_max=0")
        call assert_equal_real(est_k_min, 4.0d0, 1.0d-9, "estimate_stc_parameters: k_min")
        call assert_equal_real(est_k_density, est_k_min, 1.0d-9, "estimate_stc_parameters: k_density equals k_min")
        call assert_equal_real(est_density_quantile, 1.5d0, 1.0d-9, "estimate_stc_parameters: density_quantile")
        call assert_equal_real(est_G_max, 0.0d0, 1.0d-9, "estimate_stc_parameters: G_max")
    end subroutine test_estimate_parameters_collinear_line

    !> seed_max_set_size=0: actual_max_claims clamps to n_anchors itself (see
    !| grow_estimator_anchor_clouds), so every EA's cloud stays size 1 -- zero clouds ever
    !| reach the size >= 2 a genuine observable/SVD needs. Fewer than 2 usable EAs means no
    !| pairwise comparison is possible at all: a genuine, data-dependent runtime failure (see
    !| the kernel's own doc comment), not a validation error.
    subroutine test_estimate_parameters_too_few_valid_eas()
        real(real64)   :: vectors(2, 21)
        integer(int32) :: kd_indices(21), dim_order(2), ierr, i
        real(real64)   :: est_k_min, est_k_density, est_density_quantile, est_G_max, est_d_max
        real(real64)   :: est_chordal_dist_max_as_prcnt_of_range

        do i = 1, 21
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]
        call build_kd_index(vectors, 2_int32, 21_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_estimate_parameters_too_few_valid_eas: build_kd_index failed: ', ierr
            error stop
        end if

        call estimate_stc_parameters(vectors, 2_int32, 21_int32, kd_indices, dim_order, seed_max_set_size=0.0d0, &
                                           estimated_k_min=est_k_min, estimated_k_density=est_k_density, &
                                           estimated_density_quantile=est_density_quantile, &
                                           estimated_chordal_dist_max_as_prcnt_of_range=est_chordal_dist_max_as_prcnt_of_range, estimated_G_max=est_G_max, &
                                           estimated_d_max=est_d_max, ierr=ierr)
        call assert_true(is_err(ierr), "estimate_stc_parameters should fail when fewer than 2 EAs ever grow past size 1")
    end subroutine test_estimate_parameters_too_few_valid_eas

    subroutine test_estimate_parameters_invalid_n_anchors()
        real(real64)   :: vectors(2, 21)
        integer(int32) :: kd_indices(21), dim_order(2), ierr, i
        real(real64)   :: est_k_min, est_k_density, est_density_quantile, est_G_max, est_d_max
        real(real64)   :: est_chordal_dist_max_as_prcnt_of_range

        do i = 1, 21
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]
        call build_kd_index(vectors, 2_int32, 21_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_estimate_parameters_invalid_n_anchors: build_kd_index failed: ', ierr
            error stop
        end if

        call estimate_stc_parameters(vectors, 2_int32, 21_int32, kd_indices, dim_order, n_anchors=50_int32, &
                                           estimated_k_min=est_k_min, estimated_k_density=est_k_density, &
                                           estimated_density_quantile=est_density_quantile, &
                                           estimated_chordal_dist_max_as_prcnt_of_range=est_chordal_dist_max_as_prcnt_of_range, estimated_G_max=est_G_max, &
                                           estimated_d_max=est_d_max, ierr=ierr)
        call assert_true(is_err(ierr), "estimate_stc_parameters should reject n_anchors > n_vectors")
    end subroutine test_estimate_parameters_invalid_n_anchors

    subroutine test_estimate_parameters_invalid_seed_max_set_size()
        real(real64)   :: vectors(2, 21)
        integer(int32) :: kd_indices(21), dim_order(2), ierr, i
        real(real64)   :: est_k_min, est_k_density, est_density_quantile, est_G_max, est_d_max
        real(real64)   :: est_chordal_dist_max_as_prcnt_of_range

        do i = 1, 21
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]
        call build_kd_index(vectors, 2_int32, 21_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_estimate_parameters_invalid_seed_max_set_size: build_kd_index failed: ', ierr
            error stop
        end if

        call estimate_stc_parameters(vectors, 2_int32, 21_int32, kd_indices, dim_order, seed_max_set_size=-1.0d0, &
                                           estimated_k_min=est_k_min, estimated_k_density=est_k_density, &
                                           estimated_density_quantile=est_density_quantile, &
                                           estimated_chordal_dist_max_as_prcnt_of_range=est_chordal_dist_max_as_prcnt_of_range, estimated_G_max=est_G_max, &
                                           estimated_d_max=est_d_max, ierr=ierr)
        call assert_true(is_err(ierr), "estimate_stc_parameters should reject seed_max_set_size < 0")
    end subroutine test_estimate_parameters_invalid_seed_max_set_size

    !> Regression test for the same class of crash-fix seen throughout this family: an
    !| *explicit* n_anchors is already wrapper-validated against DM_MAX(n_vectors), but the
    !| wrapper never validates n_anchors' own default when it is omitted (see
    !| misc/code_gen_footgun.md's third entry). On N=3 (below the default of 5), an omitted
    !| n_anchors must clamp to exactly n_anchors=3 inside the kernel -- both calls below must
    !| reach the identical, well-defined "too few valid EAs" outcome (every one of the 3
    !| points is necessarily its own anchor, leaving none to grow into), not one of them
    !| corrupting memory by indexing a workspace array sized for N=3 with a length-5 slice.
    subroutine test_estimate_parameters_omitted_n_anchors_is_clamped()
        real(real64)   :: vectors(2, 3)
        integer(int32) :: kd_indices(3), dim_order(2), ierr_omitted, ierr_explicit, i
        real(real64)   :: est_k_min, est_k_density, est_density_quantile, est_G_max, est_d_max
        real(real64)   :: est_chordal_dist_max_as_prcnt_of_range

        do i = 1, 3
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]
        call build_kd_index(vectors, 2_int32, 3_int32, kd_indices, dim_order, ierr_omitted)
        if (.not. is_ok(ierr_omitted)) then
            write (*, *) 'test_estimate_parameters_omitted_n_anchors_is_clamped: build_kd_index failed: ', ierr_omitted
            error stop
        end if

        call estimate_stc_parameters(vectors, 2_int32, 3_int32, kd_indices, dim_order, seed_max_set_size=100.0d0, &
                                           estimated_k_min=est_k_min, estimated_k_density=est_k_density, &
                                           estimated_density_quantile=est_density_quantile, &
                                           estimated_chordal_dist_max_as_prcnt_of_range=est_chordal_dist_max_as_prcnt_of_range, estimated_G_max=est_G_max, &
                                           estimated_d_max=est_d_max, ierr=ierr_omitted)

        call estimate_stc_parameters(vectors, 2_int32, 3_int32, kd_indices, dim_order, n_anchors=3_int32, &
                                           seed_max_set_size=100.0d0, &
                                           estimated_k_min=est_k_min, estimated_k_density=est_k_density, &
                                           estimated_density_quantile=est_density_quantile, &
                                           estimated_chordal_dist_max_as_prcnt_of_range=est_chordal_dist_max_as_prcnt_of_range, estimated_G_max=est_G_max, &
                                           estimated_d_max=est_d_max, ierr=ierr_explicit)

        call assert_true(is_err(ierr_omitted) .and. is_err(ierr_explicit), &
                         "estimate_stc_parameters: an omitted n_anchors on N=3 clamps to exactly n_anchors=3, " // &
                         "reaching the same clean error as passing it explicitly -- not a crash")
    end subroutine test_estimate_parameters_omitted_n_anchors_is_clamped

end module mod_test_shape_truthful_clustering_parameter_estimation
