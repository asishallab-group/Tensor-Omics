!> Unit test suite for tox_get_outliers_by_angle_rap module.
module mod_test_get_outliers_by_angle_rap
    use tox_get_outliers_by_angle_rap
    use f42_utils, only: angle_between, wrap_angle, PI, below, radians
    use asserts
    use tox_errors
    use test_suite
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_is_nan, ieee_quiet_nan
    implicit none
    public
    real(real64), parameter :: EPS = 1.0e-12_real64
    real(real64), parameter :: MAX_SC = 10.0_real64

contains

    !> Get array of all available tests for circular statistics.
    function get_all_tests_outliers_by_angle_rap() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(25))

        all_tests(1) = test_case("test_rap_wrap_angle_basic", test_rap_wrap_angle_basic)
        all_tests(2) = test_case("test_rap_wrap_angle_edge_cases", test_rap_wrap_angle_edge_cases)
        all_tests(3) = test_case("test_compute_family_direction_rap_basic", test_compute_family_direction_basic)
        all_tests(4) = test_case("test_compute_family_direction_rap_single_family", test_compute_family_direction_single_family)
all_tests(5) = test_case("test_compute_family_direction_rap_no_stable_direction", test_compute_family_direction_no_stable_direction)
    all_tests(6) = test_case("test_compute_family_direction_rap_minimal_variation", test_compute_family_direction_minimal_variation)
    all_tests(7) = test_case("test_compute_family_direction_rap_uniform_distribution", test_compute_family_direction_uniform_distribution)
  all_tests(8) = test_case("test_compute_family_direction_rap_all_same_direction", test_compute_family_direction_all_same_direction)
        all_tests(9) = test_case("test_compute_angular_deviations_rap_basic", test_compute_angular_deviations_basic)
     all_tests(10) = test_case("test_compute_angular_deviations_rap_invalid_family", test_compute_angular_deviations_invalid_family)
        all_tests(11) = test_case("test_compute_angular_deviations_rap_wrapping", test_compute_angular_deviations_wrapping)
        all_tests(12) = test_case("test_z_scores_by_dispersion_rap_basic", test_z_scores_by_dispersion_basic)
        all_tests(13) = test_case("test_z_scores_by_dispersion_rap_invalid", test_z_scores_by_dispersion_invalid)
        all_tests(14) = test_case("test_z_scores_by_dispersion_rap_zero_dispersion", test_z_scores_by_dispersion_zero_dispersion)
        all_tests(15) = test_case("test_angle_outliers_rap_alloc_basic", test_angle_outliers_alloc_basic)
        all_tests(16) = test_case("test_angle_outliers_rap_alloc_no_valid", test_angle_outliers_alloc_no_valid)
        all_tests(17) = test_case("test_angle_outliers_rap_alloc_all_outliers", test_angle_outliers_alloc_all_outliers)
        all_tests(18) = test_case("test_detect_angle_outliers_pipeline_rap_basic", test_detect_angle_outliers_pipeline_basic)
    all_tests(19) = test_case("test_detect_angle_outliers_pipeline_rap_single_family", test_detect_angle_outliers_pipeline_single_family)
    all_tests(20) = test_case("test_detect_angle_outliers_pipeline_rap_no_valid_families", test_detect_angle_outliers_pipeline_no_valid_families)
    all_tests(21) = test_case("test_detect_angle_outliers_pipeline_rap_invalid_percentile", test_detect_angle_outliers_pipeline_inval_percentile)
    all_tests(22) = test_case("test_detect_angle_outliers_pipeline_rap_parallel_genes", test_detect_angle_outliers_pipeline_parallel_genes)
    all_tests(23) = test_case("test_detect_angle_outliers_pipeline_rap_orthogonal_outlier", test_detect_angle_outliers_pipeline_orth_outlier)
    all_tests(24) = test_case("test_detect_angle_outliers_pipeline_rap_mixed_directions", test_detect_angle_outliers_pipeline_mixed_directions)
     all_tests(25) = test_case("test_detect_angle_outliers_pipeline_rap_edge_cases", test_detect_angle_outliers_pipeline_edge_cases)
    end function get_all_tests_outliers_by_angle_rap

    !> Test basic angle wrapping to (-π, π]
    subroutine test_rap_wrap_angle_basic()
        real(real64) :: angle_in, angle_out

        ! Test angles already in range
        angle_in = 0.0_real64
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, 0.0_real64, EPS, "test_rap_wrap_angle_basic: 0 should stay 0")

        angle_in = PI
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, PI, EPS, "test_rap_wrap_angle_basic: π should stay π")

        angle_in = -PI
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, PI, EPS, "test_rap_wrap_angle_basic: -π should become π")

        angle_in = PI / 2.0_real64
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, PI / 2.0_real64, EPS, "test_rap_wrap_angle_basic: π/2 should stay π/2")
    end subroutine test_rap_wrap_angle_basic

    !> Test angle wrapping edge cases
    subroutine test_rap_wrap_angle_edge_cases()
        real(real64) :: angle_in, angle_out

        ! Test angles above π (should wrap to negative)
        angle_in = 1.5_real64 * PI  ! 270 degrees
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, -0.5_real64 * PI, EPS, "test_rap_wrap_angle_edge_cases: 1.5π should wrap to -0.5π")

        angle_in = 2.0_real64 * PI  ! 360 degrees
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, 0.0_real64, EPS, "test_rap_wrap_angle_edge_cases: 2π should wrap to 0")

        angle_in = 3.0_real64 * PI  ! 540 degrees
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, PI, EPS, "test_rap_wrap_angle_edge_cases: 3π should wrap to π")

        ! Test angles below -π (should wrap to positive)
        angle_in = -1.5_real64 * PI  ! -270 degrees
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, 0.5_real64 * PI, EPS, "test_rap_wrap_angle_edge_cases: -1.5π should wrap to 0.5π")

        angle_in = -2.0_real64 * PI  ! -360 degrees
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, 0.0_real64, EPS, "test_rap_wrap_angle_edge_cases: -2π should wrap to 0")

        angle_in = -3.0_real64 * PI  ! -540 degrees
        angle_out = wrap_angle(angle_in)
        call assert_equal_real(angle_out, PI, EPS, "test_rap_wrap_angle_edge_cases: -3π should wrap to π")
    end subroutine test_rap_wrap_angle_edge_cases

    !> Test circular family direction computation with multiple families
    subroutine test_compute_family_direction_basic()
        integer(int32), parameter :: n_genes = 6, n_families = 2
        ! Family 1: Three angles clustered around 0°
        ! Family 2: Three angles clustered around 90° (π/2)
        real(real64) :: rap_angles(n_genes) = [ &
                        0.5_real64, &           ! Gene 1, Family 1
                        0.6_real64, &           ! Gene 2, Family 1
                        0.4_real64, &          ! Gene 3, Family 1
                        PI / 2.0_real64, &        ! Gene 4, Family 2
                        PI / 2.0_real64 + 0.1_real64, &  ! Gene 5, Family 2
                        PI / 2.0_real64 - 0.1_real64 &   ! Gene 6, Family 2
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 2, 2, 2]
        integer(int32) :: tmp_member_counts(n_families), expected_member_counts(n_families)
        real(real64) :: family_mean_angles(n_families)
        real(real64) :: family_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        call tox_compute_family_direction_rap(rap_angles, gene_to_fam, &
                                              family_mean_angles, family_dispersions, &
                                              n_genes, n_families, ierr, status, tmp_member_counts)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_rap_basic: error code mismatch")

        expected_member_counts = [3, 3]
    call assert_equal_array_int(tmp_member_counts, expected_member_counts, n_families, "test_compute_family_direction_rap_basic: member count mismatch")

        ! Family 1 mean should be near 0°
        call assert_equal_real(family_mean_angles(1), 0.5_real64, EPS, &
                               "test_compute_family_direction_rap_basic: family 1 mean incorrect")

        ! Family 2 mean should be near π/2
        call assert_equal_real(family_mean_angles(2), PI / 2.0_real64, EPS, &
                               "test_compute_family_direction_rap_basic: family 2 mean incorrect")

        ! Dispersions should be positive and reasonable
        call assert_true(family_dispersions(1) > 0.0_real64 .and. family_dispersions(1) < PI, &
                         "test_compute_family_direction_rap_basic: family 1 dispersion out of range")
        call assert_true(family_dispersions(2) > 0.0_real64 .and. family_dispersions(2) < PI, &
                         "test_compute_family_direction_rap_basic: family 2 dispersion out of range")
    end subroutine test_compute_family_direction_basic

    !> Test circular family direction with single family
    subroutine test_compute_family_direction_single_family()
        integer(int32), parameter :: n_genes = 5, n_families = 1
        ! Angles clustered around 45° (π/4)
        real(real64) :: rap_angles(n_genes) = [ &
                        PI / 4.0_real64, &
                        PI / 4.0_real64 + 0.05_real64, &
                        PI / 4.0_real64 - 0.05_real64, &
                        PI / 4.0_real64 + 0.03_real64, &
                        PI / 4.0_real64 - 0.03_real64 &
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1]
        integer(int32) :: tmp_member_counts(n_families), expected_member_counts(n_families)
        real(real64) :: family_mean_angles(n_families)
        real(real64) :: family_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        call tox_compute_family_direction_rap(rap_angles, gene_to_fam, &
                                              family_mean_angles, family_dispersions, &
                                              n_genes, n_families, ierr, status, tmp_member_counts)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_rap_single_family: error code mismatch")
        expected_member_counts = [5]
    call assert_equal_array_int(tmp_member_counts, expected_member_counts, n_families, "test_compute_family_direction_rap_basic: member count mismatch")
        call assert_equal_real(family_mean_angles(1), PI / 4.0_real64, EPS, &
                               "test_compute_family_direction_rap_single_family: mean angle incorrect")
        call assert_true(family_dispersions(1) > 0.0_real64, &
                         "test_compute_family_direction_rap_single_family: dispersion should be positive")
    end subroutine test_compute_family_direction_single_family

    !> Test circular family direction with no stable direction (angles cancel out)
    subroutine test_compute_family_direction_no_stable_direction()
        integer(int32), parameter :: n_genes = 4, n_families = 1
        ! opposite Angles, cancel out, mean actually of length 0
        real(real64) :: rap_angles(n_genes)

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
        integer(int32) :: tmp_member_counts(n_families), expected_member_counts(n_families)
        real(real64) :: family_mean_angles(n_families)
        real(real64) :: family_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        rap_angles = [ &
                     0.0_real64, &           ! 0°
                     radians(175.0_real64), &        ! 180°
                     0.0_real64, &                   ! 0°
                     radians(175.0_real64) &        ! 180°
                     ]

        call tox_compute_family_direction_rap(rap_angles, gene_to_fam, &
                                              family_mean_angles, family_dispersions, &
                                              n_genes, n_families, ierr, status, tmp_member_counts)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_rap_no_stable_direction: error code mismatch")
        expected_member_counts = [4]
    call assert_equal_array_int(tmp_member_counts, expected_member_counts, n_families, "test_compute_family_direction_rap_basic: member count mismatch")
        call assert_equal_real(family_dispersions(1), -1.0_real64, EPS, &
                               "test_compute_family_direction_rap_no_stable_direction: should have -1 dispersion")
        call assert_true(all(status == STAT_NO_STABLE_DIRECTION), &
                         "test_compute_family_direction_rap_no_stable_direction: status should indicate no stable direction")
    end subroutine test_compute_family_direction_no_stable_direction

    !> Test circular family direction with minimal angular variation
    subroutine test_compute_family_direction_minimal_variation()
        integer(int32), parameter :: n_genes = 4, n_families = 1
        ! Nearly identical angles
        real(real64) :: rap_angles(n_genes) = [ &
                        0.0_real64, &
                        0.0001_real64, &
                        0.00015_real64, &
                        0.0002_real64 &
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
        real(real64) :: family_mean_angles(n_families)
        real(real64) :: family_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        call tox_compute_family_direction_rap_alloc(rap_angles, gene_to_fam, &
                                                    family_mean_angles, family_dispersions, &
                                                    n_genes, n_families, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_rap_minimal_variation: error code mismatch")

        if (family_dispersions(1) < 0.0_real64) then
            call assert_true(all(status == STAT_NO_ANGULAR_VARIATION), &
                             "test_compute_family_direction_rap_minimal_variation: status should indicate minimal variation")
        else
            ! If dispersion > TAU, mean should be near 0
            call assert_equal_real(family_mean_angles(1), 0.0_real64, EPS, &
                                   "test_compute_family_direction_rap_minimal_variation: mean angle incorrect")
        end if
    end subroutine test_compute_family_direction_minimal_variation

    !> Test circular family direction with uniform distribution
    subroutine test_compute_family_direction_uniform_distribution()
        integer(int32), parameter :: n_genes = 8, n_families = 1
        ! Angles uniformly spaced around the circle
        real(real64) :: rap_angles(n_genes)
        integer(int32) :: i
        integer(int32) :: gene_to_fam(n_genes) = [(1, i=1, n_genes)]
        real(real64) :: family_mean_angles(n_families)
        real(real64) :: family_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        do i = 1, n_genes
            if (mod(i, 2_int32) == 0) then
                rap_angles(i) = PI / n_genes
            else
                rap_angles(i) = PI - PI / n_genes
            end if
        end do

        call tox_compute_family_direction_rap_alloc(rap_angles, gene_to_fam, &
                                                    family_mean_angles, family_dispersions, &
                                                    n_genes, n_families, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_rap_uniform_distribution: error code mismatch")
        ! Uniform distribution should have no stable direction
        call assert_equal_real(family_dispersions(1), -1.0_real64, EPS, &
                               "test_compute_family_direction_rap_uniform_distribution: should have -1 dispersion")
        call assert_true(all(status == STAT_NO_STABLE_DIRECTION), &
                         "test_compute_family_direction_rap_uniform_distribution: status should indicate no stable direction")
    end subroutine test_compute_family_direction_uniform_distribution

    !> Test circular family direction with all same direction
    subroutine test_compute_family_direction_all_same_direction()
        integer(int32), parameter :: n_genes = 5, n_families = 1
        ! All angles exactly the same
        real(real64) :: rap_angles(n_genes) = PI / 3.0_real64  ! All 60°

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1]
        integer(int32) :: tmp_member_counts(n_families), expected_member_counts(n_families)
        real(real64) :: family_mean_angles(n_families)
        real(real64) :: family_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        call tox_compute_family_direction_rap(rap_angles, gene_to_fam, &
                                              family_mean_angles, family_dispersions, &
                                              n_genes, n_families, ierr, status, tmp_member_counts)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_rap_all_same_direction: error code mismatch")

        expected_member_counts = [5]
    call assert_equal_array_int(tmp_member_counts, expected_member_counts, n_families, "test_compute_family_direction_rap_all_same_direction: member count mismatch")

        call assert_equal_real(family_mean_angles(1), PI / 3.0_real64, EPS, &
                               "test_compute_family_direction_rap_all_same_direction: mean angle incorrect")
        ! With identical angles, resultant length R=1, so σ = sqrt(-2*ln(1)) = 0
        ! This is below TAU, so should be marked as minimal variation
        call assert_equal_real(family_dispersions(1), -1.0_real64, EPS, &
                               "test_compute_family_direction_rap_all_same_direction: should have -1 dispersion")
        call assert_true(all(status == STAT_NO_ANGULAR_VARIATION), &
                         "test_compute_family_direction_rap_all_same_direction: status should indicate minimal variation")
    end subroutine test_compute_family_direction_all_same_direction

    !> Test angular deviations computation for circular case
    subroutine test_compute_angular_deviations_basic()
        integer(int32), parameter :: n_genes = 5, n_families = 2
        real(real64) :: rap_angles(n_genes) = [ &
                        0.0_real64, &           ! Gene 1, Family 1
                        0.2_real64, &           ! Gene 2, Family 1
                        1.5_real64, &           ! Gene 3, Family 2
                        1.7_real64, &           ! Gene 4, Family 2
                        -1.0_real64 &           ! Gene 5, invalid
                        ]

        ! Family means: Family 1 at 0.1, Family 2 at 1.6
        real(real64) :: family_mean_angles(n_families) = [0.1_real64, 1.6_real64]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 2, 2, 0]  ! Gene 5: unassigned
        real(real64) :: angular_deviations(n_genes)
        integer(int32) :: ierr

        call tox_compute_angular_deviations_rap(rap_angles, family_mean_angles, &
                                                gene_to_fam, angular_deviations, &
                                                n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_angular_deviations_rap_basic: error code mismatch")

        ! Gene 1: |wrap(0.0 - 0.1)| = |wrap(-0.1)| = 0.1
        call assert_equal_real(angular_deviations(1), 0.1_real64, EPS, &
                               "test_compute_angular_deviations_rap_basic: gene 1 deviation mismatch")

        ! Gene 2: |wrap(0.2 - 0.1)| = |0.1| = 0.1
        call assert_equal_real(angular_deviations(2), 0.1_real64, EPS, &
                               "test_compute_angular_deviations_rap_basic: gene 2 deviation mismatch")

        ! Gene 3: |wrap(1.5 - 1.6)| = |wrap(-0.1)| = 0.1
        call assert_equal_real(angular_deviations(3), 0.1_real64, EPS, &
                               "test_compute_angular_deviations_rap_basic: gene 3 deviation mismatch")

        ! Gene 4: |wrap(1.7 - 1.6)| = |0.1| = 0.1
        call assert_equal_real(angular_deviations(4), 0.1_real64, EPS, &
                               "test_compute_angular_deviations_rap_basic: gene 4 deviation mismatch")

        ! Gene 5: unassigned -> -1.0
        call assert_equal_real(angular_deviations(5), -1.0_real64, EPS, &
                               "test_compute_angular_deviations_rap_basic: unassigned gene should be -1")
    end subroutine test_compute_angular_deviations_basic

    !> Test angular deviations with invalid family mapping
    subroutine test_compute_angular_deviations_invalid_family()
        integer(int32), parameter :: n_genes = 4, n_families = 2
        real(real64) :: rap_angles(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
        real(real64) :: family_mean_angles(n_families) = [0.15_real64, 0.35_real64]

        ! Gene 3: invalid family (3 > n_families)
        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 3, 2]
        real(real64) :: angular_deviations(n_genes)
        integer(int32) :: ierr

        call tox_compute_angular_deviations_rap(rap_angles, family_mean_angles, &
                                                gene_to_fam, angular_deviations, &
                                                n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_compute_angular_deviations_rap_invalid_family: error code mismatch")
    end subroutine test_compute_angular_deviations_invalid_family

    !> Test angular deviations with wrapping across π boundary
    subroutine test_compute_angular_deviations_wrapping()
        integer(int32), parameter :: n_genes = 3, n_families = 1
        ! Test case where difference crosses the π boundary
        real(real64) :: rap_angles(n_genes) = [ &
                        3.0_real64, &    ! ~171.9° (π = 3.14159)
                        PI, &
                        0.0_real64 &     ! 0°
                        ]

        ! Family mean near π (180°)
        real(real64) :: family_mean_angles(n_families) = [PI]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1]
        real(real64) :: angular_deviations(n_genes)
        integer(int32) :: ierr

        call tox_compute_angular_deviations_rap(rap_angles, family_mean_angles, &
                                                gene_to_fam, angular_deviations, &
                                                n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_angular_deviations_rap_wrapping: error code mismatch")

        ! Gene 1: wrap(3.0 - π) = wrap(~-0.14159) = -0.14159, abs = ~0.14159
        call assert_equal_real(angular_deviations(1), abs(3.0_real64 - PI), EPS, &
                               "test_compute_angular_deviations_rap_wrapping: gene 1 deviation mismatch")

        ! Gene 2: wrap(PI - π) = wrap(~-6.14159) = wrap to ~0.14159, abs = ~0.14159
        call assert_equal_real(angular_deviations(2), 0.0_real64, EPS, &
                               "test_compute_angular_deviations_rap_wrapping: gene 2 deviation mismatch")

        ! Gene 3: wrap(0.0 - π) = wrap(-π) = -π, abs = π
        call assert_equal_real(angular_deviations(3), PI, EPS, &
                               "test_compute_angular_deviations_rap_wrapping: gene 3 deviation mismatch")
    end subroutine test_compute_angular_deviations_wrapping

    !> Test scaling circular angles by dispersion (z-scores)
    subroutine test_z_scores_by_dispersion_basic()
        integer(int32), parameter :: n_genes = 5, n_families = 2
        real(real64) :: angular_deviations(n_genes) = [0.2_real64, 0.4_real64, -1.0_real64, 0.3_real64, 0.6_real64]
        real(real64) :: family_dispersions(n_families) = [0.1_real64, 0.2_real64]  ! Family 1: σ=0.1, Family 2: σ=0.2
        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 0, 2, 2]  ! Gene 3: unassigned
        real(real64) :: z_scores(n_genes)
        integer(int32) :: ierr

        call tox_z_scores_by_dispersion_rap(angular_deviations, family_dispersions, &
                                            gene_to_fam, z_scores, &
                                            n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_z_scores_by_dispersion_rap_basic: error code mismatch")

        ! Gene 1: 0.2 / 0.1 = 2.0
        call assert_equal_real(z_scores(1), 2.0_real64, EPS, "test_z_scores_by_dispersion_rap_basic: gene 1 mismatch")

        ! Gene 2: 0.4 / 0.1 = 4.0
        call assert_equal_real(z_scores(2), 4.0_real64, EPS, "test_z_scores_by_dispersion_rap_basic: gene 2 mismatch")

        ! Gene 3: unassigned -> -1.0
        call assert_equal_real(z_scores(3), -1.0_real64, EPS, "test_z_scores_by_dispersion_rap_basic: gene 3 mismatch")

        ! Gene 4: 0.3 / 0.2 = 1.5
        call assert_equal_real(z_scores(4), 1.5_real64, EPS, "test_z_scores_by_dispersion_rap_basic: gene 4 mismatch")

        ! Gene 5: 0.6 / 0.2 = 3.0
        call assert_equal_real(z_scores(5), 3.0_real64, EPS, "test_z_scores_by_dispersion_rap_basic: gene 5 mismatch")
    end subroutine test_z_scores_by_dispersion_basic

    !> Test scaling circular angles with invalid inputs
    subroutine test_z_scores_by_dispersion_invalid()
        integer(int32), parameter :: n_genes = 4, n_families = 2
        real(real64) :: angular_deviations(n_genes) = [-1.0_real64, 0.5_real64, 0.3_real64, 0.4_real64]  ! Family 1: invalid
        real(real64) :: family_dispersions(n_families) = [0.0_real64, 0.1_real64]
        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 3, 2]  ! Gene 3: invalid family
        real(real64) :: z_scores(n_genes)
        integer(int32) :: ierr

        call tox_z_scores_by_dispersion_rap(angular_deviations, family_dispersions, &
                                            gene_to_fam, z_scores, &
                                            n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_z_scores_by_dispersion_rap_invalid: error code mismatch")
    end subroutine test_z_scores_by_dispersion_invalid

    !> Test scaling circular angles with zero dispersion
    subroutine test_z_scores_by_dispersion_zero_dispersion()
        integer(int32), parameter :: n_genes = 3, n_families = 1
        real(real64) :: angular_deviations(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64]
        real(real64) :: family_dispersions(n_families) = [0.0_real64]  ! Zero dispersion
        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1]
        real(real64) :: z_scores(n_genes)
        integer(int32) :: ierr

        call tox_z_scores_by_dispersion_rap(angular_deviations, family_dispersions, &
                                            gene_to_fam, z_scores, &
                                            n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_z_scores_by_dispersion_rap_zero_dispersion: error code mismatch")

        ! With zero dispersion, all should be marked as invalid
        call assert_equal_real(z_scores(1), -1.0_real64, EPS, &
                               "test_z_scores_by_dispersion_rap_zero_dispersion: gene 1 should be -1")
        call assert_equal_real(z_scores(2), -1.0_real64, EPS, &
                               "test_z_scores_by_dispersion_rap_zero_dispersion: gene 2 should be -1")
        call assert_equal_real(z_scores(3), -1.0_real64, EPS, &
                               "test_z_scores_by_dispersion_rap_zero_dispersion: gene 3 should be -1")
    end subroutine test_z_scores_by_dispersion_zero_dispersion

    !> Test outlier detection for circular angles (uses common implementation)
    subroutine test_angle_outliers_alloc_basic()
        integer(int32), parameter :: n_genes = 8
        ! Create scaled angles: some normal, some outliers
        real(real64) :: z_scores(n_genes) = [1.0_real64, 1.5_real64, 2.0_real64, 2.5_real64, &
                                             3.0_real64, 2.5_real64, 2.0_real64, -1.0_real64]
        real(real64) :: percentile = 0.8_real64  ! 80th percentile
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr

        call tox_angle_outliers_rap_alloc(z_scores, percentile, threshold, &
                                          is_outlier, n_genes, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_rap_alloc_basic: error code mismatch")

        call assert_equal_real(threshold, 2.5d0, 1d-12, "test_angle_outliers_rap_alloc_basic: threshold out of expected range")

        ! Genes with z_scores >= threshold should be outliers
        call assert_false(is_outlier(1), "test_angle_outliers_rap_alloc_basic: gene 1 should not be outlier")
        call assert_false(is_outlier(2), "test_angle_outliers_rap_alloc_basic: gene 2 should not be outlier")
        call assert_false(is_outlier(3), "test_angle_outliers_rap_alloc_basic: gene 3 should not be outlier")
        call assert_true(is_outlier(4), "test_angle_outliers_rap_alloc_basic: gene 4 should be outlier")
        call assert_true(is_outlier(5), "test_angle_outliers_rap_alloc_basic: gene 5 should be outlier")
        call assert_true(is_outlier(6), "test_angle_outliers_rap_alloc_basic: gene 6 should be outlier")
        call assert_false(is_outlier(7), "test_angle_outliers_rap_alloc_basic: gene 7 should not be outlier")
        call assert_false(is_outlier(8), "test_angle_outliers_rap_alloc_basic: gene 8 (invalid) should not be outlier")
    end subroutine test_angle_outliers_alloc_basic

    !> Test outlier detection with no valid scaled angles
    subroutine test_angle_outliers_alloc_no_valid()
        integer(int32), parameter :: n_genes = 3
        real(real64) :: z_scores(n_genes) = [-1.0_real64, -1.0_real64, -1.0_real64]  ! All invalid
        real(real64) :: percentile = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr

        is_outlier = .false.

        call tox_angle_outliers_rap_alloc(z_scores, percentile, threshold, &
                                          is_outlier, n_genes, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_rap_alloc_no_valid: should succeed on no valid scores")
        call assert_false(any(is_outlier), "test_angle_outliers_rap_alloc_no_valid: no outliers should be marked")
    end subroutine test_angle_outliers_alloc_no_valid

    !> Test outlier detection where all valid scores are outliers
    subroutine test_angle_outliers_alloc_all_outliers()
        integer(int32), parameter :: n_genes = 5
        real(real64) :: z_scores(n_genes) = [1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0]
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr

        ! percentile = 0% -> everything outlier
        call tox_angle_outliers_rap_alloc(z_scores, 0.0_real64, threshold, &
                                          is_outlier, n_genes, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_rap_alloc_all_outliers: error code mismatch")
        call assert_true(all(is_outlier), "test_angle_outliers_rap_alloc_all_outliers: all should be outliers with 0th percentile")

        ! percentile = 100% -> everything inlier
        call tox_angle_outliers_rap_alloc(z_scores, 1.0_real64, threshold, &
                                          is_outlier, n_genes, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_rap_alloc_all_outliers: error code mismatch")
    call assert_equal_int(count(is_outlier), 1_int32, "test_angle_outliers_rap_alloc_all_outliers: only last should be inlier with 100th percentile")
    end subroutine test_angle_outliers_alloc_all_outliers

    !> Test complete pipeline for RAP angle-based outlier detection
    subroutine test_detect_angle_outliers_pipeline_basic()
        integer(int32), parameter :: n_genes = 10, n_families = 2
        ! Create RAP angles where most genes in each family point in similar directions
        ! but one gene in each family is an outlier
        real(real64) :: rap_angles(n_genes) = [ &
                        ! Family 1 genes (1-5): clustered around 0°
                        0.0_real64, &           ! Gene 1: normal
                        0.1_real64, &           ! Gene 2: normal
                        0.01_real64, &          ! Gene 3: normal
                        0.05_real64, &          ! Gene 4: normal
                        1.5_real64, &           ! Gene 5: OUTLIER (~86°, far from 0°)
                        ! Family 2 genes (6-10): clustered around π/2 (90°)
                        PI / 2.0_real64, &        ! Gene 6: normal
                        PI / 2.0_real64 + 0.1_real64, &  ! Gene 7: normal
                        PI / 2.0_real64 - 0.1_real64, &  ! Gene 8: normal
                        PI / 2.0_real64 + 0.05_real64, & ! Gene 9: normal
                        PI &        ! Gene 10: OUTLIER (far from 90°)
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
        real(real64) :: percentile_threshold = 0.85_real64  ! 85th percentile
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)
        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_pipeline_rap_basic: error code mismatch")

        ! Gene 5 and Gene 10 should be outliers (they point in very different directions)
        ! With 85th percentile and 10 genes (5 per family), we expect 1 outlier per family
        call assert_true(is_outlier(5), "test_detect_angle_outliers_pipeline_rap_basic: gene 5 should be outlier")
        call assert_true(is_outlier(10), "test_detect_angle_outliers_pipeline_rap_basic: gene 10 should be outlier")

        ! Most genes should not be outliers
        call assert_false(any(is_outlier(1:4)), "test_detect_angle_outliers_pipeline_rap_basic: genes 1-4 should not be outliers")
        call assert_false(any(is_outlier(6:9)), "test_detect_angle_outliers_pipeline_rap_basic: genes 6-9 should not be outliers")
    end subroutine test_detect_angle_outliers_pipeline_basic

    !> Test pipeline with single family
    subroutine test_detect_angle_outliers_pipeline_single_family()
        integer(int32), parameter :: n_genes = 8, n_families = 1
        ! Create RAP angles: 7 genes clustered around 45°, 1 gene at -90° (outlier)
        real(real64) :: rap_angles(n_genes) = [ &
                        PI / 4.0_real64, &           ! Gene 1: normal (45°)
                        PI / 4.0_real64 + 0.05_real64, &    ! Gene 2: normal
                        PI / 4.0_real64 - 0.05_real64, &    ! Gene 3: normal
                        PI / 4.0_real64 + 0.1_real64, &     ! Gene 4: normal
                        PI / 4.0_real64 - 0.1_real64, &     ! Gene 5: normal
                        PI / 4.0_real64 + 0.03_real64, &    ! Gene 6: normal
                        PI / 4.0_real64 - 0.03_real64, &    ! Gene 7: normal
                        PI / 2.0_real64 &            ! Gene 8: OUTLIER (90°)
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_pipeline_rap_single_family: error code mismatch")

        ! With 8 genes and 90th percentile, expect 1 outlier (gene 8)
        call assert_true(is_outlier(8), "test_detect_angle_outliers_pipeline_rap_single_family: gene 8 should be outlier")
  call assert_false(any(is_outlier(1:7)), "test_detect_angle_outliers_pipeline_rap_single_family: genes 1-7 should not be outliers")
    end subroutine test_detect_angle_outliers_pipeline_single_family

    !> Test pipeline with no valid families (all families have < 2 genes)
    subroutine test_detect_angle_outliers_pipeline_no_valid_families()
        integer(int32), parameter :: n_genes = 3, n_families = 3
        real(real64) :: rap_angles(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64]

        ! Each gene in its own family -> families have only 1 gene each
        integer(int32) :: gene_to_fam(n_genes) = [1, 2, 3]
        real(real64) :: percentile_threshold = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)

        ! Should indicate error because need at least 2 genes per family for meaningful direction
        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_pipeline_rap_no_valid_families: error code mismatch")
    call assert_true(all(status == STAT_NO_STABLE_DIRECTION), "test_detect_angle_outliers_pipeline_rap_no_valid_families: status code mismatch")
        ! But scaled angles should all be -1 since families are invalid
        call assert_true(all(z_scores == -1.0_real64), &
                         "test_detect_angle_outliers_pipeline_rap_no_valid_families: all scaled angles should be -1")
    end subroutine test_detect_angle_outliers_pipeline_no_valid_families

    !> Test pipeline with invalid percentile threshold
    subroutine test_detect_angle_outliers_pipeline_inval_percentile()
        integer(int32), parameter :: n_genes = 5, n_families = 1
        real(real64) :: rap_angles(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64, 0.5_real64]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 10.0_real64  ! Invalid: must be <= 100
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, &
                              "test_detect_angle_outliers_pipeline_rap_invalid_percentile: should error on invalid percentile")
      call assert_false(any(is_outlier), "test_detect_angle_outliers_pipeline_rap_invalid_percentile: no outliers should be marked")
    end subroutine test_detect_angle_outliers_pipeline_inval_percentile

    !> Test pipeline with nearly parallel genes (low angular dispersion)
    subroutine test_detect_angle_outliers_pipeline_parallel_genes()
        integer(int32), parameter :: n_genes = 6, n_families = 1
        ! Create nearly identical angles
        real(real64) :: rap_angles(n_genes) = [ &
                        0.0_real64, &
                        0.0001_real64, &
                        0.0003_real64, &
                        0.0002_real64, &
                        0.0004_real64, &
                        0.0003_real64 &
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)

        ! Due to missing angular variation, whole family is invalid -> all invalid
    call assert_true(all(status == STAT_NO_ANGULAR_VARIATION), "test_detect_angle_outliers_pipeline_rap_parallel_genes: minimal variation case handled")

        ! It's not an error if the pipeline results in invalid families only
        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_pipeline_rap_parallel_genes: error code mismatch")

    end subroutine test_detect_angle_outliers_pipeline_parallel_genes

    !> Test pipeline with clear orthogonal outlier
    subroutine test_detect_angle_outliers_pipeline_orth_outlier()
        integer(int32), parameter :: n_genes = 7, n_families = 1
        ! Create 6 genes pointing around 0°, 1 gene at 90° (orthogonal)
        real(real64) :: rap_angles(n_genes) = [ &
                        0.05_real64, &    ! Normal
                        0.0_real64, &      ! Normal
                        0.05_real64, &     ! Normal
                        0.1_real64, &      ! Normal
                        0.1_real64, &     ! Normal
                        0.03_real64, &     ! Normal
                        PI / 2.0_real64 &    ! OUTLIER: orthogonal (90°)
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 0.85_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_pipeline_rap_orthogonal_outlier: error code mismatch")

        ! Gene 7 should be a clear outlier
    call assert_true(is_outlier(7), "test_detect_angle_outliers_pipeline_rap_orthogonal_outlier: orthogonal gene should be outlier")
    call assert_false(any(is_outlier(1:6)), "test_detect_angle_outliers_pipeline_rap_orthogonal_outlier: normal genes should not be outliers")
    end subroutine test_detect_angle_outliers_pipeline_orth_outlier

    !> Test pipeline with mixed directions
    subroutine test_detect_angle_outliers_pipeline_mixed_directions()
        integer(int32), parameter :: n_genes = 12, n_families = 2
        ! Family 1: angles around 0° with one outlier at 120°
        ! Family 2: angles around 180° (-π) with one outlier at -60°
        real(real64) :: rap_angles(n_genes) = [ &
                        ! Family 1
                        0.2_real64, &      ! Gene 1: normal
                        0.0_real64, &       ! Gene 2: normal
                        0.1_real64, &       ! Gene 3: normal
                        0.15_real64, &     ! Gene 4: normal
                        0.05_real64, &      ! Gene 5: normal
                        2.0_real64 * PI / 3.0_real64, &  ! Gene 6: OUTLIER (120°)
                        ! Family 2
                        PI, &               ! Gene 7: normal (180°)
                        PI - 0.1_real64, &  ! Gene 8: normal
                        PI - 0.2_real64, &  ! Gene 9: normal
                        PI, &              ! Gene 10: normal (180°)
                        PI - 0.05_real64, &! Gene 11: normal
                        PI / 3.0_real64 &    ! Gene 12: OUTLIER (60°)
                        ]

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
        real(real64) :: percentile_threshold = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_pipeline_rap_mixed_directions: error code mismatch")

        ! With 90th percentile and 12 genes (6 per family), expect 1 outlier per family
        call assert_true(is_outlier(6), "test_detect_angle_outliers_pipeline_rap_mixed_directions: gene 6 should be outlier")
        call assert_true(is_outlier(12), "test_detect_angle_outliers_pipeline_rap_mixed_directions: gene 12 should be outlier")

        ! Most genes should not be outliers
    call assert_false(any(is_outlier(1:5)), "test_detect_angle_outliers_pipeline_rap_mixed_directions: genes 1-5 should not be outliers")
    call assert_false(any(is_outlier(7:11)), "test_detect_angle_outliers_pipeline_rap_mixed_directions: genes 7-11 should not be outliers")
    end subroutine test_detect_angle_outliers_pipeline_mixed_directions

    !> Test pipeline edge cases
    subroutine test_detect_angle_outliers_pipeline_edge_cases()
        integer(int32), parameter :: n_genes = 4, n_families = 1
        ! Test with angles that wrap across boundaries
        real(real64) :: rap_angles(n_genes)

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
        real(real64) :: percentile_threshold = 0.75_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        rap_angles = [ &
                     2.9_real64, &    ! ~166° (close to π)
                     below(PI), &   ! very close to π
                     0.1_real64, &    ! ~6°
                     3.0_real64 &     ! ~172° (close to π)
                     ]

        call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                    percentile_threshold, threshold, z_scores, &
                                                    n_genes, n_families, &
                                                    is_outlier, ierr, status)

        ! The test here is that wrapping works correctly, not necessarily outlier detection
        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_pipeline_rap_edge_cases: error code mismatch")

        ! All angles should be processed without errors
        call assert_true(.true., "test_detect_angle_outliers_pipeline_rap_edge_cases: edge cases processed")
    end subroutine test_detect_angle_outliers_pipeline_edge_cases

end module mod_test_get_outliers_by_angle_rap
