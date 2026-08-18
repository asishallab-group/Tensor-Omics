!> Unit test suite for tox_get_outliers_by_angle module.
module mod_test_get_outliers_by_angle
    use tox_get_outliers_by_angle
    use f42_utils, only: angle_between, wrap_angle, PI
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

    !> Get array of all available tests.
    function get_all_tests_outliers_by_angle() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(23))

        all_tests(1) = test_case("test_angle_between_vectors_basic", test_angle_between_vectors_basic)
        all_tests(2) = test_case("test_angle_between_vectors_parallel", test_angle_between_vectors_parallel)
        all_tests(3) = test_case("test_angle_between_vectors_orthogonal", test_angle_between_vectors_orthogonal)
        all_tests(4) = test_case("test_angle_between_vectors_opposite", test_angle_between_vectors_opposite)
        all_tests(5) = test_case("test_angle_between_vectors_zero_vector", test_angle_between_vectors_zero_vector)
        all_tests(6) = test_case("test_angle_between_vectors_empty", test_angle_between_vectors_empty)
        all_tests(7) = test_case("test_compute_family_direction_basic", test_compute_family_direction_basic)
        all_tests(8) = test_case("test_compute_family_direction_single_family", test_compute_family_direction_single_family)
    all_tests(9) = test_case("test_compute_family_direction_no_stable_direction", test_compute_family_direction_no_stable_direction)
       all_tests(10) = test_case("test_compute_family_direction_minimal_variation", test_compute_family_direction_minimal_variation)
        all_tests(11) = test_case("test_compute_angles_to_direction_basic", test_compute_angles_to_direction_basic)
       all_tests(12) = test_case("test_compute_angles_to_direction_invalid_family", test_compute_angles_to_direction_invalid_family)
        all_tests(13) = test_case("test_z_scores_by_dispersion_basic", test_z_scores_by_dispersion_basic)
        all_tests(14) = test_case("test_z_scores_by_dispersion_invalid", test_z_scores_by_dispersion_invalid)
        all_tests(15) = test_case("test_angle_outliers_alloc_basic", test_angle_outliers_alloc_basic)
        all_tests(16) = test_case("test_angle_outliers_alloc_no_valid", test_angle_outliers_alloc_no_valid)
        all_tests(17) = test_case("test_angle_outliers_alloc_all_outliers", test_angle_outliers_alloc_all_outliers)
        all_tests(18) = test_case("test_detect_angle_outliers_basic", test_detect_angle_outliers_basic)
        all_tests(19) = test_case("test_detect_angle_outliers_single_family", test_detect_angle_outliers_single_family)
        all_tests(20) = test_case("test_detect_angle_outliers_no_valid_families", test_detect_angle_outliers_no_valid_families)
        all_tests(21) = test_case("test_detect_angle_outliers_invalid_percentile", test_detect_angle_outliers_invalid_percentile)
        all_tests(22) = test_case("test_detect_angle_outliers_parallel_genes", test_detect_angle_outliers_parallel_genes)
        all_tests(23) = test_case("test_detect_angle_outliers_orthogonal_outlier", test_detect_angle_outliers_orthogonal_outlier)
    end function get_all_tests_outliers_by_angle

    !> Test basic angle calculation between vectors
    subroutine test_angle_between_vectors_basic()
        real(real64), dimension(3) :: vec1 = [1.0d0, 0.0d0, 0.0d0]
        real(real64), dimension(3) :: vec2 = [0.0d0, 1.0d0, 0.0d0]
        real(real64) :: angle
        integer(int32) :: ierr

        call angle_between(vec1, vec2, 3, angle, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_between_vectors_basic: error code mismatch")
        call assert_equal_real(angle, PI / 2.0d0, 1d-12, "test_angle_between_vectors_basic: angle mismatch")
    end subroutine test_angle_between_vectors_basic

    !> Test angle calculation for parallel vectors
    subroutine test_angle_between_vectors_parallel()
        real(real64), dimension(3) :: vec1 = [1.0d0, 2.0d0, 3.0d0]
        real(real64), dimension(3) :: vec2 = [2.0d0, 4.0d0, 6.0d0]  ! Parallel to vec1
        real(real64) :: angle
        integer(int32) :: ierr

        call angle_between(vec1, vec2, 3, angle, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_between_vectors_parallel: error code mismatch")
        call assert_equal_real(angle, 0.0d0, 1d-12, "test_angle_between_vectors_parallel: angle mismatch")
    end subroutine test_angle_between_vectors_parallel

    !> Test angle calculation for orthogonal vectors
    subroutine test_angle_between_vectors_orthogonal()
        real(real64), dimension(3) :: vec1 = [1.0d0, 0.0d0, 0.0d0]
        real(real64), dimension(3) :: vec2 = [0.0d0, 1.0d0, 0.0d0]
        real(real64) :: angle
        integer(int32) :: ierr

        call angle_between(vec1, vec2, 3, angle, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_between_vectors_orthogonal: error code mismatch")
        call assert_equal_real(angle, PI / 2.0d0, 1d-12, "test_angle_between_vectors_orthogonal: angle mismatch")
    end subroutine test_angle_between_vectors_orthogonal

    !> Test angle calculation for opposite vectors
    subroutine test_angle_between_vectors_opposite()
        real(real64), dimension(3) :: vec1 = [1.0d0, 0.0d0, 0.0d0]
        real(real64), dimension(3) :: vec2 = [-1.0d0, 0.0d0, 0.0d0]
        real(real64) :: angle
        integer(int32) :: ierr

        call angle_between(vec1, vec2, 3, angle, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_between_vectors_opposite: error code mismatch")
        call assert_equal_real(angle, PI, 1d-12, "test_angle_between_vectors_opposite: angle mismatch")
    end subroutine test_angle_between_vectors_opposite

    !> Test angle calculation with zero vector
    subroutine test_angle_between_vectors_zero_vector()
        real(real64), dimension(3) :: vec1 = [1.0d0, 2.0d0, 3.0d0]
        real(real64), dimension(3) :: vec2 = [0.0d0, 0.0d0, 0.0d0]
        real(real64) :: angle
        integer(int32) :: ierr

        call angle_between(vec1, vec2, 3, angle, ierr)

        call assert_equal_int(ierr, ERR_DIVISION_BY_ZERO, "test_angle_between_vectors_zero_vector: error code mismatch")
        call assert_equal_real(angle, 0.0d0, 1d-12, "test_angle_between_vectors_zero_vector: angle mismatch")
    end subroutine test_angle_between_vectors_zero_vector

    !> Test angle calculation with empty vectors
    subroutine test_angle_between_vectors_empty()
        real(real64), dimension(0) :: vec1, vec2
        real(real64) :: angle
        integer(int32) :: ierr

        call angle_between(vec1, vec2, 0, angle, ierr)

        call assert_equal_int(ierr, ERR_DIVISION_BY_ZERO, "test_angle_between_vectors_empty: error code mismatch")
        call assert_equal_real(angle, 0.0d0, 1d-12, "test_angle_between_vectors_empty: angle mismatch")
    end subroutine test_angle_between_vectors_empty

    !> Test family direction computation with multiple families
    subroutine test_compute_family_direction_basic()
        integer(int32), parameter :: n_samples = 3, n_genes = 6, n_families = 2
        real(real64) :: unit_vectors(n_samples, n_genes) = reshape([ &
                                                                   ! Family 1: Three vectors pointing in similar direction
                                                                   1.0d0, 0.0d0, 0.0d0, &    ! Gene 1
                                                                   0.9d0, 0.1d0, 0.0d0, &    ! Gene 2
                                                                   0.8d0, 0.2d0, 0.0d0, &    ! Gene 3
                                                                   ! Family 2: Three vectors pointing in another direction
                                                                   0.0d0, 1.0d0, 0.0d0, &    ! Gene 4
                                                                   0.1d0, 0.9d0, 0.0d0, &    ! Gene 5
                                                                   0.2d0, 0.8d0, 0.0d0 &     ! Gene 6
                                                                   ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 2, 2, 2]
        integer(int32) :: tmp_member_counts(n_families), expected_member_counts(n_families)
        real(real64) :: family_directions(n_samples, n_families)
        real(real64) :: angular_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: norm_fam1, norm_fam2

        call tox_compute_family_direction(unit_vectors, gene_to_fam, &
                                          family_directions, angular_dispersions, &
                                          n_samples, n_genes, n_families, ierr, status, tmp_member_counts)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_basic: error code mismatch")
        call assert_true(all(status == ERR_OK), "test_compute_family_direction_basic: status code mismatch")

        expected_member_counts = [3, 3]
    call assert_equal_array_int(tmp_member_counts, expected_member_counts, n_families, "test_compute_family_direction_basic: member count mismatch")

        ! Check family directions are unit vectors
        norm_fam1 = norm2(family_directions(:, 1))
        norm_fam2 = norm2(family_directions(:, 2))
        call assert_equal_real(norm_fam1, 1.0d0, 1d-12, "test_compute_family_direction_basic: family 1 direction not normalized")
        call assert_equal_real(norm_fam2, 1.0d0, 1d-12, "test_compute_family_direction_basic: family 2 direction not normalized")

        ! Check angular dispersions are positive and reasonable
        call assert_true(angular_dispersions(1) > 0.0d0 .and. angular_dispersions(1) < PI, &
                         "test_compute_family_direction_basic: family 1 dispersion out of range")
        call assert_true(angular_dispersions(2) > 0.0d0 .and. angular_dispersions(2) < PI, &
                         "test_compute_family_direction_basic: family 2 dispersion out of range")

        ! Family 1 should point roughly in [1,0,0] direction
        call assert_true(family_directions(1, 1) > 0.9d0, "test_compute_family_direction_basic: family 1 direction incorrect")
        call assert_true(abs(family_directions(2, 1)) < 0.2d0, "test_compute_family_direction_basic: family 1 direction incorrect")
    end subroutine test_compute_family_direction_basic

    !> Test family direction with single family
    subroutine test_compute_family_direction_single_family()
        integer(int32), parameter :: n_samples = 3, n_genes = 4, n_families = 1
        real(real64) :: unit_vectors(n_samples, n_genes) = reshape([ &
                                                                   1.0d0, 0.0d0, 0.0d0, &
                                                                   0.9d0, 0.1d0, 0.0d0, &
                                                                   0.8d0, 0.2d0, 0.0d0, &
                                                                   0.7d0, 0.3d0, 0.0d0 &
                                                                   ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
        integer(int32) :: tmp_member_counts(n_families), expected_member_counts(n_families)
        real(real64) :: family_directions(n_samples, n_families)
        real(real64) :: angular_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        call tox_compute_family_direction(unit_vectors, gene_to_fam, &
                                          family_directions, angular_dispersions, &
                                          n_samples, n_genes, n_families, ierr, status, tmp_member_counts)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_single_family: error code mismatch")
      call assert_true(angular_dispersions(1) > 0.0d0, "test_compute_family_direction_single_family: dispersion should be positive")
        expected_member_counts = [4]
    call assert_equal_array_int(tmp_member_counts, expected_member_counts, n_families, "test_compute_family_direction_single_family: member count mismatch")
    end subroutine test_compute_family_direction_single_family

    !> Test family direction with no stable direction (vectors cancel out)
    subroutine test_compute_family_direction_no_stable_direction()
        integer(int32), parameter :: n_samples = 3, n_genes = 4, n_families = 1
        ! Create vectors that point in opposite directions, cancelling out
        real(real64) :: unit_vectors(n_samples, n_genes) = reshape([ &
                                                                   1.0d0, 0.0d0, 0.0d0, &
                                                                   -1.0d0, 0.0d0, 0.0d0, &
                                                                   0.0d0, 1.0d0, 0.0d0, &
                                                                   0.0d0, -1.0d0, 0.0d0 &
                                                                   ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
        integer(int32) :: tmp_member_counts(n_families), expected_member_counts(n_families)
        real(real64) :: family_directions(n_samples, n_families)
        real(real64) :: angular_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        call tox_compute_family_direction(unit_vectors, gene_to_fam, &
                                          family_directions, angular_dispersions, &
                                          n_samples, n_genes, n_families, ierr, status, tmp_member_counts)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_no_stable_direction: error code mismatch")
        expected_member_counts = [4]
    call assert_equal_array_int(tmp_member_counts, expected_member_counts, n_families, "test_compute_family_direction_no_stable_direction: member count mismatch")
        call assert_equal_real(angular_dispersions(1), -1.0d0, 1d-12, &
                               "test_compute_family_direction_no_stable_direction: should have -1 dispersion")
        ! Status should indicate no stable direction
        call assert_true(all(status == STAT_NO_STABLE_DIRECTION), &
                         "test_compute_family_direction_no_stable_direction: status should indicate no stable direction")
    end subroutine test_compute_family_direction_no_stable_direction

    !> Test family direction with minimal angular variation
    subroutine test_compute_family_direction_minimal_variation()
        integer(int32), parameter :: n_samples = 3, n_genes = 4, n_families = 1
        ! Create nearly identical vectors (very small angular dispersion)
        real(real64) :: unit_vectors(n_samples, n_genes) = reshape([ &
                                                                   1.0d0, 0.0d0, 0.0d0, &
                                                                   0.9999999d0, 0.0001d0, 0.0d0, &
                                                                   0.9999999d0, -0.0001d0, 0.0d0, &
                                                                   0.9999998d0, 0.0002d0, 0.0d0 &
                                                                   ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
        real(real64) :: family_directions(n_samples, n_families)
        real(real64) :: angular_dispersions(n_families)
        integer(int32) :: ierr, status(n_families)

        call tox_compute_family_direction_alloc(unit_vectors, gene_to_fam, &
                                                family_directions, angular_dispersions, &
                                                n_samples, n_genes, n_families, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_compute_family_direction_minimal_variation: error code mismatch")
        ! These nearly identical vectors have a dispersion well below the minimum, so the
        ! family is deterministically flagged as having no angular variation (sentinel dispersion).
        call assert_true(angular_dispersions(1) < 0.0d0, &
                         "test_compute_family_direction_minimal_variation: dispersion should be the sentinel")
        call assert_true(all(status == STAT_NO_ANGULAR_VARIATION), &
                         "test_compute_family_direction_minimal_variation: status should indicate minimal variation")
    end subroutine test_compute_family_direction_minimal_variation

    !> Test angle computation to family directions
    subroutine test_compute_angles_to_direction_basic()
        integer(int32), parameter :: n_samples = 3, n_genes = 4, n_families = 2
        real(real64) :: unit_vectors(n_samples, n_genes) = reshape([ &
                                                                   1.0d0, 0.0d0, 0.0d0, &  ! Gene 1, Family 1
                                                                   0.8d0, 0.6d0, 0.0d0, &  ! Gene 2, Family 1
                                                                   0.0d0, 1.0d0, 0.0d0, &  ! Gene 3, Family 2
                                                                   0.6d0, 0.8d0, 0.0d0 &   ! Gene 4, Family 2
                                                                   ], [n_samples, n_genes])

        ! Family directions: Family 1 points to [1,0,0], Family 2 points to [0,1,0]
        real(real64) :: family_directions(n_samples, n_families) = reshape([ &
                                                                           1.0d0, 0.0d0, 0.0d0, &
                                                                           0.0d0, 1.0d0, 0.0d0 &
                                                                           ], [n_samples, n_families])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 2, 2]
        real(real64) :: angles(n_genes)
        integer(int32) :: ierr

        call tox_compute_angles_to_direction(unit_vectors, family_directions, &
                                             gene_to_fam, angles, &
                                             n_samples, n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_angles_to_direction_basic: error code mismatch")

        ! Gene 1: parallel to family direction -> angle = 0
        call assert_equal_real(angles(1), 0.0d0, 1d-12, "test_compute_angles_to_direction_basic: gene 1 angle mismatch")

        ! Gene 2: angle between [0.8,0.6,0] and [1,0,0] = arccos(0.8) ≈ 0.6435 rad
        call assert_equal_real(angles(2), acos(0.8d0), 1d-12, "test_compute_angles_to_direction_basic: gene 2 angle mismatch")

        ! Gene 3: parallel to family direction -> angle = 0
        call assert_equal_real(angles(3), 0.0d0, 1d-12, "test_compute_angles_to_direction_basic: gene 3 angle mismatch")

        ! Gene 4: angle between [0.6,0.8,0] and [0,1,0] = arccos(0.8) ≈ 0.6435 rad
        call assert_equal_real(angles(4), acos(0.8d0), 1d-12, "test_compute_angles_to_direction_basic: gene 4 angle mismatch")
    end subroutine test_compute_angles_to_direction_basic

    !> Test angle computation with invalid family mapping
    subroutine test_compute_angles_to_direction_invalid_family()
        integer(int32), parameter :: n_samples = 3, n_genes = 4, n_families = 2
        real(real64) :: unit_vectors(n_samples, n_genes) = reshape([ &
                                                                   1.0d0, 0.0d0, 0.0d0, &
                                                                   0.0d0, 1.0d0, 0.0d0, &
                                                                   0.0d0, 0.0d0, 1.0d0, &
                                                                   1.0d0 / sqrt(3.0d0), 1.0d0 / sqrt(3.0d0), 1.0d0 / sqrt(3.0d0) &
                                                                   ], [n_samples, n_genes])

        real(real64) :: family_directions(n_samples, n_families) = reshape([ &
                                                                           1.0d0, 0.0d0, 0.0d0, &
                                                                           0.0d0, 1.0d0, 0.0d0 &
                                                                           ], [n_samples, n_families])

        ! Gene 2: unassigned (0), Gene 3: invalid family (3 > n_families)
        integer(int32) :: gene_to_fam(n_genes) = [1, 0, 3, 2]
        real(real64) :: angles(n_genes)
        integer(int32) :: ierr

        call tox_compute_angles_to_direction(unit_vectors, family_directions, &
                                             gene_to_fam, angles, &
                                             n_samples, n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_compute_angles_to_direction_invalid_family: error code mismatch")
    end subroutine test_compute_angles_to_direction_invalid_family

    !> Test scaling angles by dispersion (z-scores)
    subroutine test_z_scores_by_dispersion_basic()
        integer(int32), parameter :: n_genes = 5, n_families = 2
        real(real64) :: angles(n_genes) = [0.2d0, 0.4d0, -1.0d0, 0.3d0, 0.6d0]
        real(real64) :: angular_dispersions(n_families) = [0.1d0, 0.2d0]  ! Family 1: σ=0.1, Family 2: σ=0.2
        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 0, 2, 2]  ! Gene 3: unassigned
        real(real64) :: z_scores(n_genes)
        integer(int32) :: ierr

        call tox_z_scores_by_dispersion(angles, angular_dispersions, &
                                        gene_to_fam, z_scores, &
                                        n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_z_scores_by_dispersion_basic: error code mismatch")

        ! Gene 1: 0.2 / 0.1 = 2.0
        call assert_equal_real(z_scores(1), 2.0d0, 1d-12, "test_z_scores_by_dispersion_basic: gene 1 mismatch")

        ! Gene 2: 0.4 / 0.1 = 4.0
        call assert_equal_real(z_scores(2), 4.0d0, 1d-12, "test_z_scores_by_dispersion_basic: gene 2 mismatch")

        ! Gene 3: unassigned -> -1.0
        call assert_equal_real(z_scores(3), -1.0d0, 1d-12, "test_z_scores_by_dispersion_basic: gene 3 mismatch")

        ! Gene 4: 0.3 / 0.2 = 1.5
        call assert_equal_real(z_scores(4), 1.5d0, 1d-12, "test_z_scores_by_dispersion_basic: gene 4 mismatch")

        ! Gene 5: 0.6 / 0.2 = 3.0
        call assert_equal_real(z_scores(5), 3.0d0, 1d-12, "test_z_scores_by_dispersion_basic: gene 5 mismatch")
    end subroutine test_z_scores_by_dispersion_basic

    !> Test scaling angles with invalid inputs
    subroutine test_z_scores_by_dispersion_invalid()
        integer(int32), parameter :: n_genes = 4, n_families = 2
        real(real64) :: angles(n_genes) = [-1.0d0, 0.5d0, 0.3d0, 0.4d0]
        real(real64) :: angular_dispersions(n_families) = [-1.0d0, 0.0d0]  ! Family 1: invalid, Family 2: zero
        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 3, 2]  ! Gene 3: invalid family
        real(real64) :: z_scores(n_genes)
        integer(int32) :: ierr

        call tox_z_scores_by_dispersion(angles, angular_dispersions, &
                                        gene_to_fam, z_scores, &
                                        n_genes, n_families, ierr)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_z_scores_by_dispersion_invalid: error code mismatch")
    end subroutine test_z_scores_by_dispersion_invalid

    !> Test outlier detection based on scaled angles
    subroutine test_angle_outliers_alloc_basic()
        integer(int32), parameter :: n_genes = 8
        ! Create scaled angles: some normal, some outliers
        real(real64) :: z_scores(n_genes) = [10.0d0, 10.5d0, 20.0d0, 20.5d0, 30.0d0, 20.5d0, 20.0d0, -1.0d0]
        real(real64) :: percentile = 0.8_real64  ! 80th percentile
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr

        call tox_angle_outliers_alloc(z_scores, percentile, threshold, &
                                      is_outlier, n_genes, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_alloc_basic: error code mismatch")

        ! With 7 valid scores, 80th percentile should be value at index 5 + 80% of diff to next larger one (0.8*(7-1)+1 = 5.8)
        ! Sorted: 1.0, 1.5, 2.0, 2.0, 2.5, 2.5, 3.0
        call assert_equal_real(threshold, 20.5d0, 1d-12, "test_angle_outliers_alloc_basic: threshold out of expected range")

        ! Genes with z_scores >= threshold should be outliers
        call assert_false(is_outlier(1), "test_angle_outliers_alloc_basic: gene 1 should not be outlier")
        call assert_false(is_outlier(2), "test_angle_outliers_alloc_basic: gene 2 should not be outlier")
        call assert_false(is_outlier(3), "test_angle_outliers_alloc_basic: gene 3 should not be outlier")
        call assert_true(is_outlier(4), "test_angle_outliers_alloc_basic: gene 4 should be outlier")
        call assert_true(is_outlier(5), "test_angle_outliers_alloc_basic: gene 5 should be outlier")
        call assert_true(is_outlier(6), "test_angle_outliers_alloc_basic: gene 6 should be outlier")
        call assert_false(is_outlier(7), "test_angle_outliers_alloc_basic: gene 7 should not be outlier")
        call assert_false(is_outlier(8), "test_angle_outliers_alloc_basic: gene 8 (invalid) should not be outlier")
    end subroutine test_angle_outliers_alloc_basic

    !> Test outlier detection with no valid scaled angles
    subroutine test_angle_outliers_alloc_no_valid()
        integer(int32), parameter :: n_genes = 3
        real(real64) :: z_scores(n_genes) = [-1.0d0, -1.0d0, -1.0d0]  ! All invalid
        real(real64) :: percentile = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr

        is_outlier = .false.

        call tox_angle_outliers_alloc(z_scores, percentile, threshold, &
                                      is_outlier, n_genes, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_alloc_no_valid: should succeed on no valid scores")
        call assert_false(any(is_outlier), "test_angle_outliers_alloc_no_valid: no outliers should be marked")
    end subroutine test_angle_outliers_alloc_no_valid

    !> Test outlier detection where all valid scores are outliers
    subroutine test_angle_outliers_alloc_all_outliers()
        integer(int32), parameter :: n_genes = 5
        real(real64) :: z_scores(n_genes) = [1.0d0, 1.1d0, 1.2d0, 1.3d0, 1.4d0]
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr

        ! percentile = 0% -> everything outlier
        call tox_angle_outliers_alloc(z_scores, 0.0_real64, threshold, &
                                      is_outlier, n_genes, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_alloc_all_outliers: error code mismatch")
        call assert_true(all(is_outlier), "test_angle_outliers_alloc_all_outliers: all should be outliers with 0th percentile")

        ! percentile = 100% -> everything inlier
        call tox_angle_outliers_alloc(z_scores, 1.0_real64, threshold, &
                                      is_outlier, n_genes, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_angle_outliers_alloc_all_outliers: error code mismatch")
        call assert_false(all(is_outlier), "test_angle_outliers_alloc_all_outliers: all should be inliers with 100th percentile")
    end subroutine test_angle_outliers_alloc_all_outliers

    !> Test complete pipeline for angle-based outlier detection
    subroutine test_detect_angle_outliers_basic()
        integer(int32), parameter :: n_samples = 3, n_genes = 10, n_families = 2
        ! Create expression vectors where most genes in each family point in similar directions
        ! but one gene in each family is an outlier
        real(real64) :: expression_vectors(n_samples, n_genes) = reshape([ &
                                                                         ! Family 1 genes (1-5)
                                                                         1.0d0, 0.1d0, 0.0d0, &    ! Gene 1: normal
                                                                         0.9d0, 0.2d0, 0.0d0, &    ! Gene 2: normal
                                                                         0.8d0, 0.3d0, 0.0d0, &    ! Gene 3: normal
                                                                         0.7d0, 0.4d0, 0.0d0, &    ! Gene 4: normal
                                                                         0.0d0, 1.0d0, 0.0d0, &    ! Gene 5: OUTLIER (orthogonal to others)
                                                                         ! Family 2 genes (6-10)
                                                                         0.1d0, 1.0d0, 0.1d0, &    ! Gene 6: normal
                                                                         0.2d0, 0.9d0, 0.1d0, &    ! Gene 7: normal
                                                                         0.3d0, 0.8d0, 0.1d0, &    ! Gene 8: normal
                                                                         0.4d0, 0.7d0, 0.1d0, &    ! Gene 9: normal
                                                                         1.0d0, 0.0d0, 0.0d0 &     ! Gene 10: OUTLIER (orthogonal to others)
                                                                         ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
        real(real64) :: percentile_threshold = 0.85_real64  ! 85th percentile
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                                percentile_threshold, threshold, z_scores, &
                                                n_samples, n_genes, n_families, &
                                                is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_basic: error code mismatch")

        ! Gene 5 and Gene 10 should be outliers (they point in very different directions)
        ! With 85th percentile and 10 genes (5 per family), we expect 1 outlier per family
        call assert_true(is_outlier(5), "test_detect_angle_outliers_basic: gene 5 should be outlier")
        call assert_true(is_outlier(10), "test_detect_angle_outliers_basic: gene 10 should be outlier")

        ! Most genes should not be outliers
        call assert_false(any(is_outlier(1:4)), "test_detect_angle_outliers_basic: genes 1-4 should not be outliers")
        call assert_false(any(is_outlier(6:9)), "test_detect_angle_outliers_basic: genes 6-9 should not be outliers")
    end subroutine test_detect_angle_outliers_basic

    !> Test pipeline with single family
    subroutine test_detect_angle_outliers_single_family()
        integer(int32), parameter :: n_samples = 3, n_genes = 8, n_families = 1
        real(real64) :: expression_vectors(n_samples, n_genes) = reshape([ &
                                                                         1.0d0, 0.1d0, 0.0d0, &
                                                                         0.9d0, 0.2d0, 0.0d0, &
                                                                         0.8d0, 0.3d0, 0.0d0, &
                                                                         0.7d0, 0.4d0, 0.0d0, &
                                                                         0.6d0, 0.5d0, 0.0d0, &
                                                                         0.5d0, 0.6d0, 0.0d0, &
                                                                         0.4d0, 0.7d0, 0.0d0, &
                                                                         0.0d0, 1.0d0, 0.0d0 &  ! OUTLIER: orthogonal to main direction
                                                                         ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                                percentile_threshold, threshold, z_scores, &
                                                n_samples, n_genes, n_families, &
                                                is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_single_family: error code mismatch")

        ! With 8 genes and 90th percentile, expect 1 outlier (gene 8)
        call assert_true(is_outlier(8), "test_detect_angle_outliers_single_family: gene 8 should be outlier")
        call assert_false(any(is_outlier(1:7)), "test_detect_angle_outliers_single_family: genes 1-7 should not be outliers")
    end subroutine test_detect_angle_outliers_single_family

    !> Test pipeline with no valid families (all families have < 2 genes)
    subroutine test_detect_angle_outliers_no_valid_families()
        integer(int32), parameter :: n_samples = 3, n_genes = 3, n_families = 3
        real(real64) :: expression_vectors(n_samples, n_genes) = reshape([ &
                                                                         1.0d0, 0.0d0, 0.0d0, &
                                                                         0.0d0, 1.0d0, 0.0d0, &
                                                                         0.0d0, 0.0d0, 1.0d0 &
                                                                         ], [n_samples, n_genes])

        ! Each gene in its own family -> families have only 1 gene each
        integer(int32) :: gene_to_fam(n_genes) = [1, 2, 3]
        real(real64) :: percentile_threshold = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                                percentile_threshold, threshold, z_scores, &
                                                n_samples, n_genes, n_families, &
                                                is_outlier, ierr, status)

        ! Should indicate error because need at least 2 genes per family for meaningful direction (error occurs in tox_angle_outliers_alloc)
        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_no_valid_families: error code mismatch")

        call assert_false(any(is_outlier), "test_detect_angle_outliers_no_valid_families: no outliers should be marked")
    end subroutine test_detect_angle_outliers_no_valid_families

    !> Test pipeline with invalid percentile threshold
    subroutine test_detect_angle_outliers_invalid_percentile()
        integer(int32), parameter :: n_samples = 3, n_genes = 5, n_families = 1
        real(real64) :: expression_vectors(n_samples, n_genes) = reshape([ &
                                                                         1.0d0, 0.0d0, 0.0d0, &
                                                                         0.9d0, 0.1d0, 0.0d0, &
                                                                         0.8d0, 0.2d0, 0.0d0, &
                                                                         0.7d0, 0.3d0, 0.0d0, &
                                                                         0.6d0, 0.4d0, 0.0d0 &
                                                                         ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 1.01_real64  ! Invalid: must be 0 <= x <= 100
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                                percentile_threshold, threshold, z_scores, &
                                                n_samples, n_genes, n_families, &
                                                is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, &
                              "test_detect_angle_outliers_invalid_percentile: should error on invalid percentile")
        call assert_false(any(is_outlier), "test_detect_angle_outliers_invalid_percentile: no outliers should be marked")
    end subroutine test_detect_angle_outliers_invalid_percentile

    !> Test pipeline with nearly parallel genes (low angular dispersion)
    subroutine test_detect_angle_outliers_parallel_genes()
        integer(int32), parameter :: n_samples = 3, n_genes = 6, n_families = 1
        ! Create nearly identical vectors
        real(real64) :: expression_vectors(n_samples, n_genes) = reshape([ &
                                                                         1.0d0, 0.001d0, 0.0d0, &
                                                                         0.999d0, 0.002d0, 0.0d0, &
                                                                         0.998d0, 0.001d0, 0.0d0, &
                                                                         0.999d0, 0.0d0, 0.001d0, &
                                                                         0.999d0, -0.001d0, 0.0d0, &
                                                                         1.0d0, -0.002d0, 0.0d0 &
                                                                         ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 0.9_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                                percentile_threshold, threshold, z_scores, &
                                                n_samples, n_genes, n_families, &
                                                is_outlier, ierr, status)

        ! Minimal-variation families are reported via `status`, not `ierr`, so the pipeline
        ! still succeeds. With such parallel genes the family is flagged as having no angular
        ! variation and no gene can be a directional outlier.
        call assert_equal_int(ierr, ERR_OK, &
                              "test_detect_angle_outliers_parallel_genes: error code should be OK")
        call assert_true(all(status == STAT_NO_ANGULAR_VARIATION), &
                         "test_detect_angle_outliers_parallel_genes: family should be flagged as minimal variation")
        call assert_false(any(is_outlier), &
                          "test_detect_angle_outliers_parallel_genes: with parallel genes, no outliers expected")
    end subroutine test_detect_angle_outliers_parallel_genes

    !> Test pipeline with clear orthogonal outlier
    subroutine test_detect_angle_outliers_orthogonal_outlier()
        integer(int32), parameter :: n_samples = 3, n_genes = 7, n_families = 1
        ! Create 6 genes pointing in similar direction, 1 gene orthogonal
        real(real64) :: expression_vectors(n_samples, n_genes) = reshape([ &
                                                                         1.0d0, 0.1d0, 0.0d0, &    ! Normal
                                                                         0.9d0, 0.2d0, 0.0d0, &    ! Normal
                                                                         0.8d0, 0.3d0, 0.0d0, &    ! Normal
                                                                         0.7d0, 0.4d0, 0.0d0, &    ! Normal
                                                                         0.6d0, 0.5d0, 0.0d0, &    ! Normal
                                                                         0.5d0, 0.6d0, 0.0d0, &    ! Normal
                                                                         0.0d0, 1.0d0, 0.0d0 &     ! OUTLIER: orthogonal
                                                                         ], [n_samples, n_genes])

        integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1, 1]
        real(real64) :: percentile_threshold = 0.85_real64
        real(real64) :: threshold
        logical :: is_outlier(n_genes)
        integer(int32) :: ierr, status(n_families)
        real(real64) :: z_scores(n_genes)

        call tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                                percentile_threshold, threshold, z_scores, &
                                                n_samples, n_genes, n_families, &
                                                is_outlier, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_detect_angle_outliers_orthogonal_outlier: error code mismatch")

        ! Gene 7 should be a clear outlier
        call assert_true(is_outlier(7), "test_detect_angle_outliers_orthogonal_outlier: orthogonal gene should be outlier")
       call assert_false(any(is_outlier(1:6)), "test_detect_angle_outliers_orthogonal_outlier: normal genes should not be outliers")
    end subroutine test_detect_angle_outliers_orthogonal_outlier

end module mod_test_get_outliers_by_angle
