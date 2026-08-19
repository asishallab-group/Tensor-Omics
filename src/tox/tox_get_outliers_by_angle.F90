#include <src/macros.h>

module tox_get_outliers_by_angle
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use f42_utils, only: sort_real_heapsort, calc_percentile_helper, PI, is_close, add_vector, angle_between, wrap_angle
    use tox_errors, only: STAT_NO_ANGULAR_VARIATION, STAT_NO_STABLE_DIRECTION, set_ok, set_err, validate_dimension_size, ERR_ALLOC_FAIL, &
                        validate_all_in_range_real, validate_all_in_range_int, validate_in_range_real, is_err, validate_in_range_int
    use tox_normalization, only: tox_normalize_vectors_unit_length_helper
    implicit none

    private
    public :: tox_compute_family_direction, &
              tox_compute_angles_to_direction, tox_z_scores_by_dispersion, &
              tox_angle_outliers, tox_detect_angle_outliers_pipeline
    public :: tox_compute_family_direction_alloc, tox_angle_outliers_alloc

#define CM_MIN_ANGULAR_DISPERSION_DEFAULT 1.0e-02_real64
#define CM_MAX_ANGULAR_DISPERSION_DEFAULT sqrt(-2.0_real64 * log(0.5_real64))

#define CM_ANGLES_SENTINEL -1.0_real64
#define CM_ANGULAR_DISPERSION_SENTINEL -1.0_real64
#define CM_Z_SCORES_SENTINEL -1.0_real64

contains

    !> Compute family reference direction (spherical mean) and angular dispersion.
    !| Some directions might result in sentinel values. In such case the respective family is being marked by the `status(i_family)` set to the respective status code.
    pure subroutine tox_compute_family_direction(unit_vectors, gene_to_fam, &
                                                 family_directions, angular_dispersions, &
                    n_samples, n_genes, n_families, ierr, status, tmp_member_counts, min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_samples
            !! Number of samples (dimension of vectors)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_samples, n_genes), intent(in) :: unit_vectors
            !! Unit vectors for each gene, shape (n_samples, n_genes)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_samples, n_families), intent(out) :: family_directions
            !! Family reference directions (spherical means)
        real(real64), dimension(n_families), intent(out) :: angular_dispersions
            !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_ANGULAR_DISPERSION_SENTINEL` if there is no stable direction or too less angular variation (`status`)
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                               |
            !! |-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                   |
            !! | STAT_NO_STABLE_DIRECTION  |             The `unit_vectors` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `unit_vectors` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        integer(int32), dimension(n_families), intent(out) :: tmp_member_counts
            !! Number of members per family, determined by `gene_to_fam`
        real(real64), intent(in), optional :: min_angular_dispersion
            !! The minimum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MIN_ANGULAR_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_angular_dispersion
            !! The maximum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MAX_ANGULAR_DISPERSION_DEFAULT`

        ! Input validation
        call set_ok(ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(unit_vectors, size(unit_vectors, kind=int32), ierr, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        call tox_compute_family_direction_helper(unit_vectors, gene_to_fam, &
                                                 family_directions, angular_dispersions, &
                          n_samples, n_genes, n_families, status, tmp_member_counts, min_angular_dispersion, max_angular_dispersion)
    end subroutine tox_compute_family_direction

    !> (no input validation) Compute family reference direction (spherical mean) and angular dispersion.
    !| Some directions might result in sentinel values. In such case the respective family is being marked by the `status(i_family)` set to the respective status code.
    pure subroutine tox_compute_family_direction_helper(unit_vectors, gene_to_fam, &
                                                        family_directions, angular_dispersions, &
                          n_samples, n_genes, n_families, status, tmp_member_counts, min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_samples
            !! Number of samples (dimension of vectors)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_samples, n_genes), intent(in) :: unit_vectors
            !! Unit vectors for each gene, shape (n_samples, n_genes)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_samples, n_families), intent(out) :: family_directions
            !! Family reference directions (spherical means)
        real(real64), dimension(n_families), intent(out) :: angular_dispersions
            !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_ANGULAR_DISPERSION_SENTINEL` if there is no stable direction or too less angular variation (`status`)
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                               |
            !! |-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                   |
            !! | STAT_NO_STABLE_DIRECTION  |             The `unit_vectors` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `unit_vectors` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        integer(int32), dimension(n_families), intent(out) :: tmp_member_counts
            !! Number of members per family, determined by `gene_to_fam`
        real(real64), intent(in), optional :: min_angular_dispersion
            !! The minimum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MIN_ANGULAR_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_angular_dispersion
            !! The maximum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MAX_ANGULAR_DISPERSION_DEFAULT`

        integer(int32) :: i_gene, i_family, family_idx
        real(real64) :: mean_length, angular_dispersion, fam_norm, actual_min_angular_dispersion, actual_max_angular_dispersion

        M_DEFAULT_VAL(min_angular_dispersion, actual_min_angular_dispersion, CM_MIN_ANGULAR_DISPERSION_DEFAULT)
        M_DEFAULT_VAL(max_angular_dispersion, actual_max_angular_dispersion, CM_MAX_ANGULAR_DISPERSION_DEFAULT)

        call set_ok(status)
        family_directions = 0.0_real64
        angular_dispersions = CM_ANGULAR_DISPERSION_SENTINEL
        tmp_member_counts = 0_int32

        do i_gene = 1, n_genes
            family_idx = gene_to_fam(i_gene)
            if (family_idx /= M_GENE_TO_FAM_SENTINEL) then
                tmp_member_counts(family_idx) = tmp_member_counts(family_idx) + 1
                call add_vector(family_directions(:, family_idx), unit_vectors(:, i_gene))
            end if
        end do

        ! Compute family directions and dispersions
        do concurrent(i_family=1:n_families) &
            local(fam_norm, mean_length, angular_dispersion) &
            shared (tmp_member_counts, angular_dispersions, status, family_directions, n_samples, actual_min_angular_dispersion, actual_max_angular_dispersion)

            if (tmp_member_counts(i_family) < 2) then
                angular_dispersions(i_family) = CM_ANGULAR_DISPERSION_SENTINEL
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                cycle
            end if

            ! Compute R_F = ‖sum‖ / count
            fam_norm = norm2(family_directions(:, i_family))

            ! A (near) zero resultant length means the unit vectors cancel out, so there
            ! is no stable mean direction. Skip before dividing by `fam_norm` to avoid
            ! polluting `family_directions` with NaN/huge values.
            if (is_close(fam_norm, 0.0_real64)) then
                angular_dispersions(i_family) = CM_ANGULAR_DISPERSION_SENTINEL
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                cycle
            end if

            mean_length = fam_norm / real(tmp_member_counts(i_family), real64)

            ! Compute family direction (spherical mean). `fam_norm > 0` is guaranteed by the
            ! zero-resultant check above, so the division is safe.
            family_directions(:, i_family) = family_directions(:, i_family) / fam_norm

            ! Edge case handling, will always exceed `max_angular_dispersion`, thus no stable direction
            if (mean_length == 0.0_real64) then
                angular_dispersions(i_family) = CM_ANGULAR_DISPERSION_SENTINEL
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
            else
                angular_dispersion = sqrt(-2.0_real64 * log(mean_length))

                ! Check for minimal angular variation
                if (angular_dispersion < actual_min_angular_dispersion) then
                    angular_dispersions(i_family) = CM_ANGULAR_DISPERSION_SENTINEL
                    call set_err(status(i_family), STAT_NO_ANGULAR_VARIATION)

                    ! Check for maximum angular variation
                else if (angular_dispersion > actual_max_angular_dispersion) then
                    angular_dispersions(i_family) = CM_ANGULAR_DISPERSION_SENTINEL
                    call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                else
                    angular_dispersions(i_family) = angular_dispersion
                end if
            end if
        end do

    end subroutine tox_compute_family_direction_helper

    !> Compute family reference direction (spherical mean) and angular dispersion.
    !| Some directions might result in sentinel values. In such case the respective family is being marked by the `status(i_family)` set to the respective status code.
    pure subroutine tox_compute_family_direction_alloc(unit_vectors, gene_to_fam, &
                                                       family_directions, angular_dispersions, &
                                       n_samples, n_genes, n_families, ierr, status, min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_samples
            !! Number of samples (dimension of vectors)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_samples, n_genes), intent(in) :: unit_vectors
            !! Unit vectors for each gene
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_samples, n_families), intent(out) :: family_directions
            !! Family reference directions
        real(real64), dimension(n_families), intent(out) :: angular_dispersions
            !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_ANGULAR_DISPERSION_SENTINEL` if there is no stable direction or too less angular variation (`status`)
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                               |
            !! |-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                   |
            !! | STAT_NO_STABLE_DIRECTION  |             The `unit_vectors` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `unit_vectors` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        real(real64), intent(in), optional :: min_angular_dispersion
            !! The minimum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MIN_ANGULAR_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_angular_dispersion
            !! The maximum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MAX_ANGULAR_DISPERSION_DEFAULT`

        integer(int32), dimension(:), allocatable :: tmp_member_counts

        ! Input validation
        call set_ok(ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(unit_vectors, size(unit_vectors, kind=int32), ierr, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_member_counts(n_families))

        call tox_compute_family_direction_helper(unit_vectors, gene_to_fam, &
                                                 family_directions, angular_dispersions, &
                          n_samples, n_genes, n_families, status, tmp_member_counts, min_angular_dispersion, max_angular_dispersion)
    end subroutine tox_compute_family_direction_alloc

    !> Compute angles between each gene and its family reference direction
    pure subroutine tox_compute_angles_to_direction(unit_vectors, family_directions, &
                                                    gene_to_fam, angles, &
                                                    n_samples, n_genes, n_families, ierr)
        integer(int32), intent(in) :: n_samples
            !! Number of samples (dimension of vectors)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_samples, n_genes), intent(in) :: unit_vectors
            !! Unit vectors for each gene
        real(real64), dimension(n_samples, n_families), intent(in) :: family_directions
            !! Family reference directions
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: angles
            !! Output angles in radians [0, π], `CM_ANGLES_SENTINEL` for invalid
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success

        call set_ok(ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
      call validate_all_in_range_real(family_directions, size(family_directions, kind=int32), ierr, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        call tox_compute_angles_to_direction_helper(unit_vectors, family_directions, &
                                                    gene_to_fam, angles, &
                                                    n_samples, n_genes, n_families)

    end subroutine tox_compute_angles_to_direction

    !> (no input validation) Compute angles between each gene and its family reference direction
    pure subroutine tox_compute_angles_to_direction_helper(unit_vectors, family_directions, &
                                                           gene_to_fam, angles, &
                                                           n_samples, n_genes, n_families)
        integer(int32), intent(in) :: n_samples
            !! Number of samples (dimension of vectors)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_samples, n_genes), intent(in) :: unit_vectors
            !! Unit vectors for each gene
        real(real64), dimension(n_samples, n_families), intent(in) :: family_directions
            !! Family reference directions
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: angles
            !! Output angles in radians [0, π], `CM_ANGLES_SENTINEL` for invalid (trying to compute angle for zero vectors)

        integer(int32) :: tmp_ierr

        integer(int32) :: i_gene, fam_idx

        do concurrent (i_gene = 1:n_genes) local(tmp_ierr, fam_idx) shared(gene_to_fam, n_families, angles, unit_vectors, family_directions, n_samples)
            call set_ok(tmp_ierr)

            fam_idx = gene_to_fam(i_gene)

            if (fam_idx == M_GENE_TO_FAM_SENTINEL) then
                angles(i_gene) = CM_ANGLES_SENTINEL
                cycle
            end if

            ! Use the reliable angle_between function
            call angle_between(unit_vectors(:, i_gene), &
                               family_directions(:, fam_idx), &
                               n_samples, angles(i_gene), tmp_ierr)

            ! If angle_between returns error, mark as invalid
            if (is_err(tmp_ierr)) then
                angles(i_gene) = CM_ANGLES_SENTINEL
            end if
        end do
    end subroutine tox_compute_angles_to_direction_helper

    !> Compute scaled angles (z-scores) by using angular dispersion (angular std dev)
    !| \[z\_scores(i\_gene) = \frac{angles(i\_gene)}{angular\_dispersions(gene\_to\_fam(i\_gene))}\]
    pure subroutine tox_z_scores_by_dispersion(angles, angular_dispersions, &
                                               gene_to_fam, z_scores, &
                                               n_genes, n_families, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: angles
            !! Angles in radians [0, π], `CM_ANGLES_SENTINEL` for invalid
        real(real64), dimension(n_families), intent(in) :: angular_dispersions
            !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\) from [[tox_get_outliers_by_angle(module):tox_compute_family_direction(subroutine)]], in range `[min_angular_dispersion, max_angular_dispersion]` with `CM_ANGULAR_DISPERSION_SENTINEL` for invalid families
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(angles, size(angles, kind=int32), ierr, min=0.0_real64, max=PI, sentinel=CM_ANGLES_SENTINEL)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        call validate_all_in_range_real(angular_dispersions, size(angular_dispersions, kind=int32), ierr, min=0.0_real64, sentinel=CM_ANGULAR_DISPERSION_SENTINEL)
        if (is_err(ierr)) return

        call tox_z_scores_by_dispersion_helper(angles, angular_dispersions, &
                                               gene_to_fam, z_scores, &
                                               n_genes, n_families)

    end subroutine tox_z_scores_by_dispersion

    !> (no input validation) Compute scaled angles (z-scores) by using angular dispersion (angular std dev)
    !| \[z\_scores(i\_gene) = \frac{angles(i\_gene)}{angular\_dispersions(gene\_to\_fam(i\_gene))}\]
    pure subroutine tox_z_scores_by_dispersion_helper(angles, angular_dispersions, &
                                                      gene_to_fam, z_scores, &
                                                      n_genes, n_families)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: angles
            !! Angles in radians [0, π], `CM_ANGLES_SENTINEL` for invalid
        real(real64), dimension(n_families), intent(in) :: angular_dispersions
            !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\) from [[tox_get_outliers_by_angle(module):tox_compute_family_direction(subroutine)]], in range `[min_angular_dispersion, max_angular_dispersion]` with `CM_ANGULAR_DISPERSION_SENTINEL` for invalid families
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid

        integer(int32) :: i_gene, fam_idx

        do concurrent(i_gene=1:n_genes) local(fam_idx) shared(gene_to_fam, n_families, z_scores, angles, angular_dispersions)
            fam_idx = gene_to_fam(i_gene)
            if (fam_idx < 1 .or. fam_idx > n_families) then
                z_scores(i_gene) = CM_Z_SCORES_SENTINEL
                cycle
            end if

            if (angles(i_gene) == CM_ANGLES_SENTINEL .or. angular_dispersions(fam_idx) == CM_ANGULAR_DISPERSION_SENTINEL) then
                z_scores(i_gene) = CM_Z_SCORES_SENTINEL
                cycle
            end if

            if (is_close(angular_dispersions(fam_idx), 0.0_real64)) then
                z_scores(i_gene) = CM_Z_SCORES_SENTINEL
            else
                z_scores(i_gene) = angles(i_gene) / angular_dispersions(fam_idx)
            end if
        end do
    end subroutine tox_z_scores_by_dispersion_helper

    !> Identify directional outliers based on scaled angles using percentile threshold
    !| \[is\_outlier(i\_gene) = z\_scores(i\_gene) \ge threshold\]
    pure subroutine tox_angle_outliers(z_scores, threshold, is_outlier, n_genes, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        real(real64), intent(in) :: threshold
            !! `z_scores` exceeding this value will be identified as outliers
        real(real64), dimension(n_genes), intent(in) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output boolean array indicating outlier status per gene
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success

        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_in_range_real(threshold, ierr)
        call validate_all_in_range_real(z_scores, size(z_scores, kind=int32), ierr, min=0.0_real64, sentinel=CM_Z_SCORES_SENTINEL)
        if (is_err(ierr)) return

        call tox_angle_outliers_helper(z_scores, threshold, is_outlier, n_genes)
    end subroutine tox_angle_outliers

    !> (no input validation) Identify directional outliers based on scaled angles using percentile threshold
    !| \[is\_outlier(i\_gene) = z\_scores(i\_gene) \ge threshold\]
    pure subroutine tox_angle_outliers_helper(z_scores, threshold, is_outlier, n_genes)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        real(real64), intent(in) :: threshold
            !! `z_scores>=threshold` will be identified as outliers
        real(real64), dimension(n_genes), intent(in) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output boolean array indicating outlier status per gene

        integer(int32) :: i_gene

        ! Identify outliers
        do concurrent(i_gene=1:n_genes) shared(z_scores, is_outlier, threshold)
            is_outlier(i_gene) = z_scores(i_gene) >= threshold
        end do
    end subroutine tox_angle_outliers_helper

    !> Identify directional outliers based on scaled angles using percentile threshold
    !| \[is\_outlier(i\_gene) = z\_scores(i\_gene) \ge threshold\]
    pure subroutine tox_angle_outliers_alloc(z_scores, percentile, threshold, is_outlier, n_genes, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        real(real64), dimension(n_genes), intent(in) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
        real(real64), intent(in) :: percentile
            !! Percentile threshold (`0.0 <= p <= 1.0`)
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output boolean array indicating outlier status per gene
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success
        real(real64), intent(out) :: threshold
            !! Threshold an outlier's angle exceeds, `CM_Z_SCORES_SENTINEL` if no threshold could be computed (fewer than two valid z-scores or invalid input)

        integer(int32) :: i_gene, count_valid, idx
        integer(int32), dimension(:), allocatable :: perm

        call set_ok(ierr)
        ! Defined default so `threshold` is never returned uninitialized on an early return
        threshold = CM_Z_SCORES_SENTINEL

        call validate_dimension_size(n_genes, ierr)
        call validate_in_range_real(percentile, ierr, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(z_scores, size(z_scores, kind=int32), ierr, min=0.0_real64, sentinel=CM_Z_SCORES_SENTINEL)
        if (is_err(ierr)) return

        ! Count valid scaled angles
        count_valid = 0
        do i_gene = 1, n_genes
            if (z_scores(i_gene) /= CM_Z_SCORES_SENTINEL) count_valid = count_valid + 1
        end do

        if (count_valid < 2) then
            threshold = maxval(z_scores)
            is_outlier = .false.
            return
        end if

        M_ALLOCATE(perm(count_valid))

        ! Collect valid scores
        idx = 1
        do i_gene = 1, n_genes
            if (z_scores(i_gene) >= 0.0_real64) then
                perm(idx) = i_gene
                idx = idx + 1
            end if
        end do

        ! Sort and compute percentile
        call sort_real_heapsort(z_scores, perm)

        call calc_percentile_helper(z_scores, perm, percentile, threshold)
        call tox_angle_outliers_helper(z_scores, threshold, is_outlier, n_genes)
    end subroutine tox_angle_outliers_alloc

    !> Complete pipeline for angle-based outlier detection using directional statistics
    !| Full pipeline:
    !|
    !| 1. [[tox_normalization(module):tox_normalize_vectors_unit_length(subroutine)]]
    !| 2. [[tox_get_outliers_by_angle(module):tox_compute_family_direction(subroutine)]]
    !| 3. [[tox_get_outliers_by_angle(module):tox_compute_angles_to_direction(subroutine)]]
    !| 4. [[tox_get_outliers_by_angle(module):tox_z_scores_by_dispersion(subroutine)]]
    !| 5. [[tox_get_outliers_by_angle(module):tox_angle_outliers_alloc(subroutine)]]
    !|
    pure subroutine tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                                       percentile, threshold, z_scores, &
                                                       n_samples, n_genes, n_families, &
                                                       is_outlier, ierr, status, min_angular_dispersion, max_angular_dispersion)
        integer(int32), intent(in) :: n_samples
            !! Number of samples (dimension of expression vectors)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_samples, n_genes), intent(in) :: expression_vectors
            !! Gene expression vectors, shape (n_samples, n_genes)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
        real(real64), intent(in) :: percentile
            !! Percentile for outlier detection (0.0-100.0)
        real(real64), intent(out) :: threshold
            !! Threshold value calculated with the percentile
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output outlier status for each gene
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                               |
            !! |-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                   |
            !! | STAT_NO_STABLE_DIRECTION  |             The `unit_vectors` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `unit_vectors` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        real(real64), intent(in), optional :: min_angular_dispersion
            !! The minimum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MIN_ANGULAR_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_angular_dispersion
            !! The maximum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. default: `CM_MAX_ANGULAR_DISPERSION_DEFAULT`

        ! Local variables
        real(real64), dimension(:, :), allocatable :: unit_vectors
        real(real64), dimension(:, :), allocatable :: family_directions
        real(real64), dimension(:), allocatable :: angular_dispersions
        real(real64), dimension(:), allocatable :: angles
        integer(int32), dimension(:), allocatable :: tmp_member_counts

        call set_ok(ierr)
        call set_ok(status)
        is_outlier = .false.

        ! Input validation
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_in_range_real(percentile, ierr, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(expression_vectors, size(expression_vectors, kind=int32), ierr)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        ! Allocate work arrays
        M_ALLOCATE(unit_vectors(n_samples, n_genes))
        M_ALLOCATE(family_directions(n_samples, n_families))
        M_ALLOCATE(angular_dispersions(n_families))
        M_ALLOCATE(angles(n_genes))
        M_ALLOCATE(tmp_member_counts(n_families))

        ! Step 1: Normalize vectors to unit length
        call tox_normalize_vectors_unit_length_helper(expression_vectors, unit_vectors, &
                                                      n_samples, n_genes, ierr)
        if (is_err(ierr)) return

        ! Step 2: Compute family reference directions and angular dispersions
        call tox_compute_family_direction_helper(unit_vectors, gene_to_fam, &
                                                 family_directions, angular_dispersions, &
                          n_samples, n_genes, n_families, status, tmp_member_counts, min_angular_dispersion, max_angular_dispersion)

        ! Step 3: Compute angles between genes and family directions
        call tox_compute_angles_to_direction_helper(unit_vectors, family_directions, &
                                                    gene_to_fam, angles, &
                                                    n_samples, n_genes, n_families)

        ! Step 4: Scale angles by angular dispersion
        call tox_z_scores_by_dispersion_helper(angles, angular_dispersions, &
                                               gene_to_fam, z_scores, &
                                               n_genes, n_families)

        call tox_angle_outliers_alloc(z_scores, percentile, threshold, is_outlier, n_genes, ierr)
    end subroutine tox_detect_angle_outliers_pipeline

end module tox_get_outliers_by_angle

!> C wrapper for tox_detect_angle_outliers_pipeline
pure subroutine tox_detect_angle_outliers_pipeline_c(n_samples, n_genes, n_families, &
                                                     expression_vectors, gene_to_fam, &
                                                     percentile, threshold, is_outlier, z_scores, &
                                                     ierr, status, min_angular_dispersion, max_angular_dispersion) &
    bind(C, name="tox_detect_angle_outliers_pipeline_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle, only: tox_detect_angle_outliers_pipeline
    use tox_conversions, only: logical_as_c_int
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_samples
        !! Number of samples
    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_samples, n_genes), intent(in), target :: expression_vectors
        !! Gene expression vectors
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping
    real(c_double), intent(in), target :: percentile
        !! Percentile for outlier detection
    real(c_double), intent(out), target :: threshold
        !! Threshold value calculated with the percentile
    integer(c_int), dimension(n_genes), intent(out), target :: is_outlier
        !! Output outlier status per gene
    real(c_double), dimension(n_genes), intent(out), target :: z_scores
        !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), dimension(n_families), intent(out), target :: status
        !! Informs about special families
        !!
        !! |            Code             |                                                               Meaning                                                               |
        !! |-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
        !! |          ERR_OK           |                                                   Normal family, no special case.                                                   |
        !! | STAT_NO_STABLE_DIRECTION  |             The `unit_vectors` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
        !! | STAT_NO_ANGULAR_VARIATION | The `unit_vectors` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
        !!
    real(c_double), intent(in), target :: min_angular_dispersion
        !! The minimum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. recommended: `CM_MIN_ANGULAR_DISPERSION_DEFAULT`
    real(c_double), intent(in), target :: max_angular_dispersion
        !! The maximum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. recommended: `CM_MAX_ANGULAR_DISPERSION_DEFAULT`

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(expression_vectors)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(is_outlier)
    M_CHECK_NON_NULL(z_scores)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(status)
    M_CHECK_NON_NULL(min_angular_dispersion)
    M_CHECK_NON_NULL(max_angular_dispersion)

    M_ALLOCATE(is_outlier_f(n_genes))

    call tox_detect_angle_outliers_pipeline(expression_vectors, gene_to_fam, &
                                            percentile, threshold, z_scores, &
                                            n_samples, n_genes, n_families, &
                                            is_outlier_f, ierr, status, min_angular_dispersion, max_angular_dispersion)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine tox_detect_angle_outliers_pipeline_c

!> C wrapper for tox_compute_family_direction
pure subroutine tox_compute_family_direction_expert_c(unit_vectors, gene_to_fam, &
                                                      family_directions, angular_dispersions, &
                                                      n_samples, n_genes, n_families, &
                                                  ierr, status, tmp_member_counts, min_angular_dispersion, max_angular_dispersion) &
    bind(C, name="tox_compute_family_direction_expert_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle, only: tox_compute_family_direction
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_samples
        !! Number of samples
    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_samples, n_genes), intent(in), target :: unit_vectors
        !! Unit vectors for each gene
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    integer(c_int), dimension(n_families), intent(out), target :: tmp_member_counts
        !! Number of members per family, determined by `gene_to_fam`
    real(c_double), dimension(n_samples, n_families), intent(out), target :: family_directions
        !! Family reference directions
    real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
        !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\)
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), dimension(n_families), intent(out), target :: status
        !! Informs about special families
        !!
        !! |            Code             |                                                               Meaning                                                               |
        !! |-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
        !! |          ERR_OK           |                                                   Normal family, no special case.                                                   |
        !! | STAT_NO_STABLE_DIRECTION  |             The `unit_vectors` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
        !! | STAT_NO_ANGULAR_VARIATION | The `unit_vectors` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
        !!
    real(c_double), intent(in), target :: min_angular_dispersion
        !! The minimum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. recommended: `CM_MIN_ANGULAR_DISPERSION_DEFAULT`
    real(c_double), intent(in), target :: max_angular_dispersion
        !! The maximum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. recommended: `CM_MAX_ANGULAR_DISPERSION_DEFAULT`

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(unit_vectors)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(tmp_member_counts)
    M_CHECK_NON_NULL(family_directions)
    M_CHECK_NON_NULL(angular_dispersions)
    M_CHECK_NON_NULL(status)
    M_CHECK_NON_NULL(min_angular_dispersion)
    M_CHECK_NON_NULL(max_angular_dispersion)

    call tox_compute_family_direction(unit_vectors, gene_to_fam, &
                                      family_directions, angular_dispersions, &
                                      n_samples, n_genes, n_families, &
                                      ierr, status, tmp_member_counts, min_angular_dispersion, max_angular_dispersion)
end subroutine tox_compute_family_direction_expert_c

!> C wrapper for tox_compute_family_direction_alloc
pure subroutine tox_compute_family_direction_c(unit_vectors, n_samples, n_genes, n_families, &
                                               gene_to_fam, family_directions, angular_dispersions, &
                                               ierr, status, min_angular_dispersion, max_angular_dispersion) &
    bind(C, name="tox_compute_family_direction_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle, only: tox_compute_family_direction_alloc
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_samples
        !! Number of samples
    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_samples, n_genes), intent(in), target :: unit_vectors
        !! Unit vectors for each gene
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    real(c_double), dimension(n_samples, n_families), intent(out), target :: family_directions
        !! Family reference directions
    real(c_double), dimension(n_families), intent(out), target :: angular_dispersions
        !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\)
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), dimension(n_families), intent(out), target :: status
        !! Informs about special families
        !!
        !! |            Code             |                                                               Meaning                                                               |
        !! |-----------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
        !! |          ERR_OK           |                                                   Normal family, no special case.                                                   |
        !! | STAT_NO_STABLE_DIRECTION  |             The `unit_vectors` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
        !! | STAT_NO_ANGULAR_VARIATION | The `unit_vectors` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
        !!
    real(c_double), intent(in), target :: min_angular_dispersion
        !! The minimum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. recommended: `CM_MIN_ANGULAR_DISPERSION_DEFAULT`
    real(c_double), intent(in), target :: max_angular_dispersion
        !! The maximum angular dispersion, otherwise `CM_ANGULAR_DISPERSION_SENTINEL` is used. recommended: `CM_MAX_ANGULAR_DISPERSION_DEFAULT`

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(unit_vectors)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(family_directions)
    M_CHECK_NON_NULL(angular_dispersions)
    M_CHECK_NON_NULL(status)
    M_CHECK_NON_NULL(min_angular_dispersion)
    M_CHECK_NON_NULL(max_angular_dispersion)

    call tox_compute_family_direction_alloc(unit_vectors, gene_to_fam, &
                                            family_directions, angular_dispersions, &
                                            n_samples, n_genes, n_families, &
                                            ierr, status, min_angular_dispersion, max_angular_dispersion)
end subroutine tox_compute_family_direction_c

!> C wrapper for tox_compute_angles_to_direction
pure subroutine tox_compute_angles_to_direction_c(unit_vectors, n_samples, n_genes, n_families, &
                                                  gene_to_fam, family_directions, angles, ierr) &
    bind(C, name="tox_compute_angles_to_direction_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle, only: tox_compute_angles_to_direction
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_samples
        !! Number of samples
    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_samples, n_genes), intent(in), target :: unit_vectors
        !! Unit vectors for each gene
    real(c_double), dimension(n_samples, n_families), intent(in), target :: family_directions
        !! Family reference directions
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    real(c_double), dimension(n_genes), intent(out), target :: angles
        !! Output angles in radians [0, π], `CM_ANGLES_SENTINEL` for invalid (trying to compute angle for zero vectors)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(unit_vectors)
    M_CHECK_NON_NULL(family_directions)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(angles)

    call tox_compute_angles_to_direction(unit_vectors, family_directions, &
                                         gene_to_fam, angles, &
                                         n_samples, n_genes, n_families, ierr)
end subroutine tox_compute_angles_to_direction_c

!> C wrapper for tox_z_scores_by_dispersion
pure subroutine tox_z_scores_by_dispersion_c(angles, n_genes, n_families, &
                                             gene_to_fam, angular_dispersions, &
                                             z_scores, ierr) &
    bind(C, name="tox_z_scores_by_dispersion_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle, only: tox_z_scores_by_dispersion
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_genes), intent(in), target :: angles
        !! Angles in radians [0, π], `CM_ANGLES_SENTINEL` for invalid
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    real(c_double), dimension(n_families), intent(in), target :: angular_dispersions
        !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\) from [[tox_get_outliers_by_angle(module):tox_compute_family_direction(subroutine)]], in range `[min_angular_dispersion, max_angular_dispersion]` with `CM_ANGULAR_DISPERSION_SENTINEL` for invalid families
    real(c_double), dimension(n_genes), intent(out), target :: z_scores
        !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(angles)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(angular_dispersions)
    M_CHECK_NON_NULL(z_scores)

    call tox_z_scores_by_dispersion(angles, angular_dispersions, &
                                    gene_to_fam, z_scores, &
                                    n_genes, n_families, ierr)
end subroutine tox_z_scores_by_dispersion_c

!> C wrapper for tox_angle_outliers_alloc
pure subroutine tox_angle_outliers_c(z_scores, n_genes, percentile, &
                                     threshold, is_outlier, ierr) &
    bind(C, name="tox_angle_outliers_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle, only: tox_angle_outliers_alloc
    use tox_conversions, only: logical_as_c_int
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    real(c_double), dimension(n_genes), intent(in), target :: z_scores
        !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
    real(c_double), intent(in), target :: percentile
        !! Percentile threshold (`0.0 < p < 100.0`)
    real(c_double), intent(out), target :: threshold
        !! Computed threshold value for outlier detection
    integer(c_int), dimension(n_genes), intent(out), target :: is_outlier
        !! Output boolean array indicating outlier status per gene
    integer(c_int), intent(out), target :: ierr
        !! Error code

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(z_scores)
    M_CHECK_NON_NULL(is_outlier)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(threshold)

    M_ALLOCATE(is_outlier_f(n_genes))

    call tox_angle_outliers_alloc(z_scores, percentile, threshold, &
                                  is_outlier_f, n_genes, ierr)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine tox_angle_outliers_c

pure subroutine tox_angle_outliers_expert_c(z_scores, threshold, is_outlier, n_genes, ierr) &
    bind(C, name="tox_angle_outliers_expert_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle, only: tox_angle_outliers
    use tox_conversions, only: logical_as_c_int
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    real(c_double), intent(in), target :: threshold
        !! `z_scores` exceeding this value will be identified as outliers
    real(c_double), dimension(n_genes), intent(in), target :: z_scores
        !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
    integer(c_int), dimension(n_genes), intent(out), target :: is_outlier
        !! Output boolean array indicating outlier status per gene
    integer(c_int), intent(out), target :: ierr
        !! Error code

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(z_scores)
    M_CHECK_NON_NULL(threshold)
    M_CHECK_NON_NULL(is_outlier)

    M_ALLOCATE(is_outlier_f(n_genes))

    call tox_angle_outliers(z_scores, threshold, is_outlier_f, n_genes, ierr)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine tox_angle_outliers_expert_c
