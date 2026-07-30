#include <src/macros.h>

module tox_get_outliers_by_angle_rap
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use f42_utils, only: PI, is_close, wrap_angle
    use tox_errors, only: STAT_NO_ANGULAR_VARIATION, STAT_NO_STABLE_DIRECTION, set_ok, validate_dimension_size, ERR_ALLOC_FAIL, &
                          validate_all_in_range_real, validate_all_in_range_int, &
                          validate_in_range_real, is_err, set_err
    use tox_get_outliers_by_angle, only: tox_angle_outliers, tox_angle_outliers_alloc
    implicit none

    private
    public :: tox_compute_family_direction_rap, &
              tox_compute_angular_deviations_rap, &
              tox_z_scores_by_dispersion_rap, &
              tox_angle_outliers_rap, &
              tox_detect_angle_outliers_pipeline_rap

    public :: tox_compute_family_direction_rap_alloc, &
              tox_angle_outliers_rap_alloc

#define CM_MIN_FAMILY_DISPERSION_DEFAULT 1.0e-02_real64
#define CM_MAX_FAMILY_DISPERSION_DEFAULT sqrt(-2.0_real64 * log(0.5_real64))

#define CM_RAP_ANGLES_SENTINEL -1.0_real64
#define CM_FAMILY_DISPERSION_SENTINEL -1.0_real64
#define CM_Z_SCORES_SENTINEL -1.0_real64
#define CM_ANGULAR_DEVIATIONS_SENTINEL -1.0_real64

contains

    !> Compute circular mean and dispersion for each family
    !| Some angles might result in sentinel values. In such case the respective family is being marked by the `status(i_family)` set to the respective status code.
    pure subroutine tox_compute_family_direction_rap(rap_angles, gene_to_fam, &
                                                     family_mean_angles, family_dispersions, &
                                 n_genes, n_families, ierr, status, tmp_member_counts, min_family_dispersion, max_family_dispersion)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: rap_angles
            !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_families), intent(out) :: family_mean_angles
            !! Family mean angles in radians, `0<=x<=PI`. The circular mean is kept for every family with >=2 members and a non-zero resultant length (even when the family is flagged via the dispersion bounds); it is `0` only for families with fewer than 2 members or a zero resultant length
        real(real64), dimension(n_families), intent(out) :: family_dispersions
            !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_FAMILY_DISPERSION_SENTINEL` if there is no stabe direction or too less angular variation (`status`)
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                             |
            !! |-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                 |
            !! | STAT_NO_STABLE_DIRECTION  |             The `rap_angles` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `rap_angles` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        integer(int32), dimension(n_families), intent(out) :: tmp_member_counts
            !! Permutation of gene_to_fam sorted by family index
        real(real64), intent(in), optional :: min_family_dispersion
            !! The minimum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MIN_FAMILY_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_family_dispersion
            !! The maximum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MAX_FAMILY_DISPERSION_DEFAULT`

        ! Input validation
        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(rap_angles, size(rap_angles, kind=int32), ierr, min=0.0_real64, max=PI, sentinel=CM_RAP_ANGLES_SENTINEL)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        call tox_compute_family_direction_rap_helper(rap_angles, gene_to_fam, &
                                                     family_mean_angles, family_dispersions, &
                                       n_genes, n_families, status, tmp_member_counts, min_family_dispersion, max_family_dispersion)
    end subroutine tox_compute_family_direction_rap

    !> (no input validation) Compute circular mean and dispersion for each family
    !| Some angles might result in sentinel values. In such case the respective family is being marked by the `status(i_family)` set to the respective status code.
    pure subroutine tox_compute_family_direction_rap_helper(rap_angles, gene_to_fam, &
                                                            family_mean_angles, family_dispersions, &
                                       n_genes, n_families, status, tmp_member_counts, min_family_dispersion, max_family_dispersion)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: rap_angles
            !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_families), intent(out), target :: family_mean_angles
            !! Family mean angles in radians, `0<=x<=PI`. The circular mean is kept for every family with >=2 members and a non-zero resultant length (even when the family is flagged via the dispersion bounds); it is `0` only for families with fewer than 2 members or a zero resultant length
        real(real64), dimension(n_families), intent(out), target :: family_dispersions
            !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_FAMILY_DISPERSION_SENTINEL` if there is no stabe direction or too less angular variation (`status`)
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                             |
            !! |-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                 |
            !! | STAT_NO_STABLE_DIRECTION  |             The `rap_angles` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `rap_angles` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        integer(int32), dimension(n_families), intent(out) :: tmp_member_counts
            !! Number of members per family, determined by `gene_to_fam`
        real(real64), intent(in), optional :: min_family_dispersion
            !! The minimum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MIN_FAMILY_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_family_dispersion
            !! The maximum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MAX_FAMILY_DISPERSION_DEFAULT`

        integer(int32) :: i_gene, i_family, family_idx
        real(real64) :: mean_length, angular_dispersion, actual_min_family_dispersion, actual_max_family_dispersion
        real(real64), dimension(:), pointer :: sum_sin_ptr, sum_cos_ptr

        M_DEFAULT_VAL(min_family_dispersion, actual_min_family_dispersion, CM_MIN_FAMILY_DISPERSION_DEFAULT)
        M_DEFAULT_VAL(max_family_dispersion, actual_max_family_dispersion, CM_MAX_FAMILY_DISPERSION_DEFAULT)

        call set_ok(status)

        family_mean_angles = 0.0_real64
        family_dispersions = 0.0_real64
        tmp_member_counts = 0_int32

        ! Reusing arrays for precalculations
        sum_sin_ptr => family_mean_angles
        sum_cos_ptr => family_dispersions
        do i_gene = 1, n_genes
            family_idx = gene_to_fam(i_gene)
            if (family_idx /= M_GENE_TO_FAM_SENTINEL .and. rap_angles(i_gene) /= CM_RAP_ANGLES_SENTINEL) then
                tmp_member_counts(family_idx) = tmp_member_counts(family_idx) + 1
                sum_sin_ptr(family_idx) = sum_sin_ptr(family_idx) + sin(rap_angles(i_gene))
                sum_cos_ptr(family_idx) = sum_cos_ptr(family_idx) + cos(rap_angles(i_gene))
            end if
        end do

        ! Compute circular statistics per family
        do concurrent(i_family=1:n_families) &
            local(mean_length, angular_dispersion) &
            shared(tmp_member_counts, family_dispersions, status, sum_cos_ptr, sum_sin_ptr, family_mean_angles)

            ! Need at least 2 genes for meaningful statistics. With fewer members the
            ! circular mean is undefined, so reset it to 0 instead of leaving the leftover
            ! `sum_sin` accumulator (a sine value, not an angle) in `family_mean_angles`.
            if (tmp_member_counts(i_family) < 2) then
                family_mean_angles(i_family) = 0.0_real64
                family_dispersions(i_family) = CM_FAMILY_DISPERSION_SENTINEL
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                cycle
            end if

            ! Compute resultant length R_F
            mean_length = sqrt(sum_cos_ptr(i_family)**2 + sum_sin_ptr(i_family)**2) / real(tmp_member_counts(i_family), real64)

            ! Edge case handling, will always exceed `max_family_dispersion`, thus no stable direction.
            ! A zero resultant length means there is no mean direction, so reset it to 0.
            if (mean_length == 0.0_real64) then
                family_mean_angles(i_family) = 0.0_real64
                family_dispersions(i_family) = CM_FAMILY_DISPERSION_SENTINEL
                call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
            else
                ! IMPORTANT: `sum_sin_ptr` and `sum_cos_ptr` point to `family_mean_angles` and `family_dispersions` to reuse them as work space.
                ! The circular mean is well defined here (>=2 members, non-zero resultant), so it is kept even when the
                ! family is subsequently flagged via the dispersion bounds below.
                family_mean_angles(i_family) = atan2(sum_sin_ptr(i_family), sum_cos_ptr(i_family))

                angular_dispersion = sqrt(-2.0_real64 * log(mean_length))

                ! Check for minimal angular variation
                if (angular_dispersion < actual_min_family_dispersion) then
                    family_dispersions(i_family) = CM_FAMILY_DISPERSION_SENTINEL
                    call set_err(status(i_family), STAT_NO_ANGULAR_VARIATION)

                    ! Check for maximum angular variation
                else if (angular_dispersion > actual_max_family_dispersion) then
                    family_dispersions(i_family) = CM_FAMILY_DISPERSION_SENTINEL
                    call set_err(status(i_family), STAT_NO_STABLE_DIRECTION)
                else
                    family_dispersions(i_family) = angular_dispersion
                end if
            end if
        end do

    end subroutine tox_compute_family_direction_rap_helper

    !> Compute circular mean and dispersion for each family
    !| Some angles might result in sentinel values. In such case the respective family is being marked by the `status(i_family)` set to the respective status code.
    pure subroutine tox_compute_family_direction_rap_alloc(rap_angles, gene_to_fam, &
                                                           family_mean_angles, family_dispersions, &
                                                    n_genes, n_families, ierr, status, min_family_dispersion, max_family_dispersion)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: rap_angles
            !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_families), intent(out) :: family_mean_angles
            !! Family mean angles in radians, `0<=x<=PI`. The circular mean is kept for every family with >=2 members and a non-zero resultant length (even when the family is flagged via the dispersion bounds); it is `0` only for families with fewer than 2 members or a zero resultant length
        real(real64), dimension(n_families), intent(out) :: family_dispersions
            !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_FAMILY_DISPERSION_SENTINEL` if there is no stabe direction or too less angular variation (`status`)
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                             |
            !! |-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                 |
            !! | STAT_NO_STABLE_DIRECTION  |             The `rap_angles` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `rap_angles` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        real(real64), intent(in), optional :: min_family_dispersion
            !! The minimum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MIN_FAMILY_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_family_dispersion
            !! The maximum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MAX_FAMILY_DISPERSION_DEFAULT`

        integer(int32), allocatable :: tmp_member_counts(:)

        ! Input validation
        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(rap_angles, size(rap_angles, kind=int32), ierr, min=0.0_real64, max=PI, sentinel=CM_RAP_ANGLES_SENTINEL)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_member_counts(n_families))

        ! Call helper pure subroutine
        call tox_compute_family_direction_rap_helper(rap_angles, gene_to_fam, &
                                                     family_mean_angles, family_dispersions, &
                                       n_genes, n_families, status, tmp_member_counts, min_family_dispersion, max_family_dispersion)

    end subroutine tox_compute_family_direction_rap_alloc

    !> Compute the absoute angular deviation of a `rap_angle` from its corresponding `family_mean_angle`
    pure subroutine tox_compute_angular_deviations_rap(rap_angles, family_mean_angles, &
                                                       gene_to_fam, angular_deviations, &
                                                       n_genes, n_families, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: rap_angles
            !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
        real(real64), dimension(n_families), intent(in) :: family_mean_angles
            !! Family mean angles in radians, `0<=x<=PI`
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: angular_deviations
            !! Absolute wrapped angular deviations [0, π], `CM_ANGULAR_DEVIATIONS_SENTINEL` for invalid
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(rap_angles, size(rap_angles, kind=int32), ierr, min=0.0_real64, max=PI, sentinel=CM_RAP_ANGLES_SENTINEL)
        call validate_all_in_range_real(family_mean_angles, size(family_mean_angles, kind=int32), ierr, min=0.0_real64, max=PI)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        call tox_compute_angular_deviations_rap_helper(rap_angles, family_mean_angles, &
                                                       gene_to_fam, angular_deviations, &
                                                       n_genes, n_families)
    end subroutine tox_compute_angular_deviations_rap

    !> (no input validation) Compute the absoute angular deviation of a `rap_angle` from its corresponding `family_mean_angle`
    pure subroutine tox_compute_angular_deviations_rap_helper(rap_angles, family_mean_angles, &
                                                              gene_to_fam, angular_deviations, &
                                                              n_genes, n_families)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: rap_angles
            !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
        real(real64), dimension(n_families), intent(in) :: family_mean_angles
            !! Family mean angles in radians, `0<=x<=PI`
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: angular_deviations
            !! Absolute wrapped angular deviations [0, π], `CM_ANGULAR_DEVIATIONS_SENTINEL` for invalid

        integer(int32) :: i_gene, fam_idx
        real(real64) :: delta

        do concurrent (i_gene = 1:n_genes) local(fam_idx, delta) shared(gene_to_fam, n_families, angular_deviations, rap_angles, family_mean_angles)
            fam_idx = gene_to_fam(i_gene)

            ! Check for invalid/unassigned gene
            if (fam_idx < 1 .or. fam_idx > n_families) then
                angular_deviations(i_gene) = CM_ANGULAR_DEVIATIONS_SENTINEL
                cycle
            end if

            ! Compute wrapped angular difference
            delta = rap_angles(i_gene) - family_mean_angles(fam_idx)
            angular_deviations(i_gene) = abs(wrap_angle(delta))
        end do

    end subroutine tox_compute_angular_deviations_rap_helper

    !> Compute scaled angular deviations (z-scores) by using angular dispersion (angular std dev)
    !| \[z\_scores(i\_gene) = \frac{angular\_deviations(i\_gene)}{family\_dispersions(gene\_to\_fam(i\_gene))}\]
    pure subroutine tox_z_scores_by_dispersion_rap(angular_deviations, family_dispersions, &
                                                   gene_to_fam, z_scores, &
                                                   n_genes, n_families, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: angular_deviations
            !! Angular deviations in radians `0<=x<=PI` and `CM_ANGULAR_DEVIATIONS_SENTINEL` for invalid/unassigned genes
        real(real64), dimension(n_families), intent(in) :: family_dispersions
            !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\) from [[tox_get_outliers_by_angle_rap(module):tox_compute_family_direction_rap(subroutine)]], in range `[min_family_dispersion, max_family_dispersion]` with `CM_FAMILY_DISPERSION_SENTINEL` for invalid families
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(angular_deviations, size(angular_deviations, kind=int32), ierr, min=0.0_real64, max=PI, sentinel=CM_ANGULAR_DEVIATIONS_SENTINEL)
        call validate_all_in_range_real(family_dispersions, size(family_dispersions, kind=int32), ierr, min=0.0_real64, sentinel=CM_FAMILY_DISPERSION_SENTINEL)
        call validate_all_in_range_int(gene_to_fam, size(gene_to_fam, kind=int32), ierr, min=1_int32, max=n_families, sentinel=M_GENE_TO_FAM_SENTINEL)
        if (is_err(ierr)) return

        call tox_z_scores_by_dispersion_rap_helper(angular_deviations, family_dispersions, &
                                                   gene_to_fam, z_scores, &
                                                   n_genes, n_families)
    end subroutine tox_z_scores_by_dispersion_rap

    !> (no input validation) Compute scaled angular deviations (z-scores) by using angular dispersion (angular std dev)
    !| \[z\_scores(i\_gene) = \frac{angular\_deviations(i\_gene)}{family\_dispersions(gene\_to\_fam(i\_gene))}\]
    pure subroutine tox_z_scores_by_dispersion_rap_helper(angular_deviations, family_dispersions, &
                                                          gene_to_fam, z_scores, &
                                                          n_genes, n_families)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        integer(int32), intent(in) :: n_families
            !! Number of gene families
        real(real64), dimension(n_genes), intent(in) :: angular_deviations
            !! Angular deviations in radians `0<=x<=PI` and `CM_ANGULAR_DEVIATIONS_SENTINEL` for invalid/unassigned genes
        real(real64), dimension(n_families), intent(in) :: family_dispersions
            !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\) from [[tox_get_outliers_by_angle_rap(module):tox_compute_family_direction_rap(subroutine)]], in range `[min_family_dispersion, max_family_dispersion]` with `CM_FAMILY_DISPERSION_SENTINEL` for invalid families
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), dimension(n_genes), intent(out) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes

        integer(int32) :: i_gene, fam_idx

        z_scores = CM_Z_SCORES_SENTINEL

    do concurrent(i_gene=1:n_genes) local(fam_idx) shared(n_families, gene_to_fam, angular_deviations, family_dispersions, z_scores)
            fam_idx = gene_to_fam(i_gene)

            if (fam_idx < 1 .or. fam_idx > n_families) then
                cycle
            end if

            if (angular_deviations(i_gene) < 0.0_real64 .or. family_dispersions(fam_idx) < 0.0_real64) then
                cycle
            end if

            if (is_close(family_dispersions(fam_idx), 0.0_real64)) then
                cycle
            end if

            z_scores(i_gene) = angular_deviations(i_gene) / family_dispersions(fam_idx)
        end do

    end subroutine tox_z_scores_by_dispersion_rap_helper

    !> Identify directional outliers based on scaled angles using percentile threshold
    !| \[is\_outlier(i\_gene) = z\_scores(i\_gene) \ge threshold\]
    pure subroutine tox_angle_outliers_rap(z_scores, threshold, is_outlier, n_genes, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        real(real64), intent(in) :: threshold
            !! `z_scores` exceeding this value will be identified as outliers
        real(real64), dimension(n_genes), intent(in) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output outlier status
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success

        ! Just call the common implementation
        call tox_angle_outliers(z_scores, threshold, is_outlier, n_genes, ierr)
    end subroutine tox_angle_outliers_rap

    !> Identify directional outliers based on scaled angles using percentile threshold
    !| \[is\_outlier(i\_gene) = z\_scores(i\_gene) \ge threshold\]
    pure subroutine tox_angle_outliers_rap_alloc(z_scores, percentile, threshold, is_outlier, n_genes, ierr)
        integer(int32), intent(in) :: n_genes
            !! Number of genes
        real(real64), dimension(n_genes), intent(in) :: z_scores
            !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
        real(real64), intent(in) :: percentile
            !! Percentile threshold (`0.0 < p < 100.0`)
        real(real64), intent(out) :: threshold
            !! Computed threshold value
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output outlier status
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success

        ! Just call the common implementation
        call tox_angle_outliers_alloc(z_scores, percentile, threshold, is_outlier, n_genes, ierr)
    end subroutine tox_angle_outliers_rap_alloc

    !> Complete pipeline for RAP angle-based outlier detection
    !| Full pipeline:
    !|
    !| 1. [[tox_get_outliers_by_angle_rap(module):tox_compute_family_direction_rap(subroutine)]]
    !| 2. [[tox_get_outliers_by_angle_rap(module):tox_compute_angular_deviations_rap(subroutine)]]
    !| 3. [[tox_get_outliers_by_angle_rap(module):tox_z_scores_by_dispersion_rap(subroutine)]]
    !| 4. [[tox_get_outliers_by_angle_rap(module):tox_angle_outliers_rap_alloc(subroutine)]]
    !|
    pure subroutine tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                           percentile, threshold, z_scores, &
                                                           n_genes, n_families, &
                                                           is_outlier, ierr, status, min_family_dispersion, max_family_dispersion)
        integer(int32), intent(in) :: n_genes
            !! Number of Genes
        integer(int32), intent(in) :: n_families
            !! Number of Families
        real(real64), dimension(n_genes), intent(in) :: rap_angles
            !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
        real(real64), intent(in) :: percentile
            !! Percentile threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
            !! Threshold value calculated with the percentile
        real(real64), dimension(n_genes), intent(out) :: z_scores
            !! Output scaled angles
        logical, dimension(n_genes), intent(out) :: is_outlier
            !! Output outlier status
        integer(int32), intent(out) :: ierr
            !! Error code, 0 on success
        integer(int32), dimension(n_families), intent(out) :: status
            !! Informs about special families
            !!
            !! |            Code             |                                                               Meaning                                                             |
            !! |-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
            !! |          ERR_OK           |                                                   Normal family, no special case.                                                 |
            !! | STAT_NO_STABLE_DIRECTION  |             The `rap_angles` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
            !! | STAT_NO_ANGULAR_VARIATION | The `rap_angles` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
            !!
        real(real64), intent(in), optional :: min_family_dispersion
            !! The minimum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MIN_FAMILY_DISPERSION_DEFAULT`
        real(real64), intent(in), optional :: max_family_dispersion
            !! The maximum family dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. default: `CM_MAX_FAMILY_DISPERSION_DEFAULT`

        ! Local variables
        real(real64), allocatable :: family_mean_angles(:)
        real(real64), allocatable :: family_dispersions(:)
        real(real64), allocatable :: angular_deviations(:)

        call set_ok(ierr)
        call set_ok(status)
        is_outlier = .false.

        ! Input validation
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_in_range_real(percentile, ierr, min=0.0_real64, max=100.0_real64)
        if (is_err(ierr)) return

        ! Allocate arrays
        M_ALLOCATE(family_mean_angles(n_families))
        M_ALLOCATE(family_dispersions(n_families))
        M_ALLOCATE(angular_deviations(n_genes))

        ! Step 1: Compute family circular statistics
        call tox_compute_family_direction_rap_alloc(rap_angles, gene_to_fam, &
                                                    family_mean_angles, family_dispersions, &
                                                    n_genes, n_families, ierr, status, min_family_dispersion, max_family_dispersion)
        if (is_err(ierr)) return

        ! Step 2: Compute angular deviations
        call tox_compute_angular_deviations_rap_helper(rap_angles, family_mean_angles, &
                                                       gene_to_fam, angular_deviations, &
                                                       n_genes, n_families)

        ! Step 3: Scale by dispersion
        call tox_z_scores_by_dispersion_rap_helper(angular_deviations, family_dispersions, &
                                                   gene_to_fam, z_scores, &
                                                   n_genes, n_families)

        ! Step 4: Detect outliers using the common function
        call tox_angle_outliers_rap_alloc(z_scores, percentile, threshold, is_outlier, n_genes, ierr)
    end subroutine tox_detect_angle_outliers_pipeline_rap

end module tox_get_outliers_by_angle_rap

!> C wrapper for the complete RAP pipeline
pure subroutine tox_detect_angle_outliers_pipeline_rap_c(rap_angles, gene_to_fam, &
                                                         percentile, threshold, z_scores, &
                                                         n_genes, n_families, &
                                                         is_outlier, ierr, status, min_family_dispersion, max_family_dispersion) &
    bind(C, name="tox_detect_angle_outliers_pipeline_rap_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle_rap, only: tox_detect_angle_outliers_pipeline_rap
    use tox_conversions, only: logical_as_c_int
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of families
    real(c_double), dimension(n_genes), intent(in), target :: rap_angles
        !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    real(c_double), intent(in), target :: percentile
        !! Percentile threshold
    real(c_double), intent(out), target :: threshold
        !! Threshold value calculated with the percentile
    real(c_double), dimension(n_genes), intent(out), target :: z_scores
        !! Output scaled angles
    integer(c_int), dimension(n_genes), intent(out), target :: is_outlier
        !! Output outlier status
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), dimension(n_families), intent(out), target :: status
        !! Informs about special families
        !!
        !! |            Code             |                                                               Meaning                                                             |
        !! |-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
        !! |          ERR_OK           |                                                   Normal family, no special case.                                                 |
        !! | STAT_NO_STABLE_DIRECTION  |             The `rap_angles` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
        !! | STAT_NO_ANGULAR_VARIATION | The `rap_angles` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
        !!
    real(c_double), intent(in), target :: min_family_dispersion
        !! The minimum angular dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. recommended: `CM_MIN_FAMILY_DISPERSION_DEFAULT`
    real(c_double), intent(in), target :: max_family_dispersion
        !! The maximum angular dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. recommended: `CM_MAX_FAMILY_DISPERSION_DEFAULT`

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(rap_angles)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(is_outlier)
    M_CHECK_NON_NULL(z_scores)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(threshold)
    M_CHECK_NON_NULL(status)
    M_CHECK_NON_NULL(min_family_dispersion)
    M_CHECK_NON_NULL(max_family_dispersion)

    M_ALLOCATE(is_outlier_f(n_genes))

    call tox_detect_angle_outliers_pipeline_rap(rap_angles, gene_to_fam, &
                                                percentile, threshold, z_scores, &
                                                n_genes, n_families, &
                                                is_outlier_f, ierr, status, min_family_dispersion, max_family_dispersion)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine tox_detect_angle_outliers_pipeline_rap_c

!> C wrapper for tox_compute_family_direction_rap
pure subroutine tox_compute_family_direction_rap_c(rap_angles, gene_to_fam, &
                                                   family_mean_angles, family_dispersions, &
                                                   n_genes, n_families, &
                                                   ierr, status, min_family_dispersion, max_family_dispersion) &
    bind(C, name="tox_compute_family_direction_rap_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle_rap, only: tox_compute_family_direction_rap_alloc
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_genes), intent(in), target :: rap_angles
        !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    real(c_double), dimension(n_families), intent(out), target :: family_mean_angles
        !! Family mean angles in radians, `0<=x<=PI`. The circular mean is kept for every family with >=2 members and a non-zero resultant length (even when the family is flagged via the dispersion bounds); it is `0` only for families with fewer than 2 members or a zero resultant length
    real(c_double), dimension(n_families), intent(out), target :: family_dispersions
        !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_FAMILY_DISPERSION_SENTINEL` if there is no stabe direction or too less angular variation (`status`)
    integer(c_int), intent(out), target :: ierr
        !! Error code, 0 on success
    integer(c_int), dimension(n_families), intent(out), target :: status
        !! Informs about special families
        !!
        !! |            Code             |                                                               Meaning                                                             |
        !! |-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
        !! |          ERR_OK           |                                                   Normal family, no special case.                                                 |
        !! | STAT_NO_STABLE_DIRECTION  |             The `rap_angles` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
        !! | STAT_NO_ANGULAR_VARIATION | The `rap_angles` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
        !!
    real(c_double), intent(in), target :: min_family_dispersion
        !! The minimum angular dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. recommended: `CM_MIN_FAMILY_DISPERSION_DEFAULT`
    real(c_double), intent(in), target :: max_family_dispersion
        !! The maximum angular dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. recommended: `CM_MAX_FAMILY_DISPERSION_DEFAULT`

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(rap_angles)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(family_mean_angles)
    M_CHECK_NON_NULL(family_dispersions)
    M_CHECK_NON_NULL(status)
    M_CHECK_NON_NULL(min_family_dispersion)
    M_CHECK_NON_NULL(max_family_dispersion)

    call tox_compute_family_direction_rap_alloc(rap_angles, gene_to_fam, &
                                                family_mean_angles, family_dispersions, &
                                                n_genes, n_families, &
                                                ierr, status, min_family_dispersion, max_family_dispersion)
end subroutine tox_compute_family_direction_rap_c

!> Expert C wrapper for tox_compute_family_direction_rap
pure subroutine tox_compute_family_direction_rap_expert_c(rap_angles, gene_to_fam, &
                                                          family_mean_angles, family_dispersions, &
                                                          n_genes, n_families, &
                                                    ierr, status, tmp_member_counts, min_family_dispersion, max_family_dispersion) &
    bind(C, name="tox_compute_family_direction_rap_expert_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle_rap, only: tox_compute_family_direction_rap
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_genes), intent(in), target :: rap_angles
        !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    real(c_double), dimension(n_families), intent(out), target :: family_mean_angles
        !! Family mean angles in radians, `0<=x<=PI`. The circular mean is kept for every family with >=2 members and a non-zero resultant length (even when the family is flagged via the dispersion bounds); it is `0` only for families with fewer than 2 members or a zero resultant length
    real(c_double), dimension(n_families), intent(out), target :: family_dispersions
        !! Angular dispersion: \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\), `CM_FAMILY_DISPERSION_SENTINEL` if there is no stabe direction or too less angular variation (`status`)
    integer(c_int), intent(out), target :: ierr
        !! Error code, 0 on success
    integer(c_int), dimension(n_families), intent(out), target :: status
        !! Informs about special families
        !!
        !! |            Code             |                                                               Meaning                                                             |
        !! |-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
        !! |          ERR_OK           |                                                   Normal family, no special case.                                                 |
        !! | STAT_NO_STABLE_DIRECTION  |             The `rap_angles` have a too large spread, meaning dispersion being greater than `max_angular_dispersion`              |
        !! | STAT_NO_ANGULAR_VARIATION | The `rap_angles` have a too small spread, meaning being too similar and thus dispersion being lower than `min_angular_dispersion` |
        !!
    integer(c_int), dimension(n_families), intent(out), target :: tmp_member_counts
        !! Number of members per family, determined by `gene_to_fam`
    real(c_double), intent(in), target :: min_family_dispersion
        !! The minimum angular dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. recommended: `CM_MIN_FAMILY_DISPERSION_DEFAULT`
    real(c_double), intent(in), target :: max_family_dispersion
        !! The maximum angular dispersion, otherwise `CM_FAMILY_DISPERSION_SENTINEL` is used. recommended: `CM_MAX_FAMILY_DISPERSION_DEFAULT`

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(rap_angles)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(tmp_member_counts)
    M_CHECK_NON_NULL(family_mean_angles)
    M_CHECK_NON_NULL(family_dispersions)
    M_CHECK_NON_NULL(status)
    M_CHECK_NON_NULL(min_family_dispersion)
    M_CHECK_NON_NULL(max_family_dispersion)

    call tox_compute_family_direction_rap(rap_angles, gene_to_fam, &
                                          family_mean_angles, family_dispersions, &
                                          n_genes, n_families, &
                                          ierr, status, tmp_member_counts, min_family_dispersion, max_family_dispersion)
end subroutine tox_compute_family_direction_rap_expert_c

!> C wrapper for tox_compute_angular_deviations_rap
pure subroutine tox_compute_angular_deviations_rap_c(rap_angles, family_mean_angles, &
                                                     gene_to_fam, angular_deviations, &
                                                     n_genes, n_families, ierr) &
    bind(C, name="tox_compute_angular_deviations_rap_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle_rap, only: tox_compute_angular_deviations_rap
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_genes), intent(in), target :: rap_angles
        !! RAP angles in radians, `0<=x<=PI` and `CM_RAP_ANGLES_SENTINEL` for invalid/unassigned genes
    real(c_double), dimension(n_families), intent(in), target :: family_mean_angles
        !! Family mean angles in radians, `0<=x<=PI`
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping, `M_GENE_TO_FAM_SENTINEL` for unassigned
    real(c_double), dimension(n_genes), intent(out), target :: angular_deviations
        !! Absolute wrapped angular deviations [0, π], `CM_ANGULAR_DEVIATIONS_SENTINEL` for invalid
    integer(c_int), intent(out), target :: ierr
        !! Error code, 0 on success

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(rap_angles)
    M_CHECK_NON_NULL(family_mean_angles)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(angular_deviations)

    call tox_compute_angular_deviations_rap(rap_angles, family_mean_angles, &
                                            gene_to_fam, angular_deviations, &
                                            n_genes, n_families, ierr)
end subroutine tox_compute_angular_deviations_rap_c

!> C wrapper for tox_z_scores_by_dispersion_rap
pure subroutine tox_z_scores_by_dispersion_rap_c(angular_deviations, family_dispersions, &
                                                 gene_to_fam, z_scores, &
                                                 n_genes, n_families, ierr) &
    bind(C, name="tox_z_scores_by_dispersion_rap_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle_rap, only: tox_z_scores_by_dispersion_rap
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes
    integer(c_int), intent(in), target :: n_families
        !! Number of gene families
    real(c_double), dimension(n_genes), intent(in), target :: angular_deviations
        !! Angular deviations in radians `0<=x<=PI` and `CM_ANGULAR_DEVIATIONS_SENTINEL` for invalid/unassigned genes
    real(c_double), dimension(n_families), intent(in), target :: family_dispersions
        !! Angular dispersions \(\sigma_{\phi, F} = \sqrt{-2 \ln(R_F)}\) from [[tox_get_outliers_by_angle_rap(module):tox_compute_family_direction_rap(subroutine)]], in range `[min_family_dispersion, max_family_dispersion]` with `CM_FAMILY_DISPERSION_SENTINEL` for invalid families
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! Gene to family mapping
    real(c_double), dimension(n_genes), intent(out), target :: z_scores
        !! Scaled angles (z-scores), `x>=0` and `CM_Z_SCORES_SENTINEL` for invalid/unassigned genes
    integer(c_int), intent(out), target :: ierr
        !! Error code, 0 on success

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(angular_deviations)
    M_CHECK_NON_NULL(family_dispersions)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(z_scores)

    call tox_z_scores_by_dispersion_rap(angular_deviations, family_dispersions, &
                                        gene_to_fam, z_scores, &
                                        n_genes, n_families, ierr)
end subroutine tox_z_scores_by_dispersion_rap_c

!> C wrapper for tox_angle_outliers_rap_alloc
pure subroutine tox_angle_outliers_rap_c(z_scores, n_genes, percentile, &
                                         threshold, is_outlier, ierr) &
    bind(C, name="tox_angle_outliers_rap_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle_rap, only: tox_angle_outliers_rap_alloc
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
        !! Computed threshold value
    integer(c_int), dimension(n_genes), intent(out), target :: is_outlier
        !! Output outlier status
    integer(c_int), intent(out), target :: ierr
        !! Error code, 0 on success

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(z_scores)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(threshold)
    M_CHECK_NON_NULL(is_outlier)

    M_ALLOCATE(is_outlier_f(n_genes))

    call tox_angle_outliers_rap_alloc(z_scores, percentile, threshold, &
                                      is_outlier_f, n_genes, ierr)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine tox_angle_outliers_rap_c

!> Expert C wrapper for tox_angle_outliers_rap
pure subroutine tox_angle_outliers_rap_expert_c(z_scores, threshold, is_outlier, n_genes, ierr) &
    bind(C, name="tox_angle_outliers_rap_expert_c")
    use iso_c_binding, only: c_int, c_double
    use tox_get_outliers_by_angle_rap, only: tox_angle_outliers_rap
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
        !! Output outlier status
    integer(c_int), intent(out), target :: ierr
        !! Error code, 0 on success

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(z_scores)
    M_CHECK_NON_NULL(threshold)
    M_CHECK_NON_NULL(is_outlier)

    M_ALLOCATE(is_outlier_f(n_genes))

    call tox_angle_outliers_rap(z_scores, threshold, is_outlier_f, n_genes, ierr)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine tox_angle_outliers_rap_expert_c
