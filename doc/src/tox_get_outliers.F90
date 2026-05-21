#include "macros.h"

!> Module to identify gene outliers based on their distances to family centroids.
module tox_get_outliers
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use f42_utils, only: sort_array, calc_percentile, logx, is_close, compute_empirical_p_values
    use tox_errors, only: ERR_INVALID_INPUT, set_ok, set_err_once, check_alloc_stat, is_err
    use tox_loess, only: tox_loess_required_workspace, loess_fit_robust, loess_fit_plain, EPS_LOESS, loess_evaluation
    implicit none

contains

    !> Compute family scaling factors (dscale) to normalize distances.
    !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
    subroutine compute_family_scaling( &
        n_genes, n_families, distances, gene_to_fam, dscale, &
        loess_x, loess_y, indices_used, tmp_perm, tmp_stack_left, tmp_stack_right, &
        tmp_iv, liv, tmp_wv, lv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat, &
        span, degree, mode, n_iters, low_sd_cutoff, excluded_low_sd, tmp_means_aux, ierr)

        use, intrinsic :: iso_fortran_env, only: real64, int32
        implicit none

        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! Mapping of each gene to its family (1-based)

        real(real64), intent(out) :: dscale(n_families)
            !! Output: array of scaling factors per family

        ! Buffers (reused)
        real(real64), intent(out) :: loess_x(n_families)
            !! Reference x-coordinates for LOESS smoothing
        real(real64), intent(out) :: loess_y(n_families)
            !! Reference y-coordinates for LOESS smoothing
        integer(int32), intent(inout) :: indices_used(n_families)
            !! Indices of reference points used for smoothing
        integer(int32), intent(out) :: tmp_perm(n_genes)
            !! Permutation array for sorting gene distances
        integer(int32), intent(out) :: tmp_stack_left(n_genes)
            !! Stack array for left indices during sorting
        integer(int32), intent(out) :: tmp_stack_right(n_genes)
            !! Stack array for right indices during sorting
        real(real64), intent(inout) :: tmp_means_aux(n_families)
            !! Work array for saving raw means
        integer(int32), intent(out) :: excluded_low_sd(n_families)
            !! Mask to save those families that have low sd

        ! LOESS workspace
        integer(int32), intent(in)    :: liv
            !! Length of integer workspace
        integer(int32), intent(inout) :: tmp_iv(liv)
            !! Integer workspace array
        integer(int32), intent(in)    :: lv
            !! Length of real workspace
        real(real64), intent(inout) :: tmp_wv(lv)
            !! Real workspace array

        real(real64), intent(inout) :: tmp_diagl(n_families)
            !! Diagonal elements of the weight matrix
        real(real64), intent(inout) :: tmp_w_init(n_families)
            !! Initial weights for LOESS
        real(real64), intent(inout) :: tmp_z_mat(n_families, 1)
            !! Z matrix for LOESS fitting
        real(real64), intent(inout) :: tmp_rw(n_families)
            !! Residuals for robust LOESS fitting
        real(real64), dimension(n_families), intent(out), target :: tmp_ww
            !! Working weights array
        real(real64), dimension(n_families), intent(out), target :: tmp_res
            !! Residuals array
        integer(int32), intent(inout):: tmp_pi(:)
            !! Permutation indices for robust LOESS fitting
        real(real64), intent(out)    :: tmp_yhat(:)
            !! Output array for LOESS predictions

        real(real64), intent(in)     :: span
            !! Span parameter for LOESS smoothing
        integer(int32), intent(in)   :: degree
            !! Degree of the LOESS polynomial
        integer(int32), intent(in)   :: mode
            !! Mode for LOESS fitting (0=plain, 1=robust)
        integer(int32), intent(in)   :: n_iters
            !! Number of iterations for robust LOESS fitting
        real(real64), intent(out) :: low_sd_cutoff
            !! cutoff used to filter families with low std
        integer(int32), intent(out)  :: ierr
            !! Error code

        ! Local variables
        integer(int32) :: i_gene, i_family, i_valid, family_idx, n_in_family, n_valid, k, tmp_ierr
        real(real64)   :: stddev_dist, mean_dist, sumsq, dist_val
        real(real64) :: xmin, xmax, eps_mean, eps_sd, std_median

        ! Initialize error code and output arrays
        call set_ok(ierr)
        dscale  = 0.0_real64
        loess_x = 0.0_real64
        loess_y = 0.0_real64
        n_valid = 0

        ! Validate family indices
        do i_gene = 1, n_genes
            if (gene_to_fam(i_gene) < 1 .or. gene_to_fam(i_gene) > n_families) then
                dscale = -1.0_real64
                call set_err_once(ierr, ERR_INVALID_INPUT)
                return
            end if
        end do
        tmp_means_aux = -1.0_real64

        ! ------------------------------------------------------------
        ! PASS 1: compute (mean, stddev) per family
        ! ------------------------------------------------------------

        tmp_w_init = 0.0_real64
        tmp_rw = 0.0_real64
        tmp_pi = 0

        do i_gene = 1, n_genes
            family_idx = gene_to_fam(i_gene)
            dist_val = abs(distances(i_gene))

            tmp_pi(family_idx) = tmp_pi(family_idx) + 1
            tmp_w_init(family_idx) = tmp_w_init(family_idx) + dist_val
            tmp_rw(family_idx) = tmp_rw(family_idx) + (dist_val**2)
        end do

        n_valid = 0
        do i_family = 1, n_families
            n_in_family = tmp_pi(i_family)

            if (n_in_family <= 1) cycle

            n_valid = n_valid + 1

            mean_dist = tmp_w_init(i_family)/real(n_in_family, real64)

            ! Var = (SumSq - (Sum^2)/N) / (N-1)
            sumsq = max(0.0_real64, tmp_rw(i_family) - (tmp_w_init(i_family)**2/real(n_in_family, real64)))
            stddev_dist = sqrt(sumsq/real(n_in_family - 1, real64))

            loess_x(n_valid) = mean_dist
            tmp_means_aux(i_family) = mean_dist
            loess_y(n_valid) = stddev_dist
            indices_used(n_valid) = i_family
        end do

        tmp_w_init = 0.0_real64
        tmp_rw = 0.0_real64
        tmp_pi = 0

        if (n_valid <= 1) then
            low_sd_cutoff = 0.0_real64
            return
        end if

        do i_valid = 1, n_valid
            tmp_perm(i_valid) = i_valid
        end do

        call sort_array(loess_x(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        call calc_percentile(loess_x(1:n_valid), tmp_perm(1:n_valid), 5.0_real64, eps_mean, ierr)
        if (is_err(ierr)) return

        eps_mean = max(eps_mean, EPS_LOESS)

        call sort_array(loess_y(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        if (mod(n_valid, 2) == 0) then
            std_median = 0.5_real64*( &
                         loess_y(tmp_perm(n_valid/2)) + &
                         loess_y(tmp_perm(n_valid/2 + 1)))
        else
            std_median = loess_y(tmp_perm((n_valid + 1)/2))
        end if

        eps_sd = max(1.0e-13_real64*std_median, EPS_LOESS)

        do concurrent (i_valid = 1:n_valid) local(tmp_ierr) shared(loess_x, eps_mean, loess_y, eps_sd, ierr)
            call logx(loess_x(i_valid) + eps_mean, 2.0_real64, loess_x(i_valid), tmp_ierr)
            if (is_err(tmp_ierr)) ierr = tmp_ierr
            call logx(loess_y(i_valid) + eps_sd, 2.0_real64, loess_y(i_valid), tmp_ierr)
            if (is_err(tmp_ierr)) ierr = tmp_ierr
        end do

        if (is_err(ierr)) return

        call sort_array(loess_y(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        call calc_percentile(loess_y(1:n_valid), tmp_perm(1:n_valid), 1.0_real64, low_sd_cutoff, ierr)

        if (is_err(ierr)) return

        excluded_low_sd = 1_int32
        k = 0

        do i_valid = 1, n_valid
            if (loess_y(i_valid) > low_sd_cutoff .or. is_close(loess_y(i_valid), low_sd_cutoff)) then
                family_idx = indices_used(i_valid)
                excluded_low_sd(family_idx) = 0_int32
                k = k + 1
                loess_x(k) = loess_x(i_valid)
                loess_y(k) = loess_y(i_valid)
                indices_used(k) = indices_used(i_valid)
            end if
        end do

        n_valid = k

        ! Trigger fallback case when having too few points
        if (n_valid <= 1) then
            xmin = 0.0_real64
            xmax = xmin
        else
            xmin = minval(loess_x(1:n_valid))
            xmax = maxval(loess_x(1:n_valid))
        end if

        ! Fallback: for constant prediction or not enough points for loess
        ! we use the global median for all families.
        if (xmax == xmin) then
            do concurrent (i_family = 1:n_families) shared(dscale, std_median)
                if (tmp_means_aux(i_family) >= 0.0_real64) then
                    dscale(i_family) = std_median
                else
                    dscale(i_family) = 0.0_real64
                end if
            end do

            low_sd_cutoff = max(2.0_real64**low_sd_cutoff - eps_sd, 0.0_real64)
            return
        end if

        ! ------------------------------------------------------------
        ! LOESS GLOBAL: smooth y_ref as function of x_ref (once)
        ! ------------------------------------------------------------
        tmp_w_init(1:n_valid) = 1.0_real64
        tmp_z_mat(1:n_valid, 1) = loess_x(1:n_valid)

        if (mode == 0) then
            ! If you have a plain routine, call it; otherwise keep robust always.
            call loess_fit_plain( &
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_w_init(1:n_valid), tmp_z_mat(1:n_valid, 1:1), &
                span, degree, n_valid, .false., .false., tmp_iv, liv, tmp_wv, lv, tmp_diagl(1:n_valid), tmp_yhat(1:n_valid), ierr)
        else
            call loess_fit_robust( &
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_w_init(1:n_valid), tmp_z_mat(1:n_valid, 1:1), &
                span, degree, n_valid, .false., .false., n_iters, tmp_iv, liv, tmp_wv, lv, tmp_diagl(1:n_valid), &
                tmp_rw(1:n_valid), tmp_ww(1:n_valid), tmp_res(1:n_valid), tmp_pi(1:n_valid), tmp_yhat(1:n_valid), ierr)
        end if

        if (is_err(ierr)) return

        do concurrent (i_family = 1:n_families) local(tmp_ierr) shared(tmp_means_aux, eps_mean, tmp_z_mat, xmin, xmax, ierr)
            if (tmp_means_aux(i_family) >= 0.0_real64) then
                call logx(tmp_means_aux(i_family) + eps_mean, 2.0_real64, tmp_z_mat(i_family, 1), tmp_ierr)
                if (is_err(tmp_ierr)) then
                    ierr = tmp_ierr
                else
                    if (tmp_z_mat(i_family, 1) < xmin) tmp_z_mat(i_family, 1) = xmin
                    if (tmp_z_mat(i_family, 1) > xmax) tmp_z_mat(i_family, 1) = xmax
                end if
            else
                tmp_z_mat(i_family, 1) = xmin
            end if
        end do

        if (is_err(ierr)) return

        call loess_evaluation(tmp_iv, liv, lv, tmp_wv, n_families, tmp_z_mat(1:n_families, 1:1), tmp_yhat(1:n_families))

        if (is_err(ierr)) return

        ! ------------------------------------------------------------
        ! Map compact tmp_results back to full dscale
        ! ------------------------------------------------------------

        do concurrent (i_family = 1:n_families) shared(tmp_means_aux, dscale, tmp_yhat, eps_sd)
            if (tmp_means_aux(i_family) < 0.0_real64) then
                dscale(i_family) = 0.0_real64
            else
                dscale(i_family) = max(2.0_real64**tmp_yhat(i_family) - eps_sd, 0.0_real64)
            end if
        end do

        do concurrent (i_valid = 1:n_valid) shared(loess_x, loess_y, eps_mean, eps_sd)
            ! linear scale for return
            loess_x(i_valid) = max(2.0_real64**loess_x(i_valid) - eps_mean, 0.0_real64)
            loess_y(i_valid) = max(2.0_real64**loess_y(i_valid) - eps_sd, 0.0_real64)
        end do

        if (n_valid < n_families) then
            loess_x(n_valid + 1:) = M_NAN
            loess_y(n_valid + 1:) = M_NAN
            indices_used(n_valid + 1:n_families) = 0_int32
        end if

        low_sd_cutoff = max(2.0_real64**low_sd_cutoff - eps_sd, 0.0_real64)

    end subroutine compute_family_scaling

    !> Helper routine that allocates internal arrays and calls compute_family_scaling.
    !| This makes usage easier since users don't need to care about internal array requirements.
    subroutine compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                            loess_x, loess_y, indices_used, ierr)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! Mapping of each gene to its family (1-based)
        real(real64), intent(out) :: dscale(n_families)
            !! Output: array of scaling factors per family
        real(real64), intent(inout) :: loess_x(n_families)
            !! Reference x-coordinates.
        real(real64), intent(inout) :: loess_y(n_families)
            !! Reference y-coordinates (length n_total).
        integer(int32), intent(inout) :: indices_used(n_families)
            !! Indices of reference points used for smoothing.
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Local work arrays
        integer(int32), allocatable :: tmp_perm(:)
        integer(int32), allocatable :: tmp_stack_left(:)
        integer(int32), allocatable :: tmp_stack_right(:)
        real(real64) :: low_sd_cutoff
        integer(int32), allocatable :: excluded_low_sd(:)
        real(real64), allocatable :: tmp_means_aux(:)

        ! LOESS workspace
        integer(int32) :: liv, lv, istat
        integer(int32), allocatable :: tmp_iv(:), tmp_pi(:)
        real(real64), allocatable :: tmp_wv(:), tmp_diagl(:), tmp_w_init(:), tmp_z_mat(:, :), tmp_rw(:), tmp_ww(:), tmp_res(:), tmp_yhat(:)

        ! LOESS params (defaults)
        real(real64), parameter :: span = 0.7_real64
        integer(int32), parameter :: degree = 2_int32
        integer(int32), parameter :: mode = 1_int32
        integer(int32), parameter :: n_iters = 3_int32
        logical, parameter :: setlf = .false.

        call set_ok(ierr)
        call set_ok(istat)

        ! Workspace sizes
        call tox_loess_required_workspace(1_int32, n_families, liv, lv, setlf)

        allocate (tmp_iv(liv), tmp_wv(lv), stat=istat)
        call check_alloc_stat(istat, ierr)
        if (is_err(ierr)) return

        ! For robust we also need arrays sized to n_valid (<= n_families).
        allocate (tmp_diagl(n_families), tmp_w_init(n_families), tmp_z_mat(n_families, 1), &
                  tmp_rw(n_families), tmp_ww(n_families), tmp_res(n_families), tmp_pi(n_families), &
                  tmp_yhat(n_families), tmp_perm(n_genes), tmp_stack_left(n_genes), &
                  tmp_stack_right(n_genes), excluded_low_sd(n_families), tmp_means_aux(n_families), stat=istat)

        call check_alloc_stat(istat, ierr)
        if (is_err(ierr)) return

        ! Initialize (important for netlib)
        tmp_iv = 1_int32
        tmp_wv = 0.0_real64
        tmp_rw = 1.0_real64
        tmp_pi = 0_int32

        call compute_family_scaling( &
            n_genes, n_families, distances, gene_to_fam, dscale, &
            loess_x, loess_y, indices_used, tmp_perm, tmp_stack_left, tmp_stack_right, &
            tmp_iv, liv, tmp_wv, lv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat, &
            span, degree, mode, n_iters, low_sd_cutoff, excluded_low_sd, tmp_means_aux, ierr)

        deallocate (tmp_iv, tmp_wv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat)

    end subroutine compute_family_scaling_alloc

    !> Compute the hybrid RDI for each gene.
    !| RDI = Euclidean distance / family scaling factor
    pure subroutine compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, &
                                stack_left, stack_right)
        implicit none
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene to its centroid
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! Gene-to-family mapping (1-based indexing)
        real(real64), intent(in) :: dscale(:)
            !! Array of scaling factors for each family
        real(real64), intent(out) :: rdi(n_genes)
            !! Output array of RDI values for each gene
        real(real64), intent(inout) :: sorted_rdi(n_genes)
            !! Work array for sorting (dimension n_genes)
        integer(int32), intent(inout) :: perm(n_genes)
            !! Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
        integer(int32), intent(inout) :: stack_left(n_genes)
            !! Stack array for sorting (dimension n_genes)
        integer(int32), intent(inout) :: stack_right(n_genes)
            !! Stack array for sorting (dimension n_genes)

        integer(int32) :: i, family_idx
        real(real64), parameter :: tol = epsilon(1.0_real64)

        ! Calculate RDI for each gene
        do i = 1, n_genes
            family_idx = gene_to_fam(i)

            ! Handle invalid family indices
            if (family_idx < 1 .or. family_idx > size(dscale)) then
                rdi(i) = -1.0_real64  ! Error indicator
                cycle
            end if

            ! Detect NaN input (portable)
            if (ieee_is_nan(distances(i))) then
                rdi(i) = distances(i)
            else if (abs(dscale(family_idx)) < tol) then
                rdi(i) = 0.0_real64  ! If scaling is zero, set RDI to zero (not outlier)
            else
                ! Calculate RDI
                rdi(i) = abs(distances(i))/dscale(family_idx)
            end if
        end do

        ! Create a copy of RDI for sorting (excluding error values)
        sorted_rdi = rdi

        ! Filter out error values (negative RDIs)
        where (sorted_rdi < 0.0_real64)
            sorted_rdi = 0.0_real64
        end where

        ! Sort RDI values using the tox_sorting module
        call sort_array(sorted_rdi, perm, stack_left, stack_right)

    end subroutine compute_rdi

    !> Identify gene outliers based on the top percentile of RDI values.
    !| Expects sorted_rdi to be filtered (no negative values) and perm should be sorted in ascending order before calling.
    !| If sorted_rdi contains negatives or perm is not sorted, tmp_results may be invalid.
    pure subroutine identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, percentile)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        real(real64), intent(in) :: rdi(n_genes)
            !! Array of RDI values for each gene
        real(real64), intent(in) :: sorted_rdi(n_genes)
            !! Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
        integer(int32), intent(in) :: perm(n_genes)
            !! Permutation array with sorted indices
        logical, intent(out) :: is_outlier(n_genes)
            !! Output boolean array indicating outliers
        real(real64), intent(out) :: threshold
            !! Output threshold value used for detection
        real(real64), intent(in), optional :: percentile
            !! (optional) Percentile threshold (default: 95 for top 5%)
        real(real64), intent(out) :: p_values(n_genes)
            !! Empirical one-sided upper-tail p-values for each gene. Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided upper-tail empirical p-value is used.

        integer(int32) :: i, idx
        real(real64) :: perc_pos, percentile_val

        ! Set default percentile if not present
        if (present(percentile)) then
            percentile_val = percentile
        else
            percentile_val = 95.0_real64
        end if

        ! Initialize output
        is_outlier = .false.

        ! Calculate the position corresponding to the desired percentile
        perc_pos = (n_genes*percentile_val)/100.0_real64
        idx = ceiling(perc_pos)
        ! Clamp idx to valid range
        if (idx < 1) idx = 1
        if (idx > n_genes) idx = n_genes

        ! Get the threshold value from the sorted array (sorted_rdi must be ascending)
        threshold = sorted_rdi(perm(idx))

        ! Mark genes as outliers if their RDI exceeds the threshold (and is positive)
        do i = 1, n_genes
            is_outlier(i) = (rdi(i) >= threshold .and. rdi(i) > 0.0_real64)
        end do

        call compute_empirical_p_values(n_genes, rdi, sorted_rdi, perm, p_values, 1.0_real64)

    end subroutine identify_outliers

    !> Main routine to detect outliers using RDI and LOESS-based scaling.
    subroutine detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                               tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right, &
                               is_outlier, loess_x, loess_y, loess_n, p_values, ierr, &
                               percentile)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene to its centroid
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! Gene-to-family mapping (1-based indexing)
        real(real64), intent(out) :: tmp_work_array(n_genes)
            !! Work array for sorting (dimension n_genes)
        integer(int32), intent(out) :: tmp_perm(n_genes)
            !! Permutation array for sorting (dimension n_genes)
        integer(int32), intent(out) :: tmp_stack_left(n_genes)
            !! Stack array for left indices during sorting
        integer(int32), intent(out) :: tmp_stack_right(n_genes)
            !! Stack array for right indices during sorting
        logical, intent(out) :: is_outlier(n_genes)
            !! Output boolean array indicating outliers
        real(real64), intent(out) :: loess_x(n_families)
            !! Reference x-coordinates.
        real(real64), intent(out) :: loess_y(n_families)
            !! Reference y-coordinates (length n_total).
        integer(int32), intent(out) :: loess_n(n_families)
            !! Indices of reference points used for smoothing.
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), intent(in), optional :: percentile
            !! (optional) Percentile threshold for outlier detection (default: 95)
        real(real64), intent(out) :: p_values(n_genes)
            !! Empirical one-sided upper-tail p-values for each gene. Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided upper-tail empirical p-value is used.

        ! Local variables
        real(real64) :: dscale(n_families)
        real(real64) :: rdi(n_genes)
        real(real64) :: threshold
        integer(int32) :: i
        real(real64) :: percentile_val

        ! Set default percentile if not present
        if (present(percentile)) then
            percentile_val = percentile
        else
            percentile_val = 95.0_real64
        end if

        ! Always initialize permutation array
        do i = 1, n_genes
            tmp_perm(i) = i
        end do

        call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                          loess_x, loess_y, loess_n, ierr)
        if (is_err(ierr)) return
        call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right)
        call identify_outliers(n_genes, rdi, tmp_work_array, tmp_perm, is_outlier, threshold, p_values, percentile_val)
    end subroutine detect_outliers
end module tox_get_outliers

!> C wrapper for compute_family_scaling (main version with automatic allocation).
!| Calls compute_family_scaling_alloc with C-compatible types for external interface.
!| This is the recommended version for most users as it handles memory allocation automatically.
subroutine compute_family_scaling_c(n_genes, n_families, distances, gene_to_fam, dscale, &
                                    loess_x, loess_y, indices_used, ierr) bind(C, name="compute_family_scaling_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_get_outliers, only: compute_family_scaling_alloc
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes
    integer(c_int), intent(in), target :: n_families
        !! Total number of families
    real(c_double), intent(in), target :: distances(n_genes)
        !! Array of Euclidean distances for each gene
    integer(c_int), intent(in), target :: gene_to_fam(n_genes)
        !! Mapping of each gene to its family (1-based)
    real(c_double), intent(out), target :: dscale(n_families)
        !! Output: array of scaling factors per family
    real(c_double), intent(inout), target :: loess_x(n_families)
        !! Reference x-coordinates for LOESS
    real(c_double), intent(inout), target :: loess_y(n_families)
        !! Reference y-coordinates for LOESS
    integer(c_int), intent(inout), target :: indices_used(n_families)
        !! Indices of reference points used for smoothing
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(distances)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(dscale)
    M_CHECK_NON_NULL(loess_x)
    M_CHECK_NON_NULL(loess_y)
    M_CHECK_NON_NULL(indices_used)

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                      loess_x, loess_y, indices_used, ierr)
end subroutine compute_family_scaling_c

!> C wrapper for compute_rdi.
!| Calls compute_rdi with C-compatible types for external interface.
!| Outputs both unsorted and sorted RDI, permutation, and sorting workspace arrays for downstream use.
subroutine compute_rdi_c(n_genes, n_families, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right, ierr) bind(C, name="compute_rdi_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_get_outliers, only: compute_rdi
    use tox_errors, only: set_ok
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes
    integer(c_int), intent(in), target :: n_families
        !! Total number of families
    real(c_double), intent(in), target :: distances(n_genes)
        !! Array of Euclidean distances for each gene to its centroid
    integer(c_int), intent(in), target :: gene_to_fam(n_genes)
        !! Gene-to-family mapping (1-based indexing)
    real(c_double), intent(in), target :: dscale(n_families)
        !! Array of scaling factors for each family
    real(c_double), intent(out), target :: rdi(n_genes)
        !! Output array of RDI values for each gene (unsorted)
    real(c_double), intent(out), target :: sorted_rdi(n_genes)
        !! Output array of sorted RDI values (filtered, sorted)
    integer(c_int), intent(out), target :: perm(n_genes)
        !! Output permutation array for sorting (dimension n_genes)
    integer(c_int), intent(out), target :: stack_left(n_genes)
        !! Output stack array for left indices during sorting
    integer(c_int), intent(out), target :: stack_right(n_genes)
        !! Output stack array for right indices during sorting
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(distances)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(dscale)
    M_CHECK_NON_NULL(rdi)
    M_CHECK_NON_NULL(sorted_rdi)
    M_CHECK_NON_NULL(perm)
    M_CHECK_NON_NULL(stack_left)
    M_CHECK_NON_NULL(stack_right)

    call set_ok(ierr)
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
end subroutine compute_rdi_c

!> C wrapper for identify_outliers.
!| Calls identify_outliers with C-compatible types for external interface.
subroutine identify_outliers_c(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, p_values, percentile, ierr) &
    bind(C, name="identify_outliers_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_get_outliers, only: identify_outliers
    use tox_conversions, only: logical_as_c_int
    use tox_errors, only: set_ok
    M_USE_ALLOCATION
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes
    real(c_double), intent(in), target :: rdi(n_genes)
        !! Array of RDI values for each gene
    real(c_double), intent(in), target :: sorted_rdi(n_genes)
        !! Filtered RDI array (no negatives, no NaNs)
    integer(c_int), intent(inout), target :: perm(n_genes)
        !! Permutation array with sorted indices
    integer(c_int), intent(out), target :: is_outlier(n_genes)
        !! Output integer array indicating outliers (1=outlier, 0=not)
    real(c_double), intent(out), target :: threshold
        !! Output threshold value used for detection
    real(c_double), intent(in), target :: percentile
        !! Percentile threshold for outlier detection
    real(c_double), intent(out), target :: p_values(n_genes)
        !! Empirical one-sided upper-tail p-values for each gene. Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided upper-tail empirical p-value is used.
    integer(c_int), intent(out), target :: ierr
        !! Error code

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(rdi)
    M_CHECK_NON_NULL(sorted_rdi)
    M_CHECK_NON_NULL(perm)
    M_CHECK_NON_NULL(is_outlier)
    M_CHECK_NON_NULL(threshold)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(p_values)

    M_ALLOCATE(is_outlier_f(n_genes))

    call set_ok(ierr)

    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier_f, threshold, p_values, percentile)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine identify_outliers_c

!> C wrapper for detect_outliers.
!| Calls detect_outliers with C-compatible types for external interface.
subroutine detect_outliers_c(n_genes, n_families, distances, gene_to_fam, &
                             tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right, &
                             is_outlier, loess_x, loess_y, loess_n, p_values, ierr, &
                             percentile) bind(C, name="detect_outliers_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_get_outliers, only: detect_outliers
    use tox_conversions, only: logical_as_c_int
    use tox_errors, only: is_ok
    M_USE_ALLOCATION
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes
    integer(c_int), intent(in), target :: n_families
        !! Total number of families
    real(c_double), intent(in), target :: distances(n_genes)
        !! Array of Euclidean distances for each gene to its centroid
    integer(c_int), intent(in), target :: gene_to_fam(n_genes)
        !! Gene-to-family mapping (1-based indexing)
    real(c_double), intent(out), target :: tmp_work_array(n_genes)
        !! Work array for sorting (dimension n_genes)
    integer(c_int), intent(out), target :: tmp_perm(n_genes)
        !! Permutation array for sorting (dimension n_genes)
    integer(c_int), intent(out), target :: tmp_stack_left(n_genes)
        !! Stack array for left indices during sorting
    integer(c_int), intent(out), target :: tmp_stack_right(n_genes)
        !! Stack array for right indices during sorting
    integer(c_int), intent(out), target :: is_outlier(n_genes)
        !! Output integer array indicating outliers (1=outlier, 0=not)
    real(c_double), intent(out), target :: loess_x(n_families)
        !! Reference x-coordinates for LOESS
    real(c_double), intent(out), target :: loess_y(n_families)
        !! Reference y-coordinates for LOESS
    integer(c_int), intent(out), target :: loess_n(n_families)
        !! Indices of reference points used for smoothing
    real(c_double), intent(out), target :: p_values(n_genes)
        !! Empirical one-sided upper-tail p-values for each gene. Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided upper-tail empirical p-value is used.
    integer(c_int), intent(out), target :: ierr
        !! Error code
    real(c_double), intent(in), target :: percentile
        !! Percentile threshold for outlier detection

    logical, dimension(:), allocatable :: is_outlier_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(distances)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(tmp_work_array)
    M_CHECK_NON_NULL(tmp_perm)
    M_CHECK_NON_NULL(tmp_stack_left)
    M_CHECK_NON_NULL(tmp_stack_right)
    M_CHECK_NON_NULL(is_outlier)
    M_CHECK_NON_NULL(loess_x)
    M_CHECK_NON_NULL(loess_y)
    M_CHECK_NON_NULL(loess_n)
    M_CHECK_NON_NULL(p_values)
    M_CHECK_NON_NULL(percentile)

    M_ALLOCATE(is_outlier_f(n_genes))

    call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                         tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right, &
                         is_outlier_f, loess_x, loess_y, loess_n, p_values, ierr, &
                         percentile)

    if (is_ok(ierr)) then
        call logical_as_c_int(is_outlier_f, is_outlier)
    end if
end subroutine detect_outliers_c

!> C wrapper for compute_family_scaling expert version.
!| Calls compute_family_scaling with C-compatible types for external interface.
!| This wrapper is designed for external use, providing additional arguments for advanced configurations.
subroutine compute_family_scaling_expert_c(n_genes, n_families, distances, gene_to_fam, dscale, &
                                           loess_x, loess_y, indices_used, tmp_perm, tmp_stack_left, tmp_stack_right, &
                                           tmp_iv, liv, tmp_wv, lv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat, &
                                           span, degree, mode, n_iters, low_sd_cutoff, excluded_low_sd, tmp_means_aux, ierr) bind(C, name="compute_family_scaling_expert_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_get_outliers, only: compute_family_scaling
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Total number of genes
    integer(c_int), intent(in), target :: n_families
        !! Total number of families
    real(c_double), intent(in), target :: distances(n_genes)
        !! Array of Euclidean distances for each gene
    integer(c_int), intent(in), target :: gene_to_fam(n_genes)
        !! Mapping of each gene to its family (1-based)
    real(c_double), intent(out), target :: dscale(n_families)
        !! Output: array of scaling factors per family
    real(c_double), intent(inout), target :: loess_x(n_families)
        !! Reference x-coordinates for LOESS
    real(c_double), intent(inout), target :: loess_y(n_families)
        !! Reference y-coordinates for LOESS
    integer(c_int), intent(inout), target :: indices_used(n_families)
        !! Indices of reference points used for smoothing
    integer(c_int), intent(inout), target :: tmp_perm(n_genes)
        !! Temporary array for permutation
    integer(c_int), intent(inout), target :: tmp_stack_left(n_genes)
        !! Temporary array for left stack
    integer(c_int), intent(inout), target :: tmp_stack_right(n_genes)
        !! Temporary array for right stack
    integer(c_int), intent(inout), target :: tmp_iv(liv)
        !! Integer workspace array for LOESS
    integer(c_int), intent(in), target :: liv
        !! Length of integer workspace array
    integer(c_int), intent(in), target :: lv
        !! Length of real workspace array
    real(c_double), intent(inout), target :: tmp_wv(lv)
        !! Real workspace array for LOESS
    real(c_double), intent(inout), target :: tmp_diagl(n_genes)
        !! Diagonal elements for LOESS
    real(c_double), intent(inout), target :: tmp_w_init(n_genes)
        !! Initial weights for LOESS
    real(c_double), intent(inout), target :: tmp_z_mat(n_genes, 1)
        !! Z matrix for LOESS
    real(c_double), intent(inout), target :: tmp_rw(n_genes)
        !! Residual weights for LOESS
    real(c_double), intent(inout), target :: tmp_ww(n_genes)
        !! Working weights for LOESS
    real(c_double), intent(inout), target :: tmp_res(n_genes)
        !! Residuals for LOESS
    integer(c_int), intent(inout), target :: tmp_pi(n_genes)
        !! Pi values for LOESS
    real(c_double), intent(inout), target :: tmp_yhat(n_genes)
        !! Temporary array for predicted values
    real(c_double), intent(in), target :: span
        !! Span parameter for LOESS
    integer(c_int), intent(in), target :: degree
        !! Degree of polynomial for LOESS
    integer(c_int), intent(in), target :: mode
        !! Mode for LOESS
    integer(c_int), intent(in), target :: n_iters
        !! Number of iterations for LOESS
    real(c_double), intent(out), target :: low_sd_cutoff
        !! cutoff used to filter families with low std
    integer(c_int), intent(out), target :: excluded_low_sd(n_families)
        !! Mask to save those families that have low sd
    real(c_double), intent(inout), target :: tmp_means_aux(n_families)
        !! Work array for saving raw means
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(distances)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(dscale)
    M_CHECK_NON_NULL(loess_x)
    M_CHECK_NON_NULL(loess_y)
    M_CHECK_NON_NULL(indices_used)
    M_CHECK_NON_NULL(tmp_perm)
    M_CHECK_NON_NULL(tmp_stack_left)
    M_CHECK_NON_NULL(tmp_stack_right)
    M_CHECK_NON_NULL(tmp_iv)
    M_CHECK_NON_NULL(liv)
    M_CHECK_NON_NULL(lv)
    M_CHECK_NON_NULL(tmp_wv)
    M_CHECK_NON_NULL(tmp_diagl)
    M_CHECK_NON_NULL(tmp_w_init)
    M_CHECK_NON_NULL(tmp_z_mat)
    M_CHECK_NON_NULL(tmp_rw)
    M_CHECK_NON_NULL(tmp_ww)
    M_CHECK_NON_NULL(tmp_res)
    M_CHECK_NON_NULL(tmp_pi)
    M_CHECK_NON_NULL(tmp_yhat)
    M_CHECK_NON_NULL(span)
    M_CHECK_NON_NULL(degree)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(n_iters)
    M_CHECK_NON_NULL(low_sd_cutoff)
    M_CHECK_NON_NULL(excluded_low_sd)
    M_CHECK_NON_NULL(tmp_means_aux)

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, tmp_perm, tmp_stack_left, tmp_stack_right, &
                                tmp_iv, liv, tmp_wv, lv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat, &
                                span, degree, mode, n_iters, low_sd_cutoff, excluded_low_sd, tmp_means_aux, ierr)
end subroutine compute_family_scaling_expert_c
