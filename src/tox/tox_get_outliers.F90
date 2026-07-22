#include <src/macros.h>

!> Module to identify gene outliers based on their distances to family centroids.
module tox_get_outliers
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use f42_utils, only: sort_array, calc_percentile, logx, is_close, compute_empirical_p_values, init_perm
    use tox_errors, only: ERR_INVALID_INPUT, ERR_ALLOC_FAIL, set_ok, set_err, set_err_once, is_err
    use tox_loess, only: tox_loess_required_workspace, loess_fit_robust, loess_fit_plain, EPS_LOESS, loess_evaluation
    M_IMPLICIT_NONE

contains

    !> M_EXPORT_C
    !| summary: Compute family scaling factors (dscale) to normalize distances (expert entry point)
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
    subroutine compute_family_scaling( &
        n_genes, n_families, distances, gene_to_fam, dscale, &
        loess_x, loess_y, indices_used, tmp_perm, tmp_stack_left, tmp_stack_right, &
        tmp_iv, liv, tmp_wv, lv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat, &
        span, degree, mode, n_iters, low_sd_cutoff, excluded_low_sd, tmp_means_aux, ierr)

        use, intrinsic :: iso_fortran_env, only: real64, int32
        M_IMPLICIT_NONE

        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! Mapping of each gene to its family (1-based)

        real(real64), intent(out) :: dscale(n_families)
            !! Array of scaling factors per family (output)

        ! Buffers (reused)
        real(real64), intent(out) :: loess_x(n_families)
            !! Reference x-coordinates for LOESS smoothing
        real(real64), intent(out) :: loess_y(n_families)
            !! Reference y-coordinates for LOESS smoothing
        integer(int32), intent(out) :: indices_used(n_families)
            !! Indices of reference points used for smoothing
        integer(int32), intent(out) :: tmp_perm(n_genes)
            !! Permutation array for sorting gene distances
        integer(int32), intent(out) :: tmp_stack_left(n_genes)
            !! Stack array for left indices during sorting
        integer(int32), intent(out) :: tmp_stack_right(n_genes)
            !! Stack array for right indices during sorting
        real(real64), intent(out) :: tmp_means_aux(n_families)
            !! Work array for saving raw means
        integer(int32), intent(out) :: excluded_low_sd(n_families)
            !! Mask to save those families that have low sd

        ! LOESS workspace
        integer(int32), intent(in)    :: liv
            !! Length of integer workspace.
            !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |-----------------------|------------|
            !! | n_dim                 | 1_int32    |
            !! | max_neighborhood_size | n_families |
            !! | save_factorization    | .false.    |
        integer(int32), intent(out) :: tmp_iv(liv)
            !! Integer workspace array
        integer(int32), intent(in)    :: lv
            !! Length of real workspace.
            !! DM_OUTPUT_FROM(real_workspace_size, tox_loess_required_workspace, tox_loess, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |-----------------------|------------|
            !! | n_dim                 | 1_int32    |
            !! | max_neighborhood_size | n_families |
            !! | save_factorization    | .false.    |
        real(real64), intent(out) :: tmp_wv(lv)
            !! Real workspace array

        real(real64), intent(inout) :: tmp_diagl(n_families)
            !! Diagonal elements of the weight matrix
        real(real64), intent(out) :: tmp_w_init(n_families)
            !! Initial weights for LOESS
        real(real64), intent(inout) :: tmp_z_mat(n_families, 1)
            !! Z matrix for LOESS fitting
        real(real64), intent(out) :: tmp_rw(n_families)
            !! Residuals for robust LOESS fitting
        real(real64), dimension(n_families), intent(out), target :: tmp_ww
            !! Working weights array
        real(real64), dimension(n_families), intent(out), target :: tmp_res
            !! Residuals array
        integer(int32), intent(out)  :: tmp_pi(n_families)
            !! Permutation indices for robust LOESS fitting
        real(real64), intent(out)    :: tmp_yhat(n_families)
            !! Output array for LOESS predictions

        real(real64), intent(in)     :: span
            !! Span parameter for LOESS smoothing
        integer(int32), intent(in)   :: degree
            !! Degree of the LOESS polynomial
        integer(int32), intent(in)   :: mode
            !! Mode for LOESS fitting
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Plain LOESS fitting | [[tox_loess(module):MODE_PLAIN(variable)]] |
            !! | Robust LOESS fitting | [[tox_loess(module):MODE_ROBUST(variable)]] |
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

        ! netlib reads iv and wv as workspace it has already seeded, so they are
        ! initialised here rather than by the caller: this routine is exported, and an
        ! interfacing language hands it freshly allocated -- uninitialised -- memory.
        tmp_iv = 1_int32
        tmp_wv = 0.0_real64

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

        call init_perm(tmp_perm)

        call sort_array(loess_x(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        ! Use the 5th percentile of the family means as a data-driven pseudo-count instead of a fixed
        ! constant, so log2(mean + eps_mean) below stays well-scaled across datasets with very
        ! different absolute expression ranges.
        call calc_percentile(loess_x(1:n_valid), tmp_perm(1:n_valid), 5.0_real64, eps_mean, ierr)
        if (is_err(ierr)) return

        eps_mean = max(eps_mean, EPS_LOESS)

        call init_perm(tmp_perm)
        call sort_array(loess_y(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        if (mod(n_valid, 2) == 0) then
            std_median = 0.5_real64*( &
                         loess_y(tmp_perm(n_valid/2)) + &
                         loess_y(tmp_perm(n_valid/2 + 1)))
        else
            std_median = loess_y(tmp_perm((n_valid + 1)/2))
        end if

        ! Same idea as eps_mean above, scaled relative to the median stddev (1e-13 is a near-machine-
        ! precision fraction) so the pseudo-count added to loess_y before log2 stays negligible except
        ! when it is needed to keep near-zero stddevs away from log2(0).
        eps_sd = max(1.0e-13_real64*std_median, EPS_LOESS)

        ! Fit the mean-vs-stddev trend in log2 space: intra-family distance stddev scales roughly
        ! multiplicatively with the mean distance, so a log/log relationship is closer to linear
        ! (homoscedastic) than the raw scale, which is what LOESS assumes.
        ! NOTE: kept as a plain sequential loop (not `do concurrent`) because `ierr` is a shared
        ! scalar conditionally written on the (rare/exceptional) error path -- writing it from
        ! concurrent iterations would be an unsynchronized data race.
        do i_valid = 1, n_valid
            call logx(loess_x(i_valid) + eps_mean, 2.0_real64, loess_x(i_valid), tmp_ierr)
            if (is_err(tmp_ierr)) ierr = tmp_ierr
            call logx(loess_y(i_valid) + eps_sd, 2.0_real64, loess_y(i_valid), tmp_ierr)
            if (is_err(tmp_ierr)) ierr = tmp_ierr
        end do

        if (is_err(ierr)) return

        call init_perm(tmp_perm)
        call sort_array(loess_y(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        ! Bottom 1% of (log2) stddevs is treated as "too flat to trust": these families would otherwise
        ! anchor the mean-vs-stddev LOESS curve with near-degenerate (close to zero-variance) points.
        call calc_percentile(loess_y(1:n_valid), tmp_perm(1:n_valid), 1.0_real64, low_sd_cutoff, ierr)

        if (is_err(ierr)) return

        ! Compact loess_x/loess_y/indices_used in place, keeping only families at or above the low-sd
        ! cutoff, so the subsequent global LOESS fit below is trained on the retained (k <= n_valid)
        ! subset only. excluded_low_sd flags which families were dropped from the fit (they still get a
        ! dscale prediction from the resulting curve further below).
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
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_w_init(1:n_valid), tmp_z_mat, &
                span, degree, n_valid, .false., .false., tmp_iv, liv, tmp_wv, lv, tmp_diagl(1:n_valid), tmp_yhat(1:n_valid), ierr)
        else
            call loess_fit_robust( &
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_w_init(1:n_valid), tmp_z_mat, &
                span, degree, n_valid, .false., .false., n_iters, tmp_iv, liv, tmp_wv, lv, tmp_diagl(1:n_valid), &
                tmp_rw(1:n_valid), tmp_ww(1:n_valid), tmp_res(1:n_valid), tmp_pi(1:n_valid), tmp_yhat(1:n_valid), ierr)
        end if

        if (is_err(ierr)) return

        ! NOTE: kept as a plain sequential loop (not `do concurrent`) because `ierr` is a shared
        ! scalar conditionally written on the (rare/exceptional) error path -- writing it from
        ! concurrent iterations would be an unsynchronized data race.
        do i_family = 1, n_families
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

        call loess_evaluation(tmp_iv, liv, lv, tmp_wv, n_families, tmp_z_mat, tmp_yhat(1:n_families))

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

    !> M_EXPORT_C
    !| summary: Allocates internal LOESS work arrays and computes per-family scaling factors
    !| AUTHOR_FRANZ_ERIC_SILL
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
        real(real64), intent(out) :: loess_x(n_families)
            !! Reference x-coordinates.
        real(real64), intent(out) :: loess_y(n_families)
            !! Reference y-coordinates (length n_total).
        integer(int32), intent(out) :: indices_used(n_families)
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
        integer(int32) :: liv, lv
        integer(int32), allocatable :: tmp_iv(:), tmp_pi(:)
        real(real64), allocatable :: tmp_wv(:), tmp_diagl(:), tmp_w_init(:), tmp_z_mat(:, :), tmp_rw(:), tmp_ww(:), tmp_res(:), tmp_yhat(:)

        ! LOESS params (defaults)
        real(real64), parameter :: span = 0.7_real64
        integer(int32), parameter :: degree = 2_int32
        integer(int32), parameter :: mode = 1_int32
        integer(int32), parameter :: n_iters = 3_int32
        logical, parameter :: setlf = .false.

        call set_ok(ierr)

        ! Workspace sizes
        call tox_loess_required_workspace(1_int32, n_families, liv, lv, setlf)

        M_ALLOCATE(tmp_iv(liv))
        M_ALLOCATE(tmp_wv(lv))

        ! For robust we also need arrays sized to n_valid (<= n_families).
        M_ALLOCATE(tmp_diagl(n_families))
        M_ALLOCATE(tmp_w_init(n_families))
        M_ALLOCATE(tmp_z_mat(n_families, 1))
        M_ALLOCATE(tmp_rw(n_families))
        M_ALLOCATE(tmp_ww(n_families))
        M_ALLOCATE(tmp_res(n_families))
        M_ALLOCATE(tmp_pi(n_families))
        M_ALLOCATE(tmp_yhat(n_families))
        M_ALLOCATE(tmp_perm(n_genes))
        M_ALLOCATE(tmp_stack_left(n_genes))
        M_ALLOCATE(tmp_stack_right(n_genes))
        M_ALLOCATE(excluded_low_sd(n_families))
        M_ALLOCATE(tmp_means_aux(n_families))

        call compute_family_scaling( &
            n_genes, n_families, distances, gene_to_fam, dscale, &
            loess_x, loess_y, indices_used, tmp_perm, tmp_stack_left, tmp_stack_right, &
            tmp_iv, liv, tmp_wv, lv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat, &
            span, degree, mode, n_iters, low_sd_cutoff, excluded_low_sd, tmp_means_aux, ierr)

        deallocate (tmp_iv, tmp_wv, tmp_diagl, tmp_w_init, tmp_z_mat, tmp_rw, tmp_ww, tmp_res, tmp_pi, tmp_yhat)

    end subroutine compute_family_scaling_alloc

    !> M_EXPORT_C
    !| summary: Compute the hybrid RDI (Relative Distance Index) for each gene
    !| AUTHOR_FRANZ_ERIC_SILL
    !| RDI = Euclidean distance / family scaling factor
    pure subroutine compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, &
                                tmp_stack_left, tmp_stack_right)
        M_IMPLICIT_NONE
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
        real(real64), intent(out) :: sorted_rdi(n_genes)
            !! Work array for sorting (dimension n_genes)
        integer(int32), intent(out) :: perm(n_genes)
            !! Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
        integer(int32), intent(out) :: tmp_stack_left(n_genes)
            !! Stack array for sorting (dimension n_genes)
        integer(int32), intent(out) :: tmp_stack_right(n_genes)
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

        ! the quicksort permutes `perm`, so it has to start as the identity: seeded here
        ! rather than by the caller, since this routine is exported and an interfacing
        ! language hands it freshly allocated memory
        call init_perm(perm)

        ! Sort RDI values using the tox_sorting module
        call sort_array(sorted_rdi, perm, tmp_stack_left, tmp_stack_right)

    end subroutine compute_rdi

    !> M_EXPORT_C
    !| summary: Identify gene outliers based on the top percentile of RDI values
    !| AUTHOR_FRANZ_ERIC_SILL
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
            !! Percentile threshold (top 5% for the default).
            !! DM_DEFAULT(95.0_real64)
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

        ! Guard against n_genes < 1: idx would otherwise be clamped to 0 below (idx<1 -> 1,
        ! then idx>n_genes==0 -> 0), causing an out-of-bounds access at sorted_rdi(perm(idx)).
        if (n_genes < 1) then
            threshold = 0.0_real64
            return
        end if

        ! Nearest-rank percentile: round the fractional rank up to the next integer index into the
        ! ascending-sorted array, so `threshold` is always an observed RDI value rather than an
        ! interpolated one.
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

    !> M_EXPORT_C
    !| summary: Main routine to detect outliers using RDI and LOESS-based scaling
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Orchestrates the full pipeline: computes per-family scaling factors via
    !| [[tox_get_outliers(module):compute_family_scaling_alloc(subroutine)]], derives the RDI per gene via
    !| [[tox_get_outliers(module):compute_rdi(subroutine)]], then flags outliers via
    !| [[tox_get_outliers(module):identify_outliers(subroutine)]].
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
            !! Percentile threshold for outlier detection.
            !! DM_DEFAULT(95.0_real64)
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
        call init_perm(tmp_perm)

        call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                          loess_x, loess_y, loess_n, ierr)
        if (is_err(ierr)) return
        call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right)
        call identify_outliers(n_genes, rdi, tmp_work_array, tmp_perm, is_outlier, threshold, p_values, percentile_val)
    end subroutine detect_outliers
end module tox_get_outliers
