#include <src/macros.h>

!> Gene outliers, from how far each gene sits from its family's centroid.
!|
!| The pipeline is three steps, each callable on its own. `compute_rdi` turns raw distances into
!| a relative distance index, scaled per family so families of different spread are comparable.
!| `compute_family_scaling` fits that scaling with LOESS against family size. `identify_outliers`
!| applies the threshold and reports which genes exceed it.
!|
!| `detect_outliers` runs all three in one call, and is the entry point to reach for first.
module tox_get_outliers_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use f42_math_impl, only: logx_helper, above, is_close
    use f42_sort_impl, only: sort_array, init_perm
    use f42_stats_impl, only: calc_percentile_impl, compute_scaled_distance_quantile_impl
    use tox_errors, only: ERR_INVALID_INPUT, ERR_ALLOC_FAIL, set_ok, set_err, set_err_once, is_err, &
                          validate_all_in_range_real
    use tox_loess_impl, only: tox_loess_required_workspace, EPS_LOESS, loess_evaluation, loess_fit_plain_impl, loess_fit_robust_impl
    M_IMPLICIT_NONE

#define CM_FAMILY_SPAN_DEFAULT 0.7_real64
#define CM_FAMILY_DEGREE_DEFAULT 2_int32
#define CM_FAMILY_MODE_DEFAULT 1_int32
#define CM_FAMILY_N_ITERS_DEFAULT 3_int32
#define CM_OUTLIER_PERCENTILE_DEFAULT 0.95_real64

contains

    !> summary: Compute family scaling factors (dscale) to normalize distances
    !| AUTHOR_VIVIAN_BASS
    !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
    subroutine compute_family_scaling_impl( &
        n_genes, n_families, distances, gene_to_fam, dscale, &
        loess_x, loess_y, indices_used, tmp_perm, tmp_stack_left, tmp_stack_right, &
        tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, tmp_diagl, tmp_weights, tmp_eval_points, tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, tmp_fitted_values, &
        span, degree, mode, n_iters, low_sd_cutoff, excluded_low_sd, tmp_means_aux, ierr)

        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
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
        integer(int32), intent(in)    :: int_workspace_size
            !! Length of integer workspace.
            !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |-----------------------|------------|
            !! | n_dim                 | 1_int32    |
            !! | max_neighborhood_size | n_families |
            !! | save_factorization    | .false.    |
        integer(int32), intent(out) :: tmp_int_workspace(int_workspace_size)
            !! Integer workspace array
        integer(int32), intent(in)    :: real_workspace_size
            !! Length of real workspace.
            !! DM_OUTPUT_FROM(real_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |-----------------------|------------|
            !! | n_dim                 | 1_int32    |
            !! | max_neighborhood_size | n_families |
            !! | save_factorization    | .false.    |
        real(real64), intent(out) :: tmp_real_workspace(real_workspace_size)
            !! Real workspace array

        real(real64), intent(inout) :: tmp_diagl(n_families)
            !! Diagonal elements of the weight matrix
        real(real64), intent(out) :: tmp_weights(n_families)
            !! Initial weights for LOESS
        real(real64), intent(inout) :: tmp_eval_points(n_families, 1)
            !! Z matrix for LOESS fitting
        real(real64), intent(out) :: tmp_robust_weights(n_families)
            !! Residuals for robust LOESS fitting
        real(real64), dimension(n_families), intent(out), target :: tmp_combined_weights
            !! Working weights array
        real(real64), dimension(n_families), intent(out), target :: tmp_residuals
            !! Residuals array
        integer(int32), intent(out)  :: tmp_permutation_indices(n_families)
            !! Permutation indices for robust LOESS fitting
        real(real64), intent(out)    :: tmp_fitted_values(n_families)
            !! Output array for LOESS predictions

        real(real64), intent(in), optional     :: span
            !! Span parameter for LOESS smoothing, passed straight to
            !! [[tox_loess_impl(module):loess_fit_plain_impl(subroutine)]], so it is held to that
            !! procedure's own range rather than to the NaN tolerance the distance data carries.
            !! DM_DEFAULT(CM_FAMILY_SPAN_DEFAULT)
            !! DM_MIN(EPS_LOESS)
            !! DM_MAX(1.0_real64)
        integer(int32), intent(in), optional   :: degree
            !! Degree of the LOESS polynomial
            !! DM_DEFAULT(CM_FAMILY_DEGREE_DEFAULT)
        integer(int32), intent(in), optional   :: mode
            !! Mode for LOESS fitting
            !! DM_DEFAULT(CM_FAMILY_MODE_DEFAULT)
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Plain LOESS fitting | [[tox_loess_impl(module):MODE_PLAIN(variable)]] |
            !! | Robust LOESS fitting | [[tox_loess_impl(module):MODE_ROBUST(variable)]] |
        integer(int32), intent(in), optional   :: n_iters
            !! Number of iterations for robust LOESS fitting
            !! DM_DEFAULT(CM_FAMILY_N_ITERS_DEFAULT)
        real(real64), intent(out) :: low_sd_cutoff
            !! cutoff used to filter families with low std
        integer(int32), intent(out)  :: ierr
            !! Error code

        ! Local variables
        integer(int32) :: i_gene, i_family, i_valid, family_idx, n_in_family, n_valid, k
        real(real64)   :: stddev_dist, mean_dist, sumsq, dist_val
        real(real64) :: xmin, xmax, eps_mean, eps_sd, std_median
        real(real64) :: actual_span
        integer(int32) :: actual_degree, actual_mode, actual_n_iters

        M_DEFAULT_VAL(span, actual_span, CM_FAMILY_SPAN_DEFAULT)
        M_DEFAULT_VAL(degree, actual_degree, CM_FAMILY_DEGREE_DEFAULT)
        M_DEFAULT_VAL(mode, actual_mode, CM_FAMILY_MODE_DEFAULT)
        M_DEFAULT_VAL(n_iters, actual_n_iters, CM_FAMILY_N_ITERS_DEFAULT)

        ! Initialize error code and output arrays
        call set_ok(ierr)
        dscale  = 0.0_real64
        loess_x = 0.0_real64
        loess_y = 0.0_real64
        n_valid = 0

        ! Validate family indices; the -1 sentinel output on failure is part of the contract, so this
        ! stays in the implementation rather than becoming a wrapper range check.
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

        tmp_weights = 0.0_real64
        tmp_robust_weights = 0.0_real64
        tmp_permutation_indices = 0

        do i_gene = 1, n_genes
            family_idx = gene_to_fam(i_gene)
            dist_val = abs(distances(i_gene))

            tmp_permutation_indices(family_idx) = tmp_permutation_indices(family_idx) + 1
            tmp_weights(family_idx) = tmp_weights(family_idx) + dist_val
            tmp_robust_weights(family_idx) = tmp_robust_weights(family_idx) + (dist_val**2)
        end do

        n_valid = 0
        do i_family = 1, n_families
            n_in_family = tmp_permutation_indices(i_family)

            if (n_in_family <= 1) cycle

            n_valid = n_valid + 1

            mean_dist = tmp_weights(i_family)/real(n_in_family, real64)

            ! Var = (SumSq - (Sum^2)/N) / (N-1)
            sumsq = max(0.0_real64, tmp_robust_weights(i_family) - (tmp_weights(i_family)**2/real(n_in_family, real64)))
            stddev_dist = sqrt(sumsq/real(n_in_family - 1, real64))

            loess_x(n_valid) = mean_dist
            tmp_means_aux(i_family) = mean_dist
            loess_y(n_valid) = stddev_dist
            indices_used(n_valid) = i_family
        end do

        tmp_weights = 0.0_real64
        tmp_robust_weights = 0.0_real64
        tmp_permutation_indices = 0

        if (n_valid <= 1) then
            low_sd_cutoff = 0.0_real64
            return
        end if

        call init_perm(tmp_perm)

        call sort_array(loess_x(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        ! Use the 5th percentile of the family means as a data-driven pseudo-count instead of a fixed
        ! constant, so log2(mean + eps_mean) below stays well-scaled across datasets with very
        ! different absolute expression ranges.
        call calc_percentile_impl(loess_x, n_valid, tmp_perm, 0.05_real64, eps_mean)

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
        !
        ! Validate the log2 arguments up front (sequentially, as `ierr` is a shared scalar), exactly as
        ! `logx` would: `loess_x + eps_mean > 0` and `loess_y + eps_sd > 0`, rejecting NaN/Inf. This
        ! does not assume the caller's `distances` are well-formed -- a NaN/Inf distance propagates into
        ! `loess_x`/`loess_y` and is caught here. These same family means feed the second log2 loop
        ! below via `tmp_means_aux` (identical values), so validating them here covers both loops.
        ! With the arguments guaranteed valid, the transform runs as a race-free `do concurrent`
        ! calling the non-validating `logx_helper` -- no per-iteration write to the shared `ierr`.
        call validate_all_in_range_real(loess_x(1:n_valid) + eps_mean, n_valid, ierr, min=above(0.0_real64))
        call validate_all_in_range_real(loess_y(1:n_valid) + eps_sd, n_valid, ierr, min=above(0.0_real64))
        if (is_err(ierr)) return

        do concurrent(i_valid=1:n_valid) shared(loess_x, loess_y, eps_mean, eps_sd)
            call logx_helper(loess_x(i_valid) + eps_mean, 2.0_real64, loess_x(i_valid))
            call logx_helper(loess_y(i_valid) + eps_sd, 2.0_real64, loess_y(i_valid))
        end do

        call init_perm(tmp_perm)
        call sort_array(loess_y(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        ! Bottom 1% of (log2) stddevs is treated as "too flat to trust": these families would otherwise
        ! anchor the mean-vs-stddev LOESS curve with near-degenerate (close to zero-variance) points.
        call calc_percentile_impl(loess_y, n_valid, tmp_perm, 0.01_real64, low_sd_cutoff)

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
        ! we use the global median for all families. Closeness is judged on LOESS's own terms
        ! -- EPS_LOESS, the smoothing floor -- since a spread finer than that is no spread to
        ! fit a curve through, however different the two ends are as floating point numbers.
        if (is_close(xmax, xmin, EPS_LOESS)) then
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
        tmp_weights(1:n_valid) = 1.0_real64
        tmp_eval_points(1:n_valid, 1) = loess_x(1:n_valid)

        if (actual_mode == 0) then
            ! If you have a plain routine, call it; otherwise keep robust always.
            call loess_fit_plain_impl( &
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_weights(1:n_valid), tmp_eval_points, &
                actual_span, actual_degree, n_valid, .false._c_bool, .false._c_bool, tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, tmp_diagl(1:n_valid), tmp_fitted_values(1:n_valid), ierr)
        else
            call loess_fit_robust_impl( &
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_weights(1:n_valid), tmp_eval_points, &
                actual_span, actual_degree, n_valid, .false._c_bool, .false._c_bool, actual_n_iters, tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, tmp_diagl(1:n_valid), &
                tmp_robust_weights(1:n_valid), tmp_combined_weights(1:n_valid), tmp_residuals(1:n_valid), tmp_permutation_indices(1:n_valid), tmp_fitted_values(1:n_valid), ierr)
        end if

        if (is_err(ierr)) return

        ! The `>= 0` branch feeds `tmp_means_aux(i_family) + eps_mean` to log2. Every such family had
        ! `n_in_family > 1` and therefore contributed exactly this mean to `loess_x`, whose
        ! `+ eps_mean` argument was already validated (> 0, finite) before the first log2 loop above;
        ! skipped families keep the initial `tmp_means_aux = -1` and take the `else` branch. So the
        ! argument here is guaranteed valid by that earlier validation -- not by assumption -- and the
        ! non-validating `logx_helper` runs in a race-free `do concurrent` with no shared `ierr`
        ! (each iteration writes only its own `tmp_eval_points` row).
        do concurrent(i_family=1:n_families) shared(tmp_means_aux, tmp_eval_points, eps_mean, xmin, xmax)
            if (tmp_means_aux(i_family) >= 0.0_real64) then
                call logx_helper(tmp_means_aux(i_family) + eps_mean, 2.0_real64, tmp_eval_points(i_family, 1))
                if (tmp_eval_points(i_family, 1) < xmin) tmp_eval_points(i_family, 1) = xmin
                if (tmp_eval_points(i_family, 1) > xmax) tmp_eval_points(i_family, 1) = xmax
            else
                tmp_eval_points(i_family, 1) = xmin
            end if
        end do

        call loess_evaluation(tmp_int_workspace, int_workspace_size, real_workspace_size, tmp_real_workspace, n_families, tmp_eval_points, tmp_fitted_values(1:n_families))

        if (is_err(ierr)) return

        ! ------------------------------------------------------------
        ! Map compact tmp_residuals back to full dscale
        ! ------------------------------------------------------------

        do concurrent (i_family = 1:n_families) shared(tmp_means_aux, dscale, tmp_fitted_values, eps_sd)
            if (tmp_means_aux(i_family) < 0.0_real64) then
                dscale(i_family) = 0.0_real64
            else
                dscale(i_family) = max(2.0_real64**tmp_fitted_values(i_family) - eps_sd, 0.0_real64)
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

    end subroutine compute_family_scaling_impl

    !> summary: Compute the hybrid RDI (Relative Distance Index) for each gene
    !| AUTHOR_VIVIAN_BASS
    !| RDI = Euclidean distance / family scaling factor
    pure subroutine compute_rdi_impl(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, &
                                tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene to its centroid
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! Gene-to-family mapping (1-based indexing)
        real(real64), intent(in) :: dscale(:)
            !! Array of scaling factors for each family
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
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
        ! rather than by the caller, since this routine is exported and an binding
        ! language hands it freshly allocated memory
        call init_perm(perm)

        ! Sort RDI values using the tox_sorting module
        call sort_array(sorted_rdi, perm, tmp_stack_left, tmp_stack_right)

    end subroutine compute_rdi_impl

    !> summary: Identify gene outliers based on the top percentile of RDI values
    !| AUTHOR_VIVIAN_BASS
    !| Expects sorted_rdi to be filtered (no negative values) and perm should be sorted in ascending order before calling.
    !| If sorted_rdi contains negatives or perm is not sorted, tmp_results may be invalid.
    pure subroutine identify_outliers_impl(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, quantile, percentile)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        real(real64), intent(in) :: rdi(n_genes)
            !! Array of RDI values for each gene
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), intent(in) :: sorted_rdi(n_genes)
            !! Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        integer(int32), intent(in) :: perm(n_genes)
            !! Permutation array with sorted indices
        logical(c_bool), intent(out) :: is_outlier(n_genes)
            !! Output boolean array indicating outliers
        real(real64), intent(out) :: threshold
            !! Output threshold value used for detection
        real(real64), intent(in), optional :: percentile
            !! Percentile threshold as a fraction in [0,1] (top 5% for the default).
            !! DM_DEFAULT(CM_OUTLIER_PERCENTILE_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), intent(out) :: quantile(n_genes)
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.

        integer(int32) :: i, idx
        real(real64) :: perc_pos, percentile_val

        M_DEFAULT_VAL(percentile, percentile_val, CM_OUTLIER_PERCENTILE_DEFAULT)

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
        ! interpolated one. `percentile_val` is a fraction in [0,1].
        perc_pos = n_genes*percentile_val
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

        call compute_scaled_distance_quantile_impl(n_genes, rdi, sorted_rdi, perm, quantile, 1.0_real64)

    end subroutine identify_outliers_impl

    !> summary: Main routine to detect outliers using RDI and LOESS-based scaling
    !| AUTHOR_VIVIAN_BASS
    !| Orchestrates the full pipeline: per-family scaling via
    !| [[tox_get_outliers_impl(module):compute_family_scaling_impl(subroutine)]], the RDI per gene via
    !| [[tox_get_outliers_impl(module):compute_rdi_impl(subroutine)]], then flags outliers via
    !| [[tox_get_outliers_impl(module):identify_outliers_impl(subroutine)]].
    subroutine detect_outliers_impl(n_genes, n_families, distances, gene_to_fam, &
                               tmp_perm, tmp_stack_left, tmp_stack_right, &
                               tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, &
                               tmp_diagl, tmp_weights, tmp_eval_points, tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, tmp_fitted_values, tmp_means_aux, &
                               tmp_dscale, tmp_excluded_low_sd, tmp_low_sd_cutoff, tmp_rdi, tmp_sorted_rdi, tmp_threshold, &
                               is_outlier, loess_x, loess_y, loess_n, quantile, ierr, percentile)
        integer(int32), intent(in) :: n_genes
            !! Total number of genes
        integer(int32), intent(in) :: n_families
            !! Total number of gene families
        real(real64), intent(in) :: distances(n_genes)
            !! Array of Euclidean distances for each gene to its centroid
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        integer(int32), intent(in) :: gene_to_fam(n_genes)
            !! Gene-to-family mapping (1-based indexing)

        ! Shared sort scratch (n_genes), reused by the family-scaling and RDI phases
        integer(int32), intent(out) :: tmp_perm(n_genes)
            !! Permutation array for sorting
        integer(int32), intent(out) :: tmp_stack_left(n_genes)
            !! Stack array for left indices during sorting
        integer(int32), intent(out) :: tmp_stack_right(n_genes)
            !! Stack array for right indices during sorting

        ! LOESS workspace for compute_family_scaling
        integer(int32), intent(in)    :: int_workspace_size
            !! Length of integer workspace.
            !! DM_OUTPUT_FROM(int_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |-----------------------|------------|
            !! | n_dim                 | 1_int32    |
            !! | max_neighborhood_size | n_families |
            !! | save_factorization    | .false.    |
        integer(int32), intent(out) :: tmp_int_workspace(int_workspace_size)
            !! Integer workspace array
        integer(int32), intent(in)    :: real_workspace_size
            !! Length of real workspace.
            !! DM_OUTPUT_FROM(real_workspace_size, tox_loess_required_workspace, tox_loess_impl, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |-----------------------|------------|
            !! | n_dim                 | 1_int32    |
            !! | max_neighborhood_size | n_families |
            !! | save_factorization    | .false.    |
        real(real64), intent(out) :: tmp_real_workspace(real_workspace_size)
            !! Real workspace array
        real(real64), intent(out) :: tmp_diagl(n_families)
            !! Diagonal elements of the weight matrix
        real(real64), intent(out) :: tmp_weights(n_families)
            !! Initial weights for LOESS
        real(real64), intent(out) :: tmp_eval_points(n_families, 1)
            !! Z matrix for LOESS fitting
        real(real64), intent(out) :: tmp_robust_weights(n_families)
            !! Residuals for robust LOESS fitting
        real(real64), intent(out) :: tmp_combined_weights(n_families)
            !! Working weights array
        real(real64), intent(out) :: tmp_residuals(n_families)
            !! Residuals array
        integer(int32), intent(out) :: tmp_permutation_indices(n_families)
            !! Permutation indices for robust LOESS fitting
        real(real64), intent(out) :: tmp_fitted_values(n_families)
            !! Output array for LOESS predictions
        real(real64), intent(out) :: tmp_means_aux(n_families)
            !! Work array for saving raw means
        real(real64), intent(out) :: tmp_dscale(n_families)
            !! Per-family scaling factors (intermediate, consumed by the RDI step)
        integer(int32), intent(out) :: tmp_excluded_low_sd(n_families)
            !! Low-sd family mask (intermediate, discarded)
        real(real64), intent(out) :: tmp_low_sd_cutoff
            !! Low-sd cutoff (intermediate, discarded)
        real(real64), intent(out) :: tmp_rdi(n_genes)
            !! RDI per gene (intermediate, consumed by the outlier step)
        real(real64), intent(out) :: tmp_sorted_rdi(n_genes)
            !! Sorted RDI (intermediate, consumed by the outlier step)
        real(real64), intent(out) :: tmp_threshold
            !! Detection threshold (intermediate, discarded)

        logical(c_bool), intent(out) :: is_outlier(n_genes)
            !! Output boolean array indicating outliers
        real(real64), intent(out) :: loess_x(n_families)
            !! Reference x-coordinates.
        real(real64), intent(out) :: loess_y(n_families)
            !! Reference y-coordinates (length n_total).
        integer(int32), intent(out) :: loess_n(n_families)
            !! Indices of reference points used for smoothing.
        real(real64), intent(out) :: quantile(n_genes)
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), intent(in), optional :: percentile
            !! Percentile threshold as a fraction in [0,1] for outlier detection.
            !! DM_DEFAULT(CM_OUTLIER_PERCENTILE_DEFAULT)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)

        call set_ok(ierr)

        call compute_family_scaling_impl( &
            n_genes, n_families, distances, gene_to_fam, tmp_dscale, &
            loess_x, loess_y, loess_n, tmp_perm, tmp_stack_left, tmp_stack_right, &
            tmp_int_workspace, int_workspace_size, tmp_real_workspace, real_workspace_size, tmp_diagl, tmp_weights, tmp_eval_points, tmp_robust_weights, tmp_combined_weights, tmp_residuals, tmp_permutation_indices, tmp_fitted_values, &
            span=CM_FAMILY_SPAN_DEFAULT, degree=CM_FAMILY_DEGREE_DEFAULT, mode=CM_FAMILY_MODE_DEFAULT, n_iters=CM_FAMILY_N_ITERS_DEFAULT, &
            low_sd_cutoff=tmp_low_sd_cutoff, excluded_low_sd=tmp_excluded_low_sd, tmp_means_aux=tmp_means_aux, ierr=ierr)
        if (is_err(ierr)) return

        call compute_rdi_impl(n_genes, distances, gene_to_fam, tmp_dscale, tmp_rdi, tmp_sorted_rdi, tmp_perm, tmp_stack_left, tmp_stack_right)

        call identify_outliers_impl(n_genes, tmp_rdi, tmp_sorted_rdi, tmp_perm, is_outlier, tmp_threshold, quantile, percentile)
    end subroutine detect_outliers_impl
end module tox_get_outliers_impl
