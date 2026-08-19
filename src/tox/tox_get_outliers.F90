#include <src/macros.h>

!> Module to identify gene outliers based on their distances to family centroids.
module tox_get_outliers
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use f42_utils, only: sort_array, calc_percentile, logx_helper, above, is_close, compute_scaled_distance_quantile, init_perm
    use tox_errors, only: ERR_INVALID_INPUT, ERR_ALLOC_FAIL, set_ok, set_err, set_err_once, is_err, validate_all_in_range_real
    use tox_loess, only: tox_loess_required_workspace, loess_fit_robust, loess_fit_plain, EPS_LOESS, loess_evaluation
    implicit none

contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Compute family scaling factors (dscale) to normalize distances.
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
        integer(int32) :: i_gene, i_family, i_valid, family_idx, n_in_family, n_valid, k
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

        call init_perm(tmp_perm)

        call sort_array(loess_x(1:n_valid), tmp_perm(1:n_valid), tmp_stack_left(1:n_valid), tmp_stack_right(1:n_valid))
        ! Use the 5th percentile of the family means as a data-driven pseudo-count instead of a fixed
        ! constant, so log2(mean + eps_mean) below stays well-scaled across datasets with very
        ! different absolute expression ranges.
        call calc_percentile(loess_x(1:n_valid), tmp_perm(1:n_valid), 0.05_real64, eps_mean, ierr)
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
        call calc_percentile(loess_y(1:n_valid), tmp_perm(1:n_valid), 0.01_real64, low_sd_cutoff, ierr)

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
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_w_init(1:n_valid), tmp_z_mat(1:n_valid, 1:1), &
                span, degree, n_valid, .false., .false., tmp_iv, liv, tmp_wv, lv, tmp_diagl(1:n_valid), tmp_yhat(1:n_valid), ierr)
        else
            call loess_fit_robust( &
                n_valid, loess_x(1:n_valid), loess_y(1:n_valid), tmp_w_init(1:n_valid), tmp_z_mat(1:n_valid, 1:1), &
                span, degree, n_valid, .false., .false., n_iters, tmp_iv, liv, tmp_wv, lv, tmp_diagl(1:n_valid), &
                tmp_rw(1:n_valid), tmp_ww(1:n_valid), tmp_res(1:n_valid), tmp_pi(1:n_valid), tmp_yhat(1:n_valid), ierr)
        end if

        if (is_err(ierr)) return

        ! The `>= 0` branch feeds `tmp_means_aux(i_family) + eps_mean` to log2. Every such family had
        ! `n_in_family > 1` and therefore contributed exactly this mean to `loess_x`, whose
        ! `+ eps_mean` argument was already validated (> 0, finite) before the first log2 loop above;
        ! skipped families keep the initial `tmp_means_aux = -1` and take the `else` branch. So the
        ! argument here is guaranteed valid by that earlier validation -- not by assumption -- and the
        ! non-validating `logx_helper` runs in a race-free `do concurrent` with no shared `ierr`
        ! (each iteration writes only its own `tmp_z_mat` row).
        do concurrent(i_family=1:n_families) shared(tmp_means_aux, tmp_z_mat, eps_mean, xmin, xmax)
            if (tmp_means_aux(i_family) >= 0.0_real64) then
                call logx_helper(tmp_means_aux(i_family) + eps_mean, 2.0_real64, tmp_z_mat(i_family, 1))
                if (tmp_z_mat(i_family, 1) < xmin) tmp_z_mat(i_family, 1) = xmin
                if (tmp_z_mat(i_family, 1) > xmax) tmp_z_mat(i_family, 1) = xmax
            else
                tmp_z_mat(i_family, 1) = xmin
            end if
        end do

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

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper routine that allocates internal arrays and calls [[tox_get_outliers(module):compute_family_scaling(subroutine)]].
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

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Compute the hybrid RDI (Relative Distance Index) for each gene.
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
        call init_perm(perm)
        call sort_array(sorted_rdi, perm, stack_left, stack_right)

    end subroutine compute_rdi

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Identify gene outliers based on the top percentile of RDI values.
    !| Expects sorted_rdi to be filtered (no negative values) and perm should be sorted in ascending order before calling.
    !| If sorted_rdi contains negatives or perm is not sorted, tmp_results may be invalid.
    pure subroutine identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, quantile, percentile)
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
            !! (optional) Percentile threshold as a fraction in [0,1] (default: 0.95 for top 5%)
        real(real64), intent(out) :: quantile(n_genes)
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.

        integer(int32) :: i, idx
        real(real64) :: perc_pos, percentile_val

        ! Set default percentile if not present
        if (present(percentile)) then
            percentile_val = percentile
        else
            percentile_val = 0.95_real64
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

        call compute_scaled_distance_quantile(n_genes, rdi, sorted_rdi, perm, quantile, 1.0_real64)

    end subroutine identify_outliers

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Main routine to detect outliers using RDI and LOESS-based scaling.
    !| Orchestrates the full pipeline: computes per-family scaling factors via
    !| [[tox_get_outliers(module):compute_family_scaling_alloc(subroutine)]], derives the RDI per gene via
    !| [[tox_get_outliers(module):compute_rdi(subroutine)]], then flags outliers via
    !| [[tox_get_outliers(module):identify_outliers(subroutine)]].
    subroutine detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                               tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right, &
                               is_outlier, loess_x, loess_y, loess_n, quantile, ierr, &
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
            !! (optional) Percentile threshold in [0,1] for outlier detection (default: 0.95)
        real(real64), intent(out) :: quantile(n_genes)
            !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
            !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
            !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
            !! upper-tail quantile is used.

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
            percentile_val = 0.95_real64
        end if

        ! Always initialize permutation array
        call init_perm(tmp_perm)

        call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                          loess_x, loess_y, loess_n, ierr)
        if (is_err(ierr)) return
        call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right)
        call identify_outliers(n_genes, rdi, tmp_work_array, tmp_perm, is_outlier, threshold, quantile, percentile_val)
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
subroutine identify_outliers_c(n_genes, rdi, sorted_rdi, perm, is_outlier, threshold, quantile, percentile, ierr) &
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
        !! Percentile threshold as a fraction in [0,1] for outlier detection
    real(c_double), intent(out), target :: quantile(n_genes)
        !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
        !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
        !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
        !! upper-tail quantile is used.
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
    M_CHECK_NON_NULL(quantile)

    M_ALLOCATE(is_outlier_f(n_genes))

    call set_ok(ierr)

    call identify_outliers(n_genes, rdi, sorted_rdi, perm, is_outlier_f, threshold, quantile, percentile)

    call logical_as_c_int(is_outlier_f, is_outlier)
end subroutine identify_outliers_c

!> C wrapper for detect_outliers.
!| Calls detect_outliers with C-compatible types for external interface.
subroutine detect_outliers_c(n_genes, n_families, distances, gene_to_fam, &
                             tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right, &
                             is_outlier, loess_x, loess_y, loess_n, quantile, ierr, &
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
    real(c_double), intent(out), target :: quantile(n_genes)
        !! Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
        !! observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
        !! Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
        !! upper-tail quantile is used.
    integer(c_int), intent(out), target :: ierr
        !! Error code
    real(c_double), intent(in), target :: percentile
        !! Percentile threshold as a fraction in [0,1] for outlier detection

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
    M_CHECK_NON_NULL(quantile)
    M_CHECK_NON_NULL(percentile)

    M_ALLOCATE(is_outlier_f(n_genes))

    call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                         tmp_work_array, tmp_perm, tmp_stack_left, tmp_stack_right, &
                         is_outlier_f, loess_x, loess_y, loess_n, quantile, ierr, &
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
