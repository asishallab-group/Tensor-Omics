#include <src/macros.h>

!> Variance-matched noise model: the exact model with a variance-band neighbourhood.
!|
!| Implements implementation notes A ("vband") and B ("vshrunk") of the benchmark
!| repo. The plain models build each gene's null from the residuals of its nearest
!| neighbours in MEAN expression only; when genes at the same mean differ in noise
!| level, that pool is a scale mixture -- too wide for quiet genes (lost power), too
!| narrow for noisy ones (false positives). This module keeps the empirical null of
!| `noise_model_exact` UNCHANGED (pairwise differences of the neighbours' raw
!| residuals, per-side `1/sqrt(n)` scaling, exact tail count) and changes only WHICH
!| genes become neighbours: a candidate is accepted only if its matching variable
!| `u` lies within a two-sided band `|u_j - u_g| <= tau_var` around the target's.
!|
!| The two notes differ ONLY in how `u` and `tau_var` are built, so they share this
!| module and are selected by `variance_mode`:
!|
!|   - `VARIANCE_MODE_BAND` (1, note A): `u` = raw detrended log-variance,
!|     `tau_var = c_tau * sqrt(trigamma(d/2))`.
!|   - `VARIANCE_MODE_SHRUNK` (2, note B): `u` = log of the limma-style moderated
!|     variance over the trend prior, `tau_var = c_tau * w * sqrt(trigamma(d/2))`.
!|     Degrades to the plain fixed-k model when genes show no noise heterogeneity
!|     beyond sampling error (`d0 = Inf`, `u = 0`, every candidate eligible).
!|
!| The matching variable is computed HERE, from the same replicate matrices the
!| null is built from (`compute_variance_match`); the caller passes no `u`.
!|
!| **Equal-variance assumption in SELECTION only.** One `u` per gene (pooled
!| within-group variance, `d = n_case + n_control - 2` df) selects both the case
!| and the control pool. The null itself keeps separate per-side pools and per-side
!| scaling, so it makes no equal-variance assumption.
!|
!| **Reuse, not copy.** Sorted structures, neighbour ordering, the exact tail count,
!| the statistic and log2fc helpers come from `noise_model_exact`; the LOESS is
!| `tox_loess::loess_alloc`. Only the band gather, the matching variable and the
!| (inverse) polygamma helpers are new.
module noise_model_variance
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_is_finite
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_DIM_MISMATCH, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real
    use tox_loess, only: loess_alloc
    use noise_model_exact, only: sorted_data_t, prepare_sorted_data, &
                                 find_closest_helper, choose_index, add_residuals_to_pool_helper, &
                                 compute_pvalue, compute_gene_centres_helper, compute_log2fc_helper, &
                                 median_of_finite_helper, NOISE_LOG_OFFSET, &
                                 BETA_CENTRE_NONE, BETA_CENTRE_MEDIAN
    implicit none

    integer(int32), parameter :: VARIANCE_MODE_BAND = 1_int32
        !! `variance_mode`: note A -- band on the raw detrended log-variance.
    integer(int32), parameter :: VARIANCE_MODE_SHRUNK = 2_int32
        !! `variance_mode`: note B -- band on the shrunk (moderated) log-variance.

    integer(int32), parameter :: BAND_MAX_SCAN_FACTOR = 20_int32
        !! The band gather examines at most `BAND_MAX_SCAN_FACTOR * k` candidate genes.
        !! Bounds how far in mean the neighbourhood may drift when the band is strict;
        !! a pool that is still short then is KEPT short (tau is never relaxed).
    integer(int32), parameter :: BAND_MIN_GENES = 10_int32
        !! A pool with fewer accepted genes than this makes the gene not computable
        !! (p = -1). The residual gate (`< 10` residuals) applies on top.
    integer(int32), parameter :: MIN_POOL_RESIDUALS = 10_int32
        !! Residual gate, as in the scalar modules.
    real(real64), parameter :: VAR_FLOOR = 1.0e-8_real64
        !! Numerical guard in `l = log(s2 + VAR_FLOOR)` (log2-expression units). The
        !! number of genes it touches (`s2 < VAR_FLOOR`) is reported as `n_floored`.
    integer(int32), parameter :: LOESS_ROBUST_ITERS = 4_int32
        !! Robustness iterations of the bisquare-reweighted LOESS; `stats::loess`
        !! with `family = "symmetric"` uses 4.
    real(real64), parameter :: MAD_CONSTANT = 1.4826_real64
        !! Normal-consistency constant of the MAD, applied ONCE (R's `mad()`).
    real(real64), parameter :: TRIGAMMA_INV_LO = 1.0e-6_real64
        !! Lower end of the inverse-trigamma search interval (note B, step 4).
    real(real64), parameter :: TRIGAMMA_INV_HI = 1.0e6_real64
        !! Upper end of the inverse-trigamma search interval.

contains

    ! =========================================================================
    ! polygamma helpers
    ! =========================================================================

    !> Digamma `psi(x)` for `x > 0`: recurrence up to `x >= 10`, then the asymptotic
    !| series (first omitted term ~2e-14 there).
    pure function digamma_helper(x_in) result(psi)
        real(real64), intent(in) :: x_in
        !! Argument, `> 0`
        real(real64) :: psi
        !! `psi(x_in)`

        real(real64) :: x, f

        x = x_in
        psi = 0.0_real64
        do while (x < 10.0_real64)
            psi = psi - 1.0_real64 / x
            x = x + 1.0_real64
        end do
        f = 1.0_real64 / (x * x)
        psi = psi + log(x) - 0.5_real64 / x &
              - f * (1.0_real64 / 12.0_real64 - f * (1.0_real64 / 120.0_real64 - f * (1.0_real64 / 252.0_real64 &
              - f * (1.0_real64 / 240.0_real64 - f / 132.0_real64))))
    end function digamma_helper

    !> Trigamma `psi1(x)` for `x > 0`: recurrence up to `x >= 10`, then the asymptotic
    !| series (first omitted term `5/(66 x^11)` ~ 8e-13 there; at `x >= 6` it would be
    !| 2e-10). `psi1(d/2)` is the sampling variance of `log(s2)` at `d` df.
    pure function trigamma_helper(x_in) result(psi1)
        real(real64), intent(in) :: x_in
        !! Argument, `> 0`
        real(real64) :: psi1
        !! `psi1(x_in)`

        real(real64) :: x, f

        x = x_in
        psi1 = 0.0_real64
        do while (x < 10.0_real64)
            psi1 = psi1 + 1.0_real64 / (x * x)
            x = x + 1.0_real64
        end do
        f = 1.0_real64 / (x * x)
        psi1 = psi1 + 1.0_real64 / x + 0.5_real64 * f &
               + f / x * (1.0_real64 / 6.0_real64 - f * (1.0_real64 / 30.0_real64 &
               - f * (1.0_real64 / 42.0_real64 - f / 30.0_real64)))
    end function trigamma_helper

    !> Inverse trigamma: the `y` in `[TRIGAMMA_INV_LO, TRIGAMMA_INV_HI]` with
    !| `psi1(y) = target`, by bisection on `log(y)` (`psi1` is decreasing). Values
    !| outside the interval's image are clamped to its ends, as `uniroot` on that
    !| interval would refuse; 200 halvings resolve the 1e12 ratio far below rounding.
    pure function trigamma_inverse_helper(target) result(y)
        real(real64), intent(in) :: target
        !! Target value, `> 0`
        real(real64) :: y
        !! Solution of `psi1(y) = target`

        real(real64) :: lo, hi, mid
        integer(int32) :: it

        lo = log(TRIGAMMA_INV_LO)
        hi = log(TRIGAMMA_INV_HI)
        if (target >= trigamma_helper(TRIGAMMA_INV_LO)) then
            y = TRIGAMMA_INV_LO
            return
        end if
        if (target <= trigamma_helper(TRIGAMMA_INV_HI)) then
            y = TRIGAMMA_INV_HI
            return
        end if
        do it = 1, 200
            mid = 0.5_real64 * (lo + hi)
            if (trigamma_helper(exp(mid)) > target) then
                lo = mid
            else
                hi = mid
            end if
        end do
        y = exp(0.5_real64 * (lo + hi))
    end function trigamma_inverse_helper

    ! =========================================================================
    ! matching variable
    ! =========================================================================

    !> Per-gene pooled within-group log-variance and the trend axis (notes A/B, §2).
    !|
    !| On `L = log2(max(x,0) + 1)` of the SAME input the null is built from
    !| (independent of `norm_method`: the matching variable is always on the log2
    !| scale):
    !|
    !|   s2_g = [ sum_case (L - mean_case)^2 + sum_ctrl (L - mean_ctrl)^2 ] / d
    !|   l_g  = log(s2_g + VAR_FLOOR)
    !|   m_g  = log2( mean over ALL samples of the linear input + 1 )
    pure subroutine compute_logvar_helper(replicates_case, n_rep_case, replicates_control, n_rep_control, &
                                          n_genes, s2, l, m, n_floored)
        integer(int32), intent(in) :: n_rep_case
        !! Case replicates per gene
        integer(int32), intent(in) :: n_rep_control
        !! Control replicates per gene
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        real(real64), dimension(n_rep_case, n_genes), intent(in) :: replicates_case
        !! Case replicates, linear (pre-log)
        real(real64), dimension(n_rep_control, n_genes), intent(in) :: replicates_control
        !! Control replicates, linear (pre-log)
        real(real64), dimension(n_genes), intent(out) :: s2
        !! Output: pooled within-group variance of `L`
        real(real64), dimension(n_genes), intent(out) :: l
        !! Output: `log(s2 + VAR_FLOOR)`
        real(real64), dimension(n_genes), intent(out) :: m
        !! Output: trend axis
        integer(int32), intent(out) :: n_floored
        !! Output: genes with `s2 < VAR_FLOOR`

        integer(int32) :: g
        real(real64) :: log2_factor, lc(n_rep_case), lt(n_rep_control), d

        log2_factor = 1.0_real64 / log(2.0_real64)
        d = real(n_rep_case + n_rep_control - 2, real64)
        n_floored = 0
        do g = 1, n_genes
            lc = log(max(replicates_case(:, g), 0.0_real64) + NOISE_LOG_OFFSET) * log2_factor
            lt = log(max(replicates_control(:, g), 0.0_real64) + NOISE_LOG_OFFSET) * log2_factor
            s2(g) = (sum((lc - sum(lc) / real(n_rep_case, real64))**2) &
                     + sum((lt - sum(lt) / real(n_rep_control, real64))**2)) / d
            if (s2(g) < VAR_FLOOR) n_floored = n_floored + 1
            l(g) = log(s2(g) + VAR_FLOOR)
            m(g) = log((sum(replicates_case(:, g)) + sum(replicates_control(:, g))) &
                       / real(n_rep_case + n_rep_control, real64) + NOISE_LOG_OFFSET) * log2_factor
        end do
    end subroutine compute_logvar_helper

    !> Build the matching variable `u` and the band half-width `tau_var`.
    !|
    !| **Mode A** (`VARIANCE_MODE_BAND`): `S(m)` = LOESS of `l` on `m`;
    !| `u = l - S(m)`; `tau_var = c_tau * sqrt(psi1(d/2))`. LOESS DETRENDS -- a
    !| trend over mean expression carries no gene-specific information; the
    !| gene-specific signal is the deviation from it. `d0`, `w` are not used and
    !| returned as -1; `scatter` is the MAD^2 of `u` (diagnostic).
    !|
    !| **Mode B** (`VARIANCE_MODE_SHRUNK`), note B §2 step by step:
    !|   1. `e = l - psi(d/2) + log(d/2)` (bias-corrected log variance)
    !|   2. `T(m)` = LOESS of `e` on `m`; `r = e - T(m)`
    !|   3. `V = mad(r)^2` (robust, default; MAD constant applied once, as R's `mad()`
    !|      and `estimate_dispersion_trend.R` do) or the sample variance of `r`;
    !|      `ev = V - psi1(d/2)`
    !|   4. `d0 = Inf` if `ev <= 0`, else `d0 = 2 y` with `psi1(y) = ev`
    !|   5. `s0sq = exp(T + psi(d0/2) - log(d0/2))`, or `exp(T)` when `d0 = Inf`
    !|   6. `s2m = (d0 s0sq + d s2) / (d0 + d)` -- on the VARIANCE scale, limma-style
    !|   7. `u = log(s2m / s0sq)` (identically 0 when `d0 = Inf`)
    !|   8. `w = d / (d + d0)` (0 when `d0 = Inf`)
    !|   10. `tau_var = c_tau * w * sqrt(psi1(d/2))`; when `d0 = Inf`, `u = 0` and
    !|       `tau_var = 1` so the band code path still runs and accepts everything.
    !|
    !| LOESS: `tox_loess::loess_alloc` (netlib dloess). Measured against
    !| `stats::loess(l ~ m, span = 0.3, degree = 1, family = "symmetric")` at R's
    !| default settings, mode A's `u` agrees to 5e-15. The only parametric element
    !| anywhere is the `d0` estimate, and it sets only the SELECTION weight, never
    !| the null.
    subroutine compute_variance_match(replicates_case, n_rep_case, replicates_control, n_rep_control, &
                                      n_genes, variance_mode, c_tau, loess_span, loess_degree, &
                                      loess_robust, robust_var, &
                                      match_u, tau_var, d0, w, scatter, n_floored, ierr)
        integer(int32), intent(in) :: n_rep_case
        !! Case replicates per gene (>= 2)
        integer(int32), intent(in) :: n_rep_control
        !! Control replicates per gene (>= 2)
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        real(real64), dimension(n_rep_case, n_genes), intent(in) :: replicates_case
        !! Case replicates, linear (pre-log)
        real(real64), dimension(n_rep_control, n_genes), intent(in) :: replicates_control
        !! Control replicates, linear (pre-log)
        integer(int32), intent(in) :: variance_mode
        !! `VARIANCE_MODE_BAND` (1) or `VARIANCE_MODE_SHRUNK` (2)
        real(real64), intent(in) :: c_tau
        !! Band constant `c` (> 0)
        real(real64), intent(in) :: loess_span
        !! LOESS span, `(0, 1]`
        integer(int32), intent(in) :: loess_degree
        !! LOESS degree, 0..2
        integer(int32), intent(in) :: loess_robust
        !! 1 = bisquare-robust LOESS (`family = "symmetric"`), 0 = plain
        integer(int32), intent(in) :: robust_var
        !! Mode B only: 1 = MAD-based scatter (default), 0 = sample variance
        real(real64), dimension(n_genes), intent(out) :: match_u
        !! Output: matching variable `u`, original gene order
        real(real64), intent(out) :: tau_var
        !! Output: band half-width
        real(real64), intent(out) :: d0
        !! Output: prior df (mode B; +Inf when no excess scatter); -1 in mode A
        real(real64), intent(out) :: w
        !! Output: weight on the gene's own variance (mode B); -1 in mode A
        real(real64), intent(out) :: scatter
        !! Output: scatter of the detrended variable (`V` of note B step 3; mode A:
        !! the same robust estimate of `u`'s spread, for comparison with `psi1(d/2)`)
        integer(int32), intent(out) :: n_floored
        !! Output: genes whose variance hit `VAR_FLOOR`
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64), allocatable :: s2(:), l(:), m(:), trend(:), r(:)
        real(real64), allocatable :: med_work(:)
        integer(int32), allocatable :: med_perm(:), med_left(:), med_right(:)
        real(real64) :: d, psi1_d, med, mad_raw, ev, y, mean_r
        integer(int32) :: n_finite, g

        call set_ok(ierr)
        tau_var = 0.0_real64
        d0 = -1.0_real64
        w = -1.0_real64
        scatter = 0.0_real64
        n_floored = 0
        match_u = 0.0_real64

        M_ALLOCATE(s2(n_genes))
        M_ALLOCATE(l(n_genes))
        M_ALLOCATE(m(n_genes))
        M_ALLOCATE(trend(n_genes))
        M_ALLOCATE(r(n_genes))
        M_ALLOCATE(med_work(n_genes))
        M_ALLOCATE(med_perm(n_genes))
        M_ALLOCATE(med_left(n_genes))
        M_ALLOCATE(med_right(n_genes))

        d = real(n_rep_case + n_rep_control - 2, real64)
        psi1_d = trigamma_helper(0.5_real64 * d)

        call compute_logvar_helper(replicates_case, n_rep_case, replicates_control, n_rep_control, &
                                   n_genes, s2, l, m, n_floored)

        ! Mode A detrends l; mode B detrends the bias-corrected e (step 1).
        if (variance_mode == VARIANCE_MODE_SHRUNK) &
            l = l - digamma_helper(0.5_real64 * d) + log(0.5_real64 * d)

        call loess_alloc(m, l, loess_span, loess_degree, trend, loess_robust, LOESS_ROBUST_ITERS, ierr)
        if (is_err(ierr)) return
        r = l - trend

        ! Robust scatter of the detrended variable: (1.4826 * median|r - median r|)^2.
        call median_of_finite_helper(r, n_genes, med_work, med_perm, med_left, med_right, med, n_finite)
        do g = 1, n_genes
            med_work(g) = abs(r(g) - med)
        end do
        r = med_work
        call median_of_finite_helper(r, n_genes, med_work, med_perm, med_left, med_right, mad_raw, n_finite)
        r = l - trend
        scatter = (MAD_CONSTANT * mad_raw)**2

        if (variance_mode == VARIANCE_MODE_BAND) then
            match_u = r
            tau_var = c_tau * sqrt(psi1_d)
            return
        end if

        ! ---- mode B --------------------------------------------------------------
        if (robust_var == 0) then
            mean_r = sum(r) / real(n_genes, real64)
            scatter = sum((r - mean_r)**2) / real(max(n_genes - 1, 1), real64)
        end if
        ev = scatter - psi1_d

        if (ev <= 0.0_real64) then
            d0 = ieee_value(1.0_real64, ieee_positive_inf)
            w = 0.0_real64
            match_u = 0.0_real64
            tau_var = 1.0_real64
            return
        end if

        y = trigamma_inverse_helper(ev)
        d0 = 2.0_real64 * y
        w = d / (d + d0)
        ! trend now holds T(m); s0sq = exp(T + psi(d0/2) - log(d0/2)).
        do g = 1, n_genes
            med = exp(trend(g) + digamma_helper(y) - log(y))          ! s0sq at gene g
            match_u(g) = log(((d0 * med + d * s2(g)) / (d0 + d)) / med)
        end do
        tau_var = c_tau * w * sqrt(psi1_d)
    end subroutine compute_variance_match

    ! =========================================================================
    ! band gather
    ! =========================================================================

    !> Gather a residual pool from the `k` nearest-in-mean genes whose matching
    !| variable lies inside the band (note A §3.2).
    !|
    !| Same contract as `noise_model_exact::gather_residuals_helper` (pure, no
    !| allocation, writes into the caller's pool, hard cap `max_pool_size`).
    !| Candidates are visited in mean-closeness order -- the closest slot first,
    !| then `choose_index`'s left/right expansion -- and each visited slot counts
    !| towards `n_scanned`. A candidate `j` is accepted iff `u_j` is finite and
    !| `|u_j - target_u| <= tau_var`. The walk stops when `k` genes are accepted,
    !| `scan_max` slots were examined, both sides are exhausted, or the pool is full.
    !| A short pool is returned as is; the band is never relaxed.
    !|
    !| **Self-inclusion:** the target is (normally) the closest gene to its own mean
    !| and always passes the band (distance 0), so it is in its own pool, as in the
    !| plain models. The Python prototypes excluded it; the effect is about 1/k.
    !|
    !| With every candidate eligible (huge `tau_var`) this visits exactly the genes
    !| the plain gather visits with `k_start = k_step = k_max = k`, in the same order.
    pure subroutine gather_residuals_band_helper(target_mean, target_u, sorted_data, match_u, &
                                                 k, tau_var, scan_max, &
                                                 pooled_residuals, n_pooled, n_genes_pool, n_scanned, &
                                                 max_pool_size)
        real(real64), intent(in) :: target_mean
        !! Mean of the target gene on this side (the kNN coordinate)
        real(real64), intent(in) :: target_u
        !! Matching variable of the target gene
        type(sorted_data_t), intent(in) :: sorted_data
        !! Sorted structure of this side
        real(real64), dimension(:), intent(in) :: match_u
        !! Matching variable, indexed by ORIGINAL gene index
        integer(int32), intent(in) :: k
        !! Number of accepted genes wanted
        real(real64), intent(in) :: tau_var
        !! Band half-width
        integer(int32), intent(in) :: scan_max
        !! Maximum number of candidate slots examined
        real(real64), intent(out) :: pooled_residuals(:)
        !! Output residual pool (at least `max_pool_size` long)
        integer(int32), intent(out) :: n_pooled
        !! Residuals written
        integer(int32), intent(out) :: n_genes_pool
        !! Genes accepted
        integer(int32), intent(out) :: n_scanned
        !! Candidate slots examined
        integer(int32), intent(in) :: max_pool_size
        !! Hard cap on the residual count

        integer(int32) :: pos, idx, left_cand, right_cand, new_size
        real(real64) :: u_j

        n_pooled = 0
        n_genes_pool = 0
        n_scanned = 0

        pos = find_closest_helper(target_mean, sorted_data%means_sorted, sorted_data%n_genes)
        if (pos == 0) return
        left_cand = pos - 1
        right_cand = pos + 1
        idx = pos

        do
            n_scanned = n_scanned + 1
            u_j = match_u(sorted_data%original_indices(idx))
            if (ieee_is_finite(u_j)) then
                if (abs(u_j - target_u) <= tau_var) then
                    new_size = min(n_pooled + sorted_data%n_residuals(idx), max_pool_size)
                    call add_residuals_to_pool_helper(pooled_residuals, n_pooled, &
                                                      sorted_data%residuals_packed(:, idx), &
                                                      sorted_data%n_residuals(idx), new_size)
                    n_pooled = new_size
                    n_genes_pool = n_genes_pool + 1
                end if
            end if

            if (n_genes_pool >= k .or. n_scanned >= scan_max .or. n_pooled >= max_pool_size) exit
            if (left_cand < 1 .and. right_cand > sorted_data%n_genes) exit
            call choose_index(sorted_data%means_sorted, sorted_data%n_genes, target_mean, idx, &
                              left_cand, right_cand)
        end do
    end subroutine gather_residuals_band_helper

    ! =========================================================================
    ! pipeline
    ! =========================================================================

    !> Validate, build the matching variable and the sorted structures, and run the
    !| per-gene loop of the variance-matched exact model.
    !|
    !| Inputs and the statistic follow the scalar modules exactly: the model computes
    !| the means (kNN coordinate), the observed statistic (`observed_statistic_own`,
    !| centred on request via `beta_centre`), `log2fc` and `delta_hat` from the
    !| replicates. The null is `noise_model_exact`'s: per-side `1/sqrt(n)` scaled
    !| pools, exact pairwise tail count. Only the neighbourhood differs: a FIXED
    !| `k` genes accepted by the variance band (no adaptive growth, no `tau`).
    !|
    !| A gene is not computable (p = -1) when its `u` is non-finite, when either
    !| side accepts fewer than `BAND_MIN_GENES` genes, or when either pool holds
    !| fewer than `MIN_POOL_RESIDUALS` residuals.
    subroutine compute_noise_pvalue_pipeline_variance( &
        replicates_case, n_genes_case, n_replicates_case, &
        replicates_control, n_genes_control, n_replicates_control, &
        compute_pvalue_own, beta_centre, &
        n_genes, norm_method, variance_mode, k, c_tau, &
        loess_span, loess_degree, loess_robust, robust_var, max_pool_size, &
        pvalues_own, observed_statistic_own, log2fc, delta_hat, n_genes_with_pvalue, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        n_scanned_own_case, n_scanned_own_control, &
        match_u, tau_var, d0, w, scatter, n_floored, ierr)

        integer(int32), intent(in) :: n_genes_case
        !! Number of genes in the case group
        integer(int32), intent(in) :: n_replicates_case
        !! Number of case replicates (>= 2)
        integer(int32), intent(in) :: n_genes_control
        !! Number of genes in the control group
        integer(int32), intent(in) :: n_replicates_control
        !! Number of control replicates (>= 2)
        integer(int32), intent(in) :: n_genes
        !! Number of genes; must equal both group gene counts
        real(real64), dimension(n_replicates_case, n_genes_case), intent(in) :: replicates_case
        !! Case replicate matrix (samples x genes), pre-log
        real(real64), dimension(n_replicates_control, n_genes_control), intent(in) :: replicates_control
        !! Control replicate matrix (samples x genes), pre-log
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_own
        !! 1 if the gene should be tested, 0 otherwise
        integer(int32), intent(in) :: beta_centre
        !! `BETA_CENTRE_NONE` (0) or `BETA_CENTRE_MEDIAN` (1), as in the scalar modules
        integer(int32), intent(in) :: norm_method
        !! 0 = linear residuals/statistic; non-zero = log2. The matching variable is
        !! always on the log2 scale.
        integer(int32), intent(in) :: variance_mode
        !! `VARIANCE_MODE_BAND` (1, note A) or `VARIANCE_MODE_SHRUNK` (2, note B)
        integer(int32), intent(in) :: k
        !! Fixed number of accepted neighbour GENES per pool
        real(real64), intent(in) :: c_tau
        !! Band constant `c` (> 0); `tau_var` follows from it and the mode
        real(real64), intent(in) :: loess_span
        !! Trend LOESS span, `(0, 1]` (notes: 0.3)
        integer(int32), intent(in) :: loess_degree
        !! Trend LOESS degree, 0..2 (notes: 1)
        integer(int32), intent(in) :: loess_robust
        !! 1 = robust (`family = "symmetric"`, notes' default), 0 = plain
        integer(int32), intent(in) :: robust_var
        !! Mode B: 1 = MAD scatter (default), 0 = sample variance. Ignored in mode A.
        integer(int32), intent(in) :: max_pool_size
        !! Hard cap on residuals per pool
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: observed_statistic_own
        !! Output: observed statistic (centred when `beta_centre = 1`)
        real(real64), dimension(n_genes), intent(out) :: log2fc
        !! Output: effect size, log2 fold change of the linear means; never centred
        real(real64), intent(out) :: delta_hat
        !! Output: median observed statistic across genes
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Output: genes that received a p-value
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Output: case pool size in residuals (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Output: control pool size in residuals (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: n_scanned_own_case
        !! Output: case candidate slots examined (-1 if not gathered)
        integer(int32), dimension(n_genes), intent(out) :: n_scanned_own_control
        !! Output: control candidate slots examined (-1 if not gathered)
        real(real64), dimension(n_genes), intent(out) :: match_u
        !! Output: matching variable `u`
        real(real64), intent(out) :: tau_var
        !! Output: band half-width actually used
        real(real64), intent(out) :: d0
        !! Output: prior df (mode B; +Inf = no excess scatter); -1 in mode A
        real(real64), intent(out) :: w
        !! Output: shrinkage weight `d/(d+d0)` (mode B); -1 in mode A
        real(real64), intent(out) :: scatter
        !! Output: robust (or, mode B with `robust_var = 0`, classical) scatter of the
        !! detrended log-variance
        integer(int32), intent(out) :: n_floored
        !! Output: genes whose variance hit `VAR_FLOOR`
        integer(int32), intent(out) :: ierr
        !! Error code

        type(sorted_data_t) :: sorted_case, sorted_control
        real(real64), allocatable :: means_case(:), means_control(:), centres_control(:)
        real(real64), allocatable :: med_work(:), pool_case(:), pool_control(:)
        integer(int32), allocatable :: med_perm(:), med_left(:), med_right(:)
        integer(int32) :: i_gene, n_finite, sort_ierr, scan_max
        integer(int32) :: n_pool_case, n_pool_control, n_gp_case, n_gp_control
        real(real64) :: scale_case, scale_control

        call set_ok(ierr)
        pvalues_own = -1.0_real64
        observed_statistic_own = 0.0_real64
        log2fc = 0.0_real64
        delta_hat = 0.0_real64
        n_genes_with_pvalue = 0
        neighborhood_size_own_case = -1
        neighborhood_size_own_control = -1
        n_scanned_own_case = -1
        n_scanned_own_control = -1
        match_u = 0.0_real64
        tau_var = 0.0_real64
        d0 = -1.0_real64
        w = -1.0_real64
        scatter = 0.0_real64
        n_floored = 0

        call validate_dimension_size(n_genes_case, ierr)
        call validate_dimension_size(n_replicates_case, ierr)
        call validate_dimension_size(n_genes_control, ierr)
        call validate_dimension_size(n_replicates_control, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(k, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        if (variance_mode /= VARIANCE_MODE_BAND .and. variance_mode /= VARIANCE_MODE_SHRUNK) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (beta_centre /= BETA_CENTRE_NONE .and. beta_centre /= BETA_CENTRE_MEDIAN) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (.not. (c_tau > 0.0_real64)) call set_err(ierr, ERR_INVALID_INPUT)
        if (loess_robust /= 0 .and. loess_robust /= 1) call set_err(ierr, ERR_INVALID_INPUT)
        if (robust_var /= 0 .and. robust_var /= 1) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return
        ! One gene index for both matrices and every per-gene array (see eval doc §18).
        if (n_genes_case /= n_genes .or. n_genes_control /= n_genes) then
            call set_err(ierr, ERR_DIM_MISMATCH)
            return
        end if
        call validate_all_in_range_real(replicates_case, n_replicates_case * n_genes_case, ierr)
        call validate_all_in_range_real(replicates_control, n_replicates_control * n_genes_control, ierr)
        if (is_err(ierr)) return

        ! --- means, statistic, effect size (as in the scalar modules) -------------
        M_ALLOCATE(means_case(n_genes))
        M_ALLOCATE(means_control(n_genes))
        M_ALLOCATE(centres_control(n_genes))
        M_ALLOCATE(med_work(n_genes))
        M_ALLOCATE(med_perm(n_genes))
        M_ALLOCATE(med_left(n_genes))
        M_ALLOCATE(med_right(n_genes))

        call compute_gene_centres_helper(replicates_case, n_replicates_case, n_genes, norm_method, &
                                         means_case, observed_statistic_own)
        call compute_gene_centres_helper(replicates_control, n_replicates_control, n_genes, norm_method, &
                                         means_control, centres_control)
        observed_statistic_own = observed_statistic_own - centres_control
        call compute_log2fc_helper(means_case, means_control, n_genes, log2fc)
        call median_of_finite_helper(observed_statistic_own, n_genes, med_work, med_perm, &
                                     med_left, med_right, delta_hat, n_finite)
        if (beta_centre == BETA_CENTRE_MEDIAN) observed_statistic_own = observed_statistic_own - delta_hat

        ! --- matching variable ------------------------------------------------------
        ! prepare_sorted_data rejects n < 2 per side; check it here first so the
        ! variance step never divides by d = 0.
        if (n_replicates_case < 2 .or. n_replicates_control < 2) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        call compute_variance_match(replicates_case, n_replicates_case, replicates_control, n_replicates_control, &
                                    n_genes, variance_mode, c_tau, loess_span, loess_degree, &
                                    loess_robust, robust_var, &
                                    match_u, tau_var, d0, w, scatter, n_floored, ierr)
        if (is_err(ierr)) return

        ! --- sorted structures (the exact module's) -------------------------------
        call prepare_sorted_data(means_case, replicates_case, n_replicates_case, n_genes, norm_method, &
                                 sorted_case, sort_ierr)
        call set_err(ierr, sort_ierr)
        call prepare_sorted_data(means_control, replicates_control, n_replicates_control, n_genes, norm_method, &
                                 sorted_control, sort_ierr)
        call set_err(ierr, sort_ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(pool_case(max_pool_size * 2))
        M_ALLOCATE(pool_control(max_pool_size * 2))

        scan_max = int(min(int(BAND_MAX_SCAN_FACTOR, int64) * int(k, int64), int(huge(1_int32), int64)), int32)
        ! Per-side 1/sqrt(n): the noise_model_exact null, unchanged.
        scale_case = 1.0_real64 / sqrt(real(n_replicates_case, real64))
        scale_control = 1.0_real64 / sqrt(real(n_replicates_control, real64))

        do i_gene = 1, n_genes
            if (.not. ieee_is_finite(match_u(i_gene))) cycle

            call gather_residuals_band_helper(means_case(i_gene), match_u(i_gene), sorted_case, match_u, &
                                              k, tau_var, scan_max, &
                                              pool_case, n_pool_case, n_gp_case, n_scanned_own_case(i_gene), &
                                              max_pool_size)
            if (n_gp_case < BAND_MIN_GENES .or. n_pool_case < MIN_POOL_RESIDUALS) cycle

            call gather_residuals_band_helper(means_control(i_gene), match_u(i_gene), sorted_control, match_u, &
                                              k, tau_var, scan_max, &
                                              pool_control, n_pool_control, n_gp_control, &
                                              n_scanned_own_control(i_gene), max_pool_size)
            if (n_gp_control < BAND_MIN_GENES .or. n_pool_control < MIN_POOL_RESIDUALS) cycle

            if (compute_pvalue_own(i_gene) /= 1) cycle

            pool_case(1:n_pool_case) = pool_case(1:n_pool_case) * scale_case
            pool_control(1:n_pool_control) = pool_control(1:n_pool_control) * scale_control
            call compute_pvalue(pool_case(1:n_pool_case), n_pool_case, &
                                pool_control(1:n_pool_control), n_pool_control, &
                                observed_statistic_own(i_gene), pvalues_own(i_gene), ierr)
            if (is_err(ierr)) return

            neighborhood_size_own_case(i_gene) = n_pool_case
            neighborhood_size_own_control(i_gene) = n_pool_control
            n_genes_with_pvalue = n_genes_with_pvalue + 1
        end do
    end subroutine compute_noise_pvalue_pipeline_variance

end module noise_model_variance

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

!> C-interoperable wrapper for `compute_noise_pvalue_pipeline_variance`.
!|
!| A separate entry point (not ABI-compatible with the plain scalar models): it
!| carries the variance-matching controls and diagnostics. Null-pointer checks only;
!| everything else is validated in Fortran.
subroutine compute_noise_pvalues_pipeline_variance_c( &
    replicates_case, n_genes_case, n_replicates_case, &
    replicates_control, n_genes_control, n_replicates_control, &
    compute_pvalue_own, beta_centre, &
    n_genes, norm_method, variance_mode, k, c_tau, &
    loess_span, loess_degree, loess_robust, robust_var, max_pool_size, &
    pvalues_own, observed_statistic_own, log2fc, delta_hat, n_genes_with_pvalue, &
    neighborhood_size_own_case, neighborhood_size_own_control, &
    n_scanned_own_case, n_scanned_own_control, &
    match_u, tau_var, d0, w, scatter, n_floored, ierr) &
    bind(C, name="compute_noise_pvalues_pipeline_variance_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use noise_model_variance, only: compute_noise_pvalue_pipeline_variance
    use safeguard
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes_case
    !! Number of genes in the case group
    integer(c_int), intent(in), target :: n_replicates_case
    !! Number of case replicates
    integer(c_int), intent(in), target :: n_genes_control
    !! Number of genes in the control group
    integer(c_int), intent(in), target :: n_replicates_control
    !! Number of control replicates
    integer(c_int), intent(in), target :: n_genes
    !! Number of genes
    real(c_double), dimension(n_replicates_case, n_genes_case), intent(in), target :: replicates_case
    !! Case replicates (samples x genes)
    real(c_double), dimension(n_replicates_control, n_genes_control), intent(in), target :: replicates_control
    !! Control replicates (samples x genes)
    integer(c_int), dimension(n_genes), intent(in), target :: compute_pvalue_own
    !! 1 = test this gene
    integer(c_int), intent(in), target :: beta_centre
    !! 0 = no centring; 1 = median centring of the statistic
    integer(c_int), intent(in), target :: norm_method
    !! 0 = linear; non-zero = log2
    integer(c_int), intent(in), target :: variance_mode
    !! 1 = variance band (note A); 2 = shrunk variance band (note B)
    integer(c_int), intent(in), target :: k
    !! Fixed neighbour-gene count
    real(c_double), intent(in), target :: c_tau
    !! Band constant
    real(c_double), intent(in), target :: loess_span
    !! Trend LOESS span
    integer(c_int), intent(in), target :: loess_degree
    !! Trend LOESS degree
    integer(c_int), intent(in), target :: loess_robust
    !! 1 = robust trend LOESS
    integer(c_int), intent(in), target :: robust_var
    !! Mode 2: 1 = MAD scatter, 0 = sample variance
    integer(c_int), intent(in), target :: max_pool_size
    !! Hard cap on residuals per pool
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_own
    !! Output: p-values (-1 if not computed)
    real(c_double), dimension(n_genes), intent(out), target :: observed_statistic_own
    !! Output: observed statistic
    real(c_double), dimension(n_genes), intent(out), target :: log2fc
    !! Output: effect size
    real(c_double), intent(out), target :: delta_hat
    !! Output: median statistic
    integer(c_int), intent(out), target :: n_genes_with_pvalue
    !! Output: genes with a p-value
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_case
    !! Output: case pool residuals
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_control
    !! Output: control pool residuals
    integer(c_int), dimension(n_genes), intent(out), target :: n_scanned_own_case
    !! Output: case candidates examined
    integer(c_int), dimension(n_genes), intent(out), target :: n_scanned_own_control
    !! Output: control candidates examined
    real(c_double), dimension(n_genes), intent(out), target :: match_u
    !! Output: matching variable
    real(c_double), intent(out), target :: tau_var
    !! Output: band half-width
    real(c_double), intent(out), target :: d0
    !! Output: prior df (mode 2)
    real(c_double), intent(out), target :: w
    !! Output: shrinkage weight (mode 2)
    real(c_double), intent(out), target :: scatter
    !! Output: scatter of the detrended log-variance
    integer(c_int), intent(out), target :: n_floored
    !! Output: genes at the variance floor
    integer(c_int), intent(out), target :: ierr
    !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(replicates_case)
    M_CHECK_NON_NULL(n_genes_case)
    M_CHECK_NON_NULL(n_replicates_case)
    M_CHECK_NON_NULL(replicates_control)
    M_CHECK_NON_NULL(n_genes_control)
    M_CHECK_NON_NULL(n_replicates_control)
    M_CHECK_NON_NULL(compute_pvalue_own)
    M_CHECK_NON_NULL(beta_centre)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(norm_method)
    M_CHECK_NON_NULL(variance_mode)
    M_CHECK_NON_NULL(k)
    M_CHECK_NON_NULL(c_tau)
    M_CHECK_NON_NULL(loess_span)
    M_CHECK_NON_NULL(loess_degree)
    M_CHECK_NON_NULL(loess_robust)
    M_CHECK_NON_NULL(robust_var)
    M_CHECK_NON_NULL(max_pool_size)
    M_CHECK_NON_NULL(pvalues_own)
    M_CHECK_NON_NULL(observed_statistic_own)
    M_CHECK_NON_NULL(log2fc)
    M_CHECK_NON_NULL(delta_hat)
    M_CHECK_NON_NULL(n_genes_with_pvalue)
    M_CHECK_NON_NULL(neighborhood_size_own_case)
    M_CHECK_NON_NULL(neighborhood_size_own_control)
    M_CHECK_NON_NULL(n_scanned_own_case)
    M_CHECK_NON_NULL(n_scanned_own_control)
    M_CHECK_NON_NULL(match_u)
    M_CHECK_NON_NULL(tau_var)
    M_CHECK_NON_NULL(d0)
    M_CHECK_NON_NULL(w)
    M_CHECK_NON_NULL(scatter)
    M_CHECK_NON_NULL(n_floored)

    call compute_noise_pvalue_pipeline_variance( &
        replicates_case, n_genes_case, n_replicates_case, &
        replicates_control, n_genes_control, n_replicates_control, &
        compute_pvalue_own, beta_centre, &
        n_genes, norm_method, variance_mode, k, c_tau, &
        loess_span, loess_degree, loess_robust, robust_var, max_pool_size, &
        pvalues_own, observed_statistic_own, log2fc, delta_hat, n_genes_with_pvalue, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        n_scanned_own_case, n_scanned_own_control, &
        match_u, tau_var, d0, w, scatter, n_floored, ierr)
end subroutine compute_noise_pvalues_pipeline_variance_c
