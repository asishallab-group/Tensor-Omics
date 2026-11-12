#include <src/macros.h>

!> Variance-aware noise model: the exact model with variance-aware neighbourhoods
!| and an optional cross-fitted, James-Stein-shrunk scaling of the null.
!|
!| Implements implementation notes A (variance band) and F (v2: James-Stein
!| shrinkage, cross-fitting, adaptive growth) of the benchmark repo; note F
!| supersedes note B (the parametric `d0` shrinkage, removed) and the
!| implementation part of note E v2 (whose estimator it uses). The null keeps the
!| construction of `noise_model_exact`: pairwise differences of the neighbours'
!| Bessel-corrected residuals, per-side `1/sqrt(n)` scaling, exact tail count.
!| Three orthogonal switches change it:
!|
!|   `variance_mode` (which genes enter the pools)
!|     - `VARIANCE_MODE_MEAN` (0): the plain mean-matched kNN of the exact model
!|       (`k_start`, `k_step`, `k_max`, `tau`); bit-identical to
!|       `noise_model_exact` with `scaling = SCALING_NONE`.
!|     - `VARIANCE_MODE_BAND` (1, note A): the `k` nearest-in-mean genes whose
!|       matching variable lies in `|u_j - u_g| <= tau_var`.
!|     - `VARIANCE_MODE_ADAPTIVE` (3, note F §5): block-wise growth inside a
!|       growth radius, with a variance band, a two-sided stability stop and an
!|       acceptance-rate stop.
!|     - `VARIANCE_MODE_ELLIPSE` (4, note F §5 "open alternative"): the `k` genes
!|       with the smallest joint distance `D = sqrt((dm / r_max)^2 + (du / tau_var)^2)`,
!|       `dm` the log2 mean distance and `du` the matching-variable distance, among
!|       those inside the ellipse `D <= 1`.
!|   `match_variable` (the `u` of band / adaptive)
!|     - `MATCH_RAW` (0, note A): raw detrended log2-variance,
!|       `tau_var = c_tau * sqrt(psi1(d/2))`.
!|     - `MATCH_JS` (1, note F §4): `u = log(s2s_pool / T)` of the James-Stein-shrunk
!|       variance, `tau_var = c_js * 1.4826 * MAD(u)` (no trigamma).
!|   `scaling` (null construction, note F §3)
!|     - `SCALING_NONE` (0): as the exact model.
!|     - `SCALING_XFIT_RAW` (1): cross-fitted with the raw group variances.
!|     - `SCALING_XFIT_JS` (2): cross-fitted with the James-Stein-shrunk ones.
!|     - `SCALING_SAME_GROUP` (3): benchmark-only negative control (note E v1):
!|       both pools scaled by the neighbour's own pooled shrunk variance.
!|
!| Cross-fitting divides each neighbour's CASE residuals by the sqrt of its CONTROL
!| variance and its CONTROL residuals by its CASE variance, and the tested
!| statistic by the gene's pooled variance. The scale applied to a residual is then
!| independent of that residual under H0 (equal group variances), so the null is
!| not internally studentised. It is applied once, to the packed residuals of the
!| sorted structures, so every gather and `compute_pvalue` are reused unchanged.
!|
!| **Shrinkage** (note F §2, note E v2 §2): moment-based, no distribution. On the
!| scale the model tests on (`y = log2(x + 1)` or `y = x`): unbiased group
!| variances `s2_c`, `s2_t`, pooled `s2`; the cross-group contrast
!| `delta = s2_c - s2_t` gives unbiased estimates of the sampling variances
!| (`Vhat`); the intensity `lambda = clip(sum w Vhat / sum w (s2 - T)^2, 0, 1)`
!| with `w = 1/T^2` (relative) or 1; `s2s = lambda T + (1 - lambda) s2`. The SAME
!| target `T` (global median, or LOESS of `s2` on `m`) serves the pooled and both
!| group-wise shrinkages (owner decision 28/10/2026).
!|
!| **Reuse, not copy.** Sorted structures, neighbour ordering, the plain gather,
!| the exact tail count, the statistic, log2fc and median helpers come from
!| `noise_model_exact`; the LOESS is `tox_loess::loess_alloc`.
module noise_model_variance
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_DIM_MISMATCH, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real
    use tox_loess, only: loess_alloc
    use noise_model_exact, only: sorted_data_t, prepare_sorted_data, gather_residuals_helper, &
                                 find_closest_helper, choose_index, add_residuals_to_pool_helper, &
                                 compute_pvalue, compute_gene_centres_helper, compute_log2fc_helper, &
                                 median_of_finite_helper, NOISE_LOG_OFFSET, &
                                 BETA_CENTRE_NONE, BETA_CENTRE_MEDIAN
    ! The Monte Carlo difference-of-means null of the bootstrap model, reused as is.
    use noise_model, only: compute_pvalue_bootstrap_mean_helper, RNG_BUFFER_MAX, N_BOOTSTRAP_DRAWS
    use f42_utils, only: init_random, sort_real_heapsort_expl_size
    implicit none

    integer(int32), parameter :: NULL_EXACT = 0_int32
        !! `null_construction`: every case x control residual pair, scaled 1/sqrt(n) per
        !! side, exact tail count (the exact model's null; default).
    integer(int32), parameter :: NULL_BOOTSTRAP = 1_int32
        !! `null_construction`: Monte Carlo difference of means -- `n_boot` draws of
        !! n residuals per side from the (scaled) pools, averaged (the bootstrap
        !! model's pooled null); p floor 1/(n_boot + 1).
    integer(int32), parameter :: VARIANCE_MODE_MEAN = 0_int32
        !! `variance_mode`: plain mean-matched kNN of the exact model.
    integer(int32), parameter :: VARIANCE_MODE_BAND = 1_int32
        !! `variance_mode`: note A -- fixed-k variance band.
    integer(int32), parameter :: VARIANCE_MODE_ADAPTIVE = 3_int32
        !! `variance_mode`: note F §5 -- adaptive variance-aware growth. (2 was note
        !! B's shrunk band, removed; it is rejected.)
    integer(int32), parameter :: VARIANCE_MODE_ELLIPSE = 4_int32
        !! `variance_mode`: note F §5 open alternative -- the `k` nearest genes in the
        !! joint (log2 mean, matching variable) distance, inside an ellipse.

    integer(int32), parameter :: MATCH_RAW = 0_int32
        !! `match_variable`: raw detrended log2-variance (note A).
    integer(int32), parameter :: MATCH_JS = 1_int32
        !! `match_variable`: James-Stein-shrunk variance (note F §4).

    integer(int32), parameter :: SCALING_NONE = 0_int32
        !! `scaling`: no scaling (the exact model's null).
    integer(int32), parameter :: SCALING_XFIT_RAW = 1_int32
        !! `scaling`: cross-fitted, raw group variances (lambda = 0).
    integer(int32), parameter :: SCALING_XFIT_JS = 2_int32
        !! `scaling`: cross-fitted, James-Stein-shrunk group variances.
    integer(int32), parameter :: SCALING_SAME_GROUP = 3_int32
        !! `scaling`: benchmark-only negative control, both pools scaled by the
        !! neighbour's own pooled shrunk variance (internally studentised).

    integer(int32), parameter :: SHRINK_TARGET_MEDIAN = 0_int32
        !! Shrinkage target: global median of the pooled variances.
    integer(int32), parameter :: SHRINK_TARGET_LOESS = 1_int32
        !! Shrinkage target: LOESS of the pooled variance on `m` (variance scale).

    integer(int32), parameter :: LAMBDA_RELATIVE = 0_int32
        !! Lambda weights `w_g = 1 / T_g^2` (relative squared error).
    integer(int32), parameter :: LAMBDA_ABSOLUTE = 1_int32
        !! Lambda weights `w_g = 1`.

    integer(int32), parameter :: STOP_K_MAX = 1_int32
        !! Adaptive stop reason: `k` genes accepted (or the pool is full).
    integer(int32), parameter :: STOP_RADIUS = 2_int32
        !! Adaptive stop reason: next candidate outside the growth radius.
    integer(int32), parameter :: STOP_STABILITY = 3_int32
        !! Adaptive stop reason: the last block changed the pool scale too much.
    integer(int32), parameter :: STOP_ACCEPTANCE = 4_int32
        !! Adaptive stop reason: acceptance rate of the last blocks below `a_min`.
    integer(int32), parameter :: STOP_SCAN_CAP = 5_int32
        !! Adaptive stop reason: scan cap reached.
    integer(int32), parameter :: STOP_EXHAUSTED = 6_int32
        !! Adaptive stop reason: candidates exhausted.

    integer(int32), parameter :: N_SHRINK_DIAG = 8_int32
        !! Length of `shrink_diag`: lambda_pool, lambda_c, lambda_t (used), the three
        !! unclipped ratios, effective dimension `p_hat`, participation ratio.

    integer(int32), parameter :: BAND_MAX_SCAN_FACTOR = 20_int32
        !! Band / adaptive gathers examine at most `BAND_MAX_SCAN_FACTOR * k` slots.
        !! A pool that is still short then is KEPT short (the band is never relaxed).
    integer(int32), parameter :: BAND_MIN_GENES = 10_int32
        !! A band / adaptive pool with fewer accepted genes makes the gene not
        !! computable (p = -1). The residual gate (`< 10` residuals) applies on top.
    integer(int32), parameter :: MIN_POOL_RESIDUALS = 10_int32
        !! Residual gate, as in the scalar modules.
    real(real64), parameter :: VAR_FLOOR = 1.0e-8_real64
        !! Floor for variances before a log or a division; genes touched are counted
        !! in `n_floored`.
    integer(int32), parameter :: LOESS_ROBUST_ITERS = 4_int32
        !! Robustness iterations of the bisquare-reweighted LOESS; `stats::loess`
        !! with `family = "symmetric"` uses 4.
    real(real64), parameter :: MAD_CONSTANT = 1.4826_real64
        !! Normal-consistency constant of the MAD, applied ONCE (R's `mad()`).

contains

    ! =========================================================================
    ! trigamma (note A's band width)
    ! =========================================================================

    !> Trigamma `psi1(x)` for `x > 0`: recurrence up to `x >= 10`, then the asymptotic
    !| series (first omitted term `5/(66 x^11)` ~ 8e-13 there). `psi1(d/2)` is the
    !| sampling variance of `log(s2)` at `d` df (note A's band width).
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

    ! =========================================================================
    ! note A: raw matching variable
    ! =========================================================================

    !> Per-gene pooled within-group log-variance and the trend axis (note A §2).
    !|
    !| On `L = log2(max(x,0) + 1)` of the SAME input the null is built from
    !| (independent of `norm_method`: note A's matching variable is always on the
    !| log2 scale):
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

    !> MAD (with `MAD_CONSTANT`) of the finite entries of `x(1:n)`.
    subroutine mad_helper(x, n, mad, work, perm, left, right)
        integer(int32), intent(in) :: n
        !! Number of values
        real(real64), dimension(n), intent(in) :: x
        !! Values
        real(real64), intent(out) :: mad
        !! Output: `MAD_CONSTANT * median|x - median(x)|`
        real(real64), dimension(n), intent(inout) :: work
        !! Work array
        integer(int32), dimension(n), intent(inout) :: perm, left, right
        !! Sort work arrays

        real(real64) :: med, dev(n)
        integer(int32) :: n_finite

        call median_of_finite_helper(x, n, work, perm, left, right, med, n_finite)
        dev = abs(x - med)
        call median_of_finite_helper(dev, n, work, perm, left, right, mad, n_finite)
        mad = MAD_CONSTANT * mad
    end subroutine mad_helper

    !> Note A's matching variable `u` and band half-width `tau_var`.
    !|
    !| `S(m)` = LOESS of `l` on `m`; `u = l - S(m)`; `tau_var = c_tau * sqrt(psi1(d/2))`.
    !| LOESS DETRENDS -- a trend over mean expression carries no gene-specific
    !| information; the gene-specific signal is the deviation from it. `scatter` is
    !| the MAD^2 of `u` (diagnostic, to compare with `psi1(d/2)`). Measured against
    !| `stats::loess(l ~ m, span = 0.3, degree = 1, family = "symmetric")` at R's
    !| default settings, `u` agrees to 5e-15.
    subroutine compute_variance_match(replicates_case, n_rep_case, replicates_control, n_rep_control, &
                                      n_genes, c_tau, loess_span, loess_degree, loess_robust, &
                                      match_u, tau_var, scatter, n_floored, ierr)
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
        real(real64), intent(in) :: c_tau
        !! Band constant `c` (> 0)
        real(real64), intent(in) :: loess_span
        !! LOESS span, `(0, 1]`
        integer(int32), intent(in) :: loess_degree
        !! LOESS degree, 0..2
        integer(int32), intent(in) :: loess_robust
        !! 1 = bisquare-robust LOESS (`family = "symmetric"`), 0 = plain
        real(real64), dimension(n_genes), intent(out) :: match_u
        !! Output: matching variable `u`, original gene order
        real(real64), intent(out) :: tau_var
        !! Output: band half-width
        real(real64), intent(out) :: scatter
        !! Output: MAD^2 of `u`
        integer(int32), intent(out) :: n_floored
        !! Output: genes whose variance hit `VAR_FLOOR`
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64), allocatable :: s2(:), l(:), m(:), trend(:), work(:)
        integer(int32), allocatable :: perm(:), left(:), right(:)
        real(real64) :: d, mad

        call set_ok(ierr)
        tau_var = 0.0_real64
        scatter = 0.0_real64
        n_floored = 0
        match_u = 0.0_real64

        M_ALLOCATE(s2(n_genes))
        M_ALLOCATE(l(n_genes))
        M_ALLOCATE(m(n_genes))
        M_ALLOCATE(trend(n_genes))
        M_ALLOCATE(work(n_genes))
        M_ALLOCATE(perm(n_genes))
        M_ALLOCATE(left(n_genes))
        M_ALLOCATE(right(n_genes))

        d = real(n_rep_case + n_rep_control - 2, real64)
        call compute_logvar_helper(replicates_case, n_rep_case, replicates_control, n_rep_control, &
                                   n_genes, s2, l, m, n_floored)
        call loess_alloc(m, l, loess_span, loess_degree, trend, loess_robust, LOESS_ROBUST_ITERS, ierr)
        if (is_err(ierr)) return
        match_u = l - trend
        call mad_helper(match_u, n_genes, mad, work, perm, left, right)
        scatter = mad**2
        tau_var = c_tau * sqrt(trigamma_helper(0.5_real64 * d))
    end subroutine compute_variance_match

    ! =========================================================================
    ! note F §2: James-Stein (moment) shrinkage
    ! =========================================================================

    !> Unbiased group variances, the pooled variance and the trend axis, on the
    !| scale the model tests on (`norm_method /= 0`: `y = log2(x + 1)`, else `y = x`).
    !| `m = log2(mean over all samples of the linear input + 1)` (note A's axis).
    pure subroutine compute_group_variances_helper(replicates_case, n_rep_case, replicates_control, &
                                                   n_rep_control, n_genes, norm_method, s2_c, s2_t, s2, m)
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
        integer(int32), intent(in) :: norm_method
        !! 0 = linear, non-zero = log2
        real(real64), dimension(n_genes), intent(out) :: s2_c
        !! Output: unbiased case variance
        real(real64), dimension(n_genes), intent(out) :: s2_t
        !! Output: unbiased control variance
        real(real64), dimension(n_genes), intent(out) :: s2
        !! Output: pooled variance, `d = n_c + n_t - 2` df
        real(real64), dimension(n_genes), intent(out) :: m
        !! Output: trend axis

        integer(int32) :: g
        real(real64) :: log2_factor, yc(n_rep_case), yt(n_rep_control), nc, nt

        log2_factor = 1.0_real64 / log(2.0_real64)
        nc = real(n_rep_case, real64)
        nt = real(n_rep_control, real64)
        do g = 1, n_genes
            if (norm_method /= 0) then
                yc = log(max(replicates_case(:, g), 0.0_real64) + NOISE_LOG_OFFSET) * log2_factor
                yt = log(max(replicates_control(:, g), 0.0_real64) + NOISE_LOG_OFFSET) * log2_factor
            else
                yc = replicates_case(:, g)
                yt = replicates_control(:, g)
            end if
            s2_c(g) = sum((yc - sum(yc) / nc)**2) / (nc - 1.0_real64)
            s2_t(g) = sum((yt - sum(yt) / nt)**2) / (nt - 1.0_real64)
            s2(g) = ((nc - 1.0_real64) * s2_c(g) + (nt - 1.0_real64) * s2_t(g)) / (nc + nt - 2.0_real64)
            m(g) = log((sum(replicates_case(:, g)) + sum(replicates_control(:, g))) / (nc + nt) &
                       + NOISE_LOG_OFFSET) * log2_factor
        end do
    end subroutine compute_group_variances_helper

    !> Moment-based (James-Stein type) shrinkage of the pooled and group variances.
    !|
    !|   delta       = s2_c - s2_t
    !|   Vhat_c      = ((n_t-1)/d) delta^2,  Vhat_t = ((n_c-1)/d) delta^2,
    !|   Vhat_pool   = ((n_c-1)(n_t-1)/d^2) delta^2          (balanced: exact; else approximate)
    !|   T           = median(s2) | LOESS of s2 on m;  T >= VAR_FLOOR
    !|   lambda_x    = clip( sum w Vhat_x / sum w (s2_x - T)^2, 0, 1 ),  w = 1/T^2 | 1
    !|   s2s_x       = lambda_x T + (1 - lambda_x) s2_x,     floored at VAR_FLOOR
    !|
    !| The same `T` is the target of all three. A zero denominator (every variance
    !| exactly on the target) gives lambda = 1. `shrink_diag` = [lambda_pool,
    !| lambda_c, lambda_t (as used), the three unclipped ratios, p_hat, participation
    !| ratio]; `p_hat = sum w Vhat_pool / max w Vhat_pool` (Bock's effective
    !| dimension for a diagonal Q), participation = (sum wV)^2 / sum (wV)^2.
    !|
    !| `lambda_fixed` / `target_fixed` (optional, tests only) replace the estimated
    !| lambdas / the target.
    subroutine compute_moment_shrinkage(s2_c, s2_t, s2, m, n_genes, n_rep_case, n_rep_control, &
                                        shrink_target, lambda_scale, loess_span, loess_degree, loess_robust, &
                                        target, s2s_c, s2s_t, s2s_pool, shrink_diag, n_floored, ierr, &
                                        lambda_fixed, target_fixed)
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        real(real64), dimension(n_genes), intent(in) :: s2_c
        !! Unbiased case variances
        real(real64), dimension(n_genes), intent(in) :: s2_t
        !! Unbiased control variances
        real(real64), dimension(n_genes), intent(in) :: s2
        !! Pooled variances
        real(real64), dimension(n_genes), intent(in) :: m
        !! Trend axis
        integer(int32), intent(in) :: n_rep_case
        !! Case replicates (>= 2)
        integer(int32), intent(in) :: n_rep_control
        !! Control replicates (>= 2)
        integer(int32), intent(in) :: shrink_target
        !! `SHRINK_TARGET_MEDIAN` or `SHRINK_TARGET_LOESS`
        integer(int32), intent(in) :: lambda_scale
        !! `LAMBDA_RELATIVE` or `LAMBDA_ABSOLUTE`
        real(real64), intent(in) :: loess_span
        !! LOESS span (LOESS target)
        integer(int32), intent(in) :: loess_degree
        !! LOESS degree (LOESS target)
        integer(int32), intent(in) :: loess_robust
        !! 1 = robust LOESS (LOESS target)
        real(real64), dimension(n_genes), intent(out) :: target
        !! Output: target `T_g`
        real(real64), dimension(n_genes), intent(out) :: s2s_c
        !! Output: shrunk case variances
        real(real64), dimension(n_genes), intent(out) :: s2s_t
        !! Output: shrunk control variances
        real(real64), dimension(n_genes), intent(out) :: s2s_pool
        !! Output: shrunk pooled variances
        real(real64), dimension(N_SHRINK_DIAG), intent(out) :: shrink_diag
        !! Output: lambdas, unclipped ratios, p_hat, participation ratio
        integer(int32), intent(out) :: n_floored
        !! Output: genes with any shrunk variance at `VAR_FLOOR`
        integer(int32), intent(out) :: ierr
        !! Error code
        real(real64), intent(in), optional :: lambda_fixed
        !! Tests only: use this lambda for all three shrinkages
        real(real64), intent(in), optional :: target_fixed
        !! Tests only: use this constant target

        real(real64), allocatable :: wv(:), work(:)
        integer(int32), allocatable :: perm(:), left(:), right(:)
        real(real64) :: d, nc1, nt1, med, raw(3), lam(3), wsum, wmax
        integer(int32) :: n_finite, g

        call set_ok(ierr)
        shrink_diag = 0.0_real64
        n_floored = 0

        M_ALLOCATE(wv(n_genes))
        M_ALLOCATE(work(n_genes))
        M_ALLOCATE(perm(n_genes))
        M_ALLOCATE(left(n_genes))
        M_ALLOCATE(right(n_genes))

        if (present(target_fixed)) then
            target = target_fixed
        else if (shrink_target == SHRINK_TARGET_LOESS) then
            call loess_alloc(m, s2, loess_span, loess_degree, target, loess_robust, LOESS_ROBUST_ITERS, ierr)
            if (is_err(ierr)) return
        else
            call median_of_finite_helper(s2, n_genes, work, perm, left, right, med, n_finite)
            target = med
        end if
        target = max(target, VAR_FLOOR)

        nc1 = real(n_rep_case - 1, real64)
        nt1 = real(n_rep_control - 1, real64)
        d = nc1 + nt1
        ! weights (relative: 1/T^2), held in wv for the three ratios
        if (lambda_scale == LAMBDA_RELATIVE) then
            wv = 1.0_real64 / target**2
        else
            wv = 1.0_real64
        end if
        ! work = delta^2
        work = (s2_c - s2_t)**2
        raw(1) = ratio_helper(sum(wv * (nc1 * nt1 / d**2) * work), sum(wv * (s2 - target)**2))
        raw(2) = ratio_helper(sum(wv * (nt1 / d) * work), sum(wv * (s2_c - target)**2))
        raw(3) = ratio_helper(sum(wv * (nc1 / d) * work), sum(wv * (s2_t - target)**2))
        lam = min(max(raw, 0.0_real64), 1.0_real64)
        if (present(lambda_fixed)) lam = lambda_fixed

        ! effective dimension and participation ratio of the pooled numerator
        wv = wv * (nc1 * nt1 / d**2) * work
        wsum = sum(wv)
        wmax = maxval(wv)
        if (wmax > 0.0_real64) then
            shrink_diag(7) = wsum / wmax
            shrink_diag(8) = wsum**2 / sum(wv**2)
        end if
        shrink_diag(1:3) = lam
        shrink_diag(4:6) = raw

        s2s_pool = lam(1) * target + (1.0_real64 - lam(1)) * s2
        s2s_c = lam(2) * target + (1.0_real64 - lam(2)) * s2_c
        s2s_t = lam(3) * target + (1.0_real64 - lam(3)) * s2_t
        do g = 1, n_genes
            if (s2s_pool(g) < VAR_FLOOR .or. s2s_c(g) < VAR_FLOOR .or. s2s_t(g) < VAR_FLOOR) &
                n_floored = n_floored + 1
        end do
        s2s_pool = max(s2s_pool, VAR_FLOOR)
        s2s_c = max(s2s_c, VAR_FLOOR)
        s2s_t = max(s2s_t, VAR_FLOOR)
    contains
        pure function ratio_helper(num, den) result(r)
            real(real64), intent(in) :: num, den
            real(real64) :: r
            if (den > 0.0_real64) then
                r = num / den
            else
                r = 1.0_real64
            end if
        end function ratio_helper
    end subroutine compute_moment_shrinkage

    ! =========================================================================
    ! band gather (note A)
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
    !| plain models.
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
    ! adaptive gather (note F §5)
    ! =========================================================================

    !> Adaptive, variance-aware neighbourhood growth (note F §5).
    !|
    !| Candidates are visited in mean-closeness order as in the band gather, in
    !| blocks of `block_size` slots. A candidate whose log2 mean distance
    !| `|log2(m_j + 1) - log2(m_g + 1)|` exceeds `r_max` is not visited (radius
    !| stop). Within a block, candidates passing the band (`|u_j - u_g| <= tau_var`)
    !| are accepted. After each block:
    !|   1. Stability (once at least `k_min` genes were accepted BEFORE the block, and
    !|      the block accepted any): if the pool scale `S` = mean |residual| changes
    !|      by more than `tau_stab` in either direction (`|S_new/S_old - 1|`), the
    !|      block is discarded and the walk stops.
    !|   2. A stop already met inside the block: `k` genes accepted or the pool full,
    !|      radius, scan cap (`scan_max` slots), candidates exhausted.
    !|   3. Acceptance rate: once `w_blocks` blocks were completed, if the accepted /
    !|      visited share over the last `w_blocks` blocks is below `a_min`, stop.
    !|
    !| Outputs as the band gather plus the final radius (largest log2 mean distance
    !| among accepted genes) and the stop reason (`STOP_*`). With `block_size = 1`,
    !| `r_max = tau_stab = Inf`, `a_min = 0` and `k_min = k` this visits and accepts
    !| exactly what the band gather does.
    pure subroutine gather_residuals_adaptive_helper(target_mean, target_u, sorted_data, match_u, &
                                                     tau_var, block_size, k_min, k, r_max, tau_stab, &
                                                     w_blocks, a_min, scan_max, &
                                                     pooled_residuals, n_pooled, n_genes_pool, n_scanned, &
                                                     radius, stop_reason, max_pool_size)
        real(real64), intent(in) :: target_mean
        !! Mean of the target gene on this side (the kNN coordinate)
        real(real64), intent(in) :: target_u
        !! Matching variable of the target gene
        type(sorted_data_t), intent(in) :: sorted_data
        !! Sorted structure of this side
        real(real64), dimension(:), intent(in) :: match_u
        !! Matching variable, indexed by ORIGINAL gene index
        real(real64), intent(in) :: tau_var
        !! Band half-width
        integer(int32), intent(in) :: block_size
        !! Candidate slots per block (>= 1)
        integer(int32), intent(in) :: k_min
        !! Accepted genes before the stability stop applies
        integer(int32), intent(in) :: k
        !! Maximum accepted genes (`k_max`)
        real(real64), intent(in) :: r_max
        !! Growth radius in log2 mean units (may be +Inf)
        real(real64), intent(in) :: tau_stab
        !! Stability tolerance on `|S_new/S_old - 1|` (may be +Inf)
        integer(int32), intent(in) :: w_blocks
        !! Blocks in the acceptance-rate window (>= 1)
        real(real64), intent(in) :: a_min
        !! Minimum acceptance rate over the window
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
        real(real64), intent(out) :: radius
        !! Largest log2 mean distance among accepted genes
        integer(int32), intent(out) :: stop_reason
        !! `STOP_*` code
        integer(int32), intent(in) :: max_pool_size
        !! Hard cap on the residual count

        integer(int32) :: pos, idx, left_cand, right_cand, new_size, j
        integer(int32) :: pool0, ngp0, block_acc, block_vis, pending, n_blocks, slot
        integer(int32) :: hist_acc(w_blocks), hist_vis(w_blocks)
        real(real64) :: u_j, lt, dist, abs_sum, abs0, rad0, s_old, s_new, rel
        logical :: have_cand

        n_pooled = 0
        n_genes_pool = 0
        n_scanned = 0
        radius = 0.0_real64
        stop_reason = STOP_EXHAUSTED
        abs_sum = 0.0_real64
        hist_acc = 0
        hist_vis = 0
        n_blocks = 0

        pos = find_closest_helper(target_mean, sorted_data%means_sorted, sorted_data%n_genes)
        if (pos == 0) return
        left_cand = pos - 1
        right_cand = pos + 1
        idx = pos
        have_cand = .true.
        lt = log(max(target_mean, 0.0_real64) + 1.0_real64) / log(2.0_real64)

        do
            pool0 = n_pooled
            ngp0 = n_genes_pool
            abs0 = abs_sum
            rad0 = radius
            block_acc = 0
            block_vis = 0
            pending = 0

            do while (block_vis < block_size)
                if (.not. have_cand) then
                    pending = STOP_EXHAUSTED
                    exit
                end if
                dist = abs(log(max(sorted_data%means_sorted(idx), 0.0_real64) + 1.0_real64) / log(2.0_real64) - lt)
                if (dist > r_max) then
                    pending = STOP_RADIUS
                    exit
                end if
                if (n_scanned >= scan_max) then
                    pending = STOP_SCAN_CAP
                    exit
                end if
                n_scanned = n_scanned + 1
                block_vis = block_vis + 1
                u_j = match_u(sorted_data%original_indices(idx))
                if (ieee_is_finite(u_j)) then
                    if (abs(u_j - target_u) <= tau_var) then
                        new_size = min(n_pooled + sorted_data%n_residuals(idx), max_pool_size)
                        do j = n_pooled + 1, new_size
                            abs_sum = abs_sum + abs(sorted_data%residuals_packed(j - n_pooled, idx))
                        end do
                        call add_residuals_to_pool_helper(pooled_residuals, n_pooled, &
                                                          sorted_data%residuals_packed(:, idx), &
                                                          sorted_data%n_residuals(idx), new_size)
                        n_pooled = new_size
                        n_genes_pool = n_genes_pool + 1
                        block_acc = block_acc + 1
                        radius = max(radius, dist)
                    end if
                end if
                if (left_cand < 1 .and. right_cand > sorted_data%n_genes) then
                    have_cand = .false.
                else
                    call choose_index(sorted_data%means_sorted, sorted_data%n_genes, target_mean, idx, &
                                      left_cand, right_cand)
                end if
                if (n_genes_pool >= k .or. n_pooled >= max_pool_size) then
                    pending = STOP_K_MAX
                    exit
                end if
            end do

            ! 1. two-sided stability of the pool scale
            if (ngp0 >= k_min .and. block_acc > 0 .and. pool0 > 0) then
                s_old = abs0 / real(pool0, real64)
                s_new = abs_sum / real(n_pooled, real64)
                if (s_old > 0.0_real64) then
                    rel = abs(s_new / s_old - 1.0_real64)
                else if (s_new > 0.0_real64) then
                    rel = huge(1.0_real64)
                else
                    rel = 0.0_real64
                end if
                if (rel > tau_stab) then
                    n_pooled = pool0
                    n_genes_pool = ngp0
                    abs_sum = abs0
                    radius = rad0
                    stop_reason = STOP_STABILITY
                    exit
                end if
            end if
            ! 2. a stop met inside the block
            if (pending /= 0) then
                stop_reason = pending
                exit
            end if
            ! 3. acceptance rate over the last w_blocks blocks
            n_blocks = n_blocks + 1
            slot = mod(n_blocks - 1, w_blocks) + 1
            hist_acc(slot) = block_acc
            hist_vis(slot) = block_vis
            if (n_blocks >= w_blocks .and. sum(hist_vis) > 0) then
                if (real(sum(hist_acc), real64) / real(sum(hist_vis), real64) < a_min) then
                    stop_reason = STOP_ACCEPTANCE
                    exit
                end if
            end if
        end do
    end subroutine gather_residuals_adaptive_helper

    ! =========================================================================
    ! 2-D ellipse gather (note F §5 open alternative)
    ! =========================================================================

    !> The `k` nearest genes in the joint distance
    !| `D_j = sqrt((dm_j / r_max)^2 + (du_j / tau_var)^2)` with
    !| `dm_j = |log2(m_j + 1) - log2(m_g + 1)|` (this side's linear means) and
    !| `du_j = |u_j - u_g|`, among the genes inside the ellipse `D_j <= 1`.
    !|
    !| Exact without scanning every gene: candidates are visited outward in log2 mean
    !| (two pointers on the sorted means, the nearer side first), so `dm` never
    !| decreases. Every unvisited gene has `D >= dm / r_max`, so the walk stops once
    !| `dm / r_max` exceeds the `k`-th best `D` (the `k` nearest are found,
    !| `STOP_K_MAX`) or exceeds 1 (nothing left inside the ellipse; `STOP_RADIUS`,
    !| a short pool when fewer than `k` were found), at the scan cap
    !| (`STOP_SCAN_CAP`) or when both sides are exhausted (`STOP_EXHAUSTED`). The
    !| selected genes' residuals are pooled in increasing `D`; the target itself
    !| (`D = 0`) comes first. `radius` = the largest `dm` among the selected genes.
    pure subroutine gather_residuals_ellipse_helper(target_mean, target_u, sorted_data, match_u, &
                                                    tau_var, r_max, k, scan_max, &
                                                    pooled_residuals, n_pooled, n_genes_pool, n_scanned, &
                                                    radius, stop_reason, max_pool_size)
        real(real64), intent(in) :: target_mean
        !! Mean of the target gene on this side (linear)
        real(real64), intent(in) :: target_u
        !! Matching variable of the target gene
        type(sorted_data_t), intent(in) :: sorted_data
        !! Sorted structure of this side
        real(real64), dimension(:), intent(in) :: match_u
        !! Matching variable, indexed by ORIGINAL gene index
        real(real64), intent(in) :: tau_var
        !! Semi-axis in the matching variable (the band half-width)
        real(real64), intent(in) :: r_max
        !! Semi-axis in log2 mean units
        integer(int32), intent(in) :: k
        !! Genes wanted
        integer(int32), intent(in) :: scan_max
        !! Maximum number of candidate slots examined
        real(real64), intent(out) :: pooled_residuals(:)
        !! Output residual pool (at least `max_pool_size` long)
        integer(int32), intent(out) :: n_pooled
        !! Residuals written
        integer(int32), intent(out) :: n_genes_pool
        !! Genes pooled
        integer(int32), intent(out) :: n_scanned
        !! Candidate slots examined
        real(real64), intent(out) :: radius
        !! Largest log2 mean distance among the selected genes
        integer(int32), intent(out) :: stop_reason
        !! `STOP_*` code
        integer(int32), intent(in) :: max_pool_size
        !! Hard cap on the residual count

        real(real64) :: best_d(k), best_dm(k), lt, dm_l, dm_r, dm, du, d, u_j
        integer(int32) :: best_idx(k), n_best, left, right, idx, pos, j, new_size
        real(real64), parameter :: LOG2F = 1.0_real64 / log(2.0_real64)

        n_pooled = 0
        n_genes_pool = 0
        n_scanned = 0
        radius = 0.0_real64
        stop_reason = STOP_EXHAUSTED
        n_best = 0

        pos = find_closest_helper(target_mean, sorted_data%means_sorted, sorted_data%n_genes)
        if (pos == 0) return
        lt = log(max(target_mean, 0.0_real64) + 1.0_real64) * LOG2F
        left = pos
        right = pos + 1

        do
            if (left < 1 .and. right > sorted_data%n_genes) then
                stop_reason = STOP_EXHAUSTED
                exit
            end if
            dm_l = huge(1.0_real64)
            dm_r = huge(1.0_real64)
            if (left >= 1) dm_l = abs(log(max(sorted_data%means_sorted(left), 0.0_real64) + 1.0_real64) * LOG2F - lt)
            if (right <= sorted_data%n_genes) &
                dm_r = abs(log(max(sorted_data%means_sorted(right), 0.0_real64) + 1.0_real64) * LOG2F - lt)
            if (dm_l <= dm_r) then
                idx = left
                dm = dm_l
                left = left - 1
            else
                idx = right
                dm = dm_r
                right = right + 1
            end if
            ! every unvisited gene has D >= dm / r_max
            if (dm / r_max > 1.0_real64) then
                stop_reason = merge(STOP_K_MAX, STOP_RADIUS, n_best >= k)
                exit
            end if
            if (n_best >= k) then
                if (dm / r_max >= best_d(k)) then
                    stop_reason = STOP_K_MAX
                    exit
                end if
            end if
            if (n_scanned >= scan_max) then
                stop_reason = STOP_SCAN_CAP
                exit
            end if
            n_scanned = n_scanned + 1
            u_j = match_u(sorted_data%original_indices(idx))
            if (.not. ieee_is_finite(u_j)) cycle
            du = abs(u_j - target_u)
            d = sqrt((dm / r_max)**2 + (du / tau_var)**2)
            if (d > 1.0_real64) cycle
            if (n_best >= k) then
                if (d >= best_d(k)) cycle
            else
                n_best = n_best + 1
            end if
            ! insert into the ascending list (the last slot is free or the worst)
            j = n_best
            do while (j > 1)
                if (best_d(j - 1) <= d) exit
                best_d(j) = best_d(j - 1)
                best_dm(j) = best_dm(j - 1)
                best_idx(j) = best_idx(j - 1)
                j = j - 1
            end do
            best_d(j) = d
            best_dm(j) = dm
            best_idx(j) = idx
        end do

        do j = 1, n_best
            if (n_pooled >= max_pool_size) exit
            idx = best_idx(j)
            new_size = min(n_pooled + sorted_data%n_residuals(idx), max_pool_size)
            call add_residuals_to_pool_helper(pooled_residuals, n_pooled, sorted_data%residuals_packed(:, idx), &
                                              sorted_data%n_residuals(idx), new_size)
            n_pooled = new_size
            n_genes_pool = n_genes_pool + 1
            radius = max(radius, best_dm(j))
        end do
    end subroutine gather_residuals_ellipse_helper

    ! =========================================================================
    ! diagnostics: neighbour genes and expression spread
    ! =========================================================================

    !> The genes of a mean-mode pool, in the order `gather_residuals_helper` took
    !| them (nearest mean first, own gene first). Replays its walk with the same
    !| `find_closest_helper` / `choose_index`; the number of genes follows from the
    !| pool size (every gene contributes `max_resid_per_gene` residuals; a pool cut
    !| by `max_pool_size` counts its last gene). Writes original gene indices to
    !| `ids(1:min(n, size(ids)))`; the rest is left untouched.
    pure subroutine neighbour_genes_helper(target_mean, sorted_data, n_pool, ids)
        real(real64), intent(in) :: target_mean
        !! Mean the pool was gathered around
        type(sorted_data_t), intent(in) :: sorted_data
        !! Sorted structure of the side
        integer(int32), intent(in) :: n_pool
        !! Residuals in the pool
        integer(int32), intent(inout) :: ids(:)
        !! Output: original gene indices
        integer(int32) :: n, j, pos, idx, left_cand, right_cand

        n = min(int(size(ids), int32), (n_pool + sorted_data%max_resid_per_gene - 1) / sorted_data%max_resid_per_gene)
        pos = find_closest_helper(target_mean, sorted_data%means_sorted, sorted_data%n_genes)
        if (pos == 0 .or. n < 1) return
        ids(1) = sorted_data%original_indices(pos)
        left_cand = pos - 1
        right_cand = pos + 1
        do j = 2, n
            call choose_index(sorted_data%means_sorted, sorted_data%n_genes, target_mean, idx, left_cand, right_cand)
            ids(j) = sorted_data%original_indices(idx)
        end do
    end subroutine neighbour_genes_helper

    !> Expression spread of a gene relative to its neighbourhood:
    !| `(max(own) - min(own)) / (Q95(pool) - Q05(pool))`, nearest-rank quantiles,
    !| -1 when the pool has no spread. `own` and `pool` must be on the same scale.
    subroutine spread_helper(own, n_own, pool, n_pool, perm, spread)
        integer(int32), intent(in) :: n_own
        !! Own residuals
        real(real64), intent(in) :: own(n_own)
        !! The gene's residuals (both sides)
        integer(int32), intent(in) :: n_pool
        !! Pool residuals
        real(real64), intent(in) :: pool(n_pool)
        !! Neighbourhood residuals (both sides)
        integer(int32), intent(inout) :: perm(:)
        !! Work array, size >= n_pool
        real(real64), intent(out) :: spread
        !! Output: spread, or -1
        integer(int32) :: i
        real(real64) :: q05, q95

        spread = -1.0_real64
        if (n_pool < 2 .or. n_own < 2) return
        perm(1:n_pool) = [(i, i = 1, n_pool)]
        call sort_real_heapsort_expl_size(pool, perm(1:n_pool), n_pool)
        q05 = pool(perm(1 + nint(0.05_real64 * real(n_pool - 1, real64))))
        q95 = pool(perm(1 + nint(0.95_real64 * real(n_pool - 1, real64))))
        if (q95 > q05) spread = (maxval(own) - minval(own)) / (q95 - q05)
    end subroutine spread_helper

    ! =========================================================================
    ! pipeline
    ! =========================================================================

    !> Validate, build the variance quantities and the sorted structures, apply the
    !| scaling, and run the per-gene loop of the variance-aware exact model.
    !|
    !| Means, the observed statistic (`observed_statistic_own`, centred on request via
    !| `beta_centre`), `log2fc` and `delta_hat` are computed from the replicates as in
    !| the scalar modules; `observed_statistic_own` is returned UNSCALED (the tested
    !| statistic under a scaling is `observed_statistic_own / sqrt(s2s_pool)`).
    !|
    !| A gene is not computable (p = -1) when either pool holds fewer than
    !| `MIN_POOL_RESIDUALS` residuals; for band / adaptive also when its `u` is
    !| non-finite or either side accepts fewer than `BAND_MIN_GENES` genes; for mean
    !| mode also when its statistic is non-finite (as `noise_model_exact`).
    !|
    !| Per-gene outputs that a mode does not produce are -1 (`n_scanned`,
    !| `n_accepted`, `radius`, `stop`); `s2s_pool` is 0 when no variance quantity is
    !| computed (`scaling = none` with mean mode or a raw band); `shrink_diag` is 0
    !| then. `match_u` is 0 in mean mode.
    !|
    !| Optional (all absent = the exact null, no diagnostics): `null_construction`
    !| (`NULL_EXACT` / `NULL_BOOTSTRAP`) with `n_boot` draws (default
    !| `N_BOOTSTRAP_DRAWS`; the RNG is seeded with 42 per call, the gene loop is
    !| sequential); `expr_spread_own` (per gene, see `spread_helper`; own residuals
    !| and pools as the test uses them, before the per-side 1/sqrt(n); -1 if not
    !| computed); `neighbours_own_case` / `_control` (k_cap x n_genes, original gene
    !| indices, 0-padded; MEAN MODE ONLY, zeros in the other modes).
    subroutine compute_noise_pvalue_pipeline_variance( &
        replicates_case, n_genes_case, n_replicates_case, &
        replicates_control, n_genes_control, n_replicates_control, &
        compute_pvalue_own, beta_centre, n_genes, norm_method, &
        variance_mode, match_variable, scaling, shrink_target, lambda_scale, &
        k_start, k_step, k_max, tau, k, c_tau, c_js, &
        block_size, k_min, r_max, tau_stab, w_blocks, a_min, &
        loess_span, loess_degree, loess_robust, shrink_loess_robust, max_pool_size, &
        pvalues_own, observed_statistic_own, log2fc, delta_hat, n_genes_with_pvalue, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        n_scanned_own_case, n_scanned_own_control, &
        n_accepted_own_case, n_accepted_own_control, &
        radius_own_case, radius_own_control, stop_own_case, stop_own_control, &
        match_u, s2s_pool, tau_var, scatter, shrink_diag, n_floored, ierr, &
        lambda_fixed, target_fixed, null_construction, n_boot, expr_spread_own, &
        neighbours_own_case, neighbours_own_control)

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
        !! 0 = linear residuals/statistic/variances; non-zero = log2. Note A's raw
        !! matching variable is always on the log2 scale.
        integer(int32), intent(in) :: variance_mode
        !! `VARIANCE_MODE_MEAN` (0), `_BAND` (1), `_ADAPTIVE` (3) or `_ELLIPSE` (4)
        integer(int32), intent(in) :: match_variable
        !! `MATCH_RAW` (0) or `MATCH_JS` (1); used by band and adaptive
        integer(int32), intent(in) :: scaling
        !! `SCALING_NONE` (0), `_XFIT_RAW` (1), `_XFIT_JS` (2), `_SAME_GROUP` (3)
        integer(int32), intent(in) :: shrink_target
        !! `SHRINK_TARGET_MEDIAN` (0) or `_LOESS` (1)
        integer(int32), intent(in) :: lambda_scale
        !! `LAMBDA_RELATIVE` (0) or `LAMBDA_ABSOLUTE` (1)
        integer(int32), intent(in) :: k_start
        !! Mean mode: initial neighbour genes (exact model)
        integer(int32), intent(in) :: k_step
        !! Mean mode: genes added per adaptive step
        integer(int32), intent(in) :: k_max
        !! Mean mode: maximum neighbour genes
        real(real64), intent(in) :: tau
        !! Mean mode: relative pool-scale increase that stops the growth
        integer(int32), intent(in) :: k
        !! Band: fixed accepted neighbour genes; adaptive: maximum accepted genes
        real(real64), intent(in) :: c_tau
        !! Raw band constant (> 0): `tau_var = c_tau * sqrt(psi1(d/2))`
        real(real64), intent(in) :: c_js
        !! JS band constant (> 0): `tau_var = c_js * MAD(u)`
        integer(int32), intent(in) :: block_size
        !! Adaptive: slots per block (>= 1)
        integer(int32), intent(in) :: k_min
        !! Adaptive: accepted genes before the stability stop applies (1..k)
        real(real64), intent(in) :: r_max
        !! Adaptive: growth radius; ellipse: semi-axis; log2 mean units (> 0, may be +Inf)
        real(real64), intent(in) :: tau_stab
        !! Adaptive: stability tolerance (> 0, may be +Inf)
        integer(int32), intent(in) :: w_blocks
        !! Adaptive: blocks in the acceptance-rate window (>= 1)
        real(real64), intent(in) :: a_min
        !! Adaptive: minimum acceptance rate, `[0, 1]`
        real(real64), intent(in) :: loess_span
        !! LOESS span, `(0, 1]` (note A trend and LOESS target)
        integer(int32), intent(in) :: loess_degree
        !! LOESS degree, 0..2
        integer(int32), intent(in) :: loess_robust
        !! Note A's trend: 1 = robust (`family = "symmetric"`), 0 = plain
        integer(int32), intent(in) :: shrink_loess_robust
        !! LOESS shrinkage target: 1 = robust, 0 = plain (note E v2 default)
        integer(int32), intent(in) :: max_pool_size
        !! Hard cap on residuals per pool
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: observed_statistic_own
        !! Output: observed statistic (centred when `beta_centre = 1`), unscaled
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
        !! Output: case candidate slots examined (band/adaptive; -1 otherwise)
        integer(int32), dimension(n_genes), intent(out) :: n_scanned_own_control
        !! Output: control candidate slots examined (band/adaptive; -1 otherwise)
        integer(int32), dimension(n_genes), intent(out) :: n_accepted_own_case
        !! Output: case genes accepted (band/adaptive; -1 otherwise)
        integer(int32), dimension(n_genes), intent(out) :: n_accepted_own_control
        !! Output: control genes accepted (band/adaptive; -1 otherwise)
        real(real64), dimension(n_genes), intent(out) :: radius_own_case
        !! Output: adaptive / ellipse case radius (log2 mean units; -1 otherwise)
        real(real64), dimension(n_genes), intent(out) :: radius_own_control
        !! Output: adaptive control radius (-1 otherwise)
        integer(int32), dimension(n_genes), intent(out) :: stop_own_case
        !! Output: adaptive / ellipse case stop reason `STOP_*` (-1 otherwise)
        integer(int32), dimension(n_genes), intent(out) :: stop_own_control
        !! Output: adaptive control stop reason (-1 otherwise)
        real(real64), dimension(n_genes), intent(out) :: match_u
        !! Output: matching variable `u` (band/adaptive)
        real(real64), dimension(n_genes), intent(out) :: s2s_pool
        !! Output: pooled variance used to scale the statistic (shrunk for JS)
        real(real64), intent(out) :: tau_var
        !! Output: band half-width actually used
        real(real64), intent(out) :: scatter
        !! Output: MAD^2 of `u` (band/adaptive)
        real(real64), dimension(N_SHRINK_DIAG), intent(out) :: shrink_diag
        !! Output: shrinkage diagnostics (see `compute_moment_shrinkage`)
        integer(int32), intent(out) :: n_floored
        !! Output: genes whose variance hit `VAR_FLOOR` (raw band or shrinkage)
        integer(int32), intent(out) :: ierr
        !! Error code
        real(real64), intent(in), optional :: lambda_fixed
        !! Tests only: fixed lambda for all shrinkages
        real(real64), intent(in), optional :: target_fixed
        !! Tests only: fixed constant shrinkage target
        integer(int32), intent(in), optional :: null_construction
        !! `NULL_EXACT` (default) or `NULL_BOOTSTRAP`
        integer(int32), intent(in), optional :: n_boot
        !! Bootstrap draws per gene (>= 1; default `N_BOOTSTRAP_DRAWS`)
        real(real64), dimension(:), intent(out), optional :: expr_spread_own
        !! Output: expression spread per gene (size n_genes; -1 if not computed)
        integer(int32), dimension(:, :), intent(out), optional :: neighbours_own_case
        !! Output: case neighbour genes (k_cap x n_genes; mean mode; 0-padded)
        integer(int32), dimension(:, :), intent(out), optional :: neighbours_own_control
        !! Output: control neighbour genes (k_cap x n_genes; mean mode; 0-padded)

        type(sorted_data_t) :: sorted_case, sorted_control
        real(real64), allocatable :: means_case(:), means_control(:), centres_control(:)
        real(real64), allocatable :: med_work(:), pool_case(:), pool_control(:)
        real(real64), allocatable :: s2_c(:), s2_t(:), s2(:), m(:), target(:), s2s_c(:), s2s_t(:)
        real(real64), allocatable :: scale_case_gene(:), scale_control_gene(:), obs_test(:)
        integer(int32), allocatable :: med_perm(:), med_left(:), med_right(:)
        integer(int32) :: i_gene, n_finite, sort_ierr, scan_max, idx, nf_shrink, nf_raw
        integer(int32) :: n_pool_case, n_pool_control, n_gp_case, n_gp_control
        real(real64) :: scale_case, scale_control, mad
        logical :: need_vars, band_like
        integer(int32) :: null_c, nb, per_iter
        real(real64), allocatable :: rbuf(:), own_work(:), spread_work(:)
        integer(int32), allocatable :: pos_case(:), pos_control(:), spread_perm(:)

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
        n_accepted_own_case = -1
        n_accepted_own_control = -1
        radius_own_case = -1.0_real64
        radius_own_control = -1.0_real64
        stop_own_case = -1
        stop_own_control = -1
        match_u = 0.0_real64
        s2s_pool = 0.0_real64
        tau_var = 0.0_real64
        scatter = 0.0_real64
        shrink_diag = 0.0_real64
        n_floored = 0
        null_c = NULL_EXACT
        if (present(null_construction)) null_c = null_construction
        nb = N_BOOTSTRAP_DRAWS
        if (present(n_boot)) nb = n_boot
        if (present(expr_spread_own)) expr_spread_own = -1.0_real64
        if (present(neighbours_own_case)) neighbours_own_case = 0
        if (present(neighbours_own_control)) neighbours_own_control = 0

        call validate_dimension_size(n_genes_case, ierr)
        call validate_dimension_size(n_replicates_case, ierr)
        call validate_dimension_size(n_genes_control, ierr)
        call validate_dimension_size(n_replicates_control, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(k, ierr)
        call validate_dimension_size(k_start, ierr)
        call validate_dimension_size(k_step, ierr)
        call validate_dimension_size(k_max, ierr)
        call validate_dimension_size(block_size, ierr)
        call validate_dimension_size(k_min, ierr)
        call validate_dimension_size(w_blocks, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        if (variance_mode /= VARIANCE_MODE_MEAN .and. variance_mode /= VARIANCE_MODE_BAND .and. &
            variance_mode /= VARIANCE_MODE_ADAPTIVE .and. variance_mode /= VARIANCE_MODE_ELLIPSE) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (match_variable /= MATCH_RAW .and. match_variable /= MATCH_JS) call set_err(ierr, ERR_INVALID_INPUT)
        if (scaling < SCALING_NONE .or. scaling > SCALING_SAME_GROUP) call set_err(ierr, ERR_INVALID_INPUT)
        if (shrink_target /= SHRINK_TARGET_MEDIAN .and. shrink_target /= SHRINK_TARGET_LOESS) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (lambda_scale /= LAMBDA_RELATIVE .and. lambda_scale /= LAMBDA_ABSOLUTE) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (beta_centre /= BETA_CENTRE_NONE .and. beta_centre /= BETA_CENTRE_MEDIAN) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (.not. (c_tau > 0.0_real64) .or. .not. (c_js > 0.0_real64)) call set_err(ierr, ERR_INVALID_INPUT)
        if (.not. (tau >= 0.0_real64)) call set_err(ierr, ERR_INVALID_INPUT)
        if (k_min > k) call set_err(ierr, ERR_INVALID_INPUT)
        if (.not. (r_max > 0.0_real64) .or. .not. (tau_stab > 0.0_real64)) call set_err(ierr, ERR_INVALID_INPUT)
        if (.not. (a_min >= 0.0_real64 .and. a_min <= 1.0_real64)) call set_err(ierr, ERR_INVALID_INPUT)
        if (loess_robust /= 0 .and. loess_robust /= 1) call set_err(ierr, ERR_INVALID_INPUT)
        if (shrink_loess_robust /= 0 .and. shrink_loess_robust /= 1) call set_err(ierr, ERR_INVALID_INPUT)
        if (null_c /= NULL_EXACT .and. null_c /= NULL_BOOTSTRAP) call set_err(ierr, ERR_INVALID_INPUT)
        call validate_dimension_size(nb, ierr)
        if (present(expr_spread_own)) then
            if (size(expr_spread_own) /= n_genes) call set_err(ierr, ERR_DIM_MISMATCH)
        end if
        if (present(neighbours_own_case)) then
            if (size(neighbours_own_case, 2) /= n_genes) call set_err(ierr, ERR_DIM_MISMATCH)
        end if
        if (present(neighbours_own_control)) then
            if (size(neighbours_own_control, 2) /= n_genes) call set_err(ierr, ERR_DIM_MISMATCH)
        end if
        if (is_err(ierr)) return
        ! One gene index for both matrices and every per-gene array (see eval doc §18).
        if (n_genes_case /= n_genes .or. n_genes_control /= n_genes) then
            call set_err(ierr, ERR_DIM_MISMATCH)
            return
        end if
        call validate_all_in_range_real(replicates_case, n_replicates_case * n_genes_case, ierr)
        call validate_all_in_range_real(replicates_control, n_replicates_control * n_genes_control, ierr)
        if (is_err(ierr)) return
        ! prepare_sorted_data rejects n < 2 per side; check it here first so no
        ! variance step divides by zero df.
        if (n_replicates_case < 2 .or. n_replicates_control < 2) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        ! --- means, statistic, effect size (as in the scalar modules) -------------
        M_ALLOCATE(means_case(n_genes))
        M_ALLOCATE(means_control(n_genes))
        M_ALLOCATE(centres_control(n_genes))
        M_ALLOCATE(med_work(n_genes))
        M_ALLOCATE(med_perm(n_genes))
        M_ALLOCATE(med_left(n_genes))
        M_ALLOCATE(med_right(n_genes))
        M_ALLOCATE(obs_test(n_genes))

        call compute_gene_centres_helper(replicates_case, n_replicates_case, n_genes, norm_method, &
                                         means_case, observed_statistic_own)
        call compute_gene_centres_helper(replicates_control, n_replicates_control, n_genes, norm_method, &
                                         means_control, centres_control)
        observed_statistic_own = observed_statistic_own - centres_control
        call compute_log2fc_helper(means_case, means_control, n_genes, log2fc)
        call median_of_finite_helper(observed_statistic_own, n_genes, med_work, med_perm, &
                                     med_left, med_right, delta_hat, n_finite)
        if (beta_centre == BETA_CENTRE_MEDIAN) observed_statistic_own = observed_statistic_own - delta_hat
        obs_test = observed_statistic_own

        band_like = (variance_mode /= VARIANCE_MODE_MEAN)
        need_vars = (scaling /= SCALING_NONE) .or. (band_like .and. match_variable == MATCH_JS)

        ! --- variances and shrinkage (note F §2) ----------------------------------
        if (need_vars) then
            M_ALLOCATE(s2_c(n_genes))
            M_ALLOCATE(s2_t(n_genes))
            M_ALLOCATE(s2(n_genes))
            M_ALLOCATE(m(n_genes))
            M_ALLOCATE(target(n_genes))
            M_ALLOCATE(s2s_c(n_genes))
            M_ALLOCATE(s2s_t(n_genes))
            call compute_group_variances_helper(replicates_case, n_replicates_case, replicates_control, &
                                                n_replicates_control, n_genes, norm_method, s2_c, s2_t, s2, m)
            call compute_moment_shrinkage(s2_c, s2_t, s2, m, n_genes, n_replicates_case, n_replicates_control, &
                                          shrink_target, lambda_scale, loess_span, loess_degree, &
                                          shrink_loess_robust, target, s2s_c, s2s_t, s2s_pool, shrink_diag, &
                                          nf_shrink, ierr, lambda_fixed, target_fixed)
            if (is_err(ierr)) return
            n_floored = nf_shrink
            ! JS matching variable from the SHRUNK pooled variance (before the
            ! xfit_raw override below, which only changes the scaling).
            if (band_like .and. match_variable == MATCH_JS) then
                match_u = log(s2s_pool / target)
                call mad_helper(match_u, n_genes, mad, med_work, med_perm, med_left, med_right)
                scatter = mad**2
                tau_var = c_js * mad
            end if
            if (scaling == SCALING_XFIT_RAW) then
                ! lambda = 0: the raw variances, floored (shrink_diag keeps the estimates)
                n_floored = count(s2 < VAR_FLOOR .or. s2_c < VAR_FLOOR .or. s2_t < VAR_FLOOR)
                s2s_pool = max(s2, VAR_FLOOR)
                s2s_c = max(s2_c, VAR_FLOOR)
                s2s_t = max(s2_t, VAR_FLOOR)
            end if
        end if

        ! --- raw matching variable (note A; the JS one is built above) -------------
        if (band_like .and. match_variable == MATCH_RAW) then
            call compute_variance_match(replicates_case, n_replicates_case, replicates_control, &
                                        n_replicates_control, n_genes, c_tau, loess_span, loess_degree, &
                                        loess_robust, match_u, tau_var, scatter, nf_raw, ierr)
            if (is_err(ierr)) return
            n_floored = max(n_floored, nf_raw)
        end if

        ! --- sorted structures (the exact module's) -------------------------------
        call prepare_sorted_data(means_case, replicates_case, n_replicates_case, n_genes, norm_method, &
                                 sorted_case, sort_ierr)
        call set_err(ierr, sort_ierr)
        call prepare_sorted_data(means_control, replicates_control, n_replicates_control, n_genes, norm_method, &
                                 sorted_control, sort_ierr)
        call set_err(ierr, sort_ierr)
        if (is_err(ierr)) return

        ! --- scaling of the null and of the statistic (note F §3) ------------------
        ! Applied once to the packed residuals: the case residuals of gene j are
        ! divided by sqrt of its CONTROL variance and vice versa (cross-fitting), or
        ! both by its own pooled shrunk variance (same_group, negative control).
        if (scaling /= SCALING_NONE) then
            M_ALLOCATE(scale_case_gene(n_genes))
            M_ALLOCATE(scale_control_gene(n_genes))
            if (scaling == SCALING_SAME_GROUP) then
                scale_case_gene = 1.0_real64 / sqrt(s2s_pool)
                scale_control_gene = scale_case_gene
            else
                scale_case_gene = 1.0_real64 / sqrt(s2s_t)
                scale_control_gene = 1.0_real64 / sqrt(s2s_c)
            end if
            do idx = 1, n_genes
                sorted_case%residuals_packed(:, idx) = sorted_case%residuals_packed(:, idx) &
                    * scale_case_gene(sorted_case%original_indices(idx))
                sorted_control%residuals_packed(:, idx) = sorted_control%residuals_packed(:, idx) &
                    * scale_control_gene(sorted_control%original_indices(idx))
            end do
            obs_test = observed_statistic_own / sqrt(s2s_pool)
        end if

        M_ALLOCATE(pool_case(max_pool_size * 2))
        M_ALLOCATE(pool_control(max_pool_size * 2))

        if (null_c == NULL_BOOTSTRAP) then
            ! as the bootstrap model: seeded once, one buffer reused across genes
            call init_random(42_int32)
            per_iter = n_replicates_case + n_replicates_control
            M_ALLOCATE(rbuf(per_iter * max(1_int32, min(nb, RNG_BUFFER_MAX / per_iter))))
        end if
        if (present(expr_spread_own)) then
            ! sorted position of every gene, for its own residuals
            M_ALLOCATE(pos_case(n_genes))
            M_ALLOCATE(pos_control(n_genes))
            M_ALLOCATE(own_work(n_replicates_case + n_replicates_control))
            M_ALLOCATE(spread_work(max_pool_size * 4))
            M_ALLOCATE(spread_perm(max_pool_size * 4))
            do idx = 1, n_genes
                pos_case(sorted_case%original_indices(idx)) = idx
                pos_control(sorted_control%original_indices(idx)) = idx
            end do
        end if

        scan_max = int(min(int(BAND_MAX_SCAN_FACTOR, int64) * int(k, int64), int(huge(1_int32), int64)), int32)
        ! Per-side 1/sqrt(n): the noise_model_exact null, unchanged.
        scale_case = 1.0_real64 / sqrt(real(n_replicates_case, real64))
        scale_control = 1.0_real64 / sqrt(real(n_replicates_control, real64))

        do i_gene = 1, n_genes
            select case (variance_mode)
            case (VARIANCE_MODE_MEAN)
                ! exactly noise_model_exact's loop
                call gather_residuals_helper(means_case(i_gene), sorted_case, k_start, k_step, k_max, tau, &
                                             pool_case, n_pool_case, max_pool_size)
                if (n_pool_case < MIN_POOL_RESIDUALS) cycle
                call gather_residuals_helper(means_control(i_gene), sorted_control, k_start, k_step, k_max, tau, &
                                             pool_control, n_pool_control, max_pool_size)
                if (n_pool_control < MIN_POOL_RESIDUALS) cycle
                if (obs_test(i_gene) /= obs_test(i_gene)) cycle
            case (VARIANCE_MODE_BAND)
                if (.not. ieee_is_finite(match_u(i_gene))) cycle
                call gather_residuals_band_helper(means_case(i_gene), match_u(i_gene), sorted_case, match_u, &
                                                  k, tau_var, scan_max, &
                                                  pool_case, n_pool_case, n_gp_case, n_scanned_own_case(i_gene), &
                                                  max_pool_size)
                n_accepted_own_case(i_gene) = n_gp_case
                if (n_gp_case < BAND_MIN_GENES .or. n_pool_case < MIN_POOL_RESIDUALS) cycle
                call gather_residuals_band_helper(means_control(i_gene), match_u(i_gene), sorted_control, match_u, &
                                                  k, tau_var, scan_max, &
                                                  pool_control, n_pool_control, n_gp_control, &
                                                  n_scanned_own_control(i_gene), max_pool_size)
                n_accepted_own_control(i_gene) = n_gp_control
                if (n_gp_control < BAND_MIN_GENES .or. n_pool_control < MIN_POOL_RESIDUALS) cycle
            case (VARIANCE_MODE_ELLIPSE)
                if (.not. ieee_is_finite(match_u(i_gene))) cycle
                call gather_residuals_ellipse_helper(means_case(i_gene), match_u(i_gene), sorted_case, match_u, &
                                                     tau_var, r_max, k, scan_max, pool_case, n_pool_case, n_gp_case, &
                                                     n_scanned_own_case(i_gene), radius_own_case(i_gene), &
                                                     stop_own_case(i_gene), max_pool_size)
                n_accepted_own_case(i_gene) = n_gp_case
                if (n_gp_case < BAND_MIN_GENES .or. n_pool_case < MIN_POOL_RESIDUALS) cycle
                call gather_residuals_ellipse_helper(means_control(i_gene), match_u(i_gene), sorted_control, &
                                                     match_u, tau_var, r_max, k, scan_max, pool_control, &
                                                     n_pool_control, n_gp_control, n_scanned_own_control(i_gene), &
                                                     radius_own_control(i_gene), stop_own_control(i_gene), &
                                                     max_pool_size)
                n_accepted_own_control(i_gene) = n_gp_control
                if (n_gp_control < BAND_MIN_GENES .or. n_pool_control < MIN_POOL_RESIDUALS) cycle
            case default
                if (.not. ieee_is_finite(match_u(i_gene))) cycle
                call gather_residuals_adaptive_helper(means_case(i_gene), match_u(i_gene), sorted_case, match_u, &
                                                      tau_var, block_size, k_min, k, r_max, tau_stab, w_blocks, a_min, &
                                                      scan_max, pool_case, n_pool_case, n_gp_case, &
                                                      n_scanned_own_case(i_gene), radius_own_case(i_gene), &
                                                      stop_own_case(i_gene), max_pool_size)
                n_accepted_own_case(i_gene) = n_gp_case
                if (n_gp_case < BAND_MIN_GENES .or. n_pool_case < MIN_POOL_RESIDUALS) cycle
                call gather_residuals_adaptive_helper(means_control(i_gene), match_u(i_gene), sorted_control, &
                                                      match_u, tau_var, block_size, k_min, k, r_max, tau_stab, &
                                                      w_blocks, a_min, scan_max, pool_control, n_pool_control, &
                                                      n_gp_control, n_scanned_own_control(i_gene), &
                                                      radius_own_control(i_gene), stop_own_control(i_gene), &
                                                      max_pool_size)
                n_accepted_own_control(i_gene) = n_gp_control
                if (n_gp_control < BAND_MIN_GENES .or. n_pool_control < MIN_POOL_RESIDUALS) cycle
            end select

            if (compute_pvalue_own(i_gene) /= 1) cycle

            if (present(expr_spread_own)) then
                own_work(1:n_replicates_case) = sorted_case%residuals_packed(1:n_replicates_case, pos_case(i_gene))
                own_work(n_replicates_case + 1:) = &
                    sorted_control%residuals_packed(1:n_replicates_control, pos_control(i_gene))
                spread_work(1:n_pool_case) = pool_case(1:n_pool_case)
                spread_work(n_pool_case + 1:n_pool_case + n_pool_control) = pool_control(1:n_pool_control)
                call spread_helper(own_work, size(own_work, kind=int32), &
                                   spread_work(1:n_pool_case + n_pool_control), n_pool_case + n_pool_control, &
                                   spread_perm, expr_spread_own(i_gene))
            end if
            if (variance_mode == VARIANCE_MODE_MEAN) then
                if (present(neighbours_own_case)) &
                    call neighbour_genes_helper(means_case(i_gene), sorted_case, n_pool_case, &
                                                neighbours_own_case(:, i_gene))
                if (present(neighbours_own_control)) &
                    call neighbour_genes_helper(means_control(i_gene), sorted_control, n_pool_control, &
                                                neighbours_own_control(:, i_gene))
            end if

            if (null_c == NULL_BOOTSTRAP) then
                call compute_pvalue_bootstrap_mean_helper(pool_case(1:n_pool_case), n_pool_case, &
                                                          pool_control(1:n_pool_control), n_pool_control, &
                                                          n_replicates_case, n_replicates_control, &
                                                          abs(obs_test(i_gene)), nb, rbuf, pvalues_own(i_gene))
            else
                pool_case(1:n_pool_case) = pool_case(1:n_pool_case) * scale_case
                pool_control(1:n_pool_control) = pool_control(1:n_pool_control) * scale_control
                call compute_pvalue(pool_case(1:n_pool_case), n_pool_case, &
                                    pool_control(1:n_pool_control), n_pool_control, &
                                    obs_test(i_gene), pvalues_own(i_gene), ierr)
                if (is_err(ierr)) return
            end if

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
!| carries the variance-aware controls and diagnostics. Null-pointer checks only;
!| everything else is validated in Fortran. The test-only `lambda_fixed` /
!| `target_fixed` are not exposed. `null_construction` / `n_boot` select the null;
!| `expr_spread_own` is always filled; the neighbour arrays are
!| `max(1, neighbour_cap) x n_genes` and meaningful only when `neighbour_cap > 0`.
subroutine compute_noise_pvalues_pipeline_variance_c( &
    replicates_case, n_genes_case, n_replicates_case, &
    replicates_control, n_genes_control, n_replicates_control, &
    compute_pvalue_own, beta_centre, n_genes, norm_method, &
    variance_mode, match_variable, scaling, shrink_target, lambda_scale, &
    k_start, k_step, k_max, tau, k, c_tau, c_js, &
    block_size, k_min, r_max, tau_stab, w_blocks, a_min, &
    loess_span, loess_degree, loess_robust, shrink_loess_robust, max_pool_size, &
    pvalues_own, observed_statistic_own, log2fc, delta_hat, n_genes_with_pvalue, &
    neighborhood_size_own_case, neighborhood_size_own_control, &
    n_scanned_own_case, n_scanned_own_control, &
    n_accepted_own_case, n_accepted_own_control, &
    radius_own_case, radius_own_control, stop_own_case, stop_own_control, &
    match_u, s2s_pool, tau_var, scatter, shrink_diag, n_floored, &
    null_construction, n_boot, neighbour_cap, expr_spread_own, neighbours_own_case, neighbours_own_control, ierr) &
    bind(C, name="compute_noise_pvalues_pipeline_variance_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use noise_model_variance, only: compute_noise_pvalue_pipeline_variance, N_SHRINK_DIAG
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
    !! 0 = mean kNN; 1 = variance band; 3 = adaptive growth; 4 = 2-D ellipse
    integer(c_int), intent(in), target :: match_variable
    !! 0 = raw (note A); 1 = James-Stein
    integer(c_int), intent(in), target :: scaling
    !! 0 = none; 1 = cross-fit raw; 2 = cross-fit JS; 3 = same group (control)
    integer(c_int), intent(in), target :: shrink_target
    !! 0 = median; 1 = LOESS
    integer(c_int), intent(in), target :: lambda_scale
    !! 0 = relative weights; 1 = absolute
    integer(c_int), intent(in), target :: k_start
    !! Mean mode: initial neighbour genes
    integer(c_int), intent(in), target :: k_step
    !! Mean mode: genes per step
    integer(c_int), intent(in), target :: k_max
    !! Mean mode: maximum neighbour genes
    real(c_double), intent(in), target :: tau
    !! Mean mode: growth tolerance
    integer(c_int), intent(in), target :: k
    !! Band: fixed accepted genes; adaptive: maximum accepted genes
    real(c_double), intent(in), target :: c_tau
    !! Raw band constant
    real(c_double), intent(in), target :: c_js
    !! JS band constant
    integer(c_int), intent(in), target :: block_size
    !! Adaptive: slots per block
    integer(c_int), intent(in), target :: k_min
    !! Adaptive: genes before the stability stop applies
    real(c_double), intent(in), target :: r_max
    !! Adaptive: growth radius (log2 mean units)
    real(c_double), intent(in), target :: tau_stab
    !! Adaptive: stability tolerance
    integer(c_int), intent(in), target :: w_blocks
    !! Adaptive: acceptance-rate window in blocks
    real(c_double), intent(in), target :: a_min
    !! Adaptive: minimum acceptance rate
    real(c_double), intent(in), target :: loess_span
    !! LOESS span
    integer(c_int), intent(in), target :: loess_degree
    !! LOESS degree
    integer(c_int), intent(in), target :: loess_robust
    !! 1 = robust note-A trend LOESS
    integer(c_int), intent(in), target :: shrink_loess_robust
    !! 1 = robust LOESS shrinkage target
    integer(c_int), intent(in), target :: max_pool_size
    !! Hard cap on residuals per pool
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_own
    !! Output: p-values (-1 if not computed)
    real(c_double), dimension(n_genes), intent(out), target :: observed_statistic_own
    !! Output: observed statistic (unscaled)
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
    integer(c_int), dimension(n_genes), intent(out), target :: n_accepted_own_case
    !! Output: case genes accepted
    integer(c_int), dimension(n_genes), intent(out), target :: n_accepted_own_control
    !! Output: control genes accepted
    real(c_double), dimension(n_genes), intent(out), target :: radius_own_case
    !! Output: adaptive case radius
    real(c_double), dimension(n_genes), intent(out), target :: radius_own_control
    !! Output: adaptive control radius
    integer(c_int), dimension(n_genes), intent(out), target :: stop_own_case
    !! Output: adaptive case stop reason
    integer(c_int), dimension(n_genes), intent(out), target :: stop_own_control
    !! Output: adaptive control stop reason
    real(c_double), dimension(n_genes), intent(out), target :: match_u
    !! Output: matching variable
    real(c_double), dimension(n_genes), intent(out), target :: s2s_pool
    !! Output: pooled (shrunk) variance used for the statistic
    real(c_double), intent(out), target :: tau_var
    !! Output: band half-width
    real(c_double), intent(out), target :: scatter
    !! Output: MAD^2 of u
    real(c_double), dimension(N_SHRINK_DIAG), intent(out), target :: shrink_diag
    !! Output: lambdas, unclipped ratios, p_hat, participation ratio
    integer(c_int), intent(out), target :: n_floored
    !! Output: genes at the variance floor
    integer(c_int), intent(in), target :: null_construction
    !! 0 = exact pairwise null; 1 = bootstrap difference of means
    integer(c_int), intent(in), target :: n_boot
    !! Bootstrap draws per gene
    integer(c_int), intent(in), target :: neighbour_cap
    !! Neighbour genes exported per gene and side (0 = none)
    real(c_double), dimension(n_genes), intent(out), target :: expr_spread_own
    !! Output: expression spread per gene
    integer(c_int), dimension(max(1, neighbour_cap), n_genes), intent(out), target :: neighbours_own_case
    !! Output: case neighbour genes (mean mode)
    integer(c_int), dimension(max(1, neighbour_cap), n_genes), intent(out), target :: neighbours_own_control
    !! Output: control neighbour genes (mean mode)
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
    M_CHECK_NON_NULL(match_variable)
    M_CHECK_NON_NULL(scaling)
    M_CHECK_NON_NULL(shrink_target)
    M_CHECK_NON_NULL(lambda_scale)
    M_CHECK_NON_NULL(k_start)
    M_CHECK_NON_NULL(k_step)
    M_CHECK_NON_NULL(k_max)
    M_CHECK_NON_NULL(tau)
    M_CHECK_NON_NULL(k)
    M_CHECK_NON_NULL(c_tau)
    M_CHECK_NON_NULL(c_js)
    M_CHECK_NON_NULL(block_size)
    M_CHECK_NON_NULL(k_min)
    M_CHECK_NON_NULL(r_max)
    M_CHECK_NON_NULL(tau_stab)
    M_CHECK_NON_NULL(w_blocks)
    M_CHECK_NON_NULL(a_min)
    M_CHECK_NON_NULL(loess_span)
    M_CHECK_NON_NULL(loess_degree)
    M_CHECK_NON_NULL(loess_robust)
    M_CHECK_NON_NULL(shrink_loess_robust)
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
    M_CHECK_NON_NULL(n_accepted_own_case)
    M_CHECK_NON_NULL(n_accepted_own_control)
    M_CHECK_NON_NULL(radius_own_case)
    M_CHECK_NON_NULL(radius_own_control)
    M_CHECK_NON_NULL(stop_own_case)
    M_CHECK_NON_NULL(stop_own_control)
    M_CHECK_NON_NULL(match_u)
    M_CHECK_NON_NULL(s2s_pool)
    M_CHECK_NON_NULL(tau_var)
    M_CHECK_NON_NULL(scatter)
    M_CHECK_NON_NULL(shrink_diag)
    M_CHECK_NON_NULL(n_floored)
    M_CHECK_NON_NULL(null_construction)
    M_CHECK_NON_NULL(n_boot)
    M_CHECK_NON_NULL(neighbour_cap)
    M_CHECK_NON_NULL(expr_spread_own)
    M_CHECK_NON_NULL(neighbours_own_case)
    M_CHECK_NON_NULL(neighbours_own_control)

    call compute_noise_pvalue_pipeline_variance( &
        replicates_case, n_genes_case, n_replicates_case, &
        replicates_control, n_genes_control, n_replicates_control, &
        compute_pvalue_own, beta_centre, n_genes, norm_method, &
        variance_mode, match_variable, scaling, shrink_target, lambda_scale, &
        k_start, k_step, k_max, tau, k, c_tau, c_js, &
        block_size, k_min, r_max, tau_stab, w_blocks, a_min, &
        loess_span, loess_degree, loess_robust, shrink_loess_robust, max_pool_size, &
        pvalues_own, observed_statistic_own, log2fc, delta_hat, n_genes_with_pvalue, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        n_scanned_own_case, n_scanned_own_control, &
        n_accepted_own_case, n_accepted_own_control, &
        radius_own_case, radius_own_control, stop_own_case, stop_own_control, &
        match_u, s2s_pool, tau_var, scatter, shrink_diag, n_floored, ierr, &
        null_construction=null_construction, n_boot=n_boot, expr_spread_own=expr_spread_own, &
        neighbours_own_case=neighbours_own_case, neighbours_own_control=neighbours_own_control)
end subroutine compute_noise_pvalues_pipeline_variance_c
