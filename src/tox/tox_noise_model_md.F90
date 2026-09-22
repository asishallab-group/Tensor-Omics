#include <src/macros.h>

!> Multidimensional (multi-axis) noise model — Issue #185.
!|
!| Generalises the scalar noise model from a single case/control contrast to `d`
!| INDEPENDENT axes (tissues, stages, conditions). A gene is ONE VECTOR; the axes
!| are not tested individually inside this test.
!|
!| Per gene `g` and axis `a`:
!|
!|   beta_hat(a) = mean_case(g,a) - mean_ctrl(g,a)        (#185 §1; OLS, two groups)
!|   D           = sqrt( sum_a beta_hat(a)^2 )            effect size, log2FC units
!|   D_std       = sqrt( sum_a beta_hat(a)^2 / Var_null(a) )   test statistic, unitless
!|
!| Both norms are computed and both are calibrated, because the weights are applied
!| identically to the observed statistic and to every null draw. They answer
!| different questions — how large the response is, versus how confident we are that
!| there is one — exactly as `logFC` and `t` do in limma. `D_std` drives the p-value
!| and the ranking. With independent axes the null covariance is diagonal, so
!| `D_std` IS the full Mahalanobis statistic: no covariance matrix, no Cholesky, no
!| ridge.
!|
!| **The null is a product of residual vectors, not a bootstrap of means.** One
!| draw takes ONE residual per axis per side and scales each by `1/sqrt(n)` of ITS
!| OWN side:
!|
!|   for a = 1..d:
!|       beta_null(a) = eps_case(a)/sqrt(n_case(a)) - eps_ctrl(a)/sqrt(n_ctrl(a))
!|
!| The scaling is **per side, not per axis**. With ragged designs — say
!| `n_case = [4,6,7,6,5]` against `n_ctrl = [3,5,5,9,3]` — each of the ten numbers
!| enters its own term; a single per-axis factor would be right only when both sides
!| have equal `n`. The scaling is needed at all because the observed statistic is a
!| difference of MEANS over `n` replicates while a raw residual difference is at
!| single-observation scale: without it the null is inflated by roughly `sqrt(n)` and
!| the test is strongly conservative. This is the same correction
!| `noise_model_exact` already carries in 1-D.
!|
!| `Var_null(a)` then follows in closed form from the pools, with no sampling:
!|
!|   Var_null(a) = Var(pool_case(a))/n_case(a) + Var(pool_ctrl(a))/n_ctrl(a)
!|
!| which, once the pools are scaled in place, is just the sum of the two scaled
!| pools' variances. It is returned per gene per axis (`var_null`), so the weights
!| are inspectable rather than implicit.
!|
!| **Scope assumption, enforced upstream.** Every sample belongs to exactly one
!| axis — no plant, subject, extraction batch or control group shared between axes.
!| That is precisely the condition under which the Cartesian-product null is
!| correct: because the axes share no samples, drawing each axis independently
!| reproduces the true joint null. It is also why there is no `null_coupling`
!| argument. The assumption is metadata, not something the numbers reveal, so the R
!| wrapper carries the guard: it errors if any sample identifier appears on more
!| than one axis.
!|
!| **Three paths to the p-value**, decided per gene and REPORTED per gene in
!| `method_used`:
!|
!|   - ENUMERATED (`method_used = 0`). Walk the joint space `prod_a S_a`, where
!|     `S_a = n_pool_case(a) * n_pool_ctrl(a)` is every pairing of one case-side
!|     residual with one control-side residual, and count the weighted tail exactly.
!|     No floor. Taken while `prod_a S_a <= enum_max_product`. In practice this is
!|     the `d = 1` regression path and little else: `S_a` grows QUADRATICALLY in the
!|     neighbourhood size, so larger `k` makes enumeration worse, and at the usual
!|     `k = 20+` settings the joint space is already ~1e7 at `d = 2`. Keep the
!|     branch for what it is — at `d = 1` it reproduces `noise_model_exact`'s
!|     `pvalues_own` exactly.
!|   - SYSTEMATIC (`method_used = 1`, the default). Instead of drawing random
!|     indices, walk each axis's pool with a fixed stride chosen coprime to that
!|     pool's size, so the marginals stay balanced. Deterministic and reproducible,
!|     no RNG overhead, no seed dependence, parallelises without stream splitting,
!|     and lower variance than random draws at the same `M`. Caveat for the method
!|     text: `(1+B)/(M+1)` is justified as a *valid* Monte-Carlo p-value under
!|     RANDOM draws; under systematic sampling it is instead an approximation of the
!|     exact enumerated tail, and the error becomes a coverage error rather than a
!|     sampling error. Both variants are implemented so they can be compared on the
!|     mock null before systematic becomes the default in anger.
!|   - RNG (`method_used = 2`). The classical Monte-Carlo path, `sampling_mode = 1`.
!|
!| `method_used = 3` is reserved for the convolution null (exact tail without
!| enumeration, staged as a later phase); it is never produced here.
!|
!| **Resolution.** `M` (`n_draws_max`) is set by BH over `G` genes, not by `d`: we
!| estimate a one-dimensional tail probability, `SE(p) = sqrt(p(1-p)/M)`, with no
!| dimensional surcharge. It is a user parameter, not a compile-time constant;
!| `DEFAULT_N_DRAWS_MAX` documents the recommended 2e4.
!|
!| **Reuse, not copy.** `sorted_data_t`, `prepare_sorted_data`,
!| `gather_residuals_helper` and `trim_pool_tails_helper` come from `noise_model`
!| directly; the validated scalar modules are not modified.
module noise_model_md
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_NAN_INF, &
                          ERR_DIM_MISMATCH, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real
    use f42_utils, only: sort_real, init_random
    use noise_model, only: sorted_data_t, prepare_sorted_data, &
                           gather_residuals_helper, trim_pool_tails_helper, &
                           NULL_METHOD_POOLED, RNG_BUFFER_MAX, NOISE_LOG_OFFSET
    implicit none

    integer(int32), parameter :: BETA_MODE_INTERNAL = 0_int32
        !! `beta_mode`: compute `beta_hat` inside this module from the replicates.
        !! The default and the only scale-coherent option: the observed statistic is
        !! then built on exactly the scale the residual pools carry.

    integer(int32), parameter :: BETA_MODE_SUPPLIED = 1_int32
        !! `beta_mode`: use the caller's `beta_obs`. The internal estimate is still
        !! computed and reported against it through `beta_check_cor` /
        !! `beta_check_mad` — this module does NOT refuse a mismatch, it surfaces it.
        !! Correlation ~1 with a large MAD is a units problem (edgeR's coefficients
        !! are natural log, a factor of ln 2 against log2 residuals, and shrunk by
        !! default); degraded correlation is shrinkage or non-linearity.

    integer(int32), parameter :: BETA_CENTRE_NONE = 0_int32
        !! `beta_centre`: no composition centring.

    integer(int32), parameter :: BETA_CENTRE_MEDIAN = 1_int32
        !! `beta_centre`: subtract the per-axis MEDIAN of `beta_hat` across genes —
        !! median centring of log-ratios, the TPM-side analogue of what TMM does for
        !! counts (not TMM itself: TMM's weights come from count magnitudes, which do
        !! not exist on TPM). On a compositional scale `log2FC = beta - Delta` with
        !! `Delta` common to every gene, and the null — built from WITHIN-group
        !! residuals — cannot contain information about a BETWEEN-group offset.
        !! Running uncorrected is not the assumption-free option; it is the
        !! assumption `Delta = 0`. **Must be off whenever the R side has already
        !! normalised**, or the data are centred twice.

    integer(int32), parameter :: SAMPLING_SYSTEMATIC = 0_int32
        !! `sampling_mode`: walk each pool with a fixed coprime stride (default).

    integer(int32), parameter :: SAMPLING_RNG = 1_int32
        !! `sampling_mode`: draw indices from the global RNG stream.

    integer(int32), parameter :: METHOD_ENUMERATED = 0_int32
        !! `method_used`: exact weighted tail count over the whole joint space.
    integer(int32), parameter :: METHOD_SYSTEMATIC = 1_int32
        !! `method_used`: systematic (stride) sampling.
    integer(int32), parameter :: METHOD_RNG = 2_int32
        !! `method_used`: random Monte Carlo.
    integer(int32), parameter :: METHOD_CONVOLUTION = 3_int32
        !! `method_used`: reserved for the FFT-convolution null; never produced yet.
    integer(int32), parameter :: METHOD_NOT_COMPUTED = -1_int32
        !! `method_used`: the gene was not tested.

    integer(int32), parameter :: MIN_POOL_RESIDUALS = 10_int32
        !! Minimum residuals in an axis's pool before that axis is usable. Matches
        !! the scalar modules' gate. Not a live concern in practice —
        !! `k_start`/`k_step`/`k_max` count GENES and each gene contributes
        !! `n_replicates`, so 15 genes x 3 replicates already gives 45 — but the
        !! branch stays for safety. A gene whose pool falls below this on ANY axis is
        !! dropped ENTIRELY: dropping only the thin axis would make `D` incomparable
        !! between genes and break the magnitude ranking. The failing axis stays
        !! visible through `neighborhood_size_case` / `neighborhood_size_control`,
        !! which are filled for every axis that was gathered.

    integer(int32), parameter :: DEFAULT_N_DRAWS_MAX = 20000_int32
        !! Recommended `M`. Not enforced here (`n_draws_max` is a caller argument) —
        !! it documents the default the R wrapper applies. `M` is set by BH over the
        !! gene count, not by `d`: a gene at the floor `1/(M+1)` is rejected only if
        !! its rank reaches `p_floor*G/q`, so at `G = 20000`, `q = 0.05` roughly 20
        !! genes must be tied at the floor before any is rejected. `M = 4e5` is where
        !! the floor stops binding altogether; with sequential stopping that is close
        !! to free (see `n_exceed_target`).

    integer(int64), parameter :: DEFAULT_ENUM_MAX_PRODUCT = 100000_int64
        !! Recommended `enum_max_product`. Chosen so realistic `d = 1` pools always
        !! enumerate (pools of 316 per side) while `d >= 2` essentially never does —
        !! at `k = 15`, `n_rep = 3` the joint space is already 4.1e6 at `d = 2`.

contains

    ! =========================================================================
    ! compute_beta_md
    ! =========================================================================

    !> Core implementation: per-axis observed mean difference `beta_hat`, from the
    !| packed replicate matrices.
    !|
    !| **Why the replicates and not the means.** The plan sketches this routine as
    !| taking `means_case` / `means_control`. It cannot: on the log2 branch the
    !| statistic that is scale-coherent with the residual pools is the difference of
    !| FRECHET means, `mean_i log2(x_i + c)`, and by Jensen's inequality that is not
    !| recoverable from the linear-space mean `log2(mean + c)`. The scalar pipelines
    !| get their `obs_own` from the caller, who computes exactly
    !| `colMeans(log2(case+1)) - colMeans(log2(ctrl+1))`; reproducing that here
    !| (unit test 2) requires the replicates. On the linear branch the two
    !| definitions coincide, so nothing is lost.
    !|
    !| The `means` arrays keep their real job — the coordinate the kNN neighbourhood
    !| is matched in — and are not involved in the statistic.
    !|
    !| Matches `prepare_sorted_data_helper`'s centring expression term for term
    !| (`sum(log(max(x,0) + NOISE_LOG_OFFSET)) * log2_factor / n`), so the observed
    !| statistic and the residuals it is scored against are built from the same
    !| arithmetic.
    !|
    !| No allocation; `beta_obs` must be pre-shaped `(n_axes, n_genes)` by the caller.
    pure subroutine compute_beta_md_helper(replicates_case_packed, n_rep_case_per_axis, &
                                           replicates_control_packed, n_rep_control_per_axis, &
                                           n_genes, n_axes, norm_method, beta_obs)
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        integer(int32), intent(in) :: n_axes
        !! Number of axes `d`
        real(real64), dimension(:), intent(in) :: replicates_case_packed
        !! Case replicates, axis blocks back to back; block `a` is `(n_rep_case_a, n_genes)`
        integer(int32), dimension(n_axes), intent(in) :: n_rep_case_per_axis
        !! Case replicate count per axis (ragged by design)
        real(real64), dimension(:), intent(in) :: replicates_control_packed
        !! Control replicates, same packing
        integer(int32), dimension(n_axes), intent(in) :: n_rep_control_per_axis
        !! Control replicate count per axis
        integer(int32), intent(in) :: norm_method
        !! 0 = linear-scale means; non-zero = log2-space (Frechet) means
        real(real64), dimension(n_axes, n_genes), intent(out) :: beta_obs
        !! Output: observed per-axis mean difference, gene-contiguous

        integer(int32) :: i_axis, i_gene, i_rep, n_rep_c, n_rep_t, off_c, off_t, base_c, base_t
        real(real64) :: sum_c, sum_t, log2_factor
        logical :: use_log

        use_log = (norm_method /= 0)
        log2_factor = 1.0_real64 / log(2.0_real64)

        off_c = 0
        off_t = 0
        do i_axis = 1, n_axes
            n_rep_c = n_rep_case_per_axis(i_axis)
            n_rep_t = n_rep_control_per_axis(i_axis)

            do i_gene = 1, n_genes
                base_c = off_c + (i_gene - 1) * n_rep_c
                base_t = off_t + (i_gene - 1) * n_rep_t

                sum_c = 0.0_real64
                sum_t = 0.0_real64
                if (use_log) then
                    do i_rep = 1, n_rep_c
                        sum_c = sum_c + log(max(replicates_case_packed(base_c + i_rep), 0.0_real64) &
                                            + NOISE_LOG_OFFSET)
                    end do
                    do i_rep = 1, n_rep_t
                        sum_t = sum_t + log(max(replicates_control_packed(base_t + i_rep), 0.0_real64) &
                                            + NOISE_LOG_OFFSET)
                    end do
                    beta_obs(i_axis, i_gene) = sum_c * log2_factor / real(n_rep_c, real64) &
                                             - sum_t * log2_factor / real(n_rep_t, real64)
                else
                    do i_rep = 1, n_rep_c
                        sum_c = sum_c + replicates_case_packed(base_c + i_rep)
                    end do
                    do i_rep = 1, n_rep_t
                        sum_t = sum_t + replicates_control_packed(base_t + i_rep)
                    end do
                    beta_obs(i_axis, i_gene) = sum_c / real(n_rep_c, real64) &
                                             - sum_t / real(n_rep_t, real64)
                end if
            end do

            off_c = off_c + n_rep_c * n_genes
            off_t = off_t + n_rep_t * n_genes
        end do
    end subroutine compute_beta_md_helper

    !> Validate and compute the per-axis observed `beta_hat` (#185 §1).
    !|
    !| Public entry point around `compute_beta_md_helper`. Allocates nothing; it only
    !| adds argument validation, so it is safe to call standalone (the pipeline calls
    !| the helper directly, having already validated).
    subroutine compute_beta_md(replicates_case_packed, n_rep_case_per_axis, &
                               replicates_control_packed, n_rep_control_per_axis, &
                               n_genes, n_axes, norm_method, beta_obs, ierr)
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        integer(int32), intent(in) :: n_axes
        !! Number of axes `d`
        real(real64), dimension(:), intent(in) :: replicates_case_packed
        !! Case replicates, packed by axis
        integer(int32), dimension(n_axes), intent(in) :: n_rep_case_per_axis
        !! Case replicate count per axis
        real(real64), dimension(:), intent(in) :: replicates_control_packed
        !! Control replicates, packed by axis
        integer(int32), dimension(n_axes), intent(in) :: n_rep_control_per_axis
        !! Control replicate count per axis
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2 (Frechet) scale
        real(real64), dimension(n_axes, n_genes), intent(out) :: beta_obs
        !! Output: `(n_axes, n_genes)` observed mean differences
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: n_packed_case, n_packed_control

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_axes, ierr)
        if (is_err(ierr)) return

        call validate_axis_replicates_helper(n_rep_case_per_axis, n_rep_control_per_axis, &
                                             n_axes, n_genes, n_packed_case, n_packed_control, ierr)
        if (is_err(ierr)) return

        if (size(replicates_case_packed) < n_packed_case .or. &
            size(replicates_control_packed) < n_packed_control) then
            call set_err(ierr, ERR_DIM_MISMATCH)
            return
        end if

        call validate_all_in_range_real(replicates_case_packed(1:n_packed_case), n_packed_case, ierr)
        call validate_all_in_range_real(replicates_control_packed(1:n_packed_control), n_packed_control, ierr)
        if (is_err(ierr)) return

        call compute_beta_md_helper(replicates_case_packed, n_rep_case_per_axis, &
                                    replicates_control_packed, n_rep_control_per_axis, &
                                    n_genes, n_axes, norm_method, beta_obs)
    end subroutine compute_beta_md

    !> Check the per-axis replicate counts and derive the packed array lengths.
    !|
    !| Every axis needs at least 2 replicates PER SIDE: fewer leaves no residual
    !| variance to build a null from, and `prepare_sorted_data` rejects it anyway.
    !| Also guards the `int32` overflow of `sum_a n_rep_a * n_genes`, which is a real
    !| risk once `d` and `n_genes` are both large.
    pure subroutine validate_axis_replicates_helper(n_rep_case_per_axis, n_rep_control_per_axis, &
                                                     n_axes, n_genes, &
                                                     n_packed_case, n_packed_control, ierr)
        integer(int32), intent(in) :: n_axes
        !! Number of axes
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        integer(int32), dimension(n_axes), intent(in) :: n_rep_case_per_axis
        !! Case replicate count per axis
        integer(int32), dimension(n_axes), intent(in) :: n_rep_control_per_axis
        !! Control replicate count per axis
        integer(int32), intent(out) :: n_packed_case
        !! Derived total length of `replicates_case_packed`
        integer(int32), intent(out) :: n_packed_control
        !! Derived total length of `replicates_control_packed`
        integer(int32), intent(inout) :: ierr
        !! Error code

        integer(int32) :: i_axis
        integer(int64) :: total_case, total_control

        n_packed_case = 0
        n_packed_control = 0

        total_case = 0_int64
        total_control = 0_int64
        do i_axis = 1, n_axes
            if (n_rep_case_per_axis(i_axis) < 2 .or. n_rep_control_per_axis(i_axis) < 2) then
                call set_err(ierr, ERR_INVALID_INPUT)
                return
            end if
            total_case = total_case + int(n_rep_case_per_axis(i_axis), int64) * int(n_genes, int64)
            total_control = total_control + int(n_rep_control_per_axis(i_axis), int64) * int(n_genes, int64)
        end do

        if (total_case > int(huge(1_int32), int64) .or. total_control > int(huge(1_int32), int64)) then
            call set_err(ierr, ERR_DIM_MISMATCH)
            return
        end if

        n_packed_case = int(total_case, int32)
        n_packed_control = int(total_control, int32)
    end subroutine validate_axis_replicates_helper

    ! =========================================================================
    ! null-variance and sampling-stride helpers
    ! =========================================================================

    !> Variance of a uniform draw from `pool(1:n)` — the population variance.
    !|
    !| Not the `n-1` sample variance: a null draw picks one pool entry uniformly at
    !| random, so `mean(x^2) - mean(x)^2` is the EXACT variance of the distribution
    !| being drawn from, not an estimate of some wider population's. (The residuals
    !| already carry the Bessel correction from `prepare_sorted_data_helper`, which
    !| is where the `n-1` adjustment belongs.)
    !|
    !| Called on the ALREADY-SCALED pools, so the result is `Var(pool)/n_rep` and the
    !| two sides simply add to give `Var_null(a)` of the plan's closed form.
    !| Returns 0 for an empty pool.
    pure function pool_variance_helper(pool, n) result(var)
        integer(int32), intent(in) :: n
        !! Number of valid entries
        real(real64), dimension(n), intent(in) :: pool
        !! Residual pool
        real(real64) :: var
        !! Population variance of a uniform draw from the pool

        integer(int32) :: i
        real(real64) :: s, ss, rn

        if (n < 1) then
            var = 0.0_real64
            return
        end if

        s = 0.0_real64
        ss = 0.0_real64
        do i = 1, n
            s = s + pool(i)
            ss = ss + pool(i) * pool(i)
        end do
        rn = real(n, real64)
        var = max(0.0_real64, ss / rn - (s / rn) * (s / rn))
    end function pool_variance_helper

    !> Greatest common divisor, Euclid. Used only to pick sampling strides.
    pure function gcd_helper(a, b) result(g)
        integer(int32), intent(in) :: a
        !! First operand (non-negative)
        integer(int32), intent(in) :: b
        !! Second operand (non-negative)
        integer(int32) :: g
        !! Greatest common divisor

        integer(int32) :: x, y, t

        x = abs(a)
        y = abs(b)
        do while (y /= 0)
            t = modulo(x, y)
            x = y
            y = t
        end do
        g = x
    end function gcd_helper

    !> A stride for the systematic walk of a pool of `n` entries, coprime to `n`.
    !|
    !| Coprimality is what makes the walk visit every entry before repeating, so the
    !| marginal over that pool stays exactly balanced — the property that gives
    !| systematic sampling lower variance than random draws at the same `M`. The
    !| starting point is near `n/phi` (the golden-ratio fraction) so that successive
    !| draws are spread across the pool rather than adjacent, then walks up until the
    !| gcd is 1. Always terminates: 1 is coprime to everything and the search wraps
    !| within `[1, n]`.
    !|
    !| `salt` separates the strides of different axes and sides so the `d` walks do
    !| not move in lock step, which would collapse the product walk onto a diagonal.
    pure function coprime_stride_helper(n, salt) result(stride)
        integer(int32), intent(in) :: n
        !! Pool size
        integer(int32), intent(in) :: salt
        !! Per-axis/per-side offset, to decorrelate the walks
        integer(int32) :: stride
        !! A value in `[1, n]` coprime to `n`

        integer(int32) :: i

        if (n <= 1) then
            stride = 1
            return
        end if

        stride = modulo(int(real(n, real64) * 0.6180339887498949_real64, int32) + salt, n) + 1
        do i = 1, n
            if (gcd_helper(stride, n) == 1) return
            stride = modulo(stride, n) + 1
        end do
        stride = 1
    end function coprime_stride_helper

    ! =========================================================================
    ! per-gene multidimensional null
    ! =========================================================================

    !> p-value for one gene's `d`-dimensional statistic, by whichever of the three
    !| paths the joint outcome space allows.
    !|
    !| The pools handed in are ALREADY scaled by `1/sqrt(n)` of their own side, so
    !| one null draw is simply `pools_case(i, a) - pools_control(j, a)` per axis —
    !| the plan's residual-vector construction, with the per-side scaling folded into
    !| the pool once per gene instead of once per draw.
    !|
    !| The tail is counted on `D_std^2 = sum_a w_a * beta_a^2` rather than on
    !| `D_std`: `w_a >= 0` makes both sides non-negative, so `D_std_null >= D_std_obs`
    !| and its square are the same event, and it saves two square roots per draw. At
    !| `d = 1` this is `w*v*v >= w*b*b`, which is `|v| >= |b|` — the comparison
    !| `noise_model_exact` makes — since multiplying by a positive constant and
    !| squaring are both monotone under IEEE rounding; they can differ only in the
    !| measure-zero case where two distinct magnitudes square to the same double.
    !|
    !| `sum_d_sq` accumulates the UNWEIGHTED `D^2` so the caller can report
    !| `D2_adj = max(0, D^2 - mean(D^2_null))`; `D` is upward biased as a magnitude
    !| (`E[D] ~ sigma*sqrt(d) != 0` even under H0), and that is the correction.
    !|
    !| **Besag-Clifford sequential stopping** (`n_exceed_target > 0`, sampled paths
    !| only): stop as soon as `c` exceedances have been seen and report `p = c/m`
    !| instead of always paying `M` draws. Expected draws `~ min(M, c/p)`, so a 40x
    !| increase in `M` costs about 1.5x in average draws — which is what makes a
    !| large `M` affordable. It is a DIFFERENT estimator with granularity `c/m` and
    !| its interaction with BH is unverified, so it ships off by default and the
    !| fixed-`M` arm stays the reference.
    subroutine compute_pvalue_md_helper(pools_case, n_pool_case, &
                                        pools_control, n_pool_control, &
                                        n_axes, axis_w, acc_obs_std, &
                                        sampling_mode, n_draws_max, n_exceed_target, &
                                        enum_max_product, rbuf, gene_index, &
                                        p_value, d_sq_null_mean, n_draws_used, method_used)
        integer(int32), intent(in) :: n_axes
        !! Number of axes
        real(real64), dimension(:, :), intent(in) :: pools_case
        !! Case pools, column `a` holding this gene's axis-`a` pool, already scaled
        integer(int32), dimension(n_axes), intent(in) :: n_pool_case
        !! Valid length of each case pool column
        real(real64), dimension(:, :), intent(in) :: pools_control
        !! Control pools, same layout, already scaled
        integer(int32), dimension(n_axes), intent(in) :: n_pool_control
        !! Valid length of each control pool column
        real(real64), dimension(n_axes), intent(in) :: axis_w
        !! `1/Var_null(a)`, or 0 for a degenerate axis
        real(real64), intent(in) :: acc_obs_std
        !! Observed `D_std^2 = sum_a w_a * beta_a^2`
        integer(int32), intent(in) :: sampling_mode
        !! `SAMPLING_SYSTEMATIC` or `SAMPLING_RNG`; ignored on the enumerated path
        integer(int32), intent(in) :: n_draws_max
        !! `M`; p-value floor `1/(M+1)` on the sampled paths
        integer(int32), intent(in) :: n_exceed_target
        !! Besag-Clifford exceedance target `c`; 0 disables
        integer(int64), intent(in) :: enum_max_product
        !! Enumerate while `prod_a S_a <= this`; 0 forces sampling
        real(real64), dimension(:), intent(inout) :: rbuf
        !! Pre-allocated random-number buffer, reused across genes (RNG path only)
        integer(int32), intent(in) :: gene_index
        !! Gene index, used to offset the systematic walk so genes do not share a lattice
        real(real64), intent(out) :: p_value
        !! p-value
        real(real64), intent(out) :: d_sq_null_mean
        !! Mean of the UNWEIGHTED `D^2` over the draws actually taken
        integer(int32), intent(out) :: n_draws_used
        !! Draws (or enumerated atoms) actually taken
        integer(int32), intent(out) :: method_used
        !! `METHOD_ENUMERATED` / `METHOD_SYSTEMATIC` / `METHOD_RNG`

        integer(int32) :: i_axis, i_c, i_t, a, per_iter, chunk, j, buf_pos, buf_filled
        integer(int32) :: ic(n_axes), it(n_axes), stride_c(n_axes), stride_t(n_axes)
        integer(int64) :: n_total, s_axis(n_axes), idx(n_axes), count_ge, m
        real(real64) :: bsq(n_axes)
        real(real64) :: acc_std, acc_plain, v, sum_plain
        logical :: stop_early, use_rng

        ! Size of the joint outcome space: every pairing of one case-side residual
        ! with one control-side residual, per axis, multiplied out. Saturates rather
        ! than overflowing -- the product is only ever compared against a cap.
        n_total = 1_int64
        do i_axis = 1, n_axes
            s_axis(i_axis) = int(n_pool_case(i_axis), int64) * int(n_pool_control(i_axis), int64)
            if (s_axis(i_axis) <= 0_int64) then
                n_total = huge(1_int64)
            else if (n_total > huge(1_int64) / s_axis(i_axis)) then
                n_total = huge(1_int64)
            else
                n_total = n_total * s_axis(i_axis)
            end if
        end do

        count_ge = 0_int64
        m = 0_int64
        sum_plain = 0.0_real64

        if (enum_max_product > 0_int64 .and. n_total <= enum_max_product) then
            ! ---------------- ENUMERATED: exact weighted tail, no floor ----------
            ! Mixed-radix odometer over the d per-axis atom indices. Only the axes
            ! that actually rolled over are recomputed, so the walk costs O(1)
            ! amortised per atom plus the O(d) accumulation.
            idx = 0_int64
            do i_axis = 1, n_axes
                bsq(i_axis) = atom_bsq_helper(pools_case(:, i_axis), n_pool_case(i_axis), &
                                              pools_control(:, i_axis), n_pool_control(i_axis), &
                                              0_int64)
            end do

            outer: do
                acc_std = 0.0_real64
                acc_plain = 0.0_real64
                do i_axis = 1, n_axes
                    acc_std = acc_std + axis_w(i_axis) * bsq(i_axis)
                    acc_plain = acc_plain + bsq(i_axis)
                end do
                m = m + 1_int64
                sum_plain = sum_plain + acc_plain
                if (acc_std >= acc_obs_std) count_ge = count_ge + 1_int64

                a = 1
                do
                    idx(a) = idx(a) + 1_int64
                    if (idx(a) < s_axis(a)) then
                        bsq(a) = atom_bsq_helper(pools_case(:, a), n_pool_case(a), &
                                                 pools_control(:, a), n_pool_control(a), idx(a))
                        exit
                    end if
                    idx(a) = 0_int64
                    bsq(a) = atom_bsq_helper(pools_case(:, a), n_pool_case(a), &
                                             pools_control(:, a), n_pool_control(a), 0_int64)
                    a = a + 1
                    if (a > n_axes) exit outer
                end do
            end do outer

            p_value = real(count_ge + 1_int64, real64) / real(m + 1_int64, real64)
            method_used = METHOD_ENUMERATED
        else
            ! ---------------- SAMPLED: systematic stride, or RNG -----------------
            use_rng = (sampling_mode == SAMPLING_RNG)
            per_iter = max(1, 2 * n_axes)
            chunk = max(1, min(n_draws_max, int(size(rbuf), int32) / per_iter))
            buf_pos = 0
            buf_filled = 0

            if (.not. use_rng) then
                do i_axis = 1, n_axes
                    stride_c(i_axis) = coprime_stride_helper(n_pool_case(i_axis), 2 * i_axis)
                    stride_t(i_axis) = coprime_stride_helper(n_pool_control(i_axis), 2 * i_axis + 1)
                    ! Offset each gene's walk so genes do not all traverse the same
                    ! lattice positions in the same order.
                    ic(i_axis) = modulo(gene_index * (i_axis + 1), n_pool_case(i_axis))
                    it(i_axis) = modulo(gene_index * (i_axis + 2), n_pool_control(i_axis))
                end do
            end if

            stop_early = .false.
            do while (m < int(n_draws_max, int64) .and. .not. stop_early)
                if (use_rng .and. buf_pos >= buf_filled) then
                    j = min(chunk, n_draws_max - int(m, int32))
                    call random_number(rbuf(1:j * per_iter))
                    buf_filled = j * per_iter
                    buf_pos = 0
                end if

                acc_std = 0.0_real64
                acc_plain = 0.0_real64
                do i_axis = 1, n_axes
                    if (use_rng) then
                        i_c = min(int(rbuf(buf_pos + 1) * real(n_pool_case(i_axis), real64), int32) + 1, &
                                  n_pool_case(i_axis))
                        i_t = min(int(rbuf(buf_pos + 2) * real(n_pool_control(i_axis), real64), int32) + 1, &
                                  n_pool_control(i_axis))
                        buf_pos = buf_pos + 2
                    else
                        i_c = ic(i_axis) + 1
                        i_t = it(i_axis) + 1
                        ic(i_axis) = modulo(ic(i_axis) + stride_c(i_axis), n_pool_case(i_axis))
                        it(i_axis) = modulo(it(i_axis) + stride_t(i_axis), n_pool_control(i_axis))
                    end if
                    v = pools_case(i_c, i_axis) - pools_control(i_t, i_axis)
                    acc_std = acc_std + axis_w(i_axis) * v * v
                    acc_plain = acc_plain + v * v
                end do

                m = m + 1_int64
                sum_plain = sum_plain + acc_plain
                if (acc_std >= acc_obs_std) count_ge = count_ge + 1_int64
                if (n_exceed_target > 0 .and. count_ge >= int(n_exceed_target, int64)) stop_early = .true.
            end do

            if (stop_early) then
                ! Besag-Clifford: c exceedances in m draws.
                p_value = real(count_ge, real64) / real(m, real64)
            else
                p_value = real(count_ge + 1_int64, real64) / real(n_draws_max + 1, real64)
            end if
            method_used = merge(METHOD_RNG, METHOD_SYSTEMATIC, use_rng)
        end if

        n_draws_used = int(min(m, int(huge(1_int32), int64)), int32)
        if (m > 0_int64) then
            d_sq_null_mean = sum_plain / real(m, real64)
        else
            d_sq_null_mean = 0.0_real64
        end if
    end subroutine compute_pvalue_md_helper

    !> Squared null contribution of one axis at joint-space atom `t`.
    !|
    !| Atom `t` in `[0, n_case*n_control)` decodes to the pair
    !| `(t / n_control + 1, mod(t, n_control) + 1)`, so walking `t` walks every
    !| case/control residual pairing exactly once. Pools are already `1/sqrt(n)`
    !| scaled, so the difference is `beta_null(a)` directly.
    pure function atom_bsq_helper(pool_case, n_pool_case, pool_control, n_pool_control, t) result(bsq)
        integer(int32), intent(in) :: n_pool_case
        !! Case pool size
        real(real64), dimension(:), intent(in) :: pool_case
        !! Case pool column (only the first `n_pool_case` entries are valid)
        integer(int32), intent(in) :: n_pool_control
        !! Control pool size
        real(real64), dimension(:), intent(in) :: pool_control
        !! Control pool column
        integer(int64), intent(in) :: t
        !! Atom index in `[0, n_pool_case * n_pool_control)`
        real(real64) :: bsq
        !! `beta_null(a)^2` for that pairing

        integer(int32) :: i_c, i_t
        real(real64) :: v

        i_c = int(t / int(n_pool_control, int64), int32) + 1
        i_t = int(modulo(t, int(n_pool_control, int64)), int32) + 1
        v = pool_case(i_c) - pool_control(i_t)
        bsq = v * v
    end function atom_bsq_helper

    ! =========================================================================
    ! per-axis summaries: median, correlation
    ! =========================================================================

    !> Median of the finite entries of `values(1:n)`, by indirect sort.
    !|
    !| Even counts average the two central order statistics, matching R's `median`.
    !| Returns 0 and `n_finite = 0` when nothing is finite, which the caller treats as
    !| "no shift estimable" (centring by 0 is the identity).
    pure subroutine median_of_finite_helper(values, n, work, perm, stack_left, stack_right, med, n_finite)
        integer(int32), intent(in) :: n
        !! Number of candidate values
        real(real64), dimension(n), intent(in) :: values
        !! Candidate values (may contain NaN / +-Inf)
        real(real64), dimension(n), intent(inout) :: work
        !! Work array: the finite subset, compacted
        integer(int32), dimension(n), intent(inout) :: perm
        !! Work array: sort permutation
        integer(int32), dimension(n), intent(inout) :: stack_left
        !! Work array: quicksort stack
        integer(int32), dimension(n), intent(inout) :: stack_right
        !! Work array: quicksort stack
        real(real64), intent(out) :: med
        !! Median of the finite values; 0 when there are none
        integer(int32), intent(out) :: n_finite
        !! Number of finite values found

        integer(int32) :: i, half

        n_finite = 0
        do i = 1, n
            ! NaN fails the self-comparison; +-Inf fails the magnitude bound.
            if (values(i) == values(i) .and. abs(values(i)) <= huge(1.0_real64)) then
                n_finite = n_finite + 1
                work(n_finite) = values(i)
            end if
        end do

        if (n_finite == 0) then
            med = 0.0_real64
            return
        end if

        do i = 1, n_finite
            perm(i) = i
        end do
        call sort_real(work(1:n_finite), perm(1:n_finite), stack_left(1:n_finite), stack_right(1:n_finite))

        half = n_finite / 2
        if (mod(n_finite, 2) == 1) then
            med = work(perm(half + 1))
        else
            med = 0.5_real64 * (work(perm(half)) + work(perm(half + 1)))
        end if
    end subroutine median_of_finite_helper

    !> Pearson correlation of two vectors over the entries where both are finite.
    !|
    !| Returns 0 when fewer than two usable pairs exist or either side has zero
    !| variance — a degenerate correlation is reported as "no agreement measured"
    !| rather than as NaN, so the diagnostic never contaminates downstream output.
    pure function finite_correlation_helper(a, b, n) result(corr)
        integer(int32), intent(in) :: n
        !! Length of both vectors
        real(real64), dimension(n), intent(in) :: a
        !! First vector
        real(real64), dimension(n), intent(in) :: b
        !! Second vector
        real(real64) :: corr
        !! Pearson correlation over the jointly finite entries; 0 if undefined

        integer(int32) :: i, m
        real(real64) :: sa, sb, saa, sbb, sab, rm, cov, va, vb

        m = 0
        sa = 0.0_real64; sb = 0.0_real64
        saa = 0.0_real64; sbb = 0.0_real64; sab = 0.0_real64
        do i = 1, n
            if (a(i) == a(i) .and. abs(a(i)) <= huge(1.0_real64) .and. &
                b(i) == b(i) .and. abs(b(i)) <= huge(1.0_real64)) then
                m = m + 1
                sa = sa + a(i); sb = sb + b(i)
                saa = saa + a(i) * a(i); sbb = sbb + b(i) * b(i); sab = sab + a(i) * b(i)
            end if
        end do

        if (m < 2) then
            corr = 0.0_real64
            return
        end if

        rm = real(m, real64)
        cov = sab - sa * sb / rm
        va = saa - sa * sa / rm
        vb = sbb - sb * sb / rm
        if (va <= 0.0_real64 .or. vb <= 0.0_real64) then
            corr = 0.0_real64
        else
            corr = cov / sqrt(va * vb)
        end if
    end function finite_correlation_helper

    ! =========================================================================
    ! per-gene loop
    ! =========================================================================

    !> Core pipeline: per-gene multidimensional p-values from pre-built structures.
    !|
    !| For each gene, per axis and per side, gather an adaptive kNN residual pool in
    !| that axis's mean space (`gather_residuals_helper`, unchanged), optionally trim
    !| its tails, then scale it by `1/sqrt(n)` of its own side. Scaling the POOL once
    !| per gene rather than each draw is what turns the plan's per-draw
    !| `eps/sqrt(n)` into a plain pool difference, and it is exactly what
    !| `noise_model_exact` does in 1-D — which is why the `d = 1` enumerated path
    !| reproduces it.
    !|
    !| `Var_null(a)` then follows from the scaled pools with no sampling, and its
    !| reciprocal is the weight. A degenerate axis (zero pool variance on both sides)
    !| gets weight 0 rather than an infinity: it contributes nothing to `D_std`, on
    !| the observed side and the null side alike, so the statistic stays calibrated
    !| and the remaining axes decide the p-value.
    !|
    !| A gene is dropped ENTIRELY if ANY axis fails the `MIN_POOL_RESIDUALS` gate or
    !| carries a non-finite `beta`. All axes are still gathered and their sizes
    !| recorded before that decision, so the failing axis is identifiable from
    !| `neighborhood_size_case` / `neighborhood_size_control`.
    !|
    !| No allocation: every work array is owned by the alloc layer. The gene loop is
    !| sequential; under `SAMPLING_SYSTEMATIC` it has no cross-gene state at all (each
    !| gene's walk is a pure function of its index), so that path is ready to
    !| parallelise, while the RNG path draws from the single global stream and is not.
    subroutine compute_noise_pvalue_pipeline_md_helper( &
        sorted_case, sorted_control, means_case, means_control, &
        beta_obs, compute_pvalue_own, &
        n_genes, n_axes, k_start, k_step, k_max, tau, trim_frac, &
        sampling_mode, n_draws_max, n_exceed_target, enum_max_product, max_pool_size, &
        pvalues_own, d_obs, d_std_obs, d_sq_null_mean, method_used, n_draws_used, &
        neighborhood_size_case, neighborhood_size_control, var_null, &
        scale_case, scale_control, &
        tmp_pools_case, tmp_pools_control, tmp_n_pool_case, tmp_n_pool_control, &
        tmp_axis_w, rbuf, n_genes_with_pvalue, ierr)

        integer(int32), intent(in) :: n_genes
        !! Number of genes
        integer(int32), intent(in) :: n_axes
        !! Number of axes
        type(sorted_data_t), dimension(n_axes), intent(in) :: sorted_case
        !! One sorted case structure per axis
        type(sorted_data_t), dimension(n_axes), intent(in) :: sorted_control
        !! One sorted control structure per axis
        real(real64), dimension(n_genes, n_axes), intent(in) :: means_case
        !! Per-axis case means (the kNN matching coordinate)
        real(real64), dimension(n_genes, n_axes), intent(in) :: means_control
        !! Per-axis control means
        real(real64), dimension(n_axes, n_genes), intent(in) :: beta_obs
        !! Observed per-axis statistic, already centred if centring is on
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_own
        !! 1 if this gene should be tested, 0 otherwise
        integer(int32), intent(in) :: k_start
        !! Minimum neighbour genes before adaptive stopping applies
        integer(int32), intent(in) :: k_step
        !! Neighbour genes added per adaptive round
        integer(int32), intent(in) :: k_max
        !! Hard cap on neighbour genes
        real(real64), intent(in) :: tau
        !! Relative-change threshold for adaptive pool growth
        real(real64), intent(in) :: trim_frac
        !! Per-tail residual trim fraction, already norm-gated by the caller
        integer(int32), intent(in) :: sampling_mode
        !! `SAMPLING_SYSTEMATIC` or `SAMPLING_RNG`
        integer(int32), intent(in) :: n_draws_max
        !! `M`
        integer(int32), intent(in) :: n_exceed_target
        !! Besag-Clifford exceedance target; 0 disables
        integer(int64), intent(in) :: enum_max_product
        !! Enumerate while the joint space fits this
        integer(int32), intent(in) :: max_pool_size
        !! Hard cap on residuals per pool
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: d_obs
        !! Output: unweighted norm, the effect size (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: d_std_obs
        !! Output: standardised norm, the test statistic (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: d_sq_null_mean
        !! Output: mean unweighted `D^2` under the null (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: method_used
        !! Output: which path produced the p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: n_draws_used
        !! Output: draws or enumerated atoms per gene (0 if not computed)
        integer(int32), dimension(n_axes, n_genes), intent(out) :: neighborhood_size_case
        !! Output: per-axis case pool size (-1 if not gathered)
        integer(int32), dimension(n_axes, n_genes), intent(out) :: neighborhood_size_control
        !! Output: per-axis control pool size (-1 if not gathered)
        real(real64), dimension(n_axes, n_genes), intent(out) :: var_null
        !! Output: per-axis null variance; the weights, exposed (-1 if not computed)
        real(real64), dimension(n_axes), intent(in) :: scale_case
        !! `1/sqrt(n_rep_case(a))`, precomputed
        real(real64), dimension(n_axes), intent(in) :: scale_control
        !! `1/sqrt(n_rep_control(a))`, precomputed
        real(real64), dimension(:, :), intent(inout) :: tmp_pools_case
        !! Work array: case pools, `(max_pool_size * 2, n_axes)`
        real(real64), dimension(:, :), intent(inout) :: tmp_pools_control
        !! Work array: control pools, same shape
        integer(int32), dimension(n_axes), intent(inout) :: tmp_n_pool_case
        !! Work array: case pool sizes for the current gene
        integer(int32), dimension(n_axes), intent(inout) :: tmp_n_pool_control
        !! Work array: control pool sizes for the current gene
        real(real64), dimension(n_axes), intent(inout) :: tmp_axis_w
        !! Work array: the current gene's per-axis weights
        real(real64), dimension(:), intent(inout) :: rbuf
        !! Work array: pre-drawn random numbers (RNG path only)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes that received a p-value
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: i_gene, i_axis, n_c, n_t
        logical :: gene_ok
        real(real64) :: p_val, d_sq_null_val, beta_val, vn
        real(real64) :: acc_obs_std, acc_obs_plain
        integer(int32) :: n_used, meth

        call set_ok(ierr)

        pvalues_own = -1.0_real64
        d_obs = -1.0_real64
        d_std_obs = -1.0_real64
        d_sq_null_mean = -1.0_real64
        method_used = METHOD_NOT_COMPUTED
        n_draws_used = 0
        neighborhood_size_case = -1
        neighborhood_size_control = -1
        var_null = -1.0_real64
        n_genes_with_pvalue = 0

        do i_gene = 1, n_genes
            gene_ok = .true.

            ! Gather every axis before deciding, so a thin axis is reported even when
            ! it is the reason the gene is dropped.
            do i_axis = 1, n_axes
                call gather_residuals_helper(means_case(i_gene, i_axis), sorted_case(i_axis), &
                                             k_start, k_step, k_max, tau, &
                                             tmp_pools_case(:, i_axis), tmp_n_pool_case(i_axis), &
                                             max_pool_size)
                call trim_pool_tails_helper(tmp_pools_case(:, i_axis), tmp_n_pool_case(i_axis), trim_frac)

                call gather_residuals_helper(means_control(i_gene, i_axis), sorted_control(i_axis), &
                                             k_start, k_step, k_max, tau, &
                                             tmp_pools_control(:, i_axis), tmp_n_pool_control(i_axis), &
                                             max_pool_size)
                call trim_pool_tails_helper(tmp_pools_control(:, i_axis), tmp_n_pool_control(i_axis), trim_frac)

                neighborhood_size_case(i_axis, i_gene) = tmp_n_pool_case(i_axis)
                neighborhood_size_control(i_axis, i_gene) = tmp_n_pool_control(i_axis)

                if (tmp_n_pool_case(i_axis) < MIN_POOL_RESIDUALS .or. &
                    tmp_n_pool_control(i_axis) < MIN_POOL_RESIDUALS) gene_ok = .false.

                beta_val = beta_obs(i_axis, i_gene)
                ! NaN fails the self-comparison; +-Inf fails the magnitude bound.
                if (.not. (beta_val == beta_val .and. abs(beta_val) <= huge(1.0_real64))) gene_ok = .false.
            end do

            if (.not. gene_ok) cycle
            if (compute_pvalue_own(i_gene) /= 1) cycle

            ! Bring each side onto the mean-difference scale, per side, then read the
            ! null variance straight off the scaled pools.
            acc_obs_std = 0.0_real64
            acc_obs_plain = 0.0_real64
            do i_axis = 1, n_axes
                n_c = tmp_n_pool_case(i_axis)
                n_t = tmp_n_pool_control(i_axis)
                tmp_pools_case(1:n_c, i_axis) = tmp_pools_case(1:n_c, i_axis) * scale_case(i_axis)
                tmp_pools_control(1:n_t, i_axis) = tmp_pools_control(1:n_t, i_axis) * scale_control(i_axis)

                vn = pool_variance_helper(tmp_pools_case(1:n_c, i_axis), n_c) &
                   + pool_variance_helper(tmp_pools_control(1:n_t, i_axis), n_t)
                var_null(i_axis, i_gene) = vn
                if (vn > 0.0_real64) then
                    tmp_axis_w(i_axis) = 1.0_real64 / vn
                else
                    ! Degenerate axis: no noise to score against. Weight 0 removes it
                    ! from D_std on BOTH sides, rather than producing an infinity.
                    tmp_axis_w(i_axis) = 0.0_real64
                end if

                beta_val = beta_obs(i_axis, i_gene)
                acc_obs_std = acc_obs_std + tmp_axis_w(i_axis) * beta_val * beta_val
                acc_obs_plain = acc_obs_plain + beta_val * beta_val
            end do

            call compute_pvalue_md_helper(tmp_pools_case, tmp_n_pool_case, &
                                          tmp_pools_control, tmp_n_pool_control, &
                                          n_axes, tmp_axis_w, acc_obs_std, &
                                          sampling_mode, n_draws_max, n_exceed_target, &
                                          enum_max_product, rbuf, i_gene, &
                                          p_val, d_sq_null_val, n_used, meth)

            pvalues_own(i_gene) = p_val
            d_obs(i_gene) = sqrt(acc_obs_plain)
            d_std_obs(i_gene) = sqrt(acc_obs_std)
            d_sq_null_mean(i_gene) = d_sq_null_val
            method_used(i_gene) = meth
            n_draws_used(i_gene) = n_used
            n_genes_with_pvalue = n_genes_with_pvalue + 1
        end do
    end subroutine compute_noise_pvalue_pipeline_md_helper

    ! =========================================================================
    ! compute_noise_pvalue_pipeline_md
    ! =========================================================================

    !> Validate, build the `2d` sorted structures, and run the per-gene pipeline.
    !|
    !| The alloc-layer entry point. It owns every allocation and delegates the
    !| per-gene work to `compute_noise_pvalue_pipeline_md_helper`. Steps:
    !|   1. Validate dimensions, replicate counts and mode arguments.
    !|   2. Seed the RNG with `seed` (only the `SAMPLING_RNG` path consumes it).
    !|   3. Fill `beta_obs` internally (`beta_mode = 0`) or keep the caller's and
    !|      report `beta_check_cor` / `beta_check_mad` against the internal estimate.
    !|   4. Estimate `delta_hat` (per-axis median across genes) and, when
    !|      `beta_centre = 1`, subtract it from `beta_obs` IN PLACE — the caller gets
    !|      the centred vector back, which is what direction and per-axis contribution
    !|      shares are computed from.
    !|   5. Hand each axis's block of the packed buffer to `prepare_sorted_data` as an
    !|      ordinary 2-D matrix by pointer bounds remapping — no copy, and no
    !|      sequence association (these are module procedures with explicit
    !|      interfaces, where a rank mismatch is a hard error).
    !|   6. Run the gene loop.
    !|
    !| `trim_frac` is raw-normalization only, as in the scalar modules: log/voom-style
    !| transforms already stabilize the mean-variance trend, so there is nothing to
    !| trim there.
    subroutine compute_noise_pvalue_pipeline_md( &
        means_case, replicates_case_packed, n_rep_case_per_axis, &
        means_control, replicates_control_packed, n_rep_control_per_axis, &
        beta_obs, compute_pvalue_own, beta_mode, beta_centre, &
        n_genes, n_axes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
        null_method, sampling_mode, n_draws_max, n_exceed_target, enum_max_product, &
        seed, max_pool_size, &
        pvalues_own, d_obs, d_std_obs, d_sq_null_mean, &
        method_used, n_draws_used, &
        neighborhood_size_case, neighborhood_size_control, var_null, &
        beta_check_cor, beta_check_mad, delta_hat, &
        n_genes_with_pvalue, ierr)

        integer(int32), intent(in) :: n_genes
        !! Number of genes
        integer(int32), intent(in) :: n_axes
        !! Number of axes `d`
        real(real64), dimension(n_genes, n_axes), intent(in) :: means_case
        !! Per-axis case means; the kNN matching coordinate, NOT the statistic
        real(real64), dimension(:), intent(in), target, contiguous :: replicates_case_packed
        !! Case replicates, axis blocks back to back; block `a` is `(n_rep_case_a, n_genes)`,
        !! samples-fastest, matching what `prepare_sorted_data` reads
        integer(int32), dimension(n_axes), intent(in) :: n_rep_case_per_axis
        !! Case replicate count per axis (ragged by design; each must be >= 2)
        real(real64), dimension(n_genes, n_axes), intent(in) :: means_control
        !! Per-axis control means
        real(real64), dimension(:), intent(in), target, contiguous :: replicates_control_packed
        !! Control replicates, same packing
        integer(int32), dimension(n_axes), intent(in) :: n_rep_control_per_axis
        !! Control replicate count per axis (each must be >= 2)
        real(real64), dimension(n_axes, n_genes), intent(inout) :: beta_obs
        !! Observed per-axis statistic, gene-contiguous. Filled here when
        !! `beta_mode = BETA_MODE_INTERNAL`; read when `BETA_MODE_SUPPLIED`. Returned
        !! CENTRED when `beta_centre = BETA_CENTRE_MEDIAN`.
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_own
        !! 1 if the gene should be tested, 0 otherwise
        integer(int32), intent(in) :: beta_mode
        !! `BETA_MODE_INTERNAL` (0) or `BETA_MODE_SUPPLIED` (1)
        integer(int32), intent(in) :: beta_centre
        !! `BETA_CENTRE_NONE` (0) or `BETA_CENTRE_MEDIAN` (1)
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2 transform. Applied per axis.
        integer(int32), intent(in) :: k_start
        !! Minimum neighbour GENES before adaptive stopping applies
        integer(int32), intent(in) :: k_step
        !! Neighbour genes added per adaptive round
        integer(int32), intent(in) :: k_max
        !! Hard cap on neighbour genes
        real(real64), intent(in) :: tau
        !! Relative-change threshold for adaptive pool growth
        real(real64), intent(in) :: trim_frac
        !! Per-tail residual trim fraction in [0, 0.5); raw normalization only
        integer(int32), intent(in) :: null_method
        !! Kept for ABI parity with the scalar modules; only `NULL_METHOD_POOLED` (0)
        !! is accepted. Under this module's construction a draw takes ONE residual per
        !! side, and each neighbour gene contributes exactly `n_rep` residuals to the
        !! pool, so "pick a gene uniformly, then a residual within it" is the same
        !! distribution as "pick a residual uniformly from the pool". A gene-blocked
        !! argument would be a no-op that merely looked meaningful, so it is rejected
        !! rather than silently ignored — the same stance `noise_model_exact` takes.
        integer(int32), intent(in) :: sampling_mode
        !! `SAMPLING_SYSTEMATIC` (0, default) or `SAMPLING_RNG` (1)
        integer(int32), intent(in) :: n_draws_max
        !! `M`; p-value floor `1/(M + 1)` on the sampled paths. See `DEFAULT_N_DRAWS_MAX`.
        integer(int32), intent(in) :: n_exceed_target
        !! Besag-Clifford exceedance target `c`; 0 disables sequential stopping
        integer(int64), intent(in) :: enum_max_product
        !! Enumerate exactly while `prod_a S_a <= this`; 0 forces sampling.
        !! See `DEFAULT_ENUM_MAX_PRODUCT`.
        integer(int32), intent(in) :: seed
        !! RNG seed; consumed only by `SAMPLING_RNG`
        integer(int32), intent(in) :: max_pool_size
        !! Hard cap on residuals per pool
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: d_obs
        !! Output: unweighted norm `D`, the effect size in log2FC units
        real(real64), dimension(n_genes), intent(out) :: d_std_obs
        !! Output: standardised norm `D_std`, the test statistic
        real(real64), dimension(n_genes), intent(out) :: d_sq_null_mean
        !! Output: mean unweighted `D^2` under the null; `D2_adj = max(0, D^2 - this)`
        integer(int32), dimension(n_genes), intent(out) :: method_used
        !! Output: 0 enumerated, 1 systematic, 2 RNG, 3 convolution, -1 not computed
        integer(int32), dimension(n_genes), intent(out) :: n_draws_used
        !! Output: draws or enumerated atoms per gene (0 if not computed)
        integer(int32), dimension(n_axes, n_genes), intent(out) :: neighborhood_size_case
        !! Output: per-axis case pool size (-1 if not gathered)
        integer(int32), dimension(n_axes, n_genes), intent(out) :: neighborhood_size_control
        !! Output: per-axis control pool size (-1 if not gathered)
        real(real64), dimension(n_axes, n_genes), intent(out) :: var_null
        !! Output: per-axis null variance (-1 if not computed)
        real(real64), dimension(n_axes), intent(out) :: beta_check_cor
        !! Output: per-axis correlation of supplied vs internally recomputed `beta`.
        !! Exactly 1 when `beta_mode = BETA_MODE_INTERNAL` (the two are the same array).
        real(real64), dimension(n_axes), intent(out) :: beta_check_mad
        !! Output: per-axis median absolute difference, same convention (0 when internal)
        real(real64), dimension(n_axes), intent(out) :: delta_hat
        !! Output: per-axis median `beta` across genes — the estimated composition
        !! shift, reported whether or not centring is applied
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes that received a p-value
        integer(int32), intent(out) :: ierr
        !! Error code

        type(sorted_data_t), allocatable :: sorted_case(:), sorted_control(:)
        real(real64), allocatable :: beta_internal(:, :)
        real(real64), allocatable :: tmp_pools_case(:, :), tmp_pools_control(:, :)
        real(real64), allocatable :: tmp_axis_w(:), rbuf(:)
        real(real64), allocatable :: med_values(:), med_work(:)
        real(real64), allocatable :: scale_case(:), scale_control(:)
        integer(int32), allocatable :: tmp_n_pool_case(:), tmp_n_pool_control(:)
        integer(int32), allocatable :: med_perm(:), med_stack_left(:), med_stack_right(:)
        real(real64), pointer, contiguous :: axis_reps(:, :)
        integer(int32) :: i_axis, i_gene, n_rep_c, n_rep_t, off_c, off_t
        integer(int32) :: n_packed_case, n_packed_control, sort_ierr
        integer(int32) :: per_iter_max, rbuf_size, n_finite
        real(real64) :: effective_trim, med

        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_axes, ierr)
        call validate_dimension_size(k_start, ierr)
        call validate_dimension_size(k_step, ierr)
        call validate_dimension_size(k_max, ierr)
        call validate_dimension_size(n_draws_max, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        if (n_exceed_target < 0) call set_err(ierr, ERR_INVALID_INPUT)
        if (enum_max_product < 0_int64) call set_err(ierr, ERR_INVALID_INPUT)
        if (beta_mode /= BETA_MODE_INTERNAL .and. beta_mode /= BETA_MODE_SUPPLIED) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (beta_centre /= BETA_CENTRE_NONE .and. beta_centre /= BETA_CENTRE_MEDIAN) &
            call set_err(ierr, ERR_INVALID_INPUT)
        if (sampling_mode /= SAMPLING_SYSTEMATIC .and. sampling_mode /= SAMPLING_RNG) &
            call set_err(ierr, ERR_INVALID_INPUT)
        ! See the `null_method` argument doc: gene-blocking is not a distinct null
        ! under one-residual-per-side draws, so a non-zero value is rejected.
        if (null_method /= NULL_METHOD_POOLED) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        call validate_axis_replicates_helper(n_rep_case_per_axis, n_rep_control_per_axis, &
                                             n_axes, n_genes, n_packed_case, n_packed_control, ierr)
        if (is_err(ierr)) return

        if (size(replicates_case_packed) < n_packed_case .or. &
            size(replicates_control_packed) < n_packed_control) then
            call set_err(ierr, ERR_DIM_MISMATCH)
            return
        end if

        call validate_all_in_range_real(means_case, n_genes * n_axes, ierr)
        call validate_all_in_range_real(means_control, n_genes * n_axes, ierr)
        call validate_all_in_range_real(replicates_case_packed(1:n_packed_case), n_packed_case, ierr)
        call validate_all_in_range_real(replicates_control_packed(1:n_packed_control), n_packed_control, ierr)
        if (is_err(ierr)) return

        ! Only the RNG path consumes this; the systematic walk is seed-independent.
        call init_random(seed)

        M_ALLOCATE(beta_internal(n_axes, n_genes))
        M_ALLOCATE(tmp_pools_case(max_pool_size * 2, n_axes))
        M_ALLOCATE(tmp_pools_control(max_pool_size * 2, n_axes))
        M_ALLOCATE(tmp_n_pool_case(n_axes))
        M_ALLOCATE(tmp_n_pool_control(n_axes))
        M_ALLOCATE(tmp_axis_w(n_axes))
        M_ALLOCATE(scale_case(n_axes))
        M_ALLOCATE(scale_control(n_axes))
        M_ALLOCATE(med_values(n_genes))
        M_ALLOCATE(med_work(n_genes))
        M_ALLOCATE(med_perm(n_genes))
        M_ALLOCATE(med_stack_left(n_genes))
        M_ALLOCATE(med_stack_right(n_genes))
        M_ALLOCATE(sorted_case(n_axes))
        M_ALLOCATE(sorted_control(n_axes))

        ! --- observed statistic -------------------------------------------------
        call compute_beta_md_helper(replicates_case_packed, n_rep_case_per_axis, &
                                    replicates_control_packed, n_rep_control_per_axis, &
                                    n_genes, n_axes, norm_method, beta_internal)

        if (beta_mode == BETA_MODE_INTERNAL) then
            beta_obs = beta_internal
            ! Supplied and internal are the same numbers, so the diagnostic is exact
            ! by construction rather than unset.
            beta_check_cor = 1.0_real64
            beta_check_mad = 0.0_real64
        else
            do i_axis = 1, n_axes
                beta_check_cor(i_axis) = finite_correlation_helper(beta_obs(i_axis, :), &
                                                                   beta_internal(i_axis, :), n_genes)
                do i_gene = 1, n_genes
                    med_values(i_gene) = abs(beta_obs(i_axis, i_gene) - beta_internal(i_axis, i_gene))
                end do
                call median_of_finite_helper(med_values, n_genes, med_work, med_perm, &
                                             med_stack_left, med_stack_right, med, n_finite)
                beta_check_mad(i_axis) = med
            end do
        end if

        ! --- composition shift --------------------------------------------------
        do i_axis = 1, n_axes
            do i_gene = 1, n_genes
                med_values(i_gene) = beta_obs(i_axis, i_gene)
            end do
            call median_of_finite_helper(med_values, n_genes, med_work, med_perm, &
                                         med_stack_left, med_stack_right, med, n_finite)
            delta_hat(i_axis) = med
        end do

        if (beta_centre == BETA_CENTRE_MEDIAN) then
            do i_gene = 1, n_genes
                do i_axis = 1, n_axes
                    beta_obs(i_axis, i_gene) = beta_obs(i_axis, i_gene) - delta_hat(i_axis)
                end do
            end do
        end if

        ! --- sorted structures, one per (group, axis) ---------------------------
        ! Pointer bounds remapping, not a copy and not sequence association: each
        ! axis block is already (n_rep_a, n_genes) samples-fastest, which is exactly
        ! what prepare_sorted_data reads.
        off_c = 0
        off_t = 0
        do i_axis = 1, n_axes
            n_rep_c = n_rep_case_per_axis(i_axis)
            axis_reps(1:n_rep_c, 1:n_genes) => &
                replicates_case_packed(off_c + 1 : off_c + n_rep_c * n_genes)
            call prepare_sorted_data(means_case(:, i_axis), axis_reps, &
                                     n_rep_c, n_genes, norm_method, sorted_case(i_axis), sort_ierr)
            call set_err(ierr, sort_ierr)
            off_c = off_c + n_rep_c * n_genes

            n_rep_t = n_rep_control_per_axis(i_axis)
            axis_reps(1:n_rep_t, 1:n_genes) => &
                replicates_control_packed(off_t + 1 : off_t + n_rep_t * n_genes)
            call prepare_sorted_data(means_control(:, i_axis), axis_reps, &
                                     n_rep_t, n_genes, norm_method, sorted_control(i_axis), sort_ierr)
            call set_err(ierr, sort_ierr)
            off_t = off_t + n_rep_t * n_genes

            ! Per-SIDE scaling: each side divided by the root of its OWN replicate
            ! count. A single per-axis factor would be right only at equal n.
            scale_case(i_axis) = 1.0_real64 / sqrt(real(n_rep_c, real64))
            scale_control(i_axis) = 1.0_real64 / sqrt(real(n_rep_t, real64))
        end do
        if (is_err(ierr)) return

        ! --- random-number buffer (RNG path only) -------------------------------
        ! Two values per axis per draw: one index per side.
        per_iter_max = max(1, 2 * n_axes)
        if (int(n_draws_max, int64) * int(per_iter_max, int64) < int(RNG_BUFFER_MAX, int64)) then
            rbuf_size = n_draws_max * per_iter_max
        else
            rbuf_size = max(per_iter_max, (RNG_BUFFER_MAX / per_iter_max) * per_iter_max)
        end if
        M_ALLOCATE(rbuf(rbuf_size))

        effective_trim = 0.0_real64
        if (norm_method == 0) effective_trim = trim_frac

        call compute_noise_pvalue_pipeline_md_helper( &
            sorted_case, sorted_control, means_case, means_control, &
            beta_obs, compute_pvalue_own, &
            n_genes, n_axes, k_start, k_step, k_max, tau, effective_trim, &
            sampling_mode, n_draws_max, n_exceed_target, enum_max_product, max_pool_size, &
            pvalues_own, d_obs, d_std_obs, d_sq_null_mean, method_used, n_draws_used, &
            neighborhood_size_case, neighborhood_size_control, var_null, &
            scale_case, scale_control, &
            tmp_pools_case, tmp_pools_control, tmp_n_pool_case, tmp_n_pool_control, &
            tmp_axis_w, rbuf, n_genes_with_pvalue, ierr)

    end subroutine compute_noise_pvalue_pipeline_md

end module noise_model_md

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

!> C-interoperable wrapper for `compute_noise_pvalue_pipeline_md`.
!|
!| Deliberately NOT ABI-compatible with `compute_noise_pvalues_pipeline_c` /
!| `compute_noise_pvalues_pipeline_exact_c`: the multidimensional model carries an
!| extra dimension on almost every array, so it is a separate entry point rather
!| than a drop-in swap. `n_axes = 1` reduces to the scalar exact model numerically,
!| not by signature.
!|
!| Two arguments exist only at this boundary: `n_packed_case` and
!| `n_packed_control`, the lengths of the two packed replicate buffers. Fortran
!| derives them from `n_rep_*_per_axis` and `n_genes`, but C cannot be trusted to
!| have allocated what it claims, so they are passed explicitly and cross-checked
!| against the derived value — a mismatch returns `ERR_DIM_MISMATCH` instead of
!| reading past the end of the caller's buffer.
!|
!| `enum_max_product` crosses as a `c_long_long` (int64): the joint outcome space
!| passes 2^31 at `d = 2` with quite ordinary neighbourhoods, so an `int` cap would
!| be unable to express the interesting thresholds.
subroutine compute_noise_pvalues_pipeline_md_c( &
    means_case, replicates_case_packed, n_packed_case, n_rep_case_per_axis, &
    means_control, replicates_control_packed, n_packed_control, n_rep_control_per_axis, &
    beta_obs, compute_pvalue_own, beta_mode, beta_centre, &
    n_genes, n_axes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
    null_method, sampling_mode, n_draws_max, n_exceed_target, enum_max_product, &
    seed, max_pool_size, &
    pvalues_own, d_obs, d_std_obs, d_sq_null_mean, method_used, n_draws_used, &
    neighborhood_size_case, neighborhood_size_control, var_null, &
    beta_check_cor, beta_check_mad, delta_hat, &
    n_genes_with_pvalue, ierr) bind(C, name="compute_noise_pvalues_pipeline_md_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_long_long
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use noise_model_md, only: compute_noise_pvalue_pipeline_md
    use tox_errors, only: ERR_DIM_MISMATCH
    use safeguard
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
    !! Number of genes
    integer(c_int), intent(in), target :: n_axes
    !! Number of axes `d`
    integer(c_int), intent(in), target :: n_packed_case
    !! Declared length of `replicates_case_packed`; cross-checked against the derived value
    integer(c_int), intent(in), target :: n_packed_control
    !! Declared length of `replicates_control_packed`; cross-checked likewise
    real(c_double), dimension(n_genes, n_axes), intent(in), target :: means_case
    !! Per-axis case means (kNN matching coordinate)
    real(c_double), dimension(n_packed_case), intent(in), target :: replicates_case_packed
    !! Case replicates, axis blocks back to back; block `a` is `(n_rep_case_a, n_genes)`
    integer(c_int), dimension(n_axes), intent(in), target :: n_rep_case_per_axis
    !! Case replicate count per axis (each >= 2)
    real(c_double), dimension(n_genes, n_axes), intent(in), target :: means_control
    !! Per-axis control means
    real(c_double), dimension(n_packed_control), intent(in), target :: replicates_control_packed
    !! Control replicates, same packing
    integer(c_int), dimension(n_axes), intent(in), target :: n_rep_control_per_axis
    !! Control replicate count per axis (each >= 2)
    real(c_double), dimension(n_axes, n_genes), intent(inout), target :: beta_obs
    !! Observed per-axis statistic; filled here when `beta_mode = 0`, returned centred
    !! when `beta_centre = 1`
    integer(c_int), dimension(n_genes), intent(in), target :: compute_pvalue_own
    !! 1 if the gene should be tested, 0 otherwise
    integer(c_int), intent(in), target :: beta_mode
    !! 0 = compute `beta` internally (default); 1 = use the supplied `beta_obs`
    integer(c_int), intent(in), target :: beta_centre
    !! 0 = no composition centring; 1 = subtract the per-axis median
    integer(c_int), intent(in), target :: norm_method
    !! 0 = linear scale; non-zero = log2 transform
    integer(c_int), intent(in), target :: k_start
    !! Minimum neighbour genes before adaptive stopping applies
    integer(c_int), intent(in), target :: k_step
    !! Neighbour genes added per adaptive round
    integer(c_int), intent(in), target :: k_max
    !! Hard cap on neighbour genes
    real(c_double), intent(in), target :: tau
    !! Relative-change threshold for adaptive pool growth
    real(c_double), intent(in), target :: trim_frac
    !! Per-tail residual trim fraction; raw normalization only
    integer(c_int), intent(in), target :: null_method
    !! ABI parity only; must be 0 (see the Fortran argument doc)
    integer(c_int), intent(in), target :: sampling_mode
    !! 0 = systematic stride (default); 1 = RNG
    integer(c_int), intent(in), target :: n_draws_max
    !! `M`; p-value floor `1/(M + 1)` on the sampled paths
    integer(c_int), intent(in), target :: n_exceed_target
    !! Besag-Clifford exceedance target; 0 disables sequential stopping
    integer(c_long_long), intent(in), target :: enum_max_product
    !! Enumerate exactly while the joint outcome space fits this; 0 forces sampling
    integer(c_int), intent(in), target :: seed
    !! RNG seed (consumed only by the RNG sampling mode)
    integer(c_int), intent(in), target :: max_pool_size
    !! Hard cap on residuals per pool
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_own
    !! Output: p-values (-1 if not computed)
    real(c_double), dimension(n_genes), intent(out), target :: d_obs
    !! Output: unweighted norm `D`, the effect size
    real(c_double), dimension(n_genes), intent(out), target :: d_std_obs
    !! Output: standardised norm `D_std`, the test statistic
    real(c_double), dimension(n_genes), intent(out), target :: d_sq_null_mean
    !! Output: mean unweighted `D^2`; `D2_adj = max(0, D^2 - this)`
    integer(c_int), dimension(n_genes), intent(out), target :: method_used
    !! Output: 0 enumerated, 1 systematic, 2 RNG, 3 convolution, -1 not computed
    integer(c_int), dimension(n_genes), intent(out), target :: n_draws_used
    !! Output: draws or enumerated atoms per gene (0 if not computed)
    integer(c_int), dimension(n_axes, n_genes), intent(out), target :: neighborhood_size_case
    !! Output: per-axis case pool size (-1 if not gathered)
    integer(c_int), dimension(n_axes, n_genes), intent(out), target :: neighborhood_size_control
    !! Output: per-axis control pool size (-1 if not gathered)
    real(c_double), dimension(n_axes, n_genes), intent(out), target :: var_null
    !! Output: per-axis null variance; the weights, exposed (-1 if not computed)
    real(c_double), dimension(n_axes), intent(out), target :: beta_check_cor
    !! Output: per-axis correlation of supplied vs recomputed `beta` (1 when internal)
    real(c_double), dimension(n_axes), intent(out), target :: beta_check_mad
    !! Output: per-axis median absolute difference (0 when internal)
    real(c_double), dimension(n_axes), intent(out), target :: delta_hat
    !! Output: per-axis median `beta` across genes (the composition shift)
    integer(c_int), intent(out), target :: n_genes_with_pvalue
    !! Number of genes that received a p-value
    integer(c_int), intent(out), target :: ierr
    !! Error code: 0 = success

    integer(int32) :: i_axis
    integer(int64) :: total_case, total_control

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(n_packed_case)
    M_CHECK_NON_NULL(n_packed_control)
    M_CHECK_NON_NULL(means_case)
    M_CHECK_NON_NULL(replicates_case_packed)
    M_CHECK_NON_NULL(n_rep_case_per_axis)
    M_CHECK_NON_NULL(means_control)
    M_CHECK_NON_NULL(replicates_control_packed)
    M_CHECK_NON_NULL(n_rep_control_per_axis)
    M_CHECK_NON_NULL(beta_obs)
    M_CHECK_NON_NULL(compute_pvalue_own)
    M_CHECK_NON_NULL(beta_mode)
    M_CHECK_NON_NULL(beta_centre)
    M_CHECK_NON_NULL(norm_method)
    M_CHECK_NON_NULL(k_start)
    M_CHECK_NON_NULL(k_step)
    M_CHECK_NON_NULL(k_max)
    M_CHECK_NON_NULL(tau)
    M_CHECK_NON_NULL(trim_frac)
    M_CHECK_NON_NULL(null_method)
    M_CHECK_NON_NULL(sampling_mode)
    M_CHECK_NON_NULL(n_draws_max)
    M_CHECK_NON_NULL(n_exceed_target)
    M_CHECK_NON_NULL(enum_max_product)
    M_CHECK_NON_NULL(seed)
    M_CHECK_NON_NULL(max_pool_size)
    M_CHECK_NON_NULL(pvalues_own)
    M_CHECK_NON_NULL(d_obs)
    M_CHECK_NON_NULL(d_std_obs)
    M_CHECK_NON_NULL(d_sq_null_mean)
    M_CHECK_NON_NULL(method_used)
    M_CHECK_NON_NULL(n_draws_used)
    M_CHECK_NON_NULL(neighborhood_size_case)
    M_CHECK_NON_NULL(neighborhood_size_control)
    M_CHECK_NON_NULL(var_null)
    M_CHECK_NON_NULL(beta_check_cor)
    M_CHECK_NON_NULL(beta_check_mad)
    M_CHECK_NON_NULL(delta_hat)
    M_CHECK_NON_NULL(n_genes_with_pvalue)

    ! The packed-buffer lengths are the one thing the Fortran entry cannot verify:
    ! it derives them and then trusts the pointer. Check them here, where the
    ! caller's own claim is available to compare against.
    total_case = 0_int64
    total_control = 0_int64
    do i_axis = 1, n_axes
        total_case = total_case + int(n_rep_case_per_axis(i_axis), int64) * int(n_genes, int64)
        total_control = total_control + int(n_rep_control_per_axis(i_axis), int64) * int(n_genes, int64)
    end do
    if (total_case /= int(n_packed_case, int64) .or. total_control /= int(n_packed_control, int64)) then
        call set_err(ierr, ERR_DIM_MISMATCH)
        return
    end if

    call compute_noise_pvalue_pipeline_md( &
        means_case, replicates_case_packed, n_rep_case_per_axis, &
        means_control, replicates_control_packed, n_rep_control_per_axis, &
        beta_obs, compute_pvalue_own, beta_mode, beta_centre, &
        n_genes, n_axes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
        null_method, sampling_mode, n_draws_max, n_exceed_target, &
        int(enum_max_product, int64), seed, max_pool_size, &
        pvalues_own, d_obs, d_std_obs, d_sq_null_mean, method_used, n_draws_used, &
        neighborhood_size_case, neighborhood_size_control, var_null, &
        beta_check_cor, beta_check_mad, delta_hat, &
        n_genes_with_pvalue, ierr)

end subroutine compute_noise_pvalues_pipeline_md_c
