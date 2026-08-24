#include <src/macros.h>

!!! TO DO
!!! 1. Implement error code to argument number mapping
!!! 2. Find solution for family based analysis

!> Noise model for exact p-value computation comparing case and control gene expression.
!|
!| EXACT VARIANT of the `noise_model` module (see `tox_noise_model.F90`). It is
!| identical in every respect EXCEPT how the per-gene p-value is computed:
!|   - The baseline module counts the pairwise-difference null by exhaustive
!|     `O(n*m)` enumeration, falling back to Monte Carlo sampling once the pool
!|     product exceeds `MAX_EXACT_COMBINATIONS`.
!|   - This variant computes the same tail count EXACTLY at every pool size by
!|     sorting the control pool once and, for each case residual, tallying the
!|     qualifying control residuals with two binary searches — `O((n+m) log m)`,
!|     no Monte Carlo, no RNG (see `compute_pvalue_helper` below).
!| Kept as a separate module (`noise_model_exact`, C entry point
!| `compute_noise_pvalues_pipeline_exact_c`) so both can be linked and run
!| side by side to compare results and performance.
!|
!| Provides routines to:
!|   - Sort and pack gene expression residuals for efficient neighbourhood lookup
!|   - Adaptively gather residual pools from the nearest-mean neighbourhood
!|   - Compute exact p-values by counting, via sorted-pool binary search, every
!|     pairwise residual difference that meets or exceeds the observed statistic
!|     (fully deterministic, no Monte Carlo sampling, exact at all pool sizes)
module noise_model_exact
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_NAN_INF, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int
    use f42_utils, only: sort_real, sort_integer
    implicit none

    !> Sorted gene data: means and packed centred residuals in mean-ascending order.
    type :: sorted_data_t
        real(real64), allocatable :: means_sorted(:)
        !! Sorted gene means (ascending)
        integer(int32), allocatable :: original_indices(:)
        !! Original gene index at each sorted position
        integer(int32), allocatable :: n_residuals(:)
        !! Number of residuals stored for each sorted gene
        real(real64), allocatable :: residuals_packed(:, :)
        !! Centred residuals: shape (max_resid_per_gene, n_genes)
        integer(int32) :: max_resid_per_gene
        !! Maximum residuals per gene (= n_samples)
        integer(int32) :: n_genes
        !! Total number of genes in the sorted structure
    end type sorted_data_t

    ! Note: this exact variant needs neither MAX_EXACT_COMBINATIONS nor
    ! MONTE_CARLO_SAMPLES — the sorted-pool tail count is exact at every pool size,
    ! so there is no enumeration/Monte-Carlo threshold and no RNG.

    real(real64), parameter :: NOISE_LOG_OFFSET = 1.0_real64
        !! Additive constant `c` in the log-space residual
        !! `log2(r + c) - mean_i[log2(r_i + c)]` (see `prepare_sorted_data_helper`),
        !! used when `norm_method /= 0` to keep the logarithm defined at zero expression

contains

    ! =========================================================================
    ! prepare_sorted_data
    ! =========================================================================

    !> Core implementation: sort genes by mean and pack centred residuals.
    !|
    !| When `norm_method == 0`, residuals are the plain linear-space deviation
    !| `r_{i,g} - mu_g`, centred on the arithmetic mean `mu_g` of the replicates.
    !|
    !| When `norm_method /= 0`, residuals are computed directly in log2-expression
    !| space, centred on the Frechet mean (barycenter) of that space rather than
    !| on `log2(mu_g + c)`:
    !|
    !|   ghat_g   = (1/n) * sum_i log2(r_{i,g} + c)   (= log2 of the geometric mean
    !|                                                    of the pseudo-count-shifted
    !|                                                    replicates, c = NOISE_LOG_OFFSET)
    !|   epsilon_{i,g} = log2(r_{i,g} + c) - ghat_g
    !|
    !| Using `log2(mu_g + c)` (the log of the linear-space mean) instead of `ghat_g`
    !| would NOT generally leave the residuals centred at zero, since log is concave
    !| (Jensen's inequality: `log2(mu_g + c) >= ghat_g`, with equality only when all
    !| replicates are equal). `ghat_g` is the arithmetic mean of the log-transformed
    !| replicates themselves, so `sum_i epsilon_{i,g} = 0` exactly, by construction.
    !|
    !| This is the single place the log transform is applied; `compute_pvalue_helper`
    !| simply differences whatever scale the pooled residuals already carry.
    !|
    !| Fills `sorted_data` in-place. Does not allocate — all arrays inside
    !| `sorted_data` must already be allocated by the caller to the correct sizes,
    !| and `tmp_perm`, `tmp_stack_left`, `tmp_stack_right` must be pre-allocated
    !| work arrays of length `n_genes`.
    pure subroutine prepare_sorted_data_helper(means, replicates, n_samples, n_genes, &
                                               norm_method, &
                                               sorted_data, &
                                               tmp_perm, tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: n_samples
        !! Number of replicates per gene
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        real(real64), dimension(n_genes), intent(in) :: means
        !! Per-gene expression means (length n_genes)
        real(real64), dimension(n_samples, n_genes), intent(in) :: replicates
        !! Replicate expression matrix (n_samples x n_genes)
        integer(int32), intent(in) :: norm_method
        !! 0 = linear-scale residuals; non-zero = log2-space residuals
        type(sorted_data_t), intent(inout) :: sorted_data
        !! Sorted data structure to fill; all allocatable fields must be pre-allocated
        integer(int32), dimension(n_genes), intent(inout) :: tmp_perm
        !! Work array: permutation vector for indirect sort (length n_genes)
        integer(int32), dimension(n_genes), intent(inout) :: tmp_stack_left
        !! Work array: quicksort left-index stack (length n_genes)
        integer(int32), dimension(n_genes), intent(inout) :: tmp_stack_right
        !! Work array: quicksort right-index stack (length n_genes)

        integer(int32) :: i_gene, i_sample, orig_idx
        real(real64) :: gene_mean, log2_mean, log2_factor, bessel
        logical :: use_log_transform

        sorted_data%n_genes = n_genes
        sorted_data%max_resid_per_gene = n_samples

        use_log_transform = (norm_method /= 0)
        if (use_log_transform) log2_factor = 1.0_real64 / log(2.0_real64)

        ! Bessel correction. Residuals x - xbar have variance sigma^2 (n-1)/n, so the
        ! raw pool understates sigma (18% at n=3). Scale every residual by sqrt(n/(n-1))
        ! once here, at construction, so every downstream null is built from
        ! unbiased-variance residuals while the observed mean-difference statistic is
        ! left untouched. prepare_sorted_data enforces n_samples >= 2.
        bessel = sqrt(real(n_samples, real64) / real(max(n_samples - 1, 1), real64))

        do concurrent(i_gene=1:n_genes) shared(tmp_perm)
            tmp_perm(i_gene) = i_gene
        end do

        call sort_real(means, tmp_perm, tmp_stack_left, tmp_stack_right)

        do i_gene = 1, n_genes
            orig_idx = tmp_perm(i_gene)
            sorted_data%original_indices(i_gene) = orig_idx
            sorted_data%means_sorted(i_gene) = means(orig_idx)
            sorted_data%n_residuals(i_gene) = n_samples

            gene_mean = sum(replicates(:, orig_idx)) / real(n_samples, real64)
            if (use_log_transform) then
                ! Frechet mean in log2-space: the arithmetic mean of the per-replicate
                ! log2-transformed values, NOT log2(gene_mean + c) — see subroutine doc.
                log2_mean = sum(log(max(replicates(:, orig_idx), 0.0_real64) + NOISE_LOG_OFFSET)) &
                            * log2_factor / real(n_samples, real64)
            end if

            do concurrent(i_sample=1:n_samples) &
                shared(sorted_data, replicates, gene_mean, log2_mean, log2_factor, &
                       use_log_transform, i_gene, orig_idx, bessel)
                if (use_log_transform) then
                    sorted_data%residuals_packed(i_sample, i_gene) = bessel * ( &
                        log(max(replicates(i_sample, orig_idx), 0.0_real64) + NOISE_LOG_OFFSET) &
                        * log2_factor - log2_mean)
                else
                    sorted_data%residuals_packed(i_sample, i_gene) = bessel * ( &
                        replicates(i_sample, orig_idx) - gene_mean)
                end if
            end do
        end do
    end subroutine prepare_sorted_data_helper

    !> Validate inputs, allocate the `sorted_data` structure, and sort genes by mean.
    !|
    !| This is the alloc-layer entry point. It allocates all fields of `sorted_data`
    !| and the internal work arrays, then delegates to `prepare_sorted_data_helper`.
    subroutine prepare_sorted_data(means, replicates, n_samples, n_genes, norm_method, sorted_data, ierr)
        integer(int32), intent(in) :: n_samples
        !! Number of replicates per gene
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        real(real64), dimension(n_genes), intent(in) :: means
        !! Per-gene expression means (length n_genes)
        real(real64), dimension(n_samples, n_genes), intent(in) :: replicates
        !! Replicate expression matrix (n_samples x n_genes)
        integer(int32), intent(in) :: norm_method
        !! 0 = linear-scale residuals; non-zero = log2-space residuals (see `prepare_sorted_data_helper`)
        type(sorted_data_t), intent(out) :: sorted_data
        !! Sorted data structure; all fields are allocated here
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: stack_size, alloc_stat
        integer(int32), allocatable :: tmp_perm(:), tmp_stack_left(:), tmp_stack_right(:)

        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_samples, ierr)
        ! Need >= 2 replicates to form a residual variance (and the Bessel correction).
        if (n_samples < 2) call set_err(ierr, ERR_INVALID_INPUT)
        call validate_all_in_range_real(means, n_genes, ierr)
        call validate_all_in_range_real(replicates, n_samples * n_genes, ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(sorted_data%means_sorted(n_genes))
        M_ALLOCATE(sorted_data%original_indices(n_genes))
        M_ALLOCATE(sorted_data%n_residuals(n_genes))
        M_ALLOCATE(sorted_data%residuals_packed(n_samples, n_genes))

        stack_size = 2 * int(log(real(n_genes, real64)) / log(2.0_real64)) + 10
        M_ALLOCATE(tmp_perm(n_genes))
        M_ALLOCATE(tmp_stack_left(stack_size))
        M_ALLOCATE(tmp_stack_right(stack_size))

        call prepare_sorted_data_helper(means, replicates, n_samples, n_genes, norm_method, &
                                        sorted_data, tmp_perm, tmp_stack_left, tmp_stack_right)
    end subroutine prepare_sorted_data

    ! =========================================================================
    ! find_closest_helper
    ! =========================================================================

    !> Binary search: return the index in `means_sorted` closest to `target`.
    !|
    !| Returns 0 when `means_sorted` is empty. Caller must guard against this.
    pure function find_closest_helper(target, means_sorted, n) result(pos)
        real(real64), intent(in) :: target
        !! Query value
        integer(int32), intent(in) :: n
        !! Length of `means_sorted`
        real(real64), dimension(n), intent(in) :: means_sorted
        !! Sorted mean array (ascending)
        integer(int32) :: pos
        !! Index of the element closest to `target`; 0 if `n == 0`

        integer(int32) :: left, right, mid

        if (n == 0) then
            pos = 0
            return
        end if

        if (target <= means_sorted(1)) then
            pos = 1
            return
        end if

        if (target >= means_sorted(n)) then
            pos = n
            return
        end if

        left = 1
        right = n
        do while (left <= right)
            mid = (left + right) / 2
            if (means_sorted(mid) < target) then
                left = mid + 1
            else if (means_sorted(mid) > target) then
                right = mid - 1
            else
                pos = mid
                return
            end if
        end do

        if (left == 1) then
            pos = 1
        else if (left > n) then
            pos = n
        else
            if (abs(means_sorted(left) - target) < abs(means_sorted(left - 1) - target)) then
                pos = left
            else
                pos = left - 1
            end if
        end if
    end function find_closest_helper

    ! =========================================================================
    ! add_residuals_to_pool_helper
    ! =========================================================================

    !> Append residuals from one gene into a pool, respecting the capacity `new_size`.
    !|
    !| Copies at most `new_size - current_size` values from `residuals` into
    !| `pool(current_size+1 : new_size)`. No allocation; caller owns all arrays.
    pure subroutine add_residuals_to_pool_helper(pool, current_size, residuals, n_resid, new_size)
        real(real64), intent(inout) :: pool(:)
        !! Target residual pool (pre-allocated by caller)
        integer(int32), intent(in) :: current_size
        !! Number of elements already in `pool`
        real(real64), intent(in) :: residuals(:)
        !! Source residuals to append
        integer(int32), intent(in) :: n_resid
        !! Number of elements in `residuals`
        integer(int32), intent(in) :: new_size
        !! Maximum number of elements allowed in `pool` after the append

        integer(int32) :: n_to_copy

        n_to_copy = min(n_resid, new_size - current_size)
        if (n_to_copy > 0) then
            pool(current_size + 1:current_size + n_to_copy) = residuals(1:n_to_copy)
        end if
    end subroutine add_residuals_to_pool_helper

    ! =========================================================================
    ! gather_residuals
    ! =========================================================================

    !> Core implementation: adaptively gather residuals from the nearest-mean neighbourhood.
    !|
    !| Starting from the gene whose mean is closest to `target_mean`, neighbour GENES are
    !| added outward (alternating left/right in mean-sorted order) until one of:
    !|   - The neighbourhood contains at least `k_start` genes (initial phase), then
    !|   - The relative change in mean absolute residual exceeds `tau` (adaptive phase), or
    !|   - The neighbourhood reaches `k_max` genes.
    !|
    !| `k_start`, `k_step` and `k_max` are counted in GENES, not residuals, so the collected
    !| expression window is independent of how many replicates each gene has (each gene
    !| contributes `sorted_data%n_residuals` = n_replicates residuals). The resulting
    !| residual pool holds up to `k_max * n_replicates` values, hard-capped at `max_pool_size`.
    !|
    !| No allocation; `pooled_residuals` must be
    !| pre-allocated by the caller to at least `max_pool_size`. Candidate expansions
    !| for a round are staged in place within `pooled_residuals` past the committed
    !| size, so no separate staging buffer is needed.
    pure subroutine gather_residuals_helper(target_mean, sorted_data, &
                                            k_start, k_step, k_max, tau, &
                                            pooled_residuals, n_pooled, &
                                            max_pool_size)
        real(real64), intent(in) :: target_mean
        !! Mean value for which a matching residual neighbourhood is sought
        type(sorted_data_t), intent(in) :: sorted_data
        !! Pre-built sorted gene structure
        integer(int32), intent(in) :: k_start
        !! Minimum number of neighbour GENES before the adaptive stopping criterion applies
        integer(int32), intent(in) :: k_step
        !! Number of new neighbour GENES to add per adaptive round before re-evaluating
        integer(int32), intent(in) :: k_max
        !! Hard upper limit on the number of neighbour GENES (residuals still capped at max_pool_size)
        real(real64), intent(in) :: tau
        !! Relative-change threshold; expansion stops when the change exceeds this value
        real(real64), intent(out) :: pooled_residuals(:)
        !! Output residual pool (pre-allocated to at least `max_pool_size`)
        integer(int32), intent(out) :: n_pooled
        !! Number of residuals written into `pooled_residuals`
        integer(int32), intent(in) :: max_pool_size
        !! Allocated size of `pooled_residuals` (hard cap on the residual count)

        integer(int32) :: pos, left_cand, right_cand, idx, i_new
        integer(int32) :: current_size, genes_added, offset, trial_size, n_genes_pool
        integer(int32) :: pool_size, n_resid
        real(real64) :: S_old, S_new, rel_change, abs_sum, trial_abs_sum

        n_pooled = 0
        current_size = 0
        n_genes_pool = 0

        pos = find_closest_helper(target_mean, sorted_data%means_sorted, sorted_data%n_genes)
        if (pos == 0) return

        current_size = min(sorted_data%n_residuals(pos), max_pool_size)
        pooled_residuals(1:current_size) = sorted_data%residuals_packed(1:current_size, pos)
        n_genes_pool = 1

        left_cand = pos - 1
        right_cand = pos + 1

        ! Phase 1: expand until the neighbourhood reaches k_start GENES, regardless of tau.
        ! (k_start/k_step/k_max are gene counts; each gene adds n_replicates residuals, and
        !  the residual pool is separately hard-capped at max_pool_size for safety.)
        do while (n_genes_pool < k_start .and. current_size < max_pool_size .and. &
                  (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))

            ! choose closest candidate gene
            call choose_index(sorted_data%means_sorted, sorted_data%n_genes, target_mean, idx, left_cand, right_cand)

            ! Add its residuals to the pool
            pool_size = min(current_size + sorted_data%n_residuals(idx), max_pool_size)
            call add_residuals_to_pool_helper(pooled_residuals, current_size, &
                                              sorted_data%residuals_packed(:, idx), &
                                              sorted_data%n_residuals(idx), pool_size)
            current_size = pool_size
            n_genes_pool = n_genes_pool + 1
        end do

        ! If fewer than 10 residuals were collected, return (too little to model)
        if (current_size < 10) return

        ! Running sum of |residual| over the committed pool, maintained incrementally
        ! so each Phase 2 round costs O(residuals added) rather than O(whole pool).
        abs_sum = sum(abs(pooled_residuals(1:current_size)))
        S_old = abs_sum / real(current_size, real64)
        if (S_old == 0.0_real64) then
            n_pooled = current_size
            return
        end if

        ! Phase 2: adaptive expansion by GENES. A round's candidate residuals are appended
        ! in place past `current_size`; a rejected round is discarded by not advancing
        ! `current_size` / `n_genes_pool` (trailing trial data beyond `n_pooled` is never read).
        do while (n_genes_pool < k_max .and. current_size < max_pool_size .and. &
                  (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))

            genes_added = 0
            offset = 0
            trial_abs_sum = abs_sum

            do while (genes_added < k_step .and. &
                      n_genes_pool + genes_added < k_max .and. &
                      current_size + offset < max_pool_size .and. &
                      (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))

                call choose_index(sorted_data%means_sorted, sorted_data%n_genes, target_mean, idx, left_cand, right_cand)

                n_resid = sorted_data%n_residuals(idx)
                pool_size = min(current_size + offset + n_resid, max_pool_size)
                call add_residuals_to_pool_helper(pooled_residuals, current_size + offset, &
                                                  sorted_data%residuals_packed(:, idx), &
                                                  n_resid, pool_size)
                ! Accumulate |residual| over only the residuals actually written this step.
                do i_new = current_size + offset + 1, pool_size
                    trial_abs_sum = trial_abs_sum + abs(pooled_residuals(i_new))
                end do
                offset = offset + n_resid
                genes_added = genes_added + 1
            end do

            if (genes_added == 0) exit

            trial_size = min(current_size + offset, max_pool_size)
            S_new = trial_abs_sum / real(trial_size, real64)
            rel_change = (S_new - S_old) / S_old

            if (rel_change > tau) exit

            current_size = trial_size
            n_genes_pool = n_genes_pool + genes_added
            abs_sum = trial_abs_sum
            S_old = S_new
        end do

        n_pooled = current_size
    end subroutine gather_residuals_helper

    !> Symmetric quantile trim of a gathered residual pool (raw normalization only).
    !|
    !| Sorts the first `n_pool` residuals ascending and drops the `k` smallest and
    !| `k` largest, where `k = floor(n_pool * trim_frac)`, keeping the central
    !| `n_pool - 2k` residuals compacted into `pool(1:n_pool)`. Purpose: in RAW
    !| (linear) space the mean-neighbourhood still carries genuine variance
    !| heterogeneity and heavy tails, so a few extreme residuals can inflate the
    !| null; trimming both tails removes such artificial outliers. When the pool is
    !| already homogeneous, the trimmed residuals sit close to the rest, so little
    !| is lost. Under log normalization the trend is stabilized, so trimming is not
    !| applied there (the caller passes `trim_frac = 0`).
    !|
    !| No-ops (pool unchanged) when `trim_frac <= 0`, when `k` rounds down to 0
    !| (pool too small to trim — e.g. n_pool < 1/trim_frac), or when trimming would
    !| empty the pool (`n_pool - 2k < 1`, e.g. trim_frac >= 0.5). The downstream
    !| `< 10 residuals` gate still applies to the trimmed size. Kept byte-for-byte
    !| identical to the copy in `tox_noise_model.F90` (shared machinery).
    pure subroutine trim_pool_tails_helper(pool, n_pool, trim_frac)
        real(real64), intent(inout) :: pool(:)
        !! Residual pool; on return its central residuals occupy pool(1:n_pool)
        integer(int32), intent(inout) :: n_pool
        !! Number of valid residuals in `pool`; reduced to the kept count on return
        real(real64), intent(in) :: trim_frac
        !! Fraction to trim from EACH tail (e.g. 0.05 keeps the central 90%)

        integer(int32) :: k, n_keep, i
        integer(int32) :: perm(n_pool), stack_left(n_pool), stack_right(n_pool)
        real(real64) :: sorted_pool(n_pool)

        if (trim_frac <= 0.0_real64 .or. n_pool < 1) return
        k = int(real(n_pool, real64) * trim_frac, int32)   ! residuals dropped per tail
        if (k <= 0) return                                 ! pool too small to trim
        n_keep = n_pool - 2 * k
        if (n_keep < 1) return                             ! never trim the pool empty

        do concurrent(i = 1:n_pool) shared(perm)
            perm(i) = i
        end do
        call sort_real(pool(1:n_pool), perm, stack_left, stack_right)
        do concurrent(i = 1:n_pool) shared(sorted_pool, pool, perm)
            sorted_pool(i) = pool(perm(i))
        end do
        ! Keep the central n_keep residuals (drop the k smallest and k largest).
        do concurrent(i = 1:n_keep) shared(pool, sorted_pool, k)
            pool(i) = sorted_pool(k + i)
        end do
        n_pool = n_keep
    end subroutine trim_pool_tails_helper

    pure subroutine choose_index(means_sorted, n_genes, target_mean, idx, left_cand, right_cand)
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        real(real64), intent(in), dimension(n_genes) :: means_sorted
        !! Sorted means per gene
        real(real64), intent(in) :: target_mean
        !! Current mean
        integer(int32), intent(out) :: idx
        !! Index of the next mean
        integer(int32), intent(inout) :: left_cand
        !! Left candidate index
        integer(int32), intent(inout) :: right_cand
        !! Right candidate index

        if (left_cand >= 1 .and. right_cand <= n_genes) then
            if (abs(means_sorted(left_cand) - target_mean) <= &
                abs(means_sorted(right_cand) - target_mean)) then
                idx = left_cand
                left_cand = left_cand - 1
            else
                idx = right_cand
                right_cand = right_cand + 1
            end if
        else if (left_cand >= 1) then
            idx = left_cand
            left_cand = left_cand - 1
        else
            idx = right_cand
            right_cand = right_cand + 1
        end if
    end subroutine

    ! =========================================================================
    ! compute_pvalue
    ! =========================================================================

    !> Count entries of an ascending array below a threshold, via binary search.
    !|
    !| Returns the number of entries in `arr(1:n)` (assumed sorted ascending) that
    !| are `< x` when `inclusive` is `.false.`, or `<= x` when `inclusive` is `.true.`.
    !| O(log n). This is the primitive used by `compute_pvalue_helper` to tally the
    !| tail of the pairwise-difference null without enumerating the pairs.
    pure function count_below_helper(arr, n, x, inclusive) result(cnt)
        integer(int32), intent(in) :: n
        !! Number of entries in `arr`
        real(real64), dimension(n), intent(in) :: arr
        !! Ascending array to search
        real(real64), intent(in) :: x
        !! Threshold value
        logical, intent(in) :: inclusive
        !! `.true.` counts `arr <= x`; `.false.` counts `arr < x`
        integer(int32) :: cnt
        !! Number of qualifying entries

        integer(int32) :: lo, hi, mid

        lo = 1
        hi = n
        cnt = 0
        do while (lo <= hi)
            mid = (lo + hi) / 2
            if ((inclusive .and. arr(mid) <= x) .or. ((.not. inclusive) .and. arr(mid) < x)) then
                ! arr is ascending, so arr(1:mid) all satisfy the condition; keep
                ! this count and look right for more.
                cnt = mid
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function count_below_helper

    !> Core implementation: exact p-value via sorted-pool tail counting.
    !|
    !| Computes exactly the same quantity as exhaustive pairwise enumeration — the
    !| number of pairwise absolute differences `|r_case - r_control|` (over every
    !| case/control residual pair) that meet or exceed `observed_statistic_abs` —
    !| but WITHOUT materialising the `n_pool_case * n_pool_control` differences.
    !| Means are intentionally excluded: we are testing how much deviation is
    !| explained by noise alone.
    !|
    !| The control pool is sorted once. Then, for each case residual `a`, the count
    !| of control residuals `b` with `|a - b| >= t` (`t = observed_statistic_abs`) is
    !|
    !|   m - #{ b in the open window (a - t, a + t) }
    !|     = m - ( #{ b < a + t } - #{ b <= a - t } )
    !|
    !| via two binary searches (`count_below_helper`). Summing over the case pool
    !| gives the exact tail count in `O((n + m) log m)` instead of `O(n*m)`, so this
    !| variant is exact at every pool size and needs no Monte Carlo fallback.
    !|
    !| The `max(0, ...)` guard covers the degenerate `t = 0` case: the window
    !| `(a, a)` is empty, so every pair qualifies and `p = 1`, exactly as
    !| enumeration would give. The strict-vs-inclusive split of the two searches
    !| reproduces the `>= t` (tie-inclusive) semantics of the enumeration path.
    !|
    !| The pooled residuals are assumed to already be on whatever scale the observed
    !| statistic uses (linear or log2) — see `prepare_sorted_data_helper`, which is
    !| the single place the log transform is applied.
    pure subroutine compute_pvalue_helper(pool_case, n_pool_case, &
                                          pool_control, n_pool_control, &
                                          observed_statistic_abs, &
                                          p_value)
        integer(int32), intent(in) :: n_pool_case
        !! Size of the case residual pool
        real(real64), dimension(n_pool_case), intent(in) :: pool_case
        !! Case residual pool
        integer(int32), intent(in) :: n_pool_control
        !! Size of the control residual pool
        real(real64), dimension(n_pool_control), intent(in) :: pool_control
        !! Control residual pool
        real(real64), intent(in) :: observed_statistic_abs
        !! Absolute value of the observed test statistic
        real(real64), intent(out) :: p_value
        !! Exact p-value

        integer(int32) :: i, i_case, hi, lo
        integer(int32) :: perm(n_pool_control), stack_left(n_pool_control), stack_right(n_pool_control)
        real(real64) :: a, sorted_control(n_pool_control)
        integer(int64) :: count_ge, n_pairs

        ! Sort the control pool ascending (indirect sort, then materialise the order
        ! into a contiguous array for cache-friendly binary search).
        do concurrent(i=1:n_pool_control) shared(perm)
            perm(i) = i
        end do
        call sort_real(pool_control, perm, stack_left, stack_right)
        do concurrent(i=1:n_pool_control) shared(sorted_control, pool_control, perm)
            sorted_control(i) = pool_control(perm(i))
        end do

        ! Tail count = sum over case residuals of #{ b : |a - b| >= t }.
        count_ge = 0_int64
        do concurrent(i_case=1:n_pool_case) &
            shared(pool_case, sorted_control, observed_statistic_abs, n_pool_control) &
            reduce(+:count_ge) &
            local(a, hi, lo)

            a  = pool_case(i_case)
            hi = count_below_helper(sorted_control, n_pool_control, a + observed_statistic_abs, .false.)  ! #{b <  a+t}
            lo = count_below_helper(sorted_control, n_pool_control, a - observed_statistic_abs, .true.)   ! #{b <= a-t}
            count_ge = count_ge + int(n_pool_control - max(0_int32, hi - lo), int64)
        end do

        n_pairs = int(n_pool_case, int64) * int(n_pool_control, int64)
        p_value = real(count_ge + 1_int64, real64) / real(n_pairs + 1_int64, real64)
    end subroutine compute_pvalue_helper

    !> Validate inputs and compute the exact p-value.
    !|
    !| This is the validated entry point: it checks dimensions and residual-pool
    !| values, then delegates to `compute_pvalue_helper`, which computes the exact
    !| tail count at any pool size (no Monte Carlo fallback in this variant).
    pure subroutine compute_pvalue(pool_case, n_pool_case, &
                                   pool_control, n_pool_control, &
                                   observed_statistic, &
                                   p_value, ierr)
        integer(int32), intent(in) :: n_pool_case
        !! Size of the case residual pool
        real(real64), dimension(n_pool_case), intent(in) :: pool_case
        !! Case residual pool
        integer(int32), intent(in) :: n_pool_control
        !! Size of the control residual pool
        real(real64), dimension(n_pool_control), intent(in) :: pool_control
        !! Control residual pool
        real(real64), intent(in) :: observed_statistic
        !! Observed test statistic (sign is discarded; absolute value is used)
        real(real64), intent(out) :: p_value
        !! Exact p-value
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_pool_case, ierr)
        call validate_dimension_size(n_pool_control, ierr)
        call validate_all_in_range_real(pool_case, n_pool_case, ierr)
        call validate_all_in_range_real(pool_control, n_pool_control, ierr)
        if (is_err(ierr)) return

        call compute_pvalue_helper(pool_case, n_pool_case, &
                                   pool_control, n_pool_control, &
                                   abs(observed_statistic), &
                                   p_value)
    end subroutine compute_pvalue

    ! =========================================================================
    ! compute_noise_pvalue_pipeline
    ! =========================================================================

    !> Core pipeline: compute per-gene noise p-values using pre-built data structures.
    !|
    !| Iterates over all genes and computes, for each, `pvalues_own`: the gene vs. its
    !| own matched neighbourhood. For both the case and control sides an adaptive kNN
    !| residual pool is gathered in mean-expression space (`gather_residuals_helper`);
    !| each pool is then scaled by `1/sqrt(n_replicates)` to bring an individual-residual
    !| null onto the mean-difference scale, and the `own` p-value is computed EXACTLY
    !| from those scaled pools by `compute_pvalue`.
    !|
    !| Requires `sorted_case` and `sorted_control` to already be built (via
    !| `prepare_sorted_data`).
    !|
    !| No allocation: the two per-gene residual pools are pre-allocated by the caller
    !| and passed in as `intent(inout)` work arrays. This exact variant needs no RNG
    !| and no Monte Carlo buffer — the p-value is computed exactly at every pool size
    !| by `compute_pvalue`.
    subroutine compute_noise_pvalue_pipeline_helper( &
        sorted_case, sorted_control, &
        means_case, means_control, &
        observed_statistic_own, &
        compute_pvalue_own, &
        n_genes, k_start, k_step, k_max, tau, trim_frac, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        tmp_pool_case, tmp_pool_control_own, &
        ierr)

        type(sorted_data_t), intent(in) :: sorted_case
        !! Sorted case gene data
        type(sorted_data_t), intent(in) :: sorted_control
        !! Sorted control gene data
        integer(int32), intent(in) :: n_genes
        !! Total number of genes
        real(real64), dimension(n_genes), intent(in) :: means_case
        !! Per-gene case expression means
        real(real64), dimension(n_genes), intent(in) :: means_control
        !! Per-gene control expression means
        real(real64), dimension(n_genes), intent(in) :: observed_statistic_own
        !! Observed gene-vs-own statistic for each gene
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_own
        !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
        integer(int32), intent(in) :: k_start
        !! Minimum pool size before adaptive stopping is applied
        integer(int32), intent(in) :: k_step
        !! Residuals added per adaptive round before re-evaluating
        integer(int32), intent(in) :: k_max
        !! Hard upper limit on residual pool size
        real(real64), intent(in) :: tau
        !! Relative-change threshold for adaptive pool growth
        real(real64), intent(in) :: trim_frac
        !! Symmetric per-tail residual-pool trim fraction (0 = no trimming); the
        !! caller passes the norm-gated value (raw only)
        integer(int32), intent(in) :: max_pool_size
        !! Allocated size of all pool arrays
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: gene-vs-own p-values (-1 if not computed)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes for which the own p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Case kNN-pool size used for the gene-vs-own comparison (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Control kNN-pool size used for the gene-vs-own comparison (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for this gene (-1 if not computed)
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool_case
        !! Work array: output residual pool for this gene's case kNN neighbourhood
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool_control_own
        !! Work array: output residual pool for this gene's control kNN neighbourhood
        !! (the `own` comparison)
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_gene
        real(real64) :: mean_case_val, mean_control_val
        real(real64) :: observed_statistic_own_val
        integer(int32) :: n_pool_case, n_pool_control_own
        real(real64) :: own_scale_case, own_scale_control

        call set_ok(ierr)

        ! Quick variance correction for the individual-vs-mean scale mismatch of the
        ! `own` null: the observed statistic is a difference of MEANS (over
        ! n_replicates each side), but the null is built from differences of
        ! INDIVIDUAL residuals. Scaling each side's residuals by 1/sqrt(its own
        ! replicate count) makes a pairwise-difference of the scaled pools carry
        ! variance sigma_case^2/n_case + sigma_control^2/n_control — matching the
        ! mean difference, per side, with NO equal-variance assumption. This fixes
        ! the ~sqrt(n) over-width only; it does not reshape the null via the CLT.
        own_scale_case    = 1.0_real64 / sqrt(real(sorted_case%max_resid_per_gene, real64))
        own_scale_control = 1.0_real64 / sqrt(real(sorted_control%max_resid_per_gene, real64))

        pvalues_own = -1.0_real64
        neighborhood_size_own_case = -1
        neighborhood_size_own_control = -1
        neighborhood_size_case = -1
        n_genes_with_pvalue = 0

        do i_gene = 1, n_genes
            mean_case_val = means_case(i_gene)
            mean_control_val = means_control(i_gene)

            call gather_residuals_helper(mean_case_val, sorted_case, &
                                         k_start, k_step, k_max, tau, &
                                         tmp_pool_case, n_pool_case, &
                                         max_pool_size)
            ! Raw-only outlier trim of each tail (no-op when trim_frac == 0). Applied
            ! before both the < 10 gate and the sqrt-scaling below, so the exact
            ! null is built from the trimmed central residuals.
            call trim_pool_tails_helper(tmp_pool_case, n_pool_case, trim_frac)
            if (n_pool_case < 10) cycle

            call gather_residuals_helper(mean_control_val, sorted_control, &
                                         k_start, k_step, k_max, tau, &
                                         tmp_pool_control_own, n_pool_control_own, &
                                         max_pool_size)
            call trim_pool_tails_helper(tmp_pool_control_own, n_pool_control_own, trim_frac)
            if (n_pool_control_own < 10) cycle

            observed_statistic_own_val = observed_statistic_own(i_gene)

            ! Skip genes with a non-finite observed statistic
            if (observed_statistic_own_val /= observed_statistic_own_val) cycle

            ! The exact p-value is computed at every pool size, so there is no
            ! Monte Carlo path here and nothing to pre-seed per gene.

            if (compute_pvalue_own(i_gene) == 1) then
                ! Scale each gathered kNN pool to mean-difference variance (see the
                ! own_scale_* definitions above), then compute the exact `own` p-value
                ! directly from the scaled pools. Observed statistic is left as-is; both
                ! pools are already gated >= 10 above.
                tmp_pool_case(1:n_pool_case) = tmp_pool_case(1:n_pool_case) * own_scale_case
                tmp_pool_control_own(1:n_pool_control_own) = &
                    tmp_pool_control_own(1:n_pool_control_own) * own_scale_control
                call compute_pvalue(tmp_pool_case(1:n_pool_case), n_pool_case, &
                                    tmp_pool_control_own(1:n_pool_control_own), n_pool_control_own, &
                                    observed_statistic_own_val, &
                                    pvalues_own(i_gene), ierr)
                if (is_err(ierr)) return
                neighborhood_size_own_case(i_gene) = n_pool_case
                neighborhood_size_own_control(i_gene) = n_pool_control_own
            end if

            neighborhood_size_case(i_gene) = n_pool_case
            ! Count genes that received an own p-value (compute_pvalue_own can be 0).
            if (pvalues_own(i_gene) >= 0.0_real64) &
                n_genes_with_pvalue = n_genes_with_pvalue + 1
        end do
    end subroutine compute_noise_pvalue_pipeline_helper

    !> Validate inputs, build sorted structures, and run the per-gene pipeline.
    !|
    !| This is the alloc-layer entry point for the full noise-model pipeline. It owns
    !| every allocation needed by the computation (the per-gene work arrays). This
    !| exact variant allocates no Monte Carlo / RNG buffer.
    !| Internally it:
    !|   1. Validates all dimension and range arguments via `tox_errors`.
    !|   2. Calls `prepare_sorted_data` for both case and control data.
    !|   3. Allocates all per-gene work arrays used by the per-gene loop.
    !|   4. Delegates the per-gene loop to `compute_noise_pvalue_pipeline_helper`.
    subroutine compute_noise_pvalue_pipeline( &
        means_case, replicates_case, n_genes_case, n_replicates_case, &
        means_control, replicates_control, n_genes_control, n_replicates_control, &
        observed_statistic_own, compute_pvalue_own, &
        n_genes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        ierr)

        integer(int32), intent(in) :: n_genes_case
        !! Number of genes in the case group
        integer(int32), intent(in) :: n_replicates_case
        !! Number of case replicates
        integer(int32), intent(in) :: n_genes_control
        !! Number of genes in the control group
        integer(int32), intent(in) :: n_replicates_control
        !! Number of control replicates
        integer(int32), intent(in) :: n_genes
        !! Total number of genes for which p-values are computed
        real(real64), dimension(n_genes_case), intent(in) :: means_case
        !! Per-gene case expression means
        real(real64), dimension(n_replicates_case, n_genes_case), intent(in) :: replicates_case
        !! Case replicate expression matrix (n_replicates_case x n_genes_case)
        real(real64), dimension(n_genes_control), intent(in) :: means_control
        !! Per-gene control expression means
        real(real64), dimension(n_replicates_control, n_genes_control), intent(in) :: replicates_control
        !! Control replicate expression matrix (n_replicates_control x n_genes_control)
        real(real64), dimension(n_genes), intent(in) :: observed_statistic_own
        !! Observed gene-vs-own statistic for each gene
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_own
        !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        integer(int32), intent(in) :: k_start
        !! Minimum pool size before adaptive stopping is applied
        integer(int32), intent(in) :: k_step
        !! Residuals added per adaptive round before re-evaluating
        integer(int32), intent(in) :: k_max
        !! Hard upper limit on residual pool size
        real(real64), intent(in) :: tau
        !! Relative-change threshold for adaptive pool growth
        real(real64), intent(in) :: trim_frac
        !! Symmetric per-tail residual-pool trim fraction in [0, 0.5); applied ONLY
        !! for raw normalization (norm_method == 0). 0 = no trimming.
        integer(int32), intent(in) :: max_pool_size
        !! Maximum number of residuals in any pool
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: gene-vs-own p-values (-1 if not computed)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes for which the own p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Case kNN-pool size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Control kNN-pool size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for each gene (-1 if not computed)
        integer(int32), intent(out) :: ierr
        !! Error code

        type(sorted_data_t) :: sorted_case, sorted_control
        real(real64), allocatable :: tmp_pool_case(:), tmp_pool_control_own(:)
        integer(int32) :: sort_ierr
        real(real64) :: effective_trim

        call set_ok(ierr)

        call validate_dimension_size(n_genes_case, ierr)
        call validate_dimension_size(n_replicates_case, ierr)
        call validate_dimension_size(n_genes_control, ierr)
        call validate_dimension_size(n_replicates_control, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(k_start, ierr)
        call validate_dimension_size(k_step, ierr)
        call validate_dimension_size(k_max, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        call validate_all_in_range_real(means_case, n_genes_case, ierr)
        call validate_all_in_range_real(means_control, n_genes_control, ierr)
        call validate_all_in_range_real(replicates_case, n_replicates_case * n_genes_case, ierr)
        call validate_all_in_range_real(replicates_control, n_replicates_control * n_genes_control, ierr)
        if (is_err(ierr)) return

        call prepare_sorted_data(means_case, replicates_case, &
                                 n_replicates_case, n_genes_case, norm_method, sorted_case, sort_ierr)
        call set_err(ierr, sort_ierr)

        call prepare_sorted_data(means_control, replicates_control, &
                                 n_replicates_control, n_genes_control, norm_method, sorted_control, sort_ierr)
        call set_err(ierr, sort_ierr)
        if (is_err(ierr)) return

        ! Per-gene work arrays for compute_noise_pvalue_pipeline_helper, allocated
        ! once here so the helper itself performs no allocation.
        M_ALLOCATE(tmp_pool_case(max_pool_size * 2))
        M_ALLOCATE(tmp_pool_control_own(max_pool_size * 2))

        if (is_err(ierr)) return

        ! Residual-pool trimming is a raw-normalization-only knob: log/voom-style
        ! transforms already stabilize the mean-variance trend, so there is nothing
        ! to trim there. Gate it here so the helper receives an already-resolved
        ! fraction (0 disables it).
        effective_trim = 0.0_real64
        if (norm_method == 0) effective_trim = trim_frac

        call compute_noise_pvalue_pipeline_helper( &
            sorted_case, sorted_control, &
            means_case, means_control, &
            observed_statistic_own, compute_pvalue_own, &
            n_genes, k_start, k_step, k_max, tau, effective_trim, &
            pvalues_own, n_genes_with_pvalue, &
            max_pool_size, &
            neighborhood_size_own_case, neighborhood_size_own_control, &
            neighborhood_size_case, &
            tmp_pool_case, tmp_pool_control_own, &
            ierr)

    end subroutine compute_noise_pvalue_pipeline

end module noise_model_exact

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

!> C-interoperable wrapper for the EXACT `compute_noise_pvalue_pipeline`.
!|
!| Identical ABI to `compute_noise_pvalues_pipeline_c` (baseline module) — same
!| arguments in the same order — but bound under a distinct C name
!| (`compute_noise_pvalues_pipeline_exact_c`) and dispatching to `noise_model_exact`,
!| so both can be linked and called side by side for comparison.
!|
!| Performs null-pointer checks via `M_CHECK_IERR_NON_NULL` / `M_CHECK_NON_NULL`,
!| then delegates unconditionally to the validated Fortran entry point.
!| No computation is performed here.
subroutine compute_noise_pvalues_pipeline_exact_c( &
    means_case, replicates_case, n_genes_case, n_replicates_case, &
    means_control, replicates_control, n_genes_control, n_replicates_control, &
    observed_statistic_own, compute_pvalue_own, &
    n_genes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
    pvalues_own, n_genes_with_pvalue, &
    max_pool_size, &
    neighborhood_size_own_case, neighborhood_size_own_control, &
    neighborhood_size_case, &
    ierr) bind(C, name="compute_noise_pvalues_pipeline_exact_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use noise_model_exact, only: compute_noise_pvalue_pipeline
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
    !! Total number of genes for which p-values are computed
    real(c_double), dimension(n_genes_case), intent(in), target :: means_case
    !! Per-gene case expression means
    real(c_double), dimension(n_replicates_case, n_genes_case), intent(in), target :: replicates_case
    !! Case replicate expression matrix (n_replicates_case x n_genes_case)
    real(c_double), dimension(n_genes_control), intent(in), target :: means_control
    !! Per-gene control expression means
    real(c_double), dimension(n_replicates_control, n_genes_control), intent(in), target :: replicates_control
    !! Control replicate expression matrix (n_replicates_control x n_genes_control)
    real(c_double), dimension(n_genes), intent(in), target :: observed_statistic_own
    !! Observed gene-vs-own statistic for each gene
    integer(c_int), dimension(n_genes), intent(in), target :: compute_pvalue_own
    !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
    integer(c_int), intent(in), target :: norm_method
    !! 0 = linear scale; non-zero = log2(x+1) transform
    integer(c_int), intent(in), target :: k_start
    !! Minimum pool size before adaptive stopping is applied
    integer(c_int), intent(in), target :: k_step
    !! Residuals added per adaptive round before re-evaluating
    integer(c_int), intent(in), target :: k_max
    !! Hard upper limit on residual pool size
    real(c_double), intent(in), target :: tau
    !! Relative-change threshold for adaptive pool growth
    real(c_double), intent(in), target :: trim_frac
    !! Symmetric per-tail residual-pool trim fraction in [0, 0.5); raw norm only
    integer(c_int), intent(in), target :: max_pool_size
    !! Maximum number of residuals in any pool
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_own
    !! Output: gene-vs-own p-values (-1 if not computed)
    integer(c_int), intent(out), target :: n_genes_with_pvalue
    !! Number of genes for which at least one p-value was computed
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_case
    !! Case kNN-pool size used for gene-vs-own (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_control
    !! Control kNN-pool size used for gene-vs-own (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_case
    !! Case pool size used for each gene (-1 if not computed)
    integer(c_int), intent(out), target :: ierr
    !! Error code: 0 = success

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes_case)
    M_CHECK_NON_NULL(n_replicates_case)
    M_CHECK_NON_NULL(n_genes_control)
    M_CHECK_NON_NULL(n_replicates_control)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(means_case)
    M_CHECK_NON_NULL(replicates_case)
    M_CHECK_NON_NULL(means_control)
    M_CHECK_NON_NULL(replicates_control)
    M_CHECK_NON_NULL(observed_statistic_own)
    M_CHECK_NON_NULL(compute_pvalue_own)
    M_CHECK_NON_NULL(norm_method)
    M_CHECK_NON_NULL(k_start)
    M_CHECK_NON_NULL(k_step)
    M_CHECK_NON_NULL(k_max)
    M_CHECK_NON_NULL(tau)
    M_CHECK_NON_NULL(trim_frac)
    M_CHECK_NON_NULL(max_pool_size)
    M_CHECK_NON_NULL(pvalues_own)
    M_CHECK_NON_NULL(n_genes_with_pvalue)
    M_CHECK_NON_NULL(neighborhood_size_own_case)
    M_CHECK_NON_NULL(neighborhood_size_own_control)
    M_CHECK_NON_NULL(neighborhood_size_case)

    call compute_noise_pvalue_pipeline( &
        means_case, replicates_case, n_genes_case, n_replicates_case, &
        means_control, replicates_control, n_genes_control, n_replicates_control, &
        observed_statistic_own, compute_pvalue_own, &
        n_genes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        ierr)

end subroutine compute_noise_pvalues_pipeline_exact_c