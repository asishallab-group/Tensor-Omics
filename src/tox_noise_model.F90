#include "macros.h"

!!! TO DO
!!! 1. Implement error code to argument number mapping
!!! 2. Find solution for family based analysis

!> Noise model for exact p-value computation comparing case and control gene expression.
!|
!| This module and `noise_model_exact` (in `tox_noise_model_exact.F90`) are
!| identical EXCEPT for how the `own` comparison's null is corrected for the
!| individual-vs-mean scale mismatch (the observed statistic is a difference of
!| MEANS, so the null must be built at the mean level):
!|   - `noise_model_exact`: scales each `own` stratum by `1/sqrt(n_replicates)`
!|     — a fast, variance-only correction of an individual-residual null.
!|   - this module: BOOTSTRAPS the mean difference (resample n_replicates per side
!|     with replacement, average, difference) — correct scale AND shape, no
!|     equal-variance assumption (see `compute_pvalue_bootstrap_mean_helper`).
!| Both share the same C ABI so either can be swapped in without downstream
!| changes; keeping them as two modules lets one be accepted and the other
!| discarded in a later phase.
!|
!| The per-gene p-value itself is computed EXACTLY at every pool size (sorted
!| control pool + two binary searches, `O((n+m) log m)`; no exhaustive-vs-Monte
!| Carlo split). The bootstrap `own` null is the one place this module uses the RNG.
!|
!| Provides routines to:
!|   - Sort and pack gene expression residuals for efficient neighbourhood lookup
!|   - Adaptively gather residual pools from the nearest-mean neighbourhood
!|   - Variance-stratify the gathered pool into quantile intervals so that the
!|     `own` comparison draws only from the gene's own variance regime, rather
!|     than from a mixture of regimes across the kNN neighbourhood
!|   - Compute exact p-values by counting, via sorted-pool binary search, every
!|     pairwise residual difference that meets or exceeds the observed statistic
!|   - Run the full pipeline over all genes with pre-cached family / ortholog pools
module noise_model
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_NAN_INF, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int
    use f42_utils, only: sort_real, sort_integer, calc_percentile_helper, init_random
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

    integer(int32), parameter :: STRATA_N_SCHEDULE_STEPS = 6
        !! Number of candidate bin-count steps tried by `stratify_residuals_helper`
    integer(int32), parameter :: STRATA_BIN_COUNT_SCHEDULE(STRATA_N_SCHEDULE_STEPS) = &
        [20_int32, 10_int32, 5_int32, 4_int32, 2_int32, 1_int32]
        !! Candidate quantile bin counts, finest (5% per bin) to coarsest first,
        !! tried in order until a step satisfies all acceptance criteria. The last
        !! entry (1 bin = no stratification, the whole pool) is a hard floor: it is
        !! accepted unconditionally if no earlier step succeeds.
    real(real64), parameter :: STRATA_MAX_C_G_PROB_THRESHOLD = 0.1_real64
        !! Criterion 2: `Pr(c_g > 2)` must stay below this fraction
    integer(int32), parameter :: STRATA_MIN_RESIDUALS_PER_BIN = 1000_int32
        !! Criterion 3: every occupied quantile interval must contain at least this
        !! many residuals to support a stable exact null distribution
    ! Note: the per-gene p-value is computed exactly at every pool size (no
    ! MAX_EXACT_COMBINATIONS enumeration/Monte-Carlo threshold). The only RNG use is
    ! the bootstrap that builds the `own` mean-difference null, below.

    integer(int32), parameter :: N_BOOTSTRAP_DRAWS = 10000_int32
        !! Number of bootstrap resamples used to build the `own` mean-difference null

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
        ! once here, at construction, so every downstream null (bootstrap / exact) is
        ! built from unbiased-variance residuals while the observed mean-difference
        ! statistic is left untouched. prepare_sorted_data enforces n_samples >= 2.
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
    !| Alongside the pooled residuals, `gene_id_per_residual` records the sorted-gene-slot
    !| (index into `sorted_data%means_sorted` / `sorted_data%original_indices`) that each
    !| pooled residual was taken from. This bookkeeping is required by variance
    !| stratification (see `stratify_residuals_helper`), which needs to know how many
    !| distinct genes contribute residuals to a given quantile interval.
    !|
    !| No allocation; `pooled_residuals` and `gene_id_per_residual` must be
    !| pre-allocated by the caller to at least `max_pool_size`. Candidate expansions
    !| for a round are staged in place within `pooled_residuals` past the committed
    !| size, so no separate staging buffer is needed.
    pure subroutine gather_residuals_helper(target_mean, sorted_data, &
                                            k_start, k_step, k_max, tau, &
                                            pooled_residuals, gene_id_per_residual, n_pooled, &
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
        integer(int32), intent(out) :: gene_id_per_residual(:)
        !! Output: sorted-gene-slot that each entry of `pooled_residuals` was taken from
        !! (pre-allocated to at least `max_pool_size`)
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
        gene_id_per_residual(1:current_size) = pos
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
            call add_gene_id_to_pool_helper(gene_id_per_residual, current_size, idx, pool_size)
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
                call add_gene_id_to_pool_helper(gene_id_per_residual, current_size + offset, idx, pool_size)
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

    !> Append the sorted-gene-slot id of one gene into `gene_ids`, mirroring the copy
    !| performed by `add_residuals_to_pool_helper` on the parallel residual pool.
    !|
    !| Both helpers must be called together with the same `current_size` / `new_size`
    !| so the two parallel arrays (`pool` and `gene_ids`) stay index-aligned.
    pure subroutine add_gene_id_to_pool_helper(gene_ids, current_size, gene_id, new_size)
        integer(int32), intent(inout) :: gene_ids(:)
        !! Target gene-id pool (pre-allocated by caller), parallel to a residual pool
        integer(int32), intent(in) :: current_size
        !! Number of elements already in `gene_ids`
        integer(int32), intent(in) :: gene_id
        !! Sorted-gene-slot id to repeat for every newly appended residual
        integer(int32), intent(in) :: new_size
        !! Maximum number of elements allowed in `gene_ids` after the append

        integer(int32) :: n_to_fill

        n_to_fill = new_size - current_size
        if (n_to_fill > 0) gene_ids(current_size + 1:new_size) = gene_id
    end subroutine add_gene_id_to_pool_helper

    ! =========================================================================
    ! stratify_residuals
    ! =========================================================================

    !> Core implementation: assign each pooled residual to a quantile bin.
    !|
    !| Divides the pooled residuals into `n_bins` quantile intervals (equal
    !| *frequency*, not equal width): cutpoints are placed at the
    !| `0, 100/n_bins, 200/n_bins, ..., 100` percentiles of `pooled_residuals`,
    !| computed via `f42_utils`'s `calc_percentile_helper` (linear interpolation,
    !| matching `numpy.percentile`). Every bin therefore contains the same number of
    !| residuals (up to ties and interpolation at the edges), regardless of how the
    !| underlying residual values are distributed — unlike equal-width binning, a
    !| bin can never come up short on residuals purely because the gene that
    !| dominates it has a narrow value range.
    !|
    !| No allocation; `bin_index_per_residual` must be pre-allocated by the caller.
    !| The ascending sort permutation `perm` is supplied by the caller
    !| (`stratify_residuals_helper`) and computed once per pool, since it is
    !| identical for every candidate bin count — recomputing it here per schedule
    !| step would repeat the same `O(n log n)` sort up to `STRATA_N_SCHEDULE_STEPS`
    !| times on identical data.
    pure subroutine assign_residual_bins_helper(bin_values, perm, n_pooled, n_bins, &
                                                bin_index_per_residual, bin_edges)
        integer(int32), intent(in) :: n_pooled
        !! Number of valid entries in `bin_values`
        real(real64), dimension(n_pooled), intent(in) :: bin_values
        !! Per-residual quantity that DEFINES the bins. To stratify by expression /
        !! variance regime we pass the SOURCE GENE MEAN of each residual (not the
        !! residual value), so every residual of a gene lands in a single bin
        !! (c_g == 1 by construction) and the stratum for a target gene can be selected
        !! by its own mean.
        integer(int32), dimension(n_pooled), intent(in) :: perm
        !! Ascending sort permutation of `bin_values` (pre-computed by the caller)
        integer(int32), intent(in) :: n_bins
        !! Number of quantile bins to divide the values into
        integer(int32), dimension(n_pooled), intent(out) :: bin_index_per_residual
        !! Output: 1-based bin index for each residual
        real(real64), dimension(n_bins + 1), intent(out) :: bin_edges
        !! Output: the `n_bins + 1` quantile cutpoints, from the 0th to the 100th percentile

        integer(int32) :: i_resid, i_edge
        real(real64) :: percentile_step

        percentile_step = 100.0_real64 / real(n_bins, real64)
        do i_edge = 1, n_bins + 1
            call calc_percentile_helper(bin_values, perm, &
                                        percentile_step * real(i_edge - 1, real64), bin_edges(i_edge))
        end do

        if (bin_edges(n_bins + 1) <= bin_edges(1)) then
            bin_index_per_residual = 1
            return
        end if

        do concurrent(i_resid=1:n_pooled) shared(bin_index_per_residual, bin_values, bin_edges, n_bins)
            bin_index_per_residual(i_resid) = locate_bin_helper(bin_values(i_resid), bin_edges, n_bins)
        end do
    end subroutine assign_residual_bins_helper

    !> Locate the bin that `target_value` falls into, given pre-computed `bin_edges`.
    !|
    !| `bin_edges` must be the `n_bins + 1` ascending boundaries produced by
    !| `assign_residual_bins_helper` for the same stratification. Values outside
    !| `[bin_edges(1), bin_edges(n_bins+1)]` are clamped to the nearest end bin.
    pure function locate_bin_helper(target_value, bin_edges, n_bins) result(bin_idx)
        integer(int32), intent(in) :: n_bins
        !! Number of bins (i.e. `size(bin_edges) - 1`)
        real(real64), dimension(n_bins + 1), intent(in) :: bin_edges
        !! Ascending bin boundaries
        real(real64), intent(in) :: target_value
        !! Value to locate within the binning
        integer(int32) :: bin_idx
        !! 1-based index of the bin containing `target_value`

        integer(int32) :: i_bin

        if (target_value <= bin_edges(1)) then
            bin_idx = 1
            return
        end if
        if (target_value >= bin_edges(n_bins + 1)) then
            bin_idx = n_bins
            return
        end if

        bin_idx = n_bins
        do i_bin = 1, n_bins
            if (target_value < bin_edges(i_bin + 1)) then
                bin_idx = i_bin
                exit
            end if
        end do
    end function locate_bin_helper

    !> Core implementation: evaluate the three stratification acceptance criteria.
    !|
    !| Given a per-residual bin assignment and the originating sorted-gene-slot of
    !| each residual, computes, per distinct gene present in the pool:
    !|
    !|   `c_g = max_bin(gene) - min_bin(gene) + 1`
    !|
    !| (bins are contiguous in mean-sorted order, so the span between the gene's
    !| smallest and largest occupied bin equals the count of distinct bins it touches).
    !| The stratification is accepted when all of the following hold:
    !|
    !|   1. `median(c_g) == 1`              (biological coherence)
    !|   2. `Pr(c_g > 2) < 0.1`              (tail control)
    !|   3. every occupied bin has >= 50 residuals  (sampling stability)
    !|
    !| No allocation; `tmp_gene_min_bin`, `tmp_gene_max_bin`, `tmp_gene_seen`,
    !| `tmp_touched_gene_slots`, `tmp_c_g`, and `tmp_bin_counts` are work arrays
    !| pre-allocated by the caller.
    !|
    !| `tmp_gene_seen` must be all-`.false.` on entry. Rather than clear and rescan
    !| the entire `n_genes_total`-length universe on every call (this routine runs
    !| up to `STRATA_N_SCHEDULE_STEPS` times per pool, per side, per gene — so a full
    !| clear/scan makes the whole pipeline scale as `O(n_genes^2)`), it records the
    !| distinct gene slots it touches in `tmp_touched_gene_slots` (at most `n_pooled`
    !| of them) and, before returning, resets exactly those `tmp_gene_seen` entries
    !| back to `.false.`. The array is therefore left all-`.false.` again for the
    !| next call, and every operation here costs `O(n_pooled)`, never `O(n_genes_total)`.
    pure subroutine check_stratification_accepted_helper( &
        bin_index_per_residual, gene_id_per_residual, n_pooled, n_bins, n_genes_total, &
        tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, tmp_touched_gene_slots, tmp_c_g, tmp_bin_counts, &
        is_accepted)
        integer(int32), intent(in) :: n_pooled
        !! Number of pooled residuals
        integer(int32), intent(in) :: n_bins
        !! Number of bins in this stratification candidate
        integer(int32), intent(in) :: n_genes_total
        !! Total number of genes in the underlying sorted structure (upper bound on gene-slot ids)
        integer(int32), dimension(n_pooled), intent(in) :: bin_index_per_residual
        !! 1-based bin index for each pooled residual
        integer(int32), dimension(n_pooled), intent(in) :: gene_id_per_residual
        !! Sorted-gene-slot id (1..n_genes_total) that each pooled residual came from
        integer(int32), dimension(n_genes_total), intent(inout) :: tmp_gene_min_bin
        !! Work array: smallest occupied bin index seen so far, per gene slot
        integer(int32), dimension(n_genes_total), intent(inout) :: tmp_gene_max_bin
        !! Work array: largest occupied bin index seen so far, per gene slot
        logical, dimension(n_genes_total), intent(inout) :: tmp_gene_seen
        !! Work array: `.true.` once a gene slot has contributed a residual during this
        !! call; MUST be all-`.false.` on entry and is left all-`.false.` on return
        integer(int32), dimension(n_pooled), intent(inout) :: tmp_touched_gene_slots
        !! Work array: distinct gene slots touched by this pool (first `n_distinct_genes`
        !! entries valid); used to reset `tmp_gene_seen` without an `O(n_genes_total)` sweep
        integer(int32), dimension(n_pooled), intent(inout) :: tmp_c_g
        !! Work array: per-distinct-gene `c_g` values (only the first `n_distinct_genes` entries are used)
        integer(int32), dimension(n_bins), intent(inout) :: tmp_bin_counts
        !! Work array: number of residuals occupying each bin
        logical, intent(out) :: is_accepted
        !! `.true.` if all three acceptance criteria are satisfied

        integer(int32) :: i_resid, gene_id, bin_id, n_distinct_genes, i_distinct
        integer(int32) :: n_c_g_gt_2, median_c_g, n_occupied_bins, min_occupied_count

        tmp_bin_counts = 0
        n_distinct_genes = 0

        do i_resid = 1, n_pooled
            gene_id = gene_id_per_residual(i_resid)
            bin_id = bin_index_per_residual(i_resid)

            tmp_bin_counts(bin_id) = tmp_bin_counts(bin_id) + 1

            if (.not. tmp_gene_seen(gene_id)) then
                tmp_gene_seen(gene_id) = .true.
                tmp_gene_min_bin(gene_id) = bin_id
                tmp_gene_max_bin(gene_id) = bin_id
                n_distinct_genes = n_distinct_genes + 1
                tmp_touched_gene_slots(n_distinct_genes) = gene_id
            else
                tmp_gene_min_bin(gene_id) = min(tmp_gene_min_bin(gene_id), bin_id)
                tmp_gene_max_bin(gene_id) = max(tmp_gene_max_bin(gene_id), bin_id)
            end if
        end do

        ! Build c_g from the touched slots only, and reset tmp_gene_seen for reuse.
        do i_distinct = 1, n_distinct_genes
            gene_id = tmp_touched_gene_slots(i_distinct)
            tmp_c_g(i_distinct) = tmp_gene_max_bin(gene_id) - tmp_gene_min_bin(gene_id) + 1
            tmp_gene_seen(gene_id) = .false.
        end do

        ! ── Criterion 1: median(c_g) == 1 ────────────────────────────────────
        call median_of_int_array_helper(tmp_c_g(1:n_distinct_genes), n_distinct_genes, median_c_g)

        ! ── Criterion 2: Pr(c_g > 2) < 0.1 ───────────────────────────────────
        n_c_g_gt_2 = count(tmp_c_g(1:n_distinct_genes) > 2)

        ! ── Criterion 3: every occupied bin has >= STRATA_MIN_RESIDUALS_PER_BIN ──
        n_occupied_bins = count(tmp_bin_counts > 0)
        if (n_occupied_bins > 0) then
            min_occupied_count = minval(tmp_bin_counts, mask=(tmp_bin_counts > 0))
        else
            min_occupied_count = 0
        end if

        is_accepted = (median_c_g == 1) .and. &
                      (real(n_c_g_gt_2, real64) / real(max(n_distinct_genes, 1_int32), real64) &
                       < STRATA_MAX_C_G_PROB_THRESHOLD) .and. &
                      (min_occupied_count >= STRATA_MIN_RESIDUALS_PER_BIN)
    end subroutine check_stratification_accepted_helper

    !> Compute the median of an integer array using a caller-supplied permutation buffer.
    !|
    !| Sorts a local copy via the project's `f42_utils` quicksort (indirect, through a
    !! permutation) rather than re-implementing sorting logic, per the project convention
    !| of reusing `f42_utils` for general-purpose routines.
    pure subroutine median_of_int_array_helper(values, n, median_val)
        integer(int32), intent(in) :: n
        !! Number of elements in `values`
        integer(int32), dimension(n), intent(in) :: values
        !! Input integer array (unsorted)
        integer(int32), intent(out) :: median_val
        !! Output: the median value (lower of the two middle values for even `n`)

        integer(int32) :: perm(n), stack_left(n), stack_right(n), i

        if (n == 0) then
            median_val = 0
            return
        end if

        do concurrent(i=1:n) shared(perm)
            perm(i) = i
        end do

        call sort_integer(values, perm, stack_left, stack_right)
        median_val = values(perm((n + 1) / 2))
    end subroutine median_of_int_array_helper

    !> Core implementation: find the coarsest-to-finest accepted variance stratification.
    !|
    !| Tries each candidate bin count in `STRATA_BIN_COUNT_SCHEDULE` (finest/5% first)
    !| and accepts the first one that satisfies all three criteria in
    !| `check_stratification_accepted_helper`. If none of the finer steps are accepted,
    !| the coarsest schedule entry (50%, 2 bins) is used unconditionally as a hard floor.
    !|
    !| No allocation; all `tmp_*` arrays are pre-allocated work arrays sized for the
    !| largest schedule step (`max_bins = STRATA_BIN_COUNT_SCHEDULE(1)`). `tmp_gene_seen`
    !| must be all-`.false.` on entry (see `check_stratification_accepted_helper`); it is
    !| left all-`.false.` on return.
    !|
    !| The pooled residuals are sorted once, up front: the resulting permutation is
    !| shared across every candidate bin count, since the quantile cutpoints for all
    !| schedule steps are read off the same ordering.
    pure subroutine stratify_residuals_helper( &
        means_sorted, gene_id_per_residual, n_pooled, n_genes_total, &
        tmp_bin_index_per_residual, tmp_bin_edges, &
        tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, tmp_touched_gene_slots, tmp_c_g, tmp_bin_counts, &
        chosen_n_bins, chosen_bin_index_per_residual, chosen_bin_edges, criteria_met)
        integer(int32), intent(in) :: n_pooled
        !! Number of pooled residuals to stratify
        integer(int32), intent(in) :: n_genes_total
        !! Total number of genes in the underlying sorted structure
        real(real64), dimension(n_genes_total), intent(in) :: means_sorted
        !! Sorted gene means (`sorted_data%means_sorted`). Bins are quantiles of the
        !! per-residual SOURCE GENE MEAN, so residuals are stratified by expression /
        !! variance regime rather than by residual value.
        integer(int32), dimension(n_pooled), intent(in) :: gene_id_per_residual
        !! Sorted-gene-slot id that each pooled residual came from
        integer(int32), dimension(n_pooled), intent(inout) :: tmp_bin_index_per_residual
        !! Work array: bin index per residual for the candidate being evaluated
        real(real64), dimension(STRATA_BIN_COUNT_SCHEDULE(1) + 1), intent(inout) :: tmp_bin_edges
        !! Work array: bin edges for the candidate being evaluated
        integer(int32), dimension(n_genes_total), intent(inout) :: tmp_gene_min_bin
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(n_genes_total), intent(inout) :: tmp_gene_max_bin
        !! Work array, see `check_stratification_accepted_helper`
        logical, dimension(n_genes_total), intent(inout) :: tmp_gene_seen
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(n_pooled), intent(inout) :: tmp_touched_gene_slots
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(n_pooled), intent(inout) :: tmp_c_g
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(STRATA_BIN_COUNT_SCHEDULE(1)), intent(inout) :: tmp_bin_counts
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), intent(out) :: chosen_n_bins
        !! Number of bins in the accepted stratification
        integer(int32), dimension(n_pooled), intent(out) :: chosen_bin_index_per_residual
        !! Bin index per residual for the accepted stratification
        real(real64), dimension(STRATA_BIN_COUNT_SCHEDULE(1) + 1), intent(out) :: chosen_bin_edges
        !! Bin edges for the accepted stratification (only the first `chosen_n_bins + 1` are valid)
        logical, intent(out) :: criteria_met
        !! `.true.` if an earlier (finer) schedule step was accepted on its own merit;
        !! `.false.` if the coarsest step was used only as a hard floor fallback

        integer(int32) :: i_step, n_bins, i_resid
        integer(int32) :: perm(n_pooled), stack_left(n_pooled), stack_right(n_pooled)
        real(real64) :: gene_mean_per_residual(n_pooled)
        logical :: is_accepted

        criteria_met = .false.

        ! Stratify by the SOURCE GENE MEAN of each residual (its expression / variance
        ! regime), NOT by the residual value. Every residual of a gene then shares one
        ! bin (c_g == 1 by construction), so the acceptance criteria become meaningful
        ! and select_stratum_for_target -- which locates a gene MEAN in these edges -- is
        ! on the correct scale.
        do concurrent(i_resid=1:n_pooled) &
            shared(gene_mean_per_residual, means_sorted, gene_id_per_residual)
            gene_mean_per_residual(i_resid) = means_sorted(gene_id_per_residual(i_resid))
        end do

        ! Sort once by gene mean: the permutation is identical for every candidate bin count.
        do concurrent(i_resid=1:n_pooled) shared(perm)
            perm(i_resid) = i_resid
        end do
        call sort_real(gene_mean_per_residual, perm, stack_left, stack_right)

        do i_step = 1, STRATA_N_SCHEDULE_STEPS
            n_bins = STRATA_BIN_COUNT_SCHEDULE(i_step)

            call assign_residual_bins_helper(gene_mean_per_residual, perm, n_pooled, n_bins, &
                                             tmp_bin_index_per_residual, tmp_bin_edges(1:n_bins + 1))

            call check_stratification_accepted_helper( &
                tmp_bin_index_per_residual, gene_id_per_residual, n_pooled, n_bins, n_genes_total, &
                tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, tmp_touched_gene_slots, &
                tmp_c_g, tmp_bin_counts(1:n_bins), &
                is_accepted)

            if (is_accepted .or. i_step == STRATA_N_SCHEDULE_STEPS) then
                chosen_n_bins = n_bins
                chosen_bin_index_per_residual = tmp_bin_index_per_residual
                chosen_bin_edges = tmp_bin_edges
                criteria_met = is_accepted
                return
            end if
        end do
    end subroutine stratify_residuals_helper

    !> Sample only the residuals from the bin containing `target_value`.
    !|
    !| Given an accepted stratification (`bin_index_per_residual`, `bin_edges`,
    !| `n_bins`), locates which bin `target_value` falls into and copies that bin's
    !| residuals into `stratum_residuals`. No allocation; `stratum_residuals` must be
    !| pre-allocated to at least `n_pooled`.
    pure subroutine select_stratum_for_target_helper( &
        pooled_residuals, bin_index_per_residual, n_pooled, bin_edges, n_bins, target_value, &
        stratum_residuals, n_stratum)
        integer(int32), intent(in) :: n_pooled
        !! Number of pooled residuals
        integer(int32), intent(in) :: n_bins
        !! Number of bins in the stratification
        real(real64), dimension(n_pooled), intent(in) :: pooled_residuals
        !! Pooled residual values
        integer(int32), dimension(n_pooled), intent(in) :: bin_index_per_residual
        !! Bin index per residual, from `stratify_residuals_helper`
        real(real64), dimension(n_bins + 1), intent(in) :: bin_edges
        !! Bin edges, from `stratify_residuals_helper`
        real(real64), intent(in) :: target_value
        !! Value (typically the gene's own mean) used to select the stratum
        real(real64), dimension(n_pooled), intent(out) :: stratum_residuals
        !! Output: residuals belonging to the selected stratum (first `n_stratum` entries valid)
        integer(int32), intent(out) :: n_stratum
        !! Number of residuals written into `stratum_residuals`

        integer(int32) :: target_bin, i_resid

        target_bin = locate_bin_helper(target_value, bin_edges, n_bins)

        n_stratum = 0
        do i_resid = 1, n_pooled
            if (bin_index_per_residual(i_resid) == target_bin) then
                n_stratum = n_stratum + 1
                stratum_residuals(n_stratum) = pooled_residuals(i_resid)
            end if
        end do
    end subroutine select_stratum_for_target_helper

    !> Bootstrap the mean-difference null and return the p-value.
    !|
    !| Builds the `own` null at the level the observed statistic actually lives on —
    !| a difference of MEANS — rather than at the individual-residual level. For each
    !| of `n_boot` draws it resamples `n_rep_case` residuals from `pool_case` and
    !| `n_rep_control` from `pool_control` (both with replacement), averages each
    !| side, and takes the absolute difference. The p-value is the add-one-corrected
    !| fraction of those null statistics that meet or exceed `observed_statistic_abs`.
    !|
    !| This gives the correct `sigma^2 / n_rep` scaling AND the correct (CLT-thinned)
    !| tail shape, and — because each side is resampled from its own pool at its own
    !| `n_rep` — it makes no equal-variance assumption. `n_rep_case` / `n_rep_control`
    !| are the replicate counts the observed case / control means were averaged over
    !| (i.e. `sorted_*%max_resid_per_gene`). Indices are drawn a whole side at a time
    !| via one `random_number` array fill (not per-scalar) for speed; impure, so the
    !| global RNG must be seeded (the pipeline calls `init_random(42)` up front for
    !| reproducibility).
    subroutine compute_pvalue_bootstrap_mean_helper(pool_case, n_pool_case, &
                                                    pool_control, n_pool_control, &
                                                    n_rep_case, n_rep_control, &
                                                    observed_statistic_abs, n_boot, &
                                                    p_value)
        integer(int32), intent(in) :: n_pool_case
        !! Size of the case residual pool (the sampling source)
        real(real64), dimension(n_pool_case), intent(in) :: pool_case
        !! Case residual pool
        integer(int32), intent(in) :: n_pool_control
        !! Size of the control residual pool (the sampling source)
        real(real64), dimension(n_pool_control), intent(in) :: pool_control
        !! Control residual pool
        integer(int32), intent(in) :: n_rep_case
        !! Number of case replicates the observed mean averages over (resample size)
        integer(int32), intent(in) :: n_rep_control
        !! Number of control replicates the observed mean averages over (resample size)
        real(real64), intent(in) :: observed_statistic_abs
        !! Absolute value of the observed test statistic
        integer(int32), intent(in) :: n_boot
        !! Number of bootstrap resamples
        real(real64), intent(out) :: p_value
        !! Bootstrapped p-value

        integer(int32) :: i_boot, i_rep, idx, count_ge
        real(real64) :: sum_case, sum_control, null_stat
        real(real64) :: rbuf(max(n_rep_case, n_rep_control))

        ! Draw the resample indices for a whole side in a SINGLE array RNG call and
        ! map [0,1) -> [1, n_pool]. This replaces the millions of scalar `rand_range`
        ! calls per gene (n_boot * n_rep of them) with 2*n_boot array fills, which is
        ! far cheaper per value. `init_random(42)` (called once by the pipeline) still
        ! seeds this stream, so results stay reproducible.
        count_ge = 0
        do i_boot = 1, n_boot
            call random_number(rbuf(1:n_rep_case))
            sum_case = 0.0_real64
            do i_rep = 1, n_rep_case
                idx = min(int(rbuf(i_rep) * real(n_pool_case, real64), int32) + 1, n_pool_case)
                sum_case = sum_case + pool_case(idx)
            end do

            call random_number(rbuf(1:n_rep_control))
            sum_control = 0.0_real64
            do i_rep = 1, n_rep_control
                idx = min(int(rbuf(i_rep) * real(n_pool_control, real64), int32) + 1, n_pool_control)
                sum_control = sum_control + pool_control(idx)
            end do

            null_stat = abs(sum_case / real(n_rep_case, real64) - &
                            sum_control / real(n_rep_control, real64))
            if (null_stat >= observed_statistic_abs) count_ge = count_ge + 1
        end do

        p_value = real(count_ge + 1, real64) / real(n_boot + 1, real64)
    end subroutine compute_pvalue_bootstrap_mean_helper

    ! =========================================================================
    ! compute_noise_pvalue_pipeline
    ! =========================================================================

    !> Core pipeline: compute per-gene noise p-values using pre-built data structures.
    !|
    !| Iterates over all genes and computes up to three p-values each:
    !|   - `pvalues_own`:  gene vs. its own matched group B neighbourhood, using a
    !|     variance-stratified residual pool (see `stratify_residuals_helper`) for
    !|     both case and control pools
    !|   - `pvalues_family`:  gene vs. its gene-family pool in group B (unchanged kNN pooling)
    !|   - `pvalues_ortholog`: gene vs. its ortholog pool in control (unchanged kNN pooling)
    !|
    !| For the `own` comparison, both the case and control kNN neighbourhoods are
    !| independently partitioned into variance strata (quantile bins over the
    !| residual distribution) before computing the p-value. This addresses the loss
    !| of sensitivity caused by mixing residuals from genes with different variance
    !| regimes within the same kNN neighbourhood. Only the stratum containing the
    !| gene's own mean (case mean for case pool, control mean for control pool) is
    !| used. The `family` and `ortholog` comparisons retain the original whole-pool
    !| kNN sampling, since family/ortholog-level stratification handling is out of
    !| scope for this change.
    !|
    !| Requires `sorted_case` and `sorted_control` to already be built (via
    !| `prepare_sorted_data`) and `cache` to be pre-populated (all `is_cached`
    !| entries set for families with `family_sizes > 0`).
    !|
    !| No allocation: every work array (the per-gene residual pools and the
    !| variance-stratification scratch arrays) is pre-allocated by the caller and
    !| passed in as an `intent(inout)` work array. There is no Monte Carlo buffer —
    !| the p-value is computed exactly at every pool size by `compute_pvalue`. The
    !| `own` mean-difference null is bootstrapped in place from the strata (via
    !| `compute_pvalue_bootstrap_mean_helper`), which draws from the global RNG; that
    !| is the only impure operation, so this subroutine is not `pure`.
    subroutine compute_noise_pvalue_pipeline_helper( &
        sorted_case, sorted_control, &
        means_case, means_control, &
        observed_statistic_own, &
        compute_pvalue_own, &
        n_genes, k_start, k_step, k_max, tau, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        chosen_n_bins_own_case, chosen_n_bins_own_control, &
        tmp_pool_case, tmp_gene_id_pool_case, &
        tmp_pool_control_own, tmp_gene_id_pool_control_own, &
        tmp_strat_bin_index_per_residual_case, tmp_chosen_bin_index_per_residual_case, &
        tmp_strat_bin_edges_case, tmp_chosen_bin_edges_case, &
        tmp_strat_gene_min_bin_case, tmp_strat_gene_max_bin_case, &
        tmp_strat_gene_seen_case, tmp_strat_touched_gene_slots_case, tmp_strat_c_g_case, tmp_strat_bin_counts_case, &
        tmp_strat_bin_index_per_residual_control, tmp_chosen_bin_index_per_residual_control, &
        tmp_strat_bin_edges_control, tmp_chosen_bin_edges_control, &
        tmp_strat_gene_min_bin_control, tmp_strat_gene_max_bin_control, &
        tmp_strat_gene_seen_control, tmp_strat_touched_gene_slots_control, tmp_strat_c_g_control, tmp_strat_bin_counts_control, &
        tmp_own_stratum_pool_case, tmp_own_stratum_pool_control, &
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
        integer(int32), intent(in) :: max_pool_size
        !! Allocated size of all pool arrays
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: gene-vs-own p-values (-1 if not computed)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes for which the own p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Stratum size used for the gene-vs-own case pool (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Stratum size used for the gene-vs-own control pool (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for this gene (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_case
        !! Diagnostic: chosen bin count for the gene-vs-own CASE stratification, sign-encoded
        !! (+ = criteria met, - = coarse 2-bin fallback; -1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_control
        !! Diagnostic: chosen bin count for the gene-vs-own CONTROL stratification, sign-encoded
        !! (+ = criteria met, - = coarse 2-bin fallback; -1 if not computed)
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool_case
        !! Work array: output residual pool for this gene's case kNN neighbourhood
        integer(int32), dimension(max_pool_size * 2), intent(inout) :: tmp_gene_id_pool_case
        !! Work array: output gene-id pool for this gene's case kNN neighbourhood
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool_control_own
        !! Work array: output residual pool for this gene's control kNN neighbourhood
        !! (the `own` comparison, before variance stratification)
        integer(int32), dimension(max_pool_size * 2), intent(inout) :: tmp_gene_id_pool_control_own
        !! Work array: output gene-id pool for this gene's control kNN neighbourhood (the `own` comparison)
        ! Case stratification work arrays
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_strat_bin_index_per_residual_case
        !! Work array: bin index per residual for case stratification candidate
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_chosen_bin_index_per_residual_case
        !! Work array: bin index per residual for accepted case stratification
        real(real64), dimension(STRATA_BIN_COUNT_SCHEDULE(1) + 1), intent(inout) :: tmp_strat_bin_edges_case
        !! Work array: bin edges for case stratification candidate
        real(real64), dimension(STRATA_BIN_COUNT_SCHEDULE(1) + 1), intent(inout) :: tmp_chosen_bin_edges_case
        !! Work array: bin edges for accepted case stratification
        integer(int32), dimension(sorted_case%n_genes), intent(inout) :: tmp_strat_gene_min_bin_case
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(sorted_case%n_genes), intent(inout) :: tmp_strat_gene_max_bin_case
        !! Work array, see `check_stratification_accepted_helper`
        logical, dimension(sorted_case%n_genes), intent(inout) :: tmp_strat_gene_seen_case
        !! Work array, see `check_stratification_accepted_helper` (all-`.false.` on entry)
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_strat_touched_gene_slots_case
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_strat_c_g_case
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(STRATA_BIN_COUNT_SCHEDULE(1)), intent(inout) :: tmp_strat_bin_counts_case
        !! Work array, see `check_stratification_accepted_helper`
        ! Control stratification work arrays
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_strat_bin_index_per_residual_control
        !! Work array: bin index per residual for control stratification candidate
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_chosen_bin_index_per_residual_control
        !! Work array: bin index per residual for accepted control stratification
        real(real64), dimension(STRATA_BIN_COUNT_SCHEDULE(1) + 1), intent(inout) :: tmp_strat_bin_edges_control
        !! Work array: bin edges for control stratification candidate
        real(real64), dimension(STRATA_BIN_COUNT_SCHEDULE(1) + 1), intent(inout) :: tmp_chosen_bin_edges_control
        !! Work array: bin edges for accepted control stratification
        integer(int32), dimension(sorted_control%n_genes), intent(inout) :: tmp_strat_gene_min_bin_control
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(sorted_control%n_genes), intent(inout) :: tmp_strat_gene_max_bin_control
        !! Work array, see `check_stratification_accepted_helper`
        logical, dimension(sorted_control%n_genes), intent(inout) :: tmp_strat_gene_seen_control
        !! Work array, see `check_stratification_accepted_helper` (all-`.false.` on entry)
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_strat_touched_gene_slots_control
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_strat_c_g_control
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(STRATA_BIN_COUNT_SCHEDULE(1)), intent(inout) :: tmp_strat_bin_counts_control
        !! Work array, see `check_stratification_accepted_helper`
        real(real64), dimension(max_pool_size), intent(inout) :: tmp_own_stratum_pool_case
        !! Work array: case residuals restricted to the stratum containing the gene's own case mean
        real(real64), dimension(max_pool_size), intent(inout) :: tmp_own_stratum_pool_control
        !! Work array: control residuals restricted to the stratum containing the gene's own control mean
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_gene
        real(real64) :: mean_case_val, mean_control_val
        real(real64) :: observed_statistic_own_val
        integer(int32) :: n_pool_case, n_pool_control_own
        integer(int32) :: case_stratum_count, control_stratum_count
        integer(int32) :: chosen_n_bins_case, chosen_n_bins_control
        logical :: strata_criteria_met_case, strata_criteria_met_control

        call set_ok(ierr)

        pvalues_own = -1.0_real64
        neighborhood_size_own_case = -1
        neighborhood_size_own_control = -1
        neighborhood_size_case = -1
        chosen_n_bins_own_case = -1
        chosen_n_bins_own_control = -1
        n_genes_with_pvalue = 0

        ! This loop MUST stay sequential: the `own` mean-difference bootstrap
        ! (compute_pvalue_bootstrap_mean_helper) draws from the global RNG stream
        ! seeded once by init_random(42). Parallelising the gene loop would make the
        ! draws order-dependent and break reproducibility -- use a per-gene
        ! counter-based RNG (seeded by gene index) before ever doing so.
        do i_gene = 1, n_genes
            mean_case_val = means_case(i_gene)
            mean_control_val = means_control(i_gene)

            call gather_residuals_helper(mean_case_val, sorted_case, &
                                         k_start, k_step, k_max, tau, &
                                         tmp_pool_case, tmp_gene_id_pool_case, n_pool_case, &
                                         max_pool_size)
            if (n_pool_case < 10) cycle

            call gather_residuals_helper(mean_control_val, sorted_control, &
                                         k_start, k_step, k_max, tau, &
                                         tmp_pool_control_own, tmp_gene_id_pool_control_own, n_pool_control_own, &
                                         max_pool_size)
            if (n_pool_control_own < 10) cycle

            observed_statistic_own_val = observed_statistic_own(i_gene)

            ! Skip genes with a non-finite observed statistic
            if (observed_statistic_own_val /= observed_statistic_own_val) cycle

            if (compute_pvalue_own(i_gene) == 1) then
                ! Variance-stratify both case and control `own` neighbourhoods,
                ! then restrict sampling to the stratum that contains this gene's
                ! own mean (case mean for case pool, control mean for control pool).

                ! Stratify case pool (bins are quantiles of the source-gene mean).
                call stratify_residuals_helper( &
                    sorted_case%means_sorted, &
                    tmp_gene_id_pool_case(1:n_pool_case), &
                    n_pool_case, sorted_case%n_genes, &
                    tmp_strat_bin_index_per_residual_case(1:n_pool_case), tmp_strat_bin_edges_case, &
                    tmp_strat_gene_min_bin_case, tmp_strat_gene_max_bin_case, tmp_strat_gene_seen_case, &
                    tmp_strat_touched_gene_slots_case(1:n_pool_case), &
                    tmp_strat_c_g_case(1:n_pool_case), tmp_strat_bin_counts_case, &
                    chosen_n_bins_case, tmp_chosen_bin_index_per_residual_case(1:n_pool_case), &
                    tmp_chosen_bin_edges_case, strata_criteria_met_case)

                call select_stratum_for_target_helper( &
                    tmp_pool_case(1:n_pool_case), &
                    tmp_chosen_bin_index_per_residual_case(1:n_pool_case), &
                    n_pool_case, &
                    tmp_chosen_bin_edges_case(1:chosen_n_bins_case + 1), chosen_n_bins_case, mean_case_val, &
                    tmp_own_stratum_pool_case(1:n_pool_case), case_stratum_count)

                ! Stratify control pool (bins are quantiles of the source-gene mean).
                call stratify_residuals_helper( &
                    sorted_control%means_sorted, &
                    tmp_gene_id_pool_control_own(1:n_pool_control_own), &
                    n_pool_control_own, sorted_control%n_genes, &
                    tmp_strat_bin_index_per_residual_control(1:n_pool_control_own), tmp_strat_bin_edges_control, &
                    tmp_strat_gene_min_bin_control, tmp_strat_gene_max_bin_control, tmp_strat_gene_seen_control, &
                    tmp_strat_touched_gene_slots_control(1:n_pool_control_own), &
                    tmp_strat_c_g_control(1:n_pool_control_own), tmp_strat_bin_counts_control, &
                    chosen_n_bins_control, tmp_chosen_bin_index_per_residual_control(1:n_pool_control_own), &
                    tmp_chosen_bin_edges_control, strata_criteria_met_control)

                call select_stratum_for_target_helper( &
                    tmp_pool_control_own(1:n_pool_control_own), &
                    tmp_chosen_bin_index_per_residual_control(1:n_pool_control_own), &
                    n_pool_control_own, &
                    tmp_chosen_bin_edges_control(1:chosen_n_bins_control + 1), chosen_n_bins_control, mean_control_val, &
                    tmp_own_stratum_pool_control(1:n_pool_control_own), control_stratum_count)

                ! Diagnostics: record the chosen bin count per side, sign-encoding whether
                ! the stratification criteria were met (+ = met, - = coarse 2-bin fallback).
                ! Recorded even if the stratum is then too small to compute a p-value.
                chosen_n_bins_own_case(i_gene)    = merge(chosen_n_bins_case, &
                                                          -chosen_n_bins_case, strata_criteria_met_case)
                chosen_n_bins_own_control(i_gene) = merge(chosen_n_bins_control, &
                                                          -chosen_n_bins_control, strata_criteria_met_control)

                if (case_stratum_count >= 10 .and. control_stratum_count >= 10) then
                    ! Build the `own` null by BOOTSTRAPPING the mean difference: resample
                    ! n_replicates residuals per side (with replacement) from the strata,
                    ! average, difference. n_replicates is each side's per-gene replicate
                    ! count (sorted_*%max_resid_per_gene) — the count the observed means
                    ! were averaged over. Observed statistic is left as-is.
                    call compute_pvalue_bootstrap_mean_helper( &
                        tmp_own_stratum_pool_case(1:case_stratum_count), case_stratum_count, &
                        tmp_own_stratum_pool_control(1:control_stratum_count), control_stratum_count, &
                        sorted_case%max_resid_per_gene, sorted_control%max_resid_per_gene, &
                        abs(observed_statistic_own_val), N_BOOTSTRAP_DRAWS, &
                        pvalues_own(i_gene))
                    neighborhood_size_own_case(i_gene) = case_stratum_count
                    neighborhood_size_own_control(i_gene) = control_stratum_count
                end if
            end if

            neighborhood_size_case(i_gene) = n_pool_case
            ! Count genes that received an own p-value (the stratum gate can still
            ! skip the `own` computation, and compute_pvalue_own can be 0).
            if (pvalues_own(i_gene) >= 0.0_real64) &
                n_genes_with_pvalue = n_genes_with_pvalue + 1
        end do
    end subroutine compute_noise_pvalue_pipeline_helper

    !> Validate inputs, build sorted structures, and run the per-gene pipeline.
    !|
    !| This is the alloc-layer entry point for the full noise-model pipeline. It owns
    !| every allocation needed by the computation (the per-gene work arrays); the
    !| `own` bootstrap draws from the global RNG in place, so no draw buffer is allocated.
    !| Internally it:
    !|   1. Validates all dimension and range arguments via `tox_errors`.
    !|   2. Calls `prepare_sorted_data` for both case and control data.
    !|   3. Allocates all per-gene work arrays used by the per-gene loop.
    !|   4. Delegates the per-gene loop to `compute_noise_pvalue_pipeline_helper`.
    subroutine compute_noise_pvalue_pipeline( &
        means_case, replicates_case, n_genes_case, n_replicates_case, &
        means_control, replicates_control, n_genes_control, n_replicates_control, &
        observed_statistic_own, compute_pvalue_own, &
        n_genes, norm_method, k_start, k_step, k_max, tau, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        chosen_n_bins_own_case, chosen_n_bins_own_control, &
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
        integer(int32), intent(in) :: max_pool_size
        !! Maximum number of residuals in any pool
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: gene-vs-own p-values (-1 if not computed)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes for which the own p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Case stratum size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Control stratum size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for each gene (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_case
        !! Diagnostic: sign-encoded chosen bin count, gene-vs-own CASE stratification
        !! (+ = criteria met, - = coarse 2-bin fallback; -1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_control
        !! Diagnostic: sign-encoded chosen bin count, gene-vs-own CONTROL stratification
        !! (+ = criteria met, - = coarse 2-bin fallback; -1 if not computed)
        integer(int32), intent(out) :: ierr
        !! Error code

        type(sorted_data_t) :: sorted_case, sorted_control
        real(real64), allocatable :: tmp_pool_case(:), tmp_pool_control_own(:)
        integer(int32), allocatable :: tmp_gene_id_pool_case(:), tmp_gene_id_pool_control_own(:)
        ! Case stratification work arrays
        integer(int32), allocatable :: tmp_strat_bin_index_per_residual_case(:)
        integer(int32), allocatable :: tmp_chosen_bin_index_per_residual_case(:)
        real(real64), allocatable :: tmp_strat_bin_edges_case(:), tmp_chosen_bin_edges_case(:)
        integer(int32), allocatable :: tmp_strat_gene_min_bin_case(:), tmp_strat_gene_max_bin_case(:)
        logical, allocatable :: tmp_strat_gene_seen_case(:)
        integer(int32), allocatable :: tmp_strat_touched_gene_slots_case(:)
        integer(int32), allocatable :: tmp_strat_c_g_case(:)
        integer(int32), allocatable :: tmp_strat_bin_counts_case(:)
        ! Control stratification work arrays
        integer(int32), allocatable :: tmp_strat_bin_index_per_residual_control(:)
        integer(int32), allocatable :: tmp_chosen_bin_index_per_residual_control(:)
        real(real64), allocatable :: tmp_strat_bin_edges_control(:), tmp_chosen_bin_edges_control(:)
        integer(int32), allocatable :: tmp_strat_gene_min_bin_control(:), tmp_strat_gene_max_bin_control(:)
        logical, allocatable :: tmp_strat_gene_seen_control(:)
        integer(int32), allocatable :: tmp_strat_touched_gene_slots_control(:)
        integer(int32), allocatable :: tmp_strat_c_g_control(:)
        integer(int32), allocatable :: tmp_strat_bin_counts_control(:)
        real(real64), allocatable :: tmp_own_stratum_pool_case(:), tmp_own_stratum_pool_control(:)
        integer(int32) :: sort_ierr

        call set_ok(ierr)

        ! Seed the RNG (fixed state 42) at the very beginning so the `own`
        ! mean-difference bootstrap is reproducible.
        call init_random(42_int32)

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
        M_ALLOCATE(tmp_gene_id_pool_case(max_pool_size * 2))
        M_ALLOCATE(tmp_pool_control_own(max_pool_size * 2))
        M_ALLOCATE(tmp_gene_id_pool_control_own(max_pool_size * 2))
        
        ! Case stratification work arrays
        M_ALLOCATE(tmp_strat_bin_index_per_residual_case(max_pool_size))
        M_ALLOCATE(tmp_chosen_bin_index_per_residual_case(max_pool_size))
        M_ALLOCATE(tmp_strat_bin_edges_case(STRATA_BIN_COUNT_SCHEDULE(1) + 1))
        M_ALLOCATE(tmp_chosen_bin_edges_case(STRATA_BIN_COUNT_SCHEDULE(1) + 1))
        M_ALLOCATE(tmp_strat_gene_min_bin_case(sorted_case%n_genes))
        M_ALLOCATE(tmp_strat_gene_max_bin_case(sorted_case%n_genes))
        M_ALLOCATE(tmp_strat_gene_seen_case(sorted_case%n_genes))
        M_ALLOCATE(tmp_strat_touched_gene_slots_case(max_pool_size))
        M_ALLOCATE(tmp_strat_c_g_case(max_pool_size))
        M_ALLOCATE(tmp_strat_bin_counts_case(STRATA_BIN_COUNT_SCHEDULE(1)))

        ! Control stratification work arrays
        M_ALLOCATE(tmp_strat_bin_index_per_residual_control(max_pool_size))
        M_ALLOCATE(tmp_chosen_bin_index_per_residual_control(max_pool_size))
        M_ALLOCATE(tmp_strat_bin_edges_control(STRATA_BIN_COUNT_SCHEDULE(1) + 1))
        M_ALLOCATE(tmp_chosen_bin_edges_control(STRATA_BIN_COUNT_SCHEDULE(1) + 1))
        M_ALLOCATE(tmp_strat_gene_min_bin_control(sorted_control%n_genes))
        M_ALLOCATE(tmp_strat_gene_max_bin_control(sorted_control%n_genes))
        M_ALLOCATE(tmp_strat_gene_seen_control(sorted_control%n_genes))
        M_ALLOCATE(tmp_strat_touched_gene_slots_control(max_pool_size))
        M_ALLOCATE(tmp_strat_c_g_control(max_pool_size))
        M_ALLOCATE(tmp_strat_bin_counts_control(STRATA_BIN_COUNT_SCHEDULE(1)))

        M_ALLOCATE(tmp_own_stratum_pool_case(max_pool_size))
        M_ALLOCATE(tmp_own_stratum_pool_control(max_pool_size))

        if (is_err(ierr)) return

        ! `check_stratification_accepted_helper` requires an all-`.false.` "seen"
        ! array on entry and restores it on return, so it only needs zeroing once here.
        tmp_strat_gene_seen_case = .false.
        tmp_strat_gene_seen_control = .false.

        call compute_noise_pvalue_pipeline_helper( &
            sorted_case, sorted_control, &
            means_case, means_control, &
            observed_statistic_own, compute_pvalue_own, &
            n_genes, k_start, k_step, k_max, tau, &
            pvalues_own, n_genes_with_pvalue, &
            max_pool_size, &
            neighborhood_size_own_case, neighborhood_size_own_control, &
            neighborhood_size_case, &
            chosen_n_bins_own_case, chosen_n_bins_own_control, &
            tmp_pool_case, tmp_gene_id_pool_case, &
            tmp_pool_control_own, tmp_gene_id_pool_control_own, &
            tmp_strat_bin_index_per_residual_case, tmp_chosen_bin_index_per_residual_case, &
            tmp_strat_bin_edges_case, tmp_chosen_bin_edges_case, &
            tmp_strat_gene_min_bin_case, tmp_strat_gene_max_bin_case, tmp_strat_gene_seen_case, &
            tmp_strat_touched_gene_slots_case, tmp_strat_c_g_case, tmp_strat_bin_counts_case, &
            tmp_strat_bin_index_per_residual_control, tmp_chosen_bin_index_per_residual_control, &
            tmp_strat_bin_edges_control, tmp_chosen_bin_edges_control, &
            tmp_strat_gene_min_bin_control, tmp_strat_gene_max_bin_control, tmp_strat_gene_seen_control, &
            tmp_strat_touched_gene_slots_control, tmp_strat_c_g_control, tmp_strat_bin_counts_control, &
            tmp_own_stratum_pool_case, tmp_own_stratum_pool_control, &
            ierr)

    end subroutine compute_noise_pvalue_pipeline

end module noise_model

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

#include "macros.h"

!> C-interoperable wrapper for `compute_noise_pvalue_pipeline` (bootstrap-null module).
!|
!| Identical ABI to `compute_noise_pvalues_pipeline_exact_c` (the scaling-null
!| module) — same arguments in the same order — but bound as
!| `compute_noise_pvalues_pipeline_c` and dispatching to `noise_model`, so either
!| can be swapped in without downstream changes.
!|
!| Performs null-pointer checks via `M_CHECK_IERR_NON_NULL` / `M_CHECK_NON_NULL`,
!| then delegates unconditionally to the validated Fortran entry point.
!| No computation is performed here.
subroutine compute_noise_pvalues_pipeline_c( &
    means_case, replicates_case, n_genes_case, n_replicates_case, &
    means_control, replicates_control, n_genes_control, n_replicates_control, &
    observed_statistic_own, compute_pvalue_own, &
    n_genes, norm_method, k_start, k_step, k_max, tau, &
    pvalues_own, n_genes_with_pvalue, &
    max_pool_size, &
    neighborhood_size_own_case, neighborhood_size_own_control, &
    neighborhood_size_case, &
    chosen_n_bins_own_case, chosen_n_bins_own_control, &
    ierr) bind(C, name="compute_noise_pvalues_pipeline_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use noise_model, only: compute_noise_pvalue_pipeline
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
    integer(c_int), intent(in), target :: max_pool_size
    !! Maximum number of residuals in any pool
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_own
    !! Output: gene-vs-own p-values (-1 if not computed)
    integer(c_int), intent(out), target :: n_genes_with_pvalue
    !! Number of genes for which at least one p-value was computed
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_case
    !! Case stratum size used for gene-vs-own (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_control
    !! Control stratum size used for gene-vs-own (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_case
    !! Case pool size used for each gene (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: chosen_n_bins_own_case
    !! Diagnostic: sign-encoded chosen bin count, gene-vs-own CASE stratification
    integer(c_int), dimension(n_genes), intent(out), target :: chosen_n_bins_own_control
    !! Diagnostic: sign-encoded chosen bin count, gene-vs-own CONTROL stratification
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
    M_CHECK_NON_NULL(max_pool_size)
    M_CHECK_NON_NULL(pvalues_own)
    M_CHECK_NON_NULL(n_genes_with_pvalue)
    M_CHECK_NON_NULL(neighborhood_size_own_case)
    M_CHECK_NON_NULL(neighborhood_size_own_control)
    M_CHECK_NON_NULL(neighborhood_size_case)
    M_CHECK_NON_NULL(chosen_n_bins_own_case)
    M_CHECK_NON_NULL(chosen_n_bins_own_control)

    call compute_noise_pvalue_pipeline( &
        means_case, replicates_case, n_genes_case, n_replicates_case, &
        means_control, replicates_control, n_genes_control, n_replicates_control, &
        observed_statistic_own, compute_pvalue_own, &
        n_genes, norm_method, k_start, k_step, k_max, tau, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        chosen_n_bins_own_case, chosen_n_bins_own_control, &
        ierr)

end subroutine compute_noise_pvalues_pipeline_c