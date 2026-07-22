#include "macros.h"

!!! TO DO
!!! 1. Implement error code to argument number mapping
!!! 2. Find solution for family based analysis

!> LFCseq-style noise model comparing case and control gene expression.
!|
!| This is the third `own`-null variant, alongside `noise_model` (bootstraps the
!| BETWEEN-neighbourhood case-vs-control mean difference) and `noise_model_exact`
!| (exact pairwise count + `1/sqrt(n_replicates)` scaling). All three share the same
!| C ABI so any one can be swapped in without downstream changes.
!|
!| The distinguishing feature here is how the `own` null is built (see
!| `compute_pvalue_lfcseq_helper`), following the LFCseq idea of estimating the
!| null from WITHIN-condition comparisons rather than across conditions:
!|   - NO variance stratification: the full gathered case and control neighbourhood
!|     pools are used directly.
!|   - For each draw, two independent group means (of sizes `n_rep_case` and
!|     `n_rep_control`, matching the observed statistic) are drawn from the SAME
!|     pool and differenced — once within the case pool, once within the control
!|     pool. Neither comparison crosses conditions, so under H0 it carries no
!|     signal, only noise.
!|   - The case-side and control-side null statistics are UNIONED into a single
!|     null distribution, against which the observed `|mean_case - mean_control|`
!|     is tested.
!| Group means are drawn at the `n_rep` level (not the individual-residual level),
!| so the null already lives on the same `sigma^2 / n_rep` mean-difference scale as
!| the observed statistic — no post-hoc scaling is needed.
!|
!| The `family` / `ortholog` comparisons (retained for ABI compatibility) still use
!| the exact sorted-pool pairwise count (`compute_pvalue`). The within-neighbourhood
!| `own` null is the one place this module uses the RNG.
!|
!| Provides routines to:
!|   - Sort and pack gene expression residuals for efficient neighbourhood lookup
!|   - Adaptively gather residual pools from the nearest-mean neighbourhood
!|   - Build the `own` null by unioned within-neighbourhood mean comparisons (LFCseq)
!|   - Compute exact family/ortholog p-values via sorted-pool binary search
!|   - Run the full pipeline over all genes with pre-cached family / ortholog pools
module noise_model_lfcseq
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_NAN_INF, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int
    use f42_utils, only: sort_real, init_random
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

    ! No variance stratification in this LFCseq variant: the `own` null is built
    ! from the full gathered neighbourhood pools, so none of the STRATA_* schedule
    ! constants used by `noise_model` are needed here.

    integer(int32), parameter :: N_BOOTSTRAP_DRAWS = 10000_int32
        !! Number of within-neighbourhood resamples drawn PER SIDE. The case-side and
        !! control-side draws are unioned, so the `own` null has `2 * N_BOOTSTRAP_DRAWS`
        !! entries (see `compute_pvalue_lfcseq_helper`).

    real(real64), parameter :: NOISE_LOG_OFFSET = 1.0_real64
        !! Additive constant `c` in the log-space residual
        !! `log2(r + c) - mean_i[log2(r_i + c)]` (see `prepare_sorted_data_helper`),
        !! used when `norm_method /= 0` to keep the logarithm defined at zero expression

    !> Pre-computed residual pools for each gene family and its ortholog set.
    type :: family_cache_t
        real(real64), allocatable :: family_pools(:, :)
        !! Residual pools for families: shape (max_pool_size, n_families)
        real(real64), allocatable :: orth_pools(:, :)
        !! Residual pools for ortholog sets: shape (max_pool_size, n_families)
        integer(int32), allocatable :: family_pool_sizes(:)
        !! Number of valid entries in each family pool
        integer(int32), allocatable :: orth_pool_sizes(:)
        !! Number of valid entries in each ortholog pool
        logical, allocatable :: is_cached(:)
        !! `.true.` if the pool for this family has been computed
    end type family_cache_t

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
        real(real64) :: gene_mean, log2_mean, log2_factor
        logical :: use_log_transform

        sorted_data%n_genes = n_genes
        sorted_data%max_resid_per_gene = n_samples

        use_log_transform = (norm_method /= 0)
        if (use_log_transform) log2_factor = 1.0_real64 / log(2.0_real64)

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
                       use_log_transform, i_gene, orig_idx)
                if (use_log_transform) then
                    sorted_data%residuals_packed(i_sample, i_gene) = &
                        log(max(replicates(i_sample, orig_idx), 0.0_real64) + NOISE_LOG_OFFSET) &
                        * log2_factor - log2_mean
                else
                    sorted_data%residuals_packed(i_sample, i_gene) = &
                        replicates(i_sample, orig_idx) - gene_mean
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
    !| pooled residual was taken from. This model does not use that bookkeeping (it does
    !| not stratify); the argument is kept only so the gather ABI matches the other two
    !| noise-model variants, and the caller may ignore it.
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

    !> Validate inputs and gather an adaptive residual pool for `target_mean`.
    !|
    !| Delegates all computation to `gather_residuals_helper` after checking
    !| that dimension parameters and `tau` are valid.
    pure subroutine gather_residuals(target_mean, sorted_data, &
                                     k_start, k_step, k_max, tau, &
                                     pooled_residuals, gene_id_per_residual, n_pooled, &
                                     max_pool_size, ierr)
        real(real64), intent(in) :: target_mean
        !! Mean value for which a matching residual neighbourhood is sought
        type(sorted_data_t), intent(in) :: sorted_data
        !! Pre-built sorted gene structure
        integer(int32), intent(in) :: k_start
        !! Minimum pool size before the adaptive stopping criterion is applied
        integer(int32), intent(in) :: k_step
        !! Number of new residuals added per adaptive round before re-evaluating
        integer(int32), intent(in) :: k_max
        !! Hard upper limit on pool size (also capped at `max_pool_size`)
        real(real64), intent(in) :: tau
        !! Relative-change threshold; expansion stops when the change exceeds this value
        real(real64), intent(out) :: pooled_residuals(:)
        !! Output residual pool (pre-allocated to at least `max_pool_size`)
        integer(int32), intent(out) :: gene_id_per_residual(:)
        !! Output: sorted-gene-slot that each entry of `pooled_residuals` was taken from
        integer(int32), intent(out) :: n_pooled
        !! Number of residuals written into `pooled_residuals`
        integer(int32), intent(in) :: max_pool_size
        !! Allocated size of `pooled_residuals`
        integer(int32), intent(inout) :: ierr
        !! Error code

        call set_ok(ierr)
        call validate_dimension_size(k_start, ierr)
        call validate_dimension_size(k_step, ierr)
        call validate_dimension_size(k_max, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        if (is_err(ierr)) return

        call gather_residuals_helper(target_mean, sorted_data, &
                                     k_start, k_step, k_max, tau, &
                                     pooled_residuals, gene_id_per_residual, n_pooled, &
                                     max_pool_size)
    end subroutine gather_residuals

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

    !> Draw one bootstrap group mean of `n_draw` residuals (with replacement) from a pool.
    !|
    !| Fills `rbuf(1:n_draw)` with a single `random_number` array call, maps each
    !| `[0,1)` value to an index in `[1, n_pool]` (same mapping as `noise_model`), and
    !| returns the mean of the sampled residuals. `rbuf` must be at least `n_draw`
    !| long. Impure (draws from the global RNG); the pipeline seeds it once via
    !| `init_random(42)`. Drawing each group in its own call, in a fixed statement
    !| order, keeps the whole `own` null reproducible.
    subroutine draw_group_mean_helper(pool, n_pool, n_draw, rbuf, group_mean)
        integer(int32), intent(in) :: n_pool
        !! Size of the source residual pool
        real(real64), dimension(n_pool), intent(in) :: pool
        !! Source residual pool
        integer(int32), intent(in) :: n_draw
        !! Number of residuals to draw (with replacement) into the group
        real(real64), dimension(:), intent(inout) :: rbuf
        !! Scratch buffer for the random draws (length >= n_draw)
        real(real64), intent(out) :: group_mean
        !! Mean of the `n_draw` sampled residuals

        integer(int32) :: i, idx
        real(real64) :: s

        call random_number(rbuf(1:n_draw))
        s = 0.0_real64
        do i = 1, n_draw
            idx = min(int(rbuf(i) * real(n_pool, real64), int32) + 1, n_pool)
            s = s + pool(idx)
        end do
        group_mean = s / real(n_draw, real64)
    end subroutine draw_group_mean_helper

    !> Build the `own` null by unioned WITHIN-neighbourhood mean comparisons (LFCseq).
    !|
    !| The observed statistic is a difference of MEANS, `|mean_case - mean_control|`,
    !| where the two group means average over `n_rep_case` and `n_rep_control`
    !| replicates respectively. This routine estimates the null of that statistic the
    !| LFCseq way — from comparisons that never cross the case/control boundary, so
    !| that under H0 they contain only noise:
    !|
    !|   For each of `n_boot` draws:
    !|     * CASE side:    draw two independent group means (sizes `n_rep_case` and
    !|                     `n_rep_control`) BOTH from `pool_case`; null = |diff|.
    !|     * CONTROL side: draw two independent group means (same sizes) BOTH from
    !|                     `pool_control`; null = |diff|.
    !|
    !| The case-side and control-side null statistics are UNIONED into one null of
    !| `2 * n_boot` entries, and the p-value is the add-one-corrected fraction of that
    !| union which meets or exceeds `observed_statistic_abs`.
    !|
    !| Because the two groups in each comparison have the same sizes as the observed
    !| means, the null already lives on the correct `sigma^2 / n_rep` mean-difference
    !| scale — no `sqrt(n)` correction is applied. No variance stratification is used;
    !| `pool_case` / `pool_control` are the full gathered neighbourhood pools. Impure
    !| (draws from the global RNG), so the pipeline must have seeded it up front.
    subroutine compute_pvalue_lfcseq_helper(pool_case, n_pool_case, &
                                            pool_control, n_pool_control, &
                                            n_rep_case, n_rep_control, &
                                            observed_statistic_abs, n_boot, &
                                            p_value)
        integer(int32), intent(in) :: n_pool_case
        !! Size of the case neighbourhood pool (the case-side sampling source)
        real(real64), dimension(n_pool_case), intent(in) :: pool_case
        !! Case neighbourhood residual pool
        integer(int32), intent(in) :: n_pool_control
        !! Size of the control neighbourhood pool (the control-side sampling source)
        real(real64), dimension(n_pool_control), intent(in) :: pool_control
        !! Control neighbourhood residual pool
        integer(int32), intent(in) :: n_rep_case
        !! Size of the first group in each within-neighbourhood comparison
        !! (= case replicate count the observed case mean averages over)
        integer(int32), intent(in) :: n_rep_control
        !! Size of the second group in each within-neighbourhood comparison
        !! (= control replicate count the observed control mean averages over)
        real(real64), intent(in) :: observed_statistic_abs
        !! Absolute value of the observed test statistic
        integer(int32), intent(in) :: n_boot
        !! Number of within-neighbourhood draws PER SIDE (union has 2*n_boot entries)
        real(real64), intent(out) :: p_value
        !! LFCseq within-neighbourhood p-value

        integer(int32) :: i_boot, count_ge
        real(real64) :: mean_a, mean_b, null_case, null_control
        real(real64) :: rbuf(max(n_rep_case, n_rep_control))

        count_ge = 0
        do i_boot = 1, n_boot
            ! WITHIN the case neighbourhood: two same-side group means differenced.
            call draw_group_mean_helper(pool_case, n_pool_case, n_rep_case, rbuf, mean_a)
            call draw_group_mean_helper(pool_case, n_pool_case, n_rep_control, rbuf, mean_b)
            null_case = abs(mean_a - mean_b)
            if (null_case >= observed_statistic_abs) count_ge = count_ge + 1

            ! WITHIN the control neighbourhood: two same-side group means differenced.
            call draw_group_mean_helper(pool_control, n_pool_control, n_rep_case, rbuf, mean_a)
            call draw_group_mean_helper(pool_control, n_pool_control, n_rep_control, rbuf, mean_b)
            null_control = abs(mean_a - mean_b)
            if (null_control >= observed_statistic_abs) count_ge = count_ge + 1
        end do

        ! Union of the case-side and control-side null sets => 2*n_boot draws total.
        p_value = real(count_ge + 1, real64) / real(2 * n_boot + 1, real64)
    end subroutine compute_pvalue_lfcseq_helper

    ! =========================================================================
    ! =========================================================================
    ! compute_noise_pvalue_pipeline
    ! =========================================================================

    !> Core pipeline: compute per-gene noise p-values using pre-built data structures.
    !|
    !| Iterates over all genes and computes up to three p-values each:
    !|   - `pvalues_own`:  gene vs. its own matched control neighbourhood, using the
    !|     LFCseq within-neighbourhood null (see `compute_pvalue_lfcseq_helper`) built
    !|     from the FULL gathered case and control pools (NO variance stratification)
    !|   - `pvalues_family`:  gene vs. its gene-family pool in control (exact kNN count)
    !|   - `pvalues_ortholog`: gene vs. its ortholog pool in control (exact kNN count)
    !|
    !| For the `own` comparison there is no stratification: both the case and control
    !| kNN neighbourhoods are used whole. The null is formed by unioned within-case and
    !| within-control mean comparisons at the observed `n_rep_case` / `n_rep_control`
    !| group sizes, so it already sits on the mean-difference scale of the observed
    !| statistic. The `family` and `ortholog` comparisons retain the exact whole-pool
    !| kNN pairwise count.
    !|
    !| Requires `sorted_case` and `sorted_control` to already be built (via
    !| `prepare_sorted_data`) and `cache` to be pre-populated (all `is_cached`
    !| entries set for families with `family_sizes > 0`).
    !|
    !| No allocation: the per-gene residual pools are pre-allocated by the caller and
    !| passed in as `intent(inout)` work arrays. The `own` within-neighbourhood null is
    !| drawn from the global RNG (via `compute_pvalue_lfcseq_helper`), which is the only
    !| impure operation, so this subroutine is not `pure`. The gene-id pools are filled
    !| by `gather_residuals_helper` but unused here (this model does not stratify); they
    !| are kept only so the gather ABI is unchanged.
    subroutine compute_noise_pvalue_pipeline_helper( &
        sorted_case, sorted_control, &
        means_case, means_control, &
        observed_statistic_own, observed_statistic_family, observed_statistic_ortholog, &
        family_means, ortholog_means, &
        compute_pvalue_own, compute_pvalue_family, compute_pvalue_ortholog, &
        family_sizes, gene_to_family, &
        n_genes, n_families, k_start, k_step, k_max, tau, &
        pvalues_own, pvalues_family, pvalues_ortholog, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_family, &
        neighborhood_size_ortholog, neighborhood_size_case, &
        chosen_n_bins_own_case, chosen_n_bins_own_control, &
        cache, &
        tmp_pool_case, tmp_gene_id_pool_case, &
        tmp_pool_control_own, tmp_gene_id_pool_control_own, &
        ierr)

        type(sorted_data_t), intent(in) :: sorted_case
        !! Sorted case gene data
        type(sorted_data_t), intent(in) :: sorted_control
        !! Sorted control gene data
        integer(int32), intent(in) :: n_genes
        !! Total number of genes
        integer(int32), intent(in) :: n_families
        !! Total number of gene families
        real(real64), dimension(n_genes), intent(in) :: means_case
        !! Per-gene case expression means
        real(real64), dimension(n_genes), intent(in) :: means_control
        !! Per-gene control expression means
        real(real64), dimension(n_genes), intent(in) :: observed_statistic_own
        !! Observed gene-vs-own statistic for each gene
        real(real64), dimension(n_genes), intent(in) :: observed_statistic_family
        !! Observed gene-vs-family statistic for each gene
        real(real64), dimension(n_genes), intent(in) :: observed_statistic_ortholog
        !! Observed gene-vs-ortholog statistic for each gene
        real(real64), dimension(n_families), intent(in) :: family_means
        !! Mean expression of each gene family (control)
        real(real64), dimension(n_families), intent(in) :: ortholog_means
        !! Mean expression of each ortholog set (control)
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_own
        !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_family
        !! 1 if the gene-vs-family p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_ortholog
        !! 1 if the gene-vs-ortholog p-value should be computed, 0 otherwise
        integer(int32), dimension(n_families), intent(in) :: family_sizes
        !! Number of genes in each family
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
        !! Maps each gene index to its family index
        integer(int32), intent(in) :: k_start
        !! Minimum number of neighbour genes before adaptive stopping is applied
        integer(int32), intent(in) :: k_step
        !! Neighbour genes added per adaptive round before re-evaluating
        integer(int32), intent(in) :: k_max
        !! Hard upper limit on the number of neighbour genes
        real(real64), intent(in) :: tau
        !! Relative-change threshold for adaptive pool growth
        integer(int32), intent(in) :: max_pool_size
        !! Allocated size of all pool arrays
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: gene-vs-own p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: pvalues_family
        !! Output: gene-vs-family p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: pvalues_ortholog
        !! Output: gene-vs-ortholog p-values (-1 if not computed)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes for which at least one p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Case pool size used for the gene-vs-own null (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Control pool size used for the gene-vs-own null (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_family
        !! Control pool size used for gene-vs-family p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_ortholog
        !! Control pool size used for gene-vs-ortholog p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for this gene (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_case
        !! Diagnostic (kept for ABI parity with the stratified models): this model does
        !! not stratify, so it is set to 1 (= whole pool, one bin) when the `own` p-value
        !! is computed, and -1 otherwise.
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_control
        !! Diagnostic (kept for ABI parity): 1 when the `own` p-value is computed
        !! (no stratification, whole pool), -1 otherwise.
        type(family_cache_t), intent(in) :: cache
        !! Pre-computed family and ortholog residual pools
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool_case
        !! Work array: residual pool for this gene's case kNN neighbourhood
        integer(int32), dimension(max_pool_size * 2), intent(inout) :: tmp_gene_id_pool_case
        !! Work array: gene-id pool for the case neighbourhood (filled by gather, unused here)
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool_control_own
        !! Work array: residual pool for this gene's control kNN neighbourhood (the `own` comparison)
        integer(int32), dimension(max_pool_size * 2), intent(inout) :: tmp_gene_id_pool_control_own
        !! Work array: gene-id pool for the control neighbourhood (filled by gather, unused here)
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_gene, family_id
        real(real64) :: mean_case_val, mean_control_val
        real(real64) :: observed_statistic_own_val, observed_statistic_family_val, observed_statistic_ortholog_val
        integer(int32) :: n_pool_case, n_pool_control_own

        call set_ok(ierr)

        pvalues_own = -1.0_real64
        pvalues_family = -1.0_real64
        pvalues_ortholog = -1.0_real64
        neighborhood_size_own_case = -1
        neighborhood_size_own_control = -1
        neighborhood_size_family = -1
        neighborhood_size_ortholog = -1
        neighborhood_size_case = -1
        chosen_n_bins_own_case = -1
        chosen_n_bins_own_control = -1
        n_genes_with_pvalue = 0

        do i_gene = 1, n_genes
            mean_case_val = means_case(i_gene)
            mean_control_val = means_control(i_gene)
            family_id = gene_to_family(i_gene)

            if (family_id < 1 .or. family_id > n_families) cycle
            if (.not. cache%is_cached(family_id)) cycle

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
            observed_statistic_family_val = observed_statistic_family(i_gene)
            observed_statistic_ortholog_val = observed_statistic_ortholog(i_gene)

            ! Skip genes with non-finite observed statistics
            if (observed_statistic_own_val /= observed_statistic_own_val .or. &
                observed_statistic_family_val /= observed_statistic_family_val .or. &
                observed_statistic_ortholog_val /= observed_statistic_ortholog_val) cycle

            if (compute_pvalue_own(i_gene) == 1) then
                ! LFCseq within-neighbourhood null on the FULL gathered pools (no
                ! stratification): union of within-case and within-control mean
                ! comparisons at the observed replicate group sizes. Observed
                ! statistic is left as-is.
                call compute_pvalue_lfcseq_helper( &
                    tmp_pool_case(1:n_pool_case), n_pool_case, &
                    tmp_pool_control_own(1:n_pool_control_own), n_pool_control_own, &
                    sorted_case%max_resid_per_gene, sorted_control%max_resid_per_gene, &
                    abs(observed_statistic_own_val), N_BOOTSTRAP_DRAWS, &
                    pvalues_own(i_gene))
                neighborhood_size_own_case(i_gene) = n_pool_case
                neighborhood_size_own_control(i_gene) = n_pool_control_own
                ! No stratification: report 1 bin (= whole pool) on both sides.
                chosen_n_bins_own_case(i_gene) = 1
                chosen_n_bins_own_control(i_gene) = 1
            end if

            if (compute_pvalue_family(i_gene) == 1) then
                call compute_pvalue(tmp_pool_case(1:n_pool_case), n_pool_case, &
                                    cache%family_pools(1:cache%family_pool_sizes(family_id), family_id), &
                                    cache%family_pool_sizes(family_id), &
                                    observed_statistic_family_val, &
                                    pvalues_family(i_gene), ierr)
                if (is_err(ierr)) return
                neighborhood_size_family(i_gene) = cache%family_pool_sizes(family_id)
            end if

            if (compute_pvalue_ortholog(i_gene) == 1) then
                call compute_pvalue(tmp_pool_case(1:n_pool_case), n_pool_case, &
                                    cache%orth_pools(1:cache%orth_pool_sizes(family_id), family_id), &
                                    cache%orth_pool_sizes(family_id), &
                                    observed_statistic_ortholog_val, &
                                    pvalues_ortholog(i_gene), ierr)
                if (is_err(ierr)) return
                neighborhood_size_ortholog(i_gene) = cache%orth_pool_sizes(family_id)
            end if

            neighborhood_size_case(i_gene) = n_pool_case
            n_genes_with_pvalue = n_genes_with_pvalue + 1
        end do
    end subroutine compute_noise_pvalue_pipeline_helper

    !> Validate inputs, build sorted structures, pre-cache family pools, and run the pipeline.
    !|
    !| This is the alloc-layer entry point for the full noise-model pipeline. It owns
    !| every allocation needed by the computation (the per-gene residual pools); the
    !| `own` within-neighbourhood null draws from the global RNG in place, so no draw
    !| buffer is allocated. Internally it:
    !|   1. Validates all dimension and range arguments via `tox_errors`.
    !|   2. Calls `prepare_sorted_data` for both case and control data.
    !|   3. Pre-computes family and ortholog residual pools via `gather_residuals_helper`.
    !|   4. Allocates the per-gene residual pools used by the per-gene loop.
    !|   5. Delegates the per-gene loop to `compute_noise_pvalue_pipeline_helper`.
    subroutine compute_noise_pvalue_pipeline( &
        means_case, replicates_case, n_genes_case, n_replicates_case, &
        means_control, replicates_control, n_genes_control, n_replicates_control, &
        observed_statistic_own, observed_statistic_family, observed_statistic_ortholog, &
        family_means, ortholog_means, &
        compute_pvalue_own, compute_pvalue_family, compute_pvalue_ortholog, &
        family_sizes, gene_to_family, &
        n_genes, n_families, norm_method, k_start, k_step, k_max, tau, &
        pvalues_own, pvalues_family, pvalues_ortholog, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_family, &
        neighborhood_size_ortholog, neighborhood_size_case, &
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
        integer(int32), intent(in) :: n_families
        !! Total number of gene families
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
        real(real64), dimension(n_genes), intent(in) :: observed_statistic_family
        !! Observed gene-vs-family statistic for each gene
        real(real64), dimension(n_genes), intent(in) :: observed_statistic_ortholog
        !! Observed gene-vs-ortholog statistic for each gene
        real(real64), dimension(n_families), intent(in) :: family_means
        !! Mean expression of each gene family (control)
        real(real64), dimension(n_families), intent(in) :: ortholog_means
        !! Mean expression of each ortholog set (control)
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_own
        !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_family
        !! 1 if the gene-vs-family p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: compute_pvalue_ortholog
        !! 1 if the gene-vs-ortholog p-value should be computed, 0 otherwise
        integer(int32), dimension(n_families), intent(in) :: family_sizes
        !! Number of genes in each family
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
        !! Maps each gene index to its family index (1-based)
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        integer(int32), intent(in) :: k_start
        !! Minimum number of neighbour genes before adaptive stopping is applied
        integer(int32), intent(in) :: k_step
        !! Neighbour genes added per adaptive round before re-evaluating
        integer(int32), intent(in) :: k_max
        !! Hard upper limit on the number of neighbour genes
        real(real64), intent(in) :: tau
        !! Relative-change threshold for adaptive pool growth
        integer(int32), intent(in) :: max_pool_size
        !! Maximum number of residuals in any pool
        real(real64), dimension(n_genes), intent(out) :: pvalues_own
        !! Output: gene-vs-own p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: pvalues_family
        !! Output: gene-vs-family p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: pvalues_ortholog
        !! Output: gene-vs-ortholog p-values (-1 if not computed)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes for which at least one p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Case pool size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Control pool size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_family
        !! Control pool size used for gene-vs-family (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_ortholog
        !! Control pool size used for gene-vs-ortholog (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for each gene (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_case
        !! Diagnostic (ABI parity): 1 when `own` computed (no stratification), else -1
        integer(int32), dimension(n_genes), intent(out) :: chosen_n_bins_own_control
        !! Diagnostic (ABI parity): 1 when `own` computed (no stratification), else -1
        integer(int32), intent(out) :: ierr
        !! Error code

        type(sorted_data_t) :: sorted_case, sorted_control
        type(family_cache_t) :: cache
        real(real64), allocatable :: tmp_pool_case(:), tmp_pool_control_own(:)
        integer(int32), allocatable :: tmp_gene_id_pool_case(:), tmp_gene_id_pool_control_own(:)
        integer(int32), allocatable :: family_gene_id_pool(:), orth_gene_id_pool(:)
        integer(int32) :: family_id, sort_ierr

        call set_ok(ierr)

        ! Seed the RNG (fixed state 42) at the very beginning so the `own`
        ! within-neighbourhood null is reproducible.
        call init_random(42_int32)

        call validate_dimension_size(n_genes_case, ierr)
        call validate_dimension_size(n_replicates_case, ierr)
        call validate_dimension_size(n_genes_control, ierr)
        call validate_dimension_size(n_replicates_control, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_dimension_size(k_start, ierr)
        call validate_dimension_size(k_step, ierr)
        call validate_dimension_size(k_max, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        call validate_all_in_range_real(means_case, n_genes_case, ierr)
        call validate_all_in_range_real(means_control, n_genes_control, ierr)
        call validate_all_in_range_real(replicates_case, n_replicates_case * n_genes_case, ierr)
        call validate_all_in_range_real(replicates_control, n_replicates_control * n_genes_control, ierr)
        call validate_all_in_range_int(gene_to_family, n_genes, ierr, min=1, max=n_families)
        if (is_err(ierr)) return

        call prepare_sorted_data(means_case, replicates_case, &
                                 n_replicates_case, n_genes_case, norm_method, sorted_case, sort_ierr)
        call set_err(ierr, sort_ierr)

        call prepare_sorted_data(means_control, replicates_control, &
                                 n_replicates_control, n_genes_control, norm_method, sorted_control, sort_ierr)
        call set_err(ierr, sort_ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(cache%family_pools(max_pool_size, n_families))
        M_ALLOCATE(cache%orth_pools(max_pool_size, n_families))
        M_ALLOCATE(cache%family_pool_sizes(n_families))
        M_ALLOCATE(cache%orth_pool_sizes(n_families))
        M_ALLOCATE(cache%is_cached(n_families))
        M_ALLOCATE(family_gene_id_pool(max_pool_size))
        M_ALLOCATE(orth_gene_id_pool(max_pool_size))

        cache%is_cached = .false.

        do family_id = 1, n_families
            if (family_sizes(family_id) > 0) then
                call gather_residuals_helper(family_means(family_id), sorted_control, &
                                             k_start, k_step, k_max, tau, &
                                             cache%family_pools(:, family_id), family_gene_id_pool, &
                                             cache%family_pool_sizes(family_id), &
                                             max_pool_size)
                call gather_residuals_helper(ortholog_means(family_id), sorted_control, &
                                             k_start, k_step, k_max, tau, &
                                             cache%orth_pools(:, family_id), orth_gene_id_pool, &
                                             cache%orth_pool_sizes(family_id), &
                                             max_pool_size)
                cache%is_cached(family_id) = .true.
            end if
        end do

        ! Per-gene work arrays for compute_noise_pvalue_pipeline_helper, allocated
        ! once here so the helper itself performs no allocation. No stratification
        ! scratch arrays are needed in this model.
        M_ALLOCATE(tmp_pool_case(max_pool_size * 2))
        M_ALLOCATE(tmp_gene_id_pool_case(max_pool_size * 2))
        M_ALLOCATE(tmp_pool_control_own(max_pool_size * 2))
        M_ALLOCATE(tmp_gene_id_pool_control_own(max_pool_size * 2))

        if (is_err(ierr)) return

        call compute_noise_pvalue_pipeline_helper( &
            sorted_case, sorted_control, &
            means_case, means_control, &
            observed_statistic_own, observed_statistic_family, observed_statistic_ortholog, &
            family_means, ortholog_means, &
            compute_pvalue_own, compute_pvalue_family, compute_pvalue_ortholog, &
            family_sizes, gene_to_family, &
            n_genes, n_families, k_start, k_step, k_max, tau, &
            pvalues_own, pvalues_family, pvalues_ortholog, n_genes_with_pvalue, &
            max_pool_size, &
            neighborhood_size_own_case, neighborhood_size_own_control, &
            neighborhood_size_family, &
            neighborhood_size_ortholog, neighborhood_size_case, &
            chosen_n_bins_own_case, chosen_n_bins_own_control, &
            cache, &
            tmp_pool_case, tmp_gene_id_pool_case, &
            tmp_pool_control_own, tmp_gene_id_pool_control_own, &
            ierr)

    end subroutine compute_noise_pvalue_pipeline

end module noise_model_lfcseq

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

#include "macros.h"

!> C-interoperable wrapper for `compute_noise_pvalue_pipeline` (LFCseq-null module).
!|
!| Identical ABI to `compute_noise_pvalues_pipeline_c` (bootstrap null) and
!| `compute_noise_pvalues_pipeline_exact_c` (scaling null) — same arguments in the
!| same order — but bound as `compute_noise_pvalues_pipeline_lfcseq_c` and
!| dispatching to `noise_model_lfcseq`, so any of the three can be swapped in
!| without downstream changes. The `chosen_n_bins_own_*` outputs are retained for
!| ABI parity but carry no stratification meaning here (see the pipeline helper).
!|
!| Performs null-pointer checks via `M_CHECK_IERR_NON_NULL` / `M_CHECK_NON_NULL`,
!| then delegates unconditionally to the validated Fortran entry point.
!| No computation is performed here.
subroutine compute_noise_pvalues_pipeline_lfcseq_c( &
    means_case, replicates_case, n_genes_case, n_replicates_case, &
    means_control, replicates_control, n_genes_control, n_replicates_control, &
    observed_statistic_own, observed_statistic_family, observed_statistic_ortholog, &
    family_means, ortholog_means, &
    compute_pvalue_own, compute_pvalue_family, compute_pvalue_ortholog, &
    family_sizes, gene_to_family, &
    n_genes, n_families, norm_method, k_start, k_step, k_max, tau, &
    pvalues_own, pvalues_family, pvalues_ortholog, n_genes_with_pvalue, &
    max_pool_size, &
    neighborhood_size_own_case, neighborhood_size_own_control, &
    neighborhood_size_family, &
    neighborhood_size_ortholog, neighborhood_size_case, &
    chosen_n_bins_own_case, chosen_n_bins_own_control, &
    ierr) bind(C, name="compute_noise_pvalues_pipeline_lfcseq_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use noise_model_lfcseq, only: compute_noise_pvalue_pipeline
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
    integer(c_int), intent(in), target :: n_families
    !! Total number of gene families
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
    real(c_double), dimension(n_genes), intent(in), target :: observed_statistic_family
    !! Observed gene-vs-family statistic for each gene
    real(c_double), dimension(n_genes), intent(in), target :: observed_statistic_ortholog
    !! Observed gene-vs-ortholog statistic for each gene
    real(c_double), dimension(n_families), intent(in), target :: family_means
    !! Mean expression of each gene family (control)
    real(c_double), dimension(n_families), intent(in), target :: ortholog_means
    !! Mean expression of each ortholog set (control)
    integer(c_int), dimension(n_genes), intent(in), target :: compute_pvalue_own
    !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
    integer(c_int), dimension(n_genes), intent(in), target :: compute_pvalue_family
    !! 1 if the gene-vs-family p-value should be computed, 0 otherwise
    integer(c_int), dimension(n_genes), intent(in), target :: compute_pvalue_ortholog
    !! 1 if the gene-vs-ortholog p-value should be computed, 0 otherwise
    integer(c_int), dimension(n_families), intent(in), target :: family_sizes
    !! Number of genes in each family
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_family
    !! Maps each gene index to its family index (1-based)
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
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_family
    !! Output: gene-vs-family p-values (-1 if not computed)
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_ortholog
    !! Output: gene-vs-ortholog p-values (-1 if not computed)
    integer(c_int), intent(out), target :: n_genes_with_pvalue
    !! Number of genes for which at least one p-value was computed
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_case
    !! Case stratum size used for gene-vs-own (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own_control
    !! Control stratum size used for gene-vs-own (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_family
    !! Control pool size used for gene-vs-family (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_ortholog
    !! Control pool size used for gene-vs-ortholog (-1 if not computed)
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
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(means_case)
    M_CHECK_NON_NULL(replicates_case)
    M_CHECK_NON_NULL(means_control)
    M_CHECK_NON_NULL(replicates_control)
    M_CHECK_NON_NULL(observed_statistic_own)
    M_CHECK_NON_NULL(observed_statistic_family)
    M_CHECK_NON_NULL(observed_statistic_ortholog)
    M_CHECK_NON_NULL(family_means)
    M_CHECK_NON_NULL(ortholog_means)
    M_CHECK_NON_NULL(compute_pvalue_own)
    M_CHECK_NON_NULL(compute_pvalue_family)
    M_CHECK_NON_NULL(compute_pvalue_ortholog)
    M_CHECK_NON_NULL(family_sizes)
    M_CHECK_NON_NULL(gene_to_family)
    M_CHECK_NON_NULL(norm_method)
    M_CHECK_NON_NULL(k_start)
    M_CHECK_NON_NULL(k_step)
    M_CHECK_NON_NULL(k_max)
    M_CHECK_NON_NULL(tau)
    M_CHECK_NON_NULL(max_pool_size)
    M_CHECK_NON_NULL(pvalues_own)
    M_CHECK_NON_NULL(pvalues_family)
    M_CHECK_NON_NULL(pvalues_ortholog)
    M_CHECK_NON_NULL(n_genes_with_pvalue)
    M_CHECK_NON_NULL(neighborhood_size_own_case)
    M_CHECK_NON_NULL(neighborhood_size_own_control)
    M_CHECK_NON_NULL(neighborhood_size_family)
    M_CHECK_NON_NULL(neighborhood_size_ortholog)
    M_CHECK_NON_NULL(neighborhood_size_case)
    M_CHECK_NON_NULL(chosen_n_bins_own_case)
    M_CHECK_NON_NULL(chosen_n_bins_own_control)

    call compute_noise_pvalue_pipeline( &
        means_case, replicates_case, n_genes_case, n_replicates_case, &
        means_control, replicates_control, n_genes_control, n_replicates_control, &
        observed_statistic_own, observed_statistic_family, observed_statistic_ortholog, &
        family_means, ortholog_means, &
        compute_pvalue_own, compute_pvalue_family, compute_pvalue_ortholog, &
        family_sizes, gene_to_family, &
        n_genes, n_families, norm_method, k_start, k_step, k_max, tau, &
        pvalues_own, pvalues_family, pvalues_ortholog, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_family, &
        neighborhood_size_ortholog, neighborhood_size_case, &
        chosen_n_bins_own_case, chosen_n_bins_own_control, &
        ierr)

end subroutine compute_noise_pvalues_pipeline_lfcseq_c