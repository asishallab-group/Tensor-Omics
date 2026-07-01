#include "macros.h"

!> Noise model for exact p-value computation comparing case and control gene expression.
!|
!| Provides routines to:
!|   - Sort and pack gene expression residuals for efficient neighbourhood lookup
!|   - Adaptively gather residual pools from the nearest-mean neighbourhood
!|   - Variance-stratify the gathered pool into quantile intervals so that the
!|     `own` comparison draws only from the gene's own variance regime, rather
!|     than from a mixture of regimes across the kNN neighbourhood
!|   - Compute exact p-values by exhaustively enumerating every pairwise
!|     residual difference between the case and control residual pools (fully
!|     deterministic, no Monte Carlo sampling)
!|   - Run the full pipeline over all genes with pre-cached family / ortholog pools
module noise_model
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_NAN_INF, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int
    use f42_utils, only: sort_real, sort_integer, calc_percentile_helper
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

    integer(int32), parameter :: STRATA_N_SCHEDULE_STEPS = 5
        !! Number of candidate bin-count steps tried by `stratify_residuals_helper`
    integer(int32), parameter :: STRATA_BIN_COUNT_SCHEDULE(STRATA_N_SCHEDULE_STEPS) = &
        [20_int32, 10_int32, 5_int32, 4_int32, 2_int32]
        !! Candidate quantile bin counts, finest (5% per bin) to coarsest (50% per
        !! bin) first, tried in order until a step satisfies all acceptance
        !! criteria. The last entry (2 bins, 50% per bin) is a hard floor: it is
        !! accepted unconditionally if no earlier step succeeds.
    real(real64), parameter :: STRATA_MAX_C_G_PROB_THRESHOLD = 0.1_real64
        !! Criterion 2: `Pr(c_g > 2)` must stay below this fraction
    integer(int32), parameter :: STRATA_MIN_RESIDUALS_PER_BIN = 50_int32
        !! Criterion 3: every occupied quantile interval must contain at least this
        !! many residuals to support a stable exact null distribution
    
    !> Threshold for switching from exact enumeration to Monte Carlo sampling
    integer(int64), parameter :: MAX_EXACT_COMBINATIONS = 10000_int32
        !! Maximum number of pairwise combinations for exhaustive enumeration
    integer(int32), parameter :: MONTE_CARLO_SAMPLES = 10000_int32
        !! Number of samples for Monte Carlo p-value estimation

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
    !| Fills `sorted_data` in-place. Does not allocate — all arrays inside
    !| `sorted_data` must already be allocated by the caller to the correct sizes,
    !| and `tmp_perm`, `tmp_stack_left`, `tmp_stack_right` must be pre-allocated
    !| work arrays of length `n_genes`.
    pure subroutine prepare_sorted_data_helper(means, replicates, n_samples, n_genes, &
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
        type(sorted_data_t), intent(inout) :: sorted_data
        !! Sorted data structure to fill; all allocatable fields must be pre-allocated
        integer(int32), dimension(n_genes), intent(inout) :: tmp_perm
        !! Work array: permutation vector for indirect sort (length n_genes)
        integer(int32), dimension(n_genes), intent(inout) :: tmp_stack_left
        !! Work array: quicksort left-index stack (length n_genes)
        integer(int32), dimension(n_genes), intent(inout) :: tmp_stack_right
        !! Work array: quicksort right-index stack (length n_genes)

        integer(int32) :: i_gene, i_sample, orig_idx
        real(real64) :: gene_mean

        sorted_data%n_genes = n_genes
        sorted_data%max_resid_per_gene = n_samples

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
            do concurrent(i_sample=1:n_samples) shared(sorted_data, replicates, gene_mean, i_gene, orig_idx)
                sorted_data%residuals_packed(i_sample, i_gene) = &
                    replicates(i_sample, orig_idx) - gene_mean
            end do
        end do
    end subroutine prepare_sorted_data_helper

    !> Validate inputs, allocate the `sorted_data` structure, and sort genes by mean.
    !|
    !| This is the alloc-layer entry point. It allocates all fields of `sorted_data`
    !| and the internal work arrays, then delegates to `prepare_sorted_data_helper`.
    subroutine prepare_sorted_data(means, replicates, n_samples, n_genes, sorted_data, ierr)
        integer(int32), intent(in) :: n_samples
        !! Number of replicates per gene
        integer(int32), intent(in) :: n_genes
        !! Number of genes
        real(real64), dimension(n_genes), intent(in) :: means
        !! Per-gene expression means (length n_genes)
        real(real64), dimension(n_samples, n_genes), intent(in) :: replicates
        !! Replicate expression matrix (n_samples x n_genes)
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

        call prepare_sorted_data_helper(means, replicates, n_samples, n_genes, &
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
    !| Starting from the gene whose mean is closest to `target_mean`, neighbours are
    !| added outward (alternating left/right in mean-sorted order) until one of:
    !|   - The pool contains at least `k_start` residuals (initial phase), then
    !|   - The relative change in mean absolute residual exceeds `tau` (adaptive phase), or
    !|   - The pool reaches `min(k_max, max_pool_size)` residuals.
    !|
    !| Alongside the pooled residuals, `gene_id_per_residual` records the sorted-gene-slot
    !| (index into `sorted_data%means_sorted` / `sorted_data%original_indices`) that each
    !| pooled residual was taken from. This bookkeeping is required by variance
    !| stratification (see `stratify_residuals_helper`), which needs to know how many
    !| distinct genes contribute residuals to a given quantile interval.
    !|
    !| No allocation; `pooled_residuals`, `gene_id_per_residual`, and `tmp_work` /
    !| `tmp_gene_id_work` must be pre-allocated by caller.
    pure subroutine gather_residuals_helper(target_mean, sorted_data, &
                                            k_start, k_step, k_max, tau, &
                                            pooled_residuals, gene_id_per_residual, n_pooled, &
                                            max_pool_size, &
                                            tmp_work, tmp_gene_id_work)
        real(real64), intent(in) :: target_mean
        !! Mean value for which a matching residual neighbourhood is sought
        type(sorted_data_t), intent(in) :: sorted_data
        !! Pre-built sorted gene structure
        integer(int32), intent(in) :: k_start
        !! Minimum pool size before the adaptive stopping criterion is applied
        integer(int32), intent(in) :: k_step
        !! Number of new residuals to add per adaptive round before re-evaluating
        integer(int32), intent(in) :: k_max
        !! Hard upper limit on pool size (also capped at `max_pool_size`)
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
        !! Allocated size of `pooled_residuals`
        real(real64), intent(inout) :: tmp_work(:)
        !! Work array for staging candidate expansions (pre-allocated to `max_pool_size`)
        integer(int32), intent(inout) :: tmp_gene_id_work(:)
        !! Work array for staging candidate gene-id expansions (pre-allocated to `max_pool_size`)

        integer(int32) :: pos, left_cand, right_cand, idx
        integer(int32) :: current_size, added_this_round, genes_added, offset
        integer(int32) :: pool_size, n_resid
        real(real64) :: S_old, S_new, rel_change

        n_pooled = 0
        current_size = 0

        pos = find_closest_helper(target_mean, sorted_data%means_sorted, sorted_data%n_genes)
        if (pos == 0) return

        current_size = min(sorted_data%n_residuals(pos), max_pool_size)
        pooled_residuals(1:current_size) = sorted_data%residuals_packed(1:current_size, pos)
        gene_id_per_residual(1:current_size) = pos

        left_cand = pos - 1
        right_cand = pos + 1

        ! Phase 1: expand until pool reaches k_start to ensure minimum size regardless of tau
        do while (current_size < k_start .and. &
                  (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))

            ! choose closest candidate
            call choose_index(sorted_data%means_sorted, sorted_data%n_genes, target_mean, idx, left_cand, right_cand)

            ! Add residuals to pool
            pool_size = min(current_size + sorted_data%n_residuals(idx), max_pool_size)
            call add_residuals_to_pool_helper(pooled_residuals, current_size, &
                                              sorted_data%residuals_packed(:, idx), &
                                              sorted_data%n_residuals(idx), pool_size)
            call add_gene_id_to_pool_helper(gene_id_per_residual, current_size, idx, pool_size)
            current_size = pool_size
        end do

        ! If max found residuals are below 10 return (possibly error or status here?)
        if (current_size < 10) return

        ! Compute current variance measure of the residuals
        S_old = sum(abs(pooled_residuals(1:current_size))) / real(current_size, real64)
        if (S_old == 0.0_real64) then
            n_pooled = current_size
            return
        end if

        ! Phase 2: adaptive expansion 
        do while (current_size < min(k_max, max_pool_size) .and. &
                  (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))

            added_this_round = 0
            genes_added = 0
            offset = 0

            tmp_work(1:current_size) = pooled_residuals(1:current_size)
            tmp_gene_id_work(1:current_size) = gene_id_per_residual(1:current_size)

            do while (added_this_round < k_step .and. &
                      (left_cand >= 1 .or. right_cand <= sorted_data%n_genes) .and. &
                      current_size + offset < min(k_max, max_pool_size))

                call choose_index(sorted_data%means_sorted, sorted_data%n_genes, target_mean, idx, left_cand, right_cand)

                n_resid = sorted_data%n_residuals(idx)
                pool_size = min(current_size + offset + n_resid, max_pool_size)
                call add_residuals_to_pool_helper(tmp_work, current_size + offset, &
                                                  sorted_data%residuals_packed(:, idx), &
                                                  n_resid, pool_size)
                call add_gene_id_to_pool_helper(tmp_gene_id_work, current_size + offset, idx, pool_size)
                offset = offset + n_resid
                added_this_round = added_this_round + n_resid
                genes_added = genes_added + 1
            end do

            if (genes_added == 0) exit

            S_new = sum(abs(tmp_work(1:current_size + offset))) / &
                    real(current_size + offset, real64)
            rel_change = (S_new - S_old) / S_old

            if (rel_change > tau) exit

            current_size = current_size + offset
            pooled_residuals(1:current_size) = tmp_work(1:current_size)
            gene_id_per_residual(1:current_size) = tmp_gene_id_work(1:current_size)
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
                                     max_pool_size, &
                                     tmp_work, tmp_gene_id_work, ierr)
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
        real(real64), intent(inout) :: tmp_work(:)
        !! Work array for staging candidate expansions (pre-allocated to `max_pool_size`)
        integer(int32), intent(inout) :: tmp_gene_id_work(:)
        !! Work array for staging candidate gene-id expansions (pre-allocated to `max_pool_size`)
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
                                     max_pool_size, &
                                     tmp_work, tmp_gene_id_work)
    end subroutine gather_residuals

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
    !| Uses local automatic work arrays (`perm`, `stack_left`, `stack_right`) sized
    !| `n_pooled` for the one-off sort needed to compute percentiles; this keeps the
    !| routine pure without requiring the caller to pre-allocate yet another pair of
    !| `tmp_*` arrays for what is a self-contained, bounded-size sort.
    pure subroutine assign_residual_bins_helper(pooled_residuals, n_pooled, n_bins, &
                                                bin_index_per_residual, bin_edges)
        integer(int32), intent(in) :: n_pooled
        !! Number of valid entries in `pooled_residuals`
        real(real64), dimension(n_pooled), intent(in) :: pooled_residuals
        !! Pooled residual values
        integer(int32), intent(in) :: n_bins
        !! Number of quantile bins to divide the residuals into
        integer(int32), dimension(n_pooled), intent(out) :: bin_index_per_residual
        !! Output: 1-based bin index for each residual
        real(real64), dimension(n_bins + 1), intent(out) :: bin_edges
        !! Output: the `n_bins + 1` quantile cutpoints, from the 0th to the 100th percentile

        integer(int32) :: i_resid, i_edge, perm(n_pooled), stack_left(n_pooled), stack_right(n_pooled)
        real(real64) :: percentile_step

        do concurrent(i_resid=1:n_pooled) shared(perm)
            perm(i_resid) = i_resid
        end do
        call sort_real(pooled_residuals, perm, stack_left, stack_right)

        percentile_step = 100.0_real64 / real(n_bins, real64)
        do i_edge = 1, n_bins + 1
            call calc_percentile_helper(pooled_residuals, perm, &
                                        percentile_step * real(i_edge - 1, real64), bin_edges(i_edge))
        end do

        if (bin_edges(n_bins + 1) <= bin_edges(1)) then
            bin_index_per_residual = 1
            return
        end if

        do concurrent(i_resid=1:n_pooled) shared(bin_index_per_residual, pooled_residuals, bin_edges, n_bins)
            bin_index_per_residual(i_resid) = locate_bin_helper(pooled_residuals(i_resid), bin_edges, n_bins)
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
    !| `tmp_c_g`, and `tmp_bin_counts` are work arrays pre-allocated by the caller.
    pure subroutine check_stratification_accepted_helper( &
        bin_index_per_residual, gene_id_per_residual, n_pooled, n_bins, n_genes_total, &
        tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, tmp_c_g, tmp_bin_counts, &
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
        !! Work array: `.true.` once a gene slot has contributed at least one residual
        integer(int32), dimension(n_pooled), intent(inout) :: tmp_c_g
        !! Work array: per-distinct-gene `c_g` values (only the first `n_distinct_genes` entries are used)
        integer(int32), dimension(n_bins), intent(inout) :: tmp_bin_counts
        !! Work array: number of residuals occupying each bin
        logical, intent(out) :: is_accepted
        !! `.true.` if all three acceptance criteria are satisfied

        integer(int32) :: i_resid, gene_id, bin_id, n_distinct_genes, i_gene_slot
        integer(int32) :: n_c_g_gt_2, median_c_g, n_occupied_bins, min_occupied_count

        tmp_gene_seen = .false.
        tmp_bin_counts = 0

        do i_resid = 1, n_pooled
            gene_id = gene_id_per_residual(i_resid)
            bin_id = bin_index_per_residual(i_resid)

            tmp_bin_counts(bin_id) = tmp_bin_counts(bin_id) + 1

            if (.not. tmp_gene_seen(gene_id)) then
                tmp_gene_seen(gene_id) = .true.
                tmp_gene_min_bin(gene_id) = bin_id
                tmp_gene_max_bin(gene_id) = bin_id
            else
                tmp_gene_min_bin(gene_id) = min(tmp_gene_min_bin(gene_id), bin_id)
                tmp_gene_max_bin(gene_id) = max(tmp_gene_max_bin(gene_id), bin_id)
            end if
        end do

        n_distinct_genes = 0
        do i_gene_slot = 1, n_genes_total
            if (tmp_gene_seen(i_gene_slot)) then
                n_distinct_genes = n_distinct_genes + 1
                tmp_c_g(n_distinct_genes) = tmp_gene_max_bin(i_gene_slot) - tmp_gene_min_bin(i_gene_slot) + 1
            end if
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
    !| largest schedule step (`max_bins = STRATA_BIN_COUNT_SCHEDULE(1)`).
    pure subroutine stratify_residuals_helper( &
        pooled_residuals, gene_id_per_residual, n_pooled, n_genes_total, &
        tmp_bin_index_per_residual, tmp_bin_edges, &
        tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, tmp_c_g, tmp_bin_counts, &
        chosen_n_bins, chosen_bin_index_per_residual, chosen_bin_edges, criteria_met)
        integer(int32), intent(in) :: n_pooled
        !! Number of pooled residuals to stratify
        integer(int32), intent(in) :: n_genes_total
        !! Total number of genes in the underlying sorted structure
        real(real64), dimension(n_pooled), intent(in) :: pooled_residuals
        !! Pooled residual values
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

        integer(int32) :: i_step, n_bins
        logical :: is_accepted

        criteria_met = .false.

        do i_step = 1, STRATA_N_SCHEDULE_STEPS
            n_bins = STRATA_BIN_COUNT_SCHEDULE(i_step)

            call assign_residual_bins_helper(pooled_residuals, n_pooled, n_bins, &
                                             tmp_bin_index_per_residual, tmp_bin_edges(1:n_bins + 1))

            call check_stratification_accepted_helper( &
                tmp_bin_index_per_residual, gene_id_per_residual, n_pooled, n_bins, n_genes_total, &
                tmp_gene_min_bin, tmp_gene_max_bin, tmp_gene_seen, tmp_c_g, tmp_bin_counts(1:n_bins), &
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

    ! =========================================================================
    ! compute_pvalue
    ! =========================================================================

    !> Core implementation: exact p-value via exhaustive pairwise residual enumeration.
    !|
    !| Builds the null distribution deterministically: for every residual
    !| `r_case` in `pool_case` and every residual `r_control` in `pool_control`,
    !| the null statistic is `(mean_case + r_case) - (mean_control + r_control)`
    !| (or its log2-transformed form when `norm_method /= 0`). The p-value is the
    !| fraction of these `n_pool_case * n_pool_control` pairwise null statistics
    !| whose magnitude meets or exceeds `observed_statistic_abs`, with the usual
    !| add-one continuity correction.
    !|
    !| This replaces Monte Carlo resampling: since residual pools are bounded by
    !| `max_pool_size` and therefore never very large, exhaustively enumerating
    !| every pair is both fully deterministic (no RNG, no sampling variance) and
    !| cheap enough to be the more accurate choice here.
    !|
    !| When `norm_method /= 0` a log2 transform is applied:
    !| `|log2(x_case+1) - log2(x_control+1)|`. Otherwise the difference is on the
    !| original scale: `|x_case - x_control|`.
    pure subroutine compute_pvalue_helper(mean_case, pool_case, n_pool_case, &
                                          mean_control, pool_control, n_pool_control, &
                                          observed_statistic_abs, norm_method, &
                                          p_value)
        real(real64), intent(in) :: mean_case
        !! Case mean for this gene
        integer(int32), intent(in) :: n_pool_case
        !! Size of the case residual pool
        real(real64), dimension(n_pool_case), intent(in) :: pool_case
        !! Case residual pool
        real(real64), intent(in) :: mean_control
        !! Control mean for this gene
        integer(int32), intent(in) :: n_pool_control
        !! Size of the control residual pool
        real(real64), dimension(n_pool_control), intent(in) :: pool_control
        !! Control residual pool
        real(real64), intent(in) :: observed_statistic_abs
        !! Absolute value of the observed test statistic
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        real(real64), intent(out) :: p_value
        !! Exact p-value

        integer(int32) :: i_case, i_control, count_ge
        real(real64) :: log2_factor, x_case, x_control, null_dist
        logical :: use_log_transform
        integer(int64) :: n_pairs

        use_log_transform = (norm_method /= 0)
        if (use_log_transform) log2_factor = 1.0_real64 / log(2.0_real64)

        count_ge = 0
        do concurrent(i_case=1:n_pool_case, i_control=1:n_pool_control) &
            shared(pool_case, pool_control, mean_case, mean_control, use_log_transform, log2_factor, &
                   observed_statistic_abs) &
            reduce(+:count_ge) &
            local(x_case, x_control, null_dist)

            x_case = pool_case(i_case)
            x_control = pool_control(i_control)

            if (use_log_transform) then
                ! Currently NOT correct, causes errors!
                x_case = max(x_case, 0.0_real64) + 1.0_real64
                x_control = max(x_control, 0.0_real64) + 1.0_real64
                null_dist = abs((log(x_case) - log(x_control)) * log2_factor)
            else
                null_dist = abs(x_case - x_control)
            end if

            if (null_dist >= observed_statistic_abs) count_ge = count_ge + 1
        end do

        n_pairs = int(n_pool_case, int64) * int(n_pool_control, int64)
        p_value = real(count_ge + 1, real64) / real(n_pairs + 1_int64, real64)
    end subroutine compute_pvalue_helper

    !> Core implementation: Monte Carlo p-value estimation.
    !|
    !| When the number of pairwise combinations exceeds MAX_EXACT_COMBINATIONS,
    !| this routine draws MONTE_CARLO_SAMPLES random pairs from the case and control
    !| pools (with replacement) and estimates the p-value as the fraction of sampled
    !| null statistics that exceed or equal the observed statistic.
    !|
    !| Uses pre-generated random indices stored in `rand_indices_case` and
    !| `rand_indices_control` (each of length MONTE_CARLO_SAMPLES, with values in
    !| 1..n_pool_case and 1..n_pool_control respectively) so that the RNG is
    !| called only once per gene.
    pure subroutine compute_pvalue_monte_carlo_helper( &
        mean_case, pool_case, n_pool_case, &
        mean_control, pool_control, n_pool_control, &
        observed_statistic_abs, norm_method, &
        rand_indices_case, rand_indices_control, &
        p_value)
        real(real64), intent(in) :: mean_case
        !! Case mean for this gene
        integer(int32), intent(in) :: n_pool_case
        !! Size of the case residual pool
        real(real64), dimension(n_pool_case), intent(in) :: pool_case
        !! Case residual pool
        real(real64), intent(in) :: mean_control
        !! Control mean for this gene
        integer(int32), intent(in) :: n_pool_control
        !! Size of the control residual pool
        real(real64), dimension(n_pool_control), intent(in) :: pool_control
        !! Control residual pool
        real(real64), intent(in) :: observed_statistic_abs
        !! Absolute value of the observed test statistic
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        integer(int32), dimension(MONTE_CARLO_SAMPLES), intent(in) :: rand_indices_case
        !! Pre-generated random indices into pool_case
        integer(int32), dimension(MONTE_CARLO_SAMPLES), intent(in) :: rand_indices_control
        !! Pre-generated random indices into pool_control
        real(real64), intent(out) :: p_value
        !! Monte Carlo p-value estimate

        integer(int32) :: i_sample, count_ge
        real(real64) :: log2_factor, x_case, x_control, null_dist
        logical :: use_log_transform

        use_log_transform = (norm_method /= 0)
        if (use_log_transform) log2_factor = 1.0_real64 / log(2.0_real64)

        count_ge = 0
        do i_sample = 1, MONTE_CARLO_SAMPLES
            x_case = pool_case(rand_indices_case(i_sample))
            x_control = pool_control(rand_indices_control(i_sample))

            if (use_log_transform) then
                x_case = max(x_case, 0.0_real64) + 1.0_real64
                x_control = max(x_control, 0.0_real64) + 1.0_real64
                null_dist = abs((log(x_case) - log(x_control)) * log2_factor)
            else
                null_dist = abs(x_case - x_control)
            end if

            if (null_dist >= observed_statistic_abs) count_ge = count_ge + 1
        end do

        p_value = real(count_ge + 1, real64) / real(MONTE_CARLO_SAMPLES + 1, real64)
    end subroutine compute_pvalue_monte_carlo_helper

    !> Validate inputs and compute a p-value (exact or Monte Carlo) based on pool sizes.
    !|
    !| This is the validated entry point: it checks dimensions and residual-pool
    !| values, then delegates to either `compute_pvalue_helper` (exact enumeration)
    !| or `compute_pvalue_monte_carlo_helper` (Monte Carlo sampling) depending on
    !| whether `n_pool_case * n_pool_control <= MAX_EXACT_COMBINATIONS`.
    !|
    !| For Monte Carlo mode, random indices are expected to be pre-generated for
    !| the gene and passed in as `rand_indices_case` and `rand_indices_control`.
    pure subroutine compute_pvalue(mean_case, pool_case, n_pool_case, &
                                   mean_control, pool_control, n_pool_control, &
                                   observed_statistic, norm_method, &
                                   rand_indices_case, rand_indices_control, &
                                   p_value, ierr)
        real(real64), intent(in) :: mean_case
        !! Case mean for this gene
        integer(int32), intent(in) :: n_pool_case
        !! Size of the case residual pool
        real(real64), dimension(n_pool_case), intent(in) :: pool_case
        !! Case residual pool
        real(real64), intent(in) :: mean_control
        !! Control mean for this gene
        integer(int32), intent(in) :: n_pool_control
        !! Size of the control residual pool
        real(real64), dimension(n_pool_control), intent(in) :: pool_control
        !! Control residual pool
        real(real64), intent(in) :: observed_statistic
        !! Observed test statistic (sign is discarded; absolute value is used)
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        integer(int32), dimension(MONTE_CARLO_SAMPLES), intent(in), optional :: rand_indices_case
        !! Pre-generated random indices into pool_case (required for Monte Carlo mode)
        integer(int32), dimension(MONTE_CARLO_SAMPLES), intent(in), optional :: rand_indices_control
        !! Pre-generated random indices into pool_control (required for Monte Carlo mode)
        real(real64), intent(out) :: p_value
        !! p-value (exact or Monte Carlo estimated)
        integer(int32), intent(inout) :: ierr
        !! Error code

        integer(int64) :: n_combinations

        call validate_dimension_size(n_pool_case, ierr)
        call validate_dimension_size(n_pool_control, ierr)
        call validate_all_in_range_real(pool_case, n_pool_case, ierr)
        call validate_all_in_range_real(pool_control, n_pool_control, ierr)
        if (is_err(ierr)) return

        n_combinations = int(n_pool_case, int64) * int(n_pool_control, int64)

        if (n_combinations <= MAX_EXACT_COMBINATIONS) then
            call compute_pvalue_helper(mean_case, pool_case, n_pool_case, &
                                       mean_control, pool_control, n_pool_control, &
                                       abs(observed_statistic), norm_method, &
                                       p_value)
        else
            if (.not. present(rand_indices_case) .or. .not. present(rand_indices_control)) then
                call set_err(ierr, ERR_INVALID_INPUT)
                return
            end if
            call compute_pvalue_monte_carlo_helper(mean_case, pool_case, n_pool_case, &
                                                   mean_control, pool_control, n_pool_control, &
                                                   abs(observed_statistic), norm_method, &
                                                   rand_indices_case, rand_indices_control, &
                                                   p_value)
        end if
    end subroutine compute_pvalue

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
    !| passed in as an `intent(inout)` work array, as required for a pure core
    !| routine. p-values are computed by either exact pairwise residual enumeration
    !| or Monte Carlo sampling (see `compute_pvalue`), so no random-number
    !| machinery is needed except for the pre-generated random indices.
    pure subroutine compute_noise_pvalue_pipeline_helper( &
        sorted_case, sorted_control, &
        means_case, means_control, &
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
        cache, &
        tmp_pool, tmp_gene_id_pool, &
        tmp_pool_case, tmp_gene_id_pool_case, &
        tmp_pool_control_own, tmp_gene_id_pool_control_own, &
        tmp_strat_bin_index_per_residual_case, tmp_chosen_bin_index_per_residual_case, &
        tmp_strat_bin_edges_case, tmp_chosen_bin_edges_case, &
        tmp_strat_gene_min_bin_case, tmp_strat_gene_max_bin_case, &
        tmp_strat_gene_seen_case, tmp_strat_c_g_case, tmp_strat_bin_counts_case, &
        tmp_strat_bin_index_per_residual_control, tmp_chosen_bin_index_per_residual_control, &
        tmp_strat_bin_edges_control, tmp_chosen_bin_edges_control, &
        tmp_strat_gene_min_bin_control, tmp_strat_gene_max_bin_control, &
        tmp_strat_gene_seen_control, tmp_strat_c_g_control, tmp_strat_bin_counts_control, &
        tmp_own_stratum_pool_case, tmp_own_stratum_pool_control, &
        rand_indices_case, rand_indices_control)

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
        !! Stratum size used for the gene-vs-own case pool (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Stratum size used for the gene-vs-own control pool (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_family
        !! Control pool size used for gene-vs-family p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_ortholog
        !! Control pool size used for gene-vs-ortholog p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for this gene (-1 if not computed)
        type(family_cache_t), intent(in) :: cache
        !! Pre-computed family and ortholog residual pools
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool
        !! Work array: generic residual-pool expansion staging for `gather_residuals_helper`,
        !! shared across the case, control, family, and ortholog gather calls
        integer(int32), dimension(max_pool_size * 2), intent(inout) :: tmp_gene_id_pool
        !! Work array: generic gene-id expansion staging for `gather_residuals_helper`,
        !! shared across the case, control, family, and ortholog gather calls
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
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(max_pool_size), intent(inout) :: tmp_strat_c_g_control
        !! Work array, see `check_stratification_accepted_helper`
        integer(int32), dimension(STRATA_BIN_COUNT_SCHEDULE(1)), intent(inout) :: tmp_strat_bin_counts_control
        !! Work array, see `check_stratification_accepted_helper`
        real(real64), dimension(max_pool_size), intent(inout) :: tmp_own_stratum_pool_case
        !! Work array: case residuals restricted to the stratum containing the gene's own case mean
        real(real64), dimension(max_pool_size), intent(inout) :: tmp_own_stratum_pool_control
        !! Work array: control residuals restricted to the stratum containing the gene's own control mean
        integer(int32), dimension(MONTE_CARLO_SAMPLES), intent(inout) :: rand_indices_case
        !! Pre-generated random indices for case pool sampling (Monte Carlo mode)
        integer(int32), dimension(MONTE_CARLO_SAMPLES), intent(inout) :: rand_indices_control
        !! Pre-generated random indices for control pool sampling (Monte Carlo mode)

        integer(int32) :: i_gene, family_id
        real(real64) :: mean_case_val, mean_control_val
        real(real64) :: observed_statistic_own_val, observed_statistic_family_val, observed_statistic_ortholog_val
        integer(int32) :: n_pool_case, n_pool_control_own
        integer(int32) :: case_stratum_count, control_stratum_count
        integer(int32) :: chosen_n_bins_case, chosen_n_bins_control
        logical :: strata_criteria_met_case, strata_criteria_met_control
        integer(int32) :: dummy_ierr

        pvalues_own = -1.0_real64
        pvalues_family = -1.0_real64
        pvalues_ortholog = -1.0_real64
        neighborhood_size_own_case = -1
        neighborhood_size_own_control = -1
        neighborhood_size_family = -1
        neighborhood_size_ortholog = -1
        neighborhood_size_case = -1
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
                                         max_pool_size, &
                                         tmp_pool, tmp_gene_id_pool)
            if (n_pool_case < 10) cycle

            call gather_residuals_helper(mean_control_val, sorted_control, &
                                         k_start, k_step, k_max, tau, &
                                         tmp_pool_control_own, tmp_gene_id_pool_control_own, n_pool_control_own, &
                                         max_pool_size, &
                                         tmp_pool, tmp_gene_id_pool)
            if (n_pool_control_own < 10) cycle

            observed_statistic_own_val = observed_statistic_own(i_gene)
            observed_statistic_family_val = observed_statistic_family(i_gene)
            observed_statistic_ortholog_val = observed_statistic_ortholog(i_gene)

            ! Skip genes with non-finite observed statistics
            if (observed_statistic_own_val /= observed_statistic_own_val .or. &
                observed_statistic_family_val /= observed_statistic_family_val .or. &
                observed_statistic_ortholog_val /= observed_statistic_ortholog_val) cycle

            if (compute_pvalue_own(i_gene) == 1) then
                ! Variance-stratify both case and control `own` neighbourhoods,
                ! then restrict sampling to the stratum that contains this gene's
                ! own mean (case mean for case pool, control mean for control pool).
                
                ! Stratify case pool
                call stratify_residuals_helper( &
                    tmp_pool_case(1:n_pool_case), &
                    tmp_gene_id_pool_case(1:n_pool_case), &
                    n_pool_case, sorted_case%n_genes, &
                    tmp_strat_bin_index_per_residual_case(1:n_pool_case), tmp_strat_bin_edges_case, &
                    tmp_strat_gene_min_bin_case, tmp_strat_gene_max_bin_case, tmp_strat_gene_seen_case, &
                    tmp_strat_c_g_case(1:n_pool_case), tmp_strat_bin_counts_case, &
                    chosen_n_bins_case, tmp_chosen_bin_index_per_residual_case(1:n_pool_case), &
                    tmp_chosen_bin_edges_case, strata_criteria_met_case)

                call select_stratum_for_target_helper( &
                    tmp_pool_case(1:n_pool_case), &
                    tmp_chosen_bin_index_per_residual_case(1:n_pool_case), &
                    n_pool_case, &
                    tmp_chosen_bin_edges_case(1:chosen_n_bins_case + 1), chosen_n_bins_case, mean_case_val, &
                    tmp_own_stratum_pool_case(1:n_pool_case), case_stratum_count)

                ! Stratify control pool
                call stratify_residuals_helper( &
                    tmp_pool_control_own(1:n_pool_control_own), &
                    tmp_gene_id_pool_control_own(1:n_pool_control_own), &
                    n_pool_control_own, sorted_control%n_genes, &
                    tmp_strat_bin_index_per_residual_control(1:n_pool_control_own), tmp_strat_bin_edges_control, &
                    tmp_strat_gene_min_bin_control, tmp_strat_gene_max_bin_control, tmp_strat_gene_seen_control, &
                    tmp_strat_c_g_control(1:n_pool_control_own), tmp_strat_bin_counts_control, &
                    chosen_n_bins_control, tmp_chosen_bin_index_per_residual_control(1:n_pool_control_own), &
                    tmp_chosen_bin_edges_control, strata_criteria_met_control)

                call select_stratum_for_target_helper( &
                    tmp_pool_control_own(1:n_pool_control_own), &
                    tmp_chosen_bin_index_per_residual_control(1:n_pool_control_own), &
                    n_pool_control_own, &
                    tmp_chosen_bin_edges_control(1:chosen_n_bins_control + 1), chosen_n_bins_control, mean_control_val, &
                    tmp_own_stratum_pool_control(1:n_pool_control_own), control_stratum_count)

                ! Use stratified pools if both have at least 10 residuals
                if (case_stratum_count >= 10 .and. control_stratum_count >= 10) then
                    call compute_pvalue(mean_case_val, tmp_own_stratum_pool_case(1:case_stratum_count), case_stratum_count, &
                                        mean_control_val, tmp_own_stratum_pool_control(1:control_stratum_count), &
                                        control_stratum_count, &
                                        observed_statistic_own_val, norm_method, &
                                        rand_indices_case, rand_indices_control, &
                                        pvalues_own(i_gene), dummy_ierr)
                    neighborhood_size_own_case(i_gene) = case_stratum_count
                    neighborhood_size_own_control(i_gene) = control_stratum_count
                end if
            end if

            if (compute_pvalue_family(i_gene) == 1) then
                call compute_pvalue(mean_case_val, tmp_pool_case(1:n_pool_case), n_pool_case, &
                                    family_means(family_id), &
                                    cache%family_pools(1:cache%family_pool_sizes(family_id), family_id), &
                                    cache%family_pool_sizes(family_id), &
                                    observed_statistic_family_val, norm_method, &
                                    rand_indices_case, rand_indices_control, &
                                    pvalues_family(i_gene), dummy_ierr)
                neighborhood_size_family(i_gene) = cache%family_pool_sizes(family_id)
            end if

            if (compute_pvalue_ortholog(i_gene) == 1) then
                call compute_pvalue(mean_case_val, tmp_pool_case(1:n_pool_case), n_pool_case, &
                                    ortholog_means(family_id), &
                                    cache%orth_pools(1:cache%orth_pool_sizes(family_id), family_id), &
                                    cache%orth_pool_sizes(family_id), &
                                    observed_statistic_ortholog_val, norm_method, &
                                    rand_indices_case, rand_indices_control, &
                                    pvalues_ortholog(i_gene), dummy_ierr)
                neighborhood_size_ortholog(i_gene) = cache%orth_pool_sizes(family_id)
            end if

            neighborhood_size_case(i_gene) = n_pool_case
            n_genes_with_pvalue = n_genes_with_pvalue + 1
        end do
    end subroutine compute_noise_pvalue_pipeline_helper

    !> Validate inputs, build sorted structures, pre-cache family pools, and run the pipeline.
    !|
    !| This is the alloc-layer entry point for the full noise-model pipeline. It owns
    !| every allocation needed by the computation, including the per-gene work
    !| arrays, so that `compute_noise_pvalue_pipeline_helper` itself can remain pure
    !| (no allocation). Internally it:
    !|   1. Validates all dimension and range arguments via `tox_errors`.
    !|   2. Calls `prepare_sorted_data` for both case and control data.
    !|   3. Pre-computes family and ortholog residual pools via `gather_residuals_helper`.
    !|   4. Allocates all per-gene work arrays used by the per-gene loop.
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
        real(real64), dimension(n_genes), intent(out) :: pvalues_family
        !! Output: gene-vs-family p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: pvalues_ortholog
        !! Output: gene-vs-ortholog p-values (-1 if not computed)
        integer(int32), intent(out) :: n_genes_with_pvalue
        !! Number of genes for which at least one p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_case
        !! Case stratum size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own_control
        !! Control stratum size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_family
        !! Control pool size used for gene-vs-family (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_ortholog
        !! Control pool size used for gene-vs-ortholog (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_case
        !! Case pool size used for each gene (-1 if not computed)
        integer(int32), intent(out) :: ierr
        !! Error code

        type(sorted_data_t) :: sorted_case, sorted_control
        type(family_cache_t) :: cache
        real(real64), allocatable :: tmp_pool(:), tmp_pool_case(:), tmp_pool_control_own(:)
        integer(int32), allocatable :: tmp_gene_id_pool(:), tmp_gene_id_pool_case(:), tmp_gene_id_pool_control_own(:)
        integer(int32), allocatable :: family_gene_id_pool(:), orth_gene_id_pool(:)
        ! Case stratification work arrays
        integer(int32), allocatable :: tmp_strat_bin_index_per_residual_case(:)
        integer(int32), allocatable :: tmp_chosen_bin_index_per_residual_case(:)
        real(real64), allocatable :: tmp_strat_bin_edges_case(:), tmp_chosen_bin_edges_case(:)
        integer(int32), allocatable :: tmp_strat_gene_min_bin_case(:), tmp_strat_gene_max_bin_case(:)
        logical, allocatable :: tmp_strat_gene_seen_case(:)
        integer(int32), allocatable :: tmp_strat_c_g_case(:)
        integer(int32), allocatable :: tmp_strat_bin_counts_case(:)
        ! Control stratification work arrays
        integer(int32), allocatable :: tmp_strat_bin_index_per_residual_control(:)
        integer(int32), allocatable :: tmp_chosen_bin_index_per_residual_control(:)
        real(real64), allocatable :: tmp_strat_bin_edges_control(:), tmp_chosen_bin_edges_control(:)
        integer(int32), allocatable :: tmp_strat_gene_min_bin_control(:), tmp_strat_gene_max_bin_control(:)
        logical, allocatable :: tmp_strat_gene_seen_control(:)
        integer(int32), allocatable :: tmp_strat_c_g_control(:)
        integer(int32), allocatable :: tmp_strat_bin_counts_control(:)
        real(real64), allocatable :: tmp_own_stratum_pool_case(:), tmp_own_stratum_pool_control(:)
        integer(int32), allocatable :: rand_indices_case(:), rand_indices_control(:)
        integer(int32) :: family_id, sort_ierr, i

        call set_ok(ierr)

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
                                 n_replicates_case, n_genes_case, sorted_case, sort_ierr)
        call set_err(ierr, sort_ierr)

        call prepare_sorted_data(means_control, replicates_control, &
                                 n_replicates_control, n_genes_control, sorted_control, sort_ierr)
        call set_err(ierr, sort_ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(cache%family_pools(max_pool_size, n_families))
        M_ALLOCATE(cache%orth_pools(max_pool_size, n_families))
        M_ALLOCATE(cache%family_pool_sizes(n_families))
        M_ALLOCATE(cache%orth_pool_sizes(n_families))
        M_ALLOCATE(cache%is_cached(n_families))
        M_ALLOCATE(tmp_pool(max_pool_size * 2))
        M_ALLOCATE(tmp_gene_id_pool(max_pool_size * 2))
        M_ALLOCATE(family_gene_id_pool(max_pool_size))
        M_ALLOCATE(orth_gene_id_pool(max_pool_size))

        cache%is_cached = .false.

        do family_id = 1, n_families
            if (family_sizes(family_id) > 0) then
                call gather_residuals_helper(family_means(family_id), sorted_control, &
                                             k_start, k_step, k_max, tau, &
                                             cache%family_pools(:, family_id), family_gene_id_pool, &
                                             cache%family_pool_sizes(family_id), &
                                             max_pool_size, &
                                             tmp_pool, tmp_gene_id_pool)
                call gather_residuals_helper(ortholog_means(family_id), sorted_control, &
                                             k_start, k_step, k_max, tau, &
                                             cache%orth_pools(:, family_id), orth_gene_id_pool, &
                                             cache%orth_pool_sizes(family_id), &
                                             max_pool_size, &
                                             tmp_pool, tmp_gene_id_pool)
                cache%is_cached(family_id) = .true.
            end if
        end do

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
        M_ALLOCATE(tmp_strat_c_g_control(max_pool_size))
        M_ALLOCATE(tmp_strat_bin_counts_control(STRATA_BIN_COUNT_SCHEDULE(1)))
        
        M_ALLOCATE(tmp_own_stratum_pool_case(max_pool_size))
        M_ALLOCATE(tmp_own_stratum_pool_control(max_pool_size))
        
        ! Random indices for Monte Carlo sampling
        M_ALLOCATE(rand_indices_case(MONTE_CARLO_SAMPLES))
        M_ALLOCATE(rand_indices_control(MONTE_CARLO_SAMPLES))
        
        if (is_err(ierr)) return

        call compute_noise_pvalue_pipeline_helper( &
            sorted_case, sorted_control, &
            means_case, means_control, &
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
            cache, &
            tmp_pool, tmp_gene_id_pool, &
            tmp_pool_case, tmp_gene_id_pool_case, &
            tmp_pool_control_own, tmp_gene_id_pool_control_own, &
            tmp_strat_bin_index_per_residual_case, tmp_chosen_bin_index_per_residual_case, &
            tmp_strat_bin_edges_case, tmp_chosen_bin_edges_case, &
            tmp_strat_gene_min_bin_case, tmp_strat_gene_max_bin_case, tmp_strat_gene_seen_case, &
            tmp_strat_c_g_case, tmp_strat_bin_counts_case, &
            tmp_strat_bin_index_per_residual_control, tmp_chosen_bin_index_per_residual_control, &
            tmp_strat_bin_edges_control, tmp_chosen_bin_edges_control, &
            tmp_strat_gene_min_bin_control, tmp_strat_gene_max_bin_control, tmp_strat_gene_seen_control, &
            tmp_strat_c_g_control, tmp_strat_bin_counts_control, &
            tmp_own_stratum_pool_case, tmp_own_stratum_pool_control, &
            rand_indices_case, rand_indices_control)

    end subroutine compute_noise_pvalue_pipeline

end module noise_model

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

#include "macros.h"

!> C-interoperable wrapper for `compute_noise_pvalue_pipeline`.
!|
!| Performs null-pointer checks via `M_CHECK_IERR_NON_NULL` / `M_CHECK_NON_NULL`,
!| then delegates unconditionally to the validated Fortran entry point.
!| No computation is performed here.
pure subroutine compute_noise_pvalues_pipeline_c( &
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
        ierr)

end subroutine compute_noise_pvalues_pipeline_c