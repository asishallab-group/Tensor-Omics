#include "macros.h"

!> Noise model for Monte Carlo p-value computation in cancer vs. healthy gene expression.
!|
!| Provides routines to:
!|   - Sort and pack gene expression residuals for efficient neighbourhood lookup
!|   - Adaptively gather residual pools from the nearest-mean neighbourhood
!|   - Compute Monte Carlo p-values by resampling residuals under the null
!|   - Run the full pipeline over all genes with pre-cached family / ortholog pools
module tox_noise_model
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
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
        !! Maximum residuals per gene
        integer(int32) :: n_genes
        !! Total number of genes in the sorted structure
    end type sorted_data_t

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
                sorted_data%residuals_packed(i_sample, i_gene) = replicates(i_sample, orig_idx) - gene_mean
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
    pure function find_closest_helper(target, means_sorted, n_means) result(pos)
        real(real64), intent(in) :: target
        !! Query value
        integer(int32), intent(in) :: n_means
        !! Length of `means_sorted`
        real(real64), dimension(n_means), intent(in) :: means_sorted
        !! Sorted mean array (ascending)
        integer(int32) :: pos
        !! Index of the element closest to `target`; 0 if `n_means == 0`

        integer(int32) :: left, right, mid

        if (n_means == 0) then
            pos = 0
            return
        end if

        if (target <= means_sorted(1)) then
            pos = 1
            return
        end if

        if (target >= means_sorted(n_means)) then
            pos = n_means
            return
        end if

        left = 1
        right = n_means
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
        else if (left > n_means) then
            pos = n_means
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
    !| No allocation; `pooled_residuals` and `tmp_work` must be pre-allocated by caller.
    pure subroutine gather_residuals_helper(target_mean, sorted_data, &
                                            k_start, k_step, k_max, tau, &
                                            pooled_residuals, n_pooled, max_pool_size, &
                                            tmp_work)
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
        integer(int32), intent(out) :: n_pooled
        !! Number of residuals written into `pooled_residuals`
        integer(int32), intent(in) :: max_pool_size
        !! Allocated size of `pooled_residuals`
        real(real64), intent(inout) :: tmp_work(:)
        !! Work array for staging candidate expansions (pre-allocated to `max_pool_size`)

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

        left_cand = pos - 1
        right_cand = pos + 1

        ! Phase 1: expand until pool reaches k_start
        do while (current_size < k_start .and. &
                  (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))

            if (left_cand >= 1 .and. right_cand <= sorted_data%n_genes) then
                if (abs(sorted_data%means_sorted(left_cand) - target_mean) <= &
                    abs(sorted_data%means_sorted(right_cand) - target_mean)) then
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

            pool_size = min(current_size + sorted_data%n_residuals(idx), max_pool_size)
            call add_residuals_to_pool_helper(pooled_residuals, current_size, &
                                              sorted_data%residuals_packed(:, idx), &
                                              sorted_data%n_residuals(idx), pool_size)
            current_size = pool_size
        end do

        S_old = sum(abs(pooled_residuals(1:current_size))) / real(current_size, real64)
        if (S_old == 0.0_real64) then
            n_pooled = current_size
            return
        end if

        ! ── Phase 2: adaptive expansion ──────────────────────────────────────
        do while (current_size < min(k_max, max_pool_size) .and. &
                  (left_cand >= 1 .or. right_cand <= sorted_data%n_genes))

            added_this_round = 0
            genes_added = 0
            offset = 0

            tmp_work(1:current_size) = pooled_residuals(1:current_size)

            do while (added_this_round < k_step .and. &
                      (left_cand >= 1 .or. right_cand <= sorted_data%n_genes) .and. &
                      current_size + offset < min(k_max, max_pool_size))

                if (left_cand >= 1 .and. right_cand <= sorted_data%n_genes) then
                    if (abs(sorted_data%means_sorted(left_cand) - target_mean) <= &
                        abs(sorted_data%means_sorted(right_cand) - target_mean)) then
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

                n_resid = sorted_data%n_residuals(idx)
                pool_size = min(current_size + offset + n_resid, max_pool_size)
                call add_residuals_to_pool_helper(tmp_work, current_size + offset, &
                                                  sorted_data%residuals_packed(:, idx), &
                                                  n_resid, pool_size)
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
            S_old = S_new
        end do

        n_pooled = current_size
    end subroutine gather_residuals_helper

    !> Validate inputs and gather an adaptive residual pool for `target_mean`.
    !|
    !| Delegates all computation to `gather_residuals_helper` after checking
    !| that dimension parameters and `tau` are valid.
    pure subroutine gather_residuals(target_mean, sorted_data, &
                                     k_start, k_step, k_max, tau, &
                                     pooled_residuals, n_pooled, max_pool_size, &
                                     tmp_work, ierr)
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
        integer(int32), intent(out) :: n_pooled
        !! Number of residuals written into `pooled_residuals`
        integer(int32), intent(in) :: max_pool_size
        !! Allocated size of `pooled_residuals`
        real(real64), intent(inout) :: tmp_work(:)
        !! Work array for staging candidate expansions (pre-allocated to `max_pool_size`)
        integer(int32), intent(inout) :: ierr
        !! Error code

        call validate_dimension_size(k_start, ierr)
        call validate_dimension_size(k_step, ierr)
        call validate_dimension_size(k_max, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        if (is_err(ierr)) return

        call gather_residuals_helper(target_mean, sorted_data, &
                                     k_start, k_step, k_max, tau, &
                                     pooled_residuals, n_pooled, max_pool_size, &
                                     tmp_work)

    end subroutine gather_residuals

    ! =========================================================================
    ! compute_pvalue
    ! =========================================================================

    !> Core implementation: Monte Carlo p-value via residual resampling.
    !|
    !| Samples `r_c` residuals from `resid_c` and `r_h` from `resid_h`, adds them
    !| to `mu_c` and `mu_h` respectively, and counts how often the simulated absolute
    !| difference meets or exceeds `obs_abs`. The raw random matrices are passed in
    !| as pre-generated work arrays to allow batch generation outside the hot loop.
    !|
    !| When `norm_method /= 0` a log2 transform is applied: `|log2(x_c+1) - log2(x_h+1)|`.
    !| Otherwise the difference is on the original scale: `|x_c - x_h|`.
    pure subroutine compute_pvalue_helper(mu_c, r_c, resid_c, n_c, &
                                          mu_h, r_h, resid_h, n_h, &
                                          obs_abs, n_draws, norm_method, &
                                          random_c, random_h, &
                                          p)
        real(real64), intent(in) :: mu_c
        !! Cancer group mean for this gene
        integer(int32), intent(in) :: r_c
        !! Number of cancer replicates to resample per null draw
        integer(int32), intent(in) :: n_c
        !! Size of the cancer residual pool
        real(real64), dimension(n_c), intent(in) :: resid_c
        !! Cancer residual pool
        real(real64), intent(in) :: mu_h
        !! Healthy group mean for this gene
        integer(int32), intent(in) :: r_h
        !! Number of healthy replicates to resample per null draw
        integer(int32), intent(in) :: n_h
        !! Size of the healthy residual pool
        real(real64), dimension(n_h), intent(in) :: resid_h
        !! Healthy residual pool
        real(real64), intent(in) :: obs_abs
        !! Absolute value of the observed test statistic
        integer(int32), intent(in) :: n_draws
        !! Number of Monte Carlo draws
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        real(real64), dimension(r_c, n_draws), intent(in) :: random_c
        !! Pre-generated uniform random numbers for cancer resampling (r_c x n_draws)
        real(real64), dimension(r_h, n_draws), intent(in) :: random_h
        !! Pre-generated uniform random numbers for healthy resampling (r_h x n_draws)
        real(real64), intent(out) :: p
        !! Monte Carlo p-value

        integer(int32) :: i_draw, j, idx_c, idx_h, count_ge
        real(real64) :: log2_factor, eta_c, eta_h, x_c, x_h, null_dist
        logical :: use_log_transform

        use_log_transform = (norm_method /= 0)
        if (use_log_transform) log2_factor = 1.0_real64 / log(2.0_real64)

        count_ge = 0

        do i_draw = 1, n_draws
            eta_c = 0.0_real64
            do j = 1, r_c
                idx_c = min(max(int(random_c(j, i_draw) * real(n_c, real64)) + 1, 1), n_c)
                eta_c = eta_c + resid_c(idx_c)
            end do
            eta_c = eta_c / real(r_c, real64)

            eta_h = 0.0_real64
            do j = 1, r_h
                idx_h = min(max(int(random_h(j, i_draw) * real(n_h, real64)) + 1, 1), n_h)
                eta_h = eta_h + resid_h(idx_h)
            end do
            eta_h = eta_h / real(r_h, real64)

            x_c = mu_c + eta_c
            x_h = mu_h + eta_h

            if (use_log_transform) then
                x_c = max(x_c, 0.0_real64) + 1.0_real64
                x_h = max(x_h, 0.0_real64) + 1.0_real64
                null_dist = abs((log(x_c) - log(x_h)) * log2_factor)
            else
                null_dist = abs(x_c - x_h)
            end if

            if (null_dist >= obs_abs) count_ge = count_ge + 1
        end do

        p = real(count_ge + 1, real64) / real(n_draws + 1, real64)
    end subroutine compute_pvalue_helper

    !> Validate inputs, allocate random matrices, and compute a Monte Carlo p-value.
    !|
    !| Generates two batches of uniform random numbers (`r_c x n_draws` and `r_h x n_draws`),
    !| then delegates to `compute_pvalue_helper`. This is the alloc-layer entry point.
    subroutine compute_pvalue(mu_c, r_c, resid_c, n_c, &
                              mu_h, r_h, resid_h, n_h, &
                              obs, n_draws, norm_method, &
                              p, ierr)
        real(real64), intent(in) :: mu_c
        !! Cancer group mean for this gene
        integer(int32), intent(in) :: r_c
        !! Number of cancer replicates to resample per null draw
        integer(int32), intent(in) :: n_c
        !! Size of the cancer residual pool
        real(real64), dimension(n_c), intent(in) :: resid_c
        !! Cancer residual pool
        real(real64), intent(in) :: mu_h
        !! Healthy group mean for this gene
        integer(int32), intent(in) :: r_h
        !! Number of healthy replicates to resample per null draw
        integer(int32), intent(in) :: n_h
        !! Size of the healthy residual pool
        real(real64), dimension(n_h), intent(in) :: resid_h
        !! Healthy residual pool
        real(real64), intent(in) :: obs
        !! Observed test statistic (sign is discarded; absolute value is used)
        integer(int32), intent(in) :: n_draws
        !! Number of Monte Carlo draws
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        real(real64), intent(out) :: p
        !! Monte Carlo p-value
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64), allocatable :: random_c(:, :), random_h(:, :)

        call set_ok(ierr)

        call validate_dimension_size(r_c, ierr)
        call validate_dimension_size(r_h, ierr)
        call validate_dimension_size(n_c, ierr)
        call validate_dimension_size(n_h, ierr)
        call validate_dimension_size(n_draws, ierr)
        call validate_all_in_range_real(resid_c, n_c, ierr)
        call validate_all_in_range_real(resid_h, n_h, ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(random_c(r_c, n_draws))
        M_ALLOCATE(random_h(r_h, n_draws))

        call random_number(random_c)
        call random_number(random_h)

        call compute_pvalue_helper(mu_c, r_c, resid_c, n_c, &
                                   mu_h, r_h, resid_h, n_h, &
                                   abs(obs), n_draws, norm_method, &
                                   random_c, random_h, &
                                   p)
    end subroutine compute_pvalue

    ! =========================================================================
    ! compute_noise_pvalue_pipeline
    ! =========================================================================

    !> Core pipeline: compute per-gene noise p-values using pre-built data structures.
    !|
    !| Iterates over all genes and computes up to three p-values each:
    !|   - `pvalues_own`:  gene vs. its own matched healthy neighbourhood
    !|   - `pvalues_fam`:  gene vs. its gene-family healthy pool
    !|   - `pvalues_orth`: gene vs. its ortholog healthy pool
    !|
    !| Requires `cancer_sorted` and `healthy_sorted` to already be built (via
    !| `prepare_sorted_data`) and `cache` to be pre-populated (all `is_cached`
    !| entries set for families with `family_sizes > 0`).
    !|
    !| No allocation; all arrays are passed in by the caller.
    subroutine compute_noise_pvalue_pipeline_helper( &
        cancer_sorted, healthy_sorted, &
        cancer_means, healthy_means, &
        cancer_n_samples, healthy_n_samples, &
        obs_own, obs_fam, obs_orth, &
        family_means, ortholog_means, &
        valid_genes_own, valid_genes_fam, valid_genes_orth, &
        family_sizes, is_ortholog_sum, gene_to_family, &
        n_genes, n_families, norm_method, n_draws, k_start, k_step, k_max, tau, &
        max_pool_size, &
        pvalues_own, pvalues_fam, pvalues_orth, n_success, &
        neighborhood_size_own, neighborhood_size_fam, &
        neighborhood_size_orth, neighborhood_size_cancer, &
        cache, tmp_pool)

        type(sorted_data_t), intent(in) :: cancer_sorted
        !! Sorted cancer gene data
        type(sorted_data_t), intent(in) :: healthy_sorted
        !! Sorted healthy gene data
        integer(int32), intent(in) :: cancer_n_samples
        !! Number of cancer replicates
        integer(int32), intent(in) :: healthy_n_samples
        !! Number of healthy replicates
        integer(int32), intent(in) :: n_genes
        !! Total number of genes
        integer(int32), intent(in) :: n_families
        !! Total number of gene families
        real(real64), dimension(n_genes), intent(in) :: cancer_means
        !! Per-gene cancer expression means
        real(real64), dimension(n_genes), intent(in) :: healthy_means
        !! Per-gene healthy expression means
        real(real64), dimension(n_genes), intent(in) :: obs_own
        !! Observed gene-vs-own statistic for each gene
        real(real64), dimension(n_genes), intent(in) :: obs_fam
        !! Observed gene-vs-family statistic for each gene
        real(real64), dimension(n_genes), intent(in) :: obs_orth
        !! Observed gene-vs-ortholog statistic for each gene
        real(real64), dimension(n_families), intent(in) :: family_means
        !! Mean expression of each gene family (healthy)
        real(real64), dimension(n_families), intent(in) :: ortholog_means
        !! Mean expression of each ortholog set (healthy)
        integer(int32), dimension(n_genes), intent(in) :: valid_genes_own
        !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: valid_genes_fam
        !! 1 if the gene-vs-family p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: valid_genes_orth
        !! 1 if the gene-vs-ortholog p-value should be computed, 0 otherwise
        integer(int32), dimension(n_families), intent(in) :: family_sizes
        !! Number of genes in each family
        integer(int32), intent(in) :: is_ortholog_sum
        !! Total number of ortholog genes across all families
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
        !! Maps each gene index to its family index
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        integer(int32), intent(in) :: n_draws
        !! Number of Monte Carlo draws per p-value
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
        real(real64), dimension(n_genes), intent(out) :: pvalues_fam
        !! Output: gene-vs-family p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: pvalues_orth
        !! Output: gene-vs-ortholog p-values (-1 if not computed)
        integer(int32), intent(out) :: n_success
        !! Number of genes for which at least one p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own
        !! Healthy pool size used for gene-vs-own p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_fam
        !! Healthy pool size used for gene-vs-family p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_orth
        !! Healthy pool size used for gene-vs-ortholog p-value (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_cancer
        !! Cancer pool size used for this gene (-1 if not computed)
        type(family_cache_t), intent(in) :: cache
        !! Pre-computed family and ortholog residual pools
        real(real64), dimension(max_pool_size * 2), intent(inout) :: tmp_pool
        !! Shared work array for gather_residuals_helper staging

        integer(int32) :: i_gene, fam_id, r_c, r_h
        real(real64) :: mu_c, mu_h
        real(real64) :: obs_own_val, obs_fam_val, obs_orth_val
        real(real64), allocatable :: cancer_pool(:), healthy_pool_own(:)
        integer(int32) :: cancer_pool_count, healthy_pool_count_own
        integer(int32) :: dummy_ierr

        pvalues_own = -1.0_real64
        pvalues_fam = -1.0_real64
        pvalues_orth = -1.0_real64
        neighborhood_size_own = -1
        neighborhood_size_fam = -1
        neighborhood_size_orth = -1
        neighborhood_size_cancer = -1
        n_success = 0

        r_c = cancer_n_samples
        r_h = healthy_n_samples

        allocate(cancer_pool(max_pool_size))
        allocate(healthy_pool_own(max_pool_size))

        do i_gene = 1, n_genes
            mu_c = cancer_means(i_gene)
            mu_h = healthy_means(i_gene)
            fam_id = gene_to_family(i_gene)

            if (fam_id < 1 .or. fam_id > n_families) cycle
            if (.not. cache%is_cached(fam_id)) cycle

            call gather_residuals_helper(mu_c, cancer_sorted, &
                                         k_start, k_step, k_max, tau, &
                                         cancer_pool, cancer_pool_count, &
                                         max_pool_size, tmp_pool)
            if (cancer_pool_count < 10) cycle

            call gather_residuals_helper(mu_h, healthy_sorted, &
                                         k_start, k_step, k_max, tau, &
                                         healthy_pool_own, healthy_pool_count_own, &
                                         max_pool_size, tmp_pool)
            if (healthy_pool_count_own < 10) cycle

            obs_own_val = obs_own(i_gene)
            obs_fam_val = obs_fam(i_gene)
            obs_orth_val = obs_orth(i_gene)

            ! Skip genes with non-finite observed statistics
            if (obs_own_val /= obs_own_val .or. &
                obs_fam_val /= obs_fam_val .or. &
                obs_orth_val /= obs_orth_val) cycle

            if (valid_genes_own(i_gene) == 1) then
                call compute_pvalue(mu_c, r_c, cancer_pool(1:cancer_pool_count), &
                                    cancer_pool_count, &
                                    mu_h, r_h, healthy_pool_own(1:healthy_pool_count_own), &
                                    healthy_pool_count_own, &
                                    obs_own_val, n_draws, norm_method, &
                                    pvalues_own(i_gene), dummy_ierr)
                neighborhood_size_own(i_gene) = healthy_pool_count_own
            end if

            if (valid_genes_fam(i_gene) == 1) then
                call compute_pvalue(mu_c, r_c, cancer_pool(1:cancer_pool_count), &
                                    cancer_pool_count, &
                                    family_means(fam_id), r_h, &
                                    cache%family_pools(1:cache%family_pool_sizes(fam_id), fam_id), &
                                    cache%family_pool_sizes(fam_id), &
                                    obs_fam_val, n_draws, norm_method, &
                                    pvalues_fam(i_gene), dummy_ierr)
                neighborhood_size_fam(i_gene) = cache%family_pool_sizes(fam_id)
            end if

            if (valid_genes_orth(i_gene) == 1) then
                call compute_pvalue(mu_c, r_c, cancer_pool(1:cancer_pool_count), &
                                    cancer_pool_count, &
                                    ortholog_means(fam_id), r_h, &
                                    cache%orth_pools(1:cache%orth_pool_sizes(fam_id), fam_id), &
                                    cache%orth_pool_sizes(fam_id), &
                                    obs_orth_val, n_draws, norm_method, &
                                    pvalues_orth(i_gene), dummy_ierr)
                neighborhood_size_orth(i_gene) = cache%orth_pool_sizes(fam_id)
            end if

            neighborhood_size_cancer(i_gene) = cancer_pool_count
            n_success = n_success + 1
        end do

        deallocate(cancer_pool, healthy_pool_own)
    end subroutine compute_noise_pvalue_pipeline_helper

    !> Validate inputs, build sorted structures, pre-cache family pools, and run the pipeline.
    !|
    !| This is the alloc-layer entry point for the full noise-model pipeline.
    !| Internally it:
    !|   1. Validates all dimension and range arguments via `tox_errors`.
    !|   2. Calls `prepare_sorted_data` for both cancer and healthy data.
    !|   3. Pre-computes family and ortholog residual pools via `gather_residuals_helper`.
    !|   4. Delegates the per-gene loop to `compute_noise_pvalue_pipeline_helper`.
    subroutine compute_noise_pvalue_pipeline( &
        cancer_means, cancer_replicates, cancer_n_genes, cancer_n_samples, &
        healthy_means, healthy_replicates, healthy_n_genes, healthy_n_samples, &
        obs_own, obs_fam, obs_orth, &
        family_means, ortholog_means, &
        valid_genes_own, valid_genes_fam, valid_genes_orth, &
        family_sizes, is_ortholog_sum, gene_to_family, &
        n_genes, n_families, norm_method, n_draws, k_start, k_step, k_max, tau, &
        max_pool_size, &
        pvalues_own, pvalues_fam, pvalues_orth, n_success, &
        neighborhood_size_own, neighborhood_size_fam, &
        neighborhood_size_orth, neighborhood_size_cancer, &
        ierr)

        integer(int32), intent(in) :: cancer_n_genes
        !! Number of cancer genes
        integer(int32), intent(in) :: cancer_n_samples
        !! Number of cancer replicates
        integer(int32), intent(in) :: healthy_n_genes
        !! Number of healthy genes
        integer(int32), intent(in) :: healthy_n_samples
        !! Number of healthy replicates
        integer(int32), intent(in) :: n_genes
        !! Total number of genes for which p-values are computed
        integer(int32), intent(in) :: n_families
        !! Total number of gene families
        real(real64), dimension(cancer_n_genes), intent(in) :: cancer_means
        !! Per-gene cancer expression means
        real(real64), dimension(cancer_n_samples, cancer_n_genes), intent(in) :: cancer_replicates
        !! Cancer replicate expression matrix (cancer_n_samples x cancer_n_genes)
        real(real64), dimension(healthy_n_genes), intent(in) :: healthy_means
        !! Per-gene healthy expression means
        real(real64), dimension(healthy_n_samples, healthy_n_genes), intent(in) :: healthy_replicates
        !! Healthy replicate expression matrix (healthy_n_samples x healthy_n_genes)
        real(real64), dimension(n_genes), intent(in) :: obs_own
        !! Observed gene-vs-own statistic for each gene
        real(real64), dimension(n_genes), intent(in) :: obs_fam
        !! Observed gene-vs-family statistic for each gene
        real(real64), dimension(n_genes), intent(in) :: obs_orth
        !! Observed gene-vs-ortholog statistic for each gene
        real(real64), dimension(n_families), intent(in) :: family_means
        !! Mean expression of each gene family (healthy)
        real(real64), dimension(n_families), intent(in) :: ortholog_means
        !! Mean expression of each ortholog set (healthy)
        integer(int32), dimension(n_genes), intent(in) :: valid_genes_own
        !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: valid_genes_fam
        !! 1 if the gene-vs-family p-value should be computed, 0 otherwise
        integer(int32), dimension(n_genes), intent(in) :: valid_genes_orth
        !! 1 if the gene-vs-ortholog p-value should be computed, 0 otherwise
        integer(int32), dimension(n_families), intent(in) :: family_sizes
        !! Number of genes in each family
        integer(int32), intent(in) :: is_ortholog_sum
        !! Total number of ortholog genes across all families
        integer(int32), dimension(n_genes), intent(in) :: gene_to_family
        !! Maps each gene index to its family index (1-based)
        integer(int32), intent(in) :: norm_method
        !! 0 = linear scale; non-zero = log2(x+1) transform
        integer(int32), intent(in) :: n_draws
        !! Number of Monte Carlo draws per p-value
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
        real(real64), dimension(n_genes), intent(out) :: pvalues_fam
        !! Output: gene-vs-family p-values (-1 if not computed)
        real(real64), dimension(n_genes), intent(out) :: pvalues_orth
        !! Output: gene-vs-ortholog p-values (-1 if not computed)
        integer(int32), intent(out) :: n_success
        !! Number of genes for which at least one p-value was computed
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_own
        !! Healthy pool size used for gene-vs-own (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_fam
        !! Healthy pool size used for gene-vs-family (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_orth
        !! Healthy pool size used for gene-vs-ortholog (-1 if not computed)
        integer(int32), dimension(n_genes), intent(out) :: neighborhood_size_cancer
        !! Cancer pool size used for each gene (-1 if not computed)
        integer(int32), intent(out) :: ierr
        !! Error code

        type(sorted_data_t) :: cancer_sorted, healthy_sorted
        type(family_cache_t) :: cache
        real(real64), allocatable :: tmp_pool(:)
        integer(int32) :: fam, sort_ierr

        call set_ok(ierr)

        call validate_dimension_size(cancer_n_genes, ierr)
        call validate_dimension_size(cancer_n_samples, ierr)
        call validate_dimension_size(healthy_n_genes, ierr)
        call validate_dimension_size(healthy_n_samples, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_dimension_size(n_draws, ierr)
        call validate_dimension_size(k_start, ierr)
        call validate_dimension_size(k_step, ierr)
        call validate_dimension_size(k_max, ierr)
        call validate_dimension_size(max_pool_size, ierr)
        call validate_all_in_range_real(cancer_means, cancer_n_genes, ierr)
        call validate_all_in_range_real(healthy_means, healthy_n_genes, ierr)
        call validate_all_in_range_real(cancer_replicates, cancer_n_samples * cancer_n_genes, ierr)
        call validate_all_in_range_real(healthy_replicates, healthy_n_samples * healthy_n_genes, ierr)
        call validate_all_in_range_int(gene_to_family, n_genes, ierr, min=1, max=n_families)
        if (is_err(ierr)) return

        call init_random()

        call prepare_sorted_data(cancer_means, cancer_replicates, &
                                 cancer_n_samples, cancer_n_genes, cancer_sorted, sort_ierr)
        call set_err(ierr, sort_ierr)

        call prepare_sorted_data(healthy_means, healthy_replicates, &
                                 healthy_n_samples, healthy_n_genes, healthy_sorted, sort_ierr)
        call set_err(ierr, sort_ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(cache%family_pools(max_pool_size, n_families))
        M_ALLOCATE(cache%orth_pools(max_pool_size, n_families))
        M_ALLOCATE(cache%family_pool_sizes(n_families))
        M_ALLOCATE(cache%orth_pool_sizes(n_families))
        M_ALLOCATE(cache%is_cached(n_families))
        M_ALLOCATE(tmp_pool(max_pool_size * 2))

        cache%is_cached = .false.

        do fam = 1, n_families
            if (family_sizes(fam) > 0) then
                call gather_residuals_helper(family_means(fam), healthy_sorted, &
                                             k_start, k_step, k_max, tau, &
                                             cache%family_pools(:, fam), &
                                             cache%family_pool_sizes(fam), &
                                             max_pool_size, tmp_pool)
                call gather_residuals_helper(ortholog_means(fam), healthy_sorted, &
                                             k_start, k_step, k_max, tau, &
                                             cache%orth_pools(:, fam), &
                                             cache%orth_pool_sizes(fam), &
                                             max_pool_size, tmp_pool)
                cache%is_cached(fam) = .true.
            end if
        end do

        call compute_noise_pvalue_pipeline_helper( &
            cancer_sorted, healthy_sorted, &
            cancer_means, healthy_means, &
            cancer_n_samples, healthy_n_samples, &
            obs_own, obs_fam, obs_orth, &
            family_means, ortholog_means, &
            valid_genes_own, valid_genes_fam, valid_genes_orth, &
            family_sizes, is_ortholog_sum, gene_to_family, &
            n_genes, n_families, norm_method, n_draws, k_start, k_step, k_max, tau, &
            max_pool_size, &
            pvalues_own, pvalues_fam, pvalues_orth, n_success, & 
            neighborhood_size_own, neighborhood_size_fam, &
            neighborhood_size_orth, neighborhood_size_cancer, &
            cache, tmp_pool)

    end subroutine compute_noise_pvalue_pipeline

end module tox_noise_model

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

!> C-interoperable wrapper for `compute_noise_pvalue_pipeline`.
!|
!| Performs null-pointer checks via `M_CHECK_IERR_NON_NULL` / `M_CHECK_NON_NULL`,
!| then delegates unconditionally to the validated Fortran entry point.
!| No computation is performed here.
subroutine compute_noise_pvalues_pipeline_C( &
    cancer_means, cancer_replicates, cancer_n_genes, cancer_n_samples, &
    healthy_means, healthy_replicates, healthy_n_genes, healthy_n_samples, &
    obs_own, obs_fam, obs_orth, &
    family_means, ortholog_means, &
    valid_genes_own, valid_genes_fam, valid_genes_orth, &
    family_sizes, is_ortholog_sum, gene_to_family, &
    n_genes, n_families, norm_method, n_draws, k_start, k_step, k_max, tau, &
    max_pool_size, &
    pvalues_own, pvalues_fam, pvalues_orth, n_success, &
    neighborhood_size_own, neighborhood_size_fam, &
    neighborhood_size_orth, neighborhood_size_cancer, &
    ierr) bind(C, name="compute_noise_pvalues_pipeline_C")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_noise_model, only: compute_noise_pvalue_pipeline
    use safeguard
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: cancer_n_genes
    !! Number of cancer genes
    integer(c_int), intent(in), target :: cancer_n_samples
    !! Number of cancer replicates
    integer(c_int), intent(in), target :: healthy_n_genes
    !! Number of healthy genes
    integer(c_int), intent(in), target :: healthy_n_samples
    !! Number of healthy replicates
    integer(c_int), intent(in), target :: n_genes
    !! Total number of genes for which p-values are computed
    integer(c_int), intent(in), target :: n_families
    !! Total number of gene families
    real(c_double), dimension(cancer_n_genes), intent(in), target :: cancer_means
    !! Per-gene cancer expression means
    real(c_double), dimension(cancer_n_samples, cancer_n_genes), intent(in), target :: cancer_replicates
    !! Cancer replicate expression matrix (cancer_n_samples x cancer_n_genes)
    real(c_double), dimension(healthy_n_genes), intent(in), target :: healthy_means
    !! Per-gene healthy expression means
    real(c_double), dimension(healthy_n_samples, healthy_n_genes), intent(in), target :: healthy_replicates
    !! Healthy replicate expression matrix (healthy_n_samples x healthy_n_genes)
    real(c_double), dimension(n_genes), intent(in), target :: obs_own
    !! Observed gene-vs-own statistic for each gene
    real(c_double), dimension(n_genes), intent(in), target :: obs_fam
    !! Observed gene-vs-family statistic for each gene
    real(c_double), dimension(n_genes), intent(in), target :: obs_orth
    !! Observed gene-vs-ortholog statistic for each gene
    real(c_double), dimension(n_families), intent(in), target :: family_means
    !! Mean expression of each gene family (healthy)
    real(c_double), dimension(n_families), intent(in), target :: ortholog_means
    !! Mean expression of each ortholog set (healthy)
    integer(c_int), dimension(n_genes), intent(in), target :: valid_genes_own
    !! 1 if the gene-vs-own p-value should be computed, 0 otherwise
    integer(c_int), dimension(n_genes), intent(in), target :: valid_genes_fam
    !! 1 if the gene-vs-family p-value should be computed, 0 otherwise
    integer(c_int), dimension(n_genes), intent(in), target :: valid_genes_orth
    !! 1 if the gene-vs-ortholog p-value should be computed, 0 otherwise
    integer(c_int), dimension(n_families), intent(in), target :: family_sizes
    !! Number of genes in each family
    integer(c_int), intent(in), target :: is_ortholog_sum
    !! Total number of ortholog genes across all families
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_family
    !! Maps each gene index to its family index (1-based)
    integer(c_int), intent(in), target :: norm_method
    !! 0 = linear scale; non-zero = log2(x+1) transform
    integer(c_int), intent(in), target :: n_draws
    !! Number of Monte Carlo draws per p-value
    integer(c_int), intent(in), target :: k_start
    !! Minimum pool size before adaptive stopping is applied
    integer(c_int), intent(in), target :: k_step
    !! Residuals added per adaptive round before re-evaluating
    integer(c_int), intent(in), target :: k_max
    !! Hard upper limit on residual pool size
    real(c_double), intent(in), target :: tau
    !! Relative-change threshold for adaptive pool growth
    integer(c_int), intent(in), target :: max_pool_size
    !! Maximum size of residual pools
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_own
    !! Output: gene-vs-own p-values (-1 if not computed)
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_fam
    !! Output: gene-vs-family p-values (-1 if not computed)
    real(c_double), dimension(n_genes), intent(out), target :: pvalues_orth
    !! Output: gene-vs-ortholog p-values (-1 if not computed)
    integer(c_int), intent(out), target :: n_success
    !! Number of genes for which at least one p-value was computed
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_own
    !! Healthy pool size used for gene-vs-own (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_fam
    !! Healthy pool size used for gene-vs-family (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_orth
    !! Healthy pool size used for gene-vs-ortholog (-1 if not computed)
    integer(c_int), dimension(n_genes), intent(out), target :: neighborhood_size_cancer
    !! Cancer pool size used for each gene (-1 if not computed)
    integer(c_int), intent(out), target :: ierr
    !! Error code: 0 = success

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(cancer_n_genes)
    M_CHECK_NON_NULL(cancer_n_samples)
    M_CHECK_NON_NULL(healthy_n_genes)
    M_CHECK_NON_NULL(healthy_n_samples)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(cancer_means)
    M_CHECK_NON_NULL(cancer_replicates)
    M_CHECK_NON_NULL(healthy_means)
    M_CHECK_NON_NULL(healthy_replicates)
    M_CHECK_NON_NULL(obs_own)
    M_CHECK_NON_NULL(obs_fam)
    M_CHECK_NON_NULL(obs_orth)
    M_CHECK_NON_NULL(family_means)
    M_CHECK_NON_NULL(ortholog_means)
    M_CHECK_NON_NULL(valid_genes_own)
    M_CHECK_NON_NULL(valid_genes_fam)
    M_CHECK_NON_NULL(valid_genes_orth)
    M_CHECK_NON_NULL(family_sizes)
    M_CHECK_NON_NULL(is_ortholog_sum)
    M_CHECK_NON_NULL(gene_to_family)
    M_CHECK_NON_NULL(norm_method)
    M_CHECK_NON_NULL(n_draws)
    M_CHECK_NON_NULL(k_start)
    M_CHECK_NON_NULL(k_step)
    M_CHECK_NON_NULL(k_max)
    M_CHECK_NON_NULL(tau)
    M_CHECK_NON_NULL(max_pool_size)
    M_CHECK_NON_NULL(pvalues_own)
    M_CHECK_NON_NULL(pvalues_fam)
    M_CHECK_NON_NULL(pvalues_orth)
    M_CHECK_NON_NULL(n_success)
    M_CHECK_NON_NULL(neighborhood_size_own)
    M_CHECK_NON_NULL(neighborhood_size_fam)
    M_CHECK_NON_NULL(neighborhood_size_orth)
    M_CHECK_NON_NULL(neighborhood_size_cancer)

    call compute_noise_pvalue_pipeline( &
        cancer_means, cancer_replicates, cancer_n_genes, cancer_n_samples, &
        healthy_means, healthy_replicates, healthy_n_genes, healthy_n_samples, &
        obs_own, obs_fam, obs_orth, &
        family_means, ortholog_means, &
        valid_genes_own, valid_genes_fam, valid_genes_orth, &
        family_sizes, is_ortholog_sum, gene_to_family, &
        n_genes, n_families, norm_method, n_draws, k_start, k_step, k_max, tau, &
        max_pool_size, &
        pvalues_own, pvalues_fam, pvalues_orth, n_success, &
        neighborhood_size_own, neighborhood_size_fam, &
        neighborhood_size_orth, neighborhood_size_cancer, &
        ierr)

end subroutine compute_noise_pvalues_pipeline_C
