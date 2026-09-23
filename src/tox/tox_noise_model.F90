#include <src/macros.h>

!!! TO DO
!!! 1. Implement error code to argument number mapping
!!! 2. Find solution for family based analysis

!> Noise model for exact p-value computation comparing case and control gene expression.
!|
!| This module and `noise_model_exact` (in `tox_noise_model_exact.F90`) are
!| identical EXCEPT for how the `own` comparison's null is corrected for the
!| individual-vs-mean scale mismatch (the observed statistic is a difference of
!| MEANS, so the null must be built at the mean level):
!|   - `noise_model_exact`: scales each `own` pool by `1/sqrt(n_replicates)`
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
!|   - Compute exact p-values by counting, via sorted-pool binary search, every
!|     pairwise residual difference that meets or exceeds the observed statistic
!|   - Run the full pipeline over all genes with pre-cached family / ortholog pools
module noise_model
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use tox_errors, only: set_ok, set_err, is_err, &
                          ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_NAN_INF, ERR_ALLOC_FAIL, &
                          validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int
    use f42_utils, only: sort_real, sort_integer, init_random
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

    ! Note: the per-gene p-value is computed exactly at every pool size (no
    ! MAX_EXACT_COMBINATIONS enumeration/Monte-Carlo threshold). The only RNG use is
    ! the bootstrap that builds the `own` mean-difference null, below.

    integer(int32), parameter :: N_BOOTSTRAP_DRAWS = 10000_int32
        !! Number of bootstrap resamples used to build the `own` mean-difference null
        !! (sampled null constructions only; the enumerated blocked null uses none)

    integer(int32), parameter :: NULL_METHOD_POOLED = 0_int32
        !! `null_method`: resample `n_rep` residuals iid from the POOLED neighbourhood.
        !! Historical behaviour of this module. One draw can combine a residual from a
        !! quiet gene with one from a noisy gene, which no real gene's mean does.

    integer(int32), parameter :: NULL_METHOD_BLOCKED = 1_int32
        !! `null_method`: GENE-BLOCKED null — pick one neighbour gene, then resample
        !! within that gene, so every null mean carries a single coherent noise level.
        !! Computed by exact multiset enumeration when that fits under
        !! `BLOCKED_MAX_ENUM_VALUES`, otherwise sampled.

    integer(int64), parameter :: BLOCKED_MAX_ENUM_VALUES = 8000_int64
        !! Cap on enumerated values per side, `C(2n-1, n) * n_genes_pool`; above it the
        !! blocked null is sampled instead. The binding constraint is RUNTIME, not
        !! memory: the per-gene cost is one sort plus `N` binary searches over `N`
        !! values, and `N` grows like `C(2n-1, n)`, i.e. ~4x per extra replicate.
        !! Measured at 2,000 genes, `k_max = 50`, against this module's own 10,000-draw
        !! pooled bootstrap on the same data (gfortran -O2, single core):
        !!
        !!   n_rep :   3      4      5       6       7
        !!   values:  500  1,750  6,300  23,100  85,800   (= C(2n-1,n) * 50)
        !!   blocked: 0.41s  1.71s  7.34s  30.6s   131s
        !!   pooled : 1.45s  1.79s  2.18s   2.63s  3.03s
        !!   ratio  : 0.28x  0.95x  3.4x    12x     43x
        !!
        !! So at `n_rep = 3` — the case where the p-value floor actually bites —
        !! enumeration is ~3.5x FASTER than the bootstrap it replaces AND removes the
        !! floor; by `n_rep = 6` it costs 12x. The default 8,000 enumerates through
        !! `n_rep = 5` at `k_max = 50` and samples from 6 up. Raise it if a lower floor
        !! at higher replicate counts is worth the time (memory is not the issue: 8,000
        !! values is 64 KB per side). The obvious optimisation before raising it far is
        !! to sort BOTH sides once and sweep them with two pointers instead of doing
        !! `N` binary searches, which removes a `log N` factor from the dominant term.

    integer(int32), parameter :: RNG_BUFFER_MAX = 1048576_int32
        !! Upper bound on the pre-drawn random-number buffer (8 MB). Random numbers
        !! are produced one CHUNK of draws at a time by a single `random_number` array
        !! fill instead of one call per draw: at `n_rep = 3` and 10,000 draws the
        !! whole gene needs ONE RNG call instead of 10,000. The buffer is allocated
        !! once by the alloc layer and reused for every gene.

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
                                                    rbuf, p_value)
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
        real(real64), intent(inout) :: rbuf(:)
        !! Pre-allocated random-number buffer (reused across genes)
        real(real64), intent(out) :: p_value
        !! Bootstrapped p-value

        integer(int32) :: i_boot, i_rep, idx, count_ge, per_iter, chunk, j, base
        real(real64) :: sum_case, sum_control, null_stat

        ! Random numbers are drawn a CHUNK OF ITERATIONS at a time with a single
        ! `random_number` array fill, instead of one fill (or worse, one scalar call)
        ! per iteration: `rbuf` holds `chunk * per_iter` values, laid out iteration by
        ! iteration, with the case indices first and the control indices after. With
        ! the default buffer this is ONE RNG call per gene at typical replicate
        ! counts, which is where the intrinsic's per-call overhead was going. The
        ! chunk is sized from the buffer so memory stays bounded when `n_rep` is
        ! large. `init_random(42)` (called once by the pipeline) seeds this stream,
        ! so results stay reproducible for a fixed gene order.
        per_iter = n_rep_case + n_rep_control
        chunk = max(1, min(n_boot, int(size(rbuf), int32) / per_iter))

        count_ge = 0
        i_boot = 0
        do while (i_boot < n_boot)
            j = min(chunk, n_boot - i_boot)
            call random_number(rbuf(1:j * per_iter))
            do base = 0, (j - 1) * per_iter, per_iter
                sum_case = 0.0_real64
                do i_rep = 1, n_rep_case
                    idx = min(int(rbuf(base + i_rep) * real(n_pool_case, real64), int32) + 1, n_pool_case)
                    sum_case = sum_case + pool_case(idx)
                end do

                sum_control = 0.0_real64
                do i_rep = 1, n_rep_control
                    idx = min(int(rbuf(base + n_rep_case + i_rep) * real(n_pool_control, real64), int32) + 1, &
                              n_pool_control)
                    sum_control = sum_control + pool_control(idx)
                end do

                null_stat = abs(sum_case / real(n_rep_case, real64) - &
                                sum_control / real(n_rep_control, real64))
                if (null_stat >= observed_statistic_abs) count_ge = count_ge + 1
            end do
            i_boot = i_boot + j
        end do

        p_value = real(count_ge + 1, real64) / real(n_boot + 1, real64)
    end subroutine compute_pvalue_bootstrap_mean_helper

    ! =========================================================================
    ! gene-blocked mean null: multiset enumeration + weighted exact tail count
    ! =========================================================================

    !> Number of distinct MULTISETS drawn when resampling `n_rep` of `n_rep` values
    !| with replacement: `C(2*n_rep - 1, n_rep)`.
    !|
    !| Computed as the running product `C(n+j, j)`, which is an integer at every
    !| step, so the alternating multiply/divide stays exact. Returns -1 on int64
    !| overflow (unreachable for any `n_rep` this module will enumerate, but the
    !| guard is what makes the cap check below safe).
    !|
    !|   n_rep :  2   3    4    5    6    8      10       12
    !|   count :  3  10   35  126  462  6435  92378  1352078
    pure function n_multisets_helper(n_rep) result(n_ms)
        integer(int32), intent(in) :: n_rep
        !! Replicate count (= resample size = number of source values)
        integer(int64) :: n_ms
        !! C(2*n_rep - 1, n_rep), or -1 on overflow

        integer(int32) :: j

        if (n_rep < 1) then
            n_ms = 0_int64
            return
        end if
        n_ms = 1_int64
        do j = 1, n_rep - 1
            if (n_ms > huge(0_int64) / int(n_rep + j, int64)) then
                n_ms = -1_int64
                return
            end if
            n_ms = (n_ms * int(n_rep + j, int64)) / int(j, int64)
        end do
    end function n_multisets_helper

    !> Binomial coefficient `C(n, k)` as a real64, for the multinomial weights.
    !|
    !| Exact for every (n, k) reachable under `BLOCKED_MAX_ENUM_VALUES`: the largest
    !| weight enumerated is `n_rep!` and `12! = 4.79e8`, far inside the 2^53 range
    !| where real64 represents integers exactly.
    pure function binom_helper(n, k) result(c)
        integer(int32), intent(in) :: n
        !! Total count
        integer(int32), intent(in) :: k
        !! Chosen count
        real(real64) :: c
        !! C(n, k)

        integer(int32) :: j, kk

        if (k < 0 .or. k > n) then
            c = 0.0_real64
            return
        end if
        kk = min(k, n - k)
        c = 1.0_real64
        do j = 1, kk
            c = (c * real(n - kk + j, real64)) / real(j, real64)
        end do
    end function binom_helper

    !> Enumerate every distinct bootstrap MEAN of one gene's residuals, with weights.
    !|
    !| Resampling `n_rep` of a gene's `n_rep` residuals with replacement produces
    !| `n_rep**n_rep` ordered tuples but only `C(2*n_rep - 1, n_rep)` distinct
    !| multisets. This walks the non-decreasing index tuples
    !| `1 <= i_1 <= i_2 <= ... <= i_n <= n` as an odometer and emits, for each, the
    !| mean of the selected residuals together with its multinomial weight
    !| `n! / prod_j(m_j!)` (`m_j` = run lengths in the tuple), so the weighted set is
    !| exactly equivalent to the full `n_rep**n_rep` enumeration but far smaller:
    !| 10 values instead of 27 at `n_rep = 3`, 462 instead of 46,656 at `n_rep = 6`.
    !|
    !| The weights sum to `n_rep**n_rep` by construction; the caller relies on that
    !| when forming the pair-weight denominator.
    !|
    !| `means_out` / `weights_out` must each have room for at least
    !| `n_multisets_helper(n_rep)` entries.
    pure subroutine enumerate_gene_means_helper(resid, n_rep, means_out, weights_out, n_out)
        real(real64), intent(in) :: resid(:)
        !! This gene's residuals (first `n_rep` entries are used)
        integer(int32), intent(in) :: n_rep
        !! Replicate count = resample size
        real(real64), intent(out) :: means_out(:)
        !! Output: distinct bootstrap means
        real(real64), intent(out) :: weights_out(:)
        !! Output: multinomial weight of each mean
        integer(int32), intent(out) :: n_out
        !! Number of entries written

        integer(int32) :: idx(n_rep)
        integer(int32) :: i, p, run_len, rem
        real(real64) :: s, w

        n_out = 0
        if (n_rep < 1) return

        idx = 1
        do
            n_out = n_out + 1

            s = 0.0_real64
            do i = 1, n_rep
                s = s + resid(idx(i))
            end do
            means_out(n_out) = s / real(n_rep, real64)

            ! Multinomial weight from the run lengths of the (non-decreasing) tuple:
            ! n! / prod_j(m_j!) built as a product of binomials, which keeps every
            ! intermediate an exact integer.
            w = 1.0_real64
            rem = n_rep
            i = 1
            do while (i <= n_rep)
                run_len = 1
                do while (i + run_len <= n_rep)
                    if (idx(i + run_len) /= idx(i)) exit
                    run_len = run_len + 1
                end do
                w = w * binom_helper(rem, run_len)
                rem = rem - run_len
                i = i + run_len
            end do
            weights_out(n_out) = w

            ! Advance the odometer: find the rightmost position that can still be
            ! incremented, bump it, and reset everything to its right to the same
            ! value (keeping the tuple non-decreasing).
            p = n_rep
            do while (p >= 1)
                if (idx(p) /= n_rep) exit
                p = p - 1
            end do
            if (p < 1) exit
            idx(p) = idx(p) + 1
            if (p < n_rep) idx(p + 1:n_rep) = idx(p)
        end do
    end subroutine enumerate_gene_means_helper

    !> Build the complete gene-blocked mean null for one gathered pool.
    !|
    !| `gather_residuals_helper` appends whole genes, so the flat pool is exactly
    !| `n_pool / n_rep` contiguous blocks of `n_rep` residuals — one per neighbour
    !| gene. (A partial trailing block can only arise from the `max_pool_size` cap;
    !| it is dropped.) Each block is enumerated separately, which is what makes the
    !| null GENE-BLOCKED: every null mean is formed from ONE gene's residuals, the
    !| way a real gene's mean is, instead of mixing residuals from a quiet gene and
    !| a noisy one as an iid draw from the pooled residuals does.
    !|
    !| Returns `ok = .false.` (and writes nothing) when the enumeration would exceed
    !| `max_values`; the caller then falls back to the sampled blocked null.
    pure subroutine build_blocked_means_helper(pool, n_pool, n_rep, &
                                               means_out, weights_out, n_out, &
                                               max_values, ok)
        real(real64), intent(in) :: pool(:)
        !! Flat residual pool as returned by `gather_residuals_helper`
        integer(int32), intent(in) :: n_pool
        !! Number of valid residuals in `pool`
        integer(int32), intent(in) :: n_rep
        !! Residuals contributed by each gene (= replicate count)
        real(real64), intent(out) :: means_out(:)
        !! Output: all enumerated block means
        real(real64), intent(out) :: weights_out(:)
        !! Output: their multinomial weights
        integer(int32), intent(out) :: n_out
        !! Number of entries written (0 when `ok` is .false.)
        integer(int64), intent(in) :: max_values
        !! Cap on the number of enumerated values
        logical, intent(out) :: ok
        !! .true. if the enumeration was performed

        integer(int32) :: n_genes_pool, g, m, base
        integer(int64) :: n_ms

        n_out = 0
        ok = .false.
        if (n_rep < 1) return

        n_genes_pool = n_pool / n_rep
        if (n_genes_pool < 1) return

        n_ms = n_multisets_helper(n_rep)
        if (n_ms < 1_int64) return
        if (n_ms > max_values / int(n_genes_pool, int64)) return
        if (n_ms * int(n_genes_pool, int64) > int(min(size(means_out), size(weights_out)), int64)) return

        do g = 1, n_genes_pool
            base = (g - 1) * n_rep
            call enumerate_gene_means_helper(pool(base + 1:base + n_rep), n_rep, &
                                             means_out(n_out + 1:), weights_out(n_out + 1:), m)
            n_out = n_out + m
        end do
        ok = .true.
    end subroutine build_blocked_means_helper

    !> Count entries of an ascending array below a threshold, via binary search.
    !|
    !| Returns the number of entries in `arr(1:n)` (assumed ascending) that are
    !| `< x` when `inclusive` is .false., or `<= x` when .true. O(log n).
    pure function count_below_real_helper(arr, n, x, inclusive) result(cnt)
        real(real64), intent(in) :: arr(:)
        !! Ascending array
        integer(int32), intent(in) :: n
        !! Number of valid entries
        real(real64), intent(in) :: x
        !! Threshold
        logical, intent(in) :: inclusive
        !! .true. counts entries equal to `x` as well
        integer(int32) :: cnt
        !! Number of entries below (or at) the threshold

        integer(int32) :: lo, hi, mid
        logical :: take

        lo = 1
        hi = n
        cnt = 0
        do while (lo <= hi)
            mid = lo + (hi - lo) / 2
            if (inclusive) then
                take = (arr(mid) <= x)
            else
                take = (arr(mid) < x)
            end if
            if (take) then
                cnt = mid
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function count_below_real_helper

    !> Exact p-value from the enumerated gene-blocked mean-difference null.
    !|
    !| Scores `|observed|` against every case-mean vs control-mean pair, weighted by
    !| the multinomial weights, using the same sorted-array + two-binary-searches
    !| tail count the individual-residual null uses:
    !|
    !|   p = ( sum_{i,j} w_a(i) w_b(j) [ |a_i - b_j| >= |obs| ] + 1 ) / ( W_a W_b + 1 )
    !|
    !| with `W = sum(w) = n_genes_pool * n_rep**n_rep`. There is no draw count and no
    !| RNG: this is the `n_boot -> infinity` limit of the sampled blocked bootstrap,
    !| computed in closed form.
    !|
    !| Resolution: the floor is `1 / (W_a W_b + 1)`. At `k = 30` neighbour genes and
    !| `n_rep = 3` that is `810 * 810 = 656,100` -> 1.5e-6, versus 1.2e-4 for the
    !| 90x90 individual-residual pairing and 4.0e-5 for a 25,000-draw bootstrap.
    !| NOTE this removes an ARTIFICIAL floor; it does not add information. The
    !| effective sample size behind a low-replicate pool is orders of magnitude
    !| smaller than `W_a W_b`, so p-values far below ~1e-2 at `n_rep = 3` are not
    !| supported by the pool that produced them.
    !|
    !| `W_a * W_b` can reach ~1e23 at large `n_rep`, past exact integer range in
    !| real64; the accumulation is therefore relatively (not absolutely) exact, which
    !| is immaterial for a ratio.
    !|
    !| Work arrays `perm`, `stack_left`, `stack_right`, `sorted_b` must hold at least
    !| `n_b` entries and `cumw_b` at least `n_b + 1` (indexed from 0).
    pure subroutine compute_pvalue_blocked_exact_helper(means_a, weights_a, n_a, &
                                                        means_b, weights_b, n_b, &
                                                        observed_statistic_abs, &
                                                        perm, stack_left, stack_right, &
                                                        sorted_b, cumw_b, p_value)
        real(real64), intent(in) :: means_a(:)
        !! Case-side enumerated block means
        real(real64), intent(in) :: weights_a(:)
        !! Their multinomial weights
        integer(int32), intent(in) :: n_a
        !! Number of case-side entries
        real(real64), intent(in) :: means_b(:)
        !! Control-side enumerated block means
        real(real64), intent(in) :: weights_b(:)
        !! Their multinomial weights
        integer(int32), intent(in) :: n_b
        !! Number of control-side entries
        real(real64), intent(in) :: observed_statistic_abs
        !! Absolute observed statistic
        integer(int32), intent(inout) :: perm(:)
        !! Work array: sort permutation (length >= n_b)
        integer(int32), intent(out) :: stack_left(:)
        !! Work array: quicksort stack (length >= n_b)
        integer(int32), intent(out) :: stack_right(:)
        !! Work array: quicksort stack (length >= n_b)
        real(real64), intent(out) :: sorted_b(:)
        !! Work array: ascending control means (length >= n_b)
        real(real64), intent(out) :: cumw_b(0:)
        !! Work array: prefix sums of the sorted control weights (length >= n_b + 1)
        real(real64), intent(out) :: p_value
        !! Exact weighted p-value

        integer(int32) :: i, hi_idx, lo_idx
        real(real64) :: total_w_a, total_w_b, inside_w, count_ge, a

        do i = 1, n_b
            perm(i) = i
        end do
        call sort_real(means_b(1:n_b), perm(1:n_b), stack_left(1:n_b), stack_right(1:n_b))

        cumw_b(0) = 0.0_real64
        do i = 1, n_b
            sorted_b(i) = means_b(perm(i))
            cumw_b(i) = cumw_b(i - 1) + weights_b(perm(i))
        end do
        total_w_b = cumw_b(n_b)
        total_w_a = sum(weights_a(1:n_a))

        ! Tail weight = sum over case means of the control weight OUTSIDE the open
        ! interval (a - t, a + t), i.e. every pair with |a - b| >= t.
        count_ge = 0.0_real64
        do i = 1, n_a
            a = means_a(i)
            hi_idx = count_below_real_helper(sorted_b, n_b, a + observed_statistic_abs, .false.)
            lo_idx = count_below_real_helper(sorted_b, n_b, a - observed_statistic_abs, .true.)
            inside_w = cumw_b(hi_idx) - cumw_b(lo_idx)
            if (inside_w < 0.0_real64) inside_w = 0.0_real64
            count_ge = count_ge + weights_a(i) * (total_w_b - inside_w)
        end do

        p_value = (count_ge + 1.0_real64) / (total_w_a * total_w_b + 1.0_real64)
        if (p_value > 1.0_real64) p_value = 1.0_real64
    end subroutine compute_pvalue_blocked_exact_helper

    !> Sampled gene-blocked mean null (fallback when enumeration is too large).
    !|
    !| Same null as `compute_pvalue_blocked_exact_helper`, drawn instead of
    !| enumerated: each draw picks one neighbour GENE per side and then resamples
    !| `n_rep` residuals from within that gene, so the draw carries a single coherent
    !| noise level. Used only when `C(2n-1, n) * n_genes_pool` exceeds
    !| `BLOCKED_MAX_ENUM_VALUES` (roughly `n_rep > 9` at `k_max = 50`).
    !|
    !| RNG: all random numbers for a CHUNK of draws are produced by ONE
    !| `random_number` array fill, so the number of RNG calls per gene is
    !| `ceil(n_boot / chunk)` rather than `n_boot` — typically 1. `rbuf` is allocated
    !| once by the alloc layer and reused for every gene.
    subroutine compute_pvalue_blocked_bootstrap_helper(pool_case, n_pool_case, &
                                                       pool_control, n_pool_control, &
                                                       n_rep_case, n_rep_control, &
                                                       observed_statistic_abs, n_boot, &
                                                       rbuf, p_value)
        real(real64), intent(in) :: pool_case(:)
        !! Case residual pool (contiguous blocks of `n_rep_case` per gene)
        integer(int32), intent(in) :: n_pool_case
        !! Number of valid case residuals
        real(real64), intent(in) :: pool_control(:)
        !! Control residual pool (contiguous blocks of `n_rep_control` per gene)
        integer(int32), intent(in) :: n_pool_control
        !! Number of valid control residuals
        integer(int32), intent(in) :: n_rep_case
        !! Case replicate count (= resample size, = residuals per gene)
        integer(int32), intent(in) :: n_rep_control
        !! Control replicate count
        real(real64), intent(in) :: observed_statistic_abs
        !! Absolute observed statistic
        integer(int32), intent(in) :: n_boot
        !! Number of bootstrap draws
        real(real64), intent(inout) :: rbuf(:)
        !! Pre-allocated random-number buffer (reused across genes)
        real(real64), intent(out) :: p_value
        !! Bootstrapped p-value

        integer(int32) :: n_genes_case, n_genes_control, per_iter, chunk
        integer(int32) :: i_boot, i_rep, j, base, g_case, g_control, idx, count_ge
        real(real64) :: sum_case, sum_control

        n_genes_case = n_pool_case / n_rep_case
        n_genes_control = n_pool_control / n_rep_control
        if (n_genes_case < 1 .or. n_genes_control < 1) then
            p_value = 1.0_real64
            return
        end if

        ! One gene index plus n_rep within-gene indices, per side.
        per_iter = (n_rep_case + 1) + (n_rep_control + 1)
        chunk = max(1, min(n_boot, int(size(rbuf), int32) / per_iter))

        count_ge = 0
        i_boot = 0
        do while (i_boot < n_boot)
            j = min(chunk, n_boot - i_boot)
            call random_number(rbuf(1:j * per_iter))
            do base = 0, (j - 1) * per_iter, per_iter
                g_case = min(int(rbuf(base + 1) * real(n_genes_case, real64), int32) + 1, n_genes_case)
                sum_case = 0.0_real64
                do i_rep = 1, n_rep_case
                    idx = min(int(rbuf(base + 1 + i_rep) * real(n_rep_case, real64), int32) + 1, n_rep_case)
                    sum_case = sum_case + pool_case((g_case - 1) * n_rep_case + idx)
                end do

                g_control = min(int(rbuf(base + n_rep_case + 2) * real(n_genes_control, real64), int32) + 1, &
                                n_genes_control)
                sum_control = 0.0_real64
                do i_rep = 1, n_rep_control
                    idx = min(int(rbuf(base + n_rep_case + 2 + i_rep) * real(n_rep_control, real64), int32) + 1, &
                              n_rep_control)
                    sum_control = sum_control + pool_control((g_control - 1) * n_rep_control + idx)
                end do

                if (abs(sum_case / real(n_rep_case, real64) - &
                        sum_control / real(n_rep_control, real64)) >= observed_statistic_abs) &
                    count_ge = count_ge + 1
            end do
            i_boot = i_boot + j
        end do

        p_value = real(count_ge + 1, real64) / real(n_boot + 1, real64)
    end subroutine compute_pvalue_blocked_bootstrap_helper

    ! =========================================================================
    ! compute_noise_pvalue_pipeline
    ! =========================================================================

    !> Core pipeline: compute per-gene noise p-values using pre-built data structures.
    !|
    !| Iterates over all genes and computes, for each, `pvalues_own`: the gene vs. its
    !| own matched neighbourhood. For both the case and control sides an adaptive kNN
    !| residual pool is gathered in mean-expression space (`gather_residuals_helper`),
    !| and the `own` mean-difference null is built by bootstrapping directly from those
    !| two pools (`compute_pvalue_bootstrap_mean_helper`).
    !|
    !| Requires `sorted_case` and `sorted_control` to already be built (via
    !| `prepare_sorted_data`).
    !|
    !| No allocation: the two per-gene residual pools are pre-allocated by the caller
    !| and passed in as `intent(inout)` work arrays. The `own` mean-difference null is
    !| bootstrapped from those pools, which draws from the global RNG; that is the only
    !| impure operation, so this subroutine is not `pure`.
    subroutine compute_noise_pvalue_pipeline_helper( &
        sorted_case, sorted_control, &
        means_case, means_control, &
        observed_statistic_own, &
        compute_pvalue_own, &
        n_genes, k_start, k_step, k_max, tau, null_method, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        tmp_pool_case, tmp_pool_control_own, &
        blk_means_case, blk_weights_case, blk_means_control, blk_weights_control, &
        blk_perm, blk_stack_left, blk_stack_right, blk_sorted, blk_cumw, rbuf, &
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
        integer(int32), intent(in) :: null_method
        !! `NULL_METHOD_POOLED` (0) = iid resample from the pooled neighbourhood;
        !! `NULL_METHOD_BLOCKED` (1) = gene-blocked null, enumerated where it fits
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
        real(real64), dimension(:), intent(inout) :: blk_means_case
        !! Work array: enumerated gene-blocked case means (length = enumeration capacity)
        real(real64), dimension(:), intent(inout) :: blk_weights_case
        !! Work array: their multinomial weights
        real(real64), dimension(:), intent(inout) :: blk_means_control
        !! Work array: enumerated gene-blocked control means
        real(real64), dimension(:), intent(inout) :: blk_weights_control
        !! Work array: their multinomial weights
        integer(int32), dimension(:), intent(inout) :: blk_perm
        !! Work array: sort permutation for the control-side enumeration
        integer(int32), dimension(:), intent(inout) :: blk_stack_left
        !! Work array: quicksort stack
        integer(int32), dimension(:), intent(inout) :: blk_stack_right
        !! Work array: quicksort stack
        real(real64), dimension(:), intent(inout) :: blk_sorted
        !! Work array: ascending control means
        real(real64), dimension(0:), intent(inout) :: blk_cumw
        !! Work array: prefix sums of the sorted control weights (0-based)
        real(real64), dimension(:), intent(inout) :: rbuf
        !! Work array: pre-drawn random numbers for the sampled null constructions
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_gene
        real(real64) :: mean_case_val, mean_control_val
        real(real64) :: observed_statistic_own_val
        integer(int32) :: n_pool_case, n_pool_control_own
        integer(int32) :: n_rep_case, n_rep_control, n_blk_case, n_blk_control
        logical :: enum_ok_case, enum_ok_control

        call set_ok(ierr)

        pvalues_own = -1.0_real64
        neighborhood_size_own_case = -1
        neighborhood_size_own_control = -1
        neighborhood_size_case = -1
        n_genes_with_pvalue = 0

        n_rep_case = sorted_case%max_resid_per_gene
        n_rep_control = sorted_control%max_resid_per_gene

        ! This loop stays sequential whenever a SAMPLED null is in use: those draw
        ! from the global RNG stream seeded once by init_random(42), so reordering
        ! genes would change the draws and break reproducibility. Note that
        ! null_method = NULL_METHOD_BLOCKED uses NO RNG at all whenever the
        ! enumeration fits (the common case at low replicate counts), so that path
        ! could be parallelised as-is; it is left sequential here so one code path
        ! serves both, and because the enumerated null is already much cheaper than
        ! N_BOOTSTRAP_DRAWS resamples per gene.
        do i_gene = 1, n_genes
            mean_case_val = means_case(i_gene)
            mean_control_val = means_control(i_gene)

            call gather_residuals_helper(mean_case_val, sorted_case, &
                                         k_start, k_step, k_max, tau, &
                                         tmp_pool_case, n_pool_case, &
                                         max_pool_size)
            if (n_pool_case < 10) cycle

            call gather_residuals_helper(mean_control_val, sorted_control, &
                                         k_start, k_step, k_max, tau, &
                                         tmp_pool_control_own, n_pool_control_own, &
                                         max_pool_size)
            if (n_pool_control_own < 10) cycle

            observed_statistic_own_val = observed_statistic_own(i_gene)

            ! Skip genes with a non-finite observed statistic
            if (observed_statistic_own_val /= observed_statistic_own_val) cycle

            if (compute_pvalue_own(i_gene) == 1) then
                ! Build the `own` null at the level the observed statistic lives on --
                ! a difference of MEANS over n_replicates per side, where n_replicates
                ! is each side's per-gene replicate count (sorted_*%max_resid_per_gene),
                ! the count the observed means were averaged over. The observed
                ! statistic is left as-is. Both pools are already gated >= 10 above.
                if (null_method == NULL_METHOD_BLOCKED) then
                    ! GENE-BLOCKED: each null mean comes from ONE neighbour gene.
                    ! Enumerated exactly when C(2n-1, n) * n_genes_pool fits the work
                    ! arrays, which removes both the RNG and the 1/(n_boot + 1) floor;
                    ! sampled otherwise.
                    call build_blocked_means_helper(tmp_pool_case, n_pool_case, n_rep_case, &
                                                    blk_means_case, blk_weights_case, n_blk_case, &
                                                    BLOCKED_MAX_ENUM_VALUES, enum_ok_case)
                    call build_blocked_means_helper(tmp_pool_control_own, n_pool_control_own, n_rep_control, &
                                                    blk_means_control, blk_weights_control, n_blk_control, &
                                                    BLOCKED_MAX_ENUM_VALUES, enum_ok_control)
                    if (enum_ok_case .and. enum_ok_control) then
                        call compute_pvalue_blocked_exact_helper( &
                            blk_means_case, blk_weights_case, n_blk_case, &
                            blk_means_control, blk_weights_control, n_blk_control, &
                            abs(observed_statistic_own_val), &
                            blk_perm, blk_stack_left, blk_stack_right, blk_sorted, blk_cumw, &
                            pvalues_own(i_gene))
                    else
                        call compute_pvalue_blocked_bootstrap_helper( &
                            tmp_pool_case, n_pool_case, &
                            tmp_pool_control_own, n_pool_control_own, &
                            n_rep_case, n_rep_control, &
                            abs(observed_statistic_own_val), N_BOOTSTRAP_DRAWS, &
                            rbuf, pvalues_own(i_gene))
                    end if
                else
                    ! POOLED: resample n_replicates residuals iid from the whole
                    ! neighbourhood pool per side, average, difference.
                    call compute_pvalue_bootstrap_mean_helper( &
                        tmp_pool_case(1:n_pool_case), n_pool_case, &
                        tmp_pool_control_own(1:n_pool_control_own), n_pool_control_own, &
                        n_rep_case, n_rep_control, &
                        abs(observed_statistic_own_val), N_BOOTSTRAP_DRAWS, &
                        rbuf, pvalues_own(i_gene))
                end if
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
        n_genes, norm_method, k_start, k_step, k_max, tau, null_method, &
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
        integer(int32), intent(in) :: null_method
        !! Null construction: `NULL_METHOD_POOLED` (0) resamples residuals iid from
        !! the whole neighbourhood pool; `NULL_METHOD_BLOCKED` (1) draws one neighbour
        !! GENE and resamples within it, and is computed by exact enumeration wherever
        !! `C(2n-1, n) * n_genes_pool` fits `BLOCKED_MAX_ENUM_VALUES`.
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
        real(real64), allocatable :: blk_means_case(:), blk_weights_case(:)
        real(real64), allocatable :: blk_means_control(:), blk_weights_control(:)
        real(real64), allocatable :: blk_sorted(:), blk_cumw(:), rbuf(:)
        integer(int32), allocatable :: blk_perm(:), blk_stack_left(:), blk_stack_right(:)
        integer(int32) :: sort_ierr, blk_capacity, rbuf_size, per_iter_max
        integer(int32) :: k_eff_case, k_eff_control
        integer(int64) :: n_ms_case, n_ms_control, need_case, need_control

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
        M_ALLOCATE(tmp_pool_control_own(max_pool_size * 2))

        ! Enumeration capacity for the gene-blocked null: C(2n-1, n) values per gene,
        ! times the most genes a pool can hold (k_max, or fewer if max_pool_size binds
        ! first). Sized per side and allocated at the larger. If the product overflows
        ! or exceeds BLOCKED_MAX_ENUM_VALUES the arrays stay at length 1, which makes
        ! `build_blocked_means_helper` report `ok = .false.` and the per-gene loop fall
        ! back to the sampled blocked null -- no special-casing needed downstream.
        blk_capacity = 1
        if (null_method == NULL_METHOD_BLOCKED) then
            k_eff_case = min(k_max, max(1, max_pool_size / max(1, n_replicates_case)))
            k_eff_control = min(k_max, max(1, max_pool_size / max(1, n_replicates_control)))
            n_ms_case = n_multisets_helper(n_replicates_case)
            n_ms_control = n_multisets_helper(n_replicates_control)
            need_case = 0_int64
            need_control = 0_int64
            if (n_ms_case > 0_int64 .and. n_ms_case <= BLOCKED_MAX_ENUM_VALUES / int(k_eff_case, int64)) &
                need_case = n_ms_case * int(k_eff_case, int64)
            if (n_ms_control > 0_int64 .and. &
                n_ms_control <= BLOCKED_MAX_ENUM_VALUES / int(k_eff_control, int64)) &
                need_control = n_ms_control * int(k_eff_control, int64)
            blk_capacity = int(max(1_int64, max(need_case, need_control)), int32)
        end if

        M_ALLOCATE(blk_means_case(blk_capacity))
        M_ALLOCATE(blk_weights_case(blk_capacity))
        M_ALLOCATE(blk_means_control(blk_capacity))
        M_ALLOCATE(blk_weights_control(blk_capacity))
        M_ALLOCATE(blk_perm(blk_capacity))
        M_ALLOCATE(blk_stack_left(blk_capacity))
        M_ALLOCATE(blk_stack_right(blk_capacity))
        M_ALLOCATE(blk_sorted(blk_capacity))
        M_ALLOCATE(blk_cumw(blk_capacity + 1))

        ! Random-number buffer for the sampled null constructions. Sized so that one
        ! `random_number` call covers as many draws as fit, capped at RNG_BUFFER_MAX
        ! so a large replicate count cannot blow the allocation up. The blocked
        ! sampler needs one extra value per side per draw (the gene index), so size
        ! for that; the pooled sampler simply uses less of the buffer.
        per_iter_max = max(1, (n_replicates_case + 1) + (n_replicates_control + 1))
        rbuf_size = per_iter_max
        if (int(N_BOOTSTRAP_DRAWS, int64) * int(per_iter_max, int64) < int(RNG_BUFFER_MAX, int64)) then
            rbuf_size = N_BOOTSTRAP_DRAWS * per_iter_max
        else
            rbuf_size = max(per_iter_max, (RNG_BUFFER_MAX / per_iter_max) * per_iter_max)
        end if
        M_ALLOCATE(rbuf(rbuf_size))

        if (is_err(ierr)) return

        call compute_noise_pvalue_pipeline_helper( &
            sorted_case, sorted_control, &
            means_case, means_control, &
            observed_statistic_own, compute_pvalue_own, &
            n_genes, k_start, k_step, k_max, tau, null_method, &
            pvalues_own, n_genes_with_pvalue, &
            max_pool_size, &
            neighborhood_size_own_case, neighborhood_size_own_control, &
            neighborhood_size_case, &
            tmp_pool_case, tmp_pool_control_own, &
            blk_means_case, blk_weights_case, blk_means_control, blk_weights_control, &
            blk_perm, blk_stack_left, blk_stack_right, blk_sorted, blk_cumw, rbuf, &
            ierr)

    end subroutine compute_noise_pvalue_pipeline

end module noise_model

! =============================================================================
! C wrapper (outside the module, as per project convention)
! =============================================================================

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
    n_genes, norm_method, k_start, k_step, k_max, tau, null_method, &
    pvalues_own, n_genes_with_pvalue, &
    max_pool_size, &
    neighborhood_size_own_case, neighborhood_size_own_control, &
    neighborhood_size_case, &
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
    integer(c_int), intent(in), target :: null_method
    !! Null construction: 0 = iid resample from the pooled neighbourhood (historical
    !! behaviour); 1 = gene-blocked null, exactly enumerated where it fits. Callers
    !! that predate this argument should pass 0 to reproduce previous results.
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
    M_CHECK_NON_NULL(null_method)
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
        n_genes, norm_method, k_start, k_step, k_max, tau, null_method, &
        pvalues_own, n_genes_with_pvalue, &
        max_pool_size, &
        neighborhood_size_own_case, neighborhood_size_own_control, &
        neighborhood_size_case, &
        ierr)

end subroutine compute_noise_pvalues_pipeline_c