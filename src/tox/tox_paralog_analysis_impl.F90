#include <src/macros.h>

!> Paralog-subset expression patterns -- dosage effect and subfunctionalization -- relative to an ancestral ortholog.
!|
!| Candidate paralog subsets are enumerated as bitmask-encoded gene sets, built up one gene at a time
!| in index order, starting from single genes. At every extension step the candidate is scored against
!| the pattern-specific criterion (small angle plus magnitude gain for dosage effect, or bounded
!| residual distance for subfunctionalization). A candidate is pruned instead of being extended further
!| by a heuristic rule, issue #64's: for the dosage effect when its summed vector is past the angle
!| bound, for subfunctionalization when every gene after the candidate's last one (kept by the
!| filter or not) is longer than its residual. Pruning keeps the combinatorial search tractable, but it can drop a candidate whose
!| extension would have met the criterion.
!|
!| There is a routine per pattern rather than one taking a pattern flag:
!| `detect_dosage_effect` and `detect_subfunctionalization` find the subsets, and
!| `filter_paralogs_by_pattern_dosage_effect` / `_subfunctionalization` reduce a set of genes to
!| those matching. Which criterion a result was scored against is therefore visible at the call
!| site, and each routine takes only the thresholds its own pattern needs.
module tox_paralog_analysis_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH, validate_dimension_size, validate_in_range_int, &
                          clear_err_arg_pos
    use f42_math_impl, only: PI
    use f42_vector_impl, only: add_vector, subtract_vector, norm, angle_between
    use f42_bit_masks_impl, only: bitmask_set, bitmask_clear, bitmask_last, bitmask_count
    M_IMPLICIT_NONE

    integer(int32), parameter :: MODE_DOSAGE_PATTERN = 0_int32
        !! Code for detecting dosage effect in [[tox_paralog_analysis_impl(module):detect_patterns_impl(subroutine)]]
    integer(int32), parameter :: MODE_SUBFUNC_PATTERN = 1_int32
        !! Code for detecting subfunctionalization in [[tox_paralog_analysis_impl(module):detect_patterns_impl(subroutine)]]

contains

    !> summary: Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine detect_neofunctionalization_impl(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc)
        integer(int32), intent(in) :: n_axes
            !! size of `ancestors` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_families
            !! number of vectors in `ancestors`
        real(real64), dimension(n_axes, n_families), intent(in) :: ancestors
            !! RAP projected unit length expression vector of ancestral ortholog
            !! DM_MIN(-1.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), dimension(n_axes, n_genes), intent(in) :: genes
            !! RAP projected unit length expression vectors of genes
            !! DM_MIN(-1.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! M_GENE_TO_FAM_DOC(genes)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
            !! DM_SENTINEL(M_GENE_TO_FAM_SENTINEL)
        real(real64), dimension(n_axes), intent(in) :: thresholds
            !! threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
            !! DM_MIN(-1.0_real64)
            !! DM_MAX(1.0_real64)
        logical(c_bool), dimension(n_genes, n_axes), intent(out) :: neofunc
            !! `.true.` if neofunctionalization has been detected for the respective axes, always `.false.` for unassigned genes

        integer(int32) :: i_gene, i_axis, fam_idx

        neofunc = .false.

        do concurrent (i_gene = 1:n_genes) local(fam_idx) shared(gene_to_fam, n_axes, neofunc, ancestors, genes, thresholds)
            fam_idx = gene_to_fam(i_gene)
            if (fam_idx /= M_GENE_TO_FAM_SENTINEL) then
                do concurrent (i_axis = 1:n_axes) shared(fam_idx, neofunc, i_gene, ancestors, genes, thresholds)
                    neofunc(i_gene, i_axis) = abs(ancestors(i_axis, fam_idx) - genes(i_axis, i_gene)) > thresholds(i_axis)
                end do
            end if
        end do
    end subroutine detect_neofunctionalization_impl

    !> summary: Identifies subsets of paralogs matching this pattern
    !| AUTHOR_FRANZ_ERIC_SILL
    !| A subset is grown from its genes in index order: it starts from one kept gene and is extended
    !| only by kept genes after its last one. A single gene is never a result; a result has two genes
    !| at least, and is not extended further. So which subsets are found depends on the genes' order,
    !| and a result's superset is found too when the superset's own prefix did not meet the criterion
    !| (issue #244).
    pure subroutine detect_patterns_impl(ancestor, genes, n_genes, n_dims, pattern_mode, filtered_paralogs_bitmask, n_mask_chunks, n_results, max_subset_size, paralog_subsets_bitmasks, n_paralog_subsets, tmp_active_bitmask, tmp_paralog_vector, max_angle, gain_gamma, rdi_threshold, paralog_norms, paralog_norms_perm, tmp_work_array, ierr)
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes; `bitmask_n_words` computes it
            !! DM_MIN(M_BIT_MASK_N_WORDS(n_genes))
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog. Dosage effect only: it must not be all
            !! zeros, which have no direction to measure a subset's angle against, and such an
            !! ancestor is refused with `ERR_INVALID_INPUT`; subfunctionalization takes it.
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes. Dosage effect only: a subset whose summed vector has
            !! every component exactly 0.0, bit for bit, has no direction either, and is dropped
            !! from the search like one past `max_angle`. A sum that only rounds to near zero is
            !! not dropped: where its length times the ancestor's is at most 1e-12, the angle's
            !! computation refuses the whole detection with `ERR_DIVISION_BY_ZERO` (pending issue
            !! #223).
        integer(int32), intent(in) :: pattern_mode
            !! used pattern for detection
            !!
            !! |         Mode         |                                 Value                                  |          Procedure           |
            !! |----------------------|------------------------------------------------------------------------|------------------------------|
            !! |    Dosage Effect     |  [[tox_paralog_analysis_impl(module):MODE_DOSAGE_PATTERN(variable)]]   | detect_dosage_effect         |
            !! | Subfunctionalization |  [[tox_paralog_analysis_impl(module):MODE_SUBFUNC_PATTERN(variable)]]  | detect_subfunctionalization  |
            !!
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `paralog_subsets_bitmasks`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. Too large a value is capped to the
            !! maximum valid size. The bindings cap it automatically while sizing the work
            !! array; a Fortran caller caps it by calling
            !! [[tox_paralog_analysis_impl(module):calc_paralog_subsets_bitmasks_size_impl(subroutine)]] first.
            !! Zero is in range and means there is no subset to check; the sizing routine reports it
            !! when the filter kept no gene. With fewer than two kept genes there is no pair to form:
            !! the sizing then reports a work array of no subsets (and a size of 1 for one kept
            !! gene), which this routine takes, finding nothing (`n_results` 0).
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `paralog_subsets_bitmasks`. Every kept
            !! gene but the last starts a subset, so fewer slots than that are refused with
            !! `ERR_SIZE_MISMATCH` before any subset is stored, and so is a search that runs out of
            !! slots later; zero slots are enough when fewer than two genes are kept.
            !! DM_OUTPUT_FROM(work_array_size, calc_paralog_subsets_bitmasks_size, tox_paralog_analysis, AUTO)
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: paralog_subsets_bitmasks
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each column is a bit mask of the `n_genes` genes, laid out as [[f42_bit_masks_impl(module)]] describes.
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine.
            !! A bit set past `n_genes` is refused with `ERR_INVALID_INPUT`; the sizing routine
            !! ignores such bits.
        integer(int32), dimension(n_mask_chunks), intent(out) :: tmp_active_bitmask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: tmp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), intent(in), optional :: gain_gamma
            !! positive magnitude gain for dosage effect
            !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_impl, MODE_DOSAGE_PATTERN)
            !! DM_DEFAULT(0.1_real64)
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in), optional :: max_angle
            !! maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
            !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_impl, MODE_DOSAGE_PATTERN)
            !! DM_DEFAULT(4.0_real64*atan(1.0_real64))
            !! DM_MIN(0.0_real64)
            !! DM_MAX(PI)
        real(real64), dimension(n_genes), intent(in), optional :: paralog_norms
            !! euclidean norm of each gene's vector, used for subset pruning. Where the caller
            !! passes them, they are trusted to be the genes' own.
            !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_impl, MODE_SUBFUNC_PATTERN)
            !! DM_MIN(0.0_real64)
            !! DM_OUTPUT_FROM(paralog_norms, calc_paralog_norms, tox_paralog_analysis, AUTO)
        integer(int32), dimension(n_genes), intent(in), optional :: paralog_norms_perm
            !! ascending permutation of `paralog_norms`, for subset pruning: a subset is extended only
            !! while the smallest norm among all genes after its last one -- kept by the filter or
            !! not -- is at most the length of its residual. Where the caller passes it, the order
            !! is the caller's job: only the entries' range is checked, so no index is out of
            !! bounds, and a wrong order gives wrong pruning -- wrong minima for a permutation that
            !! is not ascending, and for an order with a repeated entry (such as [1, 1, 1]) minima
            !! never written, read uninitialised. The plain tier builds it itself; for the expert
            !! tier, `np.argsort(paralog_norms) + 1` in Python and `order(paralog_norms)` in R give
            !! it, counting genes from 1 as it does.
            !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_impl, MODE_SUBFUNC_PATTERN)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_genes)
        real(real64), dimension(n_genes), intent(out), optional :: tmp_work_array
            !! work array for checking the minimum value after a given index efficiently
            !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_impl, MODE_SUBFUNC_PATTERN)
        real(real64), intent(in), optional :: rdi_threshold
            !! max allowed residual distance from `ancestor`
            !! DM_REQUIRED_IF_MODE(pattern_mode, tox_paralog_analysis_impl, MODE_SUBFUNC_PATTERN)
            !! DM_MIN(0.0_real64)

        ! Locals
        integer(int32) :: i_gene, subset_size, n_active_masks, n_new_active_masks, last_filtered_paralog_idx

        call set_ok(ierr)

        ! an ancestor of zeros gives the dosage effect no direction to measure angles against
        if (pattern_mode == MODE_DOSAGE_PATTERN .and. all(ancestor == 0.0_real64)) then
            call set_err(ierr, ERR_INVALID_INPUT, arg_pos=1_int32)
            return
        end if

        if (max_subset_size == 0) then
            n_results = 0
            return
        end if

        paralog_subsets_bitmasks = 0_int32
        n_active_masks = 0_int32
        last_filtered_paralog_idx = bitmask_last(n_mask_chunks, filtered_paralogs_bitmask)

        if (last_filtered_paralog_idx > n_genes) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        ! every kept gene but the last starts a subset, and each needs its slot
        if (bitmask_count(n_mask_chunks, filtered_paralogs_bitmask) - 1 > n_paralog_subsets) then
            call set_err(ierr, ERR_SIZE_MISMATCH, arg_pos=11_int32)
            return
        end if

        ! initialize first subsets of size 1 to be extended.
        ! The subset with last gene cannot be extended
        do i_gene = 1, last_filtered_paralog_idx - 1
            if (M_BIT_MASK_TEST(filtered_paralogs_bitmask, i_gene)) then
                n_active_masks = n_active_masks + 1
                call bitmask_set(n_mask_chunks, paralog_subsets_bitmasks(:, n_active_masks), i_gene)
            end if
        end do

        n_results = 0_int32

        do subset_size = 2, max_subset_size
            n_new_active_masks = 0_int32
            do while (n_active_masks > 0)
                ! A position from either helper numbers the helper's dummy list, not this one, and
                ! names nothing this procedure's caller passed: drop it, keep the error.
                call take_active_mask_helper(paralog_subsets_bitmasks, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, tmp_active_bitmask, ierr)
                if (is_err(ierr)) then
                    call clear_err_arg_pos(ierr)
                    return
                end if

                call generate_subsets_helper(tmp_active_bitmask, filtered_paralogs_bitmask, n_mask_chunks, pattern_mode, ancestor, genes, n_genes, n_dims, tmp_paralog_vector, paralog_subsets_bitmasks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, max_angle, gain_gamma, rdi_threshold, paralog_norms, paralog_norms_perm, tmp_work_array, ierr)
                if (is_err(ierr)) then
                    call clear_err_arg_pos(ierr)
                    return
                end if
            end do

            n_active_masks = n_new_active_masks
        end do
    end subroutine detect_patterns_impl

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Generates subsets for `candidate_mask` by extending it with one valid gene.
    !| Subsets with bad gene constellation will be pruned,
    !| others will be added to the work array, either as new active subset that will be extended in coming iterations or as result.
    !|
    !| Doesn't do any input validation.
    pure subroutine generate_subsets_helper(candidate_mask, filtered_paralogs_bitmask, n_mask_chunks, pattern, ancestor, genes, n_genes, n_dims, tmp_paralog_vector, paralog_subsets_bitmasks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_paralog_norms_perm, tmp_subfunc_work_array, ierr)
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        real(real64), intent(in), optional :: subfunc_rdi_threshold
            !! max allowed residual distance from `ancestor`
        integer(int32), intent(in) :: pattern
            !! used pattern for detection
            !!
            !! |         Mode         |                                  Value                                   |
            !! |----------------------|--------------------------------------------------------------------------|
            !! |    Dosage Effect     |   [[tox_paralog_analysis_impl(module):MODE_DOSAGE_PATTERN(variable)]]    |
            !! | Subfunctionalization |   [[tox_paralog_analysis_impl(module):MODE_SUBFUNC_PATTERN(variable)]]   |
            !!
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `paralog_subsets_bitmasks`. ***USE `calc_paralog_subsets_bitmasks_size_impl` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(inout) :: paralog_subsets_bitmasks
            !! working array to hold bitmask encoded subsets for detection; the results and active
            !! subsets already in it are kept
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(int32), dimension(n_mask_chunks), intent(inout) :: candidate_mask
            !! working array to hold a subset that is a potential result candidate
        real(real64), dimension(n_dims), intent(inout) :: tmp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(inout) :: n_results
            !! number of results in `paralog_subsets_bitmasks`
        integer(int32), intent(in) :: n_active_masks
            !! number of active subsets in `paralog_subsets_bitmasks`
        integer(int32), intent(inout) :: n_new_active_masks
            !! number of new active subsets in `paralog_subsets_bitmasks`
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), intent(in), optional :: dosage_gain_gamma
            !! positive magnitude gain for dosage, default 0.1
        real(real64), intent(in), optional :: dosage_max_angle
            !! maximum angle in radians that a subset candidate must not exceed, otherwise pruned, default is Pi
        real(real64), dimension(n_genes), intent(in), optional :: subfunc_paralog_norms
            !! euclidean norms of the genes, used for subset pruning (`norm` from `f42_vector_impl` computes them).
            !! The routine trusts them to be the genes' own.
        integer(int32), dimension(n_genes), intent(in), optional :: subfunc_paralog_norms_perm
            !! ascending permutation of `subfunc_paralog_norms`, for subset pruning: a subset is extended only
            !! while the smallest norm among all genes after its last one -- kept by the filter or
            !! not -- is at most the length of its residual. Only the entries' range is checked: a
            !! permutation that is not ascending gives wrong minima, and an order with a repeated
            !! entry (such as [1, 1, 1]) leaves part of the work array of minima unwritten, so the
            !! pruning reads uninitialised values.
        real(real64), dimension(n_genes), intent(out), optional :: tmp_subfunc_work_array
            !! work array for checking the minimum value after a given index efficiently

        integer(int32) :: i_gene

        call set_ok(ierr)

        select case (pattern)
        case (MODE_DOSAGE_PATTERN)
            block
                use f42_math_impl, only: PI
                real(real64) :: subset_angle, gain, max_angle

                M_DEFAULT_VAL(dosage_gain_gamma, gain, 0.1_real64)
                M_DEFAULT_VAL(dosage_max_angle, max_angle, PI)

                !! prepare sum vector, so the extending gene just needs to be included in one operation and excluded after calculation
                !TODO optimize: this rebuilds the subset's sum vector from scratch by scanning all n_genes on every call to generate_subsets_helper (once per active mask taken from the work array), instead of carrying the running sum forward from the parent subset that already had it computed. For large gene counts/subset counts this recomputation dominates the runtime of the whole subset-extension search.
                tmp_paralog_vector = 0
                do i_gene = 1, n_genes
                    if (M_BIT_MASK_TEST(candidate_mask, i_gene)) then
                        call add_vector(tmp_paralog_vector, genes(:, i_gene))
                    end if
                end do

                ! generate extended subsets by adding succeeding genes of the last active gene if suitable.
                do i_gene = bitmask_last(n_mask_chunks, candidate_mask) + 1_int32, n_genes
                    if (M_BIT_MASK_TEST(filtered_paralogs_bitmask, i_gene)) then
                        ! extend subset by current gene
                        call bitmask_set(n_mask_chunks, candidate_mask, i_gene)

                        ! compute sum vector of all subset's genes
                        call add_vector(tmp_paralog_vector, genes(:, i_gene))

                        ! a sum of exactly zero has no direction: dropped like a sum past max_angle
                        if (all(tmp_paralog_vector == 0.0_real64)) then
                            subset_angle = huge(subset_angle)
                        else
                            call angle_between(tmp_paralog_vector, ancestor, n_dims, subset_angle, ierr)
                            if (is_err(ierr)) return
                        end if

                        ! If angle of the subset vector is close enough, there may be dosage effect
                        if (subset_angle <= max_angle) then
                            ! If norm exceeds ancestor's norm significantly, subset is a result
                            if (norm(tmp_paralog_vector) >= (1 + gain)*norm(ancestor)) then
                                call add_to_results_helper(paralog_subsets_bitmasks, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            else
                                call add_new_active_mask_helper(paralog_subsets_bitmasks, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            end if
                            if (is_err(ierr)) return
                        end if

                        ! revert extension with current gene to efficiently reuse the variables for next gene
                        call subtract_vector(tmp_paralog_vector, genes(:, i_gene))
                        call bitmask_clear(n_mask_chunks, candidate_mask, i_gene)
                    end if
                end do
            end block
        case (MODE_SUBFUNC_PATTERN)
            block
                real(real64) :: residual_norm

                if (.not. (present(subfunc_paralog_norms) .and. present(subfunc_paralog_norms_perm) .and. present(tmp_subfunc_work_array) .and. present(subfunc_rdi_threshold))) then
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
                end if

                !! initialize work array with min values, so each index i holds the min value in subarray subfunc_paralog_norms(i:n_genes)
                call fill_array_with_minvals_for_each_idx(tmp_subfunc_work_array, subfunc_paralog_norms, subfunc_paralog_norms_perm, n_genes, ierr)
                if (is_err(ierr)) return

                !! also, prepare residual, so the extending gene just needs to be included in one operation and excluded after calculation
                tmp_paralog_vector = ancestor
                do i_gene = 1, n_genes
                    if (M_BIT_MASK_TEST(candidate_mask, i_gene)) then
                        call subtract_vector(tmp_paralog_vector, genes(:, i_gene))
                    end if
                end do

                ! generate extended subsets by adding succeeding genes of the last active gene if suitable.
                do i_gene = bitmask_last(n_mask_chunks, candidate_mask) + 1_int32, n_genes
                    if (M_BIT_MASK_TEST(filtered_paralogs_bitmask, i_gene)) then
                        ! extend subset by current gene
                        call bitmask_set(n_mask_chunks, candidate_mask, i_gene)

                        ! compute residual of current subset
                        call subtract_vector(tmp_paralog_vector, genes(:, i_gene))

                        residual_norm = norm(tmp_paralog_vector)
                        if (residual_norm <= subfunc_rdi_threshold) then
                            call add_to_results_helper(paralog_subsets_bitmasks, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                        else if (i_gene < n_genes) then
                            ! tmp_subfunc_work_array(i_gene+1) is min(norm for i in i_gene+1:n_genes )
                            ! so if the minimum norm of the remaining genes is not lower the residual, prune this subset branch
                            if (tmp_subfunc_work_array(i_gene + 1) <= residual_norm) then
                                call add_new_active_mask_helper(paralog_subsets_bitmasks, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            end if
                        end if
                        if (is_err(ierr)) return

                        ! revert extension with current gene to efficiently reuse the variables for next gene
                        call add_vector(tmp_paralog_vector, genes(:, i_gene))
                        call bitmask_clear(n_mask_chunks, candidate_mask, i_gene)
                    end if
                end do
            end block
        case default
            call set_err(ierr, ERR_INVALID_INPUT, 4_int32)
            return
        end select

    end subroutine generate_subsets_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper for subfunctionalization pruning. It initializes a working array with each index i holding the min value in subarray src_arr(i:src_arr_len)
    pure subroutine fill_array_with_minvals_for_each_idx(out_arr, src_arr, sorted_src_arr_perm, src_arr_len, ierr)
        integer(int32), intent(in) :: src_arr_len
            !! number elements in `src_arr`
        real(real64), dimension(src_arr_len), intent(in), optional :: src_arr
            !! array that holds the original values
        integer(int32), dimension(src_arr_len), intent(in), optional :: sorted_src_arr_perm
            !! sorted permutation vector, so each index holds the index of the value in `src_arr` as if `src_arr` was soted ascending
        real(real64), dimension(src_arr_len), intent(out), optional :: out_arr
            !! output array, e.g. source [50, 75, 0, 100, 25] would lead to `out_arr` [0, 0, 0, 25, 25]
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_out, i_perm, last_min_index, current_min_idx

        call set_ok(ierr)

        last_min_index = 0
        do i_perm = 1, src_arr_len
            ! Take next higher value index
            current_min_idx = sorted_src_arr_perm(i_perm)
            if (current_min_idx > src_arr_len) then
                call set_err(ierr, ERR_INVALID_INPUT, arg_pos=3_int32)
                exit
            end if

            ! If the current value comes after the previous iterated values y in original array,
            ! fill the indices after the already filled ones upto the current value's index with current value.
            ! Example: sorted_perm=[2, 1, 3]
            !   1. iteration fills indices 1-2
            !   2. iteration skips because value of index 2 is lower than of index 1
            !   3. iteration fills to the end
            if (current_min_idx > last_min_index) then
                do i_out = last_min_index + 1, current_min_idx
                    out_arr(i_out) = src_arr(current_min_idx)
                end do

                if (current_min_idx == src_arr_len) exit
                last_min_index = current_min_idx
            end if
        end do
    end subroutine fill_array_with_minvals_for_each_idx

    !> AUTHOR_FRANZ_ERIC_SILL
    !| For memory efficiency this subroutine helps holding different kinds of masks in a single array.
    !| To achieve this, `subsets` has this structure: [...results, ...active_masks, ...new_active_masks]
    !| This routine removes an active mask from the `subsets` array and returns it in `active_mask`.
    pure subroutine take_active_mask_helper(subsets, n_mask_chunks, n_subsets, n_results, n_active_masks, n_new_active_masks, active_mask, ierr)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks in a mask
        integer(int32), intent(in) :: n_subsets
            !! number of subsets that can be hold
        integer(int32), intent(in) :: n_results
            !! number of results in `subsets`
        integer(int32), intent(inout) :: n_active_masks
            !! number of active masks in `subsets`
        integer(int32), intent(in) :: n_new_active_masks
            !! number of new active masks in `subsets`
        integer(int32), dimension(n_mask_chunks, n_subsets), intent(inout) :: subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! taken active mask from `subsets`, will be remove from `subsets`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_mask_chunks, ierr)
        call validate_dimension_size(n_subsets, ierr)
        call validate_in_range_int(n_active_masks, ierr, min=1_int32)
        call validate_in_range_int(n_results, ierr, min=0_int32)
        call validate_in_range_int(n_new_active_masks, ierr, min=0_int32)
        if (is_err(ierr)) return

        if (n_subsets < n_results + n_active_masks + n_new_active_masks) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        ! take one active mask, always the last of all available actives
        active_mask = subsets(:, n_results + n_active_masks)

        ! in `subsets`, replace taken mask by last new active mask
        subsets(:, n_results + n_active_masks) = subsets(:, n_results + n_active_masks + n_new_active_masks)

        ! update count
        n_active_masks = n_active_masks - 1
    end subroutine take_active_mask_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| For memory efficiency this subroutine helps holding different kinds of masks in a single array.
    !| To achieve this, `subsets` has this structure: [...results, ...active_masks, ...new_active_masks]
    !| This routine adds a result mask to `subsets` array.
    pure subroutine add_to_results_helper(subsets, n_mask_chunks, n_subsets, n_results, n_active_masks, n_new_active_masks, result, ierr)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks in a mask
        integer(int32), intent(in) :: n_subsets
            !! number of subsets that can be hold
        integer(int32), intent(inout) :: n_results
            !! number of results in `subsets`
        integer(int32), intent(in) :: n_active_masks
            !! number of active masks in `subsets`
        integer(int32), intent(in) :: n_new_active_masks
            !! number of new active masks in `subsets`
        integer(int32), dimension(n_mask_chunks, n_subsets), intent(inout) :: subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: result
            !! result to add to `subsets`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_mask_chunks, ierr)
        call validate_dimension_size(n_subsets, ierr)
        call validate_in_range_int(n_active_masks, ierr, min=0_int32)
        call validate_in_range_int(n_results, ierr, min=0_int32)
        call validate_in_range_int(n_new_active_masks, ierr, min=0_int32)
        if (is_err(ierr)) return

        if (n_subsets < n_results + n_active_masks + n_new_active_masks + 1) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        ! free the index after results to hold a new result
        ! Thus, move first new active mask to end
        subsets(:, n_results + n_active_masks + n_new_active_masks + 1) = subsets(:, n_results + n_active_masks + 1)
        ! move first active mask before new active masks
        subsets(:, n_results + n_active_masks + 1) = subsets(:, n_results + 1)
        ! store result
        n_results = n_results + 1
        subsets(:, n_results) = result
    end subroutine add_to_results_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| For memory efficiency this subroutine helps holding different kinds of masks in a single array.
    !| To achieve this, `subsets` has this structure: [...results, ...active_masks, ...new_active_masks]
    !| This routine adds a new active mask to the `subsets` array.
    pure subroutine add_new_active_mask_helper(subsets, n_mask_chunks, n_subsets, n_results, n_active_masks, n_new_active_masks, new_active_mask, ierr)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks in a mask
        integer(int32), intent(in) :: n_subsets
            !! number of subsets that can be hold
        integer(int32), intent(in) :: n_results
            !! number of results in `subsets`
        integer(int32), intent(in) :: n_active_masks
            !! number of active masks in `subsets`
        integer(int32), intent(inout) :: n_new_active_masks
            !! number of new active masks in `subsets`
        integer(int32), dimension(n_mask_chunks, n_subsets), intent(inout) :: subsets
            !! working array to hold bitmask encoded subsets for detection; the masks already in it are kept
        integer(int32), dimension(n_mask_chunks), intent(in) :: new_active_mask
            !! new active mask to add to `subsets`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_mask_chunks, ierr)
        call validate_dimension_size(n_subsets, ierr)
        call validate_in_range_int(n_active_masks, ierr, min=0_int32)
        call validate_in_range_int(n_results, ierr, min=0_int32)
        call validate_in_range_int(n_new_active_masks, ierr, min=0_int32)
        if (is_err(ierr)) return

        if (n_subsets < n_results + n_active_masks + n_new_active_masks + 1) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        ! simply append the new active mask to the end
        n_new_active_masks = n_new_active_masks + 1
        subsets(:, n_results + n_active_masks + n_new_active_masks) = new_active_mask
    end subroutine add_new_active_mask_helper

    !> summary: Prefilters the genes for a pattern, so genes that cannot match it are not tried as subset extensions
    !| AUTHOR_FRANZ_ERIC_SILL
    !| This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    pure subroutine filter_paralogs_by_pattern_impl(pattern_mode, gene_angles, threshold, n_genes, n_families, gene_to_fam, family_bitmasks, n_mask_chunks)
        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(in) :: n_families
            !! number of families
            !! DM_MIN(1_int32)
            !! DM_MAX(n_genes)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! DM_MIN(M_BIT_MASK_N_WORDS(n_genes))
        integer(int32), intent(in) :: pattern_mode
            !! used pattern for detection
            !!
            !! |         Mode         |                                Value                                |                  Procedure                     |
            !! |----------------------|---------------------------------------------------------------------|------------------------------------------------|
            !! |    Dosage Effect     |  [[tox_paralog_analysis_impl(module):MODE_DOSAGE_PATTERN(variable)]]   | filter_paralogs_by_pattern_dosage_effect       |
            !! | Subfunctionalization |  [[tox_paralog_analysis_impl(module):MODE_SUBFUNC_PATTERN(variable)]]  | filter_paralogs_by_pattern_subfunctionalization |
            !!
        real(real64), dimension(n_genes), intent(in) :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(PI)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_families)
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: family_bitmasks
            !! bit mask that will have the indices of genes kept by this pattern set to 1, else 0

        integer(int32) :: i_gene, family_idx

        family_bitmasks = 0_int32

        select case (pattern_mode)
        case (MODE_DOSAGE_PATTERN)
            ! only genes with angles at most the threshold (e.g. the gene-family median or lower five percentile) are marked active
            do i_gene = 1, n_genes
                if (gene_angles(i_gene) <= threshold) then
                    family_idx = gene_to_fam(i_gene)
                    call bitmask_set(n_mask_chunks, family_bitmasks(:, family_idx), i_gene)
                end if
            end do
        case (MODE_SUBFUNC_PATTERN)
            ! only genes with angles at least the threshold (e.g. the gene-family median angle) are marked active
            do i_gene = 1, n_genes
                if (gene_angles(i_gene) >= threshold) then
                    family_idx = gene_to_fam(i_gene)
                    call bitmask_set(n_mask_chunks, family_bitmasks(:, family_idx), i_gene)
                end if
            end do
        end select
    end subroutine filter_paralogs_by_pattern_impl

    !> summary: The euclidean norm of each gene's vector, which the subfunctionalization prunes by
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Each norm is f42's scaled `norm`, so no vector of finite components overflows or underflows
    !| on the way: the result is +Inf only where the true norm exceeds `huge`.
    pure subroutine calc_paralog_norms_impl(genes, n_dims, n_genes, paralog_norms)
        integer(int32), intent(in) :: n_dims
            !! size of the vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        real(real64), dimension(n_genes), intent(out) :: paralog_norms
            !! euclidean norm of each gene's vector

        integer(int32) :: i_gene

        do i_gene = 1, n_genes
            paralog_norms(i_gene) = norm(genes(:, i_gene))
        end do
    end subroutine calc_paralog_norms_impl

    !> summary: Calculates the needed size for the paralog-subsets work array
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The `detect_*` subroutines need a work array for the subsets they try. It holds, side by side,
    !| the results found so far, the subsets still to be grown and the ones grown from them, so it is
    !| sized for the worst case: every subset of the first F - 1 kept genes stays a candidate, and every
    !| subset that adds the last kept gene becomes a result. With F genes kept and subsets of at most m
    !| genes, that is one subset of the first F - 1 genes for each size 1 to m: the sum over
    !| k = 1 .. min(m, F - 1) of C(F - 1, k), which is 2**(F - 1) - 1 for m = F, and 0 for fewer
    !| than two kept genes.
    pure subroutine calc_paralog_subsets_bitmasks_size_impl(max_subset_size, n_genes, work_array_size, filtered_paralogs_bitmask, n_mask_chunks)
        integer(int32), intent(in) :: n_genes
            !! number of genes. It is checked before `n_mask_chunks`, so a mask of too few words, none
            !! included, is reported as `n_genes` out of range.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_mask_chunks * 32)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        integer(int32), intent(inout) :: max_subset_size
            !! maximum size that a subset must not exceed. Zero is in range and means there is
            !! nothing to size a work array for, which is reported back as a size of zero.
            !! DM_MIN(0_int32)
            !! @warning
            !! It comes back lowered to the number of kept genes if it is larger, and lowered further
            !! to the largest size whose work array holds at most 2**31 - 1 (2,147,483,647) subsets, the
            !! maximum valid size, where the size asked for would need more.
            !! @endwarning
        integer(int32), intent(out) :: work_array_size
            !! The number of subsets the work array must hold in the worst case, for the
            !! `max_subset_size` handed back
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask of the genes the filter kept. Only the first `n_genes` bits are counted;
            !! a bit past them is ignored here, while the detections refuse it.

        integer(int32) :: i_gene, subset_size, n_genes_filtered
        integer(int64) :: n_subsets, binomial

        if (max_subset_size == 0) then
            work_array_size = 0
            return
        end if

        n_genes_filtered = 0
        do concurrent (i_gene = 1:n_genes) shared(filtered_paralogs_bitmask) reduce(+:n_genes_filtered)
            if (M_BIT_MASK_TEST(filtered_paralogs_bitmask, i_gene)) then
                n_genes_filtered = n_genes_filtered + 1
            end if
        end do
        max_subset_size = min(n_genes_filtered, max_subset_size)

        ! Sum C(F - 1, k) for k = 1 .. min(m, F - 1), in int64: each binomial is the previous one
        ! times (F - k)/k, exactly, and the product stays below 2**62 as long as the sum fits int32.
        n_subsets = 0_int64
        binomial = 1_int64
        do subset_size = 1, min(max_subset_size, n_genes_filtered - 1)
            binomial = binomial*int(n_genes_filtered - subset_size, int64)/int(subset_size, int64)
            if (n_subsets + binomial > int(huge(0_int32), int64)) then
                max_subset_size = subset_size - 1
                exit
            end if
            n_subsets = n_subsets + binomial
        end do
        work_array_size = int(n_subsets, int32)
    end subroutine calc_paralog_subsets_bitmasks_size_impl
end module tox_paralog_analysis_impl
