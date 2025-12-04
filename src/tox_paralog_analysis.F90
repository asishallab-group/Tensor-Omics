#include "macros.h"

module tox_paralog_analysis
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH, validate_dimension_size, validate_in_range_int, validate_all_in_range_int, validate_in_range_real, validate_all_in_range_real
    use f42_utils, only: add_vector, subtract_vector, norm, angle_between, above
    implicit none

    integer(int32), parameter :: DOSAGE_PATTERN = 0
    integer(int32), parameter :: SUBFUNC_PATTERN = 1

contains

    !> Identifies neofunctionalization for genes by checking whether the difference of expression to its ancestor exceeds the threshold for the respective axis.
    pure subroutine detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
        integer(int32), intent(in) :: n_axes
            !! size of `ancestors` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_families
            !! number of vectors in `ancestors`
        real(real64), dimension(n_axes, n_families), intent(in) :: ancestors
            !! RAP projected unit length expression vector of ancestral ortholog
        real(real64), dimension(n_axes, n_genes), intent(in) :: genes
            !! RAP projected unit length expression vectors of genes
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! mapping of gene index to family index
        real(real64), dimension(n_axes), intent(in) :: thresholds
            !! threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
        logical, dimension(n_genes, n_axes), intent(out) :: neofunc
            !! `.true.` if neofunctionalization has been detected for the respective axes
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_gene, i_axis, fam_idx

        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_axes, ierr)
        call validate_dimension_size(n_families, ierr)
        call validate_all_in_range_real(ancestors, n_axes * n_families, ierr, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_real(genes, n_axes * n_genes, ierr, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, min=1_int32, max=n_families)
        call validate_all_in_range_real(thresholds, n_axes, ierr, min=-1.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        do i_gene = 1, n_genes
            fam_idx = gene_to_fam(i_gene)
            do i_axis = 1, n_axes
                neofunc(i_gene, i_axis) = abs(ancestors(i_axis, fam_idx) - genes(i_axis, i_gene)) > thresholds(i_axis)
            end do
        end do
    end subroutine detect_neofunctionalization

    !> Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)
    pure subroutine detect_dosage_effect(ancestor, genes, n_genes, n_dims, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, ierr, max_angle, gain_gamma)
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code
        real(real64), intent(in), optional :: gain_gamma
            !! positive magnitude gain for dosage effect, default 0.1
        real(real64), intent(in), optional :: max_angle
            !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi

        call detect_patterns(ancestor, genes, n_genes, n_dims, DOSAGE_PATTERN, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, dosage_max_angle=max_angle, dosage_gain_gamma=gain_gamma, ierr=ierr)
    end subroutine detect_dosage_effect

    !> Identifies subsets of paralogs exhibiting significant angles to the `ancestor`
    pure subroutine detect_subfunctionalization(ancestor, genes, n_genes, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        real(real64), intent(in) :: rdi_threshold
            !! max allowed residual distance from `ancestor`
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code
        real(real64), dimension(n_genes), intent(in) :: paralog_norms
            !! needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
        integer(int32), dimension(n_genes), intent(in) :: sorted_paralog_norms_perm
            !! needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
        real(real64), dimension(n_genes), intent(out) :: temp_work_array
            !! needed for efficient check of minimum value after a certain index

        call detect_patterns(ancestor, genes, n_genes, n_dims, SUBFUNC_PATTERN, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, subfunc_rdi_threshold=rdi_threshold, subfunc_paralog_norms=paralog_norms, subfunc_sorted_paralog_norms_perm=sorted_paralog_norms_perm, subfunc_temp_work_array=temp_work_array, ierr=ierr)
    end subroutine detect_subfunctionalization

    !> Identifies subsets of paralogs where dosage effect or subfunctionalization applies, depending on `pattern`
    pure subroutine detect_patterns(ancestor, genes, n_genes, n_dims, pattern, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr)
        use f42_utils, only: PI

        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        integer(int32), intent(in) :: pattern
            !! used pattern for detection
            !!
            !! |       Pattern        | Value |
            !! |----------------------|-------|
            !! |    Dosage Effect     |   0   |
            !! | Subfunctionalization |   1   |
            !!
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code
        real(real64), intent(in), optional :: dosage_gain_gamma
            !! in dosage mode required positive magnitude gain for dosage, default 0.1
        real(real64), intent(in), optional :: dosage_max_angle
            !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
        real(real64), dimension(n_genes), intent(in), optional :: subfunc_paralog_norms
            !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` from `f42_utils` function for this)
        integer(int32), dimension(n_genes), intent(in), optional :: subfunc_sorted_paralog_norms_perm
            !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
        real(real64), dimension(n_genes), intent(out), optional :: subfunc_temp_work_array
            !! in subfunctionalization mode needed for efficient check of minimum value after a certain index
        real(real64), intent(in), optional :: subfunc_rdi_threshold
            !! max allowed residual distance from `ancestor`

        ! Locals
        integer(int32) :: i_gene, subset_size, n_active_masks, n_new_active_masks, last_filtered_paralog_idx

        call set_ok(ierr)

        if (max_subset_size == 0) then
            n_results = 0
            return
        end if

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_paralog_subsets, ierr)
        call validate_dimension_size(n_dims, ierr)
        call validate_dimension_size(n_mask_chunks, ierr)
        call validate_in_range_int(max_subset_size, ierr, min=1_int32)
        call validate_all_in_range_real(ancestor, n_dims, ierr)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr)
        call validate_in_range_real(dosage_gain_gamma, ierr, min=above(0.0_real64))
        call validate_in_range_real(dosage_max_angle, ierr, min=0.0_real64, max=PI)
        call validate_all_in_range_real(subfunc_paralog_norms, n_genes, ierr, min=0.0_real64)
        call validate_all_in_range_int(subfunc_sorted_paralog_norms_perm, n_genes, ierr, min=1_int32, max=n_genes)
        call validate_in_range_real(subfunc_rdi_threshold, ierr, min=0.0_real64)
        if (n_mask_chunks * 32 < n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        work_arr_paralog_subsets = 0_int32
        n_active_masks = 0_int32
        last_filtered_paralog_idx = mask_get_first_successor_idx(filtered_paralogs_mask) - 1

        if (last_filtered_paralog_idx > n_genes) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        ! initialize first subsets of size 1 to be extended.
        ! The subset with last gene cannot be extended
        do i_gene = 1, last_filtered_paralog_idx - 1
            if (mask_check_state(filtered_paralogs_mask, i_gene)) then
                n_active_masks = n_active_masks + 1
                call mask_set_state(work_arr_paralog_subsets(:, n_active_masks), i_gene, .true., ierr)
                if (is_err(ierr)) return
            end if
        end do

        n_results = 0_int32

        do subset_size = 2, max_subset_size
            n_new_active_masks = 0_int32
            do while (n_active_masks > 0)
                call take_active_mask_helper(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, active_mask, ierr)
                if (is_err(ierr)) return

                call generate_subsets_helper(active_mask, filtered_paralogs_mask, n_mask_chunks, pattern, ancestor, genes, n_genes, n_dims, temp_paralog_vector, work_arr_paralog_subsets, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr)
                if (is_err(ierr)) return
            end do

            n_active_masks = n_new_active_masks
        end do
    end subroutine detect_patterns

    !> Generates subsets for `candidate_mask` by extending it with one valid gene.
    !| Subsets with bad gene constellation will be pruned,
    !| others will be added to the work array, either as new active subset that will be extended in coming iterations or as result.
    !| 
    !| Doesn't do any input validation.
    pure subroutine generate_subsets_helper(candidate_mask, filtered_paralogs_mask, n_mask_chunks, pattern, ancestor, genes, n_genes, n_dims, temp_paralog_vector, work_arr_paralog_subsets, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr)
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
            !! |       Pattern        | Value |
            !! |----------------------|-------|
            !! |    Dosage Effect     |   0   |
            !! | Subfunctionalization |   1   |
            !!
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(int32), dimension(n_mask_chunks), intent(inout) :: candidate_mask
            !! working array to hold a subset that is a potential result candidate
        real(real64), dimension(n_dims), intent(inout) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(inout) :: n_results
            !! number of results in `work_arr_paralog_subsets`
        integer(int32), intent(in) :: n_active_masks
            !! number of active subsets in `work_arr_paralog_subsets`
        integer(int32), intent(inout) :: n_new_active_masks
            !! number of new active subsets in `work_arr_paralog_subsets`
        integer(int32), intent(out) :: ierr
            !! error code
        real(real64), intent(in), optional :: dosage_gain_gamma
            !! in dosage mode required positive magnitude gain for dosage, default 0.1
        real(real64), intent(in), optional :: dosage_max_angle
            !! in dosage mode maximum angle in radians that a subset candidate must not exceed, otherwise pruned, default is Pi
        real(real64), dimension(n_genes), intent(in), optional :: subfunc_paralog_norms
            !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` from `f42_utils` function for this)
        integer(int32), dimension(n_genes), intent(in), optional :: subfunc_sorted_paralog_norms_perm
            !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
        real(real64), dimension(n_genes), intent(out), optional :: subfunc_temp_work_array
            !! in subfunctionalization mode needed for efficient check of minimum value after a certain index

        integer(int32) :: i_gene

        call set_ok(ierr)

        select case (pattern)
        case (DOSAGE_PATTERN)
            block
                use f42_utils, only: PI
                real(real64) :: subset_angle, gain, max_angle

                M_DEFAULT_VAL(dosage_gain_gamma, gain, 0.1_real64)
                M_DEFAULT_VAL(dosage_max_angle, max_angle, PI)

                !! prepare sum vector, so the extending gene just needs to be included in one operation and excluded after calculation
                temp_paralog_vector = 0
                do i_gene = 1, n_genes
                    if (mask_check_state(candidate_mask, i_gene)) then
                        call add_vector(temp_paralog_vector, genes(:, i_gene))
                    end if
                end do

                ! generate extended subsets by adding succeeding genes of the last active gene if suitable.
                do i_gene = mask_get_first_successor_idx(candidate_mask), n_genes
                    if (mask_check_state(filtered_paralogs_mask, i_gene)) then
                        ! extend subset by current gene
                        call mask_set_state(candidate_mask, i_gene, .true., ierr)
                        if (is_err(ierr)) return

                        ! compute sum vector of all subset's genes
                        call add_vector(temp_paralog_vector, genes(:, i_gene))

                        call angle_between(temp_paralog_vector, ancestor, n_dims, subset_angle, ierr)
                        if (is_err(ierr)) return

                        ! If angle of the subset vector is close enough, there may be dosage effect
                        if (subset_angle <= max_angle) then
                            ! If norm exceeds ancestor's norm significantly, subset is a result
                            if (norm(temp_paralog_vector) >= (1 + gain) * norm(ancestor)) then
                                call add_to_results_helper(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            else
                                call add_new_active_mask_helper(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            end if
                            if (is_err(ierr)) return
                        end if

                        ! revert extension with current gene to efficiently reuse the variables for next gene
                        call subtract_vector(temp_paralog_vector, genes(:, i_gene))
                        call mask_set_state(candidate_mask, i_gene, .false., ierr)
                        if (is_err(ierr)) return
                    end if
                end do
            end block
        case (SUBFUNC_PATTERN)
            block
                real(real64) :: residual_norm

                if (.not. (present(subfunc_paralog_norms) .and. present(subfunc_sorted_paralog_norms_perm) .and. present(subfunc_temp_work_array) .and. present(subfunc_rdi_threshold))) then
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
                end if

                !! initialize work array with min values, so each index i holds the min value in subarray subfunc_paralog_norms(i:n_genes)
                call fill_array_with_minvals_for_each_idx(subfunc_temp_work_array, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, n_genes, ierr)
                if (is_err(ierr)) return

                !! also, prepare residual, so the extending gene just needs to be included in one operation and excluded after calculation
                temp_paralog_vector = ancestor
                do i_gene = 1, n_genes
                    if (mask_check_state(candidate_mask, i_gene)) then
                        call subtract_vector(temp_paralog_vector, genes(:, i_gene))
                    end if
                end do

                ! generate extended subsets by adding succeeding genes of the last active gene if suitable.
                do i_gene = mask_get_first_successor_idx(candidate_mask), n_genes
                    if (mask_check_state(filtered_paralogs_mask, i_gene)) then
                        ! extend subset by current gene
                        call mask_set_state(candidate_mask, i_gene, .true., ierr)
                        if (is_err(ierr)) return

                        ! compute residual of current subset
                        call subtract_vector(temp_paralog_vector, genes(:, i_gene))

                        residual_norm = norm(temp_paralog_vector)
                        if (residual_norm <= subfunc_rdi_threshold) then
                            call add_to_results_helper(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                        else if (i_gene < n_genes) then
                            ! subfunc_temp_work_array(i_gene+1) is min(norm for i in i_gene+1:n_genes )
                            ! so if the minimum norm of the remaining genes is not lower the residual, prune this subset branch
                            if (subfunc_temp_work_array(i_gene + 1) <= residual_norm) then
                                call add_new_active_mask_helper(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            end if
                        end if
                        if (is_err(ierr)) return

                        ! revert extension with current gene to efficiently reuse the variables for next gene
                        call add_vector(temp_paralog_vector, genes(:, i_gene))
                        call mask_set_state(candidate_mask, i_gene, .false., ierr)
                        if (is_err(ierr)) return
                    end if
                end do
            end block
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select

    end subroutine generate_subsets_helper

    !> Helper for subfunctionalization pruning. It initializes a working array with each index i holding the min value in subarray src_arr(i:src_arr_len)
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
            !! error code

        integer(int32) :: i_out, i_perm, last_min_index, current_min_idx

        call set_ok(ierr)

        last_min_index = 0
        do i_perm = 1, src_arr_len
            ! Take next higher value index
            current_min_idx = sorted_src_arr_perm(i_perm)
            if (current_min_idx > src_arr_len) then
                call set_err(ierr, ERR_INVALID_INPUT)
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

    !> For memory efficiency this subroutine helps holding different kinds of masks in a single array.
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
            !! error code

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

    !> For memory efficiency this subroutine helps holding different kinds of masks in a single array.
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
            !! error code

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

    !> For memory efficiency this subroutine helps holding different kinds of masks in a single array.
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
        integer(int32), dimension(n_mask_chunks, n_subsets), intent(out) :: subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: new_active_mask
            !! new active mask to add to `subsets`
        integer(int32), intent(out) :: ierr
            !! error code

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

    !> This subroutine easily determines the needed chunk count for subset bit masks, as an integer has only 32 bits.
    pure subroutine mask_chunk_count(n_genes, count)
        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(out) :: count
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes

        !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
        count = (n_genes + 31) / 32
    end subroutine mask_chunk_count

    !> This subroutine prefilters the genes for subfunctionalization,
    !| as genes that are already too close in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    pure subroutine filter_paralogs_by_pattern_subfunctionalization(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(in) :: n_families
            !! number of families
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        real(real64), dimension(n_genes), intent(in) :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: masks
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! error code

        call filter_paralogs_by_pattern(SUBFUNC_PATTERN, gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
    end subroutine filter_paralogs_by_pattern_subfunctionalization

    !> This subroutine prefilters the genes for dosage effect,
    !| as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    pure subroutine filter_paralogs_by_pattern_dosage_effect(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(in) :: n_families
            !! number of families
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        real(real64), dimension(n_genes), intent(in) :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: masks
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! error code

        call filter_paralogs_by_pattern(DOSAGE_PATTERN, gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
    end subroutine filter_paralogs_by_pattern_dosage_effect

    !> This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    pure subroutine filter_paralogs_by_pattern(pattern, gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
        use f42_utils, only: PI

        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(in) :: n_families
            !! number of families
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        integer(int32), intent(in) :: pattern
            !! used pattern for detection
            !!
            !! |       Pattern        | Value |
            !! |----------------------|-------|
            !! |    Dosage Effect     |   0   |
            !! | Subfunctionalization |   1   |
            !!
        real(real64), dimension(n_genes), intent(in) :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: masks
            !! bit mask that will have indices of genes kept by pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_gene, family_idx

        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_mask_chunks, ierr)
        call validate_in_range_real(threshold, ierr)
        call validate_all_in_range_real(gene_angles, n_genes, ierr, min=0.0_real64, max=PI)
        call validate_in_range_int(n_families, ierr, min=1_int32, max=n_genes)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, min=1_int32, max=n_families)
        if (n_mask_chunks * 32 < n_genes) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        masks = 0_int32

        select case (pattern)
        case (DOSAGE_PATTERN)
            ! only genes with angles below the gene-family median or lower five percentile are marked active
            do i_gene = 1, n_genes
                if (gene_angles(i_gene) <= threshold) then
                    family_idx = gene_to_fam(i_gene)
                    call mask_set_state(masks(:, family_idx), i_gene, .true., ierr)
                    if (is_err(ierr)) return
                end if
            end do
        case (SUBFUNC_PATTERN)
            ! only genes with angles greater than the gene-family median angle are marked active
            do i_gene = 1, n_genes
                if (gene_angles(i_gene) >= threshold) then
                    family_idx = gene_to_fam(i_gene)
                    call mask_set_state(masks(:, family_idx), i_gene, .true., ierr)
                    if (is_err(ierr)) return
                end if
            end do
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine filter_paralogs_by_pattern

    !> The `detect_*` subroutines need a work array for the to be tested subsets.
    !| In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
    !| This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
    !|
    !| This subroutine calculates the needed size for the work array.
    pure subroutine calc_work_arr_paralog_subsets_size(max_subset_size, n_genes, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        integer(int32), intent(inout) :: max_subset_size
            !! maximum size that a subset must not exceed.
            !! @warning
            !! If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.
            !!
            !! Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
            !! @endwarning
        integer(int32), intent(out) :: work_array_size
            !! The calculated needed work array size in absolute worst case scenario. Look into source for details.
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! Output mask with all genes disabled that did not pass the filter
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32), parameter :: max_int32 = huge(0_int32)
        integer(int32) :: i_gene, subset_size, extensions_count, results, previous_results, n_genes_filtered

        call set_ok(ierr)

        if (max_subset_size == 0) then
            work_array_size = 0
            return
        end if

        call validate_dimension_size(n_mask_chunks, ierr)
        call validate_in_range_int(n_genes, ierr, min=1_int32)
        call validate_in_range_int(max_subset_size, ierr, min=1_int32)
        if (is_err(ierr)) return

        n_genes_filtered = 0
        do i_gene = 1, n_genes
            if (mask_check_state(filtered_paralogs_mask, i_gene)) then
                n_genes_filtered = n_genes_filtered + 1
            end if
        end do

        max_subset_size = min(n_genes_filtered, max_subset_size)
        work_array_size = 1
        extensions_count = 1
        previous_results = 0
        results = 0
        do subset_size = 1, max_subset_size
            ! results holds the number of subsets of current subset size that don't have any succeeding genes to be extended with -> result in worst case
            ! previous_results holds the number of results that come from previous subset sizes
            if (previous_results > max_int32 - results) then
                max_subset_size = subset_size - 1
                exit
            end if
            previous_results = previous_results + results
            results = extensions_count - results

            ! calculate the number of extensions of current subsets
            ! overflow check
            if (extensions_count > max_int32 / (n_genes_filtered - subset_size + 1)) then
                max_subset_size = subset_size - 1
                exit
            end if
            extensions_count = extensions_count * (n_genes_filtered - subset_size + 1)
            extensions_count = extensions_count / subset_size

            ! The current subsets will be replaced by their extensions.
            ! In worst case all extended subsets won't be pruned.
            ! Thus, the extensions count will be the work array size, plus the subsets from previous iterations that became results.

            ! if there are less extensions than before, the work array size won't grow anymore
            if (extensions_count > max_int32 - previous_results) then
                max_subset_size = subset_size - 1
                exit
            end if
            if (extensions_count + previous_results < work_array_size) exit
            work_array_size = extensions_count + previous_results
        end do

        ! all subsets with last gene enabled are counted as a result.
        ! as the subset of size 1 with last gene is not a valid subset, remove it (can not be extended, thus also not part of initialization)
        work_array_size = work_array_size - 1
    end subroutine calc_work_arr_paralog_subsets_size

    !> Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.
    pure function mask_get_first_successor_idx(bit_mask) result(idx)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active genes
        integer(int32) :: idx
            !! index of last active gene

        integer(int32) :: i_mask_chunk

        idx = size(bit_mask) * 32
        do i_mask_chunk = size(bit_mask), 1, -1
            idx = idx - leadz(bit_mask(i_mask_chunk))
            if (mod(idx, 32) /= 0) exit
        end do
        idx = idx + 1
    end function mask_get_first_successor_idx

    !> Sets the state of a bit/gene in `bit_mask`
    pure subroutine mask_set_state(bit_mask, i_gene, state, ierr)
        integer(int32), dimension(:), intent(out) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32), intent(in) :: i_gene
            !! index of paralog to be marked active
        logical, intent(in) :: state
            !! state the bit should be set to
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_mask_chunk

        call set_ok(ierr)

        call validate_in_range_int(i_gene, ierr, min=1_int32, max=size(bit_mask, kind=int32)*32_int32)
        if (is_err(ierr)) return

        i_mask_chunk = (i_gene - 1) / 32 + 1

        if (state) then
            bit_mask(i_mask_chunk) = ibset(bit_mask(i_mask_chunk), mod(i_gene - 1, 32))
        else
            bit_mask(i_mask_chunk) = ibclr(bit_mask(i_mask_chunk), mod(i_gene - 1, 32))
        end if
    end subroutine mask_set_state

    !> Checks the state of a bit/paralog in `bit_mask` -> .true. if 1 else .false.
    pure function mask_check_state(bit_mask, i_gene) result(state)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32), intent(in) :: i_gene
            !! index of paralog to be marked active
        logical :: state
            !! check result

        integer(int32) :: i_mask_chunk, ierr

        call set_ok(ierr)
        call validate_in_range_int(i_gene, ierr, min=1_int32, max=size(bit_mask, kind=int32)*32_int32)

        if (is_err(ierr)) then
            state = .false.
        else
            i_mask_chunk = (i_gene - 1) / 32 + 1
            state = btest(bit_mask(i_mask_chunk), mod(i_gene - 1, 32))
        end if

    end function mask_check_state
end module tox_paralog_analysis

!> C-compatible wrapper for detect_neofunctionalization
pure subroutine detect_neofunctionalization_c(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr) bind(C, name="detect_neofunctionalization_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_paralog_analysis, only: detect_neofunctionalization
    use tox_errors, only: set_err, is_err, is_ok, set_ok, ERR_ALLOC_FAIL, validate_dimension_size
    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    M_USE_NULL_VALIDATION

    ! Arguments mapped to C types
    integer(c_int), intent(in), target :: n_axes
    integer(c_int), intent(in), target :: n_genes
    integer(c_int), intent(in), target :: n_families
    real(c_double), dimension(n_axes, n_families), intent(in), target :: ancestors
    real(c_double), dimension(n_axes, n_genes), intent(in), target :: genes
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
    real(c_double), dimension(n_axes), intent(in), target :: thresholds
    integer(c_int), dimension(n_genes, n_axes), intent(out), target :: neofunc
    integer(c_int), intent(out), target :: ierr

    logical, dimension(:, :), allocatable :: neofunc_f

    ! Null checks
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(ancestors)
    M_CHECK_NON_NULL(genes)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(thresholds)
    M_CHECK_NON_NULL(neofunc)

    call set_ok(ierr)
    call validate_dimension_size(n_genes, ierr)
    call validate_dimension_size(n_genes, ierr)
    if (is_err(ierr)) return

    allocate(neofunc_f(n_genes, n_axes), stat=ierr)
    if (is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if

    call c_int_as_logical(neofunc, neofunc_f)

    ! Call the pure Fortran routine
    call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc_f, ierr)

    if (is_ok(ierr)) then
        call logical_as_c_int(neofunc_f, neofunc)
    end if
end subroutine detect_neofunctionalization_c

!> C-compatible wrapper for `detect_dosage_effect`
pure subroutine detect_dosage_effect_c(ancestor, genes, n_genes, n_dims, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, max_angle, gain_gamma, ierr) bind(C, name="detect_dosage_effect_c")
    use tox_paralog_analysis, only: detect_dosage_effect
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_dims
        !! size of `ancestor` vector and vectors in `genes`
    integer(c_int), intent(in), target :: n_genes
        !! number of vectors in `genes`
    integer(c_int), intent(in), target :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
    real(c_double), dimension(n_dims), intent(in), target :: ancestor
        !! expression vector of ancestral ortholog
    real(c_double), dimension(n_dims, n_genes), intent(in), target :: genes
        !! expression vectors of genes
    integer(c_int), intent(out), target :: n_results
        !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
    integer(c_int), intent(in), target :: max_subset_size
        !! maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(c_int), intent(in), target :: n_paralog_subsets
        !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(c_int), dimension(n_mask_chunks, n_paralog_subsets), intent(out), target :: work_arr_paralog_subsets
        !! working array to hold bitmask encoded subsets for detection.
        !! @note
        !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
        !! @endnote
    integer(c_int), dimension(n_mask_chunks), intent(in), target :: filtered_paralogs_mask
        !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    integer(c_int), dimension(n_mask_chunks), intent(out), target :: active_mask
        !! working array to hold the extended subsets
    real(c_double), dimension(n_dims), intent(out), target :: temp_paralog_vector
        !! vector used for pruning subsets
    real(c_double), intent(in), target :: gain_gamma
        !! positive magnitude gain for dosage effect, default 0.1
    real(c_double), intent(in), target :: max_angle
        !! in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
    integer(c_int), intent(out), target :: ierr
        !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_mask_chunks)
    M_CHECK_NON_NULL(ancestor)
    M_CHECK_NON_NULL(genes)
    M_CHECK_NON_NULL(n_results)
    M_CHECK_NON_NULL(max_subset_size)
    M_CHECK_NON_NULL(n_paralog_subsets)
    M_CHECK_NON_NULL(work_arr_paralog_subsets)
    M_CHECK_NON_NULL(filtered_paralogs_mask)
    M_CHECK_NON_NULL(active_mask)
    M_CHECK_NON_NULL(temp_paralog_vector)
    M_CHECK_NON_NULL(gain_gamma)
    M_CHECK_NON_NULL(max_angle)

    call detect_dosage_effect(ancestor, genes, n_genes, n_dims, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, ierr, max_angle, gain_gamma)
end subroutine detect_dosage_effect_c

!> C-compatible wrapper for `detect_subfunctionalization`
pure subroutine detect_subfunctionalization_c(ancestor, genes, n_genes, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr) bind(C, name="detect_subfunctionalization_c")
    use tox_paralog_analysis, only: detect_subfunctionalization
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_dims
        !! size of `ancestor` vector and vectors in `genes`
    integer(c_int), intent(in), target :: n_genes
        !! number of vectors in `genes`
    integer(c_int), intent(in), target :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
    real(c_double), dimension(n_dims), intent(in), target :: ancestor
        !! expression vector of ancestral ortholog
    real(c_double), dimension(n_dims, n_genes), intent(in), target :: genes
        !! expression vectors of genes
    real(c_double), intent(in), target :: rdi_threshold
        !! max allowed residual distance from `ancestor`
    integer(c_int), intent(out), target :: n_results
        !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
    integer(c_int), intent(in), target :: max_subset_size
        !! maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(c_int), intent(in), target :: n_paralog_subsets
        !! number of gene subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(c_int), dimension(n_mask_chunks, n_paralog_subsets), intent(out), target :: work_arr_paralog_subsets
        !! working array to hold bitmask encoded subsets for detection.
        !! @note
        !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
        !! @endnote
    integer(c_int), dimension(n_mask_chunks), intent(in), target :: filtered_paralogs_mask
        !! bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    integer(c_int), dimension(n_mask_chunks), intent(out), target :: active_mask
        !! working array to hold the extended subsets
    real(c_double), dimension(n_dims), intent(out), target :: temp_paralog_vector
        !! vector used for pruning subsets
    real(c_double), dimension(n_genes), intent(in), target :: paralog_norms
        !! needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` function from `f42_utils` function for this)
    integer(c_int), dimension(n_genes), intent(in), target :: sorted_paralog_norms_perm
        !! needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
    real(c_double), dimension(n_genes), intent(out), target :: temp_work_array
        !! needed for efficient check of minimum value after a certain index
    integer(c_int), intent(out), target :: ierr
        !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_mask_chunks)
    M_CHECK_NON_NULL(ancestor)
    M_CHECK_NON_NULL(genes)
    M_CHECK_NON_NULL(rdi_threshold)
    M_CHECK_NON_NULL(n_results)
    M_CHECK_NON_NULL(max_subset_size)
    M_CHECK_NON_NULL(n_paralog_subsets)
    M_CHECK_NON_NULL(work_arr_paralog_subsets)
    M_CHECK_NON_NULL(filtered_paralogs_mask)
    M_CHECK_NON_NULL(active_mask)
    M_CHECK_NON_NULL(temp_paralog_vector)
    M_CHECK_NON_NULL(paralog_norms)
    M_CHECK_NON_NULL(sorted_paralog_norms_perm)
    M_CHECK_NON_NULL(temp_work_array)

    call detect_subfunctionalization(ancestor, genes, n_genes, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)
end subroutine detect_subfunctionalization_c

!> C-compatible wrapper for `filter_paralogs_by_pattern_subfunctionalization`
pure subroutine filter_paralogs_by_pattern_subfunctionalization_c(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr) bind(C, name="filter_paralogs_by_pattern_subfunctionalization_c")
    use tox_paralog_analysis, only: filter_paralogs_by_pattern_subfunctionalization
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! number of genes
    integer(c_int), intent(in), target :: n_families
        !! number of families
    integer(c_int), intent(in), target :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_genes` genes
    real(c_double), dimension(n_genes), intent(in), target :: gene_angles
        !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
    real(c_double), intent(in), target :: threshold
        !! filter threshold
    integer(c_int), dimension(n_mask_chunks, n_families), intent(out), target :: masks
        !! bit mask that will have indices of genes kept by pattern set to 1, else 0
    integer(c_int), intent(out), target :: ierr
        !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_mask_chunks)
    M_CHECK_NON_NULL(gene_angles)
    M_CHECK_NON_NULL(threshold)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(masks)

    call filter_paralogs_by_pattern_subfunctionalization(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
end subroutine filter_paralogs_by_pattern_subfunctionalization_c

!> C-compatible wrapper for `filter_paralogs_by_pattern_dosage_effect`
pure subroutine filter_paralogs_by_pattern_dosage_effect_c(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr) bind(C, name="filter_paralogs_by_pattern_dosage_effect")
    use tox_paralog_analysis, only: filter_paralogs_by_pattern_dosage_effect
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! number of genes
    integer(c_int), intent(in), target :: n_families
        !! number of families
    integer(c_int), intent(in), target :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_genes` genes
    real(c_double), dimension(n_genes), intent(in), target :: gene_angles
        !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    integer(c_int), dimension(n_genes), intent(in), target :: gene_to_fam
        !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
    real(c_double), intent(in), target :: threshold
        !! filter threshold
    integer(c_int), dimension(n_mask_chunks, n_families), intent(out), target :: masks
        !! bit mask that will have indices of genes kept by pattern set to 1, else 0
    integer(c_int), intent(out), target :: ierr
        !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_mask_chunks)
    M_CHECK_NON_NULL(gene_angles)
    M_CHECK_NON_NULL(threshold)
    M_CHECK_NON_NULL(gene_to_fam)
    M_CHECK_NON_NULL(n_families)
    M_CHECK_NON_NULL(masks)

    call filter_paralogs_by_pattern_dosage_effect(gene_angles, threshold, n_genes, n_families, gene_to_fam, masks, n_mask_chunks, ierr)
end subroutine filter_paralogs_by_pattern_dosage_effect_c

!> C-compatible wrapper for `calc_work_arr_paralog_subsets_size`
pure subroutine calc_work_arr_paralog_subsets_size_c(max_subset_size, n_genes, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr) bind(C, name="calc_work_arr_paralog_subsets_size")
    use tox_paralog_analysis, only: calc_work_arr_paralog_subsets_size
    use, intrinsic :: iso_c_binding, only: c_int
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! number of paralogs
    integer(c_int), intent(in), target :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_genes` genes
    integer(c_int), intent(inout), target :: max_subset_size
        !! maximum size that a subset must not exceed.
        !! @warning
        !! If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.
        !!
        !! Also, size will be set to number of paralogs in `filtered_paralogs_mask` if larger.
        !! @endwarning
    integer(c_int), intent(out), target :: work_array_size
        !! The calculated needed work array size in absolute worst case scenario. Look into source for details.
    integer(c_int), dimension(n_mask_chunks), intent(in), target :: filtered_paralogs_mask
        !! Output mask with all paralogs disabled that did not pass the filter
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_mask_chunks)
    M_CHECK_NON_NULL(max_subset_size)
    M_CHECK_NON_NULL(work_array_size)
    M_CHECK_NON_NULL(filtered_paralogs_mask)

    call calc_work_arr_paralog_subsets_size(max_subset_size, n_genes, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
end subroutine calc_work_arr_paralog_subsets_size_c

!> C-compatible wrapper for `mask_chunk_count`
pure subroutine mask_chunk_count_c(n_genes, count, ierr) bind(C, name="mask_chunk_count_c")
    use tox_paralog_analysis, only: mask_chunk_count
    use, intrinsic :: iso_c_binding, only: c_int
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! number of paralogs
    integer(c_int), intent(out), target :: count
        !! number of 32 bit chunks a mask needs to encode `n_genes` genes
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(count)

    call mask_chunk_count(n_genes, count)
end subroutine mask_chunk_count_c

!> C-compatible wrapper for `mask_check_state`
pure subroutine mask_check_state_c(bit_mask, n_mask_chunks, i_gene, state, ierr) bind(C, name="mask_check_state_c")
    use tox_paralog_analysis, only: mask_check_state
    use tox_conversions, only: logical_as_c_int
    use, intrinsic :: iso_c_binding, only: c_int
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_genes` genes
    integer(c_int), dimension(n_mask_chunks), intent(in), target :: bit_mask
        !! chunked mask to mark active paralogs
    integer(c_int), intent(in), target :: i_gene
        !! index of paralog to be checked, starting with 0
    integer(c_int), intent(out), target :: state
        !! check result
    integer(c_int), intent(out), target :: ierr

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_mask_chunks)
    M_CHECK_NON_NULL(bit_mask)
    M_CHECK_NON_NULL(i_gene)
    M_CHECK_NON_NULL(state)

    call logical_as_c_int(mask_check_state(bit_mask, i_gene), state)
end subroutine mask_check_state_c
