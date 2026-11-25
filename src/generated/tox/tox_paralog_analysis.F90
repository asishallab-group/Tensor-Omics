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
!|
!| Generated from [[tox_paralog_analysis_impl(module)]]; do not edit -- regenerate instead.
module tox_paralog_analysis
    use f42_safeguard
    use tox_paralog_analysis_impl, only: MODE_DOSAGE_PATTERN, MODE_SUBFUNC_PATTERN, calc_paralog_norms_impl, calc_paralog_subsets_bitmasks_size_impl
    use tox_paralog_analysis_impl, only: detect_neofunctionalization_impl, detect_patterns_impl, filter_paralogs_by_pattern_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: PI, above
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: get_err_arg_pos, set_err, set_err_arg_pos, validate_all_in_range_int
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: detect_neofunctionalization
    public :: detect_dosage_effect
    public :: detect_dosage_effect_expert
    public :: detect_subfunctionalization
    public :: detect_subfunctionalization_expert
    public :: filter_paralogs_by_pattern_dosage_effect
    public :: filter_paralogs_by_pattern_subfunctionalization
    public :: calc_paralog_norms
    public :: calc_paralog_subsets_bitmasks_size

contains

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_impl(module):detect_neofunctionalization_impl]].
    pure subroutine detect_neofunctionalization(&
            ancestors,&
            n_families,&
            genes,&
            n_axes,&
            gene_to_fam,&
            n_genes,&
            thresholds,&
            neofunc,&
            ierr&
        )
        integer(int32), intent(in) :: n_families
            !! number of vectors in `ancestors`
        integer(int32), intent(in) :: n_axes
            !! size of `ancestors` vector and vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        real(real64), dimension(n_axes, n_families), intent(in) :: ancestors
            !! RAP projected unit length expression vector of ancestral ortholog
            !! The minimum valid value is `-1.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_axes, n_genes), intent(in) :: genes
            !! RAP projected unit length expression vectors of genes
            !! The minimum valid value is `-1.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0_int32` for unassigned genes
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
            !! The value `0_int32` is additionally accepted.
        real(real64), dimension(n_axes), intent(in) :: thresholds
            !! threshold per axis that defines significant change in expression, may be a percentile of all genes' changes per axis
            !! The minimum valid value is `-1.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), dimension(n_genes, n_axes), intent(out) :: neofunc
            !! `.true.` if neofunctionalization has been detected for the respective axes, always `.false.` for unassigned genes
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_axes, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=6_int32)
        call validate_all_in_range_real(ancestors, n_axes * n_families, ierr, arg_pos=1_int32, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_real(genes, n_axes * n_genes, ierr, arg_pos=3_int32, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=5_int32, min=1_int32, max=n_families, sentinel=0_int32)
        call validate_all_in_range_real(thresholds, n_axes, ierr, arg_pos=7_int32, min=-1.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        call detect_neofunctionalization_impl(&
            ancestors = ancestors,&
            n_families = n_families,&
            genes = genes,&
            n_axes = n_axes,&
            gene_to_fam = gene_to_fam,&
            n_genes = n_genes,&
            thresholds = thresholds,&
            neofunc = neofunc&
        )
    end subroutine detect_neofunctionalization

    !> summary: Validates its inputs, prepares what [[tox_paralog_analysis_impl(module):detect_patterns_impl]] needs, then calls it. The entry point to reach for first; see [[tox_paralog_analysis(module):detect_dosage_effect_expert]] to prepare it yourself.
    !| A subset is grown from its genes in index order: it starts from one kept gene and is extended
    !| only by kept genes after its last one. A single gene is never a result; a result has two genes
    !| at least, and is not extended further. So which subsets are found depends on the genes' order,
    !| and a result's superset is found too when the superset's own prefix did not meet the criterion
    !| (issue #244).
    pure subroutine detect_dosage_effect(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_bitmask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            paralog_subsets_bitmasks,&
            n_paralog_subsets,&
            max_angle,&
            gain_gamma,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes; `bitmask_n_words` computes it
            !! The minimum valid value is `((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32)))`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `paralog_subsets_bitmasks`. Every kept
            !! gene but the last starts a subset, so fewer slots than that are refused with
            !! `ERR_SIZE_MISMATCH` before any subset is stored, and so is a search that runs out of
            !! slots later; zero slots are enough when fewer than two genes are kept.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis(module):calc_paralog_subsets_bitmasks_size]].
            !! The minimum valid value is `0_int32`.
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
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine.
            !! A bit set past `n_genes` is refused with `ERR_INVALID_INPUT`; the sizing routine
            !! ignores such bits.
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
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: paralog_subsets_bitmasks
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each column is a bit mask of the `n_genes` genes, laid out as [[f42_bit_masks_impl(module)]] describes.
            !! @endnote
        real(real64), intent(in), optional :: max_angle
            !! maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
            !! The default value is `4.0_real64*atan(1.0_real64)`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `PI`.
        real(real64), intent(in), optional :: gain_gamma
            !! positive magnitude gain for dosage effect
            !! The default value is `0.1_real64`.
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), dimension(:), allocatable :: tmp_active_bitmask
        real(real64), dimension(:), allocatable :: tmp_paralog_vector

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32))))
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_in_range_int(n_paralog_subsets, ierr, arg_pos=10_int32, min=0_int32)
        call validate_in_range_real(max_angle, ierr, arg_pos=11_int32, min=0.0_real64, max=PI)
        call validate_in_range_real(gain_gamma, ierr, arg_pos=12_int32, min=above(0.0_real64))
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_active_bitmask(n_mask_chunks))
        M_ALLOCATE(tmp_paralog_vector(n_dims))

        call detect_patterns_impl(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_DOSAGE_PATTERN,&
            filtered_paralogs_bitmask = filtered_paralogs_bitmask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            paralog_subsets_bitmasks = paralog_subsets_bitmasks,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_bitmask = tmp_active_bitmask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            max_angle = max_angle,&
            gain_gamma = gain_gamma,&
            ierr = ierr&
        )
        if (is_err(ierr)) then
            select case (get_err_arg_pos(ierr))
            case (1_int32:4_int32)
                ! the same position here
            case (6_int32)
                call set_err_arg_pos(ierr, 5_int32)
            case (7_int32)
                call set_err_arg_pos(ierr, 6_int32)
            case (8_int32)
                call set_err_arg_pos(ierr, 7_int32)
            case (9_int32)
                call set_err_arg_pos(ierr, 8_int32)
            case (10_int32)
                call set_err_arg_pos(ierr, 9_int32)
            case (11_int32)
                call set_err_arg_pos(ierr, 10_int32)
            case (14_int32)
                call set_err_arg_pos(ierr, 11_int32)
            case (15_int32)
                call set_err_arg_pos(ierr, 12_int32)
            case default
                call clear_err_arg_pos(ierr)
            end select
            return
        end if
    end subroutine detect_dosage_effect

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_impl(module):detect_patterns_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_paralog_analysis(module):detect_dosage_effect]] does both.
    !| A subset is grown from its genes in index order: it starts from one kept gene and is extended
    !| only by kept genes after its last one. A single gene is never a result; a result has two genes
    !| at least, and is not extended further. So which subsets are found depends on the genes' order,
    !| and a result's superset is found too when the superset's own prefix did not meet the criterion
    !| (issue #244).
    pure subroutine detect_dosage_effect_expert(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_bitmask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            paralog_subsets_bitmasks,&
            n_paralog_subsets,&
            tmp_active_bitmask,&
            tmp_paralog_vector,&
            max_angle,&
            gain_gamma,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes; `bitmask_n_words` computes it
            !! The minimum valid value is `((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32)))`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `paralog_subsets_bitmasks`. Every kept
            !! gene but the last starts a subset, so fewer slots than that are refused with
            !! `ERR_SIZE_MISMATCH` before any subset is stored, and so is a search that runs out of
            !! slots later; zero slots are enough when fewer than two genes are kept.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis(module):calc_paralog_subsets_bitmasks_size]].
            !! The minimum valid value is `0_int32`.
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
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine.
            !! A bit set past `n_genes` is refused with `ERR_INVALID_INPUT`; the sizing routine
            !! ignores such bits.
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
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: paralog_subsets_bitmasks
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each column is a bit mask of the `n_genes` genes, laid out as [[f42_bit_masks_impl(module)]] describes.
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(out) :: tmp_active_bitmask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: tmp_paralog_vector
            !! vector used for pruning subsets
        real(real64), intent(in), optional :: max_angle
            !! maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned
            !! The default value is `4.0_real64*atan(1.0_real64)`.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `PI`.
        real(real64), intent(in), optional :: gain_gamma
            !! positive magnitude gain for dosage effect
            !! The default value is `0.1_real64`.
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32))))
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_in_range_int(n_paralog_subsets, ierr, arg_pos=10_int32, min=0_int32)
        call validate_in_range_real(max_angle, ierr, arg_pos=13_int32, min=0.0_real64, max=PI)
        call validate_in_range_real(gain_gamma, ierr, arg_pos=14_int32, min=above(0.0_real64))
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        call detect_patterns_impl(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_DOSAGE_PATTERN,&
            filtered_paralogs_bitmask = filtered_paralogs_bitmask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            paralog_subsets_bitmasks = paralog_subsets_bitmasks,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_bitmask = tmp_active_bitmask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            max_angle = max_angle,&
            gain_gamma = gain_gamma,&
            ierr = ierr&
        )
        if (is_err(ierr)) then
            select case (get_err_arg_pos(ierr))
            case (1_int32:4_int32)
                ! the same position here
            case (6_int32)
                call set_err_arg_pos(ierr, 5_int32)
            case (7_int32)
                call set_err_arg_pos(ierr, 6_int32)
            case (8_int32)
                call set_err_arg_pos(ierr, 7_int32)
            case (9_int32)
                call set_err_arg_pos(ierr, 8_int32)
            case (10_int32)
                call set_err_arg_pos(ierr, 9_int32)
            case (11_int32)
                call set_err_arg_pos(ierr, 10_int32)
            case (12_int32)
                call set_err_arg_pos(ierr, 11_int32)
            case (13_int32)
                call set_err_arg_pos(ierr, 12_int32)
            case (14_int32)
                call set_err_arg_pos(ierr, 13_int32)
            case (15_int32)
                call set_err_arg_pos(ierr, 14_int32)
            case default
                call clear_err_arg_pos(ierr)
            end select
            return
        end if
    end subroutine detect_dosage_effect_expert

    !> summary: Validates its inputs, prepares what [[tox_paralog_analysis_impl(module):detect_patterns_impl]] needs, then calls it. The entry point to reach for first; see [[tox_paralog_analysis(module):detect_subfunctionalization_expert]] to prepare it yourself.
    !| A subset is grown from its genes in index order: it starts from one kept gene and is extended
    !| only by kept genes after its last one. A single gene is never a result; a result has two genes
    !| at least, and is not extended further. So which subsets are found depends on the genes' order,
    !| and a result's superset is found too when the superset's own prefix did not meet the criterion
    !| (issue #244).
    pure subroutine detect_subfunctionalization(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_bitmask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            paralog_subsets_bitmasks,&
            n_paralog_subsets,&
            rdi_threshold,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes; `bitmask_n_words` computes it
            !! The minimum valid value is `((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32)))`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `paralog_subsets_bitmasks`. Every kept
            !! gene but the last starts a subset, so fewer slots than that are refused with
            !! `ERR_SIZE_MISMATCH` before any subset is stored, and so is a search that runs out of
            !! slots later; zero slots are enough when fewer than two genes are kept.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis(module):calc_paralog_subsets_bitmasks_size]].
            !! The minimum valid value is `0_int32`.
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
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine.
            !! A bit set past `n_genes` is refused with `ERR_INVALID_INPUT`; the sizing routine
            !! ignores such bits.
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
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: paralog_subsets_bitmasks
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each column is a bit mask of the `n_genes` genes, laid out as [[f42_bit_masks_impl(module)]] describes.
            !! @endnote
        real(real64), intent(in) :: rdi_threshold
            !! max allowed residual distance from `ancestor`
            !! The minimum valid value is `0.0_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), dimension(:), allocatable :: tmp_active_bitmask
        real(real64), dimension(:), allocatable :: tmp_paralog_vector
        real(real64), dimension(:), allocatable :: paralog_norms
        integer(int32), dimension(:), allocatable :: paralog_norms_perm
        real(real64), dimension(:), allocatable :: tmp_work_array

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32))))
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_in_range_int(n_paralog_subsets, ierr, arg_pos=10_int32, min=0_int32)
        call validate_in_range_real(rdi_threshold, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(paralog_norms(n_genes))
        call calc_paralog_norms(&
            genes = genes,&
            n_dims = n_dims,&
            n_genes = n_genes,&
            paralog_norms = paralog_norms,&
            ierr = ierr&
        )
        if (is_err(ierr)) then
            select case (get_err_arg_pos(ierr))
            case (3_int32)
                ! the same position here
            case (1_int32)
                call set_err_arg_pos(ierr, 2_int32)
            case (2_int32)
                call set_err_arg_pos(ierr, 4_int32)
            case default
                call clear_err_arg_pos(ierr)
            end select
            return
        end if
        M_ALLOCATE(tmp_active_bitmask(n_mask_chunks))
        M_ALLOCATE(tmp_paralog_vector(n_dims))
        M_ALLOCATE(paralog_norms_perm(n_genes))
        M_ALLOCATE(tmp_work_array(n_genes))
        call init_perm(paralog_norms_perm)
        call sort_array_heapsort(paralog_norms, paralog_norms_perm)

        call detect_patterns_impl(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_SUBFUNC_PATTERN,&
            filtered_paralogs_bitmask = filtered_paralogs_bitmask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            paralog_subsets_bitmasks = paralog_subsets_bitmasks,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_bitmask = tmp_active_bitmask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            rdi_threshold = rdi_threshold,&
            paralog_norms = paralog_norms,&
            paralog_norms_perm = paralog_norms_perm,&
            tmp_work_array = tmp_work_array,&
            ierr = ierr&
        )
        if (is_err(ierr)) then
            select case (get_err_arg_pos(ierr))
            case (1_int32:4_int32)
                ! the same position here
            case (6_int32)
                call set_err_arg_pos(ierr, 5_int32)
            case (7_int32)
                call set_err_arg_pos(ierr, 6_int32)
            case (8_int32)
                call set_err_arg_pos(ierr, 7_int32)
            case (9_int32)
                call set_err_arg_pos(ierr, 8_int32)
            case (10_int32)
                call set_err_arg_pos(ierr, 9_int32)
            case (11_int32)
                call set_err_arg_pos(ierr, 10_int32)
            case (16_int32)
                call set_err_arg_pos(ierr, 11_int32)
            case default
                call clear_err_arg_pos(ierr)
            end select
            return
        end if
    end subroutine detect_subfunctionalization

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_impl(module):detect_patterns_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_paralog_analysis(module):detect_subfunctionalization]] does both.
    !| A subset is grown from its genes in index order: it starts from one kept gene and is extended
    !| only by kept genes after its last one. A single gene is never a result; a result has two genes
    !| at least, and is not extended further. So which subsets are found depends on the genes' order,
    !| and a result's superset is found too when the superset's own prefix did not meet the criterion
    !| (issue #244).
    pure subroutine detect_subfunctionalization_expert(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_bitmask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            paralog_subsets_bitmasks,&
            n_paralog_subsets,&
            tmp_active_bitmask,&
            tmp_paralog_vector,&
            rdi_threshold,&
            paralog_norms,&
            paralog_norms_perm,&
            tmp_work_array,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes; `bitmask_n_words` computes it
            !! The minimum valid value is `((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32)))`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `paralog_subsets_bitmasks`. Every kept
            !! gene but the last starts a subset, so fewer slots than that are refused with
            !! `ERR_SIZE_MISMATCH` before any subset is stored, and so is a search that runs out of
            !! slots later; zero slots are enough when fewer than two genes are kept.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis(module):calc_paralog_subsets_bitmasks_size]].
            !! The minimum valid value is `0_int32`.
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
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine.
            !! A bit set past `n_genes` is refused with `ERR_INVALID_INPUT`; the sizing routine
            !! ignores such bits.
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
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: paralog_subsets_bitmasks
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each column is a bit mask of the `n_genes` genes, laid out as [[f42_bit_masks_impl(module)]] describes.
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(out) :: tmp_active_bitmask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: tmp_paralog_vector
            !! vector used for pruning subsets
        real(real64), intent(in) :: rdi_threshold
            !! max allowed residual distance from `ancestor`
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_genes), intent(in) :: paralog_norms
            !! euclidean norm of each gene's vector, used for subset pruning. Where the caller
            !! passes them, they are trusted to be the genes' own.
            !! The minimum valid value is `0.0_real64`.
            !! It is *VERY IMPORTANT* to compute this argument from the `paralog_norms` output produced by [[tox_paralog_analysis(module):calc_paralog_norms]].
        integer(int32), dimension(n_genes), intent(in) :: paralog_norms_perm
            !! ascending permutation of `paralog_norms`, for subset pruning: a subset is extended only
            !! while the smallest norm among all genes after its last one -- kept by the filter or
            !! not -- is at most the length of its residual. Where the caller passes it, the order
            !! is the caller's job: only the entries' range is checked, so no index is out of
            !! bounds, and a wrong order gives wrong pruning -- wrong minima for a permutation that
            !! is not ascending, and for an order with a repeated entry (such as [1, 1, 1]) minima
            !! never written, read uninitialised.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        real(real64), dimension(n_genes), intent(out) :: tmp_work_array
            !! work array for checking the minimum value after a given index efficiently
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32))))
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_in_range_int(n_paralog_subsets, ierr, arg_pos=10_int32, min=0_int32)
        call validate_in_range_real(rdi_threshold, ierr, arg_pos=13_int32, min=0.0_real64)
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(paralog_norms, n_genes, ierr, arg_pos=14_int32, min=0.0_real64)
        call validate_all_in_range_int(paralog_norms_perm, n_genes, ierr, arg_pos=15_int32, min=1_int32, max=n_genes)
        if (is_err(ierr)) return
#endif

        call detect_patterns_impl(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_SUBFUNC_PATTERN,&
            filtered_paralogs_bitmask = filtered_paralogs_bitmask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            paralog_subsets_bitmasks = paralog_subsets_bitmasks,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_bitmask = tmp_active_bitmask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            rdi_threshold = rdi_threshold,&
            paralog_norms = paralog_norms,&
            paralog_norms_perm = paralog_norms_perm,&
            tmp_work_array = tmp_work_array,&
            ierr = ierr&
        )
        if (is_err(ierr)) then
            select case (get_err_arg_pos(ierr))
            case (1_int32:4_int32)
                ! the same position here
            case (6_int32)
                call set_err_arg_pos(ierr, 5_int32)
            case (7_int32)
                call set_err_arg_pos(ierr, 6_int32)
            case (8_int32)
                call set_err_arg_pos(ierr, 7_int32)
            case (9_int32)
                call set_err_arg_pos(ierr, 8_int32)
            case (10_int32)
                call set_err_arg_pos(ierr, 9_int32)
            case (11_int32)
                call set_err_arg_pos(ierr, 10_int32)
            case (12_int32)
                call set_err_arg_pos(ierr, 11_int32)
            case (13_int32)
                call set_err_arg_pos(ierr, 12_int32)
            case (16_int32)
                call set_err_arg_pos(ierr, 13_int32)
            case (17_int32)
                call set_err_arg_pos(ierr, 14_int32)
            case (18_int32)
                call set_err_arg_pos(ierr, 15_int32)
            case (19_int32)
                call set_err_arg_pos(ierr, 16_int32)
            case default
                call clear_err_arg_pos(ierr)
            end select
            return
        end if
    end subroutine detect_subfunctionalization_expert

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_impl(module):filter_paralogs_by_pattern_impl]].
    !| This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    pure subroutine filter_paralogs_by_pattern_dosage_effect(&
            gene_angles,&
            threshold,&
            n_genes,&
            n_families,&
            gene_to_fam,&
            family_bitmasks,&
            n_mask_chunks,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(in) :: n_families
            !! number of families
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! The minimum valid value is `((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32)))`.
        real(real64), dimension(n_genes), intent(in) :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `PI`.
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: family_bitmasks
            !! bit mask that will have the indices of genes kept by this pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_real(threshold, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_families, ierr, arg_pos=4_int32, min=1_int32, max=n_genes)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=7_int32, min=((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32))))
        call validate_all_in_range_real(gene_angles, n_genes, ierr, arg_pos=1_int32, min=0.0_real64, max=PI)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=5_int32, min=1_int32, max=n_families)
        if (is_err(ierr)) return
#endif

        call filter_paralogs_by_pattern_impl(&
            pattern_mode = MODE_DOSAGE_PATTERN,&
            gene_angles = gene_angles,&
            threshold = threshold,&
            n_genes = n_genes,&
            n_families = n_families,&
            gene_to_fam = gene_to_fam,&
            family_bitmasks = family_bitmasks,&
            n_mask_chunks = n_mask_chunks&
        )
    end subroutine filter_paralogs_by_pattern_dosage_effect

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_impl(module):filter_paralogs_by_pattern_impl]].
    !| This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    pure subroutine filter_paralogs_by_pattern_subfunctionalization(&
            gene_angles,&
            threshold,&
            n_genes,&
            n_families,&
            gene_to_fam,&
            family_bitmasks,&
            n_mask_chunks,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of genes
        integer(int32), intent(in) :: n_families
            !! number of families
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
            !! The minimum valid value is `((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32)))`.
        real(real64), dimension(n_genes), intent(in) :: gene_angles
            !! vector, holding the angles between ancestor and genes (0<=angle<=Pi)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `PI`.
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_genes), intent(in) :: gene_to_fam
            !! a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_families`.
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: family_bitmasks
            !! bit mask that will have the indices of genes kept by this pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_real(threshold, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_families, ierr, arg_pos=4_int32, min=1_int32, max=n_genes)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=7_int32, min=((n_genes)/32_int32 + min(1_int32, mod((n_genes), 32_int32))))
        call validate_all_in_range_real(gene_angles, n_genes, ierr, arg_pos=1_int32, min=0.0_real64, max=PI)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=5_int32, min=1_int32, max=n_families)
        if (is_err(ierr)) return
#endif

        call filter_paralogs_by_pattern_impl(&
            pattern_mode = MODE_SUBFUNC_PATTERN,&
            gene_angles = gene_angles,&
            threshold = threshold,&
            n_genes = n_genes,&
            n_families = n_families,&
            gene_to_fam = gene_to_fam,&
            family_bitmasks = family_bitmasks,&
            n_mask_chunks = n_mask_chunks&
        )
    end subroutine filter_paralogs_by_pattern_subfunctionalization

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_impl(module):calc_paralog_norms_impl]].
    !| Each norm is f42's scaled `norm`, so no vector of finite components overflows or underflows
    !| on the way: the result is +Inf only where the true norm exceeds `huge`.
    pure subroutine calc_paralog_norms(&
            genes,&
            n_dims,&
            n_genes,&
            paralog_norms,&
            ierr&
        )
        integer(int32), intent(in) :: n_dims
            !! size of the vectors in `genes`
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        real(real64), dimension(n_genes), intent(out) :: paralog_norms
            !! euclidean norm of each gene's vector
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dims, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call calc_paralog_norms_impl(&
            genes = genes,&
            n_dims = n_dims,&
            n_genes = n_genes,&
            paralog_norms = paralog_norms&
        )
    end subroutine calc_paralog_norms

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_impl(module):calc_paralog_subsets_bitmasks_size_impl]].
    !| The `detect_*` subroutines need a work array for the subsets they try. It holds, side by side,
    !| the results found so far, the subsets still to be grown and the ones grown from them, so it is
    !| sized for the worst case: every subset of the first F - 1 kept genes stays a candidate, and every
    !| subset that adds the last kept gene becomes a result. With F genes kept and subsets of at most m
    !| genes, that is one subset of the first F - 1 genes for each size 1 to m: the sum over
    !| k = 1 .. min(m, F - 1) of C(F - 1, k), which is 2**(F - 1) - 1 for m = F, and 0 for fewer
    !| than two kept genes.
    pure subroutine calc_paralog_subsets_bitmasks_size(&
            max_subset_size,&
            n_genes,&
            work_array_size,&
            filtered_paralogs_bitmask,&
            n_mask_chunks,&
            ierr&
        )
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes
        integer(int32), intent(inout) :: max_subset_size
            !! maximum size that a subset must not exceed. Zero is in range and means there is
            !! nothing to size a work array for, which is reported back as a size of zero.
            !! The minimum valid value is `0_int32`.
            !! @warning
            !! It comes back lowered to the number of kept genes if it is larger, and lowered further
            !! to the largest size whose work array holds at most 2**31 - 1 (2,147,483,647) subsets, the
            !! maximum valid size, where the size asked for would need more.
            !! @endwarning
        integer(int32), intent(in) :: n_genes
            !! number of genes. It is checked before `n_mask_chunks`, so a mask of too few words, none
            !! included, is reported as `n_genes` out of range.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_mask_chunks * 32`.
        integer(int32), intent(out) :: work_array_size
            !! The number of subsets the work array must hold in the worst case, for the
            !! `max_subset_size` handed back
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_bitmask
            !! bit mask of the genes the filter kept. Only the first `n_genes` bits are counted;
            !! a bit past them is ignored here, while the detections refuse it.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(max_subset_size, ierr, arg_pos=1_int32, min=0_int32)
        call validate_in_range_int(n_genes, ierr, arg_pos=2_int32, min=1_int32, max=n_mask_chunks * 32)
        call validate_dimension_size(n_mask_chunks, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call calc_paralog_subsets_bitmasks_size_impl(&
            max_subset_size = max_subset_size,&
            n_genes = n_genes,&
            work_array_size = work_array_size,&
            filtered_paralogs_bitmask = filtered_paralogs_bitmask,&
            n_mask_chunks = n_mask_chunks&
        )
    end subroutine calc_paralog_subsets_bitmasks_size

end module tox_paralog_analysis
