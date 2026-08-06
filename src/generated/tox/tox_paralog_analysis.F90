#include <src/macros.h>

!> summary: Wrappers for [[tox_paralog_analysis_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_paralog_analysis
    use tox_paralog_analysis_kernel, only: MODE_DOSAGE_PATTERN, MODE_SUBFUNC_PATTERN, detect_neofunctionalization_kernel, detect_patterns_kernel
    use tox_paralog_analysis_kernel, only: filter_paralogs_by_pattern_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math, only: PI, above
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size
    use tox_errors, only: validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: detect_neofunctionalization
    public :: detect_dosage_effect
    public :: detect_dosage_effect_alloc
    public :: detect_subfunctionalization
    public :: detect_subfunctionalization_alloc
    public :: filter_paralogs_by_pattern_dosage_effect
    public :: filter_paralogs_by_pattern_subfunctionalization

contains

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_kernel(module):detect_neofunctionalization_kernel]].
    subroutine detect_neofunctionalization(&
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
        logical, dimension(n_genes, n_axes), intent(out) :: neofunc
            !! `.true.` if neofunctionalization has been detected for the respective axes, always `.false.` for unassigned genes
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
        call validate_dimension_size(n_families, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_axes, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=6_int32)
        call validate_all_in_range_real(ancestors, n_axes * n_families, ierr, arg_pos=1_int32, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_real(genes, n_axes * n_genes, ierr, arg_pos=3_int32, min=-1.0_real64, max=1.0_real64)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=5_int32, min=1_int32, max=n_families, sentinel=0_int32)
        call validate_all_in_range_real(thresholds, n_axes, ierr, arg_pos=7_int32, min=-1.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        call detect_neofunctionalization_kernel(&
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

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_kernel(module):detect_patterns_kernel]].
    subroutine detect_dosage_effect(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            work_arr_paralog_subsets,&
            n_paralog_subsets,&
            tmp_active_mask,&
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
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
            !! The minimum valid value is `(n_genes + 31) / 32`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size]].
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. Too large a value is capped to the
            !! maximum valid size. The bindings cap it automatically while sizing the work
            !! array; a Fortran caller caps it by calling
            !! [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
            !! Zero is in range and means there is no subset to check -- the sizing routine reports
            !! it whenever the filtered families hold a single gene each. It reports a work array
            !! of zero slots along with it, which this routine does not accept, so a caller that
            !! gets zero back has nothing to detect and should not call here at all.
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(out) :: tmp_active_mask
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
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=(n_genes + 31) / 32)
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_dimension_size(n_paralog_subsets, ierr, arg_pos=10_int32)
        call validate_in_range_real(max_angle, ierr, arg_pos=13_int32, min=0.0_real64, max=PI)
        call validate_in_range_real(gain_gamma, ierr, arg_pos=14_int32, min=above(0.0_real64))
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        call detect_patterns_kernel(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_DOSAGE_PATTERN,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            work_arr_paralog_subsets = work_arr_paralog_subsets,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_mask = tmp_active_mask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            max_angle = max_angle,&
            gain_gamma = gain_gamma,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine detect_dosage_effect

    !> summary: Allocates its work arrays, then calls [[tox_paralog_analysis_kernel(module):detect_patterns_kernel]].
    subroutine detect_dosage_effect_alloc(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            work_arr_paralog_subsets,&
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
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
            !! The minimum valid value is `(n_genes + 31) / 32`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size]].
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. Too large a value is capped to the
            !! maximum valid size. The bindings cap it automatically while sizing the work
            !! array; a Fortran caller caps it by calling
            !! [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
            !! Zero is in range and means there is no subset to check -- the sizing routine reports
            !! it whenever the filtered families hold a single gene each. It reports a work array
            !! of zero slots along with it, which this routine does not accept, so a caller that
            !! gets zero back has nothing to detect and should not call here at all.
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
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
        integer(int32), dimension(:), allocatable :: tmp_active_mask
        real(real64), dimension(:), allocatable :: tmp_paralog_vector

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=(n_genes + 31) / 32)
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_dimension_size(n_paralog_subsets, ierr, arg_pos=10_int32)
        call validate_in_range_real(max_angle, ierr, arg_pos=11_int32, min=0.0_real64, max=PI)
        call validate_in_range_real(gain_gamma, ierr, arg_pos=12_int32, min=above(0.0_real64))
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_active_mask(n_mask_chunks))
        M_ALLOCATE(tmp_paralog_vector(n_dims))

        call detect_patterns_kernel(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_DOSAGE_PATTERN,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            work_arr_paralog_subsets = work_arr_paralog_subsets,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_mask = tmp_active_mask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            max_angle = max_angle,&
            gain_gamma = gain_gamma,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine detect_dosage_effect_alloc

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_kernel(module):detect_patterns_kernel]].
    subroutine detect_subfunctionalization(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            work_arr_paralog_subsets,&
            n_paralog_subsets,&
            tmp_active_mask,&
            tmp_paralog_vector,&
            rdi_threshold,&
            paralog_norms,&
            sorted_paralog_norms_perm,&
            tmp_work_array,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
            !! The minimum valid value is `(n_genes + 31) / 32`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size]].
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. Too large a value is capped to the
            !! maximum valid size. The bindings cap it automatically while sizing the work
            !! array; a Fortran caller caps it by calling
            !! [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
            !! Zero is in range and means there is no subset to check -- the sizing routine reports
            !! it whenever the filtered families hold a single gene each. It reports a work array
            !! of zero slots along with it, which this routine does not accept, so a caller that
            !! gets zero back has nothing to detect and should not call here at all.
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(out) :: tmp_active_mask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: tmp_paralog_vector
            !! vector used for pruning subsets
        real(real64), intent(in) :: rdi_threshold
            !! max allowed residual distance from `ancestor`
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_genes), intent(in) :: paralog_norms
            !! euclidean norms of the genes, used for subset pruning (`norm` from `f42_utils` computes them)
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_genes), intent(in) :: sorted_paralog_norms_perm
            !! ascending permutation of the norms, for subset pruning: the smallest norm among the genes that could extend a subset must not fall below the subset's angle to the ancestor
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        real(real64), dimension(n_genes), intent(out) :: tmp_work_array
            !! work array for checking the minimum value after a given index efficiently
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=(n_genes + 31) / 32)
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_dimension_size(n_paralog_subsets, ierr, arg_pos=10_int32)
        call validate_in_range_real(rdi_threshold, ierr, arg_pos=13_int32, min=0.0_real64)
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(paralog_norms, n_genes, ierr, arg_pos=14_int32, min=0.0_real64)
        call validate_all_in_range_int(sorted_paralog_norms_perm, n_genes, ierr, arg_pos=15_int32, min=1_int32, max=n_genes)
        if (is_err(ierr)) return

        call detect_patterns_kernel(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_SUBFUNC_PATTERN,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            work_arr_paralog_subsets = work_arr_paralog_subsets,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_mask = tmp_active_mask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            rdi_threshold = rdi_threshold,&
            paralog_norms = paralog_norms,&
            sorted_paralog_norms_perm = sorted_paralog_norms_perm,&
            tmp_work_array = tmp_work_array,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine detect_subfunctionalization

    !> summary: Allocates its work arrays, then calls [[tox_paralog_analysis_kernel(module):detect_patterns_kernel]].
    subroutine detect_subfunctionalization_alloc(&
            ancestor,&
            genes,&
            n_genes,&
            n_dims,&
            filtered_paralogs_mask,&
            n_mask_chunks,&
            n_results,&
            max_subset_size,&
            work_arr_paralog_subsets,&
            n_paralog_subsets,&
            rdi_threshold,&
            paralog_norms,&
            sorted_paralog_norms_perm,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! number of vectors in `genes`
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `genes`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_genes` genes. Use subroutine `mask_chunk_count` for calculation
            !! The minimum valid value is `(n_genes + 31) / 32`.
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of gene subsets that can be stored in `work_arr_paralog_subsets`.
            !! It is *VERY IMPORTANT* to compute this argument from the `work_array_size` output produced by [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size]].
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_genes), intent(in) :: genes
            !! expression vectors of genes
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with the genes' indices kept by this pattern set to 1, else 0. Build it with the matching `filter_paralogs_by_pattern_*` routine
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked gene subsets. Too large a value is capped to the
            !! maximum valid size. The bindings cap it automatically while sizing the work
            !! array; a Fortran caller caps it by calling
            !! [[tox_paralog_analysis_kernel(module):calc_work_arr_paralog_subsets_size(subroutine)]] first.
            !! Zero is in range and means there is no subset to check -- the sizing routine reports
            !! it whenever the filtered families hold a single gene each. It reports a work array
            !! of zero slots along with it, which this routine does not accept, so a caller that
            !! gets zero back has nothing to detect and should not call here at all.
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0_real64)` and represents the number of chunks
            !! @endnote
        real(real64), intent(in) :: rdi_threshold
            !! max allowed residual distance from `ancestor`
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_genes), intent(in) :: paralog_norms
            !! euclidean norms of the genes, used for subset pruning (`norm` from `f42_utils` computes them)
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_genes), intent(in) :: sorted_paralog_norms_perm
            !! ascending permutation of the norms, for subset pruning: the smallest norm among the genes that could extend a subset must not fall below the subset's angle to the ancestor
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), dimension(:), allocatable :: tmp_active_mask
        real(real64), dimension(:), allocatable :: tmp_paralog_vector
        real(real64), dimension(:), allocatable :: tmp_work_array

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dims, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=6_int32, min=(n_genes + 31) / 32)
        call validate_in_range_int(max_subset_size, ierr, arg_pos=8_int32, min=0_int32)
        call validate_dimension_size(n_paralog_subsets, ierr, arg_pos=10_int32)
        call validate_in_range_real(rdi_threshold, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_all_in_range_real(ancestor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(genes, n_dims * n_genes, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(paralog_norms, n_genes, ierr, arg_pos=12_int32, min=0.0_real64)
        call validate_all_in_range_int(sorted_paralog_norms_perm, n_genes, ierr, arg_pos=13_int32, min=1_int32, max=n_genes)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_active_mask(n_mask_chunks))
        M_ALLOCATE(tmp_paralog_vector(n_dims))
        M_ALLOCATE(tmp_work_array(n_genes))

        call detect_patterns_kernel(&
            ancestor = ancestor,&
            genes = genes,&
            n_genes = n_genes,&
            n_dims = n_dims,&
            pattern_mode = MODE_SUBFUNC_PATTERN,&
            filtered_paralogs_mask = filtered_paralogs_mask,&
            n_mask_chunks = n_mask_chunks,&
            n_results = n_results,&
            max_subset_size = max_subset_size,&
            work_arr_paralog_subsets = work_arr_paralog_subsets,&
            n_paralog_subsets = n_paralog_subsets,&
            tmp_active_mask = tmp_active_mask,&
            tmp_paralog_vector = tmp_paralog_vector,&
            rdi_threshold = rdi_threshold,&
            paralog_norms = paralog_norms,&
            sorted_paralog_norms_perm = sorted_paralog_norms_perm,&
            tmp_work_array = tmp_work_array,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine detect_subfunctionalization_alloc

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_kernel(module):filter_paralogs_by_pattern_kernel]].
    !| This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    subroutine filter_paralogs_by_pattern_dosage_effect(&
            gene_angles,&
            threshold,&
            n_genes,&
            n_families,&
            gene_to_fam,&
            masks,&
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
            !! The minimum valid value is `(n_genes + 31) / 32`.
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
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: masks
            !! bit mask that will have the indices of genes kept by this pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_in_range_real(threshold, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_families, ierr, arg_pos=4_int32, min=1_int32, max=n_genes)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=7_int32, min=(n_genes + 31) / 32)
        call validate_all_in_range_real(gene_angles, n_genes, ierr, arg_pos=1_int32, min=0.0_real64, max=PI)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=5_int32, min=1_int32, max=n_families)
        if (is_err(ierr)) return

        call filter_paralogs_by_pattern_kernel(&
            pattern_mode = MODE_DOSAGE_PATTERN,&
            gene_angles = gene_angles,&
            threshold = threshold,&
            n_genes = n_genes,&
            n_families = n_families,&
            gene_to_fam = gene_to_fam,&
            masks = masks,&
            n_mask_chunks = n_mask_chunks,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine filter_paralogs_by_pattern_dosage_effect

    !> summary: Validates its inputs, then calls [[tox_paralog_analysis_kernel(module):filter_paralogs_by_pattern_kernel]].
    !| This subroutine prefilters the genes for a specific pattern to reduce detection overhead, as less subsets need to be tried.
    subroutine filter_paralogs_by_pattern_subfunctionalization(&
            gene_angles,&
            threshold,&
            n_genes,&
            n_families,&
            gene_to_fam,&
            masks,&
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
            !! The minimum valid value is `(n_genes + 31) / 32`.
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
        integer(int32), dimension(n_mask_chunks, n_families), intent(out) :: masks
            !! bit mask that will have the indices of genes kept by this pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_in_range_real(threshold, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_genes, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_families, ierr, arg_pos=4_int32, min=1_int32, max=n_genes)
        call validate_in_range_int(n_mask_chunks, ierr, arg_pos=7_int32, min=(n_genes + 31) / 32)
        call validate_all_in_range_real(gene_angles, n_genes, ierr, arg_pos=1_int32, min=0.0_real64, max=PI)
        call validate_all_in_range_int(gene_to_fam, n_genes, ierr, arg_pos=5_int32, min=1_int32, max=n_families)
        if (is_err(ierr)) return

        call filter_paralogs_by_pattern_kernel(&
            pattern_mode = MODE_SUBFUNC_PATTERN,&
            gene_angles = gene_angles,&
            threshold = threshold,&
            n_genes = n_genes,&
            n_families = n_families,&
            gene_to_fam = gene_to_fam,&
            masks = masks,&
            n_mask_chunks = n_mask_chunks,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine filter_paralogs_by_pattern_subfunctionalization

end module tox_paralog_analysis
