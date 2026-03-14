#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation
!|
!| This module implements the pipeline to obtain the JSD value from neighborhood residuals obtained from [[tox_data_integration_preprocessing(submodule)]].
module tox_data_integration_jsd
    use safeguard
    use tox_data_integration
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_positive_inf
    use, intrinsic :: iso_c_binding, only: c_ptr
    use f42_heaps, only: top_k_heap_push, bottom_k_heap_push, init_top_k_heap, init_bottom_k_heap
    use f42_random_gsl, only: random_multinomial, create_rng, random_multiv_hypergeom
    use f42_utils, only: clamp, calc_percentile_helper, calc_percentile_rank, is_close, sort_array_heapsort, LOG_2
    use tox_errors, only: ERR_INVALID_INPUT, map_err_arg_pos, set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size, validate_in_range_real, validate_all_in_range_real, validate_in_range_int, validate_all_in_range_int
    implicit none

    integer(int32), parameter :: JOIN_MIN = 0
    integer(int32), parameter :: JOIN_MAX = 1
    integer(int32), parameter :: JOIN_MEDIAN = 2
contains

    !> Computes the shared residual range [-R, R] for the computed residuals from all studies
    pure subroutine determine_shared_residual_range(residuals, residuals_perm, n_residuals, shared_residual_range, ierr, residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        real(real64), dimension(n_residuals), intent(in) :: residuals
            !! Matrix of signed residuals
        integer(int32), dimension(n_residuals), intent(in) :: residuals_perm
            !! Sorting permutation for `residuals`
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_residuals, ierr, arg_pos=3_int32)
        call validate_in_range_real(residual_range_quantile, ierr, min=0.0_real64, max=100.0_real64, arg_pos=6_int32)
        call validate_all_in_range_int(residuals_perm, n_residuals, ierr, min=1_int32, max=n_residuals, arg_pos=2_int32)

        if (is_err(ierr)) return

        call determine_shared_residual_range_helper(residuals, residuals_perm, n_residuals, shared_residual_range, residual_range_quantile)
    end subroutine determine_shared_residual_range

    !> (no input validation) Computes the shared residual range [-R, R] for the computed residuals from all studies
    pure subroutine determine_shared_residual_range_helper(residuals, residuals_perm, n_residuals, shared_residual_range, residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        real(real64), dimension(n_residuals), intent(in) :: residuals
            !! Matrix of signed residuals
        integer(int32), dimension(n_residuals), intent(in) :: residuals_perm
            !! Sorting permutation for `residuals`
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0

        integer(int32) :: i_pool, n_pool, lower_index, left, right
        real(real64) :: actual_quantile, index, fraction, upper_val

        M_DEFAULT_VAL(residual_range_quantile, actual_quantile, 95.0_real64)
        
        n_pool = find_last_non_nan(residuals, residuals_perm, size(residuals, kind=int32))

        shared_residual_range = 0.0_real64
        
        if (n_pool == 0) return

        ! Residual range: Calculate quantile for values
        index = max(1.0_real64, calc_percentile_rank(actual_quantile, n_pool))
        lower_index = floor(index)
        fraction = index - real(lower_index, real64)
        ! As residuals might have negative values, pick from top and bottom until rank is reached
        left = 1
        right = n_pool
        do i_pool = lower_index, n_pool
            if (i_pool == n_pool) upper_val = shared_residual_range
            if (abs(residuals(residuals_perm(left))) > abs(residuals(residuals_perm(right)))) then
                shared_residual_range = abs(residuals(residuals_perm(left)))
                left = left + 1
            else
                shared_residual_range = abs(residuals(residuals_perm(right)))
                right = right - 1
            end if
        end do

        ! Do interpolation only if there are higher values
        if (lower_index /= n_pool) then
            shared_residual_range = shared_residual_range + (upper_val - shared_residual_range) * fraction
        end if
    end subroutine determine_shared_residual_range_helper

    !> Computes the shared residual range [-R, R] for the computed residuals from all studies
    pure subroutine determine_shared_residual_range_alloc(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, shared_residual_range, ierr, residual_range_quantile)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), dimension(max_n_reps_all_studies * max_n_genes_all_studies * n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: n_bins

        n_bins = 1_int32

        call determine_bin_count_and_shared_residual_range_alloc_helper(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, 1_int32, shared_residual_range, n_bins, determine_bin_count=.false., determine_shared_residual_range=.true., ierr=ierr, residual_range_quantile=residual_range_quantile)

        call map_err_arg_pos(ierr, 11_int32, 5_int32)
    end subroutine determine_shared_residual_range_alloc

    !> More efficient shorthand for the pipeline:
    !|
    !| 1. [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]] 
    !| 2. [[tox_data_integration(module):estimate_bin_count_alloc(interface)]] 
    !|
    !| It is faster because it sorts the residuals only once, while the separate calls do it twice. 
    pure subroutine determine_bin_count_and_shared_residual_range_alloc(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_neighbors, shared_residual_range, n_bins, ierr, residual_range_quantile)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), dimension(max_n_reps_all_studies * max_n_genes_all_studies * n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=4_int32)
        if (is_err(ierr)) return

        call determine_bin_count_and_shared_residual_range_alloc_helper(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_neighbors, shared_residual_range, n_bins, determine_bin_count=.true., determine_shared_residual_range=.true., ierr=ierr, residual_range_quantile=residual_range_quantile)
        
        call map_err_arg_pos(ierr, 11_int32, 9_int32)
    end subroutine determine_bin_count_and_shared_residual_range_alloc

    !> (with input validation) Root helper for:
    !|
    !| - [[tox_data_integration(module):determine_bin_count_and_shared_residual_range_alloc(interface)]]
    !| - [[tox_data_integration(module):estimate_bin_count_alloc(interface)]]
    !| - [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
    !|
    pure subroutine determine_bin_count_and_shared_residual_range_alloc_helper(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_neighbors, shared_residual_range, n_bins, determine_bin_count, determine_shared_residual_range, ierr, residual_range_quantile)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), dimension(max_n_reps_all_studies * max_n_genes_all_studies * n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        real(real64), intent(inout) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(inout) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        integer(int32), intent(out) :: ierr
            !! Error code
        logical, intent(in) :: determine_bin_count
            !! Should bin count be determined? If not, `n_bins` remains unchanged
        logical, intent(in) :: determine_shared_residual_range
            !! Should shared residual range be determined? If not, `shared_residual_range` is taken as input

        integer(int32) :: i_residual, n_residuals
        integer(int32), dimension(:), allocatable :: residuals_perm

        call set_ok(ierr)

        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=2_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=4_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_real(residual_range_quantile, ierr, min=0.0_real64, max=100.0_real64, arg_pos=11_int32)
        if (.not. determine_shared_residual_range) call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64, arg_pos=6_int32)

        if (is_err(ierr)) return

        if (determine_bin_count .or. determine_shared_residual_range) then
            n_residuals = size(residuals, kind=int32)
            M_ALLOCATE(residuals_perm(n_residuals))

            ! initialize permutation
            do concurrent (i_residual = 1:n_residuals) shared(residuals_perm)
                residuals_perm(i_residual) = i_residual
            end do

            call sort_array_heapsort(residuals, residuals_perm)

            if (determine_shared_residual_range) then
                call determine_shared_residual_range_helper(residuals, residuals_perm, n_residuals, shared_residual_range, residual_range_quantile)
            end if

            if (determine_bin_count) then
                call estimate_bin_count_helper(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, n_neighbors, shared_residual_range, n_bins)
            end if
        end if
    end subroutine determine_bin_count_and_shared_residual_range_alloc_helper

    !> Estimates the number of histogram bins for [[tox_data_integration(module):build_residual_histograms(interface)]],
    !| using the maximum value returned by the Sturges rule
    !| \[ \texttt{sturges_bins} = 1 + \lfloor \frac{\ln\left(\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}\right)}{\ln\left(2\right)} \rceil \]
    !| and the Freedman-Diaconis rule
    !| \[ \texttt{freed_diac_bins} = 2 \cdot \frac{\operatorname{IQR}(\texttt{residuals})}{\sqrt[3]{\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}}} \]
    !| so finally
    !| \[\texttt{n_bins} = \max\left(\texttt{sturges_bins}, \texttt{freed_diac_bins}\right)\]
    pure subroutine estimate_bin_count_alloc(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_neighbors, shared_residual_range, n_bins, ierr)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), dimension(max_n_reps_all_studies * max_n_genes_all_studies * n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64) :: tmp_shared_residual_range

        tmp_shared_residual_range = shared_residual_range

        call determine_bin_count_and_shared_residual_range_alloc_helper(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_neighbors, tmp_shared_residual_range, n_bins, determine_bin_count=.true., determine_shared_residual_range=.false., ierr=ierr)
    end subroutine estimate_bin_count_alloc

    !> Estimates the number of histogram bins for [[tox_data_integration(module):build_residual_histograms(interface)]],
    !| using the maximum value returned by the Sturges rule
    !| \[ \texttt{sturges_bins} = 1 + \lfloor \frac{\ln\left(\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}\right)}{\ln\left(2\right)} \rceil \]
    !| and the Freedman-Diaconis rule
    !| \[ \texttt{freed_diac_bins} = 2 \cdot \frac{\operatorname{IQR}(\texttt{residuals})}{\sqrt[3]{\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}}} \]
    !| so finally
    !| \[\texttt{n_bins} = \max\left(\texttt{sturges_bins}, \texttt{freed_diac_bins}\right)\]
    pure subroutine estimate_bin_count(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, n_neighbors, shared_residual_range, n_bins, ierr)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), dimension(n_residuals), intent(in) :: residuals
            !! Matrix of signed residuals
        integer(int32), dimension(n_residuals), intent(in) :: residuals_perm
            !! Sorting permutation for `residuals`
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_residuals, ierr, arg_pos=3_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64, arg_pos=6_int32)
        call validate_all_in_range_int(residuals_perm, n_residuals, ierr, min=1_int32, max=n_residuals, arg_pos=2_int32)

        if (is_err(ierr)) return

        call estimate_bin_count_helper(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, n_neighbors, shared_residual_range, n_bins)
    end subroutine estimate_bin_count

    !> (no input validation) Estimates the number of histogram bins for [[tox_data_integration(module):build_residual_histograms(interface)]],
    !| using the maximum value returned by the Sturges rule
    !| \[ \texttt{sturges_bins} = 1 + \lfloor \frac{\ln\left(\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}\right)}{\ln\left(2\right)} \rceil \]
    !| and the Freedman-Diaconis rule
    !| \[ \texttt{freed_diac_bins} = 2 \cdot \frac{\operatorname{IQR}(\texttt{residuals})}{\sqrt[3]{\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}}} \]
    !| so finally
    !| \[\texttt{n_bins} = \max\left(\texttt{sturges_bins}, \texttt{freed_diac_bins}\right)\]
    pure subroutine estimate_bin_count_helper(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, n_neighbors, shared_residual_range, n_bins)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), dimension(n_residuals), intent(in) :: residuals
            !! Matrix of signed residuals
        integer(int32), dimension(n_residuals), intent(in) :: residuals_perm
            !! Sorting permutation for `residuals`
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for

        integer(int32) :: n_pool
        real(real64) :: half_bin_width, quartile_25, quartile_75, n_reps_neighborhood

        n_pool = find_last_non_nan(residuals, residuals_perm, size(residuals, kind=int32))
        if (n_pool == 0) then
            n_bins = 0
        else
            n_reps_neighborhood = real(max_n_reps_all_studies * n_neighbors, kind=real64)

            ! Estimate n_bins
            ! Sturges
            n_bins = 1 + nint(log(n_reps_neighborhood) / LOG_2)

            ! Freedman-Diaconis
            call calc_percentile_helper(residuals, residuals_perm(:n_pool), 25.0_real64, quartile_25)
            call calc_percentile_helper(residuals, residuals_perm(:n_pool), 75.0_real64, quartile_75)

            ! bin width as defined by Freedman-Diaconis, but without doubling the value.
            ! With doubling the value, the bin count would be calculated as `2*shared_residual_range/(2*bin_width)` -> the 2 is unnecessary
            half_bin_width = (quartile_75 - quartile_25) / (n_reps_neighborhood ** (1.0_real64 / 3.0_real64))
            if (.not. is_close(half_bin_width, 0.0_real64)) then
                n_bins = max(n_bins, nint(shared_residual_range / half_bin_width))
            end if
        end if
    end subroutine estimate_bin_count_helper

    subroutine js_comp_test(residuals, gene_means, gene_means_perms, gene_means_perm_all, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_bins, shared_residual_range, n_points_candidates, n_point_counts, max_n_points_candidate, n_neighbors_candidates, n_neighbor_counts, max_n_neighbors_candidate, neighborhood_residuals, neighborhood_ranges, x_star, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts, tmp_mean_pmf_counts, weights, js_divergences, global_js_divergence, p_values, confidence_interval, best_confidence_interval, join_method, n_bootstraps, n_bootstrapping_top_k_jsds, bootstrapping_top_k_jsds)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_points_candidate
            !! Number of reference points in the studies
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Neighborhood size
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: neighborhood_residuals
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note 
            !! All indices in range `1<=idx<=max(max_n_neighbors_candidate, n_genes_S)`. So in case `n_genes_S` is lower than `max_n_neighbors_candidate`,
            !! remaining indices will be filled with the ones from `n_genes_S+1...max_n_neighbors_candidate`
            !! @endnote
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: neighborhood_ranges
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is [1, 1]
        real(real64), dimension(max_n_points_candidate), intent(out) :: x_star
            !! Mean-expression reference points
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
        integer(int32), intent(in) :: n_bootstraps
        integer(int32), intent(in) :: n_point_counts
        integer(int32), intent(in) :: n_neighbor_counts
        integer(int32), dimension(n_point_counts), intent(in) :: n_points_candidates
        integer(int32), dimension(n_neighbor_counts), intent(in) :: n_neighbors_candidates
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means_perm_all
        real(real64), dimension(n_bins, max_n_points_candidate, n_studies), intent(out) :: pmfs
        integer(int32), dimension(n_bins, max_n_points_candidate, n_studies), intent(out) :: counts
        integer(int32), dimension(n_bins, max_n_points_candidate), intent(out) :: mean_pmf_counts
        integer(int32), dimension(n_bins, max_n_points_candidate), intent(out) :: tmp_mean_pmf_counts
        real(real64), dimension(n_bins, max_n_points_candidate), intent(out) :: mean_pmf
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: included_n_reps
        integer(int32), dimension(max_n_points_candidate), intent(out) :: mean_pmf_included_n_reps
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: js_divergences
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: weights
        real(real64), dimension(2, n_studies), intent(out), target :: confidence_interval
        real(real64), dimension(2, n_studies), intent(out) :: best_confidence_interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: bootstrapping_top_k_jsds
        real(real64), dimension(n_studies), intent(out) :: p_values
        integer(int32), intent(in) :: join_method

        integer(int32), parameter :: min_count_per_mean_bin = 5_int32
        integer(int32), parameter :: n_permutations = 1000_int32

        integer(int32) :: n_gene_means, n_pool, n_points, n_neighbors, i_study, i_neighbor_count, i_point_count, best_params_CI_i_point_count, best_params_CI_i_neighbor_count, best_params_exceeded_CI_overlap
        logical :: plateau_found, all_have_min_overlap
        type(c_ptr) :: rng
        real(real64), dimension(:), pointer :: gene_means_flat, confidence_interval_flat
        integer(int32), dimension(:), pointer :: gene_means_perm_all_flat

        rng = create_rng()

        n_gene_means = size(gene_means, kind=int32)
        gene_means_flat(1:n_gene_means) => gene_means
        gene_means_perm_all_flat(1:n_gene_means) => gene_means_perm_all
        n_pool = find_last_non_nan(gene_means_flat, gene_means_perm_all_flat, n_gene_means)

        best_confidence_interval = M_POS_INF
        best_params_CI_i_point_count = 1
        best_params_CI_i_neighbor_count = 1
        best_params_exceeded_CI_overlap = 0
        plateau_found = .false.
        do i_point_count = 1, n_point_counts
            n_points = n_points_candidates(i_point_count)

            call pool_means_n_pool_input_helper(gene_means_flat, gene_means_perm_all_flat, n_gene_means, n_points, n_pool, x_star)

            do i_neighbor_count = 1, n_neighbor_counts
                n_neighbors = n_neighbors_candidates(i_neighbor_count)

                all_have_min_overlap = .true.
                do concurrent (i_study = 1:n_studies) shared(n_points, x_star, max_n_genes_all_studies, gene_means, gene_means_perms, neighborhood_residuals, neighborhood_ranges, n_neighbors, shared_residual_range, n_bins, counts, pmfs, included_n_reps)
                    call construct_neighborhoods_helper(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), gene_means_perms(:, i_study), neighborhood_residuals(:, :, i_study), neighborhood_ranges(:, :, i_study), n_neighbors)
                    if (test_neighborhood_overlaps_helper(neighborhood_ranges(:, :, i_study), n_points)) then
                        call build_residual_histograms_helper(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study))
                    else
                        all_have_min_overlap = .false.
                    end if
                end do

                if (all_have_min_overlap) then
                    call create_mean_pmf_helper(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts)
                    if (.not. test_mean_pmf_min_counts_helper(mean_pmf_counts, n_bins, n_points, min_count_per_mean_bin)) cycle

                    do concurrent (i_study = 1:n_studies) shared(pmfs, mean_pmf, n_points, n_bins, js_divergences, included_n_reps, mean_pmf_included_n_reps, weights)
                        call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
                        call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, confidence_interval(1, i_study), weights(:, i_study))
                        confidence_interval(2, i_study) = confidence_interval(1, i_study)
                    end do

                    call bootstrap_histogram_helper(n_bootstraps, counts, pmfs, included_n_reps, n_bins, n_points, n_studies, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps, confidence_interval, js_divergences, weights, global_js_divergence, bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, rng)
                    call check_plateau_condition(confidence_interval, best_confidence_interval, n_studies, best_params_CI_i_point_count, best_params_CI_i_neighbor_count, best_params_exceeded_CI_overlap, i_point_count, i_neighbor_count, join_method, plateau_found)
                    if (plateau_found) exit
                end if
            end do
            if (plateau_found) exit
        end do
    
        n_points = n_points_candidates(best_params_CI_i_point_count)
        n_neighbors = n_neighbors_candidates(best_params_CI_i_neighbor_count)

        call pool_means_helper(gene_means_flat, gene_means_perm_all_flat, n_gene_means, n_points, n_pool, x_star)

        do concurrent (i_study = 1:n_studies) shared(n_points, x_star, max_n_genes_all_studies, gene_means, gene_means_perms, neighborhood_residuals, neighborhood_ranges, n_neighbors)
            call construct_neighborhoods_helper(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), gene_means_perms(:, i_study), neighborhood_residuals(:, :, i_study), neighborhood_ranges(:, :, i_study), n_neighbors)
            call build_residual_histograms_helper(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study))
        end do
        call create_mean_pmf_helper(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts)

        do concurrent (i_study = 1:n_studies) shared(pmfs, mean_pmf, n_points, n_bins, js_divergences, included_n_reps, mean_pmf_included_n_reps, global_js_divergence, weights)
            call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, global_js_divergence(i_study), weights(:, i_study))
        end do

        confidence_interval_flat(1:size(confidence_interval, kind=int32)) => confidence_interval
        call gjct_permutation_test_helper(n_permutations, p_values, counts, pmfs, included_n_reps, n_bins, n_points, n_studies, mean_pmf_counts, tmp_mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, global_js_divergence, js_divergences, weights, confidence_interval_flat(1:n_studies), rng)

        do concurrent (i_study = 1:n_studies) shared(counts, pmfs, included_n_reps, n_bins, n_points)
            call calc_pmf_helper(counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, global_js_divergence(i_study), weights(:, i_study))
        end do
    end subroutine js_comp_test

    pure subroutine check_plateau_condition(confidence_interval, best_confidence_interval, n_studies, best_params_CI_i_point_count, best_params_CI_i_neighbor_count, best_params_exceeded_CI_overlap, i_point_count, i_neighbor_count, join_method, plateau_found)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(2, n_studies), intent(in) :: confidence_interval
        real(real64), dimension(2, n_studies), intent(inout) :: best_confidence_interval
        integer(int32), intent(inout) :: best_params_CI_i_point_count
        integer(int32), intent(inout) :: best_params_CI_i_neighbor_count
        integer(int32), intent(inout) :: best_params_exceeded_CI_overlap
        integer(int32), intent(in) :: i_point_count
        integer(int32), intent(in) :: i_neighbor_count
        integer(int32), intent(in) :: join_method
        logical, intent(out) :: plateau_found

        real(real64), parameter :: succeeding_overlap = 0.9_real64

        integer(int32) :: exceeds_min_CI_overlap, i_study

        exceeds_min_CI_overlap = 0_int32
        do concurrent (i_study = 1:n_studies) shared(confidence_interval, best_confidence_interval) reduce(+:exceeds_min_CI_overlap)
            if (compute_relative_overlap_helper(&
                    confidence_interval(1, i_study), confidence_interval(2, i_study),&
                    best_confidence_interval(1, i_study), best_confidence_interval(2, i_study)&
                ) > succeeding_overlap)&
            then
                exceeds_min_CI_overlap = exceeds_min_CI_overlap + 1
            end if
        end do

        ! If the parameter set is not at least as good as the previous, exit and use previous pair as best
        if (exceeds_min_CI_overlap < best_params_exceeded_CI_overlap) then
            plateau_found = .true.
        else
            best_params_CI_i_point_count = i_point_count
            best_params_CI_i_neighbor_count = i_neighbor_count
            best_params_exceeded_CI_overlap = exceeds_min_CI_overlap
            best_confidence_interval = confidence_interval

            select case (join_method)
                case (JOIN_MIN)
                    ! If all overlaps exceed the threshold, there is no lower one that doesn't
                    plateau_found = exceeds_min_CI_overlap == n_studies
                case (JOIN_MAX)
                    ! If at least one overlap exceeds the threshold, the max overlap will do either
                    plateau_found = exceeds_min_CI_overlap > 0
                case (JOIN_MEDIAN)
                    ! If 50% of the studies' overlaps exceed the threshold, the median will do either
                    plateau_found = exceeds_min_CI_overlap >= n_studies / 2
            end select
        end if
    end subroutine check_plateau_condition

    subroutine gjct_permutation_test_helper(n_permutations, p_values, counts, pmfs, included_n_reps, n_bins, n_points, n_studies, mean_pmf_counts, tmp_mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, global_jsd_observed, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, rng)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_permutations
            !! Number of bootstraps to perform
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: pmfs
        integer(int32), dimension(n_bins, n_points, n_studies), intent(out) :: counts
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_mean_pmf_counts
        real(real64), dimension(n_bins, n_points), intent(in) :: mean_pmf
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
        real(real64), dimension(n_studies), intent(out) :: p_values
        type(c_ptr), intent(in) :: rng

        integer(int32) :: i_permutation, i_point, i_study

        p_values = 0.0_real64
        do i_permutation = 1, n_permutations
            tmp_mean_pmf_counts = mean_pmf_counts
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multiv_hypergeom(rng, tmp_mean_pmf_counts(:, i_point), n_bins, sum(tmp_mean_pmf_counts(:, i_point)), included_n_reps(i_point, i_study), counts(:, i_point, i_study))
                end do
                call calc_pmf_helper(counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            end do

            do concurrent (i_study = 1:n_studies) shared(p_values, pmfs, mean_pmf, n_points, n_bins, tmp_js_divergences, included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, tmp_weights, global_jsd_observed)
                call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_helper(tmp_js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                if (tmp_global_js_divergence(i_study) >= global_jsd_observed(i_study)) then
                    p_values(i_study) = anint(p_values(i_study) + 1.0_real64)
                end if
            end do
        end do

        if (n_permutations == 0) return
        do concurrent (i_study = 1:n_studies) shared(p_values, n_permutations)
            p_values(i_study) = p_values(i_study) / real(n_permutations, real64)
        end do
    end subroutine gjct_permutation_test_helper

    subroutine bootstrap_histogram_helper(n_bootstraps, counts, pmfs, included_n_reps, n_bins, n_points, n_studies, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps, confidence_interval, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, rng)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: bootstrapping_top_k_jsds
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: pmfs
        integer(int32), dimension(n_bins, n_points, n_studies), intent(out) :: counts
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
        type(c_ptr), intent(in) :: rng

        integer(int32) :: i_bootstrap, i_point, i_study

        do concurrent (i_study = 1:n_studies) shared(n_bootstrapping_top_k_jsds, bootstrapping_top_k_jsds, confidence_interval)
            call init_bottom_k_heap(bootstrapping_top_k_jsds(:, 1, i_study), n_bootstrapping_top_k_jsds)
            call bottom_k_heap_push(bootstrapping_top_k_jsds(:, 1, i_study), n_bootstrapping_top_k_jsds, confidence_interval(1, i_study))
            call init_top_k_heap(bootstrapping_top_k_jsds(:, 2, i_study), n_bootstrapping_top_k_jsds)
            call top_k_heap_push(bootstrapping_top_k_jsds(:, 2, i_study), n_bootstrapping_top_k_jsds, confidence_interval(2, i_study))
        end do

        do i_bootstrap = 1, n_bootstraps
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multinomial(rng, mean_pmf_counts(:, i_point), n_bins, mean_pmf_included_n_reps(i_point), included_n_reps(i_point, i_study), counts(:, i_point, i_study))
                end do
                call calc_pmf_helper(counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            end do

            call create_mean_pmf_only_helper(pmfs, n_bins, n_points, n_studies, mean_pmf)

            do concurrent (i_study = 1:n_studies)
                call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_helper(tmp_js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                call bottom_k_heap_push(bootstrapping_top_k_jsds(:, 1, i_study), n_bootstrapping_top_k_jsds, tmp_global_js_divergence(i_study))
                call top_k_heap_push(bootstrapping_top_k_jsds(:, 2, i_study), n_bootstrapping_top_k_jsds, tmp_global_js_divergence(i_study))
            end do
        end do

        do concurrent (i_study = 1:n_studies) shared(n_bootstrapping_top_k_jsds, bootstrapping_top_k_jsds, confidence_interval)
            confidence_interval(1, i_study) = bootstrapping_top_k_jsds(1, 1, i_study)
            confidence_interval(2, i_study) = bootstrapping_top_k_jsds(1, 2, i_study)
        end do
    end subroutine bootstrap_histogram_helper

    pure subroutine calc_pmf_helper(counts, pmf, included_n_reps, n_bins, n_points)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(n_bins, n_points), intent(inout) :: pmf
        integer(int32), dimension(n_bins, n_points), intent(in) :: counts
        integer(int32), dimension(n_points), intent(in) :: included_n_reps

        integer(int32) :: i_point, i_bin

        do concurrent (i_point = 1:n_points)
            do concurrent (i_bin = 1:n_bins) shared(pmf, i_point, included_n_reps, counts)
                if (included_n_reps(i_point) == 0) then
                    pmf(i_bin, i_point) = 0.0_real64
                else
                    pmf(i_bin, i_point) = real(counts(i_bin, i_point), real64) / real(included_n_reps(i_point), real64)
                end if
            end do
        end do
            
    end subroutine calc_pmf_helper

    pure logical function test_mean_pmf_min_counts_helper(mean_pmf_counts, n_bins, n_points, min) result(all_bins_have_min_count)
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: min
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts

        integer(int32) :: i_point, i_bin

        all_bins_have_min_count = .true.
        do concurrent (i_point = 1:n_points)
            do concurrent (i_bin = 1:n_bins) shared(mean_pmf_counts, min, i_point, all_bins_have_min_count)
                if (mean_pmf_counts(i_bin, i_point) < min) then
                    all_bins_have_min_count = .false.
                end if
            end do
        end do
    
    end function test_mean_pmf_min_counts_helper

    pure subroutine create_mean_pmf_helper(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(n_bins, n_points, n_studies), intent(in) :: pmfs
        integer(int32), dimension(n_bins, n_points, n_studies), intent(in) :: counts
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64
        mean_pmf_included_n_reps = 0_int32
        mean_pmf_counts = 0_int32

        do i_study = 1, n_studies
            do i_point = 1, n_points
                do i_bin = 1, n_bins
                    mean_pmf(i_bin, i_point) = mean_pmf(i_bin, i_point) + pmfs(i_bin, i_point, i_study) / real(n_studies, real64)
                    mean_pmf_counts(i_bin, i_point) = mean_pmf_counts(i_bin, i_point) + counts(i_bin, i_point, i_study)
                end do
                mean_pmf_included_n_reps(i_point) = mean_pmf_included_n_reps(i_point) + included_n_reps(i_point, i_study)
            end do
        end do
    
    end subroutine create_mean_pmf_helper

    pure subroutine create_mean_pmf_only_helper(pmfs, n_bins, n_points, n_studies, mean_pmf)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(n_bins, n_points, n_studies), intent(in) :: pmfs
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64

        do i_study = 1, n_studies
            do i_point = 1, n_points
                do i_bin = 1, n_bins
                    mean_pmf(i_bin, i_point) = mean_pmf(i_bin, i_point) + pmfs(i_bin, i_point, i_study) / real(n_studies, real64)
                end do
            end do
        end do
    end subroutine create_mean_pmf_only_helper

    pure logical function test_neighborhood_overlaps_helper(neighborhood_range, n_points) result(all_have_min_overlap)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), dimension(2, n_points), intent(in) :: neighborhood_range
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is [1, 1]

        integer(int32) :: i_point
        real(real64) :: overlap

        all_have_min_overlap = .true.
        do concurrent (i_point = 1:n_points - 1) local(overlap) shared(neighborhood_range, all_have_min_overlap)
            overlap = compute_relative_overlap_helper(&
                real(neighborhood_range(1, i_point), real64),&
                real(neighborhood_range(2, i_point), real64),&
                real(neighborhood_range(1, i_point + 1), real64),&
                real(neighborhood_range(2, i_point + 1), real64)&
            )
            if (overlap < 0.1_real64) all_have_min_overlap = .false.
        end do
    
    end function test_neighborhood_overlaps_helper

    pure real(real64) function compute_relative_overlap_helper(a_min, a_max, b_min, b_max) result(overlap_percent)
        real(real64), intent(in) :: a_min
        real(real64), intent(in) :: a_max
        real(real64), intent(in) :: b_min
        real(real64), intent(in) :: b_max
        if (a_min < b_min) then
            overlap_percent = (min(a_max, b_max) - b_min) / (b_max - a_min)
        else
            overlap_percent = (min(b_max, b_max) - a_min) / (a_max - b_min)
        end if
    end function compute_relative_overlap_helper

    !> Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts, pmf, included_n_reps, ierr, neighbor_mask)
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(in) :: residuals
            !! Matrix of signed residuals for one study
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_bins, n_points), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_bins, n_points), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        integer(int32), intent(out) :: ierr
            !! Error code
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        call set_ok(ierr)

        call validate_dimension_size(max_n_reps_all_studies, ierr)
        call validate_dimension_size(max_n_genes_all_studies, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_bins, ierr)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)
        call validate_all_in_range_int(neighborhood_residuals, size(neighborhood_residuals, kind=int32), ierr, min=1_int32, max=max_n_genes_all_studies)

        if (is_err(ierr)) return

        call build_residual_histograms_helper(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
    end subroutine build_residual_histograms

    !> (no input validation) Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms_helper(neighborhood_residuals, n_neighbors, n_points, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(in) :: residuals
            !! Matrix of signed residuals for one study
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_bins, n_points), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_bins, n_points), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        real(real64) :: bin_width, clamped_residual, replicate
        integer(int32) :: bin_idx, i_neighbor, i_rep, included_reps, i_point
        logical :: filter_neighbors

        bin_width = 2.0_real64 * shared_residual_range / real(n_bins, real64)
        counts = 0_int32
        pmf = 0.0_real64

        filter_neighbors = present(neighbor_mask)

        ! 1. assign the bins to the residuals (increase the respective count)
        ! outer loop cannot be concurrent, as counts of same residuals and bins but different neighbors might be changed at the same time
        do concurrent (i_point = 1:n_points) local(included_reps) shared(included_n_reps)
            included_reps = 0_int32
            do concurrent (i_neighbor = 1:n_neighbors) &
                    local(i_rep, clamped_residual, bin_idx, replicate) &
                    shared(filter_neighbors, max_n_reps_all_studies, counts, neighborhood_residuals, shared_residual_range, bin_width) &
                    reduce(+:included_reps)
                ! Exclude neighbor if desired
                if (filter_neighbors) then
                    if (.not. neighbor_mask(i_neighbor, i_point)) cycle
                end if

                ! Count non-NaNs and assign the to a bin
                do i_rep = 1, max_n_reps_all_studies
                    replicate = residuals(i_rep, neighborhood_residuals(i_neighbor, i_point))
                    if (.not. ieee_is_nan(replicate)) then
                        ! clamp residual to histogram range
                        clamped_residual = clamp(replicate, min_val=-shared_residual_range, max_val=shared_residual_range)

                        ! assign bin to residual
                        bin_idx = min(n_bins, int( (clamped_residual + shared_residual_range) / bin_width ) + 1)
                        counts(i_point, bin_idx) = counts(i_point, bin_idx) + 1

                        included_reps = included_reps + 1
                    end if
                end do
            end do
            included_n_reps(i_point) = included_reps
        end do

        call calc_pmf_helper(counts, pmf, included_n_reps, n_bins, n_points)
    end subroutine build_residual_histograms_helper

    !> Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point(pmf_S1, pmf_S2, n_points, n_bins, js_divergences, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(pmf_S1, size(pmf_S1, kind=int32), ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pmf_S2, size(pmf_S2, kind=int32), ierr, arg_pos=2_int32, min=0.0_real64, max=1.0_real64)

        if (is_err(ierr)) return

        call compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
    end subroutine compute_divergence_per_reference_point

    !> (no input validation) Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
        real(real64), dimension(n_bins, n_points), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point

        real(real64) :: S_mean, s1_val, s2_val
        integer(int32) :: i_bin, i_point

        js_divergences = 0.0_real64

        ! 1. compute the Knullback-Leibler (KL) divergences
        ! Note that the J-S divergence is defined as `0.5 * KL_S1 + 0.5 * KL_S2`, equivalent to `0.5 * (KL_S1 + KL_S2)`.
        ! Thus, instead of computing KL_S* independently, it accumulates directly in the `js_divergences` output.
        ! Another thing, switching the loops would enable both to run concurrently and the 0.5*js_divergences step could be done in one go,
        ! but cache locality still beats that, except for the case of thousands of neighbors, which might not be the common case.
        do i_bin = 1, n_bins
            do concurrent (i_point = 1:n_points) local(s1_val, s2_val, S_mean) shared(i_bin, pmf_S1, pmf_S2, js_divergences)
                s1_val = pmf_S1(i_bin, i_point)
                s2_val = pmf_S2(i_bin, i_point)
                S_mean = 0.5_real64 * (s1_val + s2_val)

                if (.not. is_close(S_mean, 0.0_real64)) then
                    if (s1_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s1_val * log(s1_val / S_mean)
                    end if

                    if (s2_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s2_val * log(s2_val / S_mean)
                    end if
                end if
            end do
        end do

        ! 2. Compute the js_divergences
        do concurrent (i_point = 1:n_points) shared(js_divergences)
            js_divergences(i_point) = (0.5_real64 * js_divergences(i_point)) / LOG_2  ! Div by ln2 to have a 0-1 scaling instead of 0-ln2
        end do
    end subroutine compute_divergence_per_reference_point_helper

    !> Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    pure subroutine compute_weighted_global_divergence(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(js_divergences, size(js_divergences, kind=int32), ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(included_n_reps_S1, size(included_n_reps_S1, kind=int32), ierr, arg_pos=3_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps_S2, size(included_n_reps_S2, kind=int32), ierr, arg_pos=4_int32, min=0_int32)

        if (is_err(ierr)) return

        call compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
    end subroutine compute_weighted_global_divergence

    !> (no input validation) Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    pure subroutine compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

        integer(int32) :: i_point, included_reps
        real(real64) :: total_sample_count

        global_js_divergence = 0.0_real64

        total_sample_count = real(sum(included_n_reps_S1) + sum(included_n_reps_S2), real64)

        if (is_close(total_sample_count, 0.0_real64)) then
            weights = 0.0_real64
        else
            ! Calculate the global Jensen-Shannon divergence
            do concurrent (i_point = 1:n_points) local(included_reps) shared(weights, included_n_reps_S1, included_n_reps_S2, total_sample_count) reduce(+:global_js_divergence)
                included_reps = included_n_reps_S1(i_point) + included_n_reps_S2(i_point)

                weights(i_point) = real(included_reps, real64) / total_sample_count

                global_js_divergence = global_js_divergence + weights(i_point) * js_divergences(i_point)
            end do
        end if
    end subroutine compute_weighted_global_divergence_helper
end module tox_data_integration_jsd