#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSCompTest
!|
!| This module implements the pipeline to obtain the JSD value from neighborhood residuals obtained from [[tox_data_integration_preprocessing(module)]].
module tox_data_integration_js_comp_test
    use safeguard
    use tox_data_integration_preprocessing, only: find_last_non_nan, pool_means_n_pool_input_helper, construct_neighborhoods_helper, pool_means_helper
    use tox_data_integration_jsd, only: build_residual_histograms_helper, calc_pmf_helper, compute_weighted_global_divergence_helper, compute_divergence_per_reference_point_helper
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_positive_inf
    use, intrinsic :: iso_c_binding, only: c_ptr
    use f42_heaps, only: top_k_heap_push, bottom_k_heap_push, init_top_k_heap, init_bottom_k_heap
    use f42_random_gsl, only: random_multinomial, create_rng, random_multiv_hypergeom, reset_rng
    use f42_utils, only: calc_percentile_helper, calc_percentile_rank, is_close, sort_array_heapsort, LOG_2
    use tox_errors, only: map_err_arg_pos, set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size, validate_in_range_real, validate_in_range_int, validate_all_in_range_int
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
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
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
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
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
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
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
            n_bins = 1 + nint(log(n_reps_neighborhood) / LOG_2, kind=int32)

            ! Freedman-Diaconis
            call calc_percentile_helper(residuals, residuals_perm(:n_pool), 25.0_real64, quartile_25)
            call calc_percentile_helper(residuals, residuals_perm(:n_pool), 75.0_real64, quartile_75)

            ! bin width as defined by Freedman-Diaconis, but without doubling the value.
            ! With doubling the value, the bin count would be calculated as `2*shared_residual_range/(2*bin_width)` -> the 2 is unnecessary
            half_bin_width = (quartile_75 - quartile_25) / (n_reps_neighborhood ** (1.0_real64 / 3.0_real64))
            if (.not. is_close(half_bin_width, 0.0_real64)) then
                n_bins = max(n_bins, nint(shared_residual_range / half_bin_width, kind=int32))
            end if
        end if
    end subroutine estimate_bin_count_helper

    subroutine js_comp_test(&
            n_points_candidates, n_point_counts, max_n_points_candidate, n_neighbors_candidates, n_neighbor_counts, max_n_neighbors_candidate, n_points, n_neighbors,&
            shared_residual_range, n_bins, gene_means, max_n_genes_all_studies, n_studies, gene_means_perms, gene_means_perm_all,&
            residuals, max_n_reps_all_studies, neighborhood_ranges, x_star, neighborhood_residuals,&
            pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
            js_divergences, weights, global_js_divergence, p_values, best_candidate_pair_confidence_interval, join_method,&
            tmp_pmf_counts, tmp_mean_pmf_counts, tmp_confidence_interval, tmp_bootstrapping_top_k_jsds, n_bootstraps, n_bootstrapping_top_k_jsds,&
            min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, n_permutations, random_seed&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of elements in the top/bottom ends of bootstrapped values, e.g. for 95% of the values as `max(1, floor(0.025 * n_bootstraps))`
        integer(int32), intent(in) :: n_point_counts
            !! Number of candidates for `n_points`
        integer(int32), intent(in) :: n_neighbor_counts
            !! Number of candidates for `n_neighbors`
        integer(int32), dimension(n_point_counts), intent(in) :: n_points_candidates
            !! Candidates for `n_points` to test
        integer(int32), intent(in) :: max_n_points_candidate
            !! Number of reference points in the studies
        integer(int32), dimension(n_neighbor_counts), intent(in) :: n_neighbors_candidates
            !! Candidates for `n_neighbors` to test
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Neighborhood size
        integer(int32), intent(out) :: n_points
            !! The finally chosen candidate from `n_points_candidates`
        integer(int32), intent(out) :: n_neighbors
            !! The finally chosen candidate from `n_neighbors_candidates`
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study permutation vectors to sort the `gene_means`
        integer(int32), dimension(max_n_genes_all_studies * n_studies), intent(in) :: gene_means_perm_all
            !! Permutation vector to sort the flattened `gene_means`
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: neighborhood_ranges
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        real(real64), dimension(max_n_points_candidate), intent(out) :: x_star
            !! Mean-expression reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: neighborhood_residuals
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note 
            !! All indices in range `1<=idx<=max(max_n_neighbors_candidate, n_genes_S)`. So in case `n_genes_S` is lower than `max_n_neighbors_candidate`,
            !! remaining indices will be filled with the ones from `n_genes_S+1...max_n_neighbors_candidate`
            !! @endnote
        real(real64), dimension(n_bins, max_n_points_candidate, n_studies), intent(out) :: pmfs
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_bins, max_n_points_candidate, n_studies), intent(out) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
        real(real64), dimension(n_bins, max_n_points_candidate), intent(out) :: mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
        integer(int32), dimension(n_bins, max_n_points_candidate), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`
        integer(int32), dimension(max_n_points_candidate), intent(out) :: mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: js_divergences
            !! The per-reference-point Jensen-Shannon-Divergence -> finally `js_divergences(1:n_points, :)`
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: weights
            !! The per-reference-point weights for the Jensen-Shannon-Divergence -> finally `weights(1:n_points, :)`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! The global Jensen-Shannon-Divergence for the best candidate pair
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! The p-values from permutation test for the best candidate pair
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The JSD Confidence Intervals from bootstrapping for the best candidate pair. `-1.0_real64` if no candidate pair succeeded and fallback to `n_points=n_points_candidates(1)` and `n_neighbors=n_neighbors_candidates(1)`
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence intervals for candidate determination
            !!
            !! 1. JOIN_MIN: take min overlap of all studies' CI overlaps -> succeeds only if `all(ci_overlaps > min_neighbor_overlap)`
            !! 2. JOIN_MAX: take max overlap of all studies' CI overlaps -> succeeds only if `any(ci_overlaps > min_neighbor_overlap)`
            !! 3. JOIN_MEDIAN: take median overlap of all studies' CI overlaps -> succeeds only if `count(ci_overlaps > min_neighbor_overlap) >= (n_studies - 1) / 2 + 1`
        integer(int32), dimension(n_bins, max_n_points_candidate), intent(out) :: tmp_mean_pmf_counts
            !! Work array for proper resampling in permutation tests
        integer(int32), dimension(n_bins, max_n_points_candidate), intent(out) :: tmp_pmf_counts
            !! Work array for proper resampling in permutation tests
        real(real64), dimension(2, n_studies), intent(out), target :: tmp_confidence_interval
            !! Work array to hold the confidence intervals for a candidate pair. And also for permutation tests to hold the global jsd value per permutation.
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Work array for bootstrapping
        integer(int32), intent(in), optional :: min_count_per_mean_bin
            !! Number of minimum residuals a bin should have in the mean pmf to make a candidate pair eligible, default: `5`
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap in genes a neighborhood have to its succeeding neighborhood to make a candidate pair eligible, default: `0.1`
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap the confidence intervals should have to the current best confidence intervals (respecting the `join_method`) to make a candidate pair eligible, default: `0.9`
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations to perform in the permutation test ([[tox_data_integration(module):gjct_permutation_test(interface)]]), default: `1000`
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`

        integer(int32) :: n_gene_means, n_pool, i_study, i_neighbor_count, i_point_count, best_params_CI_i_point_count, best_params_CI_i_neighbor_count, best_params_exceeded_CI_overlap
        logical :: plateau_found, all_have_min_neighbor_overlap
        type(c_ptr) :: rng
        real(real64), dimension(:), pointer :: tmp_global_js_divergences
        integer(int32) :: actual_min_count_per_mean_bin, actual_n_permutations, actual_random_seed
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap

        M_DEFAULT_VAL(min_count_per_mean_bin, actual_min_count_per_mean_bin, 5_int32)
        M_DEFAULT_VAL(n_permutations, actual_n_permutations, 1000_int32)
        M_DEFAULT_VAL(random_seed, actual_random_seed, 42_int32)
        M_DEFAULT_VAL(min_neighbor_overlap, actual_min_neighbor_overlap, 0.1_real64)
        M_DEFAULT_VAL(succeeding_ci_overlap, actual_succeeding_ci_overlap, 0.9_real64)

        rng = create_rng(random_seed)

        n_gene_means = size(gene_means, kind=int32)
        n_pool = find_last_non_nan(gene_means, gene_means_perm_all, n_gene_means)

        best_candidate_pair_confidence_interval = -1.0_real64
        best_params_CI_i_point_count = 1
        best_params_CI_i_neighbor_count = 1
        best_params_exceeded_CI_overlap = 0
        plateau_found = .false.
        do i_point_count = 1, n_point_counts
            n_points = n_points_candidates(i_point_count)

            call pool_means_n_pool_input_helper(gene_means, gene_means_perm_all, n_gene_means, n_points, n_pool, x_star)

            do i_neighbor_count = 1, n_neighbor_counts
                n_neighbors = n_neighbors_candidates(i_neighbor_count)

                call reset_rng(rng, random_seed)

                all_have_min_neighbor_overlap = .true.
                do concurrent (i_study = 1:n_studies) shared(all_have_min_neighbor_overlap, n_points, x_star, max_n_genes_all_studies, gene_means, gene_means_perms, neighborhood_residuals, neighborhood_ranges, n_neighbors, shared_residual_range, n_bins, counts, pmfs, included_n_reps)
                    call construct_neighborhoods_helper(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), gene_means_perms(:, i_study), neighborhood_residuals(:, :, i_study), neighborhood_ranges(:, :, i_study), n_neighbors)
                    if (test_neighborhood_overlaps_helper(neighborhood_ranges(:, :, i_study), n_points, min_neighbor_overlap)) then
                        call build_residual_histograms_helper(neighborhood_residuals(:, :, i_study), n_neighbors, n_points, residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study))
                    else
                        all_have_min_neighbor_overlap = .false.
                    end if
                end do

                if (all_have_min_neighbor_overlap) then
                    call create_mean_pmf_helper(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts)
                    if (.not. test_mean_pmf_min_counts_helper(mean_pmf_counts, n_bins, n_points, min_count_per_mean_bin)) cycle

                    do concurrent (i_study = 1:n_studies) shared(pmfs, mean_pmf, n_points, n_bins, js_divergences, included_n_reps, mean_pmf_included_n_reps, weights)
                        call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
                        call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, tmp_confidence_interval(1, i_study), weights(:, i_study))
                        tmp_confidence_interval(2, i_study) = tmp_confidence_interval(1, i_study)
                    end do

                    call bootstrap_histogram_helper(n_bootstraps, included_n_reps, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf_included_n_reps, tmp_confidence_interval, js_divergences, weights, global_js_divergence, tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, tmp_pmf_counts, pmfs, mean_pmf, rng)
                    call check_plateau_condition(tmp_confidence_interval, best_candidate_pair_confidence_interval, n_studies, best_params_CI_i_point_count, best_params_CI_i_neighbor_count, best_params_exceeded_CI_overlap, i_point_count, i_neighbor_count, join_method, succeeding_ci_overlap, plateau_found)
                    if (plateau_found) exit
                end if
            end do
            if (plateau_found) exit
        end do
    
        n_points = n_points_candidates(best_params_CI_i_point_count)
        n_neighbors = n_neighbors_candidates(best_params_CI_i_neighbor_count)

        ! 1. Run the pipeline to determine the global JSDs
        call pool_means_helper(gene_means, gene_means_perm_all, n_gene_means, n_points, n_pool, x_star)

        do concurrent (i_study = 1:n_studies) shared(n_points, x_star, max_n_genes_all_studies, gene_means, gene_means_perms, neighborhood_residuals, neighborhood_ranges, n_neighbors)
            call construct_neighborhoods_helper(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), gene_means_perms(:, i_study), neighborhood_residuals(:, :, i_study), neighborhood_ranges(:, :, i_study), n_neighbors)
            call build_residual_histograms_helper(neighborhood_residuals(:, :, i_study), n_neighbors, n_points, residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study))
        end do
        call create_mean_pmf_helper(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts)

        do concurrent (i_study = 1:n_studies) shared(pmfs, mean_pmf, n_points, n_bins, js_divergences, included_n_reps, mean_pmf_included_n_reps, global_js_divergence, weights)
            call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, global_js_divergence(i_study), weights(:, i_study))
        end do

        ! 2. Run permutation tests
        call reset_rng(rng, random_seed)

        ! Reusing the workspace for temporary JSDs
        tmp_global_js_divergences(1:size(tmp_confidence_interval, kind=int32)) => tmp_confidence_interval
        call gjct_permutation_test_helper(n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, n_points, n_bins, included_n_reps, n_studies, global_js_divergence, p_values, pmfs, tmp_pmf_counts, tmp_mean_pmf_counts, js_divergences, weights, tmp_global_js_divergences(1:n_studies), rng)

        ! 3. Rerun the pipeline to recompute the values overwritten by permutation tests in `js_divergences`, `weights`, `pmfs`
        do concurrent (i_study = 1:n_studies) shared(counts, pmfs, included_n_reps, n_bins, n_points)
            call calc_pmf_helper(counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, global_js_divergence(i_study), weights(:, i_study))
        end do
    end subroutine js_comp_test

    pure subroutine check_plateau_condition(confidence_interval, best_candidate_pair_confidence_interval, n_studies, best_params_CI_i_point_count, best_params_CI_i_neighbor_count, best_params_exceeded_CI_overlap, i_point_count, i_neighbor_count, join_method, succeeding_ci_overlap, plateau_found)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(2, n_studies), intent(in) :: confidence_interval
            !! JSD confidence intervals from bootstrapping
        real(real64), dimension(2, n_studies), intent(inout) :: best_candidate_pair_confidence_interval
            !! JSD confidence intervals for the current best n_points/n_neighbors candidate pair, may be updated to `confidence_interval`
        integer(int32), intent(inout) :: best_params_CI_i_point_count
            !! Index for the `n_points_candidates` array for the current best n_points/n_neighbors candidate pair, may be updated to `i_point_count`
        integer(int32), intent(inout) :: best_params_CI_i_neighbor_count
            !! Index for the `n_neighbors_candidates` array for the current best n_points/n_neighbors candidate pair, may be updated to `i_point_count`
        integer(int32), intent(inout) :: best_params_exceeded_CI_overlap
            !! Number of studies for the current best n_points/n_neighbors candidate pair whose overlap exceeded `succeeding_ci_overlap`
        integer(int32), intent(in) :: i_point_count
            !! Index for the `n_points_candidates` array for the current n_points/n_neighbors candidate pair that caused `confidence_interval`
        integer(int32), intent(in) :: i_neighbor_count
            !! Index for the `n_neighbors_candidates` array for the current n_points/n_neighbors candidate pair that caused `confidence_interval`
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence intervals for candidate determination
            !!
            !! 1. JOIN_MIN: take min overlap of all studies' CI overlaps -> succeeds only if `all(ci_overlaps > min_neighbor_overlap)`
            !! 2. JOIN_MAX: take max overlap of all studies' CI overlaps -> succeeds only if `any(ci_overlaps > min_neighbor_overlap)`
            !! 3. JOIN_MEDIAN: take median overlap of all studies' CI overlaps -> succeeds only if `count(ci_overlaps > min_neighbor_overlap) >= (n_studies - 1) / 2 + 1`
        logical, intent(out) :: plateau_found
            !! `.true.` if plateau condition is met (as described for `join_method`), else `.false.`
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Minimum overlap an interval in `confidence_interval` should have with its respective one in `best_candidate_pair_confidence_interval`


        integer(int32) :: exceeds_min_CI_overlap, i_study

        exceeds_min_CI_overlap = 0_int32
        do concurrent (i_study = 1:n_studies) shared(confidence_interval, best_candidate_pair_confidence_interval) reduce(+:exceeds_min_CI_overlap)
            ! Check overlap
            ! IMPORTANT: best interval needs to be second argument, as the denominator will be the range of first interval
            ! Why this?
            !    If current interval is included in best interval, current is better -> 1.0
            !    If best interval is included in current interval, current is worse (larger range) -> <1.0
            if (compute_fractional_overlap_helper(&
                    confidence_interval(1, i_study), confidence_interval(2, i_study),&
                    best_candidate_pair_confidence_interval(1, i_study), best_candidate_pair_confidence_interval(2, i_study)&
                ) >= succeeding_ci_overlap)&
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
            best_candidate_pair_confidence_interval = confidence_interval

            select case (join_method)
                case (JOIN_MIN)
                    ! If all overlaps exceed the threshold, there is no lower one that doesn't
                    plateau_found = exceeds_min_CI_overlap == n_studies
                case (JOIN_MAX)
                    ! If at least one overlap exceeds the threshold, the max overlap will do either
                    plateau_found = exceeds_min_CI_overlap > 0
                case (JOIN_MEDIAN)
                    ! If 50% of the studies' overlaps exceed the threshold, the median will do either
                    plateau_found = exceeds_min_CI_overlap >= (n_studies - 1) / 2 + 1
            end select
        end if
    end subroutine check_plateau_condition

    subroutine gjct_permutation_test_helper(&
            n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, n_points, n_bins,&
            included_n_reps, n_studies, global_jsd_observed, p_values,&
            tmp_pmfs, tmp_counts, tmp_mean_pmf_counts, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, rng&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        integer(int32), intent(in) :: n_permutations
            !! Number of bootstraps to perform
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
        real(real64), dimension(n_bins, n_points), intent(in) :: mean_pmf
            !! The mean pmf built from all studies' pmfs
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for the studies' pmfs
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
            !! Computed JSD values for studies' pmfs
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-values of the permutation test, how many `x` permutations' JSDs exceeded the observation -> `p_values(i) = x(i) / n_permutations`
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: tmp_pmfs
            !! Work array that holds the permutations' pmfs
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_counts
            !! Work array that holds the permutations' pmfs' absolute counts
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_mean_pmf_counts
            !! Work array for proper resampling per permutation (keeps track of the residuals in the pool to pick from)
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Work array that holds the permutations' per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Work array that holds the permutations' per-point weight in the global JSD values
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Work array that holds the permutations' global weighted JSD values
        type(c_ptr), intent(in) :: rng
            !! random number generator from [[f42_random_gsl(module)]] to use

        integer(int32) :: i_permutation, i_point, i_study

        p_values = 0.0_real64
        do i_permutation = 1, n_permutations
            tmp_mean_pmf_counts = mean_pmf_counts
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multiv_hypergeom(rng, tmp_mean_pmf_counts(:, i_point), n_bins, sum(tmp_mean_pmf_counts(:, i_point)), included_n_reps(i_point, i_study), tmp_counts)
                end do
                call calc_pmf_helper(tmp_counts, tmp_pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            end do

            do concurrent (i_study = 1:n_studies) shared(p_values, tmp_pmfs, mean_pmf, n_points, n_bins, tmp_js_divergences, included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, tmp_weights, global_jsd_observed)
                call compute_divergence_per_reference_point_helper(tmp_pmfs(:, :, i_study), mean_pmf, n_points, n_bins, tmp_js_divergences(:, i_study))
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

    subroutine bootstrap_histogram_helper(n_bootstraps, included_n_reps, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf_included_n_reps, confidence_interval, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, tmp_counts, tmp_pmfs, tmp_mean_pmf, rng)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of elements in the top/bottom ends of bootstrapped values, e.g. for 95% of the values as `max(1, floor(0.025 * n_bootstraps))`,
            !! then the lower bound of the confidence interval would be the 2.5%-ile and the upper the 97.5%-ile
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: bootstrapping_top_k_jsds
            !! Used as heap for efficient percentile detection in the bootstrapped values
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! Confidence interval to be bootstrapped -> Initial values are reference
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for the studies' pmfs
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: tmp_pmfs
            !! Work array that holds the bootstraps' pmfs
        real(real64), dimension(n_bins, n_points), intent(out) :: tmp_mean_pmf
            !! Work array for the mean pmf for `tmp_pmfs`
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_counts
            !! Work array that holds the bootstraps' pmfs' absolute counts
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Work array that holds the bootstraps' per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Work array that holds the bootstraps' per-point weight in the global JSD values
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Work array that holds the bootstraps' global weighted JSD values
        type(c_ptr), intent(in) :: rng
            !! random number generator from [[f42_random_gsl(module)]] to use

        integer(int32) :: i_bootstrap, i_point, i_study

        ! Initialize heaps with initial confidence interval value
        do concurrent (i_study = 1:n_studies) shared(n_bootstrapping_top_k_jsds, bootstrapping_top_k_jsds, confidence_interval)
            bootstrapping_top_k_jsds(:, 1, i_study) = confidence_interval(1, i_study)
            bootstrapping_top_k_jsds(:, 2, i_study) = confidence_interval(2, i_study)
        end do

        do i_bootstrap = 1, n_bootstraps
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multinomial(rng, mean_pmf_counts(:, i_point), n_bins, mean_pmf_included_n_reps(i_point), included_n_reps(i_point, i_study), tmp_counts)
                end do
                call calc_pmf_helper(tmp_counts, tmp_pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            end do

            call create_mean_pmf_only_helper(tmp_pmfs, n_bins, n_points, n_studies, tmp_mean_pmf)

            do concurrent (i_study = 1:n_studies)
                call compute_divergence_per_reference_point_helper(tmp_pmfs(:, :, i_study), tmp_mean_pmf, n_points, n_bins, tmp_js_divergences(:, i_study))
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

    pure logical function test_mean_pmf_min_counts_helper(mean_pmf_counts, n_bins, n_points, min) result(all_bins_have_min_count)
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: min
            !! Minimum count that should be reached per bin
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`

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

    !> Helper to create the mean pmf and its histogram counts from all studies' pmfs
    pure subroutine create_mean_pmf_helper(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(n_bins, n_points, n_studies), intent(in) :: pmfs
            !! The probabilities of bins per point, from [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_bins, n_points, n_studies), intent(in) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64
        mean_pmf_included_n_reps = 0_int32
        mean_pmf_counts = 0_int32

        do i_study = 1, n_studies
            do i_point = 1, n_points
                do concurrent (i_bin = 1:n_bins) shared(mean_pmf, pmfs, n_studies, i_point, i_study, mean_pmf_counts, counts)
                    mean_pmf(i_bin, i_point) = mean_pmf(i_bin, i_point) + pmfs(i_bin, i_point, i_study) / real(n_studies, real64)
                    mean_pmf_counts(i_bin, i_point) = mean_pmf_counts(i_bin, i_point) + counts(i_bin, i_point, i_study)
                end do
                mean_pmf_included_n_reps(i_point) = mean_pmf_included_n_reps(i_point) + included_n_reps(i_point, i_study)
            end do
        end do
    
    end subroutine create_mean_pmf_helper

    !> Helper to create only the mean pmf from all studies' pmfs, not its histogram counts (helpful for [[tox_data_integration(module):bootstrap_histogram_helper(interface)]] where mean pmf counts don't matter)
    pure subroutine create_mean_pmf_only_helper(pmfs, n_bins, n_points, n_studies, mean_pmf)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(n_bins, n_points, n_studies), intent(in) :: pmfs
            !! The probabilities of bins per point, from [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64

        do i_study = 1, n_studies
            do i_point = 1, n_points
                do concurrent (i_bin = 1:n_bins) shared(mean_pmf, pmfs, n_studies, i_point, i_study)
                    mean_pmf(i_bin, i_point) = mean_pmf(i_bin, i_point) + pmfs(i_bin, i_point, i_study) / real(n_studies, real64)
                end do
            end do
        end do
    end subroutine create_mean_pmf_only_helper

    !> Helper function for [[tox_data_integration(module):js_comp_test(interface)]] to check if all neighborhoods for a certain n_points/n_neighbors pair have a consecutive min overlap
    pure logical function test_neighborhood_overlaps_helper(neighborhood_range, n_points, min_neighbor_overlap) result(all_have_min_neighbor_overlap)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), dimension(2, n_points), intent(in) :: neighborhood_range
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        real(real64), intent(in) :: min_neighbor_overlap
            !! Minimum overlap two consecutive neighborhoods (`i_point, i_point+1`) should have in their ranges

        integer(int32) :: i_point
        real(real64) :: overlap

        all_have_min_neighbor_overlap = .true.
        do concurrent (i_point = 1:n_points - 1) local(overlap) shared(neighborhood_range, all_have_min_neighbor_overlap)
            overlap = compute_fractional_overlap_helper(&
                real(neighborhood_range(1, i_point), real64),&
                real(neighborhood_range(2, i_point), real64),&
                real(neighborhood_range(1, i_point + 1), real64),&
                real(neighborhood_range(2, i_point + 1), real64)&
            )
            if (overlap < min_neighbor_overlap) all_have_min_neighbor_overlap = .false.
        end do
    end function test_neighborhood_overlaps_helper

    !> Calculates the fractional overlap between two intervals, the fraction of the overlap and the total range of both intervals, so in case `a_min < b_min`:
    !| \[ \frac{\min(\texttt{a_max}, \texttt{b_max}) - \texttt{b_min})} {\texttt{a_max} - \texttt{a_min}} \]
    pure real(real64) function compute_fractional_overlap_helper(a_min, a_max, b_min, b_max) result(overlap_percent)
        real(real64), intent(in) :: a_min
            !! Lower bound of the first interval
        real(real64), intent(in) :: a_max
            !! Upper bound of the first interval (assumed to be greater than `a_min`)
        real(real64), intent(in) :: b_min
            !! Lower bound of the second interval
        real(real64), intent(in) :: b_max
            !! Upper bound of the second interval (assumed to be greater than `b_min`)

        real(real64) :: left_max, right_min

        if (a_max == a_min) then
            ! If zero interval included in b interval, max overlap
            if (b_min <= a_max .and. b_max >= a_max) then
                overlap_percent = 1.0_real64
            else
                overlap_percent = 0.0_real64
            end if
        else
            left_max = min(a_max, b_max)
            right_min = max(a_min, b_min)
            ! Calculate overlap, is only negative if left_max < right_min -> a_min < a_max < b_min < b_max -> no overlap
            overlap_percent = max(0.0_real64, (left_max - right_min) / (a_max - a_min))
        end if
    end function compute_fractional_overlap_helper
end module tox_data_integration_js_comp_test