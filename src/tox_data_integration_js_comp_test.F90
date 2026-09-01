#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSCompTest
!|
!| This module implements the pipeline to perform the complete JSCompTest for all studies
module tox_data_integration_js_comp_test
    use safeguard
    use tox_data_integration_preprocessing, only: find_last_non_nan, pool_means_n_pool_input_helper, construct_neighborhoods_helper, pool_means_alloc
    use tox_data_integration_jsd, only: build_residual_histograms_helper, calc_pmf_helper, compute_weighted_global_divergence_helper, compute_divergence_per_reference_point_helper
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_heaps, only: top_k_heap_push, bottom_k_heap_push
    use f42_random_gsl, only: random_multinomial, create_rng, random_multiv_hypergeom, reset_rng, rng_t, destroy_rng
    use f42_utils, only: calc_percentile_helper, calc_percentile_rank, is_close, sort_array_heapsort, sort_real_heapsort_expl_size, LOG_2, init_perm, clamp
    use tox_errors, only: map_err_arg_pos, set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size, validate_in_range_real, validate_in_range_int, validate_all_in_range_int
    implicit none

#define CM_VALIDATE_JOIN_METHOD(ARG) call validate_in_range_int(join_method, ierr, min=0_int32, max=2_int32)
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
            !! Quantile for determining the residual range, default: `95.0`

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

            associate (&
                left_residual => abs(residuals(residuals_perm(left))),&
                right_residual => abs(residuals(residuals_perm(right)))&
            )
                if (left_residual > right_residual) then
                    shared_residual_range = left_residual
                    left = left + 1
                else
                    shared_residual_range = right_residual
                    right = right - 1
                end if
            end associate
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

    !> (with input validation) Root helper for:
    !|
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

        integer(int32), dimension(:), allocatable :: residuals_perm
        integer(int32) :: n_residuals

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

            call init_perm(residuals_perm)
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
    !| \[
    !|      \texttt{sturges_bins} = 1 + \lfloor \frac{\ln\left(\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}\right)}{\ln\left(2\right)} \rceil
    !| \]
    !| and the Freedman-Diaconis rule
    !| \[
    !|      \texttt{freed_diac_bins} = 2 \cdot \frac{\operatorname{IQR}(\texttt{residuals})}{\sqrt[3]{\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}}}
    !| \]
    !| so finally
    !| \[
    !|      \texttt{n_bins} = \max\left(\texttt{sturges_bins}, \texttt{freed_diac_bins}\right)
    !| \]
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
    !| \[
    !|      \texttt{sturges_bins} = 1 + \lfloor \frac{\ln\left(\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}\right)}{\ln\left(2\right)} \rceil
    !| \]
    !| and the Freedman-Diaconis rule
    !| \[
    !|      \texttt{freed_diac_bins} = 2 \cdot \frac{\operatorname{IQR}(\texttt{residuals})}{\sqrt[3]{\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}}}
    !| \]
    !| so finally
    !| \[
    !|      \texttt{n_bins} = \max\left(\texttt{sturges_bins}, \texttt{freed_diac_bins}\right)
    !| \]
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
    !| \[
    !|      \texttt{sturges_bins} = 1 + \lfloor \frac{\ln\left(\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}\right)}{\ln\left(2\right)} \rceil
    !| \]
    !| and the Freedman-Diaconis rule
    !| \[
    !|      \texttt{freed_diac_bins} = 2 \cdot \frac{\operatorname{IQR}(\texttt{residuals})}{\sqrt[3]{\texttt{max_n_reps_all_studies} \cdot \texttt{n_neighbors}}}
    !| \]
    !| so finally
    !| \[
    !|      \texttt{n_bins} = \max\left(\texttt{sturges_bins}, \texttt{freed_diac_bins}\right)
    !| \]
    pure subroutine estimate_bin_count_helper(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, n_neighbors, shared_residual_range, n_bins, &
            diag_quartile_25, diag_quartile_75, diag_n_pool, diag_n_reps_neighborhood, diag_half_bin_width, diag_n_bins_sturges, diag_n_bins_freedman_diaconis)
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
        real(real64), intent(out), optional :: diag_quartile_25
            !! Diagnostic: 25th percentile of the *global* (pooled) residual distribution. -1.0 if n_pool == 0.
        real(real64), intent(out), optional :: diag_quartile_75
            !! Diagnostic: 75th percentile of the *global* (pooled) residual distribution. -1.0 if n_pool == 0.
        integer(int32), intent(out), optional :: diag_n_pool
            !! Diagnostic: n_pool, the global count of non-NaN residuals used for the quartiles above
        real(real64), intent(out), optional :: diag_n_reps_neighborhood
            !! Diagnostic: max_n_reps_all_studies * n_neighbors, the "n" used as the Freedman-Diaconis sample size. -1.0 if n_pool == 0.
        real(real64), intent(out), optional :: diag_half_bin_width
            !! Diagnostic: Freedman-Diaconis half bin width estimate = (quartile_75 - quartile_25) / n_reps_neighborhood^(1/3). -1.0 if n_pool == 0.
        integer(int32), intent(out), optional :: diag_n_bins_sturges
            !! Diagnostic: the Sturges-only candidate for n_bins, before being combined with the Freedman-Diaconis estimate via max()
        integer(int32), intent(out), optional :: diag_n_bins_freedman_diaconis
            !! Diagnostic: the Freedman-Diaconis-only candidate for n_bins, before being combined with Sturges via max(). -1 if half_bin_width was ~0 (FD term skipped).

        integer(int32) :: n_pool
        real(real64) :: half_bin_width, quartile_25, quartile_75, n_reps_neighborhood

        n_pool = find_last_non_nan(residuals, residuals_perm, size(residuals, kind=int32))
        if (present(diag_n_pool)) diag_n_pool = n_pool

        if (n_pool == 0) then
            n_bins = 1
            if (present(diag_quartile_25)) diag_quartile_25 = -1.0_real64
            if (present(diag_quartile_75)) diag_quartile_75 = -1.0_real64
            if (present(diag_n_reps_neighborhood)) diag_n_reps_neighborhood = -1.0_real64
            if (present(diag_half_bin_width)) diag_half_bin_width = -1.0_real64
            if (present(diag_n_bins_sturges)) diag_n_bins_sturges = -1_int32
            if (present(diag_n_bins_freedman_diaconis)) diag_n_bins_freedman_diaconis = -1_int32
        else
            n_reps_neighborhood = real(max_n_reps_all_studies * n_neighbors, kind=real64)
            if (present(diag_n_reps_neighborhood)) diag_n_reps_neighborhood = n_reps_neighborhood

            ! Estimate n_bins
            ! Sturges
            n_bins = 1 + nint(log(n_reps_neighborhood) / LOG_2, kind=int32)
            if (present(diag_n_bins_sturges)) diag_n_bins_sturges = n_bins

            ! Freedman-Diaconis
            call calc_percentile_helper(residuals, residuals_perm(:n_pool), 25.0_real64, quartile_25)
            call calc_percentile_helper(residuals, residuals_perm(:n_pool), 75.0_real64, quartile_75)
            if (present(diag_quartile_25)) diag_quartile_25 = quartile_25
            if (present(diag_quartile_75)) diag_quartile_75 = quartile_75

            ! bin width as defined by Freedman-Diaconis, but without doubling the value.
            ! With doubling the value, the bin count would be calculated as `2*shared_residual_range/(2*bin_width)` -> the 2 is unnecessary
            half_bin_width = (quartile_75 - quartile_25) / (n_reps_neighborhood ** (1.0_real64 / 3.0_real64))
            if (present(diag_half_bin_width)) diag_half_bin_width = half_bin_width
            if (.not. is_close(half_bin_width, 0.0_real64)) then
                if (present(diag_n_bins_freedman_diaconis)) diag_n_bins_freedman_diaconis = nint(shared_residual_range / half_bin_width, kind=int32)
                n_bins = max(n_bins, nint(shared_residual_range / half_bin_width, kind=int32))
            else
                if (present(diag_n_bins_freedman_diaconis)) diag_n_bins_freedman_diaconis = -1_int32
            end if
        end if
    end subroutine estimate_bin_count_helper

    subroutine determine_js_comp_test_n_points_n_neighbors_alloc(&
            n_points, n_neighbors, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins,&
            gene_means, n_studies, n_bootstraps, best_candidate_pair_confidence_interval, join_method, ierr,&
            min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, two_sided_bootstrapping_significance_level, random_seed, residual_range_quantile&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
        integer(int32), intent(out) :: n_points
            !! The finally chosen candidate from `candidates_n_points`
        integer(int32), intent(out) :: n_neighbors
            !! The finally chosen candidate from `candidates_n_neighbors`
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: n_bins
            !! The number of bins used for the finally chosen candidate from `candidates_n_points_n_neighbors`
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The JSD Confidence Intervals from bootstrapping for the best candidate pair. `-1.0_real64` if no candidate pair succeeded and fallback to `n_points=candidates_n_points(1)` and `n_neighbors=candidates_n_neighbors(1)`
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence intervals for candidate determination
            !!
            !! 1. JOIN_MIN: take min overlap of all studies' CI overlaps -> succeeds only if `all(ci_overlaps > min_neighbor_overlap)`
            !! 2. JOIN_MAX: take max overlap of all studies' CI overlaps -> succeeds only if `any(ci_overlaps > min_neighbor_overlap)`
            !! 3. JOIN_MEDIAN: take median overlap of all studies' CI overlaps -> succeeds only if `count(ci_overlaps > min_neighbor_overlap) >= (n_studies - 1) / 2 + 1`
        integer(int32), intent(in), optional :: min_count_per_mean_bin
            !! Number of minimum residuals a bin should have in the mean pmf to make a candidate pair eligible, default: `5`
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap in genes a neighborhood have to its succeeding neighborhood to make a candidate pair eligible, default: `0.1`
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap the confidence intervals should have to the current best confidence intervals (respecting the `join_method`) to make a candidate pair eligible, default: `0.9`
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! The significance level used for obtained values in bootstrapping, default: `2.5` -> `best_candidate_pair_confidence_interval` caps `95%` from all obtained values
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: `95.0`
        logical, allocatable :: trace_neighbor_overlap_ok(:)
            !! .true. if ALL studies passed test_neighborhood_overlaps_helper for this candidate
        logical, allocatable :: trace_min_bin_count_ok(:)
            !! .true. if test_mean_pmf_min_counts_helper passed (consensus PMF had enough residuals per bin)
        integer(int32), allocatable :: trace_min_mean_pmf_bin_count(:)
            !! minval(mean_pmf_counts) actually observed for this candidate
        real(real64), allocatable :: trace_confidence_interval(:, :, :)
            !! bootstrapped [lower, upper] JSD CI per study, per candidate
        real(real64), allocatable :: trace_ci_overlap(:, :)
            !! per-study fractional CI overlap vs. the running best candidate
        integer(int32), allocatable :: trace_exceeds_ci_overlap(:)
            !! number of studies whose overlap exceeded succeeding_ci_overlap, per candidate
        logical, allocatable :: trace_plateau_found(:)
            !! .true. for the candidate at which the plateau condition triggered
        real(real64), allocatable :: trace_candidates_kx(:)
            !! The KX factor for each candidate pair, used to determine the plateau condition 
        integer(int32), allocatable :: trace_mean_pmf_counts_full(:, :, :)
            !! Full per-bin gene counts of the mean pmf, per candidate (bin, point, candidate)
        real(real64), allocatable :: trace_global_js_divergence(:, :)
            !! J_{i,t}: observed (point-estimate) global JSD per study, per candidate
        real(real64), allocatable :: trace_delta(:, :)
            !! Delta_{i,t}: relative change in J vs. the previous admissible candidate, per study
        real(real64), allocatable :: trace_delta_median(:)
            !! Delta-tilde_t: median of trace_delta(:, t) across studies
        real(real64), allocatable :: trace_delta_max(:)
            !! Delta_t^max: maximum of trace_delta(:, t) across studies
        real(real64), allocatable :: trace_x_star(:, :)
            !! Mean-expression reference point x_star(point_idx) for each candidate
        integer(int32), intent(out) :: ierr
            !! Error code

        ! needs to be ascending, as it is part of the denominator, so first elements are larger than next
        real(real64), dimension(2), parameter :: KX_FACTORS = [2.0_real64, 4.0_real64]

        integer(int32), parameter :: MAX_POINTS = 1500_int32
        integer(int32), parameter :: MIN_POINTS = 300_int32
        real(real64), parameter :: GAMMA = 0.8_real64
        integer(int32), parameter :: MAX_POINT_CANDIDATES = 8_int32 ! How often can I multiply `MAX_POINTS` with `GAMMA`
        integer(int32), parameter :: MAX_CANDIDATE_PAIRS = size(KX_FACTORS, kind=int32) * MAX_POINT_CANDIDATES

        integer(int32), dimension(:, :), allocatable :: candidates_n_points_n_neighbors
        integer(int32), dimension(:, :), allocatable :: gene_means_perms
        integer(int32), dimension(:), allocatable :: gene_means_perm_all, residuals_perm
        integer(int32), dimension(:, :), allocatable :: tmp_neighborhood_residuals
        integer(int32), dimension(:, :), allocatable :: tmp_neighborhood_ranges
        real(real64), dimension(:), allocatable :: tmp_x_star
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        integer(int32), dimension(:, :, :), allocatable :: tmp_counts
        integer(int32), dimension(:, :), allocatable :: tmp_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:), allocatable :: tmp_mean_pmf_included_n_reps, n_bins_candidates
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_confidence_interval
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds

        ! --- NEW: bin-count estimation diagnostics (see Fortran-Code Bins-Diskussion) ---
        real(real64), dimension(:), allocatable :: candidates_n_reps_neighborhood
        real(real64), dimension(:), allocatable :: candidates_half_bin_width
        integer(int32), dimension(:), allocatable :: candidates_n_bins_sturges
        integer(int32), dimension(:), allocatable :: candidates_n_bins_freedman_diaconis
        real(real64) :: diag_quartile_25, diag_quartile_75
        integer(int32) :: diag_n_pool

        integer(int32) :: n_bootstrapping_top_k_jsds, n_residuals, prev_point_candidate, n_candidates, max_n_bins_all_candidates
        integer(int32) :: i_point_candidate, point_candidate, prev_neighbor_candidate, i_neighbor_candidate, neighbor_candidate, max_n_points_candidate, max_n_neighbors_candidate
        integer(int32) :: i_study
        real(real64) :: n_points_high, n_points_low

        call set_ok(ierr)

        call validate_dimension_size(n_studies, ierr, arg_pos=9_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=5_int32)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=4_int32)

        call validate_in_range_int(n_bootstraps, ierr, arg_pos=10_int32, min=1_int32)
        call validate_in_range_int(min_count_per_mean_bin, ierr, arg_pos=14_int32, min=0_int32)
        CM_VALIDATE_JOIN_METHOD(arg_pos=12_int32)

        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=16_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=17_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=15_int32, min=0.0_real64, max=1.0_real64)

        if (is_err(ierr)) return

        if (present(two_sided_bootstrapping_significance_level)) then
            n_bootstrapping_top_k_jsds = max(1_int32, floor(two_sided_bootstrapping_significance_level / 100.0_real64 * real(n_bootstraps, real64), kind=int32))
        else
            n_bootstrapping_top_k_jsds = max(1_int32, floor(0.025_real64 * real(n_bootstraps, real64), kind=int32))
        end if


        ! 1. Determine shared residual range
        n_residuals = size(residuals, kind=int32)
        M_ALLOCATE(residuals_perm(n_residuals))
        call init_perm(residuals_perm)
        call sort_real_heapsort_expl_size(residuals, residuals_perm, n_residuals)
        call determine_shared_residual_range_helper(residuals, residuals_perm, n_residuals, shared_residual_range, residual_range_quantile)

        ! 2. Determine candidates
        M_ALLOCATE(candidates_n_points_n_neighbors(2, MAX_CANDIDATE_PAIRS))
        M_ALLOCATE(n_bins_candidates(MAX_CANDIDATE_PAIRS))
        M_ALLOCATE(trace_candidates_kx(MAX_CANDIDATE_PAIRS))
        M_ALLOCATE(candidates_n_reps_neighborhood(MAX_CANDIDATE_PAIRS))
        M_ALLOCATE(candidates_half_bin_width(MAX_CANDIDATE_PAIRS))
        M_ALLOCATE(candidates_n_bins_sturges(MAX_CANDIDATE_PAIRS))
        M_ALLOCATE(candidates_n_bins_freedman_diaconis(MAX_CANDIDATE_PAIRS))
        n_points_high = real(clamp(ceiling(4.0_real64 * sqrt(real(max_n_genes_all_studies, real64))), min_val=MIN_POINTS, max_val=MAX_POINTS), kind=real64)
        n_points_low = real(max(MIN_POINTS, ceiling(0.2_real64 * n_points_high)), kind=real64)

        prev_point_candidate = -1_int32
        n_candidates = 0_int32
        max_n_bins_all_candidates = 0_int32
        max_n_neighbors_candidate = 0_int32
        diag_quartile_25 = -1.0_real64
        diag_quartile_75 = -1.0_real64
        diag_n_pool = 0_int32

        do i_point_candidate = 1, MAX_POINT_CANDIDATES
            if (n_points_high < n_points_low) exit

            point_candidate = nint(n_points_high)
            if (point_candidate /= prev_point_candidate) then
                prev_point_candidate = point_candidate
                prev_neighbor_candidate = -1_int32

                do i_neighbor_candidate = 1, size(KX_FACTORS, kind=int32)
                    ! TO DISCUSS: 1 won't have any overlap
                    ! TO DISCUSS: if `max_n_genes_all_studies` is very large, the neighborhood genes are too far away from eachother
                    neighbor_candidate = max(1_int32, floor(real(max_n_genes_all_studies, real64) / (KX_FACTORS(i_neighbor_candidate) * n_points_high)))

                    if (neighbor_candidate /= prev_neighbor_candidate) then
                        prev_neighbor_candidate = neighbor_candidate

                        n_candidates = n_candidates + 1
                        candidates_n_points_n_neighbors(1, n_candidates) = point_candidate
                        candidates_n_points_n_neighbors(2, n_candidates) = neighbor_candidate

                        trace_candidates_kx(n_candidates) = KX_FACTORS(i_neighbor_candidate)

                        ! --- NEW: capture bin-count-estimation diagnostics per candidate.
                        ! quartile_25/quartile_75/n_pool are global (same every iteration) - just kept overwriting
                        ! them is harmless and avoids extra branching.
                        call estimate_bin_count_helper(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, neighbor_candidate, shared_residual_range, n_bins_candidates(n_candidates), &
                            diag_quartile_25=diag_quartile_25, diag_quartile_75=diag_quartile_75, diag_n_pool=diag_n_pool, &
                            diag_n_reps_neighborhood=candidates_n_reps_neighborhood(n_candidates), diag_half_bin_width=candidates_half_bin_width(n_candidates), &
                            diag_n_bins_sturges=candidates_n_bins_sturges(n_candidates), diag_n_bins_freedman_diaconis=candidates_n_bins_freedman_diaconis(n_candidates))

                        max_n_bins_all_candidates = max(max_n_bins_all_candidates, n_bins_candidates(n_candidates))
                        max_n_neighbors_candidate = max(max_n_neighbors_candidate, neighbor_candidate)
                    end if
                end do
            end if
            n_points_high = n_points_high * GAMMA
        end do
        max_n_points_candidate = candidates_n_points_n_neighbors(1, 1)

        ! 3. Allocate rest
        M_ALLOCATE(gene_means_perm_all(max_n_genes_all_studies * n_studies))
        call init_perm(gene_means_perm_all)
        call sort_real_heapsort_expl_size(gene_means, gene_means_perm_all, size(gene_means, kind=int32))
        M_ALLOCATE(gene_means_perms(max_n_genes_all_studies, n_studies))
        do concurrent (i_study = 1:n_studies) shared(gene_means, gene_means_perms)
            call init_perm(gene_means_perms(:, i_study))
            call sort_array_heapsort(gene_means(:, i_study), gene_means_perms(:, i_study))
        end do

        M_ALLOCATE(tmp_neighborhood_residuals(max_n_neighbors_candidate, max_n_points_candidate))
        M_ALLOCATE(tmp_neighborhood_ranges(2, max_n_points_candidate))
        M_ALLOCATE(tmp_x_star(max_n_points_candidate))
        M_ALLOCATE(tmp_pmfs(max_n_bins_all_candidates, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_counts(max_n_bins_all_candidates, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_included_n_reps(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_mean_pmf(max_n_bins_all_candidates, max_n_points_candidate))
        M_ALLOCATE(tmp_mean_pmf_counts(max_n_bins_all_candidates, max_n_points_candidate))
        M_ALLOCATE(tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        M_ALLOCATE(tmp_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))

        ! --- NEW: Allocate trace arrays if requested ---
        M_ALLOCATE(trace_neighbor_overlap_ok(n_candidates))
        M_ALLOCATE(trace_min_bin_count_ok(n_candidates))
        M_ALLOCATE(trace_min_mean_pmf_bin_count(n_candidates))
        M_ALLOCATE(trace_confidence_interval(2, n_studies, n_candidates))
        M_ALLOCATE(trace_ci_overlap(n_studies, n_candidates))
        M_ALLOCATE(trace_exceeds_ci_overlap(n_candidates))
        M_ALLOCATE(trace_plateau_found(n_candidates))
        M_ALLOCATE(trace_mean_pmf_counts_full(max_n_bins_all_candidates, max_n_points_candidate, n_candidates))
        M_ALLOCATE(trace_global_js_divergence(n_studies, n_candidates))
        M_ALLOCATE(trace_delta(n_studies, n_candidates))
        M_ALLOCATE(trace_delta_median(n_candidates))
        M_ALLOCATE(trace_delta_max(n_candidates))
        M_ALLOCATE(trace_x_star(max_n_points_candidate, n_candidates))

        call write_js_comp_test_bin_estimation_csv(&
            "js_comp_test_bin_estimation.csv", candidates_n_points_n_neighbors, trace_candidates_kx, n_candidates,&
            diag_quartile_25, diag_quartile_75, diag_n_pool,&
            candidates_n_reps_neighborhood, candidates_half_bin_width, candidates_n_bins_sturges, candidates_n_bins_freedman_diaconis, n_bins_candidates&
        )

        call determine_js_comp_test_n_points_n_neighbors_helper(&
            candidates_n_points_n_neighbors, n_candidates, max_n_points_candidate, max_n_neighbors_candidate,&
            n_points, n_neighbors, n_bins, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins_candidates, max_n_bins_all_candidates,&
            gene_means, gene_means_perms, gene_means_perm_all, n_studies,&
            n_bootstraps, n_bootstrapping_top_k_jsds, best_candidate_pair_confidence_interval, join_method,&
            tmp_neighborhood_residuals, tmp_neighborhood_ranges, tmp_x_star, tmp_pmfs, tmp_counts, tmp_included_n_reps,&
            tmp_mean_pmf, tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, tmp_global_js_divergence,&
            tmp_confidence_interval, tmp_bootstrapping_top_k_jsds,&
            min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, random_seed,&
            trace_neighbor_overlap_ok, trace_min_bin_count_ok, trace_min_mean_pmf_bin_count,&
            trace_confidence_interval, trace_ci_overlap, trace_exceeds_ci_overlap, trace_plateau_found, trace_candidates_kx,&
            trace_mean_pmf_counts_full, trace_global_js_divergence, trace_delta, trace_delta_median, trace_delta_max, trace_x_star&
        )

        !call write_js_comp_test_trace_csv( &
        !    "/home/laszl/TH-Bingen/Tensor-Omics/js_comp_test_trace.csv", &
        !    candidates_n_points_n_neighbors, &
        !    n_bins_candidates, &
        !    n_candidates, &
        !    n_studies, &
        !    trace_neighbor_overlap_ok, &
        !    trace_min_bin_count_ok, &
        !    trace_min_mean_pmf_bin_count, &
        !    trace_confidence_interval, &
        !    trace_ci_overlap)
    end subroutine determine_js_comp_test_n_points_n_neighbors_alloc

    subroutine determine_js_comp_test_n_points_n_neighbors(&
            candidates_n_points_n_neighbors, n_candidates_n_points_n_neighbors, max_n_points_candidate, max_n_neighbors_candidate,&
            n_points, n_neighbors, n_bins, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, candidates_n_bins, max_n_bins_all_candidates,&
            gene_means, gene_means_perms, gene_means_perm_all, n_studies,&
            n_bootstraps, n_bootstrapping_top_k_jsds, best_candidate_pair_confidence_interval, join_method,&
            tmp_neighborhood_residuals, tmp_neighborhood_ranges, tmp_x_star, tmp_pmfs, tmp_counts, tmp_included_n_reps,&
            tmp_mean_pmf, tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, tmp_global_js_divergence,&
            tmp_confidence_interval, tmp_bootstrapping_top_k_jsds, ierr,&
            min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, random_seed,&
            trace_neighbor_overlap_ok, trace_min_bin_count_ok, trace_min_mean_pmf_bin_count,&
            trace_confidence_interval, trace_ci_overlap, trace_exceeds_ci_overlap, trace_plateau_found, trace_candidates_kx,&
            trace_mean_pmf_counts_full, trace_global_js_divergence, trace_delta, trace_delta_median, trace_delta_max, trace_x_star&
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
        integer(int32), intent(in) :: n_candidates_n_points_n_neighbors
            !! Number of candidate pairs for `n_points` and `n_neighbors` in `candidates_n_points_n_neighbors`
        integer(int32), dimension(2, n_candidates_n_points_n_neighbors), intent(in) :: candidates_n_points_n_neighbors
            !! Candidates for `[n_points, n_neighbors]` pairs to test
        integer(int32), intent(in) :: max_n_points_candidate
            !! Value of the largest `n_points` candidate
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Value of the largest `n_neighbors` candidate
        integer(int32), intent(out) :: n_points
            !! The finally chosen candidate for `n_points` from `candidates_n_points_n_neighbors`
        integer(int32), intent(out) :: n_neighbors
            !! The finally chosen candidate for `n_neighbors` from `candidates_n_points_n_neighbors`
        integer(int32), intent(out) :: n_bins
            !! The number of bins used for the finally chosen candidate from `candidates_n_points_n_neighbors`
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(in), dimension(n_candidates_n_points_n_neighbors) :: candidates_n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for, one for each candidate pair.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study permutation vectors to sort the `gene_means`
        integer(int32), dimension(max_n_genes_all_studies * n_studies), intent(in) :: gene_means_perm_all
            !! Permutation vector to sort the flattened `gene_means`
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The JSD Confidence Intervals from bootstrapping for the best candidate pair. `-1.0_real64` if no candidate pair succeeded and fallback to `n_points=candidates_n_points(1)` and `n_neighbors=candidates_n_neighbors(1)`
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence intervals for candidate determination
            !!
            !! 1. JOIN_MIN: take min overlap of all studies' CI overlaps -> succeeds only if `all(ci_overlaps > min_neighbor_overlap)`
            !! 2. JOIN_MAX: take max overlap of all studies' CI overlaps -> succeeds only if `any(ci_overlaps > min_neighbor_overlap)`
            !! 3. JOIN_MEDIAN: take median overlap of all studies' CI overlaps -> succeeds only if `count(ci_overlaps > min_neighbor_overlap) >= (n_studies - 1) / 2 + 1`
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate), intent(out) :: tmp_neighborhood_residuals
            !! Indices of selected neighborhood genes per reference point from [[tox_data_integration(module):construct_neighborhoods(interface)]]
        integer(int32), dimension(2, max_n_points_candidate), intent(out) :: tmp_neighborhood_ranges
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Mean-expression reference points
        integer(int32), intent(in) :: max_n_bins_all_candidates
            !! Maximum `n_bins` value in `candidates_n_bins`
        real(real64), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_studies), intent(out) :: tmp_pmfs
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_studies), intent(out) :: tmp_counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
        real(real64), dimension(max_n_bins_all_candidates, max_n_points_candidate), intent(out) :: tmp_mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate), intent(out) :: tmp_mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! The per-reference-point Jensen-Shannon-Divergence -> finally `js_divergences(1:n_points, :)`
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! The per-reference-point weights for the Jensen-Shannon-Divergence -> finally `weights(1:n_points, :)`
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! The global Jensen-Shannon-Divergence for the best candidate pair
        real(real64), dimension(2, n_studies), intent(out) :: tmp_confidence_interval
            !! Work array to hold the confidence intervals for a candidate pair. And also for permutation tests to hold the global jsd value per permutation.
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Work array for bootstrapping
        integer(int32), intent(in), optional :: min_count_per_mean_bin
            !! Number of minimum residuals a bin should have in the mean pmf to make a candidate pair eligible, default: `5`
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap in genes a neighborhood have to its succeeding neighborhood to make a candidate pair eligible, default: `0.1`
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap the confidence intervals should have to the current best confidence intervals (respecting the `join_method`) to make a candidate pair eligible, default: `0.9`
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`
        logical, dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_neighbor_overlap_ok
            !! .true. if ALL studies passed test_neighborhood_overlaps_helper for this candidate
        logical, dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_min_bin_count_ok
            !! .true. if test_mean_pmf_min_counts_helper passed (consensus PMF had enough residuals per bin)
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_min_mean_pmf_bin_count
            !! minval(mean_pmf_counts) actually observed for this candidate
        real(real64), dimension(2, n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_confidence_interval
            !! bootstrapped [lower, upper] JSD CI per study, per candidate
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_ci_overlap
            !! per-study fractional CI overlap vs. the running best candidate
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_exceeds_ci_overlap
            !! number of studies whose overlap exceeded succeeding_ci_overlap, per candidate
        logical, dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_plateau_found
            !! .true. for the candidate at which the plateau condition triggered
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in), optional :: trace_candidates_kx
            !! The KX factor for each candidate pair, as determined by the caller during candidate generation (e.g. `_alloc`). This routine never computes k_x itself; it only forwards the value for tracing/CSV output.
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_mean_pmf_counts_full
            !! Full per-bin gene counts of the mean pmf, per candidate
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_global_js_divergence
            !! J_{i,t}: observed (point-estimate) global JSD per study, per candidate
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_delta
            !! Delta_{i,t}: relative change in J vs. the previous admissible candidate, per study
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_delta_median
            !! Delta-tilde_t: median of trace_delta(:, t) across studies
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_delta_max
            !! Delta_t^max: maximum of trace_delta(:, t) across studies
        real(real64), dimension(max_n_points_candidate, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_x_star
            !! Mean-expression reference point x_star(point_idx) for each candidate
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_candidate

        call set_ok(ierr)

        call validate_dimension_size(n_studies, ierr, arg_pos=17_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=10_int32)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=8_int32)
        call validate_dimension_size(n_bootstrapping_top_k_jsds, ierr, arg_pos=19_int32)
        call validate_dimension_size(n_candidates_n_points_n_neighbors, ierr, arg_pos=2_int32)
        call validate_dimension_size(max_n_points_candidate, ierr, arg_pos=3_int32)
        call validate_dimension_size(max_n_neighbors_candidate, ierr, arg_pos=4_int32)
        call validate_dimension_size(max_n_bins_all_candidates, ierr, arg_pos=13_int32)

        call validate_in_range_int(n_bootstraps, ierr, arg_pos=18_int32, min=1_int32)
        call validate_in_range_int(min_count_per_mean_bin, ierr, arg_pos=37_int32, min=0_int32)
        CM_VALIDATE_JOIN_METHOD(arg_pos=21_int32)

        call validate_in_range_real(shared_residual_range, ierr, arg_pos=11_int32, min=0.0_real64)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=38_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=39_int32, min=0.0_real64, max=1.0_real64)

        call validate_all_in_range_int(candidates_n_bins, n_candidates_n_points_n_neighbors, ierr, arg_pos=12_int32, min=1_int32, max=max_n_bins_all_candidates)
        call validate_all_in_range_int(gene_means_perms, size(gene_means_perms, kind=int32), ierr, arg_pos=15_int32, min=1_int32, max=size(gene_means_perms, kind=int32))
        call validate_all_in_range_int(gene_means_perm_all, size(gene_means_perm_all, kind=int32), ierr, arg_pos=16_int32, min=1_int32, max=size(gene_means_perm_all, kind=int32))

        do i_candidate = 1, n_candidates_n_points_n_neighbors
            call validate_in_range_int(candidates_n_points_n_neighbors(1, i_candidate), ierr, arg_pos=1_int32, min=1_int32, max=max_n_points_candidate)
            call validate_in_range_int(candidates_n_points_n_neighbors(2, i_candidate), ierr, arg_pos=1_int32, min=1_int32, max=max_n_neighbors_candidate)
        end do

        if (is_err(ierr)) return

        call determine_js_comp_test_n_points_n_neighbors_helper(&
            candidates_n_points_n_neighbors, n_candidates_n_points_n_neighbors, max_n_points_candidate, max_n_neighbors_candidate,&
            n_points, n_neighbors, n_bins, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, candidates_n_bins, max_n_bins_all_candidates,&
            gene_means, gene_means_perms, gene_means_perm_all, n_studies,&
            n_bootstraps, n_bootstrapping_top_k_jsds, best_candidate_pair_confidence_interval, join_method,&
            tmp_neighborhood_residuals, tmp_neighborhood_ranges, tmp_x_star, tmp_pmfs, tmp_counts, tmp_included_n_reps,&
            tmp_mean_pmf, tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, tmp_global_js_divergence,&
            tmp_confidence_interval, tmp_bootstrapping_top_k_jsds,&
            min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, random_seed,&
            trace_neighbor_overlap_ok, trace_min_bin_count_ok, trace_min_mean_pmf_bin_count,&
            trace_confidence_interval, trace_ci_overlap, trace_exceeds_ci_overlap, trace_plateau_found, trace_candidates_kx,&
            trace_mean_pmf_counts_full, trace_global_js_divergence, trace_delta, trace_delta_median, trace_delta_max, trace_x_star&
        )
    end subroutine determine_js_comp_test_n_points_n_neighbors

    subroutine determine_js_comp_test_n_points_n_neighbors_helper(&
            candidates_n_points_n_neighbors, n_candidates_n_points_n_neighbors, max_n_points_candidate, max_n_neighbors_candidate,&
            n_points, n_neighbors, n_bins, residuals, max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, candidates_n_bins, max_n_bins_all_candidates,&
            gene_means, gene_means_perms, gene_means_perm_all, n_studies,&
            n_bootstraps, n_bootstrapping_top_k_jsds, best_candidate_pair_confidence_interval, join_method,&
            tmp_neighborhood_residuals, tmp_neighborhood_ranges, tmp_x_star, tmp_pmfs, tmp_counts, tmp_included_n_reps,&
            tmp_mean_pmf, tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, tmp_global_js_divergence,&
            tmp_confidence_interval, tmp_bootstrapping_top_k_jsds,&
            min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, random_seed,&
            trace_neighbor_overlap_ok, trace_min_bin_count_ok, trace_min_mean_pmf_bin_count,&
            trace_confidence_interval, trace_ci_overlap, trace_exceeds_ci_overlap, trace_plateau_found, trace_candidates_kx,&
            trace_mean_pmf_counts_full, trace_global_js_divergence, trace_delta, trace_delta_median, trace_delta_max, trace_x_star&
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
        integer(int32), intent(in) :: n_candidates_n_points_n_neighbors
            !! Number of candidate pairs for `n_points` and `n_neighbors` in `candidates_n_points_n_neighbors`
        integer(int32), dimension(2, n_candidates_n_points_n_neighbors), intent(in) :: candidates_n_points_n_neighbors
            !! Candidates for `[n_points, n_neighbors]` pairs to test
        integer(int32), intent(in) :: max_n_points_candidate
            !! Value of the largest `n_points` candidate
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Value of the largest `n_neighbors` candidate
        integer(int32), intent(out) :: n_points
            !! The finally chosen candidate for `n_points` from `candidates_n_points_n_neighbors`
        integer(int32), intent(out) :: n_neighbors
            !! The finally chosen candidate for `n_neighbors` from `candidates_n_points_n_neighbors`
        integer(int32), intent(out) :: n_bins
            !! The number of bins used for the finally chosen candidate from `candidates_n_points_n_neighbors`
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: candidates_n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for, one for each candidate pair.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study permutation vectors to sort the `gene_means`
        integer(int32), dimension(max_n_genes_all_studies * n_studies), intent(in) :: gene_means_perm_all
            !! Permutation vector to sort the flattened `gene_means`
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The JSD Confidence Intervals from bootstrapping for the best candidate pair. `-1.0_real64` if no candidate pair succeeded and fallback to `n_points=candidates_n_points(1)` and `n_neighbors=candidates_n_neighbors(1)`
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence intervals for candidate determination
            !!
            !! 1. JOIN_MIN: take min overlap of all studies' CI overlaps -> succeeds only if `all(ci_overlaps > min_neighbor_overlap)`
            !! 2. JOIN_MAX: take max overlap of all studies' CI overlaps -> succeeds only if `any(ci_overlaps > min_neighbor_overlap)`
            !! 3. JOIN_MEDIAN: take median overlap of all studies' CI overlaps -> succeeds only if `count(ci_overlaps > min_neighbor_overlap) >= (n_studies - 1) / 2 + 1`
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate), intent(out) :: tmp_neighborhood_residuals
            !! Indices of selected neighborhood genes per reference point from [[tox_data_integration(module):construct_neighborhoods(interface)]]
        integer(int32), dimension(2, max_n_points_candidate), intent(out) :: tmp_neighborhood_ranges
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Mean-expression reference points
        integer(int32), intent(in) :: max_n_bins_all_candidates
            !! Maximum `n_bins` value in `candidates_n_bins`
        real(real64), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_studies), intent(out), target :: tmp_pmfs
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_studies), intent(out), target :: tmp_counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
        real(real64), dimension(max_n_bins_all_candidates, max_n_points_candidate), intent(out), target :: tmp_mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate), intent(out), target :: tmp_mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_js_divergences
            !! The per-reference-point Jensen-Shannon-Divergence -> finally `js_divergences(1:n_points, :)`
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_weights
            !! The per-reference-point weights for the Jensen-Shannon-Divergence -> finally `weights(1:n_points, :)`
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! The global Jensen-Shannon-Divergence for the best candidate pair
        real(real64), dimension(2, n_studies), intent(out) :: tmp_confidence_interval
            !! Work array to hold the confidence intervals for a candidate pair. And also for permutation tests to hold the global jsd value per permutation.
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Work array for bootstrapping
        integer(int32), intent(in), optional :: min_count_per_mean_bin
            !! Number of minimum residuals a bin should have in the mean pmf to make a candidate pair eligible, default: `5`
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap in genes a neighborhood have to its succeeding neighborhood to make a candidate pair eligible, default: `0.1`
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap the confidence intervals should have to the current best confidence intervals (respecting the `join_method`) to make a candidate pair eligible, default: `0.9`
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`
        logical, dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_neighbor_overlap_ok
            !! .true. if ALL studies passed test_neighborhood_overlaps_helper for this candidate
        logical, dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_min_bin_count_ok
            !! .true. if test_mean_pmf_min_counts_helper passed (consensus PMF had enough residuals per bin)
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_min_mean_pmf_bin_count
            !! minval(mean_pmf_counts) actually observed for this candidate -> directly comparable to min_count_per_mean_bin
        real(real64), dimension(2, n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_confidence_interval
            !! bootstrapped [lower, upper] JSD CI per study, per candidate
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_ci_overlap
            !! per-study fractional CI overlap vs. the running best candidate
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_exceeds_ci_overlap
            !! number of studies whose overlap exceeded succeeding_ci_overlap, per candidate
        logical, dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_plateau_found
            !! .true. for the candidate at which the plateau condition actually triggered (at most one .true. entry)
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in), optional :: trace_candidates_kx
            !! The KX factor for each candidate pair, as determined by the caller during candidate generation (e.g. `_alloc`). This routine never computes k_x itself; it only forwards the value for tracing/CSV output.
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_mean_pmf_counts_full
            !! Full per-bin gene counts of the mean pmf, per candidate. Only [1:n_bins, 1:n_points, i] are meaningful per candidate.
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_global_js_divergence
            !! J_{i,t}: observed (point-estimate) global JSD per study, per candidate. Captured before bootstrapping overwrites the work array. -1.0 if this candidate was never admissible far enough to compute it.
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_delta
            !! Delta_{i,t}: relative change in J vs. the previous *admissible* candidate, per study. -1.0 if there is no previous admissible candidate to compare against (e.g. the first admissible one).
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_delta_median
            !! Delta-tilde_t: median of trace_delta(:, t) across studies. -1.0 if undefined (no previous admissible candidate).
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(out), optional :: trace_delta_max
            !! Delta_t^max: maximum of trace_delta(:, t) across studies. -1.0 if undefined (no previous admissible candidate).
        real(real64), dimension(max_n_points_candidate, n_candidates_n_points_n_neighbors), intent(out), optional :: trace_x_star
            !! Mean-expression reference point x_star(point_idx) for each candidate. Only [1:n_points, i] are meaningful per candidate.
            !! Lets callers relate bin sparsity (js_comp_test_bin_counts.csv) back to the expression level of each neighborhood.
        logical, dimension(n_candidates_n_points_n_neighbors) :: local_trace_neighbor_overlap_ok
            !! Always-populated local mirror, used internally regardless of whether the caller requested trace output
        logical, dimension(n_candidates_n_points_n_neighbors) :: local_trace_min_bin_count_ok
        integer(int32), dimension(n_candidates_n_points_n_neighbors) :: local_trace_min_mean_pmf_bin_count
        real(real64), dimension(2, n_studies, n_candidates_n_points_n_neighbors) :: local_trace_confidence_interval
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors) :: local_trace_ci_overlap
        integer(int32), dimension(n_candidates_n_points_n_neighbors) :: local_trace_exceeds_ci_overlap
        logical, dimension(n_candidates_n_points_n_neighbors) :: local_trace_plateau_found
        real(real64), dimension(n_candidates_n_points_n_neighbors) :: local_trace_candidates_kx
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_candidates_n_points_n_neighbors) :: local_trace_mean_pmf_counts_full
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors) :: local_trace_global_js_divergence
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors) :: local_trace_delta
        real(real64), dimension(n_candidates_n_points_n_neighbors) :: local_trace_delta_median
        real(real64), dimension(n_candidates_n_points_n_neighbors) :: local_trace_delta_max
        real(real64), dimension(max_n_points_candidate, n_candidates_n_points_n_neighbors) :: local_trace_x_star
        real(real64), dimension(n_studies) :: prev_admissible_j
            !! J_{i,t-1} of the last admissible candidate, used to compute Delta_{i,t}
        logical :: has_prev_admissible
            !! .true. once at least one admissible candidate has been processed
        integer(int32), dimension(n_studies) :: delta_sort_perm
            !! Scratch permutation for the small median-of-studies computation
        real(real64), parameter :: DELTA_EPSILON = 1.0e-12_real64
            !! Small constant preventing division by zero in the relative JSD change Delta_{i,t}
        integer(int32) :: prev_n_points, n_gene_means, n_pool, i_study, i_candidate_n_points_n_neighbors, best_params_CI_i_candidate_n_points_n_neighbors, best_params_exceeded_CI_overlap
        logical :: plateau_found, all_have_min_neighbor_overlap
        integer(int32) :: actual_min_count_per_mean_bin
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64), dimension(:, :, :), pointer :: tmp_pmfs_view
        real(real64), dimension(:, :), pointer :: tmp_js_divergences_view, tmp_weights_view, tmp_mean_pmf_view
        integer(int32), dimension(:, :), pointer :: tmp_included_n_reps_view, tmp_mean_pmf_counts_view
        integer(int32), dimension(:, :, :), pointer :: tmp_counts_view

        M_DEFAULT_VAL(min_count_per_mean_bin, actual_min_count_per_mean_bin, 5_int32)
        M_DEFAULT_VAL(min_neighbor_overlap, actual_min_neighbor_overlap, 0.1_real64)
        M_DEFAULT_VAL(succeeding_ci_overlap, actual_succeeding_ci_overlap, 0.9_real64)

        call write_js_comp_test_run_metadata_csv(&
            "js_comp_test_run_metadata.csv", shared_residual_range, actual_min_neighbor_overlap,&
            actual_succeeding_ci_overlap, actual_min_count_per_mean_bin, n_bootstraps, join_method&
        )

        n_gene_means = size(gene_means, kind=int32)
        n_pool = find_last_non_nan(gene_means, gene_means_perm_all, n_gene_means)

        best_candidate_pair_confidence_interval = -1.0_real64
        best_params_CI_i_candidate_n_points_n_neighbors = 1
        best_params_CI_i_candidate_n_points_n_neighbors = 1
        best_params_exceeded_CI_overlap = 0

        plateau_found = .false.

        ! --- NEW: Initialize trace arrays (always-populated locals, regardless of caller's optional args) ---
        local_trace_neighbor_overlap_ok = .false.
        local_trace_min_bin_count_ok = .false.
        local_trace_min_mean_pmf_bin_count = 0
        local_trace_confidence_interval = -1.0_real64
        local_trace_ci_overlap = -1.0_real64
        local_trace_exceeds_ci_overlap = 0
        local_trace_plateau_found = .false.
        local_trace_mean_pmf_counts_full = 0
        local_trace_global_js_divergence = -1.0_real64
        local_trace_delta = -1.0_real64
        local_trace_delta_median = -1.0_real64
        local_trace_delta_max = -1.0_real64
        local_trace_x_star = -1.0_real64
        has_prev_admissible = .false.
        ! trace_candidates_kx is intent(in): this routine never computes k_x itself, only forwards
        ! whatever the caller (e.g. _alloc) supplied. -1.0 marks "caller did not provide k_x".
        if (present(trace_candidates_kx)) then
            local_trace_candidates_kx = trace_candidates_kx
        else
            local_trace_candidates_kx = -1.0_real64
        end if

        prev_n_points = -1

        ! Test candidate pairs and find JSD plateau
        do i_candidate_n_points_n_neighbors = 1, n_candidates_n_points_n_neighbors
            n_points = candidates_n_points_n_neighbors(1, i_candidate_n_points_n_neighbors)
            n_neighbors = candidates_n_points_n_neighbors(2, i_candidate_n_points_n_neighbors)
            n_bins = candidates_n_bins(i_candidate_n_points_n_neighbors)

            tmp_pmfs_view(1:n_bins, 1:n_points, 1:n_studies) => tmp_pmfs
            tmp_counts_view(1:n_bins, 1:n_points, 1:n_studies) => tmp_counts
            tmp_mean_pmf_view(1:n_bins, 1:n_points) => tmp_mean_pmf
            tmp_mean_pmf_counts_view(1:n_bins, 1:n_points) => tmp_mean_pmf_counts

            if (prev_n_points /= n_points) then
                prev_n_points = n_points

                tmp_js_divergences_view(1:n_points, 1:n_studies) => tmp_js_divergences
                tmp_weights_view(1:n_points, 1:n_studies) => tmp_weights
                tmp_included_n_reps_view(1:n_points, 1:n_studies) => tmp_included_n_reps

                ! determine x_star
                call pool_means_n_pool_input_helper(gene_means, gene_means_perm_all, n_gene_means, n_points, n_pool, tmp_x_star)
            end if

            ! --- NEW: capture x_star for this candidate (same values whenever n_points repeats, e.g. across k_x) ---
            local_trace_x_star(1:n_points, i_candidate_n_points_n_neighbors) = tmp_x_star(1:n_points)

            ! Calculate neighborhoods and check min overlap of consecutive neighbors
            all_have_min_neighbor_overlap = .true.
            do i_study = 1, n_studies
                call construct_neighborhoods_helper(n_points, tmp_x_star, max_n_genes_all_studies, gene_means(:, i_study), gene_means_perms(:, i_study), tmp_neighborhood_residuals, tmp_neighborhood_ranges, n_neighbors)
                if (test_neighborhood_overlaps_helper(tmp_neighborhood_ranges, n_points, actual_min_neighbor_overlap)) then
                    call build_residual_histograms_helper(tmp_neighborhood_residuals, n_neighbors, n_points, residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, tmp_counts_view(:, :, i_study), tmp_pmfs_view(:, :, i_study), tmp_included_n_reps_view(:, i_study))
                else
                    all_have_min_neighbor_overlap = .false.
                    exit
                end if
            end do

            ! --- NEW: trace neighbor overlap ---
            local_trace_neighbor_overlap_ok(i_candidate_n_points_n_neighbors) = all_have_min_neighbor_overlap

            ! If min consecutive overlap is met, check plateau condition 
            if (all_have_min_neighbor_overlap) then
                ! 1. Compute mean pmf and ensure all bins have a minimum bin count
                call create_mean_pmf_helper(tmp_pmfs_view, tmp_counts_view, n_bins, n_points, n_studies, tmp_included_n_reps_view, tmp_mean_pmf_view, tmp_mean_pmf_included_n_reps, tmp_mean_pmf_counts_view)
                
                ! --- NEW: capture bin sparsity BEFORE the early cycle ---
                local_trace_min_mean_pmf_bin_count(i_candidate_n_points_n_neighbors) = minval(tmp_mean_pmf_counts_view(1:n_bins, 1:n_points))
                local_trace_mean_pmf_counts_full(1:n_bins, 1:n_points, i_candidate_n_points_n_neighbors) = tmp_mean_pmf_counts_view(1:n_bins, 1:n_points)
                block
                    logical :: min_bin_count_ok
                    min_bin_count_ok = test_mean_pmf_min_counts_helper(tmp_mean_pmf_counts_view, n_bins, n_points, actual_min_count_per_mean_bin)
                    local_trace_min_bin_count_ok(i_candidate_n_points_n_neighbors) = min_bin_count_ok
                    if (.not. min_bin_count_ok) cycle
                end block

                ! 2. If min bin count is met, compute JSD observation
                do concurrent (i_study = 1:n_studies) shared(tmp_pmfs_view, tmp_mean_pmf_view, n_points, n_bins, tmp_js_divergences_view, tmp_included_n_reps_view, tmp_mean_pmf_included_n_reps, tmp_global_js_divergence, tmp_weights_view, tmp_confidence_interval)
                    call compute_divergence_per_reference_point_helper(tmp_pmfs_view(:, :, i_study), tmp_mean_pmf_view, n_points, n_bins, tmp_js_divergences_view(:, i_study))
                    call compute_weighted_global_divergence_helper(tmp_js_divergences_view(:, i_study), n_points, tmp_included_n_reps_view(:, i_study), tmp_mean_pmf_included_n_reps, tmp_global_js_divergence(i_study), tmp_weights_view(:, i_study))
                    tmp_confidence_interval(:, i_study) = tmp_global_js_divergence(i_study)
                end do

                ! --- NEW: capture J_{i,t} (point-estimate JSD) BEFORE bootstrapping overwrites tmp_global_js_divergence ---
                local_trace_global_js_divergence(:, i_candidate_n_points_n_neighbors) = tmp_global_js_divergence

                ! --- NEW: relative effect-size change Delta_{i,t} vs. the previous *admissible* candidate ---
                if (has_prev_admissible) then
                    do i_study = 1, n_studies
                        local_trace_delta(i_study, i_candidate_n_points_n_neighbors) = &
                            abs(tmp_global_js_divergence(i_study) - prev_admissible_j(i_study)) / max(prev_admissible_j(i_study), DELTA_EPSILON)
                    end do

                    ! Median and max across studies (n_studies is small, plain heapsort-based argsort is fine)
                    call init_perm(delta_sort_perm)
                    call sort_array_heapsort(local_trace_delta(:, i_candidate_n_points_n_neighbors), delta_sort_perm)
                    if (mod(n_studies, 2) == 0) then
                        local_trace_delta_median(i_candidate_n_points_n_neighbors) = ( &
                            local_trace_delta(delta_sort_perm(n_studies / 2), i_candidate_n_points_n_neighbors) + &
                            local_trace_delta(delta_sort_perm(n_studies / 2 + 1), i_candidate_n_points_n_neighbors) &
                        ) / 2.0_real64
                    else
                        local_trace_delta_median(i_candidate_n_points_n_neighbors) = &
                            local_trace_delta(delta_sort_perm(n_studies / 2 + 1), i_candidate_n_points_n_neighbors)
                    end if
                    local_trace_delta_max(i_candidate_n_points_n_neighbors) = maxval(local_trace_delta(:, i_candidate_n_points_n_neighbors))
                end if

                prev_admissible_j = tmp_global_js_divergence
                has_prev_admissible = .true.

                ! 3. Compute a confidence interval for the JSD value by bootstrapping
                ! Each candidate pair should have same conditions for comparability -> reset RNG
                call bootstrap_histogram_helper(n_bootstraps, tmp_included_n_reps_view, n_bins, n_points, n_studies, tmp_mean_pmf_counts_view, tmp_mean_pmf_included_n_reps, tmp_confidence_interval, tmp_js_divergences_view, tmp_weights_view, tmp_global_js_divergence, tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, tmp_counts_view, tmp_pmfs_view, tmp_mean_pmf_view, random_seed)
                
                ! --- NEW: capture CI + overlap BEFORE check_plateau_condition_helper ---
                local_trace_confidence_interval(:, :, i_candidate_n_points_n_neighbors) = tmp_confidence_interval
                
                do i_study = 1, n_studies
                    local_trace_ci_overlap(i_study, i_candidate_n_points_n_neighbors) = compute_fractional_overlap_helper(&
                        tmp_confidence_interval(1, i_study), tmp_confidence_interval(2, i_study),&
                        best_candidate_pair_confidence_interval(1, i_study), best_candidate_pair_confidence_interval(2, i_study)&
                    )
                end do
                
                call check_plateau_condition_helper(tmp_confidence_interval, best_candidate_pair_confidence_interval, n_studies, best_params_CI_i_candidate_n_points_n_neighbors, best_params_exceeded_CI_overlap, i_candidate_n_points_n_neighbors, join_method, actual_succeeding_ci_overlap, plateau_found)
                
                ! --- NEW: trace exceeds_ci_overlap and plateau_found ---
                local_trace_exceeds_ci_overlap(i_candidate_n_points_n_neighbors) = best_params_exceeded_CI_overlap
                local_trace_plateau_found(i_candidate_n_points_n_neighbors) = plateau_found
                
                if (plateau_found) exit
            end if
        end do

        ! ------------------------------------------------------------
        ! Write trace information for debugging/visualization
        ! ------------------------------------------------------------

        write(*,*) "DEBUG: About to write JS comp test trace CSV"



        call write_js_comp_test_trace_csv( &
            "js_comp_test_trace.csv", &
            candidates_n_points_n_neighbors, &
            candidates_n_bins, &
            n_candidates_n_points_n_neighbors, &
            n_studies, &
            local_trace_neighbor_overlap_ok, &
            local_trace_min_bin_count_ok, &
            local_trace_min_mean_pmf_bin_count, &
            local_trace_confidence_interval, &
            local_trace_ci_overlap, &
            local_trace_exceeds_ci_overlap, &
            local_trace_plateau_found, &
            local_trace_candidates_kx &
        )

        call write_js_comp_test_bin_counts_csv(&
            "js_comp_test_bin_counts.csv", candidates_n_points_n_neighbors, candidates_n_bins,&
            n_candidates_n_points_n_neighbors, max_n_bins_all_candidates, max_n_points_candidate,&
            local_trace_mean_pmf_counts_full, local_trace_neighbor_overlap_ok, local_trace_min_bin_count_ok, local_trace_x_star&
        )

        call write_js_comp_test_effect_size_csv(&
            "js_comp_test_effect_size.csv", candidates_n_points_n_neighbors, n_candidates_n_points_n_neighbors, n_studies,&
            local_trace_global_js_divergence, local_trace_confidence_interval, local_trace_delta,&
            local_trace_delta_median, local_trace_delta_max&
        )

        ! Nur zurückschreiben, wenn der Aufrufer die jeweilige Trace-Variable tatsächlich angefordert hat
        if (present(trace_neighbor_overlap_ok)) trace_neighbor_overlap_ok = local_trace_neighbor_overlap_ok
        if (present(trace_min_bin_count_ok)) trace_min_bin_count_ok = local_trace_min_bin_count_ok
        if (present(trace_min_mean_pmf_bin_count)) trace_min_mean_pmf_bin_count = local_trace_min_mean_pmf_bin_count
        if (present(trace_confidence_interval)) trace_confidence_interval = local_trace_confidence_interval
        if (present(trace_ci_overlap)) trace_ci_overlap = local_trace_ci_overlap
        if (present(trace_exceeds_ci_overlap)) trace_exceeds_ci_overlap = local_trace_exceeds_ci_overlap
        if (present(trace_plateau_found)) trace_plateau_found = local_trace_plateau_found
        ! trace_candidates_kx is intent(in) now - nothing to copy back, it's caller-supplied input
        if (present(trace_mean_pmf_counts_full)) trace_mean_pmf_counts_full = local_trace_mean_pmf_counts_full
        if (present(trace_global_js_divergence)) trace_global_js_divergence = local_trace_global_js_divergence
        if (present(trace_delta)) trace_delta = local_trace_delta
        if (present(trace_delta_median)) trace_delta_median = local_trace_delta_median
        if (present(trace_delta_max)) trace_delta_max = local_trace_delta_max
        if (present(trace_x_star)) trace_x_star = local_trace_x_star

        ! assign final candidate pair, will be first pair if no candidate pair met the plateau condition
        ! If there is only one candidate pair, this one is the plateau
        if (plateau_found .or. n_candidates_n_points_n_neighbors < 2) then
            n_points = candidates_n_points_n_neighbors(1, best_params_CI_i_candidate_n_points_n_neighbors)
            n_neighbors = candidates_n_points_n_neighbors(2, best_params_CI_i_candidate_n_points_n_neighbors)
            n_bins = candidates_n_bins(best_params_CI_i_candidate_n_points_n_neighbors)
        else
            n_points = candidates_n_points_n_neighbors(1, 1)
            n_neighbors = candidates_n_points_n_neighbors(2, 1)
            n_bins = candidates_n_bins(1)
            best_candidate_pair_confidence_interval = -1.0_real64
        end if
    end subroutine determine_js_comp_test_n_points_n_neighbors_helper

    !> Performs the pipeline:
    !| 
    !| 1. [[tox_data_integration(module):construct_neighborhoods(interface)]]
    !| 2. [[tox_data_integration(module):build_residual_histograms(interface)]]
    !| 3. [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    !| 4. [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    !|
    subroutine js_comp_test_alloc(&
            gene_means, max_n_genes_all_studies, n_studies, residuals, shared_residual_range,&
            n_bins, max_n_reps_all_studies, x_star, n_pool, n_points, n_neighbors, neighborhood_ranges, neighborhood_residuals,&
            pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
            js_divergences, weights, global_js_divergence, p_values, ierr, n_permutations, random_seed&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in neighborhoods
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        integer(int32), dimension(2, n_points, n_studies), intent(out) :: neighborhood_ranges
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        real(real64), dimension(n_points), intent(out) :: x_star
            !! Mean-expression reference points
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        integer(int32), dimension(n_neighbors, n_points, n_studies), intent(out) :: neighborhood_residuals
            !! Indices of selected neighborhood genes per reference point from [[tox_data_integration(module):construct_neighborhoods(interface)]]
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: pmfs
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_bins, n_points, n_studies), intent(out) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! The per-reference-point Jensen-Shannon-Divergence -> finally `js_divergences(1:n_points, :)`
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! The per-reference-point weights for the Jensen-Shannon-Divergence -> finally `weights(1:n_points, :)`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! The global Jensen-Shannon-Divergence
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! The p-values from permutation test
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations to perform in the permutation test ([[tox_data_integration(module):gjct_permutation_test(interface)]]), default: `1000`
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        integer(int32), dimension(:, :), allocatable :: gene_means_perms, tmp_mean_pmf_counts, tmp_pmf_counts
        integer(int32) :: i_study

        call set_ok(ierr)

        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=2_int32)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=7_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=10_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=11_int32)

        call validate_in_range_real(shared_residual_range, ierr, arg_pos=5_int32, min=0.0_real64)

        call validate_in_range_int(n_bins, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=25_int32, min=0_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_mean_pmf_counts(n_bins, n_points))
        M_ALLOCATE(tmp_pmf_counts(n_bins, n_points))
        M_ALLOCATE(gene_means_perms(max_n_genes_all_studies, n_studies))

        call pool_means_alloc(gene_means, n_studies, max_n_genes_all_studies, n_points, n_pool, x_star, ierr)

        if (is_err(ierr)) return

        do concurrent (i_study = 1:n_studies) shared(gene_means, gene_means_perms)
            call init_perm(gene_means_perms(:, i_study))
            call sort_array_heapsort(gene_means(:, i_study), gene_means_perms(:, i_study))
        end do

        call js_comp_test_helper(&
            gene_means, max_n_genes_all_studies, n_studies, gene_means_perms, residuals, shared_residual_range,&
            n_bins, max_n_reps_all_studies, x_star, n_points, n_neighbors, neighborhood_ranges, neighborhood_residuals,&
            pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
            js_divergences, weights, global_js_divergence, p_values,&
            tmp_global_js_divergence, tmp_pmf_counts, tmp_mean_pmf_counts, n_permutations, random_seed&
        )
    end subroutine js_comp_test_alloc

    !> Performs the pipeline:
    !| 
    !| 1. [[tox_data_integration(module):construct_neighborhoods(interface)]]
    !| 2. [[tox_data_integration(module):build_residual_histograms(interface)]]
    !| 3. [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    !| 4. [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    !|
    subroutine js_comp_test(&
            gene_means, max_n_genes_all_studies, n_studies, gene_means_perms, residuals, shared_residual_range,&
            n_bins, max_n_reps_all_studies, x_star, n_points, n_neighbors, neighborhood_ranges, neighborhood_residuals,&
            pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
            js_divergences, weights, global_js_divergence, p_values,&
            tmp_global_js_divergence, tmp_pmf_counts, tmp_mean_pmf_counts, ierr, n_permutations, random_seed&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in neighborhoods
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study permutation vectors to sort the `gene_means`
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        integer(int32), dimension(2, n_points, n_studies), intent(out) :: neighborhood_ranges
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
        integer(int32), dimension(n_neighbors, n_points, n_studies), intent(out) :: neighborhood_residuals
            !! Indices of selected neighborhood genes per reference point from [[tox_data_integration(module):construct_neighborhoods(interface)]]
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: pmfs
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_bins, n_points, n_studies), intent(out) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! The per-reference-point Jensen-Shannon-Divergence -> finally `js_divergences(1:n_points, :)`
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! The per-reference-point weights for the Jensen-Shannon-Divergence -> finally `weights(1:n_points, :)`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! The global Jensen-Shannon-Divergence
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! The p-values from permutation test
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_mean_pmf_counts
            !! Work array for proper resampling in permutation tests
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_pmf_counts
            !! Work array for proper resampling in permutation tests
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Work array for the permutations' global Jensen-Shannon-Divergences
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations to perform in the permutation test ([[tox_data_integration(module):gjct_permutation_test(interface)]]), default: `1000`
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=2_int32)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=8_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=10_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=11_int32)

        call validate_in_range_real(shared_residual_range, ierr, arg_pos=6_int32, min=0.0_real64)

        call validate_in_range_int(n_bins, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=28_int32, min=0_int32)

        call validate_all_in_range_int(gene_means_perms, size(gene_means_perms, kind=int32), ierr, arg_pos=4_int32, min=1_int32, max=size(gene_means_perms, kind=int32))

        if (is_err(ierr)) return

        call js_comp_test_helper(&
            gene_means, max_n_genes_all_studies, n_studies, gene_means_perms, residuals, shared_residual_range,&
            n_bins, max_n_reps_all_studies, x_star, n_points, n_neighbors, neighborhood_ranges, neighborhood_residuals,&
            pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
            js_divergences, weights, global_js_divergence, p_values,&
            tmp_global_js_divergence, tmp_pmf_counts, tmp_mean_pmf_counts, n_permutations, random_seed&
        )
    end subroutine js_comp_test

    !> (no input validation) Performs the pipeline:
    !| 
    !| 1. [[tox_data_integration(module):construct_neighborhoods(interface)]]
    !| 2. [[tox_data_integration(module):build_residual_histograms(interface)]]
    !| 3. [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    !| 4. [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
    !|
    subroutine js_comp_test_helper(&
            gene_means, max_n_genes_all_studies, n_studies, gene_means_perms, residuals, shared_residual_range,&
            n_bins, max_n_reps_all_studies, x_star, n_points, n_neighbors, neighborhood_ranges, neighborhood_residuals,&
            pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
            js_divergences, weights, global_js_divergence, p_values,&
            tmp_global_js_divergence, tmp_pmf_counts, tmp_mean_pmf_counts, n_permutations, random_seed&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in neighborhoods
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study permutation vectors to sort the `gene_means`
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
        integer(int32), dimension(2, n_points, n_studies), intent(out) :: neighborhood_ranges
            !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
            !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
            !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
            !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
            !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
        integer(int32), dimension(n_neighbors, n_points, n_studies), intent(out) :: neighborhood_residuals
            !! Indices of selected neighborhood genes per reference point from [[tox_data_integration(module):construct_neighborhoods(interface)]]
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: pmfs
            !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
        integer(int32), dimension(n_bins, n_points, n_studies), intent(out) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! The per-reference-point Jensen-Shannon-Divergence -> finally `js_divergences(1:n_points, :)`
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! The per-reference-point weights for the Jensen-Shannon-Divergence -> finally `weights(1:n_points, :)`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! The global Jensen-Shannon-Divergence
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! The p-values from permutation test
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_mean_pmf_counts
            !! Work array for proper resampling in permutation tests
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_pmf_counts
            !! Work array for proper resampling in permutation tests
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Work array for the permutations' global Jensen-Shannon-Divergences
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations to perform in the permutation test ([[tox_data_integration(module):gjct_permutation_test(interface)]]), default: `1000`
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`

        integer(int32) :: i_study

        do concurrent (i_study = 1:n_studies) shared(n_points, x_star, max_n_genes_all_studies, gene_means, gene_means_perms, neighborhood_residuals, neighborhood_ranges, n_neighbors, residuals, max_n_reps_all_studies, shared_residual_range, n_bins, counts, pmfs, included_n_reps)
            call construct_neighborhoods_helper(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), gene_means_perms(:, i_study), neighborhood_residuals(:, :, i_study), neighborhood_ranges(:, :, i_study), n_neighbors)
            call build_residual_histograms_helper(neighborhood_residuals(:, :, i_study), n_neighbors, n_points, residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, shared_residual_range, n_bins, counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study))
        end do
        call create_mean_pmf_helper(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, mean_pmf_included_n_reps, mean_pmf_counts)

        do concurrent (i_study = 1:n_studies) shared(pmfs, mean_pmf, n_points, n_bins, js_divergences, included_n_reps, mean_pmf_included_n_reps, global_js_divergence, weights)
            call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, global_js_divergence(i_study), weights(:, i_study))
        end do

        ! 2. Run permutation tests
        call gjct_permutation_test_helper(mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, n_points, n_bins, included_n_reps, n_studies, global_js_divergence, p_values, pmfs, tmp_pmf_counts, tmp_mean_pmf_counts, js_divergences, weights, tmp_global_js_divergence, n_permutations, random_seed)

        ! 3. Rerun the pipeline to recompute the values overwritten by permutation tests in `js_divergences`, `weights`, `pmfs`
        do concurrent (i_study = 1:n_studies) shared(counts, pmfs, included_n_reps, n_bins, n_points, mean_pmf, js_divergences, mean_pmf_included_n_reps, global_js_divergence, weights)
            call calc_pmf_helper(counts(:, :, i_study), pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            call compute_divergence_per_reference_point_helper(pmfs(:, :, i_study), mean_pmf, n_points, n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_helper(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, global_js_divergence(i_study), weights(:, i_study))
        end do
    end subroutine js_comp_test_helper

    !> Helper for [[tox_data_integration(module):determine_js_comp_test_n_points_n_neighbors(interface)]] to test a candidate on the plateau condition
    subroutine check_plateau_condition_helper(confidence_interval, best_candidate_pair_confidence_interval, n_studies, best_params_CI_i_candidate_n_points_n_neighbors, best_params_exceeded_CI_overlap, i_candidate_n_points_n_neighbors, join_method, succeeding_ci_overlap, plateau_found)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(2, n_studies), intent(in) :: confidence_interval
            !! JSD confidence intervals from bootstrapping
        real(real64), dimension(2, n_studies), intent(inout) :: best_candidate_pair_confidence_interval
            !! JSD confidence intervals for the current best n_points/n_neighbors candidate pair, may be updated to `confidence_interval`
        integer(int32), intent(inout) :: best_params_CI_i_candidate_n_points_n_neighbors
            !! Index for the `candidates_n_points_n_neighbors` matrix for the current best n_points/n_neighbors candidate pair, may be updated to `i_candidate_n_points_n_neighbors`
        integer(int32), intent(inout) :: best_params_exceeded_CI_overlap
            !! Number of studies for the current best n_points/n_neighbors candidate pair whose overlap exceeded `succeeding_ci_overlap`
        integer(int32), intent(in) :: i_candidate_n_points_n_neighbors
            !! Index for the `candidates_n_points_n_neighbors` matrix for the current n_points/n_neighbors candidate pair that caused `confidence_interval`
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
        do concurrent (i_study = 1:n_studies) shared(confidence_interval, best_candidate_pair_confidence_interval, succeeding_ci_overlap) reduce(+:exceeds_min_CI_overlap)
            associate (&
                ci_min => confidence_interval(1, i_study),&
                ci_max => confidence_interval(2, i_study),&
                best_ci_min => best_candidate_pair_confidence_interval(1, i_study),&
                best_ci_max => best_candidate_pair_confidence_interval(2, i_study)&
            )
                ! Check overlap
                ! IMPORTANT: best interval needs to be second argument, as the denominator will be the range of first interval
                ! Why this?
                !    If current interval is included in best interval, current is better -> 1.0
                !    If best interval is included in current interval, current is worse (larger range) -> <1.0
                if (compute_fractional_overlap_helper(ci_min, ci_max, best_ci_min, best_ci_max) >= succeeding_ci_overlap) then
                    exceeds_min_CI_overlap = exceeds_min_CI_overlap + 1
                end if
            end associate
        end do

        ! If the parameter set is not at least as good as the previous, exit and use previous pair as best
        if (exceeds_min_CI_overlap < best_params_exceeded_CI_overlap) then
            plateau_found = .true.
        else
            best_params_CI_i_candidate_n_points_n_neighbors = i_candidate_n_points_n_neighbors
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
    end subroutine check_plateau_condition_helper

    !> (no input validation) Resample histograms without replacement and recompute the JSD to estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable.
    subroutine gjct_permutation_test_helper(&
            mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, n_points, n_bins,&
            included_n_reps, n_studies, global_jsd_observed, p_values,&
            tmp_pmfs, tmp_counts, tmp_mean_pmf_counts, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, n_permutations, random_seed&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Appropriate number of bins to do the JSD Compatibility test for
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
            !! Work array for proper resampling per permutation (keeps track of the residuals in the pool to draw from)
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Work array that holds the permutations' per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Work array that holds the permutations' per-point weight in the global JSD values
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Work array that holds the permutations' global weighted JSD values
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations to perform, default: 1000
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`

        integer(int32) :: i_permutation, i_point, i_study, actual_n_permutations
        type(rng_t) :: rng

        M_DEFAULT_VAL(n_permutations, actual_n_permutations, 1000_int32)
        rng = create_rng(random_seed)

        p_values = 0.0_real64
        do i_permutation = 1, actual_n_permutations
            ! Resample histogram, each reference point becomes a pool without replacement -> draw from all studies' residuals
            ! tmp_mean_pmf_counts will be updated in-place. It represents the pool to draw from
            tmp_mean_pmf_counts = mean_pmf_counts
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multiv_hypergeom(rng, tmp_mean_pmf_counts(:, i_point), n_bins, sum(tmp_mean_pmf_counts(:, i_point)), included_n_reps(i_point, i_study), tmp_counts(:, i_point))
                end do
                call calc_pmf_helper(tmp_counts, tmp_pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            end do

            do concurrent (i_study = 1:n_studies) shared(tmp_pmfs, mean_pmf, n_points, n_bins, tmp_js_divergences, included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, tmp_weights, global_jsd_observed, p_values)
                call compute_divergence_per_reference_point_helper(tmp_pmfs(:, :, i_study), mean_pmf, n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_helper(tmp_js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                if (tmp_global_js_divergence(i_study) >= global_jsd_observed(i_study)) then
                    p_values(i_study) = p_values(i_study) + 1.0_real64
                end if
            end do
        end do

        if (n_permutations == 0) return
        do concurrent (i_study = 1:n_studies) shared(p_values, n_permutations)
            p_values(i_study) = anint(p_values(i_study)) / real(n_permutations, real64)
        end do

        call destroy_rng(rng)
    end subroutine gjct_permutation_test_helper

    !> (no input validation) Resample histograms with replacement and recompute the JSD to get a confidence interval for it.
    subroutine bootstrap_histogram_helper(n_bootstraps, included_n_reps, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf_included_n_reps, confidence_interval, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, tmp_counts, tmp_pmfs, tmp_mean_pmf, random_seed)
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
        integer(int32), intent(in), optional :: random_seed
            !! Random seed to use for random number generation, default: `42`

        integer(int32) :: i_bootstrap, i_point, i_study
        type(rng_t) :: rng

        rng = create_rng(random_seed)

        ! Initialize heaps with initial confidence interval value
        do concurrent (i_study = 1:n_studies) shared(n_bootstrapping_top_k_jsds, bootstrapping_top_k_jsds, confidence_interval)
            bootstrapping_top_k_jsds(:, 1, i_study) = confidence_interval(1, i_study)
            bootstrapping_top_k_jsds(:, 2, i_study) = confidence_interval(2, i_study)
        end do

        do i_bootstrap = 1, n_bootstraps
            ! Resample histogram, each reference point becomes a pool with replacement -> draw from all studies' residuals
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multinomial(rng, mean_pmf_counts(:, i_point), n_bins, mean_pmf_included_n_reps(i_point), included_n_reps(i_point, i_study), tmp_counts(:, i_point))
                end do
                call calc_pmf_helper(tmp_counts, tmp_pmfs(:, :, i_study), included_n_reps(:, i_study), n_bins, n_points)
            end do

            ! As resampling was with replacement, the mean pmf changed very likely
            call create_mean_pmf_only_helper(tmp_pmfs, n_bins, n_points, n_studies, tmp_mean_pmf)

            ! For new pmfs from resampling, run the pipeline to determine the JSD
            do concurrent (i_study = 1:n_studies) shared(tmp_pmfs, tmp_mean_pmf, n_points, n_bins, tmp_js_divergences, included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, tmp_weights, bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds)
                call compute_divergence_per_reference_point_helper(tmp_pmfs(:, :, i_study), tmp_mean_pmf, n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_helper(tmp_js_divergences(:, i_study), n_points, included_n_reps(:, i_study), mean_pmf_included_n_reps, tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                ! Keep highest and lowest values
                call bottom_k_heap_push(bootstrapping_top_k_jsds(:, 1, i_study), n_bootstrapping_top_k_jsds, tmp_global_js_divergence(i_study))
                call top_k_heap_push(bootstrapping_top_k_jsds(:, 2, i_study), n_bootstrapping_top_k_jsds, tmp_global_js_divergence(i_study))
            end do
        end do

        ! Set confidence interval -> [largest small value, smallest large value]
        do concurrent (i_study = 1:n_studies) shared(bootstrapping_top_k_jsds, confidence_interval)
            confidence_interval(1, i_study) = bootstrapping_top_k_jsds(1, 1, i_study)
            confidence_interval(2, i_study) = bootstrapping_top_k_jsds(1, 2, i_study)
        end do

        call destroy_rng(rng)
    end subroutine bootstrap_histogram_helper

    !> Helper to test the histogram of a mean pmf that al bins have a minimum count
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
        do concurrent (i_point = 1:n_points, i_bin = 1:n_bins) shared(mean_pmf_counts, min, all_bins_have_min_count)
            if (mean_pmf_counts(i_bin, i_point) < min) then
                all_bins_have_min_count = .false.
            end if
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
            !! If all mean values are `NaN`, the range is `[1, min(n_genes_S, n_neighbors)]`
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
            ! Calculate overlap, is only negative if `left_max < right_min` -> `a_min < a_max < b_min < b_max` -> no overlap
            overlap_percent = max(0.0_real64, (left_max - right_min) / (a_max - a_min))
        end if
    end function compute_fractional_overlap_helper

    subroutine write_js_comp_test_trace_csv( &
        filename, &
        candidates_n_points_n_neighbors, &
        candidates_n_bins, &
        n_candidates_n_points_n_neighbors, &
        n_studies, &
        trace_neighbor_overlap_ok, &
        trace_min_bin_count_ok, &
        trace_min_mean_pmf_bin_count, &
        trace_confidence_interval, &
        trace_ci_overlap, &
        trace_exceeds_ci_overlap, &
        trace_plateau_found, &
        trace_candidates_kx &
    )

        use, intrinsic :: iso_fortran_env, only: int32, real64

        character(len=*), intent(in) :: filename

        integer(int32), intent(in) :: n_candidates_n_points_n_neighbors
        integer(int32), intent(in) :: n_studies

        integer(int32), dimension(2, n_candidates_n_points_n_neighbors), intent(in) :: &
            candidates_n_points_n_neighbors

        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: &
            candidates_n_bins

        logical, dimension(n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_neighbor_overlap_ok

        logical, dimension(n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_min_bin_count_ok

        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_min_mean_pmf_bin_count

        real(real64), dimension(2, n_studies, n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_confidence_interval

        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_ci_overlap

        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_exceeds_ci_overlap

        logical, dimension(n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_plateau_found

        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in) :: &
            trace_candidates_kx

        integer :: unit
        integer :: i_candidate
        integer :: i_study
        integer :: ios

        ! Debugging output
        write(*,*) "DEBUG: ENTERED write_js_comp_test_trace_csv"
        write(*,*) "DEBUG: filename = ", trim(filename)

        ! Open CSV file
        open( &
            newunit=unit, &
            file=filename, &
            status="replace", &
            action="write", &
            iostat=ios)

        if (ios /= 0) then
            write(*,*) "ERROR: Could not open JSCompTest trace CSV: ", trim(filename)
            return
        end if

        ! Header
        write(unit, '(A)') &
            'n_points,n_neighbors,n_bins,study,ci_lower,ci_upper,' // &
            'ci_overlap,min_bin_count,neighbor_overlap_ok,min_bin_count_ok,' // &
            'exceeds_ci_overlap,plateau_found,k_x'

        ! One row per (candidate x study)
        do i_candidate = 1, n_candidates_n_points_n_neighbors

            do i_study = 1, n_studies

                write(unit, '(I0,",",I0,",",I0,",",I0,",",ES24.16,",",ES24.16,",",ES24.16,",",I0,",",I0,",",I0,",",I0,",",I0,",",ES24.16)') &
                    candidates_n_points_n_neighbors(1, i_candidate), &
                    candidates_n_points_n_neighbors(2, i_candidate), &
                    candidates_n_bins(i_candidate), &
                    i_study, &
                    trace_confidence_interval(1, i_study, i_candidate), &
                    trace_confidence_interval(2, i_study, i_candidate), &
                    trace_ci_overlap(i_study, i_candidate), &
                    trace_min_mean_pmf_bin_count(i_candidate), &
                    merge(1, 0, trace_neighbor_overlap_ok(i_candidate)), &
                    merge(1, 0, trace_min_bin_count_ok(i_candidate)), &
                    trace_exceeds_ci_overlap(i_candidate), &
                    merge(1, 0, trace_plateau_found(i_candidate)), &
                    trace_candidates_kx(i_candidate)

            end do

        end do

        close(unit)

    end subroutine write_js_comp_test_trace_csv

    !> Writes a small one-row CSV with the run-level parameters, so the trace CSV can be interpreted standalone.
    subroutine write_js_comp_test_run_metadata_csv(filename, shared_residual_range, min_neighbor_overlap, succeeding_ci_overlap, min_count_per_mean_bin, n_bootstraps, join_method)
        character(len=*), intent(in) :: filename
        real(real64), intent(in) :: shared_residual_range
        real(real64), intent(in) :: min_neighbor_overlap
        real(real64), intent(in) :: succeeding_ci_overlap
        integer(int32), intent(in) :: min_count_per_mean_bin
        integer(int32), intent(in) :: n_bootstraps
        integer(int32), intent(in) :: join_method

        integer :: unit, ios

        open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Could not open run metadata CSV: ", trim(filename)
            return
        end if

        write(unit, '(A)') 'shared_residual_range,min_neighbor_overlap,succeeding_ci_overlap,min_count_per_mean_bin,n_bootstraps,join_method'
        write(unit, '(ES24.16,",",ES24.16,",",ES24.16,",",I0,",",I0,",",I0)') &
            shared_residual_range, min_neighbor_overlap, succeeding_ci_overlap, min_count_per_mean_bin, n_bootstraps, join_method

        close(unit)
    end subroutine write_js_comp_test_run_metadata_csv

    !> Writes the full per-bin gene counts of the mean pmf, one row per (candidate, point, bin).
    subroutine write_js_comp_test_bin_counts_csv(filename, candidates_n_points_n_neighbors, candidates_n_bins, n_candidates_n_points_n_neighbors, max_n_bins_all_candidates, max_n_points_candidate, trace_mean_pmf_counts_full, trace_neighbor_overlap_ok, trace_min_bin_count_ok, trace_x_star)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: n_candidates_n_points_n_neighbors
        integer(int32), intent(in) :: max_n_bins_all_candidates
        integer(int32), intent(in) :: max_n_points_candidate
        integer(int32), dimension(2, n_candidates_n_points_n_neighbors), intent(in) :: candidates_n_points_n_neighbors
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: candidates_n_bins
        integer(int32), dimension(max_n_bins_all_candidates, max_n_points_candidate, n_candidates_n_points_n_neighbors), intent(in) :: trace_mean_pmf_counts_full
        logical, dimension(n_candidates_n_points_n_neighbors), intent(in) :: trace_neighbor_overlap_ok
        logical, dimension(n_candidates_n_points_n_neighbors), intent(in) :: trace_min_bin_count_ok
        real(real64), dimension(max_n_points_candidate, n_candidates_n_points_n_neighbors), intent(in) :: trace_x_star
            !! Mean-expression reference point x_star(point_idx) for each candidate - lets bin sparsity be related to expression level

        integer :: unit, ios, i_candidate, i_point, i_bin

        open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Could not open bin counts CSV: ", trim(filename)
            return
        end if

        write(unit, '(A)') 'candidate,n_points,n_neighbors,point_idx,x_star,bin_idx,gene_count,neighbor_overlap_ok,min_bin_count_ok'

        do i_candidate = 1, n_candidates_n_points_n_neighbors
            do i_point = 1, candidates_n_points_n_neighbors(1, i_candidate)
                do i_bin = 1, candidates_n_bins(i_candidate)
                    write(unit, '(I0,",",I0,",",I0,",",I0,",",ES24.16,",",I0,",",I0,",",I0,",",I0)') &
                        i_candidate, candidates_n_points_n_neighbors(1, i_candidate), candidates_n_points_n_neighbors(2, i_candidate), &
                        i_point, trace_x_star(i_point, i_candidate), i_bin, trace_mean_pmf_counts_full(i_bin, i_point, i_candidate), &
                        merge(1, 0, trace_neighbor_overlap_ok(i_candidate)), &
                        merge(1, 0, trace_min_bin_count_ok(i_candidate))
                end do
            end do
        end do

        close(unit)
    end subroutine write_js_comp_test_bin_counts_csv

    !> Writes the relative effect-size change (Delta_t) diagnostics from the "Assessing JSD Plateau
    !| Stability Using Relative Effect-Size Change" section of the method description: J_{i,t}, CI width,
    !| relative CI width, Delta_{i,t}, and their per-candidate median/max across studies.
    subroutine write_js_comp_test_effect_size_csv(filename, candidates_n_points_n_neighbors, n_candidates_n_points_n_neighbors, n_studies, &
            trace_global_js_divergence, trace_confidence_interval, trace_delta, trace_delta_median, trace_delta_max)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: n_candidates_n_points_n_neighbors
        integer(int32), intent(in) :: n_studies
        integer(int32), dimension(2, n_candidates_n_points_n_neighbors), intent(in) :: candidates_n_points_n_neighbors
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(in) :: trace_global_js_divergence
        real(real64), dimension(2, n_studies, n_candidates_n_points_n_neighbors), intent(in) :: trace_confidence_interval
        real(real64), dimension(n_studies, n_candidates_n_points_n_neighbors), intent(in) :: trace_delta
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in) :: trace_delta_median
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in) :: trace_delta_max

        real(real64), parameter :: WRITE_EPSILON = 1.0e-12_real64
            !! Same role as DELTA_EPSILON in the helper: avoids division by ~0 when reporting relative CI width

        integer :: unit, ios, i_candidate, i_study
        real(real64) :: j_val, ci_width

        open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Could not open effect size CSV: ", trim(filename)
            return
        end if

        write(unit, '(A)') 'candidate,n_points,n_neighbors,study,J,ci_lower,ci_upper,ci_width,ci_width_relative,delta,delta_median,delta_max'

        do i_candidate = 1, n_candidates_n_points_n_neighbors
            do i_study = 1, n_studies
                j_val = trace_global_js_divergence(i_study, i_candidate)
                ci_width = trace_confidence_interval(2, i_study, i_candidate) - trace_confidence_interval(1, i_study, i_candidate)

                write(unit, '(I0,",",I0,",",I0,",",I0,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16)') &
                    i_candidate, candidates_n_points_n_neighbors(1, i_candidate), candidates_n_points_n_neighbors(2, i_candidate), &
                    i_study, j_val, &
                    trace_confidence_interval(1, i_study, i_candidate), trace_confidence_interval(2, i_study, i_candidate), &
                    ci_width, ci_width / max(j_val, WRITE_EPSILON), &
                    trace_delta(i_study, i_candidate), trace_delta_median(i_candidate), trace_delta_max(i_candidate)
            end do
        end do

        close(unit)
    end subroutine write_js_comp_test_effect_size_csv

    !> Writes the bin-count-estimation diagnostics (quartiles, n_pool, n_reps_neighborhood, half_bin_width,
    !| Sturges vs. Freedman-Diaconis candidate n_bins) computed during candidate generation in `_alloc`.
    !| quartile_25 / quartile_75 / n_pool are global (identical for every candidate); the rest vary per candidate.
    subroutine write_js_comp_test_bin_estimation_csv(filename, candidates_n_points_n_neighbors, candidates_kx, n_candidates_n_points_n_neighbors, &
            quartile_25, quartile_75, n_pool, n_reps_neighborhood, half_bin_width, n_bins_sturges, n_bins_freedman_diaconis, n_bins_final)
        character(len=*), intent(in) :: filename
        integer(int32), intent(in) :: n_candidates_n_points_n_neighbors
        integer(int32), dimension(2, n_candidates_n_points_n_neighbors), intent(in) :: candidates_n_points_n_neighbors
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in) :: candidates_kx
        real(real64), intent(in) :: quartile_25
        real(real64), intent(in) :: quartile_75
        integer(int32), intent(in) :: n_pool
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in) :: n_reps_neighborhood
        real(real64), dimension(n_candidates_n_points_n_neighbors), intent(in) :: half_bin_width
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: n_bins_sturges
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: n_bins_freedman_diaconis
        integer(int32), dimension(n_candidates_n_points_n_neighbors), intent(in) :: n_bins_final

        integer :: unit, ios, i_candidate

        open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: Could not open bin estimation CSV: ", trim(filename)
            return
        end if

        write(unit, '(A)') 'candidate,n_points,n_neighbors,k_x,quartile_25,quartile_75,n_pool,n_reps_neighborhood,half_bin_width,n_bins_sturges,n_bins_freedman_diaconis,n_bins_final'

        do i_candidate = 1, n_candidates_n_points_n_neighbors
            write(unit, '(I0,",",I0,",",I0,",",ES24.16,",",ES24.16,",",ES24.16,",",I0,",",ES24.16,",",ES24.16,",",I0,",",I0,",",I0)') &
                i_candidate, candidates_n_points_n_neighbors(1, i_candidate), candidates_n_points_n_neighbors(2, i_candidate), &
                candidates_kx(i_candidate), quartile_25, quartile_75, n_pool, &
                n_reps_neighborhood(i_candidate), half_bin_width(i_candidate), &
                n_bins_sturges(i_candidate), n_bins_freedman_diaconis(i_candidate), n_bins_final(i_candidate)
        end do

        close(unit)
    end subroutine write_js_comp_test_bin_estimation_csv

end module tox_data_integration_js_comp_test