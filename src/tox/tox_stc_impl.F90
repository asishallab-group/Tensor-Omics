#include <src/macros.h>

#define CM_OBSERVABLE_COUNT 5

!> Shape-truthful clustering with a global or per-seed growth radius.
!| Compute density labels, select seeds, grow their ensembles, and merge overlapping ensembles.
!| Both radius variants use the same density labels, seed ordering, growth rule and merge rule.
module tox_stc_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err_once, is_err, validate_all_in_range_real, ERR_INVALID_INPUT
    use tox_gene_centroids_impl, only: mean_vector_impl
    use tox_euclidean_distance_impl, only: euclidean_distance_impl
    use f42_sort_impl, only: sort_array_heapsort, init_perm
    use f42_stats_impl, only: calc_quantile_impl
    use f42_kd_tree_impl, only: vicinity_vectors_helper, vicinity_vectors_count_helper
    M_IMPLICIT_NONE
    private

    integer(int32), parameter, public :: MODE_GLOBAL_RADIUS = 1_int32
    integer(int32), parameter, public :: MODE_SEED_RADIUS = 2_int32

    integer(int32), parameter, public :: STC_STOP_FIXED_POINT = 1_int32
        !! Growth accepted at least one batch, then found no further candidates.
    integer(int32), parameter, public :: STC_STOP_REJECT_AFTER_ACCEPT = 2_int32
        !! A candidate batch was rejected after at least one earlier batch was accepted.
    integer(int32), parameter, public :: STC_STOP_NEVER_ACCEPTED = 3_int32
        !! The first candidate batch was rejected.
    integer(int32), parameter, public :: STC_STOP_NO_CANDIDATES = 4_int32
        !! The seed had no growth candidates.

    public :: calculate_density_radius_impl
    public :: calculate_labels_as_density_impl
    public :: compute_ambient_density_stats_impl
    public :: merge_ensembles_impl
    public :: identify_ensemble_seeds_impl
    public :: grow_ensemble_impl
    public :: compute_ensemble_observable_impl
    public :: accept_ensemble_impl
    public :: obtain_ensembles_impl
    public :: compute_ensemble_observable_helper
    public :: accept_ensemble_helper
    public :: grow_single_seed_helper

contains

    !> summary: Compute the density-label radius from distances to the data centroid.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine calculate_density_radius_impl(vectors, n_dimensions, n_vectors, &
                                                    tmp_mean_vec, tmp_distances, tmp_perm, &
                                                    mean_to_other_vecs_dist_quant, radius, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(out) :: tmp_mean_vec(n_dimensions)
        !! Preallocated workspace for the mean vector
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated workspace for calculated distances
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated workspace for permutation indices
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
        !! Quantile of the distances to the data centroid that becomes the density radius
        !! DM_MIN(0.0_real64)
        !! DM_MAX(1.0_real64)
        !! DM_DEFAULT(0.15_real64)
        real(real64), intent(out) :: radius
        !! The resulting density search radius
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64)   :: actual_quant
        integer(int32) :: i_vec

        call set_ok(ierr)

        M_DEFAULT_VAL(mean_to_other_vecs_dist_quant, actual_quant, 0.15_real64)

        ! Initializing perm
        call init_perm(tmp_perm)

        ! Calculating mean vector
        call mean_vector_impl(vectors, n_dimensions, n_vectors, tmp_perm, n_vectors, tmp_mean_vec)

        ! Calculating each vector distance to mean vector
        do concurrent(i_vec=1:n_vectors) shared(vectors, tmp_mean_vec, n_dimensions, tmp_distances)
            call euclidean_distance_impl(tmp_mean_vec, vectors(:, i_vec), &
                                         n_dimensions, tmp_distances(i_vec))
        end do

        ! Finite input coordinates can still overflow in the distance computation.
        call validate_all_in_range_real(tmp_distances, n_vectors, ierr)
        if (is_err(ierr)) return

        ! Sort the permutation by distance.
        call sort_array_heapsort(tmp_distances, tmp_perm)

        ! Extract the adaptive search radius corresponding to the specified distance percentile.
        call calc_quantile_impl(tmp_distances, n_vectors, tmp_perm, actual_quant, radius)

    end subroutine calculate_density_radius_impl

    !> summary: Count points in each Euclidean neighborhood to label local density.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine calculate_labels_as_density_impl(vectors, n_dimensions, n_vectors, r, &
                                                       dimension_order, kd_indices, n_tiles, tmp_stack, &
                                                       label_densities, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
        !! Label-sphere search radius
        !! DM_MIN(0.0_real64)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking tree split axes by variance
        !! DM_MIN(1_int32)
        !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        !! DM_MIN(1_int32)
        !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: n_tiles
        !! Number of concurrent vector tiles, at least 1. Capped at the number of points.
        !! DM_MIN(1_int32)
        integer(int32), intent(inout) :: tmp_stack(M_KD_STACK_ENTRY_SIZE, M_KD_TRAVERSAL_STACK_DEPTH, min(n_tiles, n_vectors))
        !! Preallocated traversal stack per tile, reused serially across its vectors
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output array storing generated density scalars
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: i_tile, i_vec, neighbor_count, vector_start, vector_end, tile_base, tile_rem, active_tiles

        call set_ok(ierr)

        ! Balanced contiguous tiles.
        active_tiles = min(n_tiles, n_vectors)
        tile_base = n_vectors/active_tiles
        tile_rem = mod(n_vectors, active_tiles)

        do concurrent(i_tile=1:active_tiles) &
            shared(vectors, n_dimensions, n_vectors, r, dimension_order, kd_indices, &
                   tmp_stack, label_densities, tile_base, tile_rem) &
            local(i_vec, neighbor_count, vector_start, vector_end)

            ! Avoid forming n_vectors + 1 for empty tiles at the int32 limit.
            if (tile_base == 0_int32 .and. i_tile > tile_rem) cycle
            vector_start = (i_tile - 1_int32)*tile_base + min(i_tile - 1_int32, tile_rem) + 1_int32
            vector_end = vector_start - 1_int32 + tile_base
            if (i_tile <= tile_rem) vector_end = vector_end + 1_int32

            do i_vec = vector_start, vector_end
                call vicinity_vectors_count_helper(vectors(:, i_vec), vectors, n_dimensions, n_vectors, r, &
                                                   dimension_order, kd_indices, tmp_stack(:, :, i_tile), &
                                                   neighbor_count)
                label_densities(i_vec) = real(neighbor_count, real64)
            end do
        end do

    end subroutine calculate_labels_as_density_impl

    !> summary: Select seeds in descending density order and compute their local radii.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine identify_ensemble_seeds_impl(vectors, n_dimensions, n_vectors, density_labels, &
                                                   dimension_order, kd_indices, k_seeding, &
                                                   tmp_perm, tmp_distances, tmp_stack, &
                                                   tmp_visited_mask, tmp_newly_covered_mask, &
                                                   sorted_perm, n_seeds, seed_mask, seed_radii, &
                                                   coverage_quant)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        !! DM_MIN(2_int32)
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for all ambient vectors
        !! DM_MIN(1.0_real64)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the K-D tree structure
        !! DM_MIN(1_int32)
        !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! K-D tree index sequence array
        !! DM_MIN(1_int32)
        !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: k_seeding
        !! Number of nearest neighbors used to estimate each seed's local coverage radius
        !! DM_MIN(1_int32)
        !! DM_MAX(n_vectors - 1_int32)
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated permutation workspace
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated candidate-to-vector distance workspace
        integer(int32), intent(inout) :: tmp_stack(M_KD_STACK_ENTRY_SIZE, M_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated K-D tree traversal stack
        logical(c_bool), intent(out) :: tmp_visited_mask(n_vectors)
        !! Preallocated workspace tracking vectors already covered during seeding
        logical(c_bool), intent(out) :: tmp_newly_covered_mask(n_vectors)
        !! Preallocated workspace for the current seed's coverage query
        integer(int32), intent(out) :: sorted_perm(n_vectors)
        !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
        !! Number of selected seed vectors
        logical(c_bool), intent(out) :: seed_mask(n_vectors)
        !! Logical mask identifying selected seed vectors
        real(real64), intent(out) :: seed_radii(n_vectors)
        !! Per-seed growth radius, stored at each seed's own vector index; zero at non-seeds
        real(real64), intent(in), optional :: coverage_quant
        !! Quantile of each seed's k-nearest-neighbor distances used as its coverage and growth radius
        !! DM_MIN(0.0_real64)
        !! DM_MAX(1.0_real64)
        !! DM_DEFAULT(0.5_real64)

        integer(int32) :: i_vec, i_rank, candidate_idx, self_pos
        real(real64) :: coverage_radius, actual_coverage_quant

        M_DEFAULT_VAL(coverage_quant, actual_coverage_quant, 0.5_real64)
        call init_perm(tmp_perm)

        call sort_array_heapsort(density_labels, tmp_perm)

        ! Heapsort provides ascending density order. Reverse the permutation.
        sorted_perm = tmp_perm(n_vectors:1:-1)

        seed_mask = .false.
        seed_radii = 0.0_real64
        tmp_visited_mask = .false.
        n_seeds = 0_int32

        do i_rank = 1, n_vectors
            candidate_idx = sorted_perm(i_rank)

            if (tmp_visited_mask(candidate_idx)) cycle

            ! Highest-density unvisited vector becomes the next seed.
            seed_mask(candidate_idx) = .true.
            tmp_visited_mask(candidate_idx) = .true.
            n_seeds = n_seeds + 1_int32

            ! Calculate the distance from the current seed to every ambient vector.
            do concurrent(i_vec=1:n_vectors) &
                shared(vectors, n_dimensions, candidate_idx, tmp_distances)

                call euclidean_distance_impl(vectors(:, candidate_idx), vectors(:, i_vec), &
                                             n_dimensions, tmp_distances(i_vec))
            end do

            call init_perm(tmp_perm)

            call sort_array_heapsort(tmp_distances, tmp_perm)

            ! Exclude the seed itself from its k_seeding nearest neighbors.
            self_pos = 0_int32

            do i_vec = 1, n_vectors
                if (tmp_perm(i_vec) == candidate_idx) then
                    self_pos = i_vec
                    exit
                end if
            end do

            if ((self_pos >= 1_int32) .and. (self_pos <= k_seeding)) then
                do i_vec = self_pos, k_seeding
                    tmp_perm(i_vec) = tmp_perm(i_vec + 1_int32)
                end do
            end if

            ! Use the chosen k-neighbor distance percentile as this seed's local seeding and growth radius.
            call calc_quantile_impl(tmp_distances, n_vectors, tmp_perm, &
                                      actual_coverage_quant, coverage_radius, n_considered=k_seeding)

            seed_radii(candidate_idx) = coverage_radius

            ! Mark the region represented by the current seed as visited.
            call vicinity_vectors_helper(vectors(:, candidate_idx), vectors, &
                                         n_dimensions, n_vectors, coverage_radius, &
                                         dimension_order, kd_indices, tmp_stack, &
                                         tmp_newly_covered_mask)

            tmp_visited_mask = tmp_visited_mask .or. tmp_newly_covered_mask
        end do

    end subroutine identify_ensemble_seeds_impl

    !> summary: Compute the median and median absolute deviation of the density labels.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine compute_ambient_density_stats_impl(density_labels, n_vectors, &
                                                         tmp_perm, tmp_abs_diff, &
                                                         median_ambient, mad_ambient)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors, from `calculate_labels_as_density`.
        !! DM_MIN(1.0_real64)
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(out) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array for the absolute deviations
        real(real64), intent(out) :: median_ambient
        !! Ambient median density median_{i=1..N}(rho_i)
        real(real64), intent(out) :: mad_ambient
        !! Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)

        integer(int32) :: i_vec

        ! 1. Calculate ambient median density: median_{i=1..N}(rho_i)
        call init_perm(tmp_perm)
        call sort_array_heapsort(density_labels, tmp_perm)
        call calc_quantile_impl(density_labels, n_vectors, tmp_perm, &
                                  0.5_real64, median_ambient)

        ! 2. Compute absolute deviation from ambient median for all vectors
        do concurrent(i_vec=1:n_vectors) shared(tmp_abs_diff, density_labels, median_ambient)
            tmp_abs_diff(i_vec) = abs(density_labels(i_vec) - median_ambient)
        end do

        ! 3. Compute MAD_ambient = median_{i=1..N}(|rho_i - median_ambient|)
        call init_perm(tmp_perm)
        call sort_array_heapsort(tmp_abs_diff, tmp_perm)
        call calc_quantile_impl(tmp_abs_diff, n_vectors, tmp_perm, &
                                  0.5_real64, mad_ambient)

    end subroutine compute_ambient_density_stats_impl

    !> Filter the surface by density compatibility with the current ensemble.
    pure subroutine build_ensemble_trial_helper(ensemble_mask, n_vectors, density_labels, alpha_mad, &
                                                mad_ambient, surface_mask, tmp_perm, tmp_abs_diff, trial_mask)
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Accepted membership, always retained in the trial
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors
        real(real64), intent(in) :: alpha_mad
        !! Multiplier for the ambient MAD density tolerance
        real(real64), intent(in) :: mad_ambient
        !! Precomputed ambient median absolute deviation
        logical(c_bool), intent(in) :: surface_mask(n_vectors)
        !! Spatial neighborhood union; never filtered by density compatibility
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Sorting workspace; any frontier indices must be consumed before this call
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Workspace for accepted member densities
        logical(c_bool), intent(out) :: trial_mask(n_vectors)
        !! Complete trial membership: accepted members and density-compatible spatial neighbors

        integer(int32) :: i_vec, k, n_active
        real(real64) :: ensemble_center_density, max_allowed_dev

        max_allowed_dev = alpha_mad*mad_ambient
        n_active = count(ensemble_mask)
        if (n_active == 0_int32) then
            trial_mask = .false.
            return
        end if

        k = 0
        do i_vec = 1, n_vectors
            if (ensemble_mask(i_vec)) then
                k = k + 1
                tmp_abs_diff(k) = density_labels(i_vec)
            end if
        end do
        call init_perm(tmp_perm(1:k))
        call sort_array_heapsort(tmp_abs_diff(1:k), tmp_perm(1:k))
        call calc_quantile_impl(tmp_abs_diff, n_vectors, tmp_perm, 0.5_real64, &
                                  ensemble_center_density, n_considered=k)

        do concurrent(i_vec=1:n_vectors) &
            shared(ensemble_mask, surface_mask, trial_mask, density_labels, ensemble_center_density, max_allowed_dev)
            trial_mask(i_vec) = ensemble_mask(i_vec) .or. &
                                (surface_mask(i_vec) .and. &
                                 abs(density_labels(i_vec) - ensemble_center_density) <= max_allowed_dev)
        end do

    end subroutine build_ensemble_trial_helper

    !> Update the five observables in the history window.
    !| The window is the column count of `observables`, the same width `accept_ensemble_helper` reads.
    pure subroutine compute_ensemble_observable_helper(ensemble_mask, density_labels, n_vectors, &
                                                       candidate_count, current_iter, observables)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Selection mask for active ensemble members
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors, from `calculate_labels_as_density`.
        integer(int32), intent(in) :: candidate_count
        !! Number of candidates added in current growth step
        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(inout) :: observables(:, :)
        !! Observable history; its column count is the history window

        integer(int32) :: i_vec, active_count, max_cols, zero_count
        real(real64)   :: sum_arithmetic, sum_reciprocal, rho_arith, rho_harm
        real(real64)   :: obs_1, obs_2, obs_3, obs_4, obs_5

        max_cols = size(observables, dim=2, kind=int32)

        active_count = count(ensemble_mask)

        if (active_count == 0_int32) then
            obs_1 = 0.0_real64
            obs_2 = 0.0_real64
            obs_3 = 0.0_real64
            obs_4 = 0.0_real64
        else
            sum_arithmetic = 0.0_real64
            sum_reciprocal = 0.0_real64
            zero_count = 0_int32

            do concurrent(i_vec=1:n_vectors) &
                shared(ensemble_mask, density_labels) &
                reduce(+:sum_arithmetic, sum_reciprocal, zero_count)

                if (ensemble_mask(i_vec)) then
                    sum_arithmetic = sum_arithmetic + density_labels(i_vec)

                    if (density_labels(i_vec) > 0.0_real64) then
                        sum_reciprocal = sum_reciprocal + (1.0_real64/density_labels(i_vec))
                    else
                        zero_count = zero_count + 1_int32
                    end if
                end if
            end do

            ! 1. Arithmetic Mean Density
            rho_arith = sum_arithmetic/real(active_count, real64)

            ! 2. Harmonic Mean Density (strictly 0 if any active element is <= 0).
            if (zero_count == 0_int32 .and. sum_reciprocal > 0.0_real64) then
                rho_harm = real(active_count, real64)/sum_reciprocal
            else
                rho_harm = 0.0_real64
            end if

            obs_1 = rho_arith
            obs_2 = rho_harm

            ! 3. Density Heterogeneity H_rho.
            if (rho_harm > 0.0_real64) then
                obs_3 = rho_arith/rho_harm
            else
                obs_3 = 1.0_real64
            end if

            ! 4. Active Ensemble Size
            obs_4 = real(active_count, real64)
        end if

        ! 5. Candidate Count (recorded independently of active_count)
        obs_5 = real(candidate_count, real64)

        ! Store scalar observables matching iteration index and history window
        if (current_iter <= max_cols) then
            observables(1, current_iter) = obs_1
            observables(2, current_iter) = obs_2
            observables(3, current_iter) = obs_3
            observables(4, current_iter) = obs_4
            observables(5, current_iter) = obs_5
        else
            observables(:, 1:max_cols - 1_int32) = observables(:, 2:max_cols)
            observables(1, max_cols) = obs_1
            observables(2, max_cols) = obs_2
            observables(3, max_cols) = obs_3
            observables(4, max_cols) = obs_4
            observables(5, max_cols) = obs_5
        end if

    end subroutine compute_ensemble_observable_helper

    !> Apply the log fold-change acceptance rule.
    !| Every caller passes a history of at least two columns, so two consecutive steps always exist.
    pure subroutine accept_ensemble_helper(observables, current_iter, alpha_accept, &
                                           is_accepted)

        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(in) :: observables(:, :)
        !! Observable history [5 x window], at least two columns wide
        real(real64), intent(in), optional :: alpha_accept
        !! Log2 fold-change threshold; absent means 0.5
        logical(c_bool), intent(out) :: is_accepted
        !! Output flag indicating step acceptance decision

        real(real64)   :: actual_alpha, h_prev, h_curr, log2_fold_change
        integer(int32) :: max_cols, col_curr, col_prev

        M_DEFAULT_VAL(alpha_accept, actual_alpha, 0.5_real64)

        ! Initial seed iteration (iter == 1) is always accepted
        if (current_iter <= 1_int32) then
            is_accepted = .true.
            return
        end if

        max_cols = size(observables, dim=2, kind=int32)

        ! Determine column indices for current and previous iterations
        if (current_iter <= max_cols) then
            col_curr = current_iter
            col_prev = current_iter - 1_int32
        else
            ! Sliding window shifted: last column is current, second-to-last is previous
            col_curr = max_cols
            col_prev = max_cols - 1_int32
        end if

        ! Evaluate log2 fold change of Density Heterogeneity H_rho (Row 3)
        h_prev = observables(3, col_prev)
        h_curr = observables(3, col_curr)

        if ((h_prev > 0.0_real64) .and. (h_curr > 0.0_real64)) then
            log2_fold_change = log(h_curr/h_prev)/log(2.0_real64)

            ! Reject if fold change meets or exceeds acceptance threshold alpha_accept
            if (abs(log2_fold_change) >= actual_alpha) then
                is_accepted = .false.
            else
                is_accepted = .true.
            end if
        else
            is_accepted = .true.
        end if

    end subroutine accept_ensemble_helper

    !> summary: Grow one raw ensemble from each seed, using a global radius or per-seed radii.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine obtain_ensembles_impl(vectors, n_dimensions, n_vectors, dimension_order, &
                                            kd_indices, density_labels, seed_indices, n_seeds, &
                                            growth_mode, r, seed_radii, alpha_mad, alpha_accept, t_observables, n_tiles, &
                                            tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, &
                                            tmp_abs_diff, tmp_observables, tmp_current_mask, &
                                            ensemble_matrix, stop_reasons, &
                                            mad_ambient, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        integer(int32), intent(in) :: n_seeds
        !! Number of seed vectors, from `identify_ensemble_seeds`.
        !! DM_MIN(0_int32)
        !! DM_MAX(n_vectors)
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions for KD-tree
        !! DM_MIN(1_int32)
        !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        !! DM_MIN(1_int32)
        !! DM_MAX(n_vectors)
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors, from `calculate_labels_as_density`.
        !! DM_MIN(1.0_real64)
        integer(int32), intent(in) :: seed_indices(n_seeds)
        !! Indices of the seed vectors to grow from, from `identify_ensemble_seeds`.
        !! DM_MIN(1_int32)
        !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: growth_mode
        !! Radius to use for ensemble growth.
        !!
        !! | Mode | Value | Procedure |
        !! |-------------|-------|-----------|
        !! | Global radius | [[tox_stc_impl(module):MODE_GLOBAL_RADIUS(variable)]] | obtain_ensembles_global_radius |
        !! | Per-seed radius | [[tox_stc_impl(module):MODE_SEED_RADIUS(variable)]] | obtain_ensembles_seed_radius |
        real(real64), intent(in), optional :: r
        !! Uniform growth radius for every seed.
        !! DM_REQUIRED_IF_MODE(growth_mode, tox_stc_impl, MODE_GLOBAL_RADIUS)
        !! DM_MIN(0.0_real64)
        real(real64), intent(in), optional :: seed_radii(n_seeds)
        !! Per-seed growth radius, allowing tighter growth in dense regions and wider growth in sparse regions.
        !! DM_MIN(0.0_real64)
        !! DM_REQUIRED_IF_MODE(growth_mode, tox_stc_impl, MODE_SEED_RADIUS)
        real(real64), intent(in), optional :: alpha_mad
        !! MAD multiplier for density compatibility threshold
        !! DM_MIN(0.0_real64)
        !! DM_DEFAULT(0.5_real64)
        real(real64), intent(in), optional :: alpha_accept
        !! Acceptance log2 fold-change threshold
        !! DM_MIN(0.0_real64)
        !! DM_DEFAULT(0.5_real64)
        integer(int32), intent(in) :: t_observables
        !! History window in iterations, at least two because acceptance compares two consecutive
        !! steps. Pass n_vectors to retain the full history.
        !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_tiles
        !! Number of concurrent seed tiles, at least 1. Capped at the number of seeds when nonzero.
        !! DM_MIN(1_int32)
        integer(int32), intent(inout) :: tmp_stack(M_KD_STACK_ENTRY_SIZE, M_KD_TRAVERSAL_STACK_DEPTH, &
                                                   min(n_tiles, max(1_int32, n_seeds)))
        !! Workspace stack for tree traversal per tile
        logical(c_bool), intent(out) :: tmp_vicinity_mask(n_vectors, min(n_tiles, max(1_int32, n_seeds)))
        !! Workspace mask holding one member's neighborhood query result per tile
        logical(c_bool), intent(out) :: tmp_surface_mask(n_vectors, min(n_tiles, max(1_int32, n_seeds)))
        !! Cached spatial neighborhood union per tile, reset for each seed
        integer(int32), intent(inout) :: tmp_perm(n_vectors, min(n_tiles, max(1_int32, n_seeds)))
        !! Sorting workspace, reused for newly accepted member indices between growth rounds per tile
        real(real64), intent(out) :: tmp_abs_diff(n_vectors, min(n_tiles, max(1_int32, n_seeds)))
        !! Workspace array for deviation calculations per tile
        real(real64), intent(out) :: tmp_observables(CM_OBSERVABLE_COUNT, max(2_int32, min(t_observables, n_vectors)), &
                                                     min(n_tiles, max(1_int32, n_seeds)))
        !! Observable history per tile. Growth adds at least one member per iteration, so more than
        !! n_vectors columns are never needed; at least two are kept so acceptance can compare steps.
        logical(c_bool), intent(out) :: tmp_current_mask(n_vectors, min(n_tiles, max(1_int32, n_seeds)))
        !! Workspace holding accepted seed ensemble membership per tile
        logical(c_bool), intent(out) :: ensemble_matrix(n_vectors, n_seeds)
        !! Output matrix storing unmerged grown ensemble masks
        integer(int32), intent(out) :: stop_reasons(n_seeds)
        !! Per-seed stop code. 1: fixed point, 2: rejection after acceptance,
        !! 3: first step rejected, 4: no candidates.
        real(real64), intent(out) :: mad_ambient
        !! Median absolute deviation of `density_labels`, which scales the `alpha_mad` density tolerance.
        !! Zero when there are no seeds; otherwise positive on success, because a zero value is an error.
        integer(int32), intent(out) :: ierr
        !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
        !! half the labels are equal), which leaves no density tolerance for growth.

        real(real64) :: actual_alpha_mad, actual_alpha_accept, median_ambient
        integer(int32) :: active_tiles, i_tile, i_seed, seed_start, seed_end, tile_base, tile_rem
        real(real64) :: growth_radius
        logical(c_bool) :: use_seed_radii

        call set_ok(ierr)
        if (n_seeds == 0_int32) then
            mad_ambient = 0.0_real64
            return
        end if

        call compute_ambient_density_stats_impl(density_labels, n_vectors, tmp_perm(:, 1), tmp_abs_diff(:, 1), &
                                                median_ambient, mad_ambient)
        if (mad_ambient <= 0.0_real64) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=6_int32)
            return
        end if

        M_DEFAULT_VAL(alpha_mad, actual_alpha_mad, 0.5_real64)
        M_DEFAULT_VAL(alpha_accept, actual_alpha_accept, 0.5_real64)
        active_tiles = min(n_tiles, n_seeds)
        use_seed_radii = logical(growth_mode == MODE_SEED_RADIUS, c_bool)
        ensemble_matrix = .false.

        ! Distribute seeds evenly across tiles, with the first tile_rem tiles receiving one extra seed.
        tile_base = n_seeds/active_tiles
        tile_rem = mod(n_seeds, active_tiles)

        ! Process independent seed tiles concurrently, with seeds within each tile handled serially.
        do concurrent(i_tile=1:active_tiles) &
            shared(vectors, n_dimensions, n_vectors, dimension_order, kd_indices, &
                   density_labels, seed_indices, n_seeds, use_seed_radii, r, seed_radii, actual_alpha_mad, mad_ambient, &
                   actual_alpha_accept, tile_base, tile_rem, tmp_stack, tmp_vicinity_mask, &
                   tmp_surface_mask, tmp_perm, tmp_abs_diff, tmp_observables, &
                   tmp_current_mask, ensemble_matrix, stop_reasons) &
            local(i_seed, seed_start, seed_end, growth_radius)

            if (i_tile <= tile_rem) then
                seed_start = (i_tile - 1_int32)*(tile_base + 1_int32) + 1_int32
                seed_end = seed_start + tile_base
            else
                seed_start = tile_rem*(tile_base + 1_int32) + &
                             (i_tile - tile_rem - 1_int32)*tile_base + 1_int32
                seed_end = seed_start + tile_base - 1_int32
            end if

            do i_seed = seed_start, seed_end

                if (use_seed_radii) then
                    growth_radius = seed_radii(i_seed)
                else
                    growth_radius = r
                end if

                call grow_single_seed_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                             kd_indices, density_labels, seed_indices(i_seed), &
                                             growth_radius, actual_alpha_mad, mad_ambient, actual_alpha_accept, &
                                             tmp_stack(:, :, i_tile), tmp_vicinity_mask(:, i_tile), &
                                             tmp_surface_mask(:, i_tile), tmp_perm(:, i_tile), &
                                             tmp_abs_diff(:, i_tile), tmp_observables(:, :, i_tile), &
                                             tmp_current_mask(:, i_tile), &
                                             ensemble_matrix(:, i_seed), stop_reasons(i_seed))

            end do

        end do

    end subroutine obtain_ensembles_impl

    !> Grow one seed and record why it stopped.
    pure subroutine grow_single_seed_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                            kd_indices, density_labels, seed_idx, &
                                            r, alpha_mad, mad_ambient, alpha_accept, &
                                            tmp_stack, tmp_vicinity_mask, tmp_surface_mask, &
                                            tmp_perm, tmp_abs_diff, tmp_observables, &
                                            current_mask, out_mask, stop_reason)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions for KD-tree
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors, from `calculate_labels_as_density`.
        integer(int32), intent(in) :: seed_idx
        !! Index of starting seed vector
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        real(real64), intent(in) :: alpha_mad
        !! MAD multiplier for density compatibility threshold
        real(real64), intent(in) :: mad_ambient
        !! Ambient MAD precalculated via the ambient density statistics routine.
        real(real64), intent(in) :: alpha_accept
        !! Acceptance log2 fold-change threshold
        integer(int32), intent(inout) :: tmp_stack(M_KD_STACK_ENTRY_SIZE, M_KD_TRAVERSAL_STACK_DEPTH)
        !! Workspace stack for a single tree traversal
        logical(c_bool), intent(inout) :: tmp_vicinity_mask(n_vectors)
        !! Workspace mask holding one member's neighborhood query result
        logical(c_bool), intent(inout) :: tmp_surface_mask(n_vectors)
        !! Cached geometric union for this seed; reset once and never density-filtered
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Newly accepted indices between rounds; reused for median sorting after their queries
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Workspace array holding the active member densities
        real(real64), intent(inout) :: tmp_observables(:, :)
        !! Observable history, at least two columns; its column count is the history window
        logical(c_bool), intent(inout) :: current_mask(n_vectors)
        !! Accepted ensemble membership; modified only after a trial passes acceptance
        logical(c_bool), intent(out) :: out_mask(n_vectors)
        !! Final accepted membership; reused as trial membership until growth terminates
        integer(int32), intent(out) :: stop_reason
        !! Reason growth terminated, one of the `STC_STOP_*` constants

        integer(int32) :: iter, prev_count, curr_count, candidate_count, i_vec, i_frontier, frontier_idx, n_frontier
        logical(c_bool) :: is_growing, is_accepted

        current_mask = .false.
        current_mask(seed_idx) = .true.
        tmp_observables = 0.0_real64
        tmp_surface_mask = .false.
        n_frontier = 1_int32
        tmp_perm(1) = seed_idx

        iter = 1_int32
        call compute_ensemble_observable_helper(current_mask, density_labels, &
                                                n_vectors, 0_int32, iter, &
                                                tmp_observables)

        is_growing = .true.
        do while (is_growing)
            prev_count = count(current_mask)

            ! Query each accepted vector once before median sorting reuses tmp_perm; query new candidates only after acceptance.
            do i_frontier = 1, n_frontier
                frontier_idx = tmp_perm(i_frontier)
                call vicinity_vectors_helper(vectors(:, frontier_idx), vectors, n_dimensions, n_vectors, r, &
                                             dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask)
                tmp_surface_mask = tmp_surface_mask .or. tmp_vicinity_mask
            end do

            call build_ensemble_trial_helper(current_mask, n_vectors, density_labels, alpha_mad, &
                                             mad_ambient, tmp_surface_mask, tmp_perm, tmp_abs_diff, out_mask)

            curr_count = count(out_mask)
            candidate_count = curr_count - prev_count

            if (candidate_count == 0_int32) then
                ! iter is still 1 only on the very first pass, so the seed never grew at all
                if (iter == 1_int32) then
                    stop_reason = STC_STOP_NO_CANDIDATES
                else
                    stop_reason = STC_STOP_FIXED_POINT
                end if
                is_growing = .false.
                exit
            end if

            iter = iter + 1_int32

            call compute_ensemble_observable_helper(out_mask, density_labels, &
                                                    n_vectors, candidate_count, iter, &
                                                    tmp_observables)

            call accept_ensemble_helper(tmp_observables, iter, &
                                        alpha_accept, is_accepted)

            if (.not. is_accepted) then
                ! iter counts candidate-producing growth rounds; iter == 2 means no previous batch was accepted.
                if (iter == 2_int32) then
                    stop_reason = STC_STOP_NEVER_ACCEPTED
                else
                    stop_reason = STC_STOP_REJECT_AFTER_ACCEPT
                end if

                is_growing = .false.
            else
                ! Commit and collect only new members for next round's spatial queries.
                n_frontier = 0_int32
                do i_vec = 1, n_vectors
                    if (out_mask(i_vec) .and. .not. current_mask(i_vec)) then
                        current_mask(i_vec) = .true.
                        n_frontier = n_frontier + 1_int32
                        tmp_perm(n_frontier) = i_vec
                    end if
                end do
            end if
        end do

        out_mask = current_mask

    end subroutine grow_single_seed_helper

    !> Find a union-find component root and compress the path.
    pure subroutine find_root_helper(parent, n_nodes, node, root)

        integer(int32), intent(in) :: n_nodes
        !! Number of nodes in the forest
        integer(int32), intent(inout) :: parent(n_nodes)
        !! Parent pointers; flattened towards the root as a side effect of the lookup
        integer(int32), intent(in) :: node
        !! Node whose component root is wanted
        integer(int32), intent(out) :: root
        !! Representative of the component containing `node`

        integer(int32) :: walker, next_node

        root = node
        do while (parent(root) /= root)
            root = parent(root)
        end do

        ! Path compression keeps later lookups near constant time and prevents long union chains from restoring quadratic behaviour.
        walker = node
        do while (parent(walker) /= root)
            next_node = parent(walker)
            parent(walker) = root
            walker = next_node
        end do

    end subroutine find_root_helper

    !> summary: Merge intersecting ensembles by overlap coefficient and discard final singletons.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine merge_ensembles_impl(raw_masks, n_vectors, n_seeds, &
                                           min_overlap_coefficient, merged_masks, &
                                           tmp_active_flag, tmp_parent, n_ensembles, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of grown raw ensemble columns
        !! DM_MIN(0_int32)
        !! DM_MAX(n_vectors)
        real(real64), intent(in), optional :: min_overlap_coefficient
        !! Minimum Overlap Coefficient required to merge two ensembles
        !! DM_MIN(0.0_real64)
        !! DM_MAX(1.0_real64)
        !! DM_DEFAULT(0.9_real64)
        logical(c_bool), intent(in) :: raw_masks(n_vectors, n_seeds)
        !! Matrix of raw unmerged boolean ensemble masks
        logical(c_bool), intent(out) :: merged_masks(n_vectors, n_seeds)
        !! Output matrix storing merged boolean ensemble masks
        !! DM_RESULT_SIZE_IS(n_ensembles)
        logical(c_bool), intent(out) :: tmp_active_flag(n_seeds)
        !! Preallocated workspace marking the component roots that survive as ensembles
        integer(int32), intent(inout) :: tmp_parent(n_seeds)
        !! Preallocated union-find parent pointers
        integer(int32), intent(out) :: n_ensembles
        !! Final count of merged non-singleton ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        integer(int32) :: i, j, shared_count, size_i, size_j, root_i, root_j

        real(real64) :: actual_min_overlap

        call set_ok(ierr)
        M_DEFAULT_VAL(min_overlap_coefficient, actual_min_overlap, 0.9_real64)

        if (n_seeds == 0_int32) then
            n_ensembles = 0_int32
            merged_masks = .false.
            return
        end if

        do i = 1, n_seeds
            tmp_parent(i) = i
        end do

        ! Merge connected components of the raw-mask overlap graph in one O(S^2*N) union-find pass.
        do i = 1, n_seeds
            size_i = count(raw_masks(:, i))
            if (size_i == 0_int32) cycle

            do j = i + 1, n_seeds
                shared_count = count(raw_masks(:, i) .and. raw_masks(:, j))
                ! A disjoint pair never merges, whatever the threshold; this also keeps min() below nonzero.
                if (shared_count < 1_int32) cycle

                ! Only paid for pairs that actually intersect, so the extra count stays off the common path.
                size_j = count(raw_masks(:, j))

                ! Overlap Coefficient: |E_i intersection E_j| / min(|E_i|, |E_j|)
                if (real(shared_count, real64) / real(min(size_i, size_j), real64) < actual_min_overlap) cycle

                call find_root_helper(tmp_parent, n_seeds, i, root_i)
                call find_root_helper(tmp_parent, n_seeds, j, root_j)

                ! Attach to the smaller index so each component keeps its lowest column.
                if (root_i < root_j) then
                    tmp_parent(root_j) = root_i
                else if (root_j < root_i) then
                    tmp_parent(root_i) = root_j
                end if
            end do
        end do

        merged_masks = .false.
        tmp_active_flag = .false.

        do i = 1, n_seeds
            call find_root_helper(tmp_parent, n_seeds, i, root_i)
            merged_masks(:, root_i) = merged_masks(:, root_i) .or. raw_masks(:, i)
            tmp_active_flag(root_i) = .true.
        end do

        ! Compact merged ensembles to the front, omitting singletons as unassigned/noise.
        n_ensembles = 0_int32
        do i = 1, n_seeds
            if (tmp_active_flag(i) .and. count(merged_masks(:, i)) > 1_int32) then
                n_ensembles = n_ensembles + 1_int32
                if (n_ensembles /= i) then
                    merged_masks(:, n_ensembles) = merged_masks(:, i)
                end if
            end if
        end do

        ! Clear post-compaction columns to prevent stale or duplicate ensembles from being read beyond n_ensembles.
        merged_masks(:, n_ensembles + 1_int32:n_seeds) = .false.

    end subroutine merge_ensembles_impl

    !> summary: Grow one ensemble by one density-compatible surface step.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine grow_ensemble_impl(vectors, n_dimensions, n_vectors, ensemble_mask, &
                                         r, dimension_order, kd_indices, density_labels, &
                                         alpha_mad, tmp_stack, tmp_vicinity_mask, &
                                         tmp_surface_mask, tmp_perm, tmp_abs_diff, grown_mask, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Membership before the growth step.
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        !! DM_MIN(0.0_real64)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions
        !! DM_MIN(1_int32)
        !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        !! DM_MIN(1_int32)
        !! DM_MAX(n_vectors)
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors, from `calculate_labels_as_density`.
        !! DM_MIN(1.0_real64)
        real(real64), intent(in), optional :: alpha_mad
        !! Multiplier factor for MAD density compatibility threshold
        !! DM_MIN(0.0_real64)
        !! DM_DEFAULT(0.5_real64)
        integer(int32), intent(inout) :: tmp_stack(M_KD_STACK_ENTRY_SIZE, M_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated workspace stack for a single tree traversal
        logical(c_bool), intent(out) :: tmp_vicinity_mask(n_vectors)
        !! Workspace for one member's neighborhood query
        logical(c_bool), intent(out) :: tmp_surface_mask(n_vectors)
        !! Geometric union of accepted member neighborhoods before this growth step
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(out) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array holding the active member densities

        logical(c_bool), intent(out) :: grown_mask(n_vectors)
        !! Membership after the growth step: every input member plus each density-compatible neighbor.
        integer(int32), intent(out) :: ierr
        !! Invalid input when the median absolute deviation of `density_labels` is zero (more than
        !! half the labels are equal), which leaves no density tolerance for growth.

        real(real64) :: actual_alpha_mad, median_ambient, mad_ambient
        integer(int32) :: i_vec

        call set_ok(ierr)
        call compute_ambient_density_stats_impl(density_labels, n_vectors, tmp_perm, tmp_abs_diff, &
                                                median_ambient, mad_ambient)
        if (mad_ambient <= 0.0_real64) then
            call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=8_int32)
            return
        end if
        M_DEFAULT_VAL(alpha_mad, actual_alpha_mad, 0.5_real64)

        ! The surface is the union of every member's neighborhood.
        tmp_surface_mask = .false.
        do i_vec = 1, n_vectors
            if (.not. ensemble_mask(i_vec)) cycle
            call vicinity_vectors_helper(vectors(:, i_vec), vectors, n_dimensions, n_vectors, r, &
                                         dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask)
            tmp_surface_mask = tmp_surface_mask .or. tmp_vicinity_mask
        end do

        ! The trial keeps every member, so it is the grown ensemble; this step commits unconditionally.
        call build_ensemble_trial_helper(ensemble_mask, n_vectors, density_labels, actual_alpha_mad, &
                                         mad_ambient, tmp_surface_mask, tmp_perm, tmp_abs_diff, grown_mask)
    end subroutine grow_ensemble_impl

    !> summary: Update the five ensemble observables for a growth iteration.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine compute_ensemble_observable_impl(ensemble_mask, density_labels, n_vectors, &
                                                       candidate_count, current_iter, observables, n_history, n_observables)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Selection mask for active ensemble members
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors, from `calculate_labels_as_density`.
        !! DM_MIN(1.0_real64)
        integer(int32), intent(in) :: candidate_count
        !! Number of candidates added in current growth step
        !! DM_MIN(0_int32)
        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_history
        !! Number of columns in the observable history. This is the history window: once
        !! `current_iter` exceeds it, the columns shift left and the newest step is stored last.
        !! Acceptance compares two consecutive steps, so at least two columns are required.
        !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_observables
        !! Number of observable rows, which must be exactly five.
        !! DM_MIN(CM_OBSERVABLE_COUNT)
        !! DM_MAX(CM_OBSERVABLE_COUNT)
        real(real64), intent(inout) :: observables(n_observables, n_history)
        !! Observable history, one column per growth iteration within the window

        call compute_ensemble_observable_helper(ensemble_mask, density_labels, n_vectors, &
                                                candidate_count, current_iter, observables)
    end subroutine compute_ensemble_observable_impl

    !> summary: Test whether the change in ensemble observables permits a growth step.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine accept_ensemble_impl(observables, n_history, n_observables, current_iter, alpha_accept, &
                                           is_accepted)

        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_history
        !! Number of columns in the observable history, the same window the history was written with.
        !! Acceptance compares two consecutive steps, so at least two columns are required.
        !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_observables
        !! Number of observable rows, which must be exactly five.
        !! DM_MIN(CM_OBSERVABLE_COUNT)
        !! DM_MAX(CM_OBSERVABLE_COUNT)
        real(real64), intent(in) :: observables(n_observables, n_history)
        !! Observable history written by [[tox_stc(module):compute_ensemble_observable]]
        real(real64), intent(in), optional :: alpha_accept
        !! Log2 fold-change threshold; a step whose density heterogeneity changes by at least this much is rejected
        !! DM_MIN(0.0_real64)
        !! DM_DEFAULT(0.5_real64)
        logical(c_bool), intent(out) :: is_accepted
        !! Output flag indicating step acceptance decision

        call accept_ensemble_helper(observables, current_iter, alpha_accept, is_accepted)
    end subroutine accept_ensemble_impl

end module tox_stc_impl
