#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Seeding
!|
!| Kernels for identifying seed points to grow ensembles from: `calculate_density_radius`
!| (the local-density search radius), `density_labels` (per-vector density under that
!| radius), and `seeds` (the greedy, density-ranked, coverage-based seed selection). See
!| `misc/mod_STC.md`, section "Seeding", for the full algorithm definition.
!|
!| `density_labels` and `seeds` both take an already-built k-d tree (`kd_indices`,
!| `dimension_order`, see [[f42_kd_tree(module):build_kd_index_alloc(subroutine)]]) as input
!| rather than building their own: the same tree is also needed by `calc_ensemble_growth_radius`
!| and `grow_ensemble`, so building it once in a future orchestrator (`ensemble_identification`)
!| and passing it to every consumer avoids redundant O(N log N) rebuilds. `seeds_kernel` calls
!| `calculate_density_radius_kernel` and `density_labels_kernel` directly (both `pure`, no
!| `ierr`) rather than through their own validated entries, to avoid redundant re-validation.
module tox_shape_truthful_clustering_seeding_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_utils, only: sort_real_heapsort, init_perm, calc_percentile_helper
    use f42_kd_tree, only: kd_range_query_count_helper, kd_range_query_mask_helper
    M_IMPLICIT_NONE

#define CM_DENSITY_RADIUS_QUANTILE_DEFAULT 0.15_real64

    private
    public :: calculate_density_radius_kernel
    public :: density_labels_kernel
    public :: seeds_kernel

contains

    !> summary: Local-density search radius, a percentile of mean-to-vector distances
    !| AUTHOR_ASIS_HALLAB
    !| Computes the mean vector of `vectors`, the Euclidean distance from every vector to it,
    !| and returns the `mean_to_other_vecs_dist_quant` percentile of those distances -- the
    !| radius later used by `density_labels` to measure local density around each vector.
    pure subroutine calculate_density_radius_kernel(vectors, n_dimensions, n_vectors, &
                                                     mean_to_other_vecs_dist_quant, &
                                                     tmp_mean_vec, tmp_distances, tmp_distances_perm, &
                                                     radius)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-vector distances used as the density radius
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_DENSITY_RADIUS_QUANTILE_DEFAULT)
        real(real64), intent(out) :: tmp_mean_vec(n_dimensions)
            !! Workspace: the global mean vector
        real(real64), intent(out) :: tmp_distances(n_vectors)
            !! Workspace: mean-to-vector distances
        integer(int32), intent(out) :: tmp_distances_perm(n_vectors)
            !! Workspace: ascending sort permutation of `tmp_distances`
        real(real64), intent(out) :: radius
            !! Resulting density search radius

        real(real64)   :: actual_quant
        integer(int32) :: i_vec

        M_DEFAULT_VAL(mean_to_other_vecs_dist_quant, actual_quant, CM_DENSITY_RADIUS_QUANTILE_DEFAULT)

        tmp_mean_vec = sum(vectors, dim=2)/real(n_vectors, real64)

        do i_vec = 1, n_vectors
            tmp_distances(i_vec) = sqrt(sum((vectors(:, i_vec) - tmp_mean_vec)**2))
        end do

        call init_perm(tmp_distances_perm)
        call sort_real_heapsort(tmp_distances, tmp_distances_perm)

        call calc_percentile_helper(tmp_distances, tmp_distances_perm, actual_quant*100.0_real64, radius)

    end subroutine calculate_density_radius_kernel

    !> summary: Per-vector density label, the count of vectors (including itself) within `radius`
    !| AUTHOR_ASIS_HALLAB
    !| $\rho_i = \sum_j \mathbf{1}(d(v_i, v_j) \le radius)$, via a k-d tree range-count query
    !| per vector -- see [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]].
    pure subroutine density_labels_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                          radius, tmp_range_stack, labels)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree index over `vectors`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! Dimension order used to build `kd_indices`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in) :: radius
            !! Density search radius, see `calculate_density_radius`
            !! DM_MIN(0.0_real64)
        integer(int32), intent(out) :: tmp_range_stack(3, n_vectors)
            !! Workspace: k-d tree traversal stack
        real(real64), intent(out) :: labels(n_vectors)
            !! Per-vector density label

        integer(int32) :: i_vec, neighbor_count

        do i_vec = 1, n_vectors
            call kd_range_query_count_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                             vectors(:, i_vec), radius, tmp_range_stack, neighbor_count)
            labels(i_vec) = real(neighbor_count, real64)
        end do

    end subroutine density_labels_kernel

    !> summary: Select seed points via greedy, density-ranked, coverage-based selection
    !| AUTHOR_ASIS_HALLAB
    !| Ranks vectors by density label, descending (see `density_labels`). Starting with the
    !| highest-density unvisited vector, marks it a seed, marks every vector within the
    !| density radius of it as visited, and continues with the next-highest-density
    !| unvisited vector until none remain -- so only genuinely uncovered regions can seed
    !| another ensemble.
    pure subroutine seeds_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                 mean_to_other_vecs_dist_quant, &
                                 tmp_mean_vec, tmp_distances, tmp_distances_perm, &
                                 tmp_labels, tmp_range_stack, tmp_rank_perm, &
                                 tmp_visited_mask, tmp_newly_covered_mask, &
                                 is_seed_mask)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: kd_indices(n_vectors)
            !! Pre-built k-d tree index over `vectors`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! Dimension order used to build `kd_indices`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0), see `calculate_density_radius`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_DENSITY_RADIUS_QUANTILE_DEFAULT)
        real(real64), intent(out) :: tmp_mean_vec(n_dimensions)
            !! Workspace, see `calculate_density_radius`
        real(real64), intent(out) :: tmp_distances(n_vectors)
            !! Workspace, see `calculate_density_radius`
        integer(int32), intent(out) :: tmp_distances_perm(n_vectors)
            !! Workspace, see `calculate_density_radius`
        real(real64), intent(out) :: tmp_labels(n_vectors)
            !! Workspace: per-vector density labels, see `density_labels`
        integer(int32), intent(out) :: tmp_range_stack(3, n_vectors)
            !! Workspace: k-d tree traversal stack, reused across `density_labels` and the
            !! greedy coverage loop below
        integer(int32), intent(out) :: tmp_rank_perm(n_vectors)
            !! Workspace: density-descending sort permutation
        logical, intent(out) :: tmp_visited_mask(n_vectors)
            !! Workspace: coverage tracker across the greedy loop
        logical, intent(out) :: tmp_newly_covered_mask(n_vectors)
            !! Workspace: per-candidate range-query result
        logical, intent(out) :: is_seed_mask(n_vectors)
            !! .true. for points selected as seeds

        real(real64)   :: radius
        integer(int32) :: i, rank, candidate, swap_tmp

        call calculate_density_radius_kernel(vectors, n_dimensions, n_vectors, mean_to_other_vecs_dist_quant, &
                                             tmp_mean_vec, tmp_distances, tmp_distances_perm, radius)

        call density_labels_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, radius, &
                                   tmp_range_stack, tmp_labels)

        ! Rank vectors by density, descending: heapsort gives ascending order, so reverse the
        ! resulting permutation.
        call init_perm(tmp_rank_perm)
        call sort_real_heapsort(tmp_labels, tmp_rank_perm)
        do i = 1, n_vectors/2
            swap_tmp = tmp_rank_perm(i)
            tmp_rank_perm(i) = tmp_rank_perm(n_vectors - i + 1)
            tmp_rank_perm(n_vectors - i + 1) = swap_tmp
        end do

        is_seed_mask = .false.
        tmp_visited_mask = .false.

        do rank = 1, n_vectors
            candidate = tmp_rank_perm(rank)
            if (tmp_visited_mask(candidate)) cycle
            is_seed_mask(candidate) = .true.
            call kd_range_query_mask_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                            vectors(:, candidate), radius, tmp_range_stack, tmp_newly_covered_mask)
            tmp_visited_mask = tmp_visited_mask .or. tmp_newly_covered_mask
        end do

    end subroutine seeds_kernel

end module tox_shape_truthful_clustering_seeding_kernel
