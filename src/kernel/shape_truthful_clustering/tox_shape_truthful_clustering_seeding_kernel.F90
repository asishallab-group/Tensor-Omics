#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Seeding
!|
!| Kernels for identifying seed points to grow ensembles from: `density_labels` (an
!| adaptive-bandwidth local density estimate) and `seeds` (the greedy, density-ranked,
!| coverage-based seed selection). See `misc/mod_STC.md`, section "Seeding", for the full
!| algorithm definition.
!|
!| `density_labels` and `seeds` both take an already-built k-d tree (`kd_indices`,
!| `dimension_order`, see [[f42_kd_tree(module):build_kd_index_alloc(subroutine)]]) as input
!| rather than building their own: the same tree is also needed by `calc_ensemble_growth_radius`
!| and `grow_ensemble`, so building it once in a future orchestrator (`ensemble_identification`)
!| and passing it to every consumer avoids redundant O(N log N) rebuilds. `seeds_kernel` calls
!| `density_labels_kernel` and, for each selected seed,
!| [[tox_shape_truthful_clustering_ensemble_growing_kernel(module):calc_ensemble_growth_radius_kernel]]
!| directly (all `pure`, no `ierr`) rather than through their own validated entries, to avoid
!| redundant re-validation. Reusing `calc_ensemble_growth_radius_kernel` for the seed coverage
!| radius -- rather than a second, separately-implemented median-of-k-NN-distances computation
!| -- is why this module now depends on its sibling
!| `tox_shape_truthful_clustering_ensemble_growing_kernel`, the first sibling-to-sibling
!| dependency in this family (previously only the parent module called into siblings).
module tox_shape_truthful_clustering_seeding_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_utils, only: sort_real_heapsort, init_perm, calc_percentile_helper
    use f42_kd_tree, only: kd_knn_query_helper, kd_range_query_mask_helper
    use tox_shape_truthful_clustering_ensemble_growing_kernel, only: calc_ensemble_growth_radius_kernel
    M_IMPLICIT_NONE

#define CM_DENSITY_K_DEFAULT 30_int32
#define CM_BANDWIDTH_PERCENTILE_DEFAULT 68.27_real64
#define CM_EXCLUSION_RADIUS_PERCENTILE_DEFAULT 50.0_real64

    private
    public :: density_labels_kernel
    public :: seeds_kernel

contains

    !> summary: Per-vector local density label, an adaptive-bandwidth kernel density estimate over each vector's own k_density nearest neighbors
    !| AUTHOR_ASIS_HALLAB
    !| For each vector: find its `k_density` nearest neighbors (excluding itself), take the
    !| `bandwidth_percentile` percentile of the distances to them as a per-vector local
    !| bandwidth, then sum a Gaussian kernel over those same distances at that bandwidth,
    !| normalized by `bandwidth**n_dimensions`. Unlike a single dataset-wide radius, this
    !| bandwidth shrinks in dense regions and grows in sparse ones, so the resulting labels
    !| reflect local, not global, density. See `misc/mod_STC.md`, SKG `density_labels`.
    !|
    !| `bandwidth_percentile` is a heuristic knob, not a calibrated standard deviation, and
    !| deliberately documented as one: for a genuine 1-D Gaussian, 68.27% of its mass sits
    !| within one SD, so the 68.27th percentile of *samples from that Gaussian* equals the SD
    !| exactly -- which is where the default comes from -- but our distances are norms in
    !| `n_dimensions` dimensions, not draws from a 1-D Gaussian, and the same "percentile that
    !| equals the SD" shifts with dimension (it is ~39% at 2 dimensions, ~20% at 3, following
    !| a chi distribution with `n_dimensions` degrees of freedom, not a plain half-normal).
    !| Correcting for that would need the *local intrinsic* dimension, not the ambient one --
    !| STC's whole premise is data concentrated near a lower-dimensional manifold, so the
    !| ambient `n_dimensions` is typically the wrong dimension to correct with anyway, and the
    !| local one is not yet known at this point in the pipeline (estimating it is `observable`'s
    !| job, run later, on an actual ensemble -- not something to redo per point just to
    !| calibrate a seeding bandwidth). So this is left an explicit, undisguised heuristic:
    !| `bandwidth_percentile` is exposed for exactly this reason -- to be explored
    !| empirically (see `misc/STC-experiments/README.md`) rather than settled by further
    !| first-principles argument that does not actually resolve for manifold-concentrated data.
    !|
    !| The `bandwidth**n_dimensions` normalization is not optional, independent of which
    !| percentile is chosen: a raw `sum(exp(-d**2/(2*bandwidth**2)))`, without it, is
    !| scale-invariant -- scaling every distance (and therefore the bandwidth) by the same
    !| constant leaves `d/bandwidth`, the only thing that enters the exponent, unchanged, so a
    !| tight cluster and the same cluster stretched out 20x would score identically. The
    !| standard adaptive-KDE normalization (divide by the bandwidth to the power of the
    !| ambient dimension) ties the estimate back to an absolute scale, so a genuinely tighter
    !| neighborhood outscores a genuinely looser one, not just a differently-shaped one.
    pure subroutine density_labels_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                          k_density, bandwidth_percentile, &
                                          tmp_neighbors, tmp_distances, tmp_range_stack, tmp_sort_perm, &
                                          labels)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! DM_MIN(2_int32)
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
        integer(int32), intent(in), optional :: k_density
            !! Neighborhood size the local density estimate is taken over
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors - 1_int32)
            !! DM_DEFAULT(CM_DENSITY_K_DEFAULT)
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth -- a heuristic choice, not a calibrated standard deviation,
            !! see above
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(CM_BANDWIDTH_PERCENTILE_DEFAULT)
        integer(int32), intent(out) :: tmp_neighbors(n_vectors)
            !! Workspace: k-NN query result, indices (sized for the worst case, sliced internally)
        real(real64), intent(out) :: tmp_distances(n_vectors)
            !! Workspace: k-NN query result, distances (sized as `tmp_neighbors`)
        integer(int32), intent(out) :: tmp_range_stack(3, n_vectors)
            !! Workspace: k-d tree traversal stack, see `kd_knn_query`
        integer(int32), intent(out) :: tmp_sort_perm(n_vectors)
            !! Workspace: ascending sort permutation of the k_density distances
        real(real64), intent(out) :: labels(n_vectors)
            !! Per-vector local density label

        integer(int32) :: actual_k_density, k_query, self_pos, i_vec, j
        real(real64)   :: actual_bandwidth_percentile, bandwidth

        M_DEFAULT_VAL(k_density, actual_k_density, CM_DENSITY_K_DEFAULT)
        ! See calc_ensemble_growth_radius_kernel's identical clamp and its own comment: an
        ! *explicit* k_density is already wrapper-validated against DM_MAX(n_vectors - 1); this
        ! guards CM_DENSITY_K_DEFAULT itself, which the wrapper never validates against a
        ! runtime-dependent bound when k_density is omitted (misc/code_gen_footgun.md's third
        ! entry) -- without it, a caller on fewer than 31 points who omits k_density would reach
        ! the k-NN query below asking for more neighbors than tmp_neighbors/tmp_distances hold.
        actual_k_density = min(actual_k_density, n_vectors - 1)
        k_query = actual_k_density + 1

        M_DEFAULT_VAL(bandwidth_percentile, actual_bandwidth_percentile, CM_BANDWIDTH_PERCENTILE_DEFAULT)

        do i_vec = 1, n_vectors
            call kd_knn_query_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                     vectors(:, i_vec), k_query, tmp_range_stack, &
                                     tmp_neighbors(1:k_query), tmp_distances(1:k_query))

            ! Exclude the vector itself, at distance 0, from its own k_query nearest
            ! neighbors -- see calc_ensemble_growth_radius_kernel's identical swap-with-last
            ! approach for why this is safe irrespective of kd_knn_query's result order.
            self_pos = 0
            do j = 1, k_query
                if (tmp_neighbors(j) == i_vec) then
                    self_pos = j
                    exit
                end if
            end do
            if (self_pos > 0 .and. self_pos < k_query) then
                tmp_neighbors(self_pos) = tmp_neighbors(k_query)
                tmp_distances(self_pos) = tmp_distances(k_query)
            end if

            do j = 1, actual_k_density
                tmp_sort_perm(j) = j
            end do
            call sort_real_heapsort(tmp_distances(1:actual_k_density), tmp_sort_perm(1:actual_k_density))

            call calc_percentile_helper(tmp_distances(1:actual_k_density), tmp_sort_perm(1:actual_k_density), &
                                        actual_bandwidth_percentile, bandwidth)
            ! A percentile of the raw distances is 0 only for genuinely coincident points (the
            ! chosen percentile landing exactly on one or more zero-distance duplicates) --
            ! unlike the earlier MAD-based bandwidth, an otherwise-regular neighborhood never
            ! drives this to 0 on its own, so this is a last-resort guard, not a routine one.
            bandwidth = max(bandwidth, epsilon(1.0_real64))

            labels(i_vec) = sum(exp(-tmp_distances(1:actual_k_density)**2/(2.0_real64*bandwidth**2))) &
                            /bandwidth**n_dimensions
        end do

    end subroutine density_labels_kernel

    !> summary: Select seed points via greedy, density-ranked, coverage-based selection
    !| AUTHOR_ASIS_HALLAB
    !| Ranks vectors by density label, descending (see `density_labels`). Starting with the
    !| highest-density unvisited vector, marks it a seed, marks every vector within its own
    !| coverage radius as visited, and continues with the next-highest-density unvisited
    !| vector until none remain -- so only genuinely uncovered regions can seed another
    !| ensemble. The coverage radius is
    !| [[tox_shape_truthful_clustering_ensemble_growing_kernel(module):calc_ensemble_growth_radius_kernel]]'s
    !| own computation, called on the newly-selected seed with `k_density` in place of
    !| `k_min` -- not a separate, dataset-wide radius: a fixed global radius can suppress
    !| seed placement across a region much larger than what that seed's own ensemble will
    !| ever actually grow into, leaving points "covered" by seed-exclusion but never reached
    !| by any grown ensemble, see `misc/STC-experiments/README.md`.
    !|
    !| `exclusion_radius_percentile` exposes that computation's own `radius_percentile`
    !| (default 50.0, the median -- unchanged from this SKG's original behavior) so the
    !| exclusion radius can be tuned independently of the actual growth-phase radius any
    !| other caller of `calc_ensemble_growth_radius` relies on: shrinking it here trades
    !| fewer, larger ensembles for less over-eager seed suppression around curvature extrema
    !| (peaks, troughs, kinks) that a seed's own later growth cannot actually reach, see
    !| `misc/STC-experiments/README.md`.
    pure subroutine seeds_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                 k_density, bandwidth_percentile, exclusion_radius_percentile, &
                                 tmp_neighbors, tmp_distances, tmp_range_stack, tmp_sort_perm, &
                                 tmp_labels, tmp_rank_perm, &
                                 tmp_visited_mask, tmp_newly_covered_mask, &
                                 is_seed_mask)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! DM_MIN(2_int32)
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
        integer(int32), intent(in), optional :: k_density
            !! Neighborhood size for both the density estimate and the coverage radius, see
            !! `density_labels` and `calc_ensemble_growth_radius`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors - 1_int32)
            !! DM_DEFAULT(CM_DENSITY_K_DEFAULT)
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth, see `density_labels`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(CM_BANDWIDTH_PERCENTILE_DEFAULT)
        real(real64), intent(in), optional :: exclusion_radius_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as each seed's
            !! coverage/exclusion radius, see above
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(CM_EXCLUSION_RADIUS_PERCENTILE_DEFAULT)
        integer(int32), intent(out) :: tmp_neighbors(n_vectors)
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        real(real64), intent(out) :: tmp_distances(n_vectors)
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        integer(int32), intent(out) :: tmp_range_stack(3, n_vectors)
            !! Workspace: k-d tree traversal stack, reused across `density_labels`,
            !! `calc_ensemble_growth_radius`, and the greedy coverage loop below
        integer(int32), intent(out) :: tmp_sort_perm(n_vectors)
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        real(real64), intent(out) :: tmp_labels(n_vectors)
            !! Workspace: per-vector density labels, see `density_labels`
        integer(int32), intent(out) :: tmp_rank_perm(n_vectors)
            !! Workspace: density-descending sort permutation
        logical, intent(out) :: tmp_visited_mask(n_vectors)
            !! Workspace: coverage tracker across the greedy loop
        logical, intent(out) :: tmp_newly_covered_mask(n_vectors)
            !! Workspace: per-candidate range-query result
        logical, intent(out) :: is_seed_mask(n_vectors)
            !! .true. for points selected as seeds

        real(real64)   :: coverage_radius, actual_exclusion_radius_percentile
        integer(int32) :: i, rank, candidate, swap_tmp

        M_DEFAULT_VAL(exclusion_radius_percentile, actual_exclusion_radius_percentile, CM_EXCLUSION_RADIUS_PERCENTILE_DEFAULT)

        call density_labels_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                   k_density, bandwidth_percentile, &
                                   tmp_neighbors, tmp_distances, tmp_range_stack, tmp_sort_perm, tmp_labels)

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

            call calc_ensemble_growth_radius_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                                    candidate, k_density, &
                                                    radius_percentile=actual_exclusion_radius_percentile, &
                                                    tmp_neighbors=tmp_neighbors, tmp_distances=tmp_distances, &
                                                    tmp_range_stack=tmp_range_stack, tmp_sort_perm=tmp_sort_perm, &
                                                    growth_radius=coverage_radius)

            call kd_range_query_mask_helper(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                            vectors(:, candidate), coverage_radius, tmp_range_stack, &
                                            tmp_newly_covered_mask)
            tmp_visited_mask = tmp_visited_mask .or. tmp_newly_covered_mask
        end do

    end subroutine seeds_kernel

end module tox_shape_truthful_clustering_seeding_kernel
