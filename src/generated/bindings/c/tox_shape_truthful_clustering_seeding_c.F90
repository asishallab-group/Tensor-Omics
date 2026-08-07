#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_seeding(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_seeding_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: density_labels_expert_c
    public :: density_labels_c
    public :: seeds_expert_c
    public :: seeds_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):density_labels(subroutine)]]
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
    subroutine density_labels_expert_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            tmp_neighbors,&
            tmp_distances,&
            tmp_range_stack,&
            tmp_sort_perm,&
            labels,&
            ierr&
        ) bind(C, name="density_labels_expert_c")
        use tox_shape_truthful_clustering_seeding, only: density_labels

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: k_density
            !! Neighborhood size the local density estimate is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth -- a heuristic choice, not a calibrated standard deviation,
            !! see above
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_neighbors
            !! Workspace: k-NN query result, indices (sized for the worst case, sliced internally)
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Workspace: k-NN query result, distances (sized as `tmp_neighbors`)
        integer(c_int), dimension(3, n_vectors), intent(out), target :: tmp_range_stack
            !! Workspace: k-d tree traversal stack, see `kd_knn_query`
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_sort_perm
            !! Workspace: ascending sort permutation of the k_density distances
        real(c_double), dimension(n_vectors), intent(out), target :: labels
            !! Per-vector local density label
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(k_density)
        M_CHECK_NON_NULL(bandwidth_percentile)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbors, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_sort_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(labels, n_vectors)

        call density_labels(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            k_density = k_density,&
            bandwidth_percentile = bandwidth_percentile,&
            tmp_neighbors = tmp_neighbors,&
            tmp_distances = tmp_distances,&
            tmp_range_stack = tmp_range_stack,&
            tmp_sort_perm = tmp_sort_perm,&
            labels = labels,&
            ierr = ierr&
        )
    end subroutine density_labels_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):density_labels_alloc(subroutine)]]
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
    subroutine density_labels_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            labels,&
            ierr&
        ) bind(C, name="density_labels_c")
        use tox_shape_truthful_clustering_seeding, only: density_labels_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: k_density
            !! Neighborhood size the local density estimate is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth -- a heuristic choice, not a calibrated standard deviation,
            !! see above
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        real(c_double), dimension(n_vectors), intent(out), target :: labels
            !! Per-vector local density label
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(k_density)
        M_CHECK_NON_NULL(bandwidth_percentile)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(labels, n_vectors)

        call density_labels_alloc(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            k_density = k_density,&
            bandwidth_percentile = bandwidth_percentile,&
            labels = labels,&
            ierr = ierr&
        )
    end subroutine density_labels_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):seeds(subroutine)]]
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
    subroutine seeds_expert_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            tmp_neighbors,&
            tmp_distances,&
            tmp_range_stack,&
            tmp_sort_perm,&
            tmp_labels,&
            tmp_rank_perm,&
            tmp_visited_mask,&
            tmp_newly_covered_mask,&
            is_seed_mask,&
            ierr&
        ) bind(C, name="seeds_expert_c")
        use tox_shape_truthful_clustering_seeding, only: seeds

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: k_density
            !! Neighborhood size for both the density estimate and the coverage radius, see
            !! `density_labels` and `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth, see `density_labels`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_neighbors
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        integer(c_int), dimension(3, n_vectors), intent(out), target :: tmp_range_stack
            !! `calc_ensemble_growth_radius`, and the greedy coverage loop below
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_sort_perm
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_labels
            !! Workspace: per-vector density labels, see `density_labels`
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_rank_perm
            !! Workspace: density-descending sort permutation
        logical(c_bool), dimension(n_vectors), intent(out), target :: tmp_visited_mask
            !! Workspace: coverage tracker across the greedy loop
        logical(c_bool), dimension(n_vectors), intent(out), target :: tmp_newly_covered_mask
            !! Workspace: per-candidate range-query result
        logical(c_bool), dimension(n_vectors), intent(out), target :: is_seed_mask
            !! .true. for points selected as seeds
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_vectors) :: tmp_visited_mask_f
        logical, dimension(n_vectors) :: tmp_newly_covered_mask_f
        logical, dimension(n_vectors) :: is_seed_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(k_density)
        M_CHECK_NON_NULL(bandwidth_percentile)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbors, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_sort_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_rank_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_visited_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_newly_covered_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(is_seed_mask, n_vectors)

        call seeds(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            k_density = k_density,&
            bandwidth_percentile = bandwidth_percentile,&
            tmp_neighbors = tmp_neighbors,&
            tmp_distances = tmp_distances,&
            tmp_range_stack = tmp_range_stack,&
            tmp_sort_perm = tmp_sort_perm,&
            tmp_labels = tmp_labels,&
            tmp_rank_perm = tmp_rank_perm,&
            tmp_visited_mask = tmp_visited_mask_f,&
            tmp_newly_covered_mask = tmp_newly_covered_mask_f,&
            is_seed_mask = is_seed_mask_f,&
            ierr = ierr&
        )

        tmp_visited_mask = tmp_visited_mask_f
        tmp_newly_covered_mask = tmp_newly_covered_mask_f
        is_seed_mask = is_seed_mask_f
    end subroutine seeds_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):seeds_alloc(subroutine)]]
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
    subroutine seeds_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            is_seed_mask,&
            ierr&
        ) bind(C, name="seeds_c")
        use tox_shape_truthful_clustering_seeding, only: seeds_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        integer(c_int), dimension(n_vectors), intent(in), target :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(in), target :: k_density
            !! Neighborhood size for both the density estimate and the coverage radius, see
            !! `density_labels` and `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth, see `density_labels`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        logical(c_bool), dimension(n_vectors), intent(out), target :: is_seed_mask
            !! .true. for points selected as seeds
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_vectors) :: is_seed_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(k_density)
        M_CHECK_NON_NULL(bandwidth_percentile)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(is_seed_mask, n_vectors)

        call seeds_alloc(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            k_density = k_density,&
            bandwidth_percentile = bandwidth_percentile,&
            is_seed_mask = is_seed_mask_f,&
            ierr = ierr&
        )

        is_seed_mask = is_seed_mask_f
    end subroutine seeds_c

end module tox_shape_truthful_clustering_seeding_c
#endif
