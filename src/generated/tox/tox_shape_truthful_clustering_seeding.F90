#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_seeding_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_seeding
    use tox_shape_truthful_clustering_seeding_kernel, only: density_labels_kernel, seeds_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    use tox_errors, only: validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: density_labels
    public :: density_labels_alloc
    public :: seeds
    public :: seeds_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):density_labels_kernel]].
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
    subroutine density_labels(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in), optional :: k_density
            !! Neighborhood size the local density estimate is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth -- a heuristic choice, not a calibrated standard deviation,
            !! see above
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        integer(int32), dimension(n_vectors), intent(out) :: tmp_neighbors
            !! Workspace: k-NN query result, indices (sized for the worst case, sliced internally)
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Workspace: k-NN query result, distances (sized as `tmp_neighbors`)
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_range_stack
            !! Workspace: k-d tree traversal stack, see `kd_knn_query`
        integer(int32), dimension(n_vectors), intent(out) :: tmp_sort_perm
            !! Workspace: ascending sort permutation of the k_density distances
        real(real64), dimension(n_vectors), intent(out) :: labels
            !! Per-vector local density label
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_density, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(bandwidth_percentile, ierr, arg_pos=7_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call density_labels_kernel(&
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
            labels = labels&
        )
    end subroutine density_labels

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):density_labels_kernel]].
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
    subroutine density_labels_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            labels,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in), optional :: k_density
            !! Neighborhood size the local density estimate is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth -- a heuristic choice, not a calibrated standard deviation,
            !! see above
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        real(real64), dimension(n_vectors), intent(out) :: labels
            !! Per-vector local density label
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_neighbors
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack
        integer(int32), dimension(:), allocatable :: tmp_sort_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_density, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(bandwidth_percentile, ierr, arg_pos=7_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_neighbors(n_vectors))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_range_stack(3, n_vectors))
        M_ALLOCATE(tmp_sort_perm(n_vectors))

        call density_labels_kernel(&
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
            labels = labels&
        )
    end subroutine density_labels_alloc

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):seeds_kernel]].
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
    subroutine seeds(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in), optional :: k_density
            !! Neighborhood size for both the density estimate and the coverage radius, see
            !! `density_labels` and `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth, see `density_labels`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        integer(int32), dimension(n_vectors), intent(out) :: tmp_neighbors
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_range_stack
            !! `calc_ensemble_growth_radius`, and the greedy coverage loop below
        integer(int32), dimension(n_vectors), intent(out) :: tmp_sort_perm
            !! Workspace, see `density_labels`/`calc_ensemble_growth_radius` (reused across both)
        real(real64), dimension(n_vectors), intent(out) :: tmp_labels
            !! Workspace: per-vector density labels, see `density_labels`
        integer(int32), dimension(n_vectors), intent(out) :: tmp_rank_perm
            !! Workspace: density-descending sort permutation
        logical, dimension(n_vectors), intent(out) :: tmp_visited_mask
            !! Workspace: coverage tracker across the greedy loop
        logical, dimension(n_vectors), intent(out) :: tmp_newly_covered_mask
            !! Workspace: per-candidate range-query result
        logical, dimension(n_vectors), intent(out) :: is_seed_mask
            !! .true. for points selected as seeds
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_density, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(bandwidth_percentile, ierr, arg_pos=7_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call seeds_kernel(&
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
            tmp_visited_mask = tmp_visited_mask,&
            tmp_newly_covered_mask = tmp_newly_covered_mask,&
            is_seed_mask = is_seed_mask&
        )
    end subroutine seeds

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):seeds_kernel]].
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
    subroutine seeds_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            is_seed_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a "nearest neighbor" is undefined for a
            !! single point.
            !! The minimum valid value is `2_int32`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        integer(int32), dimension(n_vectors), intent(in) :: kd_indices
            !! Pre-built k-d tree index over `vectors`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(in), optional :: k_density
            !! Neighborhood size for both the density estimate and the coverage radius, see
            !! `density_labels` and `calc_ensemble_growth_radius`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Percentile (0 to 100) of the k_density neighbor distances used as the local
            !! Gaussian bandwidth, see `density_labels`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `68.27_real64`.
        logical, dimension(n_vectors), intent(out) :: is_seed_mask
            !! .true. for points selected as seeds
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_neighbors
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack
        integer(int32), dimension(:), allocatable :: tmp_sort_perm
        real(real64), dimension(:), allocatable :: tmp_labels
        integer(int32), dimension(:), allocatable :: tmp_rank_perm
        logical, dimension(:), allocatable :: tmp_visited_mask
        logical, dimension(:), allocatable :: tmp_newly_covered_mask

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_density, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(bandwidth_percentile, ierr, arg_pos=7_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_neighbors(n_vectors))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_range_stack(3, n_vectors))
        M_ALLOCATE(tmp_sort_perm(n_vectors))
        M_ALLOCATE(tmp_labels(n_vectors))
        M_ALLOCATE(tmp_rank_perm(n_vectors))
        M_ALLOCATE(tmp_visited_mask(n_vectors))
        M_ALLOCATE(tmp_newly_covered_mask(n_vectors))

        call seeds_kernel(&
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
            tmp_visited_mask = tmp_visited_mask,&
            tmp_newly_covered_mask = tmp_newly_covered_mask,&
            is_seed_mask = is_seed_mask&
        )
    end subroutine seeds_alloc

end module tox_shape_truthful_clustering_seeding
