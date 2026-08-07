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

    public :: calculate_density_radius_expert_c
    public :: calculate_density_radius_c
    public :: density_labels_expert_c
    public :: density_labels_c
    public :: seeds_expert_c
    public :: seeds_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):calculate_density_radius(subroutine)]]
    !| Computes the mean vector of `vectors`, the Euclidean distance from every vector to it,
    !| and returns the `mean_to_other_vecs_dist_quant` percentile of those distances -- the
    !| radius later used by `density_labels` to measure local density around each vector.
    subroutine calculate_density_radius_expert_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            mean_to_other_vecs_dist_quant,&
            tmp_mean_vec,&
            tmp_distances,&
            tmp_distances_perm,&
            radius,&
            ierr&
        ) bind(C, name="calculate_density_radius_expert_c")
        use tox_shape_truthful_clustering_seeding, only: calculate_density_radius

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        real(c_double), intent(in), target :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-vector distances used as the density radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(c_double), dimension(n_dimensions), intent(out), target :: tmp_mean_vec
            !! Workspace: the global mean vector
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Workspace: mean-to-vector distances
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_distances_perm
            !! Workspace: ascending sort permutation of `tmp_distances`
        real(c_double), intent(out), target :: radius
            !! Resulting density search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(mean_to_other_vecs_dist_quant)
        M_CHECK_NON_NULL(radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_vec, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_distances_perm, n_vectors)

        call calculate_density_radius(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            tmp_mean_vec = tmp_mean_vec,&
            tmp_distances = tmp_distances,&
            tmp_distances_perm = tmp_distances_perm,&
            radius = radius,&
            ierr = ierr&
        )
    end subroutine calculate_density_radius_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):calculate_density_radius_alloc(subroutine)]]
    !| Computes the mean vector of `vectors`, the Euclidean distance from every vector to it,
    !| and returns the `mean_to_other_vecs_dist_quant` percentile of those distances -- the
    !| radius later used by `density_labels` to measure local density around each vector.
    subroutine calculate_density_radius_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            mean_to_other_vecs_dist_quant,&
            radius,&
            ierr&
        ) bind(C, name="calculate_density_radius_c")
        use tox_shape_truthful_clustering_seeding, only: calculate_density_radius_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        real(c_double), intent(in), target :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-vector distances used as the density radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(c_double), intent(out), target :: radius
            !! Resulting density search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(mean_to_other_vecs_dist_quant)
        M_CHECK_NON_NULL(radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)

        call calculate_density_radius_alloc(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            radius = radius,&
            ierr = ierr&
        )
    end subroutine calculate_density_radius_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):density_labels(subroutine)]]
    !| $\rho_i = \sum_j \mathbf{1}(d(v_i, v_j) \le radius)$, via a k-d tree range-count query
    !| per vector -- see [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]].
    subroutine density_labels_expert_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            radius,&
            tmp_range_stack,&
            labels,&
            ierr&
        ) bind(C, name="density_labels_expert_c")
        use tox_shape_truthful_clustering_seeding, only: density_labels

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
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
        real(c_double), intent(in), target :: radius
            !! Density search radius, see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(3, n_vectors), intent(out), target :: tmp_range_stack
            !! Workspace: k-d tree traversal stack
        real(c_double), dimension(n_vectors), intent(out), target :: labels
            !! Per-vector density label
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_vectors)
        M_CHECK_ARRAY_NON_NULL(labels, n_vectors)

        call density_labels(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            labels = labels,&
            ierr = ierr&
        )
    end subroutine density_labels_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):density_labels_alloc(subroutine)]]
    !| $\rho_i = \sum_j \mathbf{1}(d(v_i, v_j) \le radius)$, via a k-d tree range-count query
    !| per vector -- see [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]].
    subroutine density_labels_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            radius,&
            labels,&
            ierr&
        ) bind(C, name="density_labels_c")
        use tox_shape_truthful_clustering_seeding, only: density_labels_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
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
        real(c_double), intent(in), target :: radius
            !! Density search radius, see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_vectors), intent(out), target :: labels
            !! Per-vector density label
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(radius)
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
            radius = radius,&
            labels = labels,&
            ierr = ierr&
        )
    end subroutine density_labels_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_seeding(module):seeds(subroutine)]]
    !| Ranks vectors by density label, descending (see `density_labels`). Starting with the
    !| highest-density unvisited vector, marks it a seed, marks every vector within the
    !| density radius of it as visited, and continues with the next-highest-density
    !| unvisited vector until none remain -- so only genuinely uncovered regions can seed
    !| another ensemble.
    subroutine seeds_expert_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            mean_to_other_vecs_dist_quant,&
            tmp_mean_vec,&
            tmp_distances,&
            tmp_distances_perm,&
            tmp_labels,&
            tmp_range_stack,&
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
            !! Number of input vectors N
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
        real(c_double), intent(in), target :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0), see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(c_double), dimension(n_dimensions), intent(out), target :: tmp_mean_vec
            !! Workspace, see `calculate_density_radius`
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Workspace, see `calculate_density_radius`
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_distances_perm
            !! Workspace, see `calculate_density_radius`
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_labels
            !! Workspace: per-vector density labels, see `density_labels`
        integer(c_int), dimension(3, n_vectors), intent(out), target :: tmp_range_stack
            !! greedy coverage loop below
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
        M_CHECK_NON_NULL(mean_to_other_vecs_dist_quant)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_vec, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_distances_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_vectors)
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
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            tmp_mean_vec = tmp_mean_vec,&
            tmp_distances = tmp_distances,&
            tmp_distances_perm = tmp_distances_perm,&
            tmp_labels = tmp_labels,&
            tmp_range_stack = tmp_range_stack,&
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
    !| highest-density unvisited vector, marks it a seed, marks every vector within the
    !| density radius of it as visited, and continues with the next-highest-density
    !| unvisited vector until none remain -- so only genuinely uncovered regions can seed
    !| another ensemble.
    subroutine seeds_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            mean_to_other_vecs_dist_quant,&
            is_seed_mask,&
            ierr&
        ) bind(C, name="seeds_c")
        use tox_shape_truthful_clustering_seeding, only: seeds_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
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
        real(c_double), intent(in), target :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0), see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        logical(c_bool), dimension(n_vectors), intent(out), target :: is_seed_mask
            !! .true. for points selected as seeds
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_vectors) :: is_seed_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(mean_to_other_vecs_dist_quant)
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
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            is_seed_mask = is_seed_mask_f,&
            ierr = ierr&
        )

        is_seed_mask = is_seed_mask_f
    end subroutine seeds_c

end module tox_shape_truthful_clustering_seeding_c
#endif
