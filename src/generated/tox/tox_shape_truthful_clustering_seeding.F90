#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_seeding_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_seeding
    use tox_shape_truthful_clustering_seeding_kernel, only: calculate_density_radius_kernel, density_labels_kernel, seeds_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: calculate_density_radius
    public :: calculate_density_radius_alloc
    public :: density_labels
    public :: density_labels_alloc
    public :: seeds
    public :: seeds_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):calculate_density_radius_kernel]].
    !| Computes the mean vector of `vectors`, the Euclidean distance from every vector to it,
    !| and returns the `mean_to_other_vecs_dist_quant` percentile of those distances -- the
    !| radius later used by `density_labels` to measure local density around each vector.
    subroutine calculate_density_radius(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            mean_to_other_vecs_dist_quant,&
            tmp_mean_vec,&
            tmp_distances,&
            tmp_distances_perm,&
            radius,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-vector distances used as the density radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(real64), dimension(n_dimensions), intent(out) :: tmp_mean_vec
            !! Workspace: the global mean vector
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Workspace: mean-to-vector distances
        integer(int32), dimension(n_vectors), intent(out) :: tmp_distances_perm
            !! Workspace: ascending sort permutation of `tmp_distances`
        real(real64), intent(out) :: radius
            !! Resulting density search radius
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, arg_pos=4_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call calculate_density_radius_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            tmp_mean_vec = tmp_mean_vec,&
            tmp_distances = tmp_distances,&
            tmp_distances_perm = tmp_distances_perm,&
            radius = radius&
        )
    end subroutine calculate_density_radius

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):calculate_density_radius_kernel]].
    !| Computes the mean vector of `vectors`, the Euclidean distance from every vector to it,
    !| and returns the `mean_to_other_vecs_dist_quant` percentile of those distances -- the
    !| radius later used by `density_labels` to measure local density around each vector.
    subroutine calculate_density_radius_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            mean_to_other_vecs_dist_quant,&
            radius,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Percentile (0.0 to 1.0) of mean-to-vector distances used as the density radius
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(real64), intent(out) :: radius
            !! Resulting density search radius
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_mean_vec
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:), allocatable :: tmp_distances_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, arg_pos=4_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_mean_vec(n_dimensions))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_distances_perm(n_vectors))

        call calculate_density_radius_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            tmp_mean_vec = tmp_mean_vec,&
            tmp_distances = tmp_distances,&
            tmp_distances_perm = tmp_distances_perm,&
            radius = radius&
        )
    end subroutine calculate_density_radius_alloc

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):density_labels_kernel]].
    !| $\rho_i = \sum_j \mathbf{1}(d(v_i, v_j) \le radius)$, via a k-d tree range-count query
    !| per vector -- see [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]].
    subroutine density_labels(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            radius,&
            tmp_range_stack,&
            labels,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
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
        real(real64), intent(in) :: radius
            !! Density search radius, see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_range_stack
            !! Workspace: k-d tree traversal stack
        real(real64), dimension(n_vectors), intent(out) :: labels
            !! Per-vector density label
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=6_int32, min=0.0_real64)
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
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            labels = labels&
        )
    end subroutine density_labels

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):density_labels_kernel]].
    !| $\rho_i = \sum_j \mathbf{1}(d(v_i, v_j) \le radius)$, via a k-d tree range-count query
    !| per vector -- see [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]].
    subroutine density_labels_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            radius,&
            labels,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
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
        real(real64), intent(in) :: radius
            !! Density search radius, see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_vectors), intent(out) :: labels
            !! Per-vector density label
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=6_int32, min=0.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_range_stack(3, n_vectors))

        call density_labels_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            labels = labels&
        )
    end subroutine density_labels_alloc

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):seeds_kernel]].
    !| Ranks vectors by density label, descending (see `density_labels`). Starting with the
    !| highest-density unvisited vector, marks it a seed, marks every vector within the
    !| density radius of it as visited, and continues with the next-highest-density
    !| unvisited vector until none remain -- so only genuinely uncovered regions can seed
    !| another ensemble.
    subroutine seeds(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
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
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0), see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        real(real64), dimension(n_dimensions), intent(out) :: tmp_mean_vec
            !! Workspace, see `calculate_density_radius`
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Workspace, see `calculate_density_radius`
        integer(int32), dimension(n_vectors), intent(out) :: tmp_distances_perm
            !! Workspace, see `calculate_density_radius`
        real(real64), dimension(n_vectors), intent(out) :: tmp_labels
            !! Workspace: per-vector density labels, see `density_labels`
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_range_stack
            !! greedy coverage loop below
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
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, arg_pos=6_int32, min=0.0_real64, max=1.0_real64)
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
            mean_to_other_vecs_dist_quant = mean_to_other_vecs_dist_quant,&
            tmp_mean_vec = tmp_mean_vec,&
            tmp_distances = tmp_distances,&
            tmp_distances_perm = tmp_distances_perm,&
            tmp_labels = tmp_labels,&
            tmp_range_stack = tmp_range_stack,&
            tmp_rank_perm = tmp_rank_perm,&
            tmp_visited_mask = tmp_visited_mask,&
            tmp_newly_covered_mask = tmp_newly_covered_mask,&
            is_seed_mask = is_seed_mask&
        )
    end subroutine seeds

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_seeding_kernel(module):seeds_kernel]].
    !| Ranks vectors by density label, descending (see `density_labels`). Starting with the
    !| highest-density unvisited vector, marks it a seed, marks every vector within the
    !| density radius of it as visited, and continues with the next-highest-density
    !| unvisited vector until none remain -- so only genuinely uncovered regions can seed
    !| another ensemble.
    subroutine seeds_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            mean_to_other_vecs_dist_quant,&
            is_seed_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
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
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
            !! Density-radius percentile (0.0 to 1.0), see `calculate_density_radius`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.15_real64`.
        logical, dimension(n_vectors), intent(out) :: is_seed_mask
            !! .true. for points selected as seeds
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_mean_vec
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:), allocatable :: tmp_distances_perm
        real(real64), dimension(:), allocatable :: tmp_labels
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack
        integer(int32), dimension(:), allocatable :: tmp_rank_perm
        logical, dimension(:), allocatable :: tmp_visited_mask
        logical, dimension(:), allocatable :: tmp_newly_covered_mask

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, arg_pos=6_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_mean_vec(n_dimensions))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_distances_perm(n_vectors))
        M_ALLOCATE(tmp_labels(n_vectors))
        M_ALLOCATE(tmp_range_stack(3, n_vectors))
        M_ALLOCATE(tmp_rank_perm(n_vectors))
        M_ALLOCATE(tmp_visited_mask(n_vectors))
        M_ALLOCATE(tmp_newly_covered_mask(n_vectors))

        call seeds_kernel(&
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
            tmp_visited_mask = tmp_visited_mask,&
            tmp_newly_covered_mask = tmp_newly_covered_mask,&
            is_seed_mask = is_seed_mask&
        )
    end subroutine seeds_alloc

end module tox_shape_truthful_clustering_seeding
