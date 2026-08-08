#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_ensemble_growing_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_ensemble_growing
    use tox_shape_truthful_clustering_ensemble_growing_kernel, only: calc_ensemble_growth_radius_kernel, grow_ensemble_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    use tox_errors, only: validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: calc_ensemble_growth_radius
    public :: calc_ensemble_growth_radius_alloc
    public :: grow_ensemble
    public :: grow_ensemble_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_ensemble_growing_kernel(module):calc_ensemble_growth_radius_kernel]].
    !| Matches LoManLe's `local_scale_i` exactly at the default `radius_percentile=50.0` (the
    !| median) -- a per-seed, locally adaptive radius rather than a single dataset-wide one
    !| (see `misc/STC_for_LoManLe.md` section 2.2). Work arrays are sized for the worst case
    !| (`k_min = n_vectors - 1`) and sliced internally, since `k_min`'s resolved value is only
    !| known once its default (if any) has been applied.
    !|
    !| `radius_percentile` generalizes what used to be a hardcoded median: `seeds_kernel`
    !| reuses this same SKG for its seed-exclusion radius (see `misc/mod_STC.md`, SKG
    !| `seeds`), and a fixed dataset-wide-in-spirit median there was observed to suppress
    !| entire uncovered regions (e.g. curvature extrema on a wavy manifold) whose own growth
    !| never actually reaches that far -- exposing the percentile lets that specific call site
    !| shrink its exclusion radius independently of this kernel's own default, without
    !| touching the actual growth-radius computation any other caller relies on.
    subroutine calc_ensemble_growth_radius(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            seed_index,&
            k_min,&
            radius_percentile,&
            tmp_neighbors,&
            tmp_distances,&
            tmp_range_stack,&
            tmp_sort_perm,&
            growth_radius,&
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
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the growth radius for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size the median distance is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in), optional :: radius_percentile
            !! Percentile (0 to 100) of the k_min neighbor distances reported as the growth
            !! radius -- 50.0 (the default) is the median, matching this SKG's original,
            !! non-parameterized behavior
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `50.0_real64`.
        integer(int32), dimension(n_vectors), intent(out) :: tmp_neighbors
            !! Workspace: k-NN query result, indices (sized for the worst case, sliced internally)
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Workspace: k-NN query result, distances (sized as `tmp_neighbors`)
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_range_stack
            !! Workspace: k-d tree traversal stack, see `kd_knn_query`
        integer(int32), dimension(n_vectors), intent(out) :: tmp_sort_perm
            !! Workspace: sort-permutation scratch (sized as `tmp_neighbors`)
        real(real64), intent(out) :: growth_radius
            !! radius_percentile-th percentile of the distances among the seed's own k_min
            !! nearest neighbors
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(seed_index, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors)
        call validate_in_range_int(k_min, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(radius_percentile, ierr, arg_pos=8_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call calc_ensemble_growth_radius_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            seed_index = seed_index,&
            k_min = k_min,&
            radius_percentile = radius_percentile,&
            tmp_neighbors = tmp_neighbors,&
            tmp_distances = tmp_distances,&
            tmp_range_stack = tmp_range_stack,&
            tmp_sort_perm = tmp_sort_perm,&
            growth_radius = growth_radius&
        )
    end subroutine calc_ensemble_growth_radius

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_ensemble_growing_kernel(module):calc_ensemble_growth_radius_kernel]].
    !| Matches LoManLe's `local_scale_i` exactly at the default `radius_percentile=50.0` (the
    !| median) -- a per-seed, locally adaptive radius rather than a single dataset-wide one
    !| (see `misc/STC_for_LoManLe.md` section 2.2). Work arrays are sized for the worst case
    !| (`k_min = n_vectors - 1`) and sliced internally, since `k_min`'s resolved value is only
    !| known once its default (if any) has been applied.
    !|
    !| `radius_percentile` generalizes what used to be a hardcoded median: `seeds_kernel`
    !| reuses this same SKG for its seed-exclusion radius (see `misc/mod_STC.md`, SKG
    !| `seeds`), and a fixed dataset-wide-in-spirit median there was observed to suppress
    !| entire uncovered regions (e.g. curvature extrema on a wavy manifold) whose own growth
    !| never actually reaches that far -- exposing the percentile lets that specific call site
    !| shrink its exclusion radius independently of this kernel's own default, without
    !| touching the actual growth-radius computation any other caller relies on.
    subroutine calc_ensemble_growth_radius_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            seed_index,&
            k_min,&
            radius_percentile,&
            growth_radius,&
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
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the growth radius for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size the median distance is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(real64), intent(in), optional :: radius_percentile
            !! Percentile (0 to 100) of the k_min neighbor distances reported as the growth
            !! radius -- 50.0 (the default) is the median, matching this SKG's original,
            !! non-parameterized behavior
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `50.0_real64`.
        real(real64), intent(out) :: growth_radius
            !! radius_percentile-th percentile of the distances among the seed's own k_min
            !! nearest neighbors
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
        call validate_in_range_int(seed_index, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors)
        call validate_in_range_int(k_min, ierr, arg_pos=7_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(radius_percentile, ierr, arg_pos=8_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_neighbors(n_vectors))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_range_stack(3, n_vectors))
        M_ALLOCATE(tmp_sort_perm(n_vectors))

        call calc_ensemble_growth_radius_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            seed_index = seed_index,&
            k_min = k_min,&
            radius_percentile = radius_percentile,&
            tmp_neighbors = tmp_neighbors,&
            tmp_distances = tmp_distances,&
            tmp_range_stack = tmp_range_stack,&
            tmp_sort_perm = tmp_sort_perm,&
            growth_radius = growth_radius&
        )
    end subroutine calc_ensemble_growth_radius_alloc

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_ensemble_growing_kernel(module):grow_ensemble_kernel]].
    !| $\mathcal{C}_{\mathcal{E}_{t+1}} = \{x_k \mid \|x_k-x_i\|\le r_{\mathcal{E}} \;\exists\,
    !| x_i\in\mathcal{E}\}$. Deliberately a plain sequential loop, not parallelized: this runs
    !| once per growth iteration per ensemble, and outer-level parallelism across
    !| ensembles/seeds (in `ensemble_identification`) is the right place for that -- a single
    !| ensemble's member count is typically small, especially early in growth. An all-.false.
    !| `is_member_mask` (an empty ensemble) is a well-defined degenerate case: there is nothing
    !| to grow from, so the result is all-.false. too, not an error.
    subroutine grow_ensemble(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            is_member_mask,&
            growth_radius,&
            tmp_range_stack,&
            tmp_member_mask_buf,&
            is_member_mask_next,&
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
        logical, dimension(n_vectors), intent(in) :: is_member_mask
            !! Current ensemble membership
        real(real64), intent(in) :: growth_radius
            !! This ensemble's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_range_stack
            !! Workspace: k-d tree traversal stack
        logical, dimension(n_vectors), intent(out) :: tmp_member_mask_buf
            !! Workspace: per-member range-query result
        logical, dimension(n_vectors), intent(out) :: is_member_mask_next
            !! Grown ensemble membership (superset of `is_member_mask`)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(growth_radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call grow_ensemble_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            is_member_mask = is_member_mask,&
            growth_radius = growth_radius,&
            tmp_range_stack = tmp_range_stack,&
            tmp_member_mask_buf = tmp_member_mask_buf,&
            is_member_mask_next = is_member_mask_next&
        )
    end subroutine grow_ensemble

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_ensemble_growing_kernel(module):grow_ensemble_kernel]].
    !| $\mathcal{C}_{\mathcal{E}_{t+1}} = \{x_k \mid \|x_k-x_i\|\le r_{\mathcal{E}} \;\exists\,
    !| x_i\in\mathcal{E}\}$. Deliberately a plain sequential loop, not parallelized: this runs
    !| once per growth iteration per ensemble, and outer-level parallelism across
    !| ensembles/seeds (in `ensemble_identification`) is the right place for that -- a single
    !| ensemble's member count is typically small, especially early in growth. An all-.false.
    !| `is_member_mask` (an empty ensemble) is a well-defined degenerate case: there is nothing
    !| to grow from, so the result is all-.false. too, not an error.
    subroutine grow_ensemble_alloc(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            is_member_mask,&
            growth_radius,&
            is_member_mask_next,&
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
        logical, dimension(n_vectors), intent(in) :: is_member_mask
            !! Current ensemble membership
        real(real64), intent(in) :: growth_radius
            !! This ensemble's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `0.0_real64`.
        logical, dimension(n_vectors), intent(out) :: is_member_mask_next
            !! Grown ensemble membership (superset of `is_member_mask`)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack
        logical, dimension(:), allocatable :: tmp_member_mask_buf

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=3_int32)
        call validate_in_range_real(growth_radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_range_stack(3, n_vectors))
        M_ALLOCATE(tmp_member_mask_buf(n_vectors))

        call grow_ensemble_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            is_member_mask = is_member_mask,&
            growth_radius = growth_radius,&
            tmp_range_stack = tmp_range_stack,&
            tmp_member_mask_buf = tmp_member_mask_buf,&
            is_member_mask_next = is_member_mask_next&
        )
    end subroutine grow_ensemble_alloc

end module tox_shape_truthful_clustering_ensemble_growing
