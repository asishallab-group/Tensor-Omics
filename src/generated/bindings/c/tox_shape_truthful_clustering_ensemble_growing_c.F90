#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_ensemble_growing(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_ensemble_growing_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_ensemble_growth_radius_expert_c
    public :: calc_ensemble_growth_radius_c
    public :: grow_ensemble_expert_c
    public :: grow_ensemble_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_ensemble_growing(module):calc_ensemble_growth_radius(subroutine)]]
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
    subroutine calc_ensemble_growth_radius_expert_c(&
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
        ) bind(C, name="calc_ensemble_growth_radius_expert_c")
        use tox_shape_truthful_clustering_ensemble_growing, only: calc_ensemble_growth_radius

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
        integer(c_int), intent(in), target :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the growth radius for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: k_min
            !! Neighborhood size the median distance is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: radius_percentile
            !! Percentile (0 to 100) of the k_min neighbor distances reported as the growth
            !! radius -- 50.0 (the default) is the median, matching this SKG's original,
            !! non-parameterized behavior
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `50.0_real64`.
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_neighbors
            !! Workspace: k-NN query result, indices (sized for the worst case, sliced internally)
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Workspace: k-NN query result, distances (sized as `tmp_neighbors`)
        integer(c_int), dimension(3, n_vectors), intent(out), target :: tmp_range_stack
            !! Workspace: k-d tree traversal stack, see `kd_knn_query`
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_sort_perm
            !! Workspace: sort-permutation scratch (sized as `tmp_neighbors`)
        real(c_double), intent(out), target :: growth_radius
            !! radius_percentile-th percentile of the distances among the seed's own k_min
            !! nearest neighbors
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(seed_index)
        M_CHECK_NON_NULL(k_min)
        M_CHECK_NON_NULL(radius_percentile)
        M_CHECK_NON_NULL(growth_radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbors, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_sort_perm, n_vectors)

        call calc_ensemble_growth_radius(&
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
            growth_radius = growth_radius,&
            ierr = ierr&
        )
    end subroutine calc_ensemble_growth_radius_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_ensemble_growing(module):calc_ensemble_growth_radius_alloc(subroutine)]]
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
    subroutine calc_ensemble_growth_radius_c(&
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
        ) bind(C, name="calc_ensemble_growth_radius_c")
        use tox_shape_truthful_clustering_ensemble_growing, only: calc_ensemble_growth_radius_alloc

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
        integer(c_int), intent(in), target :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the growth radius for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        integer(c_int), intent(in), target :: k_min
            !! Neighborhood size the median distance is taken over
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
            !! The default value is `30_int32`.
        real(c_double), intent(in), target :: radius_percentile
            !! Percentile (0 to 100) of the k_min neighbor distances reported as the growth
            !! radius -- 50.0 (the default) is the median, matching this SKG's original,
            !! non-parameterized behavior
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `50.0_real64`.
        real(c_double), intent(out), target :: growth_radius
            !! radius_percentile-th percentile of the distances among the seed's own k_min
            !! nearest neighbors
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(seed_index)
        M_CHECK_NON_NULL(k_min)
        M_CHECK_NON_NULL(radius_percentile)
        M_CHECK_NON_NULL(growth_radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)

        call calc_ensemble_growth_radius_alloc(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            seed_index = seed_index,&
            k_min = k_min,&
            radius_percentile = radius_percentile,&
            growth_radius = growth_radius,&
            ierr = ierr&
        )
    end subroutine calc_ensemble_growth_radius_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_ensemble_growing(module):grow_ensemble(subroutine)]]
    !| $\mathcal{C}_{\mathcal{E}_{t+1}} = \{x_k \mid \|x_k-x_i\|\le r_{\mathcal{E}} \;\exists\,
    !| x_i\in\mathcal{E}\}$. Deliberately a plain sequential loop, not parallelized: this runs
    !| once per growth iteration per ensemble, and outer-level parallelism across
    !| ensembles/seeds (in `ensemble_identification`) is the right place for that -- a single
    !| ensemble's member count is typically small, especially early in growth. An all-.false.
    !| `is_member_mask` (an empty ensemble) is a well-defined degenerate case: there is nothing
    !| to grow from, so the result is all-.false. too, not an error.
    subroutine grow_ensemble_expert_c(&
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
        ) bind(C, name="grow_ensemble_expert_c")
        use tox_shape_truthful_clustering_ensemble_growing, only: grow_ensemble

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
        logical(c_bool), dimension(n_vectors), intent(in), target :: is_member_mask
            !! Current ensemble membership
        real(c_double), intent(in), target :: growth_radius
            !! This ensemble's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(3, n_vectors), intent(out), target :: tmp_range_stack
            !! Workspace: k-d tree traversal stack
        logical(c_bool), dimension(n_vectors), intent(out), target :: tmp_member_mask_buf
            !! Workspace: per-member range-query result
        logical(c_bool), dimension(n_vectors), intent(out), target :: is_member_mask_next
            !! Grown ensemble membership (superset of `is_member_mask`)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_vectors) :: is_member_mask_f
        logical, dimension(n_vectors) :: tmp_member_mask_buf_f
        logical, dimension(n_vectors) :: is_member_mask_next_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(growth_radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(is_member_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_member_mask_buf, n_vectors)
        M_CHECK_ARRAY_NON_NULL(is_member_mask_next, n_vectors)

        is_member_mask_f = is_member_mask

        call grow_ensemble(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            is_member_mask = is_member_mask_f,&
            growth_radius = growth_radius,&
            tmp_range_stack = tmp_range_stack,&
            tmp_member_mask_buf = tmp_member_mask_buf_f,&
            is_member_mask_next = is_member_mask_next_f,&
            ierr = ierr&
        )

        tmp_member_mask_buf = tmp_member_mask_buf_f
        is_member_mask_next = is_member_mask_next_f
    end subroutine grow_ensemble_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_ensemble_growing(module):grow_ensemble_alloc(subroutine)]]
    !| $\mathcal{C}_{\mathcal{E}_{t+1}} = \{x_k \mid \|x_k-x_i\|\le r_{\mathcal{E}} \;\exists\,
    !| x_i\in\mathcal{E}\}$. Deliberately a plain sequential loop, not parallelized: this runs
    !| once per growth iteration per ensemble, and outer-level parallelism across
    !| ensembles/seeds (in `ensemble_identification`) is the right place for that -- a single
    !| ensemble's member count is typically small, especially early in growth. An all-.false.
    !| `is_member_mask` (an empty ensemble) is a well-defined degenerate case: there is nothing
    !| to grow from, so the result is all-.false. too, not an error.
    subroutine grow_ensemble_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            is_member_mask,&
            growth_radius,&
            is_member_mask_next,&
            ierr&
        ) bind(C, name="grow_ensemble_c")
        use tox_shape_truthful_clustering_ensemble_growing, only: grow_ensemble_alloc

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
        logical(c_bool), dimension(n_vectors), intent(in), target :: is_member_mask
            !! Current ensemble membership
        real(c_double), intent(in), target :: growth_radius
            !! This ensemble's growth radius, see `calc_ensemble_growth_radius`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), dimension(n_vectors), intent(out), target :: is_member_mask_next
            !! Grown ensemble membership (superset of `is_member_mask`)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_vectors) :: is_member_mask_f
        logical, dimension(n_vectors) :: is_member_mask_next_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(growth_radius)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(is_member_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(is_member_mask_next, n_vectors)

        is_member_mask_f = is_member_mask

        call grow_ensemble_alloc(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            is_member_mask = is_member_mask_f,&
            growth_radius = growth_radius,&
            is_member_mask_next = is_member_mask_next_f,&
            ierr = ierr&
        )

        is_member_mask_next = is_member_mask_next_f
    end subroutine grow_ensemble_c

end module tox_shape_truthful_clustering_ensemble_growing_c
#endif
