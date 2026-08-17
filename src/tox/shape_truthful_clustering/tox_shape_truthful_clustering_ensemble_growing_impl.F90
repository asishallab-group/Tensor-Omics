#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Ensemble Growing
!|
!| Kernels for growing an ensemble by one step: `calc_ensemble_growth_radius` (the
!| per-seed, locally adaptive growth radius, computed once) and `grow_ensemble` (the union,
!| over every current member, of the points within that radius). See `misc/mod_STC.md`,
!| sections "Local Radius Identification" and SKG `grow_ensemble`, for the full algorithm
!| definition. Both take an already-built k-d tree (`kd_indices`, `dimension_order`, see
!| [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]]) as input, the same one shared
!| with the seeding kernels, rather than building their own.
module tox_shape_truthful_clustering_ensemble_growing_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_sort_impl, only: sort_real_heapsort
    use f42_stats_impl, only: calc_percentile_impl
    use f42_kd_tree_impl, only: kd_knn_query_impl, kd_range_query_mask_impl
    M_IMPLICIT_NONE

#define CM_GROWTH_RADIUS_K_MIN_DEFAULT 30_int32
#define CM_GROWTH_RADIUS_PERCENTILE_DEFAULT 50.0_real64

    private
    public :: calc_ensemble_growth_radius_impl
    public :: grow_ensemble_impl

contains

    !> summary: Locally adapted ensemble growth radius, a percentile of the distances among a seed's own k_min nearest neighbors
    !| AUTHOR_ASIS_HALLAB
    !| Matches LoManLe's `local_scale_i` exactly at the default `radius_percentile=50.0` (the
    !| median) -- a per-seed, locally adaptive radius rather than a single dataset-wide one
    !| (see `misc/STC_for_LoManLe.md` section 2.2). Work arrays are sized for the worst case
    !| (`k_min = n_vectors - 1`) and sliced internally, since `k_min`'s resolved value is only
    !| known once its default (if any) has been applied.
    !|
    !| `radius_percentile` generalizes what used to be a hardcoded median: `seeds_impl`
    !| reuses this same SKG for its seed-exclusion radius (see `misc/mod_STC.md`, SKG
    !| `seeds`), and a fixed dataset-wide-in-spirit median there was observed to suppress
    !| entire uncovered regions (e.g. curvature extrema on a wavy manifold) whose own growth
    !| never actually reaches that far -- exposing the percentile lets that specific call site
    !| shrink its exclusion radius independently of this kernel's own default, without
    !| touching the actual growth-radius computation any other caller relies on.
    pure subroutine calc_ensemble_growth_radius_impl(vectors, n_dimensions, n_vectors, &
                                                        kd_indices, dimension_order, seed_index, k_min, &
                                                        radius_percentile, &
                                                        tmp_neighbors, tmp_distances, tmp_range_stack, tmp_sort_perm, &
                                                        growth_radius)
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
        integer(int32), intent(in) :: seed_index
            !! Index into `vectors`/`kd_indices` of the seed to compute the growth radius for
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in), optional :: k_min
            !! Neighborhood size the median distance is taken over
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors - 1_int32)
            !! DM_DEFAULT(CM_GROWTH_RADIUS_K_MIN_DEFAULT)
        real(real64), intent(in), optional :: radius_percentile
            !! Percentile (0 to 100) of the k_min neighbor distances reported as the growth
            !! radius -- 50.0 (the default) is the median, matching this SKG's original,
            !! non-parameterized behavior
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(CM_GROWTH_RADIUS_PERCENTILE_DEFAULT)
        integer(int32), intent(out) :: tmp_neighbors(n_vectors)
            !! Workspace: k-NN query result, indices (sized for the worst case, sliced internally)
        real(real64), intent(out) :: tmp_distances(n_vectors)
            !! Workspace: k-NN query result, distances (sized as `tmp_neighbors`)
        integer(int32), intent(out) :: tmp_range_stack(3, n_vectors)
            !! Workspace: k-d tree traversal stack, see `kd_knn_query`
        integer(int32), intent(out) :: tmp_sort_perm(n_vectors)
            !! Workspace: sort-permutation scratch (sized as `tmp_neighbors`)
        real(real64), intent(out) :: growth_radius
            !! radius_percentile-th percentile of the distances among the seed's own k_min
            !! nearest neighbors

        integer(int32) :: actual_k_min, k_query, self_pos, j
        real(real64)   :: actual_radius_percentile

        M_DEFAULT_VAL(k_min, actual_k_min, CM_GROWTH_RADIUS_K_MIN_DEFAULT)
        ! An *explicit* k_min is already wrapper-validated against DM_MAX(n_vectors - 1), so
        ! this is a no-op for it -- what it actually guards is CM_GROWTH_RADIUS_K_MIN_DEFAULT
        ! itself, a fixed constant that the wrapper never validates against a runtime-dependent
        ! bound at all when k_min is omitted (see misc/code_gen_footgun.md's third entry): a
        ! caller on fewer than 31 points who omits k_min would otherwise reach the k-NN query
        ! below asking for more neighbors than exist, past the end of tmp_neighbors/tmp_distances.
        actual_k_min = min(actual_k_min, n_vectors - 1)
        k_query = actual_k_min + 1

        M_DEFAULT_VAL(radius_percentile, actual_radius_percentile, CM_GROWTH_RADIUS_PERCENTILE_DEFAULT)

        call kd_knn_query_impl(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                               vectors(:, seed_index), k_query, tmp_range_stack, &
                               tmp_neighbors(1:k_query), tmp_distances(1:k_query))

        ! Exclude the seed itself, at distance 0, from its own k_query nearest neighbors: swap
        ! it into the last slot, then only the first actual_k_min entries are used from here
        ! on. kd_knn_query does not guarantee its output is sorted by distance (max-heap based
        ! internally), so this swap-with-last is safe irrespective of result order.
        self_pos = 0
        do j = 1, k_query
            if (tmp_neighbors(j) == seed_index) then
                self_pos = j
                exit
            end if
        end do
        if (self_pos > 0 .and. self_pos < k_query) then
            tmp_neighbors(self_pos) = tmp_neighbors(k_query)
            tmp_distances(self_pos) = tmp_distances(k_query)
        end if

        ! Sort the remaining actual_k_min distances ascending so the percentile below is meaningful.
        do j = 1, actual_k_min
            tmp_sort_perm(j) = j
        end do
        call sort_real_heapsort(tmp_distances(1:actual_k_min), tmp_sort_perm(1:actual_k_min))

        ! calc_percentile_impl takes its percentile as a [0,1] fraction; this SKG's own
        ! radius_percentile is documented and validated as 0-100, so it is rescaled at the call
        ! site rather than changing what every caller of this kernel passes.
        call calc_percentile_impl(tmp_distances(1:actual_k_min), actual_k_min, tmp_sort_perm(1:actual_k_min), &
                                  actual_radius_percentile/100.0_real64, growth_radius)

    end subroutine calc_ensemble_growth_radius_impl

    !> summary: Grow an ensemble by one step, the union of every current member's growth-radius neighborhood
    !| AUTHOR_ASIS_HALLAB
    !| $\mathcal{C}_{\mathcal{E}_{t+1}} = \{x_k \mid \|x_k-x_i\|\le r_{\mathcal{E}} \;\exists\,
    !| x_i\in\mathcal{E}\}$. Deliberately a plain sequential loop, not parallelized: this runs
    !| once per growth iteration per ensemble, and outer-level parallelism across
    !| ensembles/seeds (in `ensemble_identification`) is the right place for that -- a single
    !| ensemble's member count is typically small, especially early in growth. An all-.false.
    !| `is_member_mask` (an empty ensemble) is a well-defined degenerate case: there is nothing
    !| to grow from, so the result is all-.false. too, not an error.
    pure subroutine grow_ensemble_impl(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                         is_member_mask, growth_radius, &
                                         tmp_range_stack, tmp_member_mask_buf, &
                                         is_member_mask_next)
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
        logical(c_bool), intent(in) :: is_member_mask(n_vectors)
            !! Current ensemble membership
        real(real64), intent(in) :: growth_radius
            !! This ensemble's growth radius, see `calc_ensemble_growth_radius`
            !! DM_MIN(0.0_real64)
        integer(int32), intent(out) :: tmp_range_stack(3, n_vectors)
            !! Workspace: k-d tree traversal stack
        logical(c_bool), intent(out) :: tmp_member_mask_buf(n_vectors)
            !! Workspace: per-member range-query result
        logical(c_bool), intent(out) :: is_member_mask_next(n_vectors)
            !! Grown ensemble membership (superset of `is_member_mask`)

        integer(int32) :: i

        is_member_mask_next = is_member_mask

        do i = 1, n_vectors
            if (.not. is_member_mask(i)) cycle
            call kd_range_query_mask_impl(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                          vectors(:, i), growth_radius, tmp_range_stack, tmp_member_mask_buf)
            is_member_mask_next = is_member_mask_next .or. tmp_member_mask_buf
        end do

    end subroutine grow_ensemble_impl

end module tox_shape_truthful_clustering_ensemble_growing_impl
