#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_shape_truthful_clustering_parameter_estimation(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_parameter_estimation_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: sample_estimator_anchors_expert_c
    public :: sample_estimator_anchors_c
    public :: grow_estimator_anchor_clouds_c
    public :: estimate_stc_parameters_expert_c
    public :: estimate_stc_parameters_c

contains

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_parameter_estimation(module):sample_estimator_anchors(subroutine)]]
    !| Nearest-rank selection (not `calc_percentile_helper`'s own linear interpolation): these
    !| must be genuine point indices, not interpolated values. Percentiles are
    !| $100/n_{\text{anchors}}, 200/n_{\text{anchors}}, \ldots, 100$ -- e.g. $n_{\text{anchors}}=5$
    !| gives 20/40/60/80/100%ile. Duplicate anchor indices are possible (not deduplicated) when
    !| `n_anchors` is close to `n_vectors`; harmless -- a duplicated anchor's cloud just grows
    !| redundantly, one estimator-anchor "slot" among several effectively wasted, not incorrect.
    subroutine sample_estimator_anchors_expert_c(&
            density_labels,&
            n_vectors,&
            n_anchors,&
            tmp_sort_perm,&
            anchor_indices,&
            ierr&
        ) bind(C, name="sample_estimator_anchors_expert_c")
        use tox_shape_truthful_clustering_parameter_estimation, only: sample_estimator_anchors

        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: a percentile of a single point is
            !! undefined.
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_anchors
            !! Number of estimator anchors (EAs) to pick
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Per-vector density label, see density_labels
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_sort_perm
            !! Workspace: ascending sort permutation of density_labels
        integer(c_int), dimension(n_anchors), intent(out), target :: anchor_indices
            !! Point indices of the n_anchors estimator anchors, ascending-percentile order
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_anchors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_sort_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(anchor_indices, n_anchors)

        call sample_estimator_anchors(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            n_anchors = n_anchors,&
            tmp_sort_perm = tmp_sort_perm,&
            anchor_indices = anchor_indices,&
            ierr = ierr&
        )
    end subroutine sample_estimator_anchors_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_parameter_estimation(module):sample_estimator_anchors_alloc(subroutine)]]
    !| Nearest-rank selection (not `calc_percentile_helper`'s own linear interpolation): these
    !| must be genuine point indices, not interpolated values. Percentiles are
    !| $100/n_{\text{anchors}}, 200/n_{\text{anchors}}, \ldots, 100$ -- e.g. $n_{\text{anchors}}=5$
    !| gives 20/40/60/80/100%ile. Duplicate anchor indices are possible (not deduplicated) when
    !| `n_anchors` is close to `n_vectors`; harmless -- a duplicated anchor's cloud just grows
    !| redundantly, one estimator-anchor "slot" among several effectively wasted, not incorrect.
    subroutine sample_estimator_anchors_c(&
            density_labels,&
            n_vectors,&
            n_anchors,&
            anchor_indices,&
            ierr&
        ) bind(C, name="sample_estimator_anchors_c")
        use tox_shape_truthful_clustering_parameter_estimation, only: sample_estimator_anchors_alloc

        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: a percentile of a single point is
            !! undefined.
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_anchors
            !! Number of estimator anchors (EAs) to pick
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_vectors), intent(in), target :: density_labels
            !! Per-vector density label, see density_labels
        integer(c_int), dimension(n_anchors), intent(out), target :: anchor_indices
            !! Point indices of the n_anchors estimator anchors, ascending-percentile order
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_anchors)
        M_CHECK_ARRAY_NON_NULL(density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(anchor_indices, n_anchors)

        call sample_estimator_anchors_alloc(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            n_anchors = n_anchors,&
            anchor_indices = anchor_indices,&
            ierr = ierr&
        )
    end subroutine sample_estimator_anchors_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_parameter_estimation(module):grow_estimator_anchor_clouds(subroutine)]]
    !| Each EA starts as its own single-point cloud (its anchor). Every round, the globally
    !| closest (unclaimed point, proposing EA) pair -- closest meaning nearest to *any* member
    !| of that EA's own current cloud, not just to the anchor -- is claimed, until either no
    !| unclaimed point remains reachable or the total claimed across all EAs reaches
    !| `seed_max_set_size` percent of `n_vectors`. Implemented as a brute-force rescan every
    !| round directly on `vectors`, not via the k-d tree: `f42_kd_tree` has no "nearest
    !| unclaimed point to a growing region" primitive, and with `n_anchors` small and total
    !| growth capped small by design, the rescan is cheap in absolute terms regardless -- see
    !| `misc/mod_STC.md` for the full complexity discussion and why this is a deliberate
    !| simplicity trade, not an oversight.
    subroutine grow_estimator_anchor_clouds_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            anchor_indices,&
            n_anchors,&
            seed_max_set_size,&
            cloud_masks,&
            cloud_sizes,&
            ierr&
        ) bind(C, name="grow_estimator_anchor_clouds_c")
        use tox_shape_truthful_clustering_parameter_estimation, only: grow_estimator_anchor_clouds

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_anchors
            !! Number of estimator anchors (EAs)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), dimension(n_dimensions, n_vectors), intent(in), target :: vectors
            !! Input data matrix
        integer(c_int), dimension(n_anchors), intent(in), target :: anchor_indices
            !! Point indices of the estimator anchors, see sample_estimator_anchors
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(c_double), intent(in), target :: seed_max_set_size
            !! Percent (0 to 100) of n_vectors at which total growth across all EAs stops
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `5.0_real64`.
        logical(c_bool), dimension(n_vectors, n_anchors), intent(out), target :: cloud_masks
            !! .true. for members of each EA's (column) final cloud, including its own anchor
        integer(c_int), dimension(n_anchors), intent(out), target :: cloud_sizes
            !! Final cloud size per EA
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_vectors, n_anchors) :: cloud_masks_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_anchors)
        M_CHECK_NON_NULL(seed_max_set_size)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(anchor_indices, n_anchors)
        M_CHECK_ARRAY_NON_NULL(cloud_masks, n_vectors * n_anchors)
        M_CHECK_ARRAY_NON_NULL(cloud_sizes, n_anchors)

        call grow_estimator_anchor_clouds(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            anchor_indices = anchor_indices,&
            n_anchors = n_anchors,&
            seed_max_set_size = seed_max_set_size,&
            cloud_masks = cloud_masks_f,&
            cloud_sizes = cloud_sizes,&
            ierr = ierr&
        )

        cloud_masks = cloud_masks_f
    end subroutine grow_estimator_anchor_clouds_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_parameter_estimation(module):estimate_stc_parameters(subroutine)]]
    !| Orchestrates density_labels -> sample_estimator_anchors -> grow_estimator_anchor_clouds
    !| -> observable (once per EA) -> pairwise EA comparisons -> aggregation. See
    !| `misc/mod_STC.md`, SKG `estimate_stc_parameters`, for the full definition of every
    !| output and the reasoning behind each. EAs whose final cloud has fewer than 2 members
    !| (no meaningful SVD possible -- a documented, deliberately unguarded-against possibility
    !| of `grow_estimator_anchor_clouds`'s own stop condition, see there) are excluded from
    !| every statistic below; `ierr` is set if fewer than 2 EAs remain usable (no pairwise
    !| comparison -- and therefore no G_max/d_max estimate -- is possible at all), or if every
    !| usable pair has a zero shared rank (no chordal-distance estimate possible either, even
    !| though G_max/d_max still are). Both are the one genuine, input-shape-dependent runtime
    !| failure this SKG can hit that no simple per-argument DM_* annotation could foresee (it
    !| depends on the data's own spatial distribution, not just
    !| n_vectors/n_anchors/seed_max_set_size in isolation) -- see `codegen_guide.md` section
    !| 5.14.
    subroutine estimate_stc_parameters_expert_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            n_anchors,&
            seed_max_set_size,&
            first_quartile_percentile,&
            lwork_observable,&
            iwork_size,&
            lwork_angle,&
            tmp_neighbors,&
            tmp_distances,&
            tmp_range_stack,&
            tmp_sort_perm,&
            tmp_density_labels,&
            tmp_anchor_indices,&
            tmp_cloud_masks,&
            tmp_cloud_sizes,&
            tmp_y,&
            tmp_s,&
            tmp_u_econ,&
            tmp_vt_econ,&
            tmp_work,&
            tmp_iwork,&
            tmp_angle_m,&
            tmp_angle_s,&
            tmp_angle_work,&
            estimated_k_min,&
            estimated_k_density,&
            estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max,&
            estimated_d_max,&
            ierr&
        ) bind(C, name="estimate_stc_parameters_expert_c")
        use tox_shape_truthful_clustering_parameter_estimation, only: estimate_stc_parameters

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: density_labels' own requirement. Whether
            !! there end up being enough usable estimator anchors is a genuine, data-dependent
            !! runtime condition handled via `ierr` below, not a fixed structural minimum.
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: lwork_observable
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork_observable` output produced by [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):tox_stc_estimate_parameters_svd_workspace]].
        integer(c_int), intent(in), target :: iwork_size
            !! Size of tmp_iwork
            !! It is *VERY IMPORTANT* to compute this argument from the `iwork_size` output produced by [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):tox_stc_estimate_parameters_svd_workspace]].
        integer(c_int), intent(in), target :: lwork_angle
            !! Size of tmp_angle_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork_angle` output produced by [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):tox_stc_estimate_parameters_svd_workspace]].
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
        integer(c_int), intent(in), optional :: k_density
            !! Passed through to density_labels
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        real(c_double), intent(in), optional :: bandwidth_percentile
            !! Passed through to density_labels
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
        integer(c_int), intent(in), target :: n_anchors
            !! Number of estimator anchors (EAs), see sample_estimator_anchors
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
            !! The default value is `5_int32`.
        real(c_double), intent(in), target :: seed_max_set_size
            !! Passed through to grow_estimator_anchor_clouds
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `5.0_real64`.
        real(c_double), intent(in), target :: first_quartile_percentile
            !! Percentile (0 to 100) of the pairwise-EA-comparison distributions used for
            !! chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `25.0_real64`.
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_neighbors
            !! Workspace, see density_labels
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_distances
            !! Workspace, see density_labels
        integer(c_int), dimension(3, n_vectors), intent(out), target :: tmp_range_stack
            !! Workspace, see density_labels
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_sort_perm
            !! Workspace, see density_labels/sample_estimator_anchors
        real(c_double), dimension(n_vectors), intent(out), target :: tmp_density_labels
            !! Workspace: per-vector density labels, see density_labels
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_anchor_indices
            !! Workspace: estimator anchor indices (sized for the worst case, sliced internally)
        logical(c_bool), dimension(n_vectors, n_vectors), intent(out), target :: tmp_cloud_masks
            !! Workspace: EA cloud membership (sized for the worst case, sliced internally)
        integer(c_int), dimension(n_vectors), intent(out), target :: tmp_cloud_sizes
            !! Workspace: EA cloud sizes (sized for the worst case, sliced internally)
        real(c_double), dimension(n_dimensions, n_vectors), intent(out), target :: tmp_y
            !! Workspace, see observable (sized for the worst case)
        real(c_double), dimension(min(n_dimensions,n_vectors)), intent(out), target :: tmp_s
            !! Workspace, see observable (sized for the worst case)
        real(c_double), dimension(n_dimensions, min(n_dimensions,n_vectors)), intent(out), target :: tmp_u_econ
            !! Workspace, see observable (sized for the worst case)
        real(c_double), dimension(min(n_dimensions,n_vectors), n_vectors), intent(out), target :: tmp_vt_econ
            !! Workspace, see observable (sized for the worst case)
        real(c_double), dimension(lwork_observable), intent(out), target :: tmp_work
            !! Workspace, see observable
        integer(c_int), dimension(iwork_size), intent(out), target :: tmp_iwork
            !! Workspace, see observable
        real(c_double), dimension(n_dimensions, n_dimensions), intent(out), target :: tmp_angle_m
            !! Workspace: M = U_i(:,1:d)^T U_j(:,1:d) (sized for the worst case)
        real(c_double), dimension(n_dimensions), intent(out), target :: tmp_angle_s
            !! Workspace: singular values of tmp_angle_m (sized for the worst case)
        real(c_double), dimension(lwork_angle), intent(out), target :: tmp_angle_work
            !! Workspace: LAPACK dgesvd scratch for the principal-angle SVD
        real(c_double), intent(out), target :: estimated_k_min
            !! Estimated k_min (real-valued; round for direct use as an integer argument)
        real(c_double), intent(out), target :: estimated_k_density
            !! Estimated k_density (equal to estimated_k_min, see estimate_stc_parameters)
        real(c_double), intent(out), target :: estimated_density_quantile
            !! Estimated density_quantile -- a literal radius (data units), not a percentile
        real(c_double), intent(out), target :: estimated_chordal_dist_max_as_prcnt_of_range
            !! Estimated chordal_dist_max_as_prcnt_of_range (0 to 1)
        real(c_double), intent(out), target :: estimated_G_max
            !! Estimated G_max
        real(c_double), intent(out), target :: estimated_d_max
            !! Estimated d_max (real-valued; round for direct use as an integer argument)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success
        logical, dimension(n_vectors, n_vectors) :: tmp_cloud_masks_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_anchors)
        M_CHECK_NON_NULL(seed_max_set_size)
        M_CHECK_NON_NULL(first_quartile_percentile)
        M_CHECK_NON_NULL(lwork_observable)
        M_CHECK_NON_NULL(iwork_size)
        M_CHECK_NON_NULL(lwork_angle)
        M_CHECK_NON_NULL(estimated_k_min)
        M_CHECK_NON_NULL(estimated_k_density)
        M_CHECK_NON_NULL(estimated_density_quantile)
        M_CHECK_NON_NULL(estimated_chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(estimated_G_max)
        M_CHECK_NON_NULL(estimated_d_max)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbors, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_sort_perm, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_density_labels, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_anchor_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_cloud_masks, n_vectors * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_cloud_sizes, n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_y, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_s, (min(n_dimensions,n_vectors)))
        M_CHECK_ARRAY_NON_NULL(tmp_u_econ, n_dimensions * (min(n_dimensions,n_vectors)))
        M_CHECK_ARRAY_NON_NULL(tmp_vt_econ, (min(n_dimensions,n_vectors)) * n_vectors)
        M_CHECK_ARRAY_NON_NULL(tmp_work, lwork_observable)
        M_CHECK_ARRAY_NON_NULL(tmp_iwork, iwork_size)
        M_CHECK_ARRAY_NON_NULL(tmp_angle_m, n_dimensions * n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_angle_s, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_angle_work, lwork_angle)

        call estimate_stc_parameters(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            k_density = k_density,&
            bandwidth_percentile = bandwidth_percentile,&
            n_anchors = n_anchors,&
            seed_max_set_size = seed_max_set_size,&
            first_quartile_percentile = first_quartile_percentile,&
            lwork_observable = lwork_observable,&
            iwork_size = iwork_size,&
            lwork_angle = lwork_angle,&
            tmp_neighbors = tmp_neighbors,&
            tmp_distances = tmp_distances,&
            tmp_range_stack = tmp_range_stack,&
            tmp_sort_perm = tmp_sort_perm,&
            tmp_density_labels = tmp_density_labels,&
            tmp_anchor_indices = tmp_anchor_indices,&
            tmp_cloud_masks = tmp_cloud_masks_f,&
            tmp_cloud_sizes = tmp_cloud_sizes,&
            tmp_y = tmp_y,&
            tmp_s = tmp_s,&
            tmp_u_econ = tmp_u_econ,&
            tmp_vt_econ = tmp_vt_econ,&
            tmp_work = tmp_work,&
            tmp_iwork = tmp_iwork,&
            tmp_angle_m = tmp_angle_m,&
            tmp_angle_s = tmp_angle_s,&
            tmp_angle_work = tmp_angle_work,&
            estimated_k_min = estimated_k_min,&
            estimated_k_density = estimated_k_density,&
            estimated_density_quantile = estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range = estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max = estimated_G_max,&
            estimated_d_max = estimated_d_max,&
            ierr = ierr&
        )

        tmp_cloud_masks = tmp_cloud_masks_f
    end subroutine estimate_stc_parameters_expert_c

    !> summary: C-wrapper for [[tox_shape_truthful_clustering_parameter_estimation(module):estimate_stc_parameters_alloc(subroutine)]]
    !| Orchestrates density_labels -> sample_estimator_anchors -> grow_estimator_anchor_clouds
    !| -> observable (once per EA) -> pairwise EA comparisons -> aggregation. See
    !| `misc/mod_STC.md`, SKG `estimate_stc_parameters`, for the full definition of every
    !| output and the reasoning behind each. EAs whose final cloud has fewer than 2 members
    !| (no meaningful SVD possible -- a documented, deliberately unguarded-against possibility
    !| of `grow_estimator_anchor_clouds`'s own stop condition, see there) are excluded from
    !| every statistic below; `ierr` is set if fewer than 2 EAs remain usable (no pairwise
    !| comparison -- and therefore no G_max/d_max estimate -- is possible at all), or if every
    !| usable pair has a zero shared rank (no chordal-distance estimate possible either, even
    !| though G_max/d_max still are). Both are the one genuine, input-shape-dependent runtime
    !| failure this SKG can hit that no simple per-argument DM_* annotation could foresee (it
    !| depends on the data's own spatial distribution, not just
    !| n_vectors/n_anchors/seed_max_set_size in isolation) -- see `codegen_guide.md` section
    !| 5.14.
    subroutine estimate_stc_parameters_c(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            kd_indices,&
            dimension_order,&
            k_density,&
            bandwidth_percentile,&
            n_anchors,&
            seed_max_set_size,&
            first_quartile_percentile,&
            estimated_k_min,&
            estimated_k_density,&
            estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max,&
            estimated_d_max,&
            ierr&
        ) bind(C, name="estimate_stc_parameters_c")
        use tox_shape_truthful_clustering_parameter_estimation, only: estimate_stc_parameters_alloc

        integer(c_int), intent(in), target :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(c_int), intent(in), target :: n_vectors
            !! Number of input vectors N. At least 2: density_labels' own requirement. Whether
            !! there end up being enough usable estimator anchors is a genuine, data-dependent
            !! runtime condition handled via `ierr` below, not a fixed structural minimum.
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
        integer(c_int), intent(in), optional :: k_density
            !! Passed through to density_labels
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        real(c_double), intent(in), optional :: bandwidth_percentile
            !! Passed through to density_labels
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
        integer(c_int), intent(in), target :: n_anchors
            !! Number of estimator anchors (EAs), see sample_estimator_anchors
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
            !! The default value is `5_int32`.
        real(c_double), intent(in), target :: seed_max_set_size
            !! Passed through to grow_estimator_anchor_clouds
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `5.0_real64`.
        real(c_double), intent(in), target :: first_quartile_percentile
            !! Percentile (0 to 100) of the pairwise-EA-comparison distributions used for
            !! chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `25.0_real64`.
        real(c_double), intent(out), target :: estimated_k_min
            !! Estimated k_min (real-valued; round for direct use as an integer argument)
        real(c_double), intent(out), target :: estimated_k_density
            !! Estimated k_density (equal to estimated_k_min, see estimate_stc_parameters)
        real(c_double), intent(out), target :: estimated_density_quantile
            !! Estimated density_quantile -- a literal radius (data units), not a percentile
        real(c_double), intent(out), target :: estimated_chordal_dist_max_as_prcnt_of_range
            !! Estimated chordal_dist_max_as_prcnt_of_range (0 to 1)
        real(c_double), intent(out), target :: estimated_G_max
            !! Estimated G_max
        real(c_double), intent(out), target :: estimated_d_max
            !! Estimated d_max (real-valued; round for direct use as an integer argument)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_anchors)
        M_CHECK_NON_NULL(seed_max_set_size)
        M_CHECK_NON_NULL(first_quartile_percentile)
        M_CHECK_NON_NULL(estimated_k_min)
        M_CHECK_NON_NULL(estimated_k_density)
        M_CHECK_NON_NULL(estimated_density_quantile)
        M_CHECK_NON_NULL(estimated_chordal_dist_max_as_prcnt_of_range)
        M_CHECK_NON_NULL(estimated_G_max)
        M_CHECK_NON_NULL(estimated_d_max)
        M_CHECK_ARRAY_NON_NULL(vectors, n_dimensions * n_vectors)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_vectors)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)

        call estimate_stc_parameters_alloc(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            k_density = k_density,&
            bandwidth_percentile = bandwidth_percentile,&
            n_anchors = n_anchors,&
            seed_max_set_size = seed_max_set_size,&
            first_quartile_percentile = first_quartile_percentile,&
            estimated_k_min = estimated_k_min,&
            estimated_k_density = estimated_k_density,&
            estimated_density_quantile = estimated_density_quantile,&
            estimated_chordal_dist_max_as_prcnt_of_range = estimated_chordal_dist_max_as_prcnt_of_range,&
            estimated_G_max = estimated_G_max,&
            estimated_d_max = estimated_d_max,&
            ierr = ierr&
        )
    end subroutine estimate_stc_parameters_c

end module tox_shape_truthful_clustering_parameter_estimation_c
#endif
