#include <src/macros.h>

!> summary: Wrappers for [[tox_shape_truthful_clustering_parameter_estimation_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_shape_truthful_clustering_parameter_estimation
    use tox_shape_truthful_clustering_parameter_estimation_kernel, only: estimate_stc_parameters_kernel, grow_estimator_anchor_clouds_kernel, sample_estimator_anchors_kernel, tox_stc_estimate_parameters_svd_workspace
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size
    use tox_errors, only: validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: sample_estimator_anchors
    public :: sample_estimator_anchors_alloc
    public :: grow_estimator_anchor_clouds
    public :: estimate_stc_parameters
    public :: estimate_stc_parameters_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):sample_estimator_anchors_kernel]].
    !| Nearest-rank selection (not `calc_percentile_helper`'s own linear interpolation): these
    !| must be genuine point indices, not interpolated values. Percentiles are
    !| $100/n_{\text{anchors}}, 200/n_{\text{anchors}}, \ldots, 100$ -- e.g. $n_{\text{anchors}}=5$
    !| gives 20/40/60/80/100%ile. Duplicate anchor indices are possible (not deduplicated) when
    !| `n_anchors` is close to `n_vectors`; harmless -- a duplicated anchor's cloud just grows
    !| redundantly, one estimator-anchor "slot" among several effectively wasted, not incorrect.
    subroutine sample_estimator_anchors(&
            density_labels,&
            n_vectors,&
            n_anchors,&
            tmp_sort_perm,&
            anchor_indices,&
            ierr&
        )
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a percentile of a single point is
            !! undefined.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_anchors
            !! Number of estimator anchors (EAs) to pick
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Per-vector density label, see density_labels
        integer(int32), dimension(n_vectors), intent(out) :: tmp_sort_perm
            !! Workspace: ascending sort permutation of density_labels
        integer(int32), dimension(n_anchors), intent(out) :: anchor_indices
            !! Point indices of the n_anchors estimator anchors, ascending-percentile order
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_vectors, ierr, arg_pos=2_int32, min=2_int32)
        call validate_in_range_int(n_anchors, ierr, arg_pos=3_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call sample_estimator_anchors_kernel(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            n_anchors = n_anchors,&
            tmp_sort_perm = tmp_sort_perm,&
            anchor_indices = anchor_indices&
        )
    end subroutine sample_estimator_anchors

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):sample_estimator_anchors_kernel]].
    !| Nearest-rank selection (not `calc_percentile_helper`'s own linear interpolation): these
    !| must be genuine point indices, not interpolated values. Percentiles are
    !| $100/n_{\text{anchors}}, 200/n_{\text{anchors}}, \ldots, 100$ -- e.g. $n_{\text{anchors}}=5$
    !| gives 20/40/60/80/100%ile. Duplicate anchor indices are possible (not deduplicated) when
    !| `n_anchors` is close to `n_vectors`; harmless -- a duplicated anchor's cloud just grows
    !| redundantly, one estimator-anchor "slot" among several effectively wasted, not incorrect.
    subroutine sample_estimator_anchors_alloc(&
            density_labels,&
            n_vectors,&
            n_anchors,&
            anchor_indices,&
            ierr&
        )
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a percentile of a single point is
            !! undefined.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_anchors
            !! Number of estimator anchors (EAs) to pick
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_vectors), intent(in) :: density_labels
            !! Per-vector density label, see density_labels
        integer(int32), dimension(n_anchors), intent(out) :: anchor_indices
            !! Point indices of the n_anchors estimator anchors, ascending-percentile order
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_sort_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_vectors, ierr, arg_pos=2_int32, min=2_int32)
        call validate_in_range_int(n_anchors, ierr, arg_pos=3_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_sort_perm(n_vectors))

        call sample_estimator_anchors_kernel(&
            density_labels = density_labels,&
            n_vectors = n_vectors,&
            n_anchors = n_anchors,&
            tmp_sort_perm = tmp_sort_perm,&
            anchor_indices = anchor_indices&
        )
    end subroutine sample_estimator_anchors_alloc

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):grow_estimator_anchor_clouds_kernel]].
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
    subroutine grow_estimator_anchor_clouds(&
            vectors,&
            n_dimensions,&
            n_vectors,&
            anchor_indices,&
            n_anchors,&
            seed_max_set_size,&
            cloud_masks,&
            cloud_sizes,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_anchors
            !! Number of estimator anchors (EAs)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input data matrix
        integer(int32), dimension(n_anchors), intent(in) :: anchor_indices
            !! Point indices of the estimator anchors, see sample_estimator_anchors
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors`.
        real(real64), intent(in), optional :: seed_max_set_size
            !! Percent (0 to 100) of n_vectors at which total growth across all EAs stops
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `5.0_real64`.
        logical, dimension(n_vectors, n_anchors), intent(out) :: cloud_masks
            !! .true. for members of each EA's (column) final cloud, including its own anchor
        integer(int32), dimension(n_anchors), intent(out) :: cloud_sizes
            !! Final cloud size per EA
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(n_anchors, ierr, arg_pos=5_int32, min=1_int32, max=n_vectors)
        call validate_in_range_real(seed_max_set_size, ierr, arg_pos=6_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(anchor_indices, n_anchors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return
#endif

        call grow_estimator_anchor_clouds_kernel(&
            vectors = vectors,&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            anchor_indices = anchor_indices,&
            n_anchors = n_anchors,&
            seed_max_set_size = seed_max_set_size,&
            cloud_masks = cloud_masks,&
            cloud_sizes = cloud_sizes&
        )
    end subroutine grow_estimator_anchor_clouds

    !> summary: Validates its inputs, then calls [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):estimate_stc_parameters_kernel]].
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
    subroutine estimate_stc_parameters(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: density_labels' own requirement. Whether
            !! there end up being enough usable estimator anchors is a genuine, data-dependent
            !! runtime condition handled via `ierr` below, not a fixed structural minimum.
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: lwork_observable
            !! Size of tmp_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork_observable` output produced by [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):tox_stc_estimate_parameters_svd_workspace]].
        integer(int32), intent(in) :: iwork_size
            !! Size of tmp_iwork
            !! It is *VERY IMPORTANT* to compute this argument from the `iwork_size` output produced by [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):tox_stc_estimate_parameters_svd_workspace]].
        integer(int32), intent(in) :: lwork_angle
            !! Size of tmp_angle_work
            !! It is *VERY IMPORTANT* to compute this argument from the `lwork_angle` output produced by [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):tox_stc_estimate_parameters_svd_workspace]].
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
            !! Passed through to density_labels
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Passed through to density_labels
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
        integer(int32), intent(in), optional :: n_anchors
            !! Number of estimator anchors (EAs), see sample_estimator_anchors
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
            !! The default value is `5_int32`.
        real(real64), intent(in), optional :: seed_max_set_size
            !! Passed through to grow_estimator_anchor_clouds
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `5.0_real64`.
        real(real64), intent(in), optional :: first_quartile_percentile
            !! Percentile (0 to 100) of the pairwise-EA-comparison distributions used for
            !! chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `25.0_real64`.
        integer(int32), dimension(n_vectors), intent(out) :: tmp_neighbors
            !! Workspace, see density_labels
        real(real64), dimension(n_vectors), intent(out) :: tmp_distances
            !! Workspace, see density_labels
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_range_stack
            !! Workspace, see density_labels
        integer(int32), dimension(n_vectors), intent(out) :: tmp_sort_perm
            !! Workspace, see density_labels/sample_estimator_anchors
        real(real64), dimension(n_vectors), intent(out) :: tmp_density_labels
            !! Workspace: per-vector density labels, see density_labels
        integer(int32), dimension(n_vectors), intent(out) :: tmp_anchor_indices
            !! Workspace: estimator anchor indices (sized for the worst case, sliced internally)
        logical, dimension(n_vectors, n_vectors), intent(out) :: tmp_cloud_masks
            !! Workspace: EA cloud membership (sized for the worst case, sliced internally)
        integer(int32), dimension(n_vectors), intent(out) :: tmp_cloud_sizes
            !! Workspace: EA cloud sizes (sized for the worst case, sliced internally)
        real(real64), dimension(n_dimensions, n_vectors), intent(out) :: tmp_y
            !! Workspace, see observable (sized for the worst case)
        real(real64), dimension(min(n_dimensions,n_vectors)), intent(out) :: tmp_s
            !! Workspace, see observable (sized for the worst case)
        real(real64), dimension(n_dimensions, min(n_dimensions,n_vectors)), intent(out) :: tmp_u_econ
            !! Workspace, see observable (sized for the worst case)
        real(real64), dimension(min(n_dimensions,n_vectors), n_vectors), intent(out) :: tmp_vt_econ
            !! Workspace, see observable (sized for the worst case)
        real(real64), dimension(lwork_observable), intent(out) :: tmp_work
            !! Workspace, see observable
        integer(int32), dimension(iwork_size), intent(out) :: tmp_iwork
            !! Workspace, see observable
        real(real64), dimension(n_dimensions, n_dimensions), intent(out) :: tmp_angle_m
            !! Workspace: M = U_i(:,1:d)^T U_j(:,1:d) (sized for the worst case)
        real(real64), dimension(n_dimensions), intent(out) :: tmp_angle_s
            !! Workspace: singular values of tmp_angle_m (sized for the worst case)
        real(real64), dimension(lwork_angle), intent(out) :: tmp_angle_work
            !! Workspace: LAPACK dgesvd scratch for the principal-angle SVD
        real(real64), intent(out) :: estimated_k_min
            !! Estimated k_min (real-valued; round for direct use as an integer argument)
        real(real64), intent(out) :: estimated_k_density
            !! Estimated k_density (equal to estimated_k_min, see estimate_stc_parameters)
        real(real64), intent(out) :: estimated_density_quantile
            !! Estimated density_quantile -- a literal radius (data units), not a percentile
        real(real64), intent(out) :: estimated_chordal_dist_max_as_prcnt_of_range
            !! Estimated chordal_dist_max_as_prcnt_of_range (0 to 1)
        real(real64), intent(out) :: estimated_G_max
            !! Estimated G_max
        real(real64), intent(out) :: estimated_d_max
            !! Estimated d_max (real-valued; round for direct use as an integer argument)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=2_int32, min=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_density, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(bandwidth_percentile, ierr, arg_pos=7_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_int(n_anchors, ierr, arg_pos=8_int32, min=2_int32, max=n_vectors)
        call validate_in_range_real(seed_max_set_size, ierr, arg_pos=9_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(first_quartile_percentile, ierr, arg_pos=10_int32, min=0.0_real64, max=100.0_real64)
        call validate_dimension_size(lwork_observable, ierr, arg_pos=11_int32)
        call validate_dimension_size(iwork_size, ierr, arg_pos=12_int32)
        call validate_dimension_size(lwork_angle, ierr, arg_pos=13_int32)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call estimate_stc_parameters_kernel(&
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
            tmp_cloud_masks = tmp_cloud_masks,&
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
        call clear_err_arg_pos(ierr)
    end subroutine estimate_stc_parameters

    !> summary: Allocates its work arrays, then calls [[tox_shape_truthful_clustering_parameter_estimation_kernel(module):estimate_stc_parameters_kernel]].
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
    subroutine estimate_stc_parameters_alloc(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! The minimum valid value is `2_int32`.
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: density_labels' own requirement. Whether
            !! there end up being enough usable estimator anchors is a genuine, data-dependent
            !! runtime condition handled via `ierr` below, not a fixed structural minimum.
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
            !! Passed through to density_labels
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_vectors - 1_int32`.
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Passed through to density_labels
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
        integer(int32), intent(in), optional :: n_anchors
            !! Number of estimator anchors (EAs), see sample_estimator_anchors
            !! The minimum valid value is `2_int32`.
            !! The maximum valid value is `n_vectors`.
            !! The default value is `5_int32`.
        real(real64), intent(in), optional :: seed_max_set_size
            !! Passed through to grow_estimator_anchor_clouds
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `5.0_real64`.
        real(real64), intent(in), optional :: first_quartile_percentile
            !! Percentile (0 to 100) of the pairwise-EA-comparison distributions used for
            !! chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `25.0_real64`.
        real(real64), intent(out) :: estimated_k_min
            !! Estimated k_min (real-valued; round for direct use as an integer argument)
        real(real64), intent(out) :: estimated_k_density
            !! Estimated k_density (equal to estimated_k_min, see estimate_stc_parameters)
        real(real64), intent(out) :: estimated_density_quantile
            !! Estimated density_quantile -- a literal radius (data units), not a percentile
        real(real64), intent(out) :: estimated_chordal_dist_max_as_prcnt_of_range
            !! Estimated chordal_dist_max_as_prcnt_of_range (0 to 1)
        real(real64), intent(out) :: estimated_G_max
            !! Estimated G_max
        real(real64), intent(out) :: estimated_d_max
            !! Estimated d_max (real-valued; round for direct use as an integer argument)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success
        integer(int32) :: lwork_observable
        integer(int32) :: iwork_size
        integer(int32) :: lwork_angle
        integer(int32), dimension(:), allocatable :: tmp_neighbors
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack
        integer(int32), dimension(:), allocatable :: tmp_sort_perm
        real(real64), dimension(:), allocatable :: tmp_density_labels
        integer(int32), dimension(:), allocatable :: tmp_anchor_indices
        logical, dimension(:, :), allocatable :: tmp_cloud_masks
        integer(int32), dimension(:), allocatable :: tmp_cloud_sizes
        real(real64), dimension(:, :), allocatable :: tmp_y
        real(real64), dimension(:), allocatable :: tmp_s
        real(real64), dimension(:, :), allocatable :: tmp_u_econ
        real(real64), dimension(:, :), allocatable :: tmp_vt_econ
        real(real64), dimension(:), allocatable :: tmp_work
        integer(int32), dimension(:), allocatable :: tmp_iwork
        real(real64), dimension(:, :), allocatable :: tmp_angle_m
        real(real64), dimension(:), allocatable :: tmp_angle_s
        real(real64), dimension(:), allocatable :: tmp_angle_work

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_dimensions, ierr, arg_pos=2_int32, min=2_int32)
        call validate_in_range_int(n_vectors, ierr, arg_pos=3_int32, min=2_int32)
        call validate_in_range_int(k_density, ierr, arg_pos=6_int32, min=1_int32, max=n_vectors - 1_int32)
        call validate_in_range_real(bandwidth_percentile, ierr, arg_pos=7_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_int(n_anchors, ierr, arg_pos=8_int32, min=2_int32, max=n_vectors)
        call validate_in_range_real(seed_max_set_size, ierr, arg_pos=9_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(first_quartile_percentile, ierr, arg_pos=10_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(vectors, n_dimensions * n_vectors, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, arg_pos=4_int32, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call tox_stc_estimate_parameters_svd_workspace(&
            n_dimensions = n_dimensions,&
            n_vectors = n_vectors,&
            lwork_observable = lwork_observable,&
            iwork_size = iwork_size,&
            lwork_angle = lwork_angle&
        )
        M_ALLOCATE(tmp_neighbors(n_vectors))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_range_stack(3, n_vectors))
        M_ALLOCATE(tmp_sort_perm(n_vectors))
        M_ALLOCATE(tmp_density_labels(n_vectors))
        M_ALLOCATE(tmp_anchor_indices(n_vectors))
        M_ALLOCATE(tmp_cloud_masks(n_vectors, n_vectors))
        M_ALLOCATE(tmp_cloud_sizes(n_vectors))
        M_ALLOCATE(tmp_y(n_dimensions, n_vectors))
        M_ALLOCATE(tmp_s(min(n_dimensions,n_vectors)))
        M_ALLOCATE(tmp_u_econ(n_dimensions, min(n_dimensions,n_vectors)))
        M_ALLOCATE(tmp_vt_econ(min(n_dimensions,n_vectors), n_vectors))
        M_ALLOCATE(tmp_work(lwork_observable))
        M_ALLOCATE(tmp_iwork(iwork_size))
        M_ALLOCATE(tmp_angle_m(n_dimensions, n_dimensions))
        M_ALLOCATE(tmp_angle_s(n_dimensions))
        M_ALLOCATE(tmp_angle_work(lwork_angle))

        call estimate_stc_parameters_kernel(&
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
            tmp_cloud_masks = tmp_cloud_masks,&
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
        call clear_err_arg_pos(ierr)
    end subroutine estimate_stc_parameters_alloc

end module tox_shape_truthful_clustering_parameter_estimation
