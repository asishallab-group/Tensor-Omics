#include <src/macros.h>

!> # Shape Truthful Clustering (STC): Parameter Estimation
!|
!| A separate, optional pipeline step estimating near-optimal starting values for the crucial
!| parameters (`k_min`, `k_density`, `density_quantile`,
!| `chordal_dist_max_as_prcnt_of_range`, `G_max`, `d_max`) directly from the input data, at a
!| fraction of the cost of a grid search or a
!| resampling-based scheme: grow a handful of "estimator anchors" (EAs) into small local
!| neighborhoods using the same primitives the real pipeline already has
!| (`density_labels`, `observable`), then read the parameters off simple summary statistics of
!| that pass. See `misc/mod_STC.md`, "Estimate parameters from data", for the full algorithm
!| definition and the reasoning behind every design choice below.
module tox_shape_truthful_clustering_parameter_estimation_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_utils, only: sort_real_heapsort, init_perm, calc_percentile_helper
    use tox_errors, only: set_ok, set_err_once, ERR_INTERNAL
    use tox_shape_truthful_clustering_seeding_kernel, only: density_labels_kernel
    use tox_shape_truthful_clustering_observable_kernel, only: observable_kernel, tox_stc_observable_svd_workspace
    M_IMPLICIT_NONE

    interface
        ! Own copy of the dgesvd interface, matching tox_shape_truthful_clustering_accept_kernel's
        ! identical declaration -- module-private interface bodies are not importable via `use`,
        ! so every module that calls a given LAPACK routine declares its own, same as
        ! tox_shape_truthful_clustering_observable_kernel's own dgesdd interface.
        pure subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import :: int32, real64
            character,      intent(in)    :: jobu, jobvt
            integer(int32), intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(real64),   intent(inout) :: a(lda, n)
            real(real64),   intent(out)   :: s(min(m, n))
            real(real64),   intent(out)   :: u(ldu, *)
            real(real64),   intent(out)   :: vt(ldvt, *)
            real(real64),   intent(out)   :: work(lwork)
            integer(int32), intent(out)   :: info
        end subroutine dgesvd
    end interface

#define CM_N_ANCHORS_DEFAULT 5_int32
#define CM_SEED_MAX_SET_SIZE_DEFAULT 5.0_real64
#define CM_FIRST_QUARTILE_PERCENTILE_DEFAULT 25.0_real64

    private
    public :: sample_estimator_anchors_kernel
    public :: grow_estimator_anchor_clouds_kernel
    public :: estimate_stc_parameters_kernel
    public :: tox_stc_estimate_parameters_svd_workspace

contains

    !> summary: Pick n_anchors point indices at evenly-spaced percentiles of the density-sorted order
    !| AUTHOR_ASIS_HALLAB
    !| Nearest-rank selection (not `calc_percentile_helper`'s own linear interpolation): these
    !| must be genuine point indices, not interpolated values. Percentiles are
    !| $100/n_{\text{anchors}}, 200/n_{\text{anchors}}, \ldots, 100$ -- e.g. $n_{\text{anchors}}=5$
    !| gives 20/40/60/80/100%ile. Duplicate anchor indices are possible (not deduplicated) when
    !| `n_anchors` is close to `n_vectors`; harmless -- a duplicated anchor's cloud just grows
    !| redundantly, one estimator-anchor "slot" among several effectively wasted, not incorrect.
    pure subroutine sample_estimator_anchors_kernel(density_labels, n_vectors, n_anchors, &
                                                    tmp_sort_perm, &
                                                    anchor_indices)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: a percentile of a single point is
            !! undefined.
            !! DM_MIN(2_int32)
        real(real64), intent(in) :: density_labels(n_vectors)
            !! Per-vector density label, see density_labels
        integer(int32), intent(in) :: n_anchors
            !! Number of estimator anchors (EAs) to pick
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(out) :: tmp_sort_perm(n_vectors)
            !! Workspace: ascending sort permutation of density_labels
        integer(int32), intent(out) :: anchor_indices(n_anchors)
            !! Point indices of the n_anchors estimator anchors, ascending-percentile order

        integer(int32) :: i, rank_idx
        real(real64)   :: rank_real

        call init_perm(tmp_sort_perm)
        call sort_real_heapsort(density_labels, tmp_sort_perm)

        do i = 1, n_anchors
            rank_real = real(i, real64)/real(n_anchors, real64)*real(n_vectors - 1, real64) + 1.0_real64
            rank_idx  = nint(rank_real)
            rank_idx  = max(1_int32, min(n_vectors, rank_idx))
            anchor_indices(i) = tmp_sort_perm(rank_idx)
        end do

    end subroutine sample_estimator_anchors_kernel

    !> summary: Multi-source competitive region growth of the n_anchors estimator anchors, bounded by seed_max_set_size
    !| AUTHOR_ASIS_HALLAB
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
    pure subroutine grow_estimator_anchor_clouds_kernel(vectors, n_dimensions, n_vectors, &
                                                        anchor_indices, n_anchors, seed_max_set_size, &
                                                        cloud_masks, cloud_sizes)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
            !! DM_MIN(2_int32)
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
            !! Input data matrix
        integer(int32), intent(in) :: n_anchors
            !! Number of estimator anchors (EAs)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        integer(int32), intent(in) :: anchor_indices(n_anchors)
            !! Point indices of the estimator anchors, see sample_estimator_anchors
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors)
        real(real64), intent(in), optional :: seed_max_set_size
            !! Percent (0 to 100) of n_vectors at which total growth across all EAs stops
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(CM_SEED_MAX_SET_SIZE_DEFAULT)
        logical, intent(out) :: cloud_masks(n_vectors, n_anchors)
            !! .true. for members of each EA's (column) final cloud, including its own anchor
        integer(int32), intent(out) :: cloud_sizes(n_anchors)
            !! Final cloud size per EA

        integer(int32) :: actual_max_claims, total_claimed, e, p, q, best_ea, best_point
        real(real64)   :: actual_seed_max_set_size, d2, global_best_d2
        logical        :: claimed(n_vectors)

        M_DEFAULT_VAL(seed_max_set_size, actual_seed_max_set_size, CM_SEED_MAX_SET_SIZE_DEFAULT)

        cloud_masks = .false.
        claimed     = .false.
        do e = 1, n_anchors
            cloud_masks(anchor_indices(e), e) = .true.
            claimed(anchor_indices(e))        = .true.
        end do
        cloud_sizes = 1_int32

        actual_max_claims = max(n_anchors, &
                                ceiling(actual_seed_max_set_size/100.0_real64*real(n_vectors, real64)))
        actual_max_claims = min(actual_max_claims, n_vectors)

        total_claimed = n_anchors
        do while (total_claimed < actual_max_claims)
            best_ea         = 0
            best_point       = 0
            global_best_d2  = huge(1.0_real64)
            do e = 1, n_anchors
                do p = 1, n_vectors
                    if (claimed(p)) cycle
                    do q = 1, n_vectors
                        if (.not. cloud_masks(q, e)) cycle
                        d2 = sum((vectors(:, p) - vectors(:, q))**2)
                        if (d2 < global_best_d2) then
                            global_best_d2 = d2
                            best_ea        = e
                            best_point     = p
                        end if
                    end do
                end do
            end do
            if (best_point == 0) exit
            cloud_masks(best_point, best_ea) = .true.
            claimed(best_point)              = .true.
            cloud_sizes(best_ea)             = cloud_sizes(best_ea) + 1
            total_claimed                    = total_claimed + 1
        end do

    end subroutine grow_estimator_anchor_clouds_kernel

    !> M_EXPORT_C
    !| summary: Recommend LAPACK workspace sizes for estimate_stc_parameters' SVD calls
    !| AUTHOR_ASIS_HALLAB
    !| Worst-case sizing for both `observable`'s dgesdd (an EA cloud can be as large as
    !| n_vectors) and the pairwise principal-angle dgesvd (shared rank at most n_dimensions),
    !| see `tox_stc_observable_svd_workspace` and `tox_stc_accept_ensemble_svd_workspace` for
    !| the individual formulas this combines.
    pure subroutine tox_stc_estimate_parameters_svd_workspace(n_dimensions, n_vectors, &
                                                              lwork_observable, iwork_size, lwork_angle)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N
        integer(int32), intent(out) :: lwork_observable
            !! Recommended size of observable's real LAPACK workspace (worst case)
        integer(int32), intent(out) :: iwork_size
            !! Recommended size of observable's integer LAPACK workspace (worst case)
        integer(int32), intent(out) :: lwork_angle
            !! Recommended size of the pairwise principal-angle LAPACK workspace (worst case)

        call tox_stc_observable_svd_workspace(n_dimensions, n_vectors, lwork_observable, iwork_size)
        lwork_angle = max(1_int32, 5_int32*n_dimensions)

    end subroutine tox_stc_estimate_parameters_svd_workspace

    !> summary: Estimate k_min, k_density, density_quantile, chordal_dist_max_as_prcnt_of_range, G_max, d_max from the data
    !| AUTHOR_ASIS_HALLAB
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
    pure subroutine estimate_stc_parameters_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                                   k_density, bandwidth_percentile, &
                                                   n_anchors, seed_max_set_size, first_quartile_percentile, &
                                                   lwork_observable, iwork_size, lwork_angle, &
                                                   tmp_neighbors, tmp_distances, tmp_range_stack, tmp_sort_perm, &
                                                   tmp_density_labels, tmp_anchor_indices, &
                                                   tmp_cloud_masks, tmp_cloud_sizes, &
                                                   tmp_y, tmp_s, tmp_u_econ, tmp_vt_econ, tmp_work, tmp_iwork, &
                                                   tmp_angle_m, tmp_angle_s, tmp_angle_work, &
                                                   estimated_k_min, estimated_k_density, estimated_density_quantile, &
                                                   estimated_chordal_dist_max_as_prcnt_of_range, estimated_G_max, &
                                                   estimated_d_max, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Ambient dimension D
            !! DM_MIN(2_int32)
        integer(int32), intent(in) :: n_vectors
            !! Number of input vectors N. At least 2: density_labels' own requirement. Whether
            !! there end up being enough usable estimator anchors is a genuine, data-dependent
            !! runtime condition handled via `ierr` below, not a fixed structural minimum.
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
            !! Passed through to density_labels
            !! DM_MIN(1_int32)
            !! DM_MAX(n_vectors - 1_int32)
        real(real64), intent(in), optional :: bandwidth_percentile
            !! Passed through to density_labels
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
        integer(int32), intent(in), optional :: n_anchors
            !! Number of estimator anchors (EAs), see sample_estimator_anchors
            !! DM_MIN(2_int32)
            !! DM_MAX(n_vectors)
            !! DM_DEFAULT(CM_N_ANCHORS_DEFAULT)
        real(real64), intent(in), optional :: seed_max_set_size
            !! Passed through to grow_estimator_anchor_clouds
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(CM_SEED_MAX_SET_SIZE_DEFAULT)
        real(real64), intent(in), optional :: first_quartile_percentile
            !! Percentile (0 to 100) of the pairwise-EA-comparison distributions used for
            !! chordal_dist_max_as_prcnt_of_range/G_max/d_max, see estimate_stc_parameters
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(CM_FIRST_QUARTILE_PERCENTILE_DEFAULT)
        integer(int32), intent(in) :: lwork_observable
            !! Size of tmp_work
            !! DM_OUTPUT_FROM(lwork_observable, tox_stc_estimate_parameters_svd_workspace, tox_shape_truthful_clustering_parameter_estimation_kernel, AUTO)
        integer(int32), intent(in) :: iwork_size
            !! Size of tmp_iwork
            !! DM_OUTPUT_FROM(iwork_size, tox_stc_estimate_parameters_svd_workspace, tox_shape_truthful_clustering_parameter_estimation_kernel, AUTO)
        integer(int32), intent(in) :: lwork_angle
            !! Size of tmp_angle_work
            !! DM_OUTPUT_FROM(lwork_angle, tox_stc_estimate_parameters_svd_workspace, tox_shape_truthful_clustering_parameter_estimation_kernel, AUTO)
        integer(int32), intent(out) :: tmp_neighbors(n_vectors)
            !! Workspace, see density_labels
        real(real64), intent(out) :: tmp_distances(n_vectors)
            !! Workspace, see density_labels
        integer(int32), intent(out) :: tmp_range_stack(3, n_vectors)
            !! Workspace, see density_labels
        integer(int32), intent(out) :: tmp_sort_perm(n_vectors)
            !! Workspace, see density_labels/sample_estimator_anchors
        real(real64), intent(out) :: tmp_density_labels(n_vectors)
            !! Workspace: per-vector density labels, see density_labels
        integer(int32), intent(out) :: tmp_anchor_indices(n_vectors)
            !! Workspace: estimator anchor indices (sized for the worst case, sliced internally)
        logical, intent(out) :: tmp_cloud_masks(n_vectors, n_vectors)
            !! Workspace: EA cloud membership (sized for the worst case, sliced internally)
        integer(int32), intent(out) :: tmp_cloud_sizes(n_vectors)
            !! Workspace: EA cloud sizes (sized for the worst case, sliced internally)
        real(real64), intent(out) :: tmp_y(n_dimensions, n_vectors)
            !! Workspace, see observable (sized for the worst case)
        real(real64), intent(out) :: tmp_s(min(n_dimensions, n_vectors))
            !! Workspace, see observable (sized for the worst case)
        real(real64), intent(out) :: tmp_u_econ(n_dimensions, min(n_dimensions, n_vectors))
            !! Workspace, see observable (sized for the worst case)
        real(real64), intent(out) :: tmp_vt_econ(min(n_dimensions, n_vectors), n_vectors)
            !! Workspace, see observable (sized for the worst case)
        real(real64), intent(out) :: tmp_work(lwork_observable)
            !! Workspace, see observable
        integer(int32), intent(out) :: tmp_iwork(iwork_size)
            !! Workspace, see observable
        real(real64), intent(out) :: tmp_angle_m(n_dimensions, n_dimensions)
            !! Workspace: M = U_i(:,1:d)^T U_j(:,1:d) (sized for the worst case)
        real(real64), intent(out) :: tmp_angle_s(n_dimensions)
            !! Workspace: singular values of tmp_angle_m (sized for the worst case)
        real(real64), intent(out) :: tmp_angle_work(lwork_angle)
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

        integer(int32) :: actual_n_anchors, e, i, j, k, n_valid, d_common, info
        real(real64)   :: u_dummy(1, 1), vt_dummy(1, 1)
        real(real64)   :: actual_first_quartile_percentile
        real(real64)   :: k_vals(n_vectors), dist_vals(n_vectors), median_dist_vals(n_vectors)
        integer(int32) :: d_vals(n_vectors)
        real(real64)   :: G_vals(n_vectors)
        real(real64)   :: U_vals(n_dimensions, n_dimensions, n_vectors)
        integer(int32) :: k_perm(n_vectors), dist_perm(n_vectors)
        integer(int32) :: n_pairs, n_chordal_pairs, pair_perm(n_vectors*(n_vectors - 1)/2)
        real(real64)   :: chordal_vals(n_vectors*(n_vectors - 1)/2)
        real(real64)   :: g_ratio_vals(n_vectors*(n_vectors - 1)/2)
        real(real64)   :: d_diff_vals(n_vectors*(n_vectors - 1)/2)
        real(real64)   :: mu_local(n_dimensions), normal_error_local, tangent_scales_local(n_dimensions)
        real(real64)   :: eigenvalues_local(n_dimensions)
        real(real64)   :: cos_theta, chordal_sum_sq, tmp_median
        integer(int32) :: n_members

        call set_ok(ierr)

        M_DEFAULT_VAL(n_anchors, actual_n_anchors, CM_N_ANCHORS_DEFAULT)
        ! An *explicit* n_anchors is already wrapper-validated against DM_MAX(n_vectors), so
        ! this is a no-op for it -- what it actually guards is CM_N_ANCHORS_DEFAULT itself, a
        ! fixed constant the wrapper never validates against a runtime-dependent bound when
        ! n_anchors is omitted (see misc/code_gen_footgun.md's third entry): a caller on fewer
        ! than 5 points who omits n_anchors would otherwise reach tmp_anchor_indices(1:5) below,
        ! past the end of a workspace array sized n_vectors.
        actual_n_anchors = min(actual_n_anchors, n_vectors)
        M_DEFAULT_VAL(first_quartile_percentile, actual_first_quartile_percentile, CM_FIRST_QUARTILE_PERCENTILE_DEFAULT)

        call density_labels_kernel(vectors, n_dimensions, n_vectors, kd_indices, dimension_order, &
                                   k_density, bandwidth_percentile, &
                                   tmp_neighbors, tmp_distances, tmp_range_stack, tmp_sort_perm, &
                                   tmp_density_labels)

        call sample_estimator_anchors_kernel(tmp_density_labels, n_vectors, actual_n_anchors, &
                                             tmp_sort_perm, tmp_anchor_indices(1:actual_n_anchors))

        call grow_estimator_anchor_clouds_kernel(vectors, n_dimensions, n_vectors, &
                                                 tmp_anchor_indices(1:actual_n_anchors), actual_n_anchors, &
                                                 seed_max_set_size, &
                                                 tmp_cloud_masks(:, 1:actual_n_anchors), &
                                                 tmp_cloud_sizes(1:actual_n_anchors))

        ! --- Per-EA observable, skipping any EA whose cloud never grew past size 1 ---------
        n_valid = 0
        do e = 1, actual_n_anchors
            if (tmp_cloud_sizes(e) < 2) cycle
            n_valid = n_valid + 1

            call observable_kernel(vectors, n_dimensions, n_vectors, tmp_cloud_masks(:, e), tmp_cloud_sizes(e), &
                                   lwork_observable, iwork_size, &
                                   tmp_y(:, 1:tmp_cloud_sizes(e)), tmp_s(1:min(n_dimensions, tmp_cloud_sizes(e))), &
                                   tmp_u_econ(:, 1:min(n_dimensions, tmp_cloud_sizes(e))), &
                                   tmp_vt_econ(1:min(n_dimensions, tmp_cloud_sizes(e)), 1:tmp_cloud_sizes(e)), &
                                   tmp_work, tmp_iwork, &
                                   U_vals(:, :, n_valid), eigenvalues_local, mu_local, &
                                   d_vals(n_valid), G_vals(n_valid), normal_error_local, tangent_scales_local, ierr)
            if (ierr /= 0) return

            k_vals(n_valid) = real(tmp_cloud_sizes(e), real64)

            ! Median distance from this EA's anchor to its own cloud members (excluding itself).
            n_members = 0
            do i = 1, n_vectors
                if (.not. tmp_cloud_masks(i, e)) cycle
                if (i == tmp_anchor_indices(e)) cycle
                n_members = n_members + 1
                dist_vals(n_members) = sqrt(sum((vectors(:, i) - vectors(:, tmp_anchor_indices(e)))**2))
            end do
            call init_perm(dist_perm(1:n_members))
            call sort_real_heapsort(dist_vals(1:n_members), dist_perm(1:n_members))
            call calc_percentile_helper(dist_vals(1:n_members), dist_perm(1:n_members), 50.0_real64, tmp_median)
            median_dist_vals(n_valid) = tmp_median
        end do

        if (n_valid < 2) then
            call set_err_once(ierr, ERR_INTERNAL)
            return
        end if

        ! --- Aggregation: medians over the n_valid EAs -------------------------------------
        call init_perm(k_perm(1:n_valid))
        call sort_real_heapsort(k_vals(1:n_valid), k_perm(1:n_valid))
        call calc_percentile_helper(k_vals(1:n_valid), k_perm(1:n_valid), 50.0_real64, estimated_k_min)
        estimated_k_density = estimated_k_min

        call init_perm(dist_perm(1:n_valid))
        call sort_real_heapsort(median_dist_vals(1:n_valid), dist_perm(1:n_valid))
        call calc_percentile_helper(median_dist_vals(1:n_valid), dist_perm(1:n_valid), 50.0_real64, &
                                    estimated_density_quantile)

        ! --- Aggregation: first_quartile_percentile over all pairs of the n_valid EAs -------
        n_pairs         = 0
        n_chordal_pairs = 0
        do i = 1, n_valid - 1
            do j = i + 1, n_valid
                n_pairs = n_pairs + 1

                d_common = min(d_vals(i), d_vals(j))
                if (d_common > 0) then
                    tmp_angle_m(1:d_common, 1:d_common) = &
                        matmul(transpose(U_vals(:, 1:d_common, i)), U_vals(:, 1:d_common, j))
                    call dgesvd('N', 'N', d_common, d_common, tmp_angle_m(1:d_common, 1:d_common), d_common, &
                               tmp_angle_s(1:d_common), u_dummy, 1, vt_dummy, 1, &
                               tmp_angle_work, lwork_angle, info)
                    if (info /= 0) then
                        call set_err_once(ierr, ERR_INTERNAL)
                        return
                    end if
                    ! Chordal distance across all d_common principal angles (not just the
                    ! worst one), normalized by sqrt(d_common) -- exactly
                    ! chordal_dist_max_as_prcnt_of_range's own definition, see
                    ! misc/mod_STC.md, SKG accept_ensemble, criterion (1). Pairs with
                    ! d_common=0 are excluded from this sample entirely (nothing to compare),
                    ! not counted as a zero -- unlike g_ratio_vals/d_diff_vals below, which
                    ! stay defined for every pair regardless of d_common.
                    n_chordal_pairs = n_chordal_pairs + 1
                    chordal_sum_sq  = 0.0_real64
                    do k = 1, d_common
                        cos_theta      = max(-1.0_real64, min(1.0_real64, tmp_angle_s(k)))
                        chordal_sum_sq = chordal_sum_sq + (1.0_real64 - cos_theta**2)
                    end do
                    chordal_vals(n_chordal_pairs) = sqrt(chordal_sum_sq)/sqrt(real(d_common, real64))
                end if

                g_ratio_vals(n_pairs) = abs(log(G_vals(j)/G_vals(i)))
                d_diff_vals(n_pairs)  = real(abs(d_vals(j) - d_vals(i)), real64)
            end do
        end do

        if (n_chordal_pairs < 1) then
            call set_err_once(ierr, ERR_INTERNAL)
            return
        end if

        call init_perm(pair_perm(1:n_chordal_pairs))
        call sort_real_heapsort(chordal_vals(1:n_chordal_pairs), pair_perm(1:n_chordal_pairs))
        call calc_percentile_helper(chordal_vals(1:n_chordal_pairs), pair_perm(1:n_chordal_pairs), &
                                    actual_first_quartile_percentile, estimated_chordal_dist_max_as_prcnt_of_range)

        call init_perm(pair_perm(1:n_pairs))
        call sort_real_heapsort(g_ratio_vals(1:n_pairs), pair_perm(1:n_pairs))
        call calc_percentile_helper(g_ratio_vals(1:n_pairs), pair_perm(1:n_pairs), &
                                    actual_first_quartile_percentile, estimated_G_max)

        call init_perm(pair_perm(1:n_pairs))
        call sort_real_heapsort(d_diff_vals(1:n_pairs), pair_perm(1:n_pairs))
        call calc_percentile_helper(d_diff_vals(1:n_pairs), pair_perm(1:n_pairs), &
                                    actual_first_quartile_percentile, estimated_d_max)

    end subroutine estimate_stc_parameters_kernel

end module tox_shape_truthful_clustering_parameter_estimation_kernel
