#include <src/macros.h>

!> Local manifold learning (LoManLe): builds a skeleton/backbone through a
!| point cloud by fitting a local atlas of overlapping tangent-plane anchors,
!| stitching points within each overlap by inverse-variance-weighted
!| projection, and then reading a topological backbone graph (endpoints,
!| pass-through anchors, branch/junction anchors) off the anchor adjacency
!| graph's minimum spanning tree. See misc/smoothing_experiments.md for the
!| full algorithm write-up and misc/lomanle_literature_resume.md for how this
!| compares to LTSA/Hidalgo/MMLS/Atlas/IAN/SCMS.
module lomanle_mod
    use iso_fortran_env, only: int32, real64
    use kd_tree,         only: build_kd_index, kd_knn_query
    use tox_errors,      only: set_ok, is_ok, is_err, set_err, ERR_INVALID_INPUT, ERR_ALLOC_FAIL, &
                               validate_dimension_size, validate_in_range_real, &
                               validate_in_range_int, validate_all_in_range_real
    use f42_utils,       only: sort_array, above
    implicit none

    interface
        ! Declared pure so it can be called from do concurrent / pure contexts
        ! (grow_one_point_neighborhood, compute_anchor_svd): dsyev is a
        ! deterministic numerical routine with no I/O and no global state.
        pure subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: int32, real64
            character,        intent(in)    :: jobz, uplo
            integer(int32),   intent(in)    :: n, lda, lwork
            real(real64),     intent(inout) :: a(lda, n)
            real(real64),     intent(out)   :: w(n)
            real(real64),     intent(out)   :: work(lwork)
            integer(int32),   intent(out)   :: info
        end subroutine dsyev
    end interface

    private :: lomanle_pass, grow_adaptive_neighborhoods, grow_one_point_neighborhood, &
               construct_atlas, sort_points_by_density, select_atlas_anchors, absorb_orphans, &
               compute_anchor_svd, build_membership_matrix, &
               count_anchor_intersection_edges, fill_anchor_intersection_edges, &
               count_intersection_csr_sizes, build_intersection_graph_alloc, &
               build_intersection_graph, stitch_points, &
               stitch_multi_anchor_point, stitch_single_anchor_point, &
               build_skeleton_edges_alloc, build_skeleton_edges, &
               build_anchor_mapping, mark_tier1_candidate_pairs, count_tier2_candidate_pairs, &
               build_anchor_mst_alloc, &
               build_anchor_mst, build_tier1_mst_edges, build_tier2_mst_edges, &
               classify_anchor_roles, build_member_chains, build_branch_adjacency, &
               find_root, nearest_member, emit_branch, compute_relative_conv_tol

contains

    !> Orchestrates the LoManLe pipeline for a single pass, by calling each
    !| numbered step (see the individual step subroutines below) in sequence.
    !| Not pure: Steps 6.5-9 call build_intersection_graph_alloc, which
    !| allocates. The point_in_anchor_mask/anchor_to_point locals are
    !| allocatable only because Fortran requires an allocatable actual
    !| argument for build_intersection_graph_alloc's allocatable intent(out)
    !| dummies of the same name; this routine itself never calls allocate(),
    !| so it does not carry an _alloc suffix.
    subroutine lomanle_pass(work_coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                            o_max, o_min, stability_threshold, scale_factor, &
                            tmp_kd_indices, tmp_workspace, tmp_val_buf, tmp_perm, &
                            tmp_l_stack, tmp_r_stack, tmp_rec_stack, &
                            tmp_dim_order, lwork, &
                            sphere_radii, densities, gap_values, &
                            normal_errors, stability_values, k_selected, growth_stopped_complex, &
                            is_anchor_mask, tangent_bases, tangent_scales, labels, skeleton_coords, &
                            primary_anchor_ids, secondary_anchor_ids, anchor_centers, ierr)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: k_min
            !! Minimum neighborhood size to start adaptive growth from
        integer(int32), intent(in) :: lwork
            !! Size of the LAPACK dsyev scratch workspace
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        real(real64),   intent(in) :: g_threshold
            !! Spectral-gap threshold that stops neighborhood growth
        real(real64),   intent(in) :: o_max
            !! Maximum allowed anchor-sphere overlap ratio (Step 5)
        real(real64),   intent(in) :: o_min
            !! Minimum required anchor-sphere overlap ratio (Step 5)
        real(real64),   intent(in) :: stability_threshold
            !! Minimum tangent-basis stability to keep growing a neighborhood
        real(real64),   intent(in) :: scale_factor
            !! Cap on sphere radius as a multiple of the point's local scale
        real(real64), intent(out) :: skeleton_coords(dim+1, n_points)
            !! Row 1 = anchor_count; rows 2:dim+1 = stitched position (Step 10)
        real(real64), intent(out) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (Step 6); anchor_centers(:,i) is only
            !! meaningful where is_anchor_mask(i) is .true. Exposed so callers can build a
            !! topological backbone graph over anchors (see build_skeleton_edges_alloc)
            !! without having to recompute or re-derive anchor positions.

        integer(int32), intent(inout) :: tmp_kd_indices(n_points)
            !! KD-tree scratch buffer, see kd_tree module
        integer(int32), intent(inout) :: tmp_workspace(n_points)
            !! KD-tree/sort scratch buffer, reused across steps
        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer, see f42_utils sort_array
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! KD-tree/sort_array recursion-stack scratch buffers
        integer(int32), intent(inout) :: tmp_rec_stack(3, n_points)
            !! KD-tree build recursion-stack scratch buffer
        integer(int32), intent(inout) :: tmp_dim_order(dim)
            !! KD-tree dimension-splitting order scratch buffer
        real(real64),   intent(inout) :: tmp_val_buf(n_points)
            !! KD-tree/sort scratch buffer, reused across steps

        real(real64),   intent(out) :: sphere_radii(n_points)
            !! Per-point adaptive neighborhood radius
        real(real64),   intent(out) :: densities(n_points)
            !! Per-point local density estimate
        real(real64),   intent(out) :: gap_values(n_points)
            !! Per-point spectral gap at the kept neighborhood size
        real(real64),   intent(out) :: normal_errors(n_points)
            !! Mean squared residual of the best-kept neighborhood off its tangent subspace
        real(real64),   intent(out) :: stability_values(n_points)
            !! Subspace similarity of the best-kept neighborhood's tangent basis vs. its growth history
        integer(int32), intent(out) :: k_selected(n_points)
            !! Neighborhood size (excluding the point itself) of the best-kept neighborhood
        logical,        intent(out) :: growth_stopped_complex(n_points)
            !! .true. if growth stopped early because the tangent basis became unstable
        logical,        intent(out) :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors (Step 5)
        real(real64),   intent(inout) :: tangent_bases(dim, manifold_dim, n_points)
            !! On entry: previous outer iteration's tangent bases (zero on the first call),
            !! used to warm-start the tangent-stability check for the smallest neighborhood.
        real(real64),   intent(out) :: tangent_scales(manifold_dim, n_points)
            !! Per-point extent along each tangent direction
        integer(int32), intent(out) :: labels(n_points)
            !! BFS connected-overlap-region label; 0 if not in any intersection
        integer(int32), intent(out) :: primary_anchor_ids(n_points), secondary_anchor_ids(n_points)
            !! Point index of the highest/second-highest inverse-variance-weighted anchor
            !! covering this point (Step 10); 0 if not covered by that many anchors
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        ! Local variables
        logical, allocatable :: point_in_anchor_mask(:,:)  ! (n_points, n_anchor) membership matrix
        integer(int32), allocatable :: anchor_to_point(:)
        integer(int32) :: anchor_count(n_points)      ! how many anchors cover each point
        integer(int32) :: n_anchor

        call set_ok(ierr)

        ! STEPS 0-3: KD-tree + adaptive neighborhood growth
        call grow_adaptive_neighborhoods(work_coords, n_points, dim, manifold_dim, k_min, &
                                         g_threshold, stability_threshold, scale_factor, &
                                         tmp_kd_indices, tmp_workspace, tmp_val_buf, tmp_perm, &
                                         tmp_l_stack, tmp_r_stack, tmp_rec_stack, tmp_dim_order, &
                                         lwork, tangent_bases, &
                                         sphere_radii, densities, gap_values, normal_errors, &
                                         stability_values, k_selected, growth_stopped_complex, &
                                         tangent_scales, ierr)
        if (.not. is_ok(ierr)) return

        ! STEPS 4-5b: density sort + greedy anchor selection + orphan absorption
        call construct_atlas(work_coords, n_points, dim, o_min, o_max, densities, &
                             tmp_perm, tmp_l_stack, tmp_r_stack, sphere_radii, is_anchor_mask)

        ! STEP 6: per-anchor SVD (center, tangent basis, extent, fit quality)
        call compute_anchor_svd(work_coords, n_points, dim, manifold_dim, is_anchor_mask, &
                                  sphere_radii, lwork, &
                                  anchor_centers, tangent_bases, tangent_scales, normal_errors)

        ! STEPS 6.5-9: anchor mapping + membership matrix + intersection-graph BFS labeling
        call build_intersection_graph_alloc(work_coords, n_points, dim, is_anchor_mask, sphere_radii, &
                                            tmp_workspace, tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack, &
                                            point_in_anchor_mask, anchor_to_point, anchor_count, labels, &
                                            n_anchor, ierr)
        if (.not. is_ok(ierr)) return
        if (n_anchor == 0) return

        ! STEP 10: inverse-variance-weighted stitching
        call stitch_points(work_coords, n_points, dim, manifold_dim, n_anchor, &
                           point_in_anchor_mask, anchor_to_point, anchor_count, anchor_centers, &
                           tangent_bases, normal_errors, skeleton_coords, &
                           primary_anchor_ids, secondary_anchor_ids)

    end subroutine lomanle_pass

    !> STEP 0-3: build the KD-tree, then grow an adaptive neighborhood for every
    !| point via candidate growth, keeping the best-scoring size seen (not the
    !| last one evaluated). See misc/smoothing_experiments.md section 4 for the
    !| full rationale behind the quality score and the tangent-stability check.
    pure subroutine grow_adaptive_neighborhoods(work_coords, n_points, dim, manifold_dim, k_min, &
                                           g_threshold, stability_threshold, scale_factor, &
                                           tmp_kd_indices, tmp_workspace, tmp_val_buf, tmp_perm, &
                                           tmp_l_stack, tmp_r_stack, tmp_rec_stack, tmp_dim_order, &
                                           lwork, tangent_bases, &
                                           sphere_radii, densities, gap_values, normal_errors, &
                                           stability_values, k_selected, growth_stopped_complex, &
                                           tangent_scales, ierr)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: k_min
            !! Minimum neighborhood size to start adaptive growth from
        integer(int32), intent(in) :: lwork
            !! Size of the LAPACK dsyev scratch workspace
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        real(real64),   intent(in) :: g_threshold
            !! Spectral-gap threshold that stops neighborhood growth
        real(real64),   intent(in) :: stability_threshold
            !! Minimum tangent-basis stability to keep growing a neighborhood
        real(real64),   intent(in) :: scale_factor
            !! Cap on sphere radius as a multiple of the point's local scale

        integer(int32), intent(inout) :: tmp_kd_indices(n_points)
            !! KD-tree scratch buffer, see kd_tree module
        integer(int32), intent(inout) :: tmp_workspace(n_points)
            !! KD-tree/sort scratch buffer, reused during neighbor sorting
        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer, see f42_utils sort_array
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! KD-tree/sort_array recursion-stack scratch buffers
        integer(int32), intent(inout) :: tmp_rec_stack(3, n_points)
            !! KD-tree build recursion-stack scratch buffer
        integer(int32), intent(inout) :: tmp_dim_order(dim)
            !! KD-tree dimension-splitting order scratch buffer
        real(real64),   intent(inout) :: tmp_val_buf(n_points)
            !! KD-tree/sort scratch buffer, reused during neighbor sorting

        real(real64),   intent(inout) :: tangent_bases(dim, manifold_dim, n_points)
            !! On entry: previous outer iteration's tangent bases (zero on the first call),
            !! used to warm-start the tangent-stability check for the smallest neighborhood.
        real(real64),   intent(out) :: sphere_radii(n_points)
            !! Per-point adaptive neighborhood radius (kept-neighborhood max distance)
        real(real64),   intent(out) :: densities(n_points)
            !! Per-point local density estimate
        real(real64),   intent(out) :: gap_values(n_points)
            !! Per-point spectral gap at the kept neighborhood size
        real(real64),   intent(out) :: normal_errors(n_points)
            !! Mean squared residual of the kept neighborhood off its tangent subspace
        real(real64),   intent(out) :: stability_values(n_points)
            !! Tangent-basis stability of the kept neighborhood
        integer(int32), intent(out) :: k_selected(n_points)
            !! Neighborhood size (excluding the point itself) of the kept neighborhood
        logical,        intent(out) :: growth_stopped_complex(n_points)
            !! .true. if growth stopped early because the tangent basis became unstable
        real(real64),   intent(out) :: tangent_scales(manifold_dim, n_points)
            !! Per-point extent along each tangent direction
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        ! Local variables
        integer(int32) :: i, k_limit
        integer(int32) :: ierr_per_point(n_points)

        call set_ok(ierr)
        k_limit = n_points / 4
        do i = 1, dim ; tmp_dim_order(i) = i ; end do
        tangent_bases = 0.0_real64
        tangent_scales = 0.0_real64

        ! STEP 0: Build KD-Tree
        call build_kd_index(work_coords, dim, n_points, tmp_kd_indices, tmp_dim_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack, tmp_rec_stack, ierr)
        if (.not. is_ok(ierr)) return

        ! --- STEPS 1, 2 and 3: Adaptive Neighborhoods via Candidate Growth ---
        ! Each point's neighborhood growth depends only on shared read-only
        ! inputs and writes to its own slice of the output arrays, so the
        ! iterations are independent; every KD-NN/sort/LAPACK scratch buffer
        ! is private per-point here (tmp_kd_indices/tmp_dim_order, filled once by
        ! Step 0 above, are the only ones actually shared/read-only below).
        do concurrent (i = 1:n_points) shared(work_coords, n_points, dim, manifold_dim, k_min, &
                                              g_threshold, stability_threshold, scale_factor, k_limit, &
                                              tmp_kd_indices, tmp_dim_order, lwork, tangent_bases, &
                                              sphere_radii, densities, gap_values, normal_errors, &
                                              stability_values, k_selected, growth_stopped_complex, &
                                              tangent_scales, ierr_per_point)
            block
                integer(int32) :: n_loc_i(n_points), workspace_i(n_points), perm_i(n_points)
                integer(int32) :: l_stack_i(n_points), r_stack_i(n_points)
                real(real64)   :: d_loc_i(n_points), val_buf_i(n_points), work_lapack_i(lwork)
                integer(int32) :: ierr_i

                call grow_one_point_neighborhood(i, work_coords, n_points, dim, manifold_dim, k_min, &
                                                 g_threshold, stability_threshold, scale_factor, k_limit, &
                                                 tmp_kd_indices, tmp_dim_order, lwork, work_lapack_i, n_loc_i, d_loc_i, &
                                                 workspace_i, val_buf_i, perm_i, l_stack_i, r_stack_i, &
                                                 tangent_bases(:,:,i), sphere_radii(i), densities(i), &
                                                 gap_values(i), normal_errors(i), stability_values(i), &
                                                 k_selected(i), growth_stopped_complex(i), &
                                                 tangent_scales(:,i), ierr_i)
                ierr_per_point(i) = ierr_i
            end block
        end do

        ! Surface the first point-level failure, if any (each point's own
        ! kd_knn_query failure is already handled gracefully inside
        ! grow_one_point_neighborhood; this only matters for a genuinely
        ! fatal, dataset-wide condition).
        do i = 1, n_points
            if (.not. is_ok(ierr_per_point(i))) then
                ierr = ierr_per_point(i)
                exit
            end if
        end do

    end subroutine grow_adaptive_neighborhoods

    !> STEPS 1-3 for exactly one point: grows its neighborhood from k_min via
    !| candidate growth, scoring each candidate size with a quality score
    !| (spectral gap + tangent stability - normal error - radius), and keeps
    !| the best-scoring size found so far rather than the last one evaluated,
    !| so a single bad spectral gap reading doesn't force growth to continue.
    !| Growth stops once the tangent basis is unstable for several
    !| consecutive steps. See misc/smoothing_experiments.md section 4 for
    !| the full rationale.
    pure subroutine grow_one_point_neighborhood(i, work_coords, n_points, dim, manifold_dim, k_min, &
                                           g_threshold, stability_threshold, scale_factor, k_limit, &
                                           tmp_kd_indices, tmp_dim_order, lwork, tmp_work_lapack, tmp_n_loc, tmp_d_loc, &
                                           tmp_workspace, tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack, &
                                           point_tangent_basis, sphere_radius_i, density_i, gap_i, &
                                           point_normal_error, point_stability, k_selected_i, &
                                           stopped_complex, point_tangent_scales, ierr)

        integer(int32), intent(in) :: i
            !! Index of the point to grow a neighborhood for
        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: k_min
            !! Minimum neighborhood size to start adaptive growth from
        integer(int32), intent(in) :: k_limit
            !! Hard cap on neighborhood size (n_points / 4)
        integer(int32), intent(in) :: lwork
            !! Size of the tmp_work_lapack scratch array (see dsyev)
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        real(real64),   intent(in) :: g_threshold
            !! Spectral-gap threshold that stops neighborhood growth
        real(real64),   intent(in) :: stability_threshold
            !! Minimum tangent-basis stability to keep growing a neighborhood
        real(real64),   intent(in) :: scale_factor
            !! Cap on sphere radius as a multiple of the point's local scale

        integer(int32), intent(in) :: tmp_kd_indices(n_points)
            !! KD-tree scratch buffer, see kd_tree module
        integer(int32), intent(in) :: tmp_dim_order(dim)
            !! KD-tree dimension-splitting order
        real(real64),   intent(inout) :: tmp_work_lapack(lwork)
            !! LAPACK dsyev scratch workspace
        integer(int32), intent(inout) :: tmp_n_loc(:)
            !! KD-NN query result: neighbor point indices
        real(real64),   intent(inout) :: tmp_d_loc(:)
            !! KD-NN query result: neighbor distances
        integer(int32), intent(inout) :: tmp_workspace(:)
            !! Sort scratch buffer, reused during neighbor sorting
        real(real64),   intent(inout) :: tmp_val_buf(:)
            !! Sort scratch buffer, reused during neighbor sorting
        integer(int32), intent(inout) :: tmp_perm(:)
            !! Sort-permutation scratch buffer, see f42_utils sort_array
        integer(int32), intent(inout) :: tmp_l_stack(:), tmp_r_stack(:)
            !! sort_array recursion-stack scratch buffers

        real(real64),   intent(inout) :: point_tangent_basis(dim, manifold_dim)
            !! On entry: previous outer iteration's tangent basis for this point
            !! (zero on the first call), used to warm-start the tangent-stability
            !! check for the smallest neighborhood.
        real(real64),   intent(out) :: sphere_radius_i
            !! Adaptive neighborhood radius (kept-neighborhood max distance)
        real(real64),   intent(out) :: density_i
            !! Local density estimate
        real(real64),   intent(out) :: gap_i
            !! Spectral gap at the kept neighborhood size
        real(real64),   intent(out) :: point_normal_error
            !! Mean squared residual of the kept neighborhood off its tangent subspace
        real(real64),   intent(out) :: point_stability
            !! Tangent-basis stability of the kept neighborhood
        integer(int32), intent(out) :: k_selected_i
            !! Neighborhood size (excluding the point itself) of the kept neighborhood
        logical,        intent(out) :: stopped_complex
            !! .true. if growth stopped early because the tangent basis became unstable
        real(real64),   intent(out) :: point_tangent_scales(manifold_dim)
            !! Extent along each tangent direction
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        ! Local variables
        integer(int32) :: j, row, col, info, k_curr, d_idx, base_idx
        integer(int32) :: k_query, self_pos, best_k, patience, consecutive_unstable
        logical        :: have_previous, first_step
        real(real64)   :: sigma_i, sigma2_i, rho_i, s_gap, dist_sq
        real(real64)   :: center(dim), cov(dim, dim), w_eig(dim), p_diff(dim), projection
        real(real64)   :: previous_tangent(dim, manifold_dim)
        real(real64)   :: cov_small(manifold_dim, manifold_dim), gram_small(manifold_dim, manifold_dim)
        real(real64)   :: w_eig_small(manifold_dim)
        real(real64)   :: stability_i, normal_error_i, local_scale_i, quality_i, best_quality
        real(real64)   :: noise_ratio_i
        integer(int32), parameter :: max_instability_patience = 5
            !! Safety cap on how many consecutive unstable readings we'll tolerate
            !! before giving up on a point, however noisy its neighborhood looks.
            !! Not a per-dataset tuning knob: the actual patience used is derived
            !! from each point's own noise_ratio_i, so it self-calibrates instead
            !! of needing to be re-tuned per dataset.

        ! Warm-start: this point's tangent basis from the previous outer iteration
        ! (zero on the very first call) seeds the stability check for the smallest
        ! neighborhood; later growth steps compare against the immediately
        ! preceding step instead (see `previous_tangent` update at the loop's end).
        previous_tangent = point_tangent_basis
        have_previous = any(abs(previous_tangent) > 1.0e-12_real64)
        first_step = .true.

        k_curr = k_min
        best_quality = -huge(1.0_real64)
        best_k = 0
        local_scale_i = 0.0_real64
        stopped_complex = .false.
        consecutive_unstable = 0
        patience = 1 ! recomputed from noise_ratio_i once the k_min neighborhood is known

        adaptive_k: do
            ! Query one extra neighbor so the query point itself (always
            ! returned at distance 0) can be dropped without losing k_curr
            ! real neighbors.
            k_query = min(k_curr + 1, size(tmp_n_loc))
            call kd_knn_query(work_coords, tmp_kd_indices, dim, n_points, tmp_dim_order, &
                              work_coords(:,i), k_query, tmp_n_loc(1:k_query), tmp_d_loc(1:k_query), ierr)
            if (.not. is_ok(ierr)) exit adaptive_k

            self_pos = 0
            do j = 1, k_query
                if (tmp_n_loc(j) == i) then
                    self_pos = j
                    exit
                end if
            end do
            if (self_pos > 0 .and. self_pos < k_query) then
                tmp_n_loc(self_pos) = tmp_n_loc(k_query)
                tmp_d_loc(self_pos) = tmp_d_loc(k_query)
            end if
            k_curr = k_query - 1
            if (k_curr < manifold_dim + 1) exit adaptive_k ! not enough points for a tangent

            ! Sort ascending so the radius, local scale and jump checks are meaningful.
            ! sort_array leaves its `array` argument untouched (intent(in)) and only
            ! permutes `tmp_perm` -- which must start as the identity 1..k_curr -- so we
            ! sort a throwaway identity permutation and apply it to tmp_d_loc/tmp_n_loc
            ! ourselves via tmp_workspace/tmp_val_buf (both unused elsewhere during Steps 1-3).
            do j = 1, k_curr
                tmp_perm(j) = j
            end do
            call sort_array(tmp_d_loc(1:k_curr), tmp_perm(1:k_curr), tmp_l_stack(1:k_curr), tmp_r_stack(1:k_curr))
            tmp_workspace(1:k_curr) = tmp_n_loc(1:k_curr)
            tmp_val_buf(1:k_curr) = tmp_d_loc(1:k_curr)
            do j = 1, k_curr
                tmp_n_loc(j) = tmp_workspace(tmp_perm(j))
                tmp_d_loc(j) = tmp_val_buf(tmp_perm(j))
            end do

            if (k_curr == k_min) then
                if (mod(k_min, 2) == 1) then
                    local_scale_i = tmp_d_loc((k_min + 1) / 2)
                else
                    local_scale_i = 0.5_real64 * (tmp_d_loc(k_min / 2) + tmp_d_loc(k_min / 2 + 1))
                end if
            end if

            sigma_i = tmp_d_loc(k_curr) ! true max distance, now that neighbors are sorted
            center = sum(work_coords(:, tmp_n_loc(1:k_curr)), dim=2) / real(k_curr, real64)

            ! Unconditional floor: record this neighborhood as the current
            ! best-so-far *before* attempting the eigendecomposition below,
            ! so sphere_radius_i/density_i/k_selected_i never fall below k_min
            ! even if dsyev never succeeds for this point (e.g. a degenerate
            ! covariance from near-duplicate points).
            if (best_k == 0) then
                best_k = k_curr
                sphere_radius_i = sigma_i
                gap_i = 0.0_real64
                point_normal_error = 0.0_real64
                point_stability = 0.0_real64
                rho_i = 0.0_real64
                sigma2_i = sigma_i**2
                if (sigma_i > 1.0e-12_real64) then
                    do j = 1, k_curr
                        rho_i = rho_i + exp(-(tmp_d_loc(j)**2) / (2.0_real64 * sigma2_i))
                    end do
                else
                    rho_i = real(k_curr, real64)
                end if
                density_i = rho_i / (max(1.0e-12_real64, sigma_i)**manifold_dim)
            end if

            cov = 0.0_real64
            do j = 1, k_curr
                p_diff = work_coords(:, tmp_n_loc(j)) - center
                do col = 1, dim
                    do row = 1, dim
                        cov(row, col) = cov(row, col) + p_diff(row) * p_diff(col)
                    end do
                end do
            end do
            cov = cov / max(1.0_real64, real(k_curr - 1, real64))

            call dsyev('V', 'U', dim, cov, dim, w_eig, tmp_work_lapack, lwork, info)
            if (info /= 0) exit adaptive_k

            d_idx = dim - manifold_dim + 1
            if (w_eig(d_idx - 1) > 1.0e-12_real64) then
                s_gap = sqrt(w_eig(d_idx)) / sqrt(w_eig(d_idx - 1))
            else
                s_gap = g_threshold + 1.0_real64
            end if

            ! Normal reconstruction error: mean squared distance of neighbors
            ! from the manifold_dim-dimensional tangent subspace through `center`
            normal_error_i = 0.0_real64
            do j = 1, k_curr
                p_diff = work_coords(:, tmp_n_loc(j)) - center
                dist_sq = dot_product(p_diff, p_diff)
                do base_idx = 1, manifold_dim
                    projection = dot_product(p_diff, cov(:, dim - manifold_dim + base_idx))
                    dist_sq = dist_sq - projection**2
                end do
                normal_error_i = normal_error_i + max(0.0_real64, dist_sq)
            end do
            normal_error_i = normal_error_i / real(k_curr, real64)

            if (k_curr == k_min) then
                ! How noisy this point's neighborhood looks, relative to its own
                ! scale: dimensionless (both sides are lengths), so it self-
                ! calibrates per point/dataset instead of needing a fixed constant.
                ! Clean, curve-like data -> ~0 -> patience=1 (react immediately).
                ! Noisy, blob-like data -> ~1 or more -> patience grows toward
                ! max_instability_patience, giving the small-sample tangent
                ! estimate room to settle before we trust a bad reading.
                noise_ratio_i = sqrt(normal_error_i) / max(local_scale_i, 1.0e-12_real64)
                patience = 1 + nint(min(noise_ratio_i, 1.0_real64) * real(max_instability_patience - 1, real64))
            end if

            ! Tangent stability: vs. the warm-started previous outer iteration on
            ! the first growth step, then step-to-step against the last evaluated
            ! size. manifold_dim==1 reduces to a plain dot product; manifold_dim>1
            ! uses the smallest singular value of previous_tangent^T * current_tangent,
            ! computed as sqrt(smallest eigenvalue of its Gram matrix) so the
            ! existing dsyev interface (and its tmp_work_lapack sizing) can be reused.
            if (first_step .and. .not. have_previous) then
                stability_i = 1.0_real64
            else if (manifold_dim == 1) then
                stability_i = abs(dot_product(previous_tangent(:,1), cov(:, dim)))
            else
                cov_small = matmul(transpose(previous_tangent), cov(:, dim - manifold_dim + 1:dim))
                gram_small = matmul(transpose(cov_small), cov_small)
                call dsyev('N', 'U', manifold_dim, gram_small, manifold_dim, w_eig_small, tmp_work_lapack, lwork, info)
                if (info == 0) then
                    stability_i = sqrt(max(0.0_real64, w_eig_small(1)))
                else
                    stability_i = 0.0_real64
                end if
            end if

            ! Density (unchanged formula; still used for anchor seeding later)
            rho_i = 0.0_real64
            sigma2_i = sigma_i**2
            if (sigma_i > 1.0e-12_real64) then
                do j = 1, k_curr
                    rho_i = rho_i + exp(-(tmp_d_loc(j)**2) / (2.0_real64 * sigma2_i))
                end do
            else
                rho_i = real(k_curr, real64)
            end if

            ! Quality score: reward gap + stability, penalise normal error and an
            ! oversized radius relative to the local scale established at k_min.
            quality_i = min(s_gap / g_threshold, 2.0_real64) &
                      + stability_i &
                      - normal_error_i / max(local_scale_i**2, 1.0e-12_real64) &
                      - sigma_i / max(local_scale_i, 1.0e-12_real64)

            if (quality_i > best_quality) then
                best_quality = quality_i
                best_k = k_curr
                sphere_radius_i = sigma_i
                gap_i = s_gap
                point_normal_error = normal_error_i
                point_stability = stability_i
                density_i = rho_i / (max(1.0e-12_real64, sigma_i)**manifold_dim)
                do base_idx = 1, manifold_dim
                    point_tangent_basis(:, base_idx) = cov(:, dim - manifold_dim + base_idx)
                    point_tangent_scales(base_idx) = sqrt(max(0.0_real64, w_eig(dim - manifold_dim + base_idx)))
                end do
            end if

            ! Stopping conditions: an unstable tangent halts growth once it has
            ! been confirmed `patience` times in a row (a single bad reading can
            ! just be small-sample noise settling down; a real branch/curvature
            ! transition keeps failing as we grow further into it). A stable
            ! reading resets the counter. Otherwise growth continues only while
            ! the gap is still insufficient.
            if (stability_i < stability_threshold) then
                consecutive_unstable = consecutive_unstable + 1
                if (consecutive_unstable >= patience) then
                    stopped_complex = .true.
                    exit adaptive_k
                end if
            else
                consecutive_unstable = 0
                if (s_gap >= g_threshold) exit adaptive_k
            end if
            if (k_curr >= k_limit) exit adaptive_k
            if (sigma_i > scale_factor * max(local_scale_i, 1.0e-12_real64)) exit adaptive_k

            previous_tangent = cov(:, dim - manifold_dim + 1:dim)
            first_step = .false.
            k_curr = nint(k_curr * 1.25)
            if (k_curr > size(tmp_n_loc) - 1) k_curr = size(tmp_n_loc) - 1 ! Guard against buffer overflow
        end do adaptive_k

        k_selected_i = best_k
        if (best_k == 0) then
            ! Degenerate case (e.g. n_points too small for manifold_dim): keep
            ! outputs well-defined even though no neighborhood could be scored.
            sphere_radius_i = 0.0_real64
            gap_i = 0.0_real64
            point_normal_error = 0.0_real64
            point_stability = 0.0_real64
            density_i = 0.0_real64
        end if

    end subroutine grow_one_point_neighborhood

    !> STEP 4-5b: sort points by density, then greedily grow a set of anchor
    !| spheres subject to the o_min/o_max overlap-ratio constraint (Step 5),
    !| finally growing the nearest anchor's radius to absorb any leftover
    !| uncovered "orphan" points (Step 5b).
    pure subroutine construct_atlas(work_coords, n_points, dim, o_min, o_max, densities, &
                               tmp_perm, tmp_l_stack, tmp_r_stack, sphere_radii, is_anchor_mask)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        real(real64),   intent(in) :: o_min
            !! Minimum required anchor-sphere overlap ratio
        real(real64),   intent(in) :: o_max
            !! Maximum allowed anchor-sphere overlap ratio
        real(real64),   intent(in) :: densities(n_points)
            !! Per-point local density estimate, used to order/rank candidates

        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer; ends up holding the
            !! descending-by-density point order used to pick candidates
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! sort_array recursion-stack scratch buffers

        real(real64),   intent(inout) :: sphere_radii(n_points)
            !! In: per-point adaptive radius from Steps 1-3. Out: orphan-adjacent
            !! anchor radii may have been grown (Step 5b) to absorb an orphan.
        logical,        intent(out) :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors

        ! Local variables
        integer(int32) :: n_anchor
        logical        :: is_covered_mask(n_points)

        is_anchor_mask  = .false.

        ! STEP 4: sort points by density, descending
        call sort_points_by_density(n_points, densities, tmp_perm, tmp_l_stack, tmp_r_stack)

        ! STEP 5: greedy anchor selection under the o_min/o_max overlap constraint
        call select_atlas_anchors(work_coords, n_points, dim, o_min, o_max, densities, &
                                  tmp_perm, sphere_radii, is_anchor_mask, is_covered_mask, n_anchor)

        ! STEP 5b: absorb leftover uncovered points into the nearest anchor
        call absorb_orphans(work_coords, n_points, dim, is_anchor_mask, sphere_radii, is_covered_mask)

    end subroutine construct_atlas

    !> STEP 4: sorts point indices into tmp_perm, descending by density, so
    !| select_atlas_anchors can consider the densest remaining candidate
    !| first.
    pure subroutine sort_points_by_density(n_points, densities, tmp_perm, tmp_l_stack, tmp_r_stack)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        real(real64),   intent(in) :: densities(n_points)
            !! Per-point local density estimate

        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer; ends up holding the
            !! descending-by-density point order
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! sort_array recursion-stack scratch buffers

        integer(int32) :: i, tmp_idx

        do i = 1, n_points ; tmp_perm(i) = i ; end do
        call sort_array(densities, tmp_perm, tmp_l_stack, tmp_r_stack)
        do i = 1, n_points / 2 ! Reverse to descending
            tmp_idx = tmp_perm(i) ; tmp_perm(i) = tmp_perm(n_points - i + 1) ; tmp_perm(n_points - i + 1) = tmp_idx
        end do

    end subroutine sort_points_by_density

    !> STEP 5: greedily selects anchor spheres in density order, subject to
    !| the o_min/o_max overlap-ratio constraint; opens a new component from
    !| the highest-density uncovered point whenever no candidate satisfies
    !| the constraint.
    pure subroutine select_atlas_anchors(work_coords, n_points, dim, o_min, o_max, densities, &
                                    tmp_perm, sphere_radii, is_anchor_mask, is_covered_mask, n_anchor)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        real(real64),   intent(in) :: o_min
            !! Minimum required anchor-sphere overlap ratio
        real(real64),   intent(in) :: o_max
            !! Maximum allowed anchor-sphere overlap ratio
        real(real64),   intent(in) :: densities(n_points)
            !! Per-point local density estimate, used to order/rank candidates
        integer(int32), intent(in) :: tmp_perm(n_points)
            !! Point indices in descending-density order, see sort_points_by_density
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-point adaptive radius from Steps 1-3

        logical,        intent(out) :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors
        logical,        intent(out) :: is_covered_mask(n_points)
            !! .true. for points inside at least one selected anchor's sphere
        integer(int32), intent(out) :: n_anchor
            !! Number of anchors selected

        integer(int32) :: i, k_idx, m_idx, num_overlap, total_in_sphere
        integer(int32) :: best_candidate, n_uncovered
        logical        :: found_candidate
        real(real64)   :: current_ratio, dist_sq, best_density

        is_covered_mask = .false.
        n_anchor = 0

        ! Main loop: continue until all (or nearly all) points are covered
        atlas_construction: do
            ! Count uncovered points
            n_uncovered = count(.not. is_covered_mask)
            if (n_uncovered == 0) exit atlas_construction

            ! Find best candidate that satisfies o_min <= overlap <= o_max
            best_candidate = -1
            best_density = -1.0_real64
            found_candidate = .false.

            do k_idx = 1, n_points
                i = tmp_perm(k_idx)  ! Iterate in density order

                ! Skip if already anchor or if this point is covered
                if (is_anchor_mask(i) .or. is_covered_mask(i)) cycle

                ! Calculate overlap with current atlas
                num_overlap = 0
                total_in_sphere = 0
                do m_idx = 1, n_points
                    dist_sq = sum((work_coords(:, i) - work_coords(:, m_idx))**2)
                    if (dist_sq <= sphere_radii(i)**2) then
                        total_in_sphere = total_in_sphere + 1
                        if (is_covered_mask(m_idx)) num_overlap = num_overlap + 1
                    end if
                end do

                current_ratio = 0.0_real64
                if (total_in_sphere > 0) current_ratio = real(num_overlap, real64) / real(total_in_sphere, real64)

                ! Check if this candidate satisfies connectivity constraints
                ! For the first anchor (n_anchor == 0), we don't check o_min
                if (n_anchor == 0) then
                    if (.not. found_candidate) then
                        best_candidate = i
                        best_density = densities(i)
                        found_candidate = .true.
                        exit  ! Take the first (highest density) as seed
                    end if
                else
                    ! For subsequent anchors, enforce o_min <= overlap <= o_max
                    if (current_ratio >= o_min .and. current_ratio <= o_max) then
                        if (densities(i) > best_density) then
                            best_candidate = i
                            best_density = densities(i)
                            found_candidate = .true.
                        end if
                    end if
                end if
            end do

            ! If we found a valid candidate, add it as anchor
            if (found_candidate .and. best_candidate > 0) then
                is_anchor_mask(best_candidate) = .true.
                n_anchor = n_anchor + 1

                ! Mark all points in this sphere as covered
                do m_idx = 1, n_points
                    dist_sq = sum((work_coords(:, best_candidate) - work_coords(:, m_idx))**2)
                    if (dist_sq <= sphere_radii(best_candidate)**2) then
                        is_covered_mask(m_idx) = .true.
                    end if
                end do
            else
                ! No valid candidate found - try to start new component with highest density uncovered point
                best_candidate = -1
                best_density = -1.0_real64
                do k_idx = 1, n_points
                    i = tmp_perm(k_idx)
                    if (.not. is_covered_mask(i) .and. .not. is_anchor_mask(i)) then
                        if (densities(i) > best_density) then
                            best_candidate = i
                            best_density = densities(i)
                        end if
                    end if
                end do

                if (best_candidate > 0) then
                    is_anchor_mask(best_candidate) = .true.
                    n_anchor = n_anchor + 1
                    do m_idx = 1, n_points
                        dist_sq = sum((work_coords(:, best_candidate) - work_coords(:, m_idx))**2)
                        if (dist_sq <= sphere_radii(best_candidate)**2) is_covered_mask(m_idx) = .true.
                    end do
                else
                    ! Should rarely happen: no uncovered point left to seed a new
                    ! component from, so there is nothing further this loop can do.
                    exit atlas_construction
                end if
            end if
        end do atlas_construction

    end subroutine select_atlas_anchors

    !> STEP 5b: for every point still uncovered after select_atlas_anchors,
    !| grows the nearest anchor's sphere just enough to absorb it.
    pure subroutine absorb_orphans(work_coords, n_points, dim, is_anchor_mask, sphere_radii, is_covered_mask)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        logical,        intent(in) :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors

        real(real64),   intent(inout) :: sphere_radii(n_points)
            !! In: per-point adaptive radius from Steps 1-3. Out: orphan-adjacent
            !! anchor radii may have been grown to absorb an orphan.
        logical,        intent(inout) :: is_covered_mask(n_points)
            !! In: coverage after select_atlas_anchors. Out: orphans marked covered.

        integer(int32) :: i, j, nearest_anchor, n_orphans
        real(real64)   :: dist_sq, min_dist, new_radius

        ! Only grow spheres for remaining orphans if there are very few
        n_orphans = count(.not. is_covered_mask)
        if (n_orphans > 0) then
            do i = 1, n_points
                if (is_covered_mask(i) .or. is_anchor_mask(i)) cycle

                ! Find nearest anchor
                min_dist = huge(1.0_real64)
                nearest_anchor = -1
                do j = 1, n_points
                    if (.not. is_anchor_mask(j)) cycle
                    dist_sq = sum((work_coords(:, i) - work_coords(:, j))**2)
                    if (dist_sq < min_dist) then
                        min_dist = dist_sq
                        nearest_anchor = j
                    end if
                end do

                ! Grow that anchor's sphere to include the orphan
                if (nearest_anchor > 0) then
                    new_radius = sqrt(min_dist)
                    if (new_radius > sphere_radii(nearest_anchor)) then
                        sphere_radii(nearest_anchor) = new_radius
                    end if
                    is_covered_mask(i) = .true.
                end if
            end do
        end if

    end subroutine absorb_orphans

    !> STEP 6: for each selected anchor, computes its center, tangent basis,
    !| extent, and fit-quality (normal_error) over its FINAL sphere membership
    !| (radii may have changed since Steps 1-3 during atlas/orphan handling).
    !| Non-anchor points keep whatever Steps 1-3 already wrote for them.
    pure subroutine compute_anchor_svd(work_coords, n_points, dim, manifold_dim, is_anchor_mask, &
                                  sphere_radii, lwork, &
                                  anchor_centers, tangent_bases, tangent_scales, normal_errors)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: lwork
            !! Size of the tmp_work_lapack scratch array (see dsyev)
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        logical,        intent(in) :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-point (or, for anchors, per-anchor) sphere radius

        real(real64),   intent(out)   :: anchor_centers(dim, n_points)
            !! Default: anchor centre = anchor position; overwritten below for
            !! every valid (enough-points-for-SVD) anchor.
        real(real64),   intent(inout) :: tangent_bases(dim, manifold_dim, n_points)
            !! Only entries where is_anchor_mask(i) is .true. are overwritten here.
        real(real64),   intent(inout) :: tangent_scales(manifold_dim, n_points)
            !! Only entries where is_anchor_mask(i) is .true. are overwritten here.
        real(real64),   intent(inout) :: normal_errors(n_points)
            !! Only entries where is_anchor_mask(i) is .true. are overwritten here,
            !! with the fuller-sphere fit quality Step 10 uses for weighting.

        ! Local variables
        integer(int32) :: i

        anchor_centers = work_coords

        ! --- STEP 6: Compute SVD for Anchor Spheres ---
        ! Each anchor's SVD depends only on shared read-only inputs and writes
        ! to its own slice of the output arrays, so the iterations are
        ! independent; tmp_n_loc/tmp_work_lapack are private per-anchor scratch (each
        ! anchor's own sphere is gathered fresh, never shared across anchors).
        do concurrent (i = 1:n_points) shared(work_coords, n_points, dim, manifold_dim, lwork, &
                                              sphere_radii, is_anchor_mask, anchor_centers, &
                                              tangent_bases, tangent_scales, normal_errors)
            if (.not. is_anchor_mask(i)) cycle
            block
                integer(int32) :: j, m_idx, k_curr, col, row, info, base_idx
                real(real64)   :: dist_sq, center(dim), cov(dim, dim), p_diff(dim), w_eig(dim)
                real(real64)   :: projection, min_proj, max_proj
                integer(int32) :: tmp_n_loc(n_points)
                real(real64)   :: tmp_work_lapack(lwork)

                ! Gather all points within the final sphere radius
                k_curr = 0
                do m_idx = 1, n_points
                    dist_sq = sum((work_coords(:, i) - work_coords(:, m_idx))**2)
                    if (dist_sq <= sphere_radii(i)**2) then
                        k_curr = k_curr + 1
                        if (k_curr > size(tmp_n_loc)) exit  ! Buffer overflow protection
                        tmp_n_loc(k_curr) = m_idx
                    end if
                end do

                if (k_curr >= manifold_dim + 1) then  ! Need enough points for SVD
                    ! Compute center
                    center = sum(work_coords(:, tmp_n_loc(1:k_curr)), dim=2) / real(k_curr, real64)
                    anchor_centers(:, i) = center   ! store for Step 10 stitching

                    ! Compute covariance matrix
                    cov = 0.0_real64
                    do j = 1, k_curr
                        p_diff = work_coords(:, tmp_n_loc(j)) - center
                        do col = 1, dim
                            do row = 1, dim
                                cov(row, col) = cov(row, col) + p_diff(row) * p_diff(col)
                            end do
                        end do
                    end do
                    cov = cov / max(1.0_real64, real(k_curr - 1, real64))

                    ! Compute SVD
                    call dsyev('V', 'U', dim, cov, dim, w_eig, tmp_work_lapack, lwork, info)
                    if (info == 0) then
                        ! Update tangent bases
                        do base_idx = 1, manifold_dim
                            tangent_bases(:, base_idx, i) = cov(:, dim - manifold_dim + base_idx)

                            ! Calculate scale as actual extent along this PC direction
                            min_proj = huge(1.0_real64)
                            max_proj = -huge(1.0_real64)
                            do j = 1, k_curr
                                p_diff = work_coords(:, tmp_n_loc(j)) - center
                                projection = dot_product(p_diff, tangent_bases(:, base_idx, i))
                                min_proj = min(min_proj, projection)
                                max_proj = max(max_proj, projection)
                            end do
                            tangent_scales(base_idx, i) = (max_proj - min_proj) / 2.0_real64
                        end do

                        ! Anchor-level fit quality: mean squared distance of the FULL
                        ! sphere's members from this tangent subspace (not the smaller
                        ! Steps 1-3 neighborhood normal_errors(i) was originally computed
                        ! from). Step 10 uses this to weight anchors by how well their
                        ! tangent plane actually fits, not just by how many points
                        ! happen to be nearby.
                        normal_errors(i) = 0.0_real64
                        do j = 1, k_curr
                            p_diff = work_coords(:, tmp_n_loc(j)) - center
                            dist_sq = dot_product(p_diff, p_diff)
                            do base_idx = 1, manifold_dim
                                projection = dot_product(p_diff, tangent_bases(:, base_idx, i))
                                dist_sq = dist_sq - projection**2
                            end do
                            normal_errors(i) = normal_errors(i) + max(0.0_real64, dist_sq)
                        end do
                        normal_errors(i) = normal_errors(i) / real(k_curr, real64)
                    end if
                end if
            end block
        end do

    end subroutine compute_anchor_svd

    !> STEP 7: fills the (n_points, n_anchor) membership matrix, where
    !| M(i, i_anc) = .true. iff point i falls inside anchor i_anc's sphere.
    !| This is the central data structure for intersection identification.
    pure subroutine build_membership_matrix(work_coords, n_points, dim, n_anchor, anchor_to_point, &
                                       sphere_radii, point_in_anchor_mask)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-anchor sphere radius

        logical, intent(out) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix

        integer(int32) :: i_anc, m_idx

        ! Each (m_idx, i_anc) entry depends only on shared read-only inputs
        ! and writes to its own matrix element, so the iterations are independent.
        point_in_anchor_mask = .false.
        do concurrent (i_anc = 1:n_anchor, m_idx = 1:n_points) &
                shared(work_coords, n_points, dim, anchor_to_point, sphere_radii, point_in_anchor_mask)
            block
                integer(int32) :: i
                real(real64)   :: dist_sq
                i = anchor_to_point(i_anc)
                dist_sq = sum((work_coords(:,m_idx) - work_coords(:,i))**2)
                if (dist_sq <= sphere_radii(i)**2) point_in_anchor_mask(m_idx, i_anc) = .true.
            end block
        end do

    end subroutine build_membership_matrix

    !> STEP 8, phase a: counts pairwise anchor intersections (an edge = two
    !| anchors that share at least one point), so the caller can allocate the
    !| edge-list arrays to the exact size before fill_anchor_intersection_edges
    !| fills them.
    pure subroutine count_anchor_intersection_edges(n_points, n_anchor, point_in_anchor_mask, n_edges)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        logical,        intent(in) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix

        integer(int32), intent(out) :: n_edges
            !! Number of anchor pairs that share at least one point

        integer(int32) :: i_anc, j_anc

        n_edges = 0
        do i_anc = 1, n_anchor
            do j_anc = i_anc + 1, n_anchor
                if (any(point_in_anchor_mask(:,i_anc) .and. point_in_anchor_mask(:,j_anc))) &
                    n_edges = n_edges + 1
            end do
        end do

    end subroutine count_anchor_intersection_edges

    !> STEP 8, phase b: fills the anchor-pair edge list (edge_anc1/edge_anc2)
    !| for every pair of anchors that share at least one point.
    pure subroutine fill_anchor_intersection_edges(n_points, n_anchor, point_in_anchor_mask, n_edges, &
                                              edge_anc1, edge_anc2)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        logical,        intent(in) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix
        integer(int32), intent(in) :: n_edges
            !! Number of anchor pairs that share at least one point

        integer(int32), intent(out) :: edge_anc1(n_edges), edge_anc2(n_edges)
            !! Anchor-pair edge list (compact anchor indices)

        integer(int32) :: i_anc, j_anc, e

        e = 0
        do i_anc = 1, n_anchor
            do j_anc = i_anc + 1, n_anchor
                if (any(point_in_anchor_mask(:,i_anc) .and. point_in_anchor_mask(:,j_anc))) then
                    e = e + 1
                    edge_anc1(e) = i_anc
                    edge_anc2(e) = j_anc
                end if
            end do
        end do

    end subroutine fill_anchor_intersection_edges

    !> STEP 8, phase c: counts, in a single pass, how many anchors cover each
    !| point (anchor_count, a real output the caller needs) and how many
    !| edges touch each point/anchor-pair (n_edges_of_pt/edge_pt_count, used
    !| only to size the CSR data arrays that follow).
    pure subroutine count_intersection_csr_sizes(n_points, n_anchor, n_edges, point_in_anchor_mask, &
                                            edge_anc1, edge_anc2, anchor_count, &
                                            n_edges_of_pt, edge_pt_count)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: n_edges
            !! Number of anchor pairs that share at least one point
        logical,        intent(in) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix
        integer(int32), intent(in) :: edge_anc1(n_edges), edge_anc2(n_edges)
            !! Anchor-pair edge list (compact anchor indices)

        integer(int32), intent(out) :: anchor_count(n_points)
            !! How many anchors cover each point
        integer(int32), intent(out) :: n_edges_of_pt(n_points)
            !! How many anchor-pair edges touch each point
        integer(int32), intent(out) :: edge_pt_count(n_edges)
            !! How many points touch each anchor-pair edge

        integer(int32) :: i, i_anc, e

        anchor_count  = 0
        n_edges_of_pt = 0
        edge_pt_count = 0
        do i_anc = 1, n_anchor
            do i = 1, n_points
                if (point_in_anchor_mask(i, i_anc)) anchor_count(i) = anchor_count(i) + 1
            end do
        end do
        do e = 1, n_edges
            do i = 1, n_points
                if (point_in_anchor_mask(i, edge_anc1(e)) .and. &
                    point_in_anchor_mask(i, edge_anc2(e))) then
                    n_edges_of_pt(i) = n_edges_of_pt(i) + 1
                    edge_pt_count(e) = edge_pt_count(e) + 1
                end if
            end do
        end do

    end subroutine count_intersection_csr_sizes

    !> STEP 6.5-8d: allocates the anchor index mapping, the point/anchor
    !| membership matrix, and the CSR adjacency structures over the pairwise
    !| anchor-intersection edges, sizing each buffer exactly once its size is
    !| known (edge/CSR sizes can't be known before the membership matrix is
    !| built). Then calls build_intersection_graph to fill the CSR lists and
    !| BFS-label the connected overlap regions, which needs no further
    !| allocation.
    subroutine build_intersection_graph_alloc(work_coords, n_points, dim, is_anchor_mask, sphere_radii, &
                                              tmp_workspace, tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack, &
                                              point_in_anchor_mask, anchor_to_point, anchor_count, &
                                              labels, n_anchor, ierr)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        logical,        intent(in) :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-point (or, for anchors, per-anchor) sphere radius

        integer(int32), intent(inout) :: tmp_workspace(n_points)
            !! Reused here as the BFS queue (Step 9).
        real(real64),   intent(inout) :: tmp_val_buf(n_points)
            !! Sort scratch buffer (Step 8, ordering points by intersection degree)
        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer, see f42_utils sort_array
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! sort_array recursion-stack scratch buffers

        logical,        intent(out), allocatable :: point_in_anchor_mask(:,:)
            !! (n_points, n_anchor) membership matrix.
        integer(int32), intent(out), allocatable :: anchor_to_point(:)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(out) :: anchor_count(n_points)
            !! How many anchors cover each point.
        integer(int32), intent(out) :: labels(n_points)
            !! BFS connected-overlap-region label; 0 if not in any intersection.
        integer(int32), intent(out) :: n_anchor
            !! 0 if no anchors were found (caller must skip stitching in that case).
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: i, e
        integer(int32) :: n_edges
        integer(int32), allocatable :: point_to_anchor(:)
        integer(int32), allocatable :: edge_anc1(:), edge_anc2(:), tmp_edge_label_arr(:)
        logical, allocatable :: edge_visited_mask(:), pt_visited_mask(:)
        integer(int32), allocatable :: pt_edge_start(:), pt_edge_list(:)
        integer(int32), allocatable :: edge_pt_start(:), edge_pt_list(:)
        integer(int32), allocatable :: edge_pt_count(:)
        integer(int32) :: n_edges_of_pt(n_points)

        call set_ok(ierr)
        n_anchor = count(is_anchor_mask)
        if (n_anchor == 0) then
            print *, "Error: No anchors found."
            return
        end if

        M_ALLOCATE(anchor_to_point(n_anchor))
        M_ALLOCATE(point_to_anchor(n_points))
        call build_anchor_mapping(n_points, is_anchor_mask, n_anchor, anchor_to_point, point_to_anchor)

        M_ALLOCATE(point_in_anchor_mask(n_points, n_anchor))
        call build_membership_matrix(work_coords, n_points, dim, n_anchor, anchor_to_point, &
                                     sphere_radii, point_in_anchor_mask)

        ! --- STEP 8: Build Edge List + CSR Adjacency Structures ---
        ! An edge = two anchors that share at least one point.
        ! CSR pt->edges: O(degree) lookup instead of O(n_edges) scan in BFS.
        ! CSR edge->pts: O(edge_size) expansion instead of O(n_points) scan.
        call count_anchor_intersection_edges(n_points, n_anchor, point_in_anchor_mask, n_edges)

        M_ALLOCATE(edge_anc1(n_edges))
        M_ALLOCATE(edge_anc2(n_edges))
        M_ALLOCATE(tmp_edge_label_arr(n_edges))
        M_ALLOCATE(edge_visited_mask(n_edges))
        M_ALLOCATE(pt_visited_mask(n_points))
        M_ALLOCATE(edge_pt_count(n_edges))
        tmp_edge_label_arr = 0
        edge_visited_mask   = .false.
        pt_visited_mask     = .false.

        call fill_anchor_intersection_edges(n_points, n_anchor, point_in_anchor_mask, n_edges, &
                                            edge_anc1, edge_anc2)

        call count_intersection_csr_sizes(n_points, n_anchor, n_edges, point_in_anchor_mask, &
                                          edge_anc1, edge_anc2, anchor_count, &
                                          n_edges_of_pt, edge_pt_count)

        ! Phase 8d: Build CSR prefix sums, which fixes the exact sizes of the
        ! two CSR data arrays (pt_edge_list, edge_pt_list)
        M_ALLOCATE(pt_edge_start(n_points + 1))
        pt_edge_start(1) = 1
        do i = 1, n_points
            pt_edge_start(i+1) = pt_edge_start(i) + n_edges_of_pt(i)
        end do
        M_ALLOCATE(pt_edge_list(max(1, pt_edge_start(n_points+1) - 1)))

        M_ALLOCATE(edge_pt_start(n_edges + 1))
        edge_pt_start(1) = 1
        do e = 1, n_edges
            edge_pt_start(e+1) = edge_pt_start(e) + edge_pt_count(e)
        end do
        M_ALLOCATE(edge_pt_list(max(1, edge_pt_start(max(1,n_edges)+1) - 1)))

        deallocate(edge_pt_count)

        call build_intersection_graph(n_points, n_anchor, n_edges, point_in_anchor_mask, &
                                      edge_anc1, edge_anc2, pt_edge_start, edge_pt_start, &
                                      tmp_workspace, tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack, &
                                      pt_visited_mask, edge_visited_mask, tmp_edge_label_arr, &
                                      pt_edge_list, edge_pt_list, labels)

    end subroutine build_intersection_graph_alloc

    !> STEP 8e-9: fill the CSR point<->edge data lists into their
    !| already-allocated (and already prefix-summed) buffers, then
    !| BFS-label the connected overlap regions. No allocation.
    pure subroutine build_intersection_graph(n_points, n_anchor, n_edges, point_in_anchor_mask, &
                                        edge_anc1, edge_anc2, pt_edge_start, edge_pt_start, &
                                        tmp_workspace, tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack, &
                                        pt_visited_mask, edge_visited_mask, tmp_edge_label_arr, &
                                        pt_edge_list, edge_pt_list, labels)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: n_edges
            !! Number of pairwise anchor-intersection edges
        logical,        intent(in) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix
        integer(int32), intent(in) :: edge_anc1(n_edges), edge_anc2(n_edges)
            !! Anchor-pair endpoints of each edge
        integer(int32), intent(in) :: pt_edge_start(n_points + 1)
            !! CSR prefix sums: point -> range of its edge-list entries
        integer(int32), intent(in) :: edge_pt_start(n_edges + 1)
            !! CSR prefix sums: edge -> range of its point-list entries

        integer(int32), intent(inout) :: tmp_workspace(n_points)
            !! Reused here as the BFS queue (Step 9).
        real(real64),   intent(inout) :: tmp_val_buf(n_points)
            !! Sort scratch buffer (ordering points by intersection degree)
        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer, see f42_utils sort_array
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! sort_array recursion-stack scratch buffers
        logical,        intent(inout) :: pt_visited_mask(n_points)
            !! BFS visited-point marker, pre-initialised to .false.
        logical,        intent(inout) :: edge_visited_mask(n_edges)
            !! BFS visited-edge marker, pre-initialised to .false.
        integer(int32), intent(inout) :: tmp_edge_label_arr(n_edges)
            !! BFS region label per edge, pre-initialised to 0

        integer(int32), intent(out) :: pt_edge_list(:)
            !! CSR: edge index for each point-edge incidence
        integer(int32), intent(out) :: edge_pt_list(:)
            !! CSR: point index for each edge-point incidence
        integer(int32), intent(out) :: labels(n_points)
            !! BFS connected-overlap-region label; 0 if not in any intersection.

        ! Local variables
        integer(int32) :: i, j, e, ei, ep, k_idx, tmp_idx
        integer(int32) :: head, tail, curr_idx, current_label
        integer(int32) :: n_edges_of_pt(n_points)
        integer(int32) :: edge_pt_count(n_edges)
        logical        :: found_candidate

        ! Phase 8e: Fill both CSR lists in a single pass
        n_edges_of_pt = 0   ! per-point fill offset
        edge_pt_count = 0   ! per-edge fill offset
        do e = 1, n_edges
            do i = 1, n_points
                if (point_in_anchor_mask(i, edge_anc1(e)) .and. &
                    point_in_anchor_mask(i, edge_anc2(e))) then
                    n_edges_of_pt(i) = n_edges_of_pt(i) + 1
                    pt_edge_list(pt_edge_start(i) + n_edges_of_pt(i) - 1) = e
                    edge_pt_count(e) = edge_pt_count(e) + 1
                    edge_pt_list(edge_pt_start(e) + edge_pt_count(e) - 1) = i
                end if
            end do
        end do
        ! Restore n_edges_of_pt from prefix sums (was used as fill counter)
        do i = 1, n_points
            n_edges_of_pt(i) = pt_edge_start(i+1) - pt_edge_start(i)
        end do

        ! Sort points descending by n_edges_of_pt (most-intersected first)
        do i = 1, n_points ; tmp_perm(i) = i ; end do
        do i = 1, n_points
            tmp_val_buf(i) = real(n_edges_of_pt(i), real64)
        end do
        call sort_array(tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack)
        do i = 1, n_points / 2
            tmp_idx = tmp_perm(i) ; tmp_perm(i) = tmp_perm(n_points-i+1) ; tmp_perm(n_points-i+1) = tmp_idx
        end do

        ! --- STEP 9: BFS on Edge Graph using CSR (Boss Algorithm) ---
        ! Each dequeued point looks up only its own edges: O(degree).
        ! Each visited edge expands only its own points: O(edge_size).
        ! Total BFS cost: O(n_edges + n_points) instead of O(n_points * n_edges).
        labels        = 0
        current_label = 0

        do k_idx = 1, n_points
            i = tmp_perm(k_idx)
            if (n_edges_of_pt(i) == 0) exit   ! sorted desc; rest have 0 edges
            if (pt_visited_mask(i)) cycle

            ! Skip if all edges of this point are already visited (O(degree) check)
            found_candidate = .false.
            do ei = pt_edge_start(i), pt_edge_start(i+1) - 1
                if (.not. edge_visited_mask(pt_edge_list(ei))) then
                    found_candidate = .true.
                    exit
                end if
            end do
            if (.not. found_candidate) cycle

            current_label = current_label + 1
            pt_visited_mask(i) = .true.
            labels(i)     = current_label

            ! BFS queue (tmp_workspace reused as integer queue)
            head = 1 ; tail = 1
            tmp_workspace(1) = i

            do while (head <= tail)
                curr_idx = tmp_workspace(head)
                head = head + 1

                ! Walk only edges of curr_idx via CSR: O(degree) per point
                do ei = pt_edge_start(curr_idx), pt_edge_start(curr_idx+1) - 1
                    e = pt_edge_list(ei)
                    if (edge_visited_mask(e)) cycle
                    edge_visited_mask(e)   = .true.
                    tmp_edge_label_arr(e) = current_label

                    ! Enqueue all unvisited points in this edge via CSR: O(edge_size)
                    do ep = edge_pt_start(e), edge_pt_start(e+1) - 1
                        j = edge_pt_list(ep)
                        if (pt_visited_mask(j)) cycle
                        pt_visited_mask(j)   = .true.
                        labels(j)       = current_label
                        tail            = tail + 1
                        tmp_workspace(tail) = j
                    end do
                end do
            end do
        end do

    end subroutine build_intersection_graph

    !> STEP 10: for every point, blend the tangent-subspace projections of all
    !| anchors covering it, weighted by inverse anchor fit-quality
    !| (1/normal_error). See misc/smoothing_experiments.md section 15 for the
    !| full rationale behind inverse-variance weighting.
    pure subroutine stitch_points(work_coords, n_points, dim, manifold_dim, n_anchor, &
                             point_in_anchor_mask, anchor_to_point, anchor_count, anchor_centers, &
                             tangent_bases, normal_errors, skeleton_coords, &
                             primary_anchor_ids, secondary_anchor_ids)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        logical,        intent(in) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(in) :: anchor_count(n_points)
            !! How many anchors cover each point
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: tangent_bases(dim, manifold_dim, n_points)
            !! Per-anchor tangent basis directions (only meaningful for anchor points)
        real(real64),   intent(in) :: normal_errors(n_points)
            !! Per-anchor fit quality, used as the inverse-variance stitching weight

        real(real64),   intent(out) :: skeleton_coords(dim+1, n_points)
            !! Row 1 = anchor_count; rows 2:dim+1 = stitched position
        integer(int32), intent(out) :: primary_anchor_ids(n_points), secondary_anchor_ids(n_points)
            !! Point index of the highest/second-highest inverse-variance-weighted
            !! anchor covering this point; 0 if not covered by that many anchors

        ! Local variables
        integer(int32) :: i

        ! skeleton_coords(1,   i) = anchor_count   [for R filtering: n_anchors column]
        ! skeleton_coords(2:dim+1, i) = stitched position in R^dim
        ! Each point's stitched position depends only on shared read-only
        ! inputs and writes to its own slice of the output arrays, so the
        ! iterations are independent.
        do concurrent (i = 1:n_points) shared(n_points, dim, manifold_dim, n_anchor, work_coords, &
                                              point_in_anchor_mask, anchor_to_point, anchor_count, &
                                              anchor_centers, tangent_bases, normal_errors, &
                                              skeleton_coords, primary_anchor_ids, secondary_anchor_ids)
            skeleton_coords(1, i) = real(anchor_count(i), real64)

            if (anchor_count(i) >= 2) then
                call stitch_multi_anchor_point(i, work_coords, n_points, dim, manifold_dim, n_anchor, &
                                               point_in_anchor_mask, anchor_to_point, anchor_centers, &
                                               tangent_bases, normal_errors, skeleton_coords(2:dim+1, i), &
                                               primary_anchor_ids(i), secondary_anchor_ids(i))
            else if (anchor_count(i) == 1) then
                call stitch_single_anchor_point(i, work_coords, n_points, dim, manifold_dim, n_anchor, &
                                                point_in_anchor_mask, anchor_to_point, anchor_centers, &
                                                tangent_bases, skeleton_coords(2:dim+1, i), &
                                                primary_anchor_ids(i))
                secondary_anchor_ids(i) = 0
            else
                ! Not covered by any anchor: keep original position
                skeleton_coords(2:dim+1, i) = work_coords(:, i)
                primary_anchor_ids(i) = 0
                secondary_anchor_ids(i) = 0
            end if
        end do

    end subroutine stitch_points

    !> STEP 10, intersection points (anchor_count >= 2): projects the point
    !| onto each covering anchor's tangent plane, reconstructs in R^dim, then
    !| forms an inverse-variance-weighted average -- weights = 1/normal_error,
    !| the standard way to combine several noisy estimates: trust the anchor
    !| whose tangent plane genuinely fits its neighborhood, not just the one
    !| with more nearby points. This collapses points onto a single central
    !| trend while still allowing smooth transitions and preserving
    !| bifurcations (where two anchors genuinely fit about equally well, the
    !| average naturally lands between both branches).
    pure subroutine stitch_multi_anchor_point(i, work_coords, n_points, dim, manifold_dim, n_anchor, &
                                         point_in_anchor_mask, anchor_to_point, anchor_centers, &
                                         tangent_bases, normal_errors, stitched_position, &
                                         primary_anchor_idx, secondary_anchor_idx)

        integer(int32), intent(in) :: i
            !! Index of the point to stitch
        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        logical,        intent(in) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: tangent_bases(dim, manifold_dim, n_points)
            !! Per-anchor tangent basis directions (only meaningful for anchor points)
        real(real64),   intent(in) :: normal_errors(n_points)
            !! Per-anchor fit quality, used as the inverse-variance stitching weight

        real(real64),   intent(out) :: stitched_position(dim)
            !! Inverse-variance-weighted stitched position
        integer(int32), intent(out) :: primary_anchor_idx, secondary_anchor_idx
            !! Point index of the highest/second-highest weighted covering anchor

        integer(int32) :: i_anc, k, base_idx
        real(real64)   :: p_diff(dim), center(dim), proj_val, sigma_i
        real(real64)   :: v_ij(dim), min_dist
        real(real64)   :: primary_weight, secondary_weight

        primary_anchor_idx   = 0
        secondary_anchor_idx = 0
        primary_weight       = -1.0_real64
        secondary_weight     = -1.0_real64
        v_ij     = 0.0_real64
        min_dist = 0.0_real64

        do i_anc = 1, n_anchor
            if (.not. point_in_anchor_mask(i, i_anc)) cycle
            k = anchor_to_point(i_anc)

            p_diff = work_coords(:, i) - anchor_centers(:, k)
            center = anchor_centers(:, k)
            do base_idx = 1, manifold_dim
                proj_val = dot_product(p_diff, tangent_bases(:, base_idx, k))
                center   = center + proj_val * tangent_bases(:, base_idx, k)
            end do

            ! Inverse-variance weighting: 1/normal_error, so a well-fit
            ! anchor dominates over a poorly-fit one regardless of density
            sigma_i  = 1.0_real64 / max(normal_errors(k), 1.0e-12_real64)
            if (sigma_i > primary_weight) then
                secondary_weight = primary_weight
                secondary_anchor_idx = primary_anchor_idx
                primary_weight = sigma_i
                primary_anchor_idx = k
            else if (sigma_i > secondary_weight) then
                secondary_weight = sigma_i
                secondary_anchor_idx = k
            end if
            v_ij     = v_ij     + sigma_i * center
            min_dist = min_dist + sigma_i
        end do

        if (min_dist > 0.0_real64) then
            stitched_position = v_ij / min_dist
        else
            stitched_position = work_coords(:, i)
        end if

    end subroutine stitch_multi_anchor_point

    !> STEP 10, single-anchor points (anchor_count == 1): projects the point
    !| onto the tangent subspace of the one covering anchor, so all points in
    !| the same patch collapse onto the manifold_dim-dimensional central trend.
    pure subroutine stitch_single_anchor_point(i, work_coords, n_points, dim, manifold_dim, n_anchor, &
                                          point_in_anchor_mask, anchor_to_point, anchor_centers, &
                                          tangent_bases, stitched_position, primary_anchor_idx)

        integer(int32), intent(in) :: i
            !! Index of the point to stitch
        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of work_coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        real(real64),   intent(in) :: work_coords(dim, n_points)
            !! Current point coordinates
        logical,        intent(in) :: point_in_anchor_mask(n_points, n_anchor)
            !! (n_points, n_anchor) membership matrix
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: tangent_bases(dim, manifold_dim, n_points)
            !! Per-anchor tangent basis directions (only meaningful for anchor points)

        real(real64),   intent(out) :: stitched_position(dim)
            !! Position projected onto the covering anchor's tangent subspace
        integer(int32), intent(out) :: primary_anchor_idx
            !! Point index of the covering anchor

        integer(int32) :: i_anc, k, base_idx
        real(real64)   :: p_diff(dim), proj_val

        do i_anc = 1, n_anchor
            if (.not. point_in_anchor_mask(i, i_anc)) cycle
            k = anchor_to_point(i_anc)
            primary_anchor_idx = k
            p_diff = work_coords(:, i) - anchor_centers(:, k)
            stitched_position = anchor_centers(:, k)
            do base_idx = 1, manifold_dim
                proj_val = dot_product(p_diff, tangent_bases(:, base_idx, k))
                stitched_position = stitched_position + proj_val * tangent_bases(:, base_idx, k)
            end do
            exit
        end do

    end subroutine stitch_single_anchor_point

    !> Computes the outer convergence-loop displacement tolerance as a small
    !| fraction (relative_conv_tol) of the dataset's own resolution -- the
    !| median nearest-neighbor distance in the ORIGINAL coordinates -- rather
    !| than an absolute number tied to whatever units this particular dataset
    !| happens to use, so it self-calibrates instead of needing to be
    !| re-tuned per dataset.
    pure subroutine compute_relative_conv_tol(coords, n_points, dim, tmp_kd_indices, tmp_dim_order, &
                                         relative_conv_tol, tmp_perm, tmp_l_stack, tmp_r_stack, conv_tol)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of coords
        real(real64),   intent(in) :: coords(dim, n_points)
            !! Original point coordinates
        integer(int32), intent(in) :: tmp_kd_indices(n_points)
            !! KD-tree index buffer, see kd_tree module
        integer(int32), intent(in) :: tmp_dim_order(dim)
            !! KD-tree dimension-splitting order
        real(real64),   intent(in) :: relative_conv_tol
            !! Convergence tolerance as a fraction of the dataset's own resolution

        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer, see f42_utils sort_array
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! sort_array recursion-stack scratch buffers

        real(real64),   intent(out) :: conv_tol
            !! Displacement tolerance for the outer convergence loop

        integer(int32) :: qi, qj, self_pos2, k_query2, ierr
        integer(int32) :: nn_idx(2)
        real(real64)   :: nn_d(2)
        real(real64)   :: nn_dist(n_points)

        do qi = 1, n_points
            k_query2 = min(2, n_points)
            call kd_knn_query(coords, tmp_kd_indices, dim, n_points, tmp_dim_order, &
                              coords(:, qi), k_query2, nn_idx(1:k_query2), nn_d(1:k_query2), ierr)
            if (.not. is_ok(ierr)) then
                call set_ok(ierr)
                nn_dist(qi) = 0.0_real64
                cycle
            end if
            self_pos2 = 0
            do qj = 1, k_query2
                if (nn_idx(qj) == qi) then
                    self_pos2 = qj
                    exit
                end if
            end do
            if (self_pos2 > 0 .and. self_pos2 < k_query2) then
                nn_d(self_pos2) = nn_d(k_query2)
            end if
            if (k_query2 >= 2) then
                nn_dist(qi) = nn_d(1)
            else
                nn_dist(qi) = 0.0_real64
            end if
        end do

        do qi = 1, n_points
            tmp_perm(qi) = qi
        end do
        call sort_array(nn_dist, tmp_perm, tmp_l_stack, tmp_r_stack)
        conv_tol = relative_conv_tol * max(nn_dist(tmp_perm((n_points + 1) / 2)), 1.0e-12_real64)

    end subroutine compute_relative_conv_tol

    !> Public entry point: iterates lomanle_pass until convergence.
    !| Every iteration recomputes KD-tree, adaptive radii, densities, atlas and
    !| stitching on the updated coordinates, so all quantities converge together.
    !| Callers that don't want to manage the ~15 work buffers themselves should
    !| use lomanle_compute_alloc instead.
    subroutine lomanle_compute(coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                               o_max, o_min, stability_threshold, scale_factor, &
                               max_iterations, relative_conv_tol, &
                               tmp_kd_indices, tmp_workspace, tmp_val_buf, tmp_perm, &
                               tmp_l_stack, tmp_r_stack, tmp_rec_stack, &
                               tmp_dim_order, lwork, &
                               sphere_radii, densities, gap_values, &
                               normal_errors, stability_values, k_selected, growth_stopped_complex, &
                               is_anchor_mask, tangent_bases, tangent_scales, labels, &
                               skeleton_coords, skeleton_iter1, &
                               radii_1, densities_1, gap_1, is_anchor_mask_1, labels_1, &
                               tangent_bases_1, tangent_scales_1, &
                               primary_anchor_final, secondary_anchor_final, &
                               primary_anchor_1, secondary_anchor_1, &
                               anchor_centers_1, anchor_centers_final, &
                               max_edges, edge_from_1, edge_to_1, n_edges_1, anchor_role_1, &
                               edge_from_final, edge_to_final, n_edges_final, anchor_role_final, &
                               ierr)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: k_min
            !! Minimum neighborhood size to start adaptive growth from
        integer(int32), intent(in) :: lwork
            !! Size of the LAPACK dsyev scratch workspace
        real(real64),   intent(in) :: coords(dim, n_points)
            !! Original (unstitched) point coordinates
        real(real64),   intent(in) :: g_threshold
            !! Spectral-gap threshold that stops neighborhood growth
        real(real64),   intent(in) :: o_max
            !! Maximum allowed anchor-sphere overlap ratio
        real(real64),   intent(in) :: o_min
            !! Minimum required anchor-sphere overlap ratio
        real(real64),   intent(in) :: stability_threshold
            !! Minimum tangent-basis stability to keep growing a neighborhood
        real(real64),   intent(in) :: scale_factor
            !! Cap on sphere radius as a multiple of the point's local scale
        integer(int32), intent(in) :: max_iterations
            !! Outer convergence loop iteration cap
        real(real64),   intent(in) :: relative_conv_tol
            !! Convergence tolerance as a fraction of the dataset's own resolution
            !! (median nearest-neighbor distance in the original coordinates),
            !! not an absolute number tied to any one dataset's coordinate units.

        real(real64),   intent(out) :: skeleton_coords(dim+1, n_points)
            !! Final converged stitched positions (row 1 = anchor_count)
        real(real64),   intent(out) :: skeleton_iter1(dim+1, n_points)
            !! Stitched positions after iteration 1
        real(real64),   intent(out) :: radii_1(n_points), densities_1(n_points), gap_1(n_points)
            !! Diagnostic quantities from iteration 1
        logical,        intent(out) :: is_anchor_mask_1(n_points)
            !! .true. for points selected as atlas anchors in iteration 1
        integer(int32), intent(out) :: labels_1(n_points)
            !! BFS intersection-region label from iteration 1
        real(real64),   intent(out) :: tangent_bases_1(dim, manifold_dim, n_points)
            !! Anchor tangent bases from iteration 1
        real(real64),   intent(out) :: tangent_scales_1(manifold_dim, n_points)
            !! Anchor tangent extents from iteration 1
        integer(int32), intent(out) :: primary_anchor_final(n_points), secondary_anchor_final(n_points)
            !! Primary/secondary stitching anchor ids from the converged (final) pass
        integer(int32), intent(out) :: primary_anchor_1(n_points), secondary_anchor_1(n_points)
            !! Primary/secondary stitching anchor ids from iteration 1
        real(real64),   intent(out) :: anchor_centers_1(dim, n_points), anchor_centers_final(dim, n_points)
            !! Anchor centroids from iteration 1 / from the converged (final) pass

        integer(int32), intent(in)  :: max_edges
            !! Size of the caller-provided edge_from*/edge_to* buffers. A generous
            !! bound is 3 * n_points: internal per-anchor chains contribute at most
            !! ~2 * n_points edges, MST bridges at most n_anchor - 1 more.
        integer(int32), intent(out) :: edge_from_1(max_edges), edge_to_1(max_edges)
            !! Backbone edge list (point-index pairs) from iteration 1
        integer(int32), intent(out) :: n_edges_1
            !! Number of meaningful entries in edge_from_1/edge_to_1
        integer(int32), intent(out) :: edge_from_final(max_edges), edge_to_final(max_edges)
            !! Backbone edge list (point-index pairs) from the converged (final) pass
        integer(int32), intent(out) :: n_edges_final
            !! Number of meaningful entries in edge_from_final/edge_to_final
        integer(int32), intent(out) :: anchor_role_1(n_points), anchor_role_final(n_points)
            !! 0 = not an anchor; 1 = endpoint; 2 = pass-through; 3 = branch/junction

        integer(int32), intent(inout) :: tmp_kd_indices(n_points)
            !! KD-tree scratch buffer, see kd_tree module
        integer(int32), intent(inout) :: tmp_workspace(n_points)
            !! KD-tree/sort scratch buffer, reused across steps
        integer(int32), intent(inout) :: tmp_perm(n_points)
            !! Sort-permutation scratch buffer, see f42_utils sort_array
        integer(int32), intent(inout) :: tmp_l_stack(n_points), tmp_r_stack(n_points)
            !! KD-tree/sort_array recursion-stack scratch buffers
        integer(int32), intent(inout) :: tmp_rec_stack(3, n_points)
            !! KD-tree build recursion-stack scratch buffer
        integer(int32), intent(inout) :: tmp_dim_order(dim)
            !! KD-tree dimension-splitting order scratch buffer
        real(real64),   intent(inout) :: tmp_val_buf(n_points)
            !! KD-tree/sort scratch buffer, reused across steps

        real(real64),   intent(out) :: sphere_radii(n_points)
            !! Per-point adaptive neighborhood radius, from the final pass
        real(real64),   intent(out) :: densities(n_points)
            !! Per-point local density estimate, from the final pass
        real(real64),   intent(out) :: gap_values(n_points)
            !! Per-point spectral gap, from the final pass
        real(real64),   intent(out) :: normal_errors(n_points)
            !! Per-point/anchor fit quality, from the final pass
        real(real64),   intent(out) :: stability_values(n_points)
            !! Per-point tangent-basis stability, from the final pass
        integer(int32), intent(out) :: k_selected(n_points)
            !! Per-point kept neighborhood size, from the final pass
        logical,        intent(out) :: growth_stopped_complex(n_points)
            !! Per-point instability-triggered early stop flag, from the final pass
        logical,        intent(out) :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors, from the final pass
        real(real64),   intent(out) :: tangent_bases(dim, manifold_dim, n_points)
            !! Anchor tangent bases from the final pass
        real(real64),   intent(out) :: tangent_scales(manifold_dim, n_points)
            !! Anchor tangent extents from the final pass
        integer(int32), intent(out) :: labels(n_points)
            !! BFS intersection-region label, from the final pass
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: iter, i
        logical        :: converged
        real(real64)   :: max_disp, disp, conv_tol
        real(real64)   :: work_coords(dim, n_points)
        real(real64)   :: anchor_centers(dim, n_points)

        call set_ok(ierr)
        work_coords = coords
        converged   = .false.

        ! --- Dataset-relative convergence tolerance ---
        ! conv_tol is a small fraction (relative_conv_tol) of the data's own
        ! resolution -- the median nearest-neighbor distance in the ORIGINAL
        ! coordinates -- rather than an absolute number tied to whatever units
        ! this particular dataset happens to use, so it self-calibrates
        ! instead of needing to be re-tuned per dataset.
        do i = 1, dim
            tmp_dim_order(i) = i
        end do
        call build_kd_index(coords, dim, n_points, tmp_kd_indices, tmp_dim_order, &
                            tmp_workspace, tmp_val_buf, tmp_perm, tmp_l_stack, tmp_r_stack, tmp_rec_stack, ierr)
        if (.not. is_ok(ierr)) return

        call compute_relative_conv_tol(coords, n_points, dim, tmp_kd_indices, tmp_dim_order, &
                                       relative_conv_tol, tmp_perm, tmp_l_stack, tmp_r_stack, conv_tol)
        skeleton_iter1 = 0.0_real64
        radii_1 = 0.0_real64 ; densities_1 = 0.0_real64 ; gap_1 = 0.0_real64
        is_anchor_mask_1 = .false. ; labels_1 = 0
        tangent_bases_1 = 0.0_real64 ; tangent_scales_1 = 0.0_real64
        primary_anchor_final = 0 ; secondary_anchor_final = 0
        primary_anchor_1 = 0 ; secondary_anchor_1 = 0
        ! lomanle_pass reads tangent_bases on entry to warm-start the tangent-stability
        ! check with the previous outer iteration's basis; zero it here so the very
        ! first call sees "no previous tangent".
        tangent_bases = 0.0_real64

        do iter = 1, max_iterations
            call lomanle_pass(work_coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                              o_max, o_min, stability_threshold, scale_factor, &
                              tmp_kd_indices, tmp_workspace, tmp_val_buf, tmp_perm, &
                              tmp_l_stack, tmp_r_stack, tmp_rec_stack, &
                              tmp_dim_order, lwork, &
                              sphere_radii, densities, gap_values, &
                              normal_errors, stability_values, k_selected, growth_stopped_complex, &
                              is_anchor_mask, tangent_bases, tangent_scales, labels, skeleton_coords, &
                              primary_anchor_final, secondary_anchor_final, anchor_centers, ierr)
            if (.not. is_ok(ierr)) return

            ! Save all iteration-1 diagnostics for the report
            if (iter == 1) then
                skeleton_iter1   = skeleton_coords
                radii_1          = sphere_radii
                densities_1      = densities
                gap_1            = gap_values
                is_anchor_mask_1      = is_anchor_mask
                labels_1         = labels
                tangent_bases_1  = tangent_bases
                tangent_scales_1 = tangent_scales
                primary_anchor_1 = primary_anchor_final
                secondary_anchor_1 = secondary_anchor_final
                anchor_centers_1 = anchor_centers
            end if

            ! Convergence: max displacement of any point this iteration
            max_disp = 0.0_real64
            do i = 1, n_points
                disp = sqrt(sum((skeleton_coords(2:dim+1, i) - work_coords(:, i))**2))
                if (disp > max_disp) max_disp = disp
            end do

            ! Update positions for next iteration
            do i = 1, n_points
                work_coords(:, i) = skeleton_coords(2:dim+1, i)
            end do

            if (max_disp < conv_tol) then
                converged = .true.
                exit
            end if
        end do

        if (.not. converged) then
            print '(A,I3,A,ES10.3,A,ES10.3,A)', "  WARNING: did not converge in ", max_iterations, &
                " iterations. Last max_disp=", max_disp, " (conv_tol=", conv_tol, &
                "; consider raising max_iterations or relative_conv_tol)"
        end if

        anchor_centers_final = anchor_centers

        ! --- Backbone graph: built here (not by the caller) from the anchor
        ! adjacency graph's MST, for both the iteration-1 snapshot and the
        ! converged result -- see build_skeleton_edges_alloc for the full rationale.
        call build_skeleton_edges_alloc(n_points, dim, manifold_dim, is_anchor_mask_1, anchor_centers_1, &
                                  tangent_bases_1, primary_anchor_1, secondary_anchor_1, &
                                  radii_1, skeleton_iter1, max_edges, edge_from_1, edge_to_1, &
                                  n_edges_1, anchor_role_1, ierr)
        if (.not. is_ok(ierr)) return

        call build_skeleton_edges_alloc(n_points, dim, manifold_dim, is_anchor_mask, anchor_centers_final, &
                                  tangent_bases, primary_anchor_final, secondary_anchor_final, &
                                  sphere_radii, skeleton_coords, max_edges, edge_from_final, edge_to_final, &
                                  n_edges_final, anchor_role_final, ierr)

    end subroutine lomanle_compute

    !> Convenience wrapper around lomanle_compute: allocates every work
    !| buffer (KD-tree/sort scratch, backbone-edge buffers, and the
    !| per-iteration diagnostics the caller never actually reads)
    !| internally, so callers only need to supply the problem itself and
    !| receive the outputs they actually use. See misc/Fortran_Coding_Guides.tex
    !| section "Routine Organisation: Helper / Main / Alloc Pattern".
    subroutine lomanle_compute_alloc(coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                                     o_max, o_min, stability_threshold, scale_factor, &
                                     max_iterations, relative_conv_tol, &
                                     sphere_radii, densities, gap_values, &
                                     normal_errors, stability_values, k_selected, growth_stopped_complex, &
                                     skeleton_coords, skeleton_iter1, &
                                     radii_1, densities_1, gap_1, is_anchor_mask_1, labels_1, &
                                     tangent_bases_1, tangent_scales_1, &
                                     primary_anchor_final, secondary_anchor_final, &
                                     primary_anchor_1, secondary_anchor_1, &
                                     anchor_centers_1, anchor_centers_final, &
                                     edge_from_1, edge_to_1, n_edges_1, anchor_role_1, &
                                     edge_from_final, edge_to_final, n_edges_final, anchor_role_final, &
                                     ierr)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of coords
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: k_min
            !! Minimum neighborhood size to start adaptive growth from
        real(real64),   intent(in) :: coords(dim, n_points)
            !! Original (unstitched) point coordinates
        real(real64),   intent(in) :: g_threshold
            !! Spectral-gap threshold that stops neighborhood growth
        real(real64),   intent(in) :: o_max
            !! Maximum allowed anchor-sphere overlap ratio
        real(real64),   intent(in) :: o_min
            !! Minimum required anchor-sphere overlap ratio
        real(real64),   intent(in) :: stability_threshold
            !! Minimum tangent-basis stability to keep growing a neighborhood
        real(real64),   intent(in) :: scale_factor
            !! Cap on sphere radius as a multiple of the point's local scale
        integer(int32), intent(in) :: max_iterations
            !! Outer convergence loop iteration cap
        real(real64),   intent(in) :: relative_conv_tol
            !! Convergence tolerance as a fraction of the dataset's own resolution

        real(real64),   intent(out) :: sphere_radii(n_points)
            !! Per-point adaptive neighborhood radius, from the final pass
        real(real64),   intent(out) :: densities(n_points)
            !! Per-point local density estimate, from the final pass
        real(real64),   intent(out) :: gap_values(n_points)
            !! Per-point spectral gap, from the final pass
        real(real64),   intent(out) :: normal_errors(n_points)
            !! Per-point/anchor fit quality, from the final pass
        real(real64),   intent(out) :: stability_values(n_points)
            !! Per-point tangent-basis stability, from the final pass
        integer(int32), intent(out) :: k_selected(n_points)
            !! Per-point kept neighborhood size, from the final pass
        logical,        intent(out) :: growth_stopped_complex(n_points)
            !! Per-point instability-triggered early stop flag, from the final pass

        real(real64),   intent(out) :: skeleton_coords(dim+1, n_points)
            !! Final converged stitched positions (row 1 = anchor_count)
        real(real64),   intent(out) :: skeleton_iter1(dim+1, n_points)
            !! Stitched positions after iteration 1
        real(real64),   intent(out) :: radii_1(n_points), densities_1(n_points), gap_1(n_points)
            !! Diagnostic quantities from iteration 1
        logical,        intent(out) :: is_anchor_mask_1(n_points)
            !! .true. for points selected as atlas anchors in iteration 1
        integer(int32), intent(out) :: labels_1(n_points)
            !! BFS intersection-region label from iteration 1
        real(real64),   intent(out) :: tangent_bases_1(dim, manifold_dim, n_points)
            !! Anchor tangent bases from iteration 1
        real(real64),   intent(out) :: tangent_scales_1(manifold_dim, n_points)
            !! Anchor tangent extents from iteration 1
        integer(int32), intent(out) :: primary_anchor_final(n_points), secondary_anchor_final(n_points)
            !! Primary/secondary stitching anchor ids from the converged (final) pass
        integer(int32), intent(out) :: primary_anchor_1(n_points), secondary_anchor_1(n_points)
            !! Primary/secondary stitching anchor ids from iteration 1
        real(real64),   intent(out) :: anchor_centers_1(dim, n_points), anchor_centers_final(dim, n_points)
            !! Anchor centroids from iteration 1 / from the converged (final) pass

        integer(int32), intent(out), allocatable :: edge_from_1(:), edge_to_1(:)
            !! Backbone edge list (point-index pairs) from iteration 1
        integer(int32), intent(out) :: n_edges_1
            !! Number of meaningful entries in edge_from_1/edge_to_1
        integer(int32), intent(out) :: anchor_role_1(n_points)
            !! 0 = not an anchor; 1 = endpoint; 2 = pass-through; 3 = branch/junction (iteration 1)
        integer(int32), intent(out), allocatable :: edge_from_final(:), edge_to_final(:)
            !! Backbone edge list (point-index pairs) from the converged (final) pass
        integer(int32), intent(out) :: n_edges_final
            !! Number of meaningful entries in edge_from_final/edge_to_final
        integer(int32), intent(out) :: anchor_role_final(n_points)
            !! 0 = not an anchor; 1 = endpoint; 2 = pass-through; 3 = branch/junction (final)
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        ! Work buffers: never seen by the caller.
        integer(int32) :: lwork, max_edges
        integer(int32), allocatable :: tmp_kd_indices(:), tmp_workspace(:), tmp_perm(:)
        integer(int32), allocatable :: tmp_l_stack(:), tmp_r_stack(:), tmp_rec_stack(:,:), tmp_dim_order(:)
        real(real64),   allocatable :: tmp_val_buf(:)

        ! Per-iteration outputs of lomanle_compute that the caller never reads
        ! (only the "_1"/"_final"-suffixed snapshots and the final-iteration
        ! growth diagnostics above are meaningful downstream).
        logical,        allocatable :: is_anchor_mask(:)
        real(real64),   allocatable :: tangent_bases(:,:,:), tangent_scales(:,:)
        integer(int32), allocatable :: labels(:)

        call set_ok(ierr)
        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
        call validate_dimension_size(dim, ierr, arg_pos=3_int32)
        call validate_dimension_size(manifold_dim, ierr, arg_pos=4_int32)
        call validate_dimension_size(k_min, ierr, arg_pos=5_int32)
        call validate_dimension_size(max_iterations, ierr, arg_pos=11_int32)
        if (is_err(ierr)) return

        call validate_in_range_int(manifold_dim, ierr, arg_pos=4_int32, max=dim)
        call validate_all_in_range_real(coords, dim * n_points, ierr, arg_pos=1_int32)
        call validate_in_range_real(g_threshold, ierr, arg_pos=6_int32, min=above(0.0_real64))
        call validate_in_range_real(o_max, ierr, arg_pos=7_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(o_min, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(stability_threshold, ierr, arg_pos=9_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(scale_factor, ierr, arg_pos=10_int32, min=above(0.0_real64))
        call validate_in_range_real(relative_conv_tol, ierr, arg_pos=12_int32, min=above(0.0_real64))
        if (is_err(ierr)) return

        if (o_min > o_max) then
            call set_err(ierr, ERR_INVALID_INPUT, arg_pos=8_int32)
            return
        end if

        lwork     = max(1, 3 * dim - 1)
        max_edges = max(1, 3 * n_points)

        M_ALLOCATE(tmp_kd_indices(n_points))
        M_ALLOCATE(tmp_workspace(n_points))
        M_ALLOCATE(tmp_perm(n_points))
        M_ALLOCATE(tmp_l_stack(n_points))
        M_ALLOCATE(tmp_r_stack(n_points))
        M_ALLOCATE(tmp_rec_stack(3, n_points))
        M_ALLOCATE(tmp_dim_order(dim))
        M_ALLOCATE(tmp_val_buf(n_points))
        M_ALLOCATE(is_anchor_mask(n_points))
        M_ALLOCATE(labels(n_points))
        M_ALLOCATE(tangent_bases(dim, manifold_dim, n_points))
        M_ALLOCATE(tangent_scales(manifold_dim, n_points))
        M_ALLOCATE(edge_from_1(max_edges))
        M_ALLOCATE(edge_to_1(max_edges))
        M_ALLOCATE(edge_from_final(max_edges))
        M_ALLOCATE(edge_to_final(max_edges))

        call lomanle_compute(coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                             o_max, o_min, stability_threshold, scale_factor, &
                             max_iterations, relative_conv_tol, &
                             tmp_kd_indices, tmp_workspace, tmp_val_buf, tmp_perm, &
                             tmp_l_stack, tmp_r_stack, tmp_rec_stack, &
                             tmp_dim_order, lwork, &
                             sphere_radii, densities, gap_values, &
                             normal_errors, stability_values, k_selected, growth_stopped_complex, &
                             is_anchor_mask, tangent_bases, tangent_scales, labels, &
                             skeleton_coords, skeleton_iter1, &
                             radii_1, densities_1, gap_1, is_anchor_mask_1, labels_1, &
                             tangent_bases_1, tangent_scales_1, &
                             primary_anchor_final, secondary_anchor_final, &
                             primary_anchor_1, secondary_anchor_1, &
                             anchor_centers_1, anchor_centers_final, &
                             max_edges, edge_from_1, edge_to_1, n_edges_1, anchor_role_1, &
                             edge_from_final, edge_to_final, n_edges_final, anchor_role_final, &
                             ierr)

        deallocate(tmp_kd_indices, tmp_workspace, tmp_perm, tmp_l_stack, tmp_r_stack, tmp_rec_stack, &
                   tmp_dim_order, tmp_val_buf, &
                   is_anchor_mask, labels, tangent_bases, tangent_scales)

    end subroutine lomanle_compute_alloc

    !> Builds a topological skeleton (backbone) graph over the stitched points.
    !|
    !| Rationale: connecting individual points into a line by having each point
    !| greedily pick its own "successor" (nearest compatible neighbor, direction
    !| cones, loop-avoidance walks, etc.) fights the shape of the actual problem.
    !| The topology of the manifold -- how many branches, where they meet, where
    !| they end -- is a property of the ANCHORS and how their charts overlap, not
    !| of any individual point. That anchor-level graph is exactly what Step 10's
    !| primary/secondary anchor assignment already encodes: whenever a point blends
    !| between two anchors, those two anchors are, by construction, part of the
    !| same connected piece of the manifold.
    !|
    !| So instead of a per-point next-pointer chain, this builds:
    !|   1. A candidate graph over anchors (nodes), with an edge between two
    !|      anchors whenever some point uses them as its primary+secondary pair
    !|      (build_anchor_mapping, build_anchor_mst).
    !|   2. A minimum spanning tree of that graph (Kruskal + union-find), which
    !|      prunes redundant/spurious overlaps and guarantees an acyclic backbone
    !|      with no ad-hoc loop-detection needed downstream (build_anchor_mst).
    !|   3. Anchor degree in the MST classifies each anchor structurally:
    !|      endpoint (<=1), pass-through (2), or branch/junction (>=3) -- these
    !|      facts fall out of the graph instead of being inferred heuristically
    !|      (build_anchor_mst).
    !|   4. Every point that blends into a given anchor (as primary or secondary)
    !|      is threaded into a chain ordered by its position along that anchor's
    !|      own tangent axis, giving a fine polyline within each chart
    !|      (build_member_chains).
    !|   5. Consecutive anchors in the MST are bridged through their mutually
    !|      closest member points, so the tree edge threads through real stitched
    !|      points instead of jumping straight from one anchor center to the other
    !|      (build_branch_adjacency, then the branch walk in build_skeleton_edges,
    !|      via emit_branch).
    !|
    !| The returned edge list is over POINT indices (not next_point per point),
    !| since a branch anchor genuinely needs more than one outgoing edge.
    subroutine build_skeleton_edges_alloc(n_points, dim, manifold_dim, is_anchor_mask, anchor_centers, &
                                    tangent_bases, primary_anchor_ids, secondary_anchor_ids, &
                                    sphere_radii, skeleton_coords, max_edges, edge_from, edge_to, &
                                    n_edges, anchor_role, ierr)

        integer(int32), intent(in)  :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in)  :: dim
            !! Ambient dimension of the coordinates
        integer(int32), intent(in)  :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in)  :: max_edges
            !! Size of the caller-provided edge_from/edge_to buffers
        logical,        intent(in)  :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors
        real(real64),   intent(in)  :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in)  :: tangent_bases(dim, manifold_dim, n_points)
            !! Per-anchor tangent basis directions (only meaningful for anchor points)
        integer(int32), intent(in)  :: primary_anchor_ids(n_points), secondary_anchor_ids(n_points)
            !! Point index of the highest/second-highest weighted stitching anchor
            !! (Step 10); 0 if not covered by that many anchors
        real(real64),   intent(in)  :: sphere_radii(n_points)
            !! Per-point (or, for anchors, per-anchor) sphere radius
        real(real64),   intent(in)  :: skeleton_coords(dim+1, n_points)
            !! Stitched point positions (row 1 = anchor_count, rows 2:dim+1 = position)
        integer(int32), intent(out) :: edge_from(max_edges), edge_to(max_edges)
            !! Backbone edge list, as point-index pairs
        integer(int32), intent(out) :: n_edges
            !! Number of meaningful entries in edge_from/edge_to
        integer(int32), intent(out) :: anchor_role(n_points)
            !! 0 = not an anchor; 1 = endpoint (backbone degree <= 1);
            !! 2 = pass-through (degree 2); 3 = branch/junction (degree >= 3)
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: n_anchor
        integer(int32), allocatable :: anchor_to_point(:), point_to_anchor(:)
        integer(int32), allocatable :: mst_a1(:), mst_a2(:), anchor_degree(:)
        integer(int32) :: n_mst
        integer(int32), allocatable :: member_start(:), member_list(:)
        integer(int32), allocatable :: adj_start(:), adj_list(:), adj_edge_id(:)
        integer(int32), allocatable :: tmp_branch_pts(:), tmp_branch_perm(:), tmp_branch_lstack(:), tmp_branch_rstack(:)
        real(real64),   allocatable :: tmp_branch_g(:)
        logical,        allocatable :: edge_used_mask(:)
        integer(int32), allocatable :: anchor_branch_id(:)
        real(real64),   allocatable :: branch_axis(:,:), branch_offset(:)
        integer(int32), allocatable :: tmp_path_buf(:)

        call set_ok(ierr)
        anchor_role = 0
        n_edges = 0
        n_anchor = count(is_anchor_mask)
        if (n_anchor == 0) return

        ! 1. Anchor index mapping (compact <-> original point indices)
        M_ALLOCATE(anchor_to_point(n_anchor))
        M_ALLOCATE(point_to_anchor(n_points))
        call build_anchor_mapping(n_points, is_anchor_mask, n_anchor, anchor_to_point, point_to_anchor)

        ! 2. Anchor MST (tier-1 topological + tier-2 geometric candidates, Kruskal)
        M_ALLOCATE(mst_a1(n_anchor))
        M_ALLOCATE(mst_a2(n_anchor))
        M_ALLOCATE(anchor_degree(n_anchor))
        call build_anchor_mst_alloc(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                    sphere_radii, point_to_anchor, primary_anchor_ids, &
                                    secondary_anchor_ids, mst_a1, mst_a2, n_mst, anchor_degree, &
                                    anchor_role, ierr)
        if (.not. is_ok(ierr)) return

        ! Scratch buffers shared by the special-anchor chains (build_member_chains)
        ! and by emit_branch (branch interiors): sized to n_points, the safe upper
        ! bound for either a single hub's cluster or one branch's pooled points.
        M_ALLOCATE(tmp_branch_pts(n_points))
        M_ALLOCATE(tmp_branch_g(n_points))
        M_ALLOCATE(tmp_branch_perm(n_points))
        M_ALLOCATE(tmp_branch_lstack(n_points))
        M_ALLOCATE(tmp_branch_rstack(n_points))

        ! 4. Thread every point into its primary anchor's member chain; also
        ! locally orders and emits edges for each "special" anchor's own cluster.
        ! member_list is sized to n_points: a safe, tight upper bound, since
        ! each point is threaded into at most one anchor's member chain.
        M_ALLOCATE(member_start(n_anchor + 1))
        M_ALLOCATE(member_list(n_points))
        call build_member_chains(n_points, dim, manifold_dim, n_anchor, anchor_to_point, &
                                 point_to_anchor, primary_anchor_ids, anchor_degree, &
                                 skeleton_coords, anchor_centers, tangent_bases, max_edges, &
                                 tmp_branch_pts, tmp_branch_g, tmp_branch_perm, tmp_branch_lstack, tmp_branch_rstack, &
                                 member_start, member_list, edge_from, edge_to, n_edges)

        ! CSR adjacency of the anchor MST, so build_skeleton_edges below can walk
        ! each branch in O(branch length) instead of re-scanning all MST edges.
        ! adj_list/adj_edge_id are sized to exactly 2*n_mst: each MST edge
        ! contributes exactly one adjacency-list entry per endpoint.
        M_ALLOCATE(adj_start(n_anchor + 1))
        M_ALLOCATE(adj_list(2 * max(1, n_mst)))
        M_ALLOCATE(adj_edge_id(2 * max(1, n_mst)))
        call build_branch_adjacency(n_anchor, mst_a1, mst_a2, n_mst, anchor_degree, &
                                    adj_start, adj_list, adj_edge_id)

        ! 5. Walk every branch and emit its points as ONE continuous chain --
        ! see build_skeleton_edges for the full rationale.
        M_ALLOCATE(edge_used_mask(max(1, n_mst)))
        M_ALLOCATE(anchor_branch_id(n_anchor))
        M_ALLOCATE(branch_axis(dim, n_anchor))
        M_ALLOCATE(branch_offset(n_anchor))
        M_ALLOCATE(tmp_path_buf(max(1, n_anchor)))

        call build_skeleton_edges(n_points, dim, manifold_dim, n_anchor, &
                                  anchor_to_point, point_to_anchor, tangent_bases, anchor_centers, &
                                  primary_anchor_ids, secondary_anchor_ids, skeleton_coords, max_edges, &
                                  n_mst, anchor_degree, adj_start, adj_list, adj_edge_id, &
                                  member_start, member_list, &
                                  edge_used_mask, anchor_branch_id, branch_axis, branch_offset, tmp_path_buf, &
                                  tmp_branch_pts, tmp_branch_g, tmp_branch_perm, tmp_branch_lstack, tmp_branch_rstack, &
                                  edge_from, edge_to, n_edges)

    end subroutine build_skeleton_edges_alloc

    !> Walks the anchor MST from every endpoint/junction anchor (degree /= 2),
    !| following each chain of pass-through anchors (degree == 2) until the
    !| next hub is reached, and emits every such branch as one continuous
    !| ordered chain of points via emit_branch.
    !|
    !| Rather than each anchor ordering (and drawing) only its own members
    !| and then jumping to the next anchor's own ordering via a single
    !| bridge point, every point along a whole branch gets ONE shared,
    !| monotonically increasing coordinate: each anchor's local tangential
    !| offset is added on top of a running "distance travelled so far"
    !| accumulated anchor-to-anchor along the branch (see emit_branch),
    !| instead of resetting to zero at every chart boundary. Sorting the
    !| WHOLE branch by that single coordinate -- instead of "all of anchor
    !| A, then all of anchor B" -- lets points near a chart boundary
    !| interleave in their true order even when two neighboring anchors'
    !| independently-fit local tangent lines don't perfectly agree.
    pure subroutine build_skeleton_edges(n_points, dim, manifold_dim, n_anchor, &
                                    anchor_to_point, point_to_anchor, tangent_bases, anchor_centers, &
                                    primary_anchor_ids, secondary_anchor_ids, skeleton_coords, max_edges, &
                                    n_mst, anchor_degree, adj_start, adj_list, adj_edge_id, &
                                    member_start, member_list, &
                                    edge_used_mask, anchor_branch_id, branch_axis, branch_offset, tmp_path_buf, &
                                    tmp_branch_pts, tmp_branch_g, tmp_branch_perm, tmp_branch_lstack, tmp_branch_rstack, &
                                    edge_from, edge_to, n_edges)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of the coordinates
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(in) :: point_to_anchor(n_points)
            !! Original point index -> compact anchor index (0 if not an anchor)
        real(real64),   intent(in) :: tangent_bases(dim, manifold_dim, n_points)
            !! Per-anchor tangent basis directions (only meaningful for anchor points)
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        integer(int32), intent(in) :: primary_anchor_ids(n_points), secondary_anchor_ids(n_points)
            !! Point index of the highest/second-highest weighted stitching anchor
            !! (Step 10); 0 if not covered by that many anchors
        real(real64),   intent(in) :: skeleton_coords(dim+1, n_points)
            !! Stitched point positions (row 1 = anchor_count, rows 2:dim+1 = position)
        integer(int32), intent(in) :: max_edges
            !! Size of the caller-provided edge_from/edge_to buffers
        integer(int32), intent(in) :: n_mst
            !! Number of meaningful entries in the anchor MST edge list
        integer(int32), intent(in) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST
        integer(int32), intent(in) :: adj_start(n_anchor + 1)
            !! CSR: anchor -> range of its MST adjacency-list entries
        integer(int32), intent(in) :: adj_list(:), adj_edge_id(:)
            !! CSR: neighboring anchor / MST-edge-id for each adjacency entry
        integer(int32), intent(in) :: member_start(n_anchor + 1)
            !! CSR: anchor -> range of its member-chain entries
        integer(int32), intent(in) :: member_list(:)
            !! CSR: point index for each member-chain entry

        logical,        intent(out) :: edge_used_mask(:)
            !! Per-MST-edge walked-already flag
        integer(int32), intent(out) :: anchor_branch_id(n_anchor)
            !! Branch id each anchor was last threaded into
        real(real64),   intent(out) :: branch_axis(dim, n_anchor)
            !! Per-anchor tangent axis, oriented forward along its branch
        real(real64),   intent(out) :: branch_offset(n_anchor)
            !! Per-anchor running distance-along-branch offset
        integer(int32), intent(out) :: tmp_path_buf(:)
            !! Sequence of compact anchor indices along the branch being walked
        integer(int32), intent(inout) :: tmp_branch_pts(:), tmp_branch_perm(:)
            !! Point-index / sort-permutation scratch buffers
        integer(int32), intent(inout) :: tmp_branch_lstack(:), tmp_branch_rstack(:)
            !! sort_array recursion-stack scratch buffers
        real(real64),   intent(inout) :: tmp_branch_g(:)
            !! Per-point ordering-coordinate scratch buffer
        integer(int32), intent(inout) :: edge_from(max_edges), edge_to(max_edges)
            !! Backbone edge list, as point-index pairs; each branch's edges appended here
        integer(int32), intent(inout) :: n_edges
            !! Number of meaningful entries in edge_from/edge_to

        integer(int32) :: s, ei, e, cur, prev, next_anchor, ej, e_next
        integer(int32) :: path_len, current_branch

        edge_used_mask = .false.
        anchor_branch_id = 0

        current_branch = 0
        do s = 1, n_anchor
            if (anchor_degree(s) == 2) cycle   ! only walk out from endpoints/junctions
            ! anchor_degree(s) == 0 (isolated anchor): adj_start(s) == adj_start(s+1),
            ! so the loop below naturally does nothing -- its own cluster was
            ! already threaded above, nothing more to connect it to.

            do ei = adj_start(s), adj_start(s + 1) - 1
                e = adj_edge_id(ei)
                if (edge_used_mask(e)) cycle
                edge_used_mask(e) = .true.
                current_branch = current_branch + 1

                path_len = 1
                tmp_path_buf(1) = s
                prev = s
                cur = adj_list(ei)
                do while (anchor_degree(cur) == 2)
                    path_len = path_len + 1
                    tmp_path_buf(path_len) = cur
                    ! cur has exactly two neighbors; step to the one that isn't `prev`.
                    do ej = adj_start(cur), adj_start(cur + 1) - 1
                        if (adj_list(ej) /= prev) then
                            e_next = adj_edge_id(ej)
                            next_anchor = adj_list(ej)
                            exit
                        end if
                    end do
                    edge_used_mask(e_next) = .true.
                    prev = cur
                    cur = next_anchor
                end do
                path_len = path_len + 1
                tmp_path_buf(path_len) = cur   ! cur is now the branch's other special anchor

                call emit_branch(tmp_path_buf(1:path_len), path_len, current_branch, &
                                 n_points, dim, manifold_dim, n_anchor, &
                                 anchor_to_point, point_to_anchor, tangent_bases, anchor_centers, &
                                 member_start, member_list, primary_anchor_ids, secondary_anchor_ids, &
                                 skeleton_coords, max_edges, &
                                 branch_axis, branch_offset, anchor_branch_id, &
                                 tmp_branch_pts, tmp_branch_g, tmp_branch_perm, tmp_branch_lstack, tmp_branch_rstack, &
                                 edge_from, edge_to, n_edges)
            end do
        end do

    end subroutine build_skeleton_edges

    !> Builds the index mapping between the compact 1..n_anchor anchor
    !| numbering and the original 1..n_points point numbering, in both
    !| directions.
    pure subroutine build_anchor_mapping(n_points, is_anchor_mask, n_anchor, anchor_to_point, point_to_anchor)

        integer(int32), intent(in)  :: n_points
            !! Number of points in the dataset
        logical,        intent(in)  :: is_anchor_mask(n_points)
            !! .true. for points selected as atlas anchors
        integer(int32), intent(in)  :: n_anchor
            !! Number of atlas anchors, i.e. count(is_anchor_mask)
        integer(int32), intent(out) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(out) :: point_to_anchor(n_points)
            !! Original point index -> compact anchor index (0 if not an anchor)

        integer(int32) :: i, k

        point_to_anchor = 0
        k = 0
        do i = 1, n_points
            if (is_anchor_mask(i)) then
                k = k + 1
                anchor_to_point(k) = i
                point_to_anchor(i) = k
            end if
        end do

    end subroutine build_anchor_mapping

    !> Tier 1 candidate discovery: flags an anchor pair as a backbone
    !| candidate whenever some point already blends between them as its
    !| primary and secondary anchor (Step 10) -- the same pair driving that
    !| point's stitched position, so it is the natural proxy for "these two
    !| anchors belong to the same branch".
    pure subroutine mark_tier1_candidate_pairs(n_points, n_anchor, point_to_anchor, &
                                          primary_anchor_ids, secondary_anchor_ids, pair_seen_mask)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: point_to_anchor(n_points)
            !! Original point index -> compact anchor index (0 if not an anchor)
        integer(int32), intent(in) :: primary_anchor_ids(n_points), secondary_anchor_ids(n_points)
            !! Point index of the highest/second-highest weighted stitching anchor
            !! (Step 10); 0 if not covered by that many anchors

        logical, intent(out) :: pair_seen_mask(n_anchor, n_anchor)
            !! .true. for anchor pairs (lo, hi) flagged as tier-1 candidates

        integer(int32) :: i, a1, a2, lo, hi

        pair_seen_mask = .false.
        do i = 1, n_points
            if (primary_anchor_ids(i) <= 0 .or. secondary_anchor_ids(i) <= 0) cycle
            a1 = point_to_anchor(primary_anchor_ids(i))
            a2 = point_to_anchor(secondary_anchor_ids(i))
            if (a1 <= 0 .or. a2 <= 0 .or. a1 == a2) cycle
            lo = min(a1, a2) ; hi = max(a1, a2)
            pair_seen_mask(lo, hi) = .true.
        end do

    end subroutine mark_tier1_candidate_pairs

    !> Tier 2 candidate sizing: counts anchor pairs not already flagged by
    !| tier 1 whose spheres overlap (the same geometric test
    !| build_tier2_mst_edges later uses to fill the candidate list), so the
    !| caller can allocate the tier-2 candidate arrays to the exact size.
    pure subroutine count_tier2_candidate_pairs(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                           sphere_radii, pair_seen_mask, n_cand2)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of anchor_centers
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-anchor sphere radius
        logical,        intent(in) :: pair_seen_mask(n_anchor, n_anchor)
            !! .true. for anchor pairs already flagged by tier 1

        integer(int32), intent(out) :: n_cand2
            !! Number of tier-2 candidate edges

        integer(int32) :: lo, hi
        real(real64)   :: dist_c

        n_cand2 = 0
        do hi = 2, n_anchor
            do lo = 1, hi - 1
                if (pair_seen_mask(lo, hi)) cycle
                dist_c = sqrt(sum((anchor_centers(:, anchor_to_point(hi)) - &
                                   anchor_centers(:, anchor_to_point(lo)))**2))
                if (dist_c <= sphere_radii(anchor_to_point(lo)) + sphere_radii(anchor_to_point(hi))) &
                    n_cand2 = n_cand2 + 1
            end do
        end do

    end subroutine count_tier2_candidate_pairs

    !> Allocates the candidate anchor-anchor edge buffers (tier 1: from Step
    !| 10's own primary/secondary pairing; tier 2: geometric sphere-overlap
    !| fallback for thin overlaps tier 1 misses) and delegates to
    !| build_anchor_mst, which runs Kruskal + union-find over them and
    !| classifies anchor roles. No further processing happens here.
    subroutine build_anchor_mst_alloc(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                      sphere_radii, point_to_anchor, primary_anchor_ids, &
                                      secondary_anchor_ids, mst_a1, mst_a2, n_mst, anchor_degree, &
                                      anchor_role, ierr)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of anchor_centers
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-anchor sphere radius
        integer(int32), intent(in) :: point_to_anchor(n_points)
            !! Original point index -> compact anchor index (0 if not an anchor)
        integer(int32), intent(in) :: primary_anchor_ids(n_points), secondary_anchor_ids(n_points)
            !! Point index of the highest/second-highest weighted stitching anchor
            !! (Step 10); 0 if not covered by that many anchors

        integer(int32), intent(out) :: mst_a1(n_anchor), mst_a2(n_anchor)
            !! MST edge list, as compact anchor-index pairs
        integer(int32), intent(out) :: n_mst
            !! Number of meaningful entries in mst_a1/mst_a2
        integer(int32), intent(out) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST
        integer(int32), intent(out) :: anchor_role(n_points)
            !! 0 = not an anchor; 1 = endpoint (backbone degree <= 1);
            !! 2 = pass-through (degree 2); 3 = branch/junction (degree >= 3)
        integer(int32), intent(out) :: ierr
            !! Error code: 0 = success

        integer(int32) :: n_cand, n_cand2
        logical,        allocatable :: pair_seen_mask(:,:)
        integer(int32), allocatable :: tmp_cand_a1(:), tmp_cand_a2(:), tmp_cand_perm(:)
        integer(int32), allocatable :: tmp_cand_lstack(:), tmp_cand_rstack(:)
        real(real64),   allocatable :: tmp_cand_weight(:)
        integer(int32), allocatable :: tmp_cand2_a1(:), tmp_cand2_a2(:), tmp_cand2_perm(:)
        integer(int32), allocatable :: tmp_cand2_lstack(:), tmp_cand2_rstack(:)
        real(real64),   allocatable :: tmp_cand2_weight(:)
        integer(int32), allocatable :: tmp_uf_parent(:)

        call set_ok(ierr)

        M_ALLOCATE(pair_seen_mask(n_anchor, n_anchor))
        call mark_tier1_candidate_pairs(n_points, n_anchor, point_to_anchor, &
                                        primary_anchor_ids, secondary_anchor_ids, pair_seen_mask)

        n_cand = count(pair_seen_mask)
        M_ALLOCATE(tmp_cand_a1(max(1,n_cand)))
        M_ALLOCATE(tmp_cand_a2(max(1,n_cand)))
        M_ALLOCATE(tmp_cand_weight(max(1,n_cand)))
        M_ALLOCATE(tmp_cand_perm(max(1,n_cand)))
        M_ALLOCATE(tmp_cand_lstack(max(1,n_cand)))
        M_ALLOCATE(tmp_cand_rstack(max(1,n_cand)))

        call count_tier2_candidate_pairs(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                         sphere_radii, pair_seen_mask, n_cand2)
        M_ALLOCATE(tmp_cand2_a1(max(1,n_cand2)))
        M_ALLOCATE(tmp_cand2_a2(max(1,n_cand2)))
        M_ALLOCATE(tmp_cand2_weight(max(1,n_cand2)))
        M_ALLOCATE(tmp_cand2_perm(max(1,n_cand2)))
        M_ALLOCATE(tmp_cand2_lstack(max(1,n_cand2)))
        M_ALLOCATE(tmp_cand2_rstack(max(1,n_cand2)))

        M_ALLOCATE(tmp_uf_parent(n_anchor))

        call build_anchor_mst(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                              sphere_radii, pair_seen_mask, n_cand, n_cand2, &
                              tmp_cand_a1, tmp_cand_a2, tmp_cand_weight, tmp_cand_perm, tmp_cand_lstack, tmp_cand_rstack, &
                              tmp_cand2_a1, tmp_cand2_a2, tmp_cand2_weight, tmp_cand2_perm, tmp_cand2_lstack, tmp_cand2_rstack, &
                              tmp_uf_parent, mst_a1, mst_a2, n_mst, anchor_degree, anchor_role)

    end subroutine build_anchor_mst_alloc

    !> Runs Kruskal + union-find over the pre-counted, pre-allocated tier-1
    !| and tier-2 candidate edges (see build_anchor_mst_alloc) to build the
    !| anchor MST, then classifies each anchor's structural role from its
    !| degree in the resulting tree. No allocation.
    pure subroutine build_anchor_mst(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                sphere_radii, pair_seen_mask, n_cand, n_cand2, &
                                tmp_cand_a1, tmp_cand_a2, tmp_cand_weight, tmp_cand_perm, tmp_cand_lstack, tmp_cand_rstack, &
                                tmp_cand2_a1, tmp_cand2_a2, tmp_cand2_weight, tmp_cand2_perm, tmp_cand2_lstack, tmp_cand2_rstack, &
                                tmp_uf_parent, mst_a1, mst_a2, n_mst, anchor_degree, anchor_role)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of anchor_centers
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-anchor sphere radius
        logical,        intent(in) :: pair_seen_mask(n_anchor, n_anchor)
            !! .true. for anchor pairs already flagged by tier 1
        integer(int32), intent(in) :: n_cand
            !! Number of tier-1 candidate edges (0 if none)
        integer(int32), intent(in) :: n_cand2
            !! Number of tier-2 candidate edges (0 if none)

        integer(int32), intent(inout) :: tmp_cand_a1(:), tmp_cand_a2(:), tmp_cand_perm(:)
            !! Tier-1 candidate edges (anchor-index pairs) and sort-permutation scratch
        integer(int32), intent(inout) :: tmp_cand_lstack(:), tmp_cand_rstack(:)
            !! sort_array recursion-stack scratch buffers
        real(real64),   intent(inout) :: tmp_cand_weight(:)
            !! Tier-1 candidate edge weight (anchor-center distance)
        integer(int32), intent(inout) :: tmp_cand2_a1(:), tmp_cand2_a2(:), tmp_cand2_perm(:)
            !! Tier-2 candidate edges (anchor-index pairs) and sort-permutation scratch
        integer(int32), intent(inout) :: tmp_cand2_lstack(:), tmp_cand2_rstack(:)
            !! sort_array recursion-stack scratch buffers
        real(real64),   intent(inout) :: tmp_cand2_weight(:)
            !! Tier-2 candidate edge weight (anchor-center distance)
        integer(int32), intent(inout) :: tmp_uf_parent(n_anchor)
            !! Union-find parent array

        integer(int32), intent(out) :: mst_a1(n_anchor), mst_a2(n_anchor)
            !! MST edge list, as compact anchor-index pairs
        integer(int32), intent(out) :: n_mst
            !! Number of meaningful entries in mst_a1/mst_a2
        integer(int32), intent(out) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST
        integer(int32), intent(out) :: anchor_role(n_points)
            !! 0 = not an anchor; 1 = endpoint (backbone degree <= 1);
            !! 2 = pass-through (degree 2); 3 = branch/junction (degree >= 3)

        integer(int32) :: k

        do k = 1, n_anchor ; tmp_uf_parent(k) = k ; end do
        anchor_degree = 0
        n_mst = 0

        call build_tier1_mst_edges(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                   pair_seen_mask, n_cand, tmp_cand_a1, tmp_cand_a2, tmp_cand_weight, &
                                   tmp_cand_perm, tmp_cand_lstack, tmp_cand_rstack, &
                                   tmp_uf_parent, mst_a1, mst_a2, n_mst, anchor_degree)

        call build_tier2_mst_edges(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                   sphere_radii, pair_seen_mask, n_cand2, tmp_cand2_a1, tmp_cand2_a2, tmp_cand2_weight, &
                                   tmp_cand2_perm, tmp_cand2_lstack, tmp_cand2_rstack, &
                                   tmp_uf_parent, mst_a1, mst_a2, n_mst, anchor_degree)

        call classify_anchor_roles(n_points, n_anchor, anchor_to_point, anchor_degree, anchor_role)

    end subroutine build_anchor_mst

    !> Tier 1 (topological): Kruskal + union-find over the candidate edges
    !| that come from Step 10's own primary/secondary anchor pairing. Keeps
    !| the shortest connecting edge per pair of components and drops
    !| redundant overlaps, so the backbone is acyclic by construction and no
    !| loop-detection is ever needed downstream.
    pure subroutine build_tier1_mst_edges(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                     pair_seen_mask, n_cand, tmp_cand_a1, tmp_cand_a2, tmp_cand_weight, &
                                     tmp_cand_perm, tmp_cand_lstack, tmp_cand_rstack, &
                                     tmp_uf_parent, mst_a1, mst_a2, n_mst, anchor_degree)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of anchor_centers
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        logical,        intent(in) :: pair_seen_mask(n_anchor, n_anchor)
            !! .true. for anchor pairs already flagged by tier 1
        integer(int32), intent(in) :: n_cand
            !! Number of tier-1 candidate edges (0 if none)

        integer(int32), intent(inout) :: tmp_cand_a1(:), tmp_cand_a2(:), tmp_cand_perm(:)
            !! Tier-1 candidate edges (anchor-index pairs) and sort-permutation scratch
        integer(int32), intent(inout) :: tmp_cand_lstack(:), tmp_cand_rstack(:)
            !! sort_array recursion-stack scratch buffers
        real(real64),   intent(inout) :: tmp_cand_weight(:)
            !! Tier-1 candidate edge weight (anchor-center distance)
        integer(int32), intent(inout) :: tmp_uf_parent(n_anchor)
            !! Union-find parent array

        integer(int32), intent(inout) :: mst_a1(n_anchor), mst_a2(n_anchor)
            !! MST edge list, as compact anchor-index pairs
        integer(int32), intent(inout) :: n_mst
            !! Number of meaningful entries in mst_a1/mst_a2
        integer(int32), intent(inout) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST

        integer(int32) :: a1, a2, lo, hi, e_idx, root1, root2, idx

        if (n_cand <= 0) return

        e_idx = 0
        do hi = 1, n_anchor
            do lo = 1, hi - 1
                if (.not. pair_seen_mask(lo, hi)) cycle
                e_idx = e_idx + 1
                tmp_cand_a1(e_idx) = lo
                tmp_cand_a2(e_idx) = hi
                tmp_cand_weight(e_idx) = sqrt(sum((anchor_centers(:, anchor_to_point(hi)) - &
                                               anchor_centers(:, anchor_to_point(lo)))**2))
            end do
        end do

        do e_idx = 1, n_cand ; tmp_cand_perm(e_idx) = e_idx ; end do
        call sort_array(tmp_cand_weight(1:n_cand), tmp_cand_perm(1:n_cand), tmp_cand_lstack(1:n_cand), tmp_cand_rstack(1:n_cand))

        do e_idx = 1, n_cand
            idx = tmp_cand_perm(e_idx)
            a1 = tmp_cand_a1(idx) ; a2 = tmp_cand_a2(idx)
            call find_root(tmp_uf_parent, a1, root1)
            call find_root(tmp_uf_parent, a2, root2)
            if (root1 /= root2) then
                tmp_uf_parent(root1) = root2
                anchor_degree(a1) = anchor_degree(a1) + 1
                anchor_degree(a2) = anchor_degree(a2) + 1
                n_mst = n_mst + 1
                mst_a1(n_mst) = a1
                mst_a2(n_mst) = a2
            end if
        end do

    end subroutine build_tier1_mst_edges

    !> Tier 2 (geometric refinement): continuing the SAME union-find (started
    !| by build_tier1_mst_edges), bridges whatever is still disconnected
    !| using sphere overlap. This only ever fires for anchor pairs tier 1
    !| left in separate components (a pair already joined by a tier-1 edge
    !| is skipped below via root1 == root2), which is what happens when two
    !| neighboring anchors overlap thinly enough that no single point ranked
    !| them as its own top-2 -- leaving a visible gap in the backbone
    !| otherwise.
    pure subroutine build_tier2_mst_edges(n_points, dim, n_anchor, anchor_to_point, anchor_centers, &
                                     sphere_radii, pair_seen_mask, n_cand2, tmp_cand2_a1, tmp_cand2_a2, tmp_cand2_weight, &
                                     tmp_cand2_perm, tmp_cand2_lstack, tmp_cand2_rstack, &
                                     tmp_uf_parent, mst_a1, mst_a2, n_mst, anchor_degree)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of anchor_centers
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: sphere_radii(n_points)
            !! Per-anchor sphere radius
        logical,        intent(in) :: pair_seen_mask(n_anchor, n_anchor)
            !! .true. for anchor pairs already flagged by tier 1
        integer(int32), intent(in) :: n_cand2
            !! Number of tier-2 candidate edges (0 if none)

        integer(int32), intent(inout) :: tmp_cand2_a1(:), tmp_cand2_a2(:), tmp_cand2_perm(:)
            !! Tier-2 candidate edges (anchor-index pairs) and sort-permutation scratch
        integer(int32), intent(inout) :: tmp_cand2_lstack(:), tmp_cand2_rstack(:)
            !! sort_array recursion-stack scratch buffers
        real(real64),   intent(inout) :: tmp_cand2_weight(:)
            !! Tier-2 candidate edge weight (anchor-center distance)
        integer(int32), intent(inout) :: tmp_uf_parent(n_anchor)
            !! Union-find parent array

        integer(int32), intent(inout) :: mst_a1(n_anchor), mst_a2(n_anchor)
            !! MST edge list, as compact anchor-index pairs
        integer(int32), intent(inout) :: n_mst
            !! Number of meaningful entries in mst_a1/mst_a2
        integer(int32), intent(inout) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST

        integer(int32) :: a1, a2, lo, hi, e_idx, root1, root2, idx
        real(real64)   :: dist_c

        if (n_cand2 <= 0) return

        e_idx = 0
        do hi = 2, n_anchor
            do lo = 1, hi - 1
                if (pair_seen_mask(lo, hi)) cycle
                dist_c = sqrt(sum((anchor_centers(:, anchor_to_point(hi)) - &
                                   anchor_centers(:, anchor_to_point(lo)))**2))
                if (dist_c <= sphere_radii(anchor_to_point(lo)) + sphere_radii(anchor_to_point(hi))) then
                    e_idx = e_idx + 1
                    tmp_cand2_a1(e_idx) = lo
                    tmp_cand2_a2(e_idx) = hi
                    tmp_cand2_weight(e_idx) = dist_c
                end if
            end do
        end do

        do e_idx = 1, n_cand2 ; tmp_cand2_perm(e_idx) = e_idx ; end do
        call sort_array(tmp_cand2_weight(1:n_cand2), tmp_cand2_perm(1:n_cand2), tmp_cand2_lstack(1:n_cand2), tmp_cand2_rstack(1:n_cand2))

        do e_idx = 1, n_cand2
            idx = tmp_cand2_perm(e_idx)
            a1 = tmp_cand2_a1(idx) ; a2 = tmp_cand2_a2(idx)
            call find_root(tmp_uf_parent, a1, root1)
            call find_root(tmp_uf_parent, a2, root2)
            if (root1 /= root2) then
                tmp_uf_parent(root1) = root2
                anchor_degree(a1) = anchor_degree(a1) + 1
                anchor_degree(a2) = anchor_degree(a2) + 1
                n_mst = n_mst + 1
                mst_a1(n_mst) = a1
                mst_a2(n_mst) = a2
            end if
        end do

    end subroutine build_tier2_mst_edges

    !> Classifies each anchor's structural role from its degree in the MST.
    pure subroutine classify_anchor_roles(n_points, n_anchor, anchor_to_point, anchor_degree, anchor_role)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(in) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST

        integer(int32), intent(out) :: anchor_role(n_points)
            !! 0 = not an anchor; 1 = endpoint (backbone degree <= 1);
            !! 2 = pass-through (degree 2); 3 = branch/junction (degree >= 3)

        integer(int32) :: i, k

        anchor_role = 0
        do k = 1, n_anchor
            i = anchor_to_point(k)
            if (anchor_degree(k) <= 1) then
                anchor_role(i) = 1
            else if (anchor_degree(k) == 2) then
                anchor_role(i) = 2
            else
                anchor_role(i) = 3
            end if
        end do

    end subroutine classify_anchor_roles

    !> Threads every point into the member chain of its primary (best-fit)
    !| anchor (CSR: member_start/member_list), then locally orders and
    !| emits edges for each "special" (non-pass-through) anchor's own
    !| cluster -- a hub is threaded exactly once here, regardless of how
    !| many branches later connect through it.
    pure subroutine build_member_chains(n_points, dim, manifold_dim, n_anchor, anchor_to_point, &
                                   point_to_anchor, primary_anchor_ids, anchor_degree, &
                                   skeleton_coords, anchor_centers, tangent_bases, max_edges, &
                                   tmp_branch_pts, tmp_branch_g, tmp_branch_perm, tmp_branch_lstack, tmp_branch_rstack, &
                                   member_start, member_list, edge_from, edge_to, n_edges)

        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of the coordinates
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(in) :: point_to_anchor(n_points)
            !! Original point index -> compact anchor index (0 if not an anchor)
        integer(int32), intent(in) :: primary_anchor_ids(n_points)
            !! Point index of the highest-weighted stitching anchor covering this point
        integer(int32), intent(in) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST
        real(real64),   intent(in) :: skeleton_coords(dim+1, n_points)
            !! Stitched point positions (row 1 = anchor_count, rows 2:dim+1 = position)
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        real(real64),   intent(in) :: tangent_bases(dim, manifold_dim, n_points)
            !! Per-anchor tangent basis directions (only meaningful for anchor points)
        integer(int32), intent(in) :: max_edges
            !! Size of the caller-provided edge_from/edge_to buffers

        ! Scratch, sized n_points by the caller; also reused later by emit_branch.
        integer(int32), intent(inout) :: tmp_branch_pts(:), tmp_branch_perm(:)
            !! Point-index / sort-permutation scratch buffers
        integer(int32), intent(inout) :: tmp_branch_lstack(:), tmp_branch_rstack(:)
            !! sort_array recursion-stack scratch buffers
        real(real64),   intent(inout) :: tmp_branch_g(:)
            !! Per-point ordering-coordinate scratch buffer

        integer(int32), intent(out) :: member_start(n_anchor + 1)
            !! CSR: anchor -> range of its member-chain entries
        integer(int32), intent(out) :: member_list(n_points)
            !! CSR: point index for each member-chain entry (caller-sized upper
            !! bound: at most n_points points are ever threaded in total)
        integer(int32), intent(inout) :: edge_from(max_edges), edge_to(max_edges)
            !! Backbone edge list, as point-index pairs; hub-cluster edges appended here
        integer(int32), intent(inout) :: n_edges
            !! Number of meaningful entries in edge_from/edge_to

        integer(int32) :: i, k, a1, idx, mstart_s, mend_s, cnt_s, pidx
        integer(int32) :: member_count(n_anchor)

        ! --- Per-anchor member chains: every point is threaded into exactly ONE
        ! chain, that of its primary (best-fit) anchor -- not primary AND
        ! secondary. A point in the overlap of two anchors still only "belongs"
        ! to a single position along the backbone; assigning it to both chains
        ! would order it twice, independently, against two different tangent
        ! axes that need not agree, producing duplicated/overlapping fragments
        ! right at chart boundaries instead of one continuous line. Each point
        ! is threaded once, and the MST bridge below is what stitches
        ! consecutive anchors' chains together.
        member_count = 0
        do i = 1, n_points
            if (primary_anchor_ids(i) > 0) then
                a1 = point_to_anchor(primary_anchor_ids(i))
                if (a1 > 0) member_count(a1) = member_count(a1) + 1
            end if
        end do

        member_start(1) = 1
        do k = 1, n_anchor
            member_start(k + 1) = member_start(k) + member_count(k)
        end do

        member_count = 0   ! reuse as per-anchor fill offset
        do i = 1, n_points
            if (primary_anchor_ids(i) > 0) then
                a1 = point_to_anchor(primary_anchor_ids(i))
                if (a1 > 0) then
                    member_count(a1) = member_count(a1) + 1
                    member_list(member_start(a1) + member_count(a1) - 1) = i
                end if
            end if
        end do

        ! Scratch buffers shared by the special-anchor chains below and by
        ! emit_branch (branch interiors): sized to n_points, the safe upper
        ! bound for either a single hub's cluster or one branch's pooled points.

        ! --- Special-anchor (endpoint/junction) local chains: a "special" anchor
        ! (anchor_degree /= 2) is a hub that can be shared by several branches
        ! below (a junction is literally the meeting point of 3+ of them). Its
        ! own member cluster is threaded into a simple chain here, ONCE, using
        ! its own tangent axis -- a small local cluster doesn't need a
        ! branch-wide coordinate. Branches only pool INTERIOR (pass-through)
        ! anchors' points, so a hub's cluster is never re-threaded once per
        ! incident branch.
        do k = 1, n_anchor
            if (anchor_degree(k) == 2) cycle
            mstart_s = member_start(k) ; mend_s = member_start(k + 1) - 1
            cnt_s = mend_s - mstart_s + 1
            if (cnt_s <= 1) cycle
            pidx = anchor_to_point(k)
            do idx = 1, cnt_s
                tmp_branch_g(idx)    = dot_product(skeleton_coords(2:dim+1, member_list(mstart_s + idx - 1)) - &
                                                anchor_centers(:, pidx), tangent_bases(:, 1, pidx))
                tmp_branch_perm(idx) = idx
            end do
            call sort_array(tmp_branch_g(1:cnt_s), tmp_branch_perm(1:cnt_s), tmp_branch_lstack(1:cnt_s), tmp_branch_rstack(1:cnt_s))
            tmp_branch_lstack(1:cnt_s) = member_list(mstart_s:mend_s)
            do idx = 1, cnt_s
                member_list(mstart_s + idx - 1) = tmp_branch_lstack(tmp_branch_perm(idx))
            end do
            do idx = mstart_s, mend_s - 1
                if (n_edges >= max_edges) exit
                n_edges = n_edges + 1
                edge_from(n_edges) = member_list(idx)
                edge_to(n_edges)   = member_list(idx + 1)
            end do
        end do

    end subroutine build_member_chains

    !> CSR adjacency of the anchor MST, so each branch can be walked from
    !| one special anchor (endpoint or junction) to the next in
    !| O(branch length) instead of re-scanning all MST edges.
    pure subroutine build_branch_adjacency(n_anchor, mst_a1, mst_a2, n_mst, anchor_degree, &
                                      adj_start, adj_list, adj_edge_id)

        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: mst_a1(n_anchor), mst_a2(n_anchor)
            !! MST edge list, as compact anchor-index pairs
        integer(int32), intent(in) :: n_mst
            !! Number of meaningful entries in mst_a1/mst_a2
        integer(int32), intent(in) :: anchor_degree(n_anchor)
            !! Degree of each anchor in the MST

        integer(int32), intent(out) :: adj_start(n_anchor + 1)
            !! CSR: anchor -> range of its adjacency-list entries
        integer(int32), intent(out) :: adj_list(:)
            !! CSR: adjacent anchor index for each entry (caller-sized to 2*n_mst)
        integer(int32), intent(out) :: adj_edge_id(:)
            !! CSR: MST edge index (into mst_a1/mst_a2) for each entry (caller-sized to 2*n_mst)

        integer(int32) :: k, e_idx, a1, a2
        integer(int32) :: fill_off(n_anchor)

        ! --- CSR adjacency of the MST, so each branch can be walked from one
        ! special anchor (endpoint or junction) to the next in O(branch length).
        adj_start(1) = 1
        do k = 1, n_anchor
            adj_start(k + 1) = adj_start(k) + anchor_degree(k)
        end do
        fill_off = 0
        do e_idx = 1, n_mst
            a1 = mst_a1(e_idx) ; a2 = mst_a2(e_idx)
            fill_off(a1) = fill_off(a1) + 1
            adj_list(adj_start(a1) + fill_off(a1) - 1)    = a2
            adj_edge_id(adj_start(a1) + fill_off(a1) - 1) = e_idx
            fill_off(a2) = fill_off(a2) + 1
            adj_list(adj_start(a2) + fill_off(a2) - 1)    = a1
            adj_edge_id(adj_start(a2) + fill_off(a2) - 1) = e_idx
        end do

    end subroutine build_branch_adjacency

    !> Union-find root lookup with path compression, for the Kruskal MST in
    !| build_anchor_mst. A subroutine rather than a function since a pure
    !| function's arguments must all be intent(in), but path compression
    !| needs to mutate parent in place.
    pure subroutine find_root(parent, x, r)
        integer(int32), intent(inout) :: parent(:)
            !! Union-find parent array; path-compressed in place
        integer(int32), intent(in)    :: x
            !! Node whose root to find
        integer(int32), intent(out)   :: r
            !! Root of x's component
        integer(int32) :: cur, nxt
        cur = x
        do while (parent(cur) /= cur)
            cur = parent(cur)
        end do
        r = cur
        cur = x
        do while (parent(cur) /= r)
            nxt = parent(cur)
            parent(cur) = r
            cur = nxt
        end do
    end subroutine find_root

    !> Closest (Euclidean) point among anchor k's own members to `target`;
    !| falls back to the anchor's own point if it has no members at all.
    pure integer(int32) function nearest_member(k, target, n_points, dim, n_anchor, &
                                           anchor_to_point, member_start, member_list, &
                                           skeleton_coords) result(best)
        integer(int32), intent(in) :: k
            !! Compact index of the anchor whose members to search
        real(real64),   intent(in) :: target(dim)
            !! Point to find the nearest member to
        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of the coordinates
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(in) :: member_start(n_anchor + 1)
            !! CSR: anchor -> range of its member-chain entries
        integer(int32), intent(in) :: member_list(:)
            !! CSR: point index for each member-chain entry
        real(real64),   intent(in) :: skeleton_coords(dim+1, n_points)
            !! Stitched point positions (row 1 = anchor_count, rows 2:dim+1 = position)
        integer(int32) :: jj, mstart3, mend3
        real(real64)   :: d2, best_d2
        best = anchor_to_point(k)
        best_d2 = huge(1.0_real64)
        mstart3 = member_start(k) ; mend3 = member_start(k + 1) - 1
        do jj = mstart3, mend3
            d2 = sum((skeleton_coords(2:dim+1, member_list(jj)) - target)**2)
            if (d2 < best_d2) then
                best_d2 = d2
                best = member_list(jj)
            end if
        end do
    end function nearest_member

    !> path always has path_len >= 2: path(1) and path(path_len) are the
    !| branch's two "special" (endpoint/junction) anchors, already
    !| threaded into their own chain above; only the INTERIOR
    !| (pass-through) anchors path(2:path_len-1), if any, are pooled and
    !| chained here, then bridged into each special anchor's nearest
    !| member -- so a hub's cluster is never re-threaded once per branch.
    pure subroutine emit_branch(path, path_len, branch_id, n_points, dim, manifold_dim, n_anchor, &
                           anchor_to_point, point_to_anchor, tangent_bases, anchor_centers, &
                           member_start, member_list, primary_anchor_ids, secondary_anchor_ids, &
                           skeleton_coords, max_edges, &
                           branch_axis, branch_offset, anchor_branch_id, &
                           tmp_branch_pts, tmp_branch_g, tmp_branch_perm, tmp_branch_lstack, tmp_branch_rstack, &
                           edge_from, edge_to, n_edges)

        integer(int32), intent(in) :: path(:)
            !! Sequence of compact anchor indices along this branch
        integer(int32), intent(in) :: path_len
            !! Number of meaningful entries in path
        integer(int32), intent(in) :: branch_id
            !! Identifier assigned to this branch by the caller's walk loop
        integer(int32), intent(in) :: n_points
            !! Number of points in the dataset
        integer(int32), intent(in) :: dim
            !! Ambient dimension of the coordinates
        integer(int32), intent(in) :: manifold_dim
            !! Target intrinsic dimension of the local tangent subspaces
        integer(int32), intent(in) :: n_anchor
            !! Number of atlas anchors
        integer(int32), intent(in) :: anchor_to_point(n_anchor)
            !! Compact anchor index (1..n_anchor) -> original point index
        integer(int32), intent(in) :: point_to_anchor(n_points)
            !! Original point index -> compact anchor index (0 if not an anchor)
        real(real64),   intent(in) :: tangent_bases(dim, manifold_dim, n_points)
            !! Per-anchor tangent basis directions (only meaningful for anchor points)
        real(real64),   intent(in) :: anchor_centers(dim, n_points)
            !! Centroid of each anchor's sphere (only meaningful for anchor points)
        integer(int32), intent(in) :: member_start(n_anchor + 1)
            !! CSR: anchor -> range of its member-chain entries
        integer(int32), intent(in) :: member_list(:)
            !! CSR: point index for each member-chain entry
        integer(int32), intent(in) :: primary_anchor_ids(n_points), secondary_anchor_ids(n_points)
            !! Point index of the highest/second-highest weighted stitching anchor
            !! (Step 10); 0 if not covered by that many anchors
        real(real64),   intent(in) :: skeleton_coords(dim+1, n_points)
            !! Stitched point positions (row 1 = anchor_count, rows 2:dim+1 = position)
        integer(int32), intent(in) :: max_edges
            !! Size of the caller-provided edge_from/edge_to buffers

        real(real64),   intent(inout) :: branch_axis(dim, n_anchor)
            !! Per-anchor tangent axis, oriented forward along its branch
        real(real64),   intent(inout) :: branch_offset(n_anchor)
            !! Per-anchor running distance-along-branch offset
        integer(int32), intent(inout) :: anchor_branch_id(n_anchor)
            !! Branch id each anchor was last threaded into
        integer(int32), intent(inout) :: tmp_branch_pts(:), tmp_branch_perm(:)
            !! Point-index / sort-permutation scratch buffers
        integer(int32), intent(inout) :: tmp_branch_lstack(:), tmp_branch_rstack(:)
            !! sort_array recursion-stack scratch buffers
        real(real64),   intent(inout) :: tmp_branch_g(:)
            !! Per-point ordering-coordinate scratch buffer
        integer(int32), intent(inout) :: edge_from(max_edges), edge_to(max_edges)
            !! Backbone edge list, as point-index pairs; this branch's edges appended here
        integer(int32), intent(inout) :: n_edges
            !! Number of meaningful entries in edge_from/edge_to

        integer(int32) :: ii, kk, pidx, jj, mstart2, mend2, a1b, a2b, npts, i
        integer(int32) :: s_anchor, e_anchor, first_pt, last_pt
        real(real64)   :: axis_i(dim), dir(dim), step_dist, gi, gj

        ! Orient each anchor's own tangent axis to point "forward" along the
        ! branch, and accumulate a running offset (plain anchor-to-anchor
        ! distance) so anchor k's local coordinate continues where the
        ! previous anchor's left off, instead of restarting at zero. Computed
        ! for every anchor in path (including the two special ends), since
        ! interior points need a continuous coordinate relative to them.
        do ii = 1, path_len
            kk = path(ii)
            pidx = anchor_to_point(kk)
            axis_i = tangent_bases(:, 1, pidx)
            if (ii < path_len) then
                dir = anchor_centers(:, anchor_to_point(path(ii + 1))) - anchor_centers(:, pidx)
            else
                dir = anchor_centers(:, pidx) - anchor_centers(:, anchor_to_point(path(ii - 1)))
            end if
            if (dot_product(axis_i, dir) < 0.0_real64) axis_i = -axis_i
            branch_axis(:, kk) = axis_i

            if (ii == 1) then
                branch_offset(kk) = 0.0_real64
            else
                step_dist = sqrt(sum((anchor_centers(:, pidx) - &
                                      anchor_centers(:, anchor_to_point(path(ii - 1))))**2))
                branch_offset(kk) = branch_offset(path(ii - 1)) + step_dist
            end if
            anchor_branch_id(kk) = branch_id
        end do

        s_anchor = path(1)
        e_anchor = path(path_len)

        ! Pool ONLY interior anchors' points -- the two special ends' own
        ! clusters were already threaded once, above.
        npts = 0
        do ii = 2, path_len - 1
            kk = path(ii)
            mstart2 = member_start(kk) ; mend2 = member_start(kk + 1) - 1
            do jj = mstart2, mend2
                npts = npts + 1
                tmp_branch_pts(npts) = member_list(jj)
            end do
        end do

        if (npts == 0) then
            ! No interior anchors between the two hubs: bridge them directly.
            if (n_edges < max_edges) then
                n_edges = n_edges + 1
                edge_from(n_edges) = nearest_member(s_anchor, anchor_centers(:, anchor_to_point(e_anchor)), &
                                                    n_points, dim, n_anchor, anchor_to_point, &
                                                    member_start, member_list, skeleton_coords)
                edge_to(n_edges)   = nearest_member(e_anchor, anchor_centers(:, anchor_to_point(s_anchor)), &
                                                    n_points, dim, n_anchor, anchor_to_point, &
                                                    member_start, member_list, skeleton_coords)
            end if
            return
        end if

        do jj = 1, npts
            i   = tmp_branch_pts(jj)
            a1b = point_to_anchor(primary_anchor_ids(i))
            gi  = branch_offset(a1b) + dot_product(skeleton_coords(2:dim+1, i) - &
                  anchor_centers(:, anchor_to_point(a1b)), branch_axis(:, a1b))

            ! If this point also blends into a secondary anchor that is on
            ! the SAME branch (the usual case right at a chart boundary),
            ! average both anchors' coordinates instead of using only the
            ! primary one -- this is what actually smooths the transition,
            ! the same way Step 10 already blends the point's position.
            if (secondary_anchor_ids(i) > 0) then
                a2b = point_to_anchor(secondary_anchor_ids(i))
                if (a2b > 0) then
                    if (anchor_branch_id(a2b) == branch_id) then
                        gj = branch_offset(a2b) + dot_product(skeleton_coords(2:dim+1, i) - &
                             anchor_centers(:, anchor_to_point(a2b)), branch_axis(:, a2b))
                        gi = 0.5_real64 * (gi + gj)
                    end if
                end if
            end if
            tmp_branch_g(jj)    = gi
            tmp_branch_perm(jj) = jj
        end do

        call sort_array(tmp_branch_g(1:npts), tmp_branch_perm(1:npts), tmp_branch_lstack(1:npts), tmp_branch_rstack(1:npts))

        tmp_branch_lstack(1:npts) = tmp_branch_pts(1:npts)
        do jj = 1, npts
            tmp_branch_pts(jj) = tmp_branch_lstack(tmp_branch_perm(jj))
        end do

        do jj = 1, npts - 1
            if (n_edges >= max_edges) exit
            n_edges = n_edges + 1
            edge_from(n_edges) = tmp_branch_pts(jj)
            edge_to(n_edges)   = tmp_branch_pts(jj + 1)
        end do

        ! Bridge each end of the interior chain into its special anchor's
        ! nearest member, so the branch actually connects into both hubs.
        first_pt = tmp_branch_pts(1)
        last_pt  = tmp_branch_pts(npts)
        if (n_edges < max_edges) then
            n_edges = n_edges + 1
            edge_from(n_edges) = nearest_member(s_anchor, skeleton_coords(2:dim+1, first_pt), &
                                                n_points, dim, n_anchor, anchor_to_point, &
                                                member_start, member_list, skeleton_coords)
            edge_to(n_edges)   = first_pt
        end if
        if (n_edges < max_edges) then
            n_edges = n_edges + 1
            edge_from(n_edges) = last_pt
            edge_to(n_edges)   = nearest_member(e_anchor, skeleton_coords(2:dim+1, last_pt), &
                                                n_points, dim, n_anchor, anchor_to_point, &
                                                member_start, member_list, skeleton_coords)
        end if
    end subroutine emit_branch


end module lomanle_mod