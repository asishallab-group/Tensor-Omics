module mod_test_smoothing_metrics
  use asserts
  use f42_utils, only: compute_roughness, compute_rmse, get_approx_diameter, compute_coverage, compute_smoothing_score
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

    ! ===========================================
    ! ROUGHNESS 
    ! ===========================================

    !> Constant field: y = C => roughness = 0
    subroutine test_constant_field()
        integer(int32), parameter :: n_points = 5, n_dims = 3, k_neighbors = 2
        real(real64) :: y(n_dims, n_points), sigma(n_points), roughness
        integer(int32) :: neighbors(k_neighbors, n_points)
        real(real64) :: distances(k_neighbors, n_points)

        y = 3.14159_real64
        sigma = 1.0_real64
        neighbors = reshape([2,3, 1,3, 1,2, 1,2, 1,2], [k_neighbors, n_points])
        distances = 1.0_real64

        call compute_roughness(y, sigma, neighbors, distances, n_points, n_dims, k_neighbors, roughness)
        call assert_equal_real(roughness, 0.0_real64, 1d-12, "Invariant: Constant field must be 0")
    end subroutine test_constant_field

    !> Translation invariance: y2 = y + offset => roughness(y2) == roughness(y)
    subroutine test_translation_invariance()
        integer(int32), parameter :: n_points = 3, n_dims = 2, k_neighbors = 2
        real(real64) :: y1(n_dims, n_points), y2(n_dims, n_points)
        real(real64) :: sigma(n_points), dist(k_neighbors, n_points), r1, r2
        integer(int32) :: neigh(k_neighbors, n_points)

        y1 = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [n_dims, n_points])
        y2 = y1 + 100.0_real64
        sigma = 1.5_real64
        dist = 0.5_real64
        neigh = reshape([2,3, 1,3, 1,2], [k_neighbors, n_points])

        call compute_roughness(y1, sigma, neigh, dist, n_points, n_dims, k_neighbors, r1)
        call compute_roughness(y2, sigma, neigh, dist, n_points, n_dims, k_neighbors, r2)
        call assert_equal_real(r1, r2, 1d-12, "Invariant: Translation does not affect roughness")
    end subroutine test_translation_invariance

    !> Quadratic scaling: y2 = c*y => roughness(y2) = c^2 * roughness(y)
    subroutine test_quadratic_scaling()
        integer(int32), parameter :: n_points = 3, n_dims = 1, k_neighbors = 2
        real(real64) :: y(n_dims, n_points), c, r1, r2
        real(real64) :: sigma(n_points), dist(k_neighbors, n_points)
        integer(int32) :: neigh(k_neighbors, n_points)

        c = 2.0_real64
        y(:,1)=1.0; y(:,2)=4.0; y(:,3)=9.0
        sigma = 1.0; dist = 1.0
        neigh = reshape([2,3, 1,3, 1,2], [k_neighbors, n_points])

        call compute_roughness(y, sigma, neigh, dist, n_points, n_dims, k_neighbors, r1)
        call compute_roughness(y*c, sigma, neigh, dist, n_points, n_dims, k_neighbors, r2)
        call assert_equal_real(r2, (c**2)*r1, 1d-12, "Invariant: Quadratic scaling with field amplitude")
    end subroutine test_quadratic_scaling

    !> 2 points analytical check: R = (y1-y2)^2
    subroutine test_two_points_analytical()
        integer(int32), parameter :: n_points = 2, n_dims = 1, k_neighbors = 1
        real(real64) :: y(n_dims, n_points), sigma(n_points), dist(k_neighbors, n_points), r
        integer(int32) :: neigh(k_neighbors, n_points)
        real(real64) :: expected

        y(:,1) = 10.0; y(:,2) = 20.0
        sigma = [1.0, 2.0]
        dist = 0.5
        neigh(:,1) = [2]; neigh(:,2) = [1]
        
        ! For 2 mutual points, normalization Z cancels weights R, leaving only diff2
        expected = (y(1,1) - y(1,2))**2
        call compute_roughness(y, sigma, neigh, dist, n_points, n_dims, k_neighbors, r)
        call assert_equal_real(r, expected, 1d-12, "Analytical: 2-point mutual neighbors")
    end subroutine test_two_points_analytical

    !> Distance = 0 check: weights become 1.0, roughness = simple mean of diff2
    subroutine test_zero_distance_mean()
        integer(int32), parameter :: n_points = 3, n_dims = 1, k_neighbors = 2
        real(real64) :: y(n_dims, n_points), r
        real(real64) :: sigma(n_points), dist(k_neighbors, n_points)
        integer(int32) :: neigh(k_neighbors, n_points)
        real(real64) :: d12, d13, d23, expected

        y(:,1)=1.0; y(:,2)=2.0; y(:,3)=5.0
        dist = 0.0_real64 ! exp(0) = 1.0 weights
        sigma = 1.0
        neigh = reshape([2,3, 1,3, 1,2], [k_neighbors, n_points])

        d12 = (1.0-2.0)**2; d13 = (1.0-5.0)**2; d23 = (2.0-5.0)**2
        ! Total R = 2*(d12 + d13 + d23), Total Z = 2*(1+1+1) = 6
        expected = (d12 + d13 + d12 + d23 + d13 + d23) / 6.0_real64
        
        call compute_roughness(y, sigma, neigh, dist, n_points, n_dims, k_neighbors, r)
        call assert_equal_real(r, expected, 1d-12, "Analytical: Zero distance gives arithmetic mean")
    end subroutine test_zero_distance_mean

    !> Invalid neighbors: Check if indices like 0 or > n_points are skipped
    subroutine test_invalid_indices()
        integer(int32), parameter :: n_points = 2, n_dims = 1, k_neighbors = 2
        real(real64) :: y(n_dims, n_points), sigma(n_points)
        real(real64) :: dist_dirty(k_neighbors, n_points), dist_clean(k_neighbors, n_points)
        real(real64) :: r_dirty, r_clean
        integer(int32) :: neigh_dirty(k_neighbors, n_points), neigh_clean(k_neighbors, n_points)

        y(:,1) = 0.0_real64; y(:,2) = 10.0_real64
        sigma = 1.0_real64

        ! Dirty: one valid neighbor + one invalid
        neigh_dirty(:,1) = [2_int32, 999_int32]
        neigh_dirty(:,2) = [1_int32, -5_int32]
        dist_dirty(:,1)  = [0.1_real64, 0.1_real64]   ! second entry ignored anyway
        dist_dirty(:,2)  = [0.1_real64, 0.1_real64]

        ! Clean: replace invalid neighbors by duplicating the valid neighbor
        ! (This keeps the same number of valid contributions as dirty case.)
        neigh_clean(:,1) = [2_int32, 2_int32]
        neigh_clean(:,2) = [1_int32, 1_int32]
        dist_clean = dist_dirty

        call compute_roughness(y, sigma, neigh_dirty, dist_dirty, n_points, n_dims, k_neighbors, r_dirty)
        call compute_roughness(y, sigma, neigh_clean, dist_clean, n_points, n_dims, k_neighbors, r_clean)

        call assert_equal_real(r_dirty, r_clean, 1d-12, "Invalid indices must be ignored (match clean equivalent)")
    end subroutine test_invalid_indices


    !> Self-loops: If neighbor is the point itself, it must be ignored
    subroutine test_self_loop_exclusion()
        integer(int32), parameter :: n_points = 2, n_dims = 1, k_neighbors = 1
        real(real64) :: y(n_dims, n_points), sigma(n_points)
        real(real64) :: dist(k_neighbors, n_points)
        real(real64) :: r_valid, r_self
        integer(int32) :: neigh_valid(k_neighbors, n_points), neigh_self(k_neighbors, n_points)

        y(:,1)=1.0_real64; y(:,2)=5.0_real64
        sigma = 1.0_real64
        dist  = 0.5_real64

        ! Valid mutual neighbors
        neigh_valid(:,1) = [2_int32]
        neigh_valid(:,2) = [1_int32]

        ! Self-loops only -> should produce Z=0 -> roughness = 0
        neigh_self(:,1) = [1_int32]
        neigh_self(:,2) = [2_int32]

        call compute_roughness(y, sigma, neigh_valid, dist, n_points, n_dims, k_neighbors, r_valid)
        call compute_roughness(y, sigma, neigh_self,  dist, n_points, n_dims, k_neighbors, r_self)

        call assert_true(r_valid > 0.0_real64, "Sanity: valid neighbors must give positive roughness here")
        call assert_equal_real(r_self, 0.0_real64, 1d-12, "Self-loops must be ignored -> no contributions -> roughness 0")
    end subroutine test_self_loop_exclusion


    !> All invalid: If all weights are 0 or neighbors invalid, result must be 0
    subroutine test_all_invalid_results_zero()
        integer(int32), parameter :: n_points = 2, n_dims = 1, k_neighbors = 1
        real(real64) :: y(n_dims, n_points), sigma(n_points), dist(k_neighbors, n_points), r
        integer(int32) :: neigh(k_neighbors, n_points)

        y = 1.0; dist = 1.0
        sigma = 0.0_real64 ! All sigmas invalid
        neigh(:,1) = [2]; neigh(:,2) = [1]

        call compute_roughness(y, sigma, neigh, dist, n_points, n_dims, k_neighbors, r)
        call assert_equal_real(r, 0.0_real64, 1d-12, "Robustness: Total invalid input gives 0 roughness")
    end subroutine test_all_invalid_results_zero

    subroutine test_nonnegativity_and_finite()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        integer(int32), parameter :: n_points = 3, n_dims = 2, k_neighbors = 2
        real(real64) :: y(n_dims, n_points), sigma(n_points), dist(k_neighbors, n_points), r
        integer(int32) :: neigh(k_neighbors, n_points)

        y = reshape([ 1.0_real64, 2.0_real64, &
                    3.0_real64, 4.0_real64, &
                    5.0_real64, 6.0_real64 ], [n_dims, n_points])
        sigma = 1.0_real64
        dist = 0.7_real64
        neigh = reshape([2,3, 1,3, 1,2], [k_neighbors, n_points])

        call compute_roughness(y, sigma, neigh, dist, n_points, n_dims, k_neighbors, r)

        call assert_true(.not. ieee_is_nan(r), "Roughness must not be NaN")
        call assert_true(r >= 0.0_real64, "Roughness must be non-negative")
    end subroutine test_nonnegativity_and_finite

    subroutine test_mixed_sigma_exclusion()
        integer(int32), parameter :: n_points = 3, n_dims = 1, k_neighbors = 2
        real(real64) :: y(n_dims, n_points), sigma_dirty(n_points), sigma_clean(n_points)
        real(real64) :: dist(k_neighbors, n_points)
        integer(int32) :: neigh_dirty(k_neighbors, n_points), neigh_clean(k_neighbors, n_points)
        real(real64) :: r_dirty, r_clean

        y(:,1)=1.0_real64; y(:,2)=10.0_real64; y(:,3)=2.0_real64
        dist = 0.3_real64

        ! Full graph (complete among {1,2,3})
        neigh_dirty = reshape([2,3, 1,3, 1,2], [k_neighbors, n_points])

        ! Sigma dirty: point 2 invalid
        sigma_dirty = [1.0_real64, 0.0_real64, 1.0_real64]

        ! Clean reference: remove any contributions involving point 2 by forcing self-loops for its row,
        ! and also not referencing it as a neighbor for others.
        neigh_clean(:,1) = [3_int32, 3_int32]   ! 1 connects only to 3
        neigh_clean(:,2) = [2_int32, 2_int32]   ! self only
        neigh_clean(:,3) = [1_int32, 1_int32]   ! 3 connects only to 1
        sigma_clean = [1.0_real64, 1.0_real64, 1.0_real64]  ! all valid now

        call compute_roughness(y, sigma_dirty, neigh_dirty, dist, n_points, n_dims, k_neighbors, r_dirty)
        call compute_roughness(y, sigma_clean, neigh_clean, dist, n_points, n_dims, k_neighbors, r_clean)

        call assert_equal_real(r_dirty, r_clean, 1d-12, "Sigma<=eps must exclude point contributions (match clean pruned graph)")
    end subroutine test_mixed_sigma_exclusion

    subroutine test_edge_duplication_invariance()
        integer(int32), parameter :: n_points = 3, n_dims = 1
        integer(int32), parameter :: k1 = 1, k2 = 2
        real(real64) :: y(n_dims, n_points), sigma(n_points)
        real(real64) :: dist1(k1, n_points), dist2(k2, n_points)
        integer(int32) :: neigh1(k1, n_points), neigh2(k2, n_points)
        real(real64) :: r1, r2

        y(:,1)=1.0_real64; y(:,2)=2.0_real64; y(:,3)=5.0_real64
        sigma = 1.0_real64

        ! k=1: a simple directed ring 1->2, 2->3, 3->1
        neigh1(:,1) = [2_int32]
        neigh1(:,2) = [3_int32]
        neigh1(:,3) = [1_int32]
        dist1 = 0.8_real64

        ! k=2: duplicate the same neighbor twice
        neigh2(:,1) = [2_int32, 2_int32]
        neigh2(:,2) = [3_int32, 3_int32]
        neigh2(:,3) = [1_int32, 1_int32]
        dist2(:,1) = [0.8_real64, 0.8_real64]
        dist2(:,2) = [0.8_real64, 0.8_real64]
        dist2(:,3) = [0.8_real64, 0.8_real64]

        call compute_roughness(y, sigma, neigh1, dist1, n_points, n_dims, k1, r1)
        call compute_roughness(y, sigma, neigh2, dist2, n_points, n_dims, k2, r2)

        call assert_equal_real(r1, r2, 1d-12, "Metamorphic: duplicating edges must not change roughness")
    end subroutine test_edge_duplication_invariance

    subroutine test_permutation_invariance()
        integer(int32), parameter :: n_points = 3, n_dims = 1, k_neighbors = 2
        real(real64) :: y(n_dims, n_points), yP(n_dims, n_points)
        real(real64) :: sigma(n_points), sigmaP(n_points)
        real(real64) :: dist(k_neighbors, n_points), distP(k_neighbors, n_points)
        integer(int32) :: neigh(k_neighbors, n_points), neighP(k_neighbors, n_points)
        real(real64) :: r, rP
        integer(int32) :: P(n_points), Pinv(n_points)
        integer(int32) :: i, j

        ! Original
        y(:,1)=1.0_real64; y(:,2)=2.0_real64; y(:,3)=5.0_real64
        sigma = [1.0_real64, 2.0_real64, 1.5_real64]
        neigh = reshape([2,3, 1,3, 1,2], [k_neighbors, n_points])
        dist  = reshape([0.4_real64, 0.9_real64, &
                        0.4_real64, 0.7_real64, &
                        0.9_real64, 0.7_real64], [k_neighbors, n_points])

        ! Permutation P: new_index i corresponds to old_index P(i)
        ! Example: swap 1 and 3 => P = [3,2,1]
        P    = [3_int32, 2_int32, 1_int32]
        Pinv = 0_int32
        do i=1_int32,n_points
            Pinv(P(i)) = i
        end do

        ! Apply permutation to y and sigma: new arrays are ordered by old P(i)
        do i=1_int32,n_points
            yP(:,i) = y(:,P(i))
            sigmaP(i) = sigma(P(i))
        end do

        ! Rebuild neighbors/distances under permutation:
        ! If original edge is i -> neigh(j,i), in permuted space it becomes Pinv(i_old) -> Pinv(neigh_old)
        do i=1_int32,n_points
            do j=1_int32,k_neighbors
                neighP(j, Pinv(i)) = Pinv( neigh(j,i) )
                distP(j, Pinv(i))  = dist(j,i)
            end do
        end do

        call compute_roughness(y,  sigma,  neigh,  dist,  n_points, n_dims, k_neighbors, r)
        call compute_roughness(yP, sigmaP, neighP, distP, n_points, n_dims, k_neighbors, rP)

        call assert_equal_real(r, rP, 1d-12, "Permutation: consistent reindexing must preserve roughness")
    end subroutine test_permutation_invariance

    ! ===========================================
    ! RMSE 
    ! ===========================================
    !> Identity test: y == y0 => RMSE = 0
    subroutine test_rmse_identity()
        integer(int32), parameter :: n_points = 5, n_dims = 2
        real(real64) :: y0(n_dims, n_points), y(n_dims, n_points), r
        
        y0 = 1.5_real64
        y  = y0
        call compute_rmse(y0, y, n_points, n_dims, r)
        call assert_equal_real(r, 0.0_real64,  1d-12, "Identity: Identical fields must yield 0.0")
    end subroutine test_rmse_identity

    !> Analytical check: 1 point, 2D error (3, 4) => RMSE = 5
    subroutine test_rmse_pythagoras()
        integer(int32), parameter :: n_points = 1, n_dims = 2
        real(real64) :: y0(2,1), y(2,1), r
        
        y0(:,1) = [0.0_real64, 0.0_real64]
        y(:,1)  = [3.0_real64, 4.0_real64] 
        ! S = 3^2 + 4^2 = 25. MSE = 25/1. RMSE = 5.
        call compute_rmse(y0, y, n_points, n_dims, r)
        call assert_equal_real(r, 5.0_real64, 1d-12, "Analytical: 2D Pythagorean error check")
    end subroutine test_rmse_pythagoras

    !> Joint translation invariance: RMSE(y+a, y0+a) == RMSE(y, y0)
    subroutine test_rmse_joint_translation()
        integer(int32), parameter :: n_points = 3, n_dims = 1
        real(real64) :: y0(1,3), y(1,3), offset(1,3), r1, r2
        
        y0(1,:) = [1.0, 2.0, 3.0]
        y(1,:)  = [1.1, 2.1, 3.1]
        offset  = 100.0_real64
        
        call compute_rmse(y0, y, n_points, 1_int32, r1)
        call compute_rmse(y0 + offset, y + offset, n_points, 1_int32, r2)
        
        call assert_equal_real(r1, r2, 1d-12, "Invariant: Joint translation must not change RMSE")
    end subroutine test_rmse_joint_translation

    !> Linear scaling: RMSE(c*y, c*y0) = c * RMSE(y, y0)
    subroutine test_rmse_scaling()
        integer(int32), parameter :: n_points = 2, n_dims = 1
        real(real64) :: y0(1,2), y(1,2), r1, r2, c
        
        c = 2.5_real64
        y0 = 0.0; y = 1.0
        
        call compute_rmse(y0, y, n_points, 1_int32, r1)
        call compute_rmse(y0*c, y*c, n_points, 1_int32, r2)
        
        call assert_equal_real(r2, r1 * c, 1d-12, "Invariant: RMSE must scale linearly with error magnitude")
    end subroutine test_rmse_scaling

    !> Robustness: Ensure result is always non-negative and finite
    subroutine test_rmse_non_negativity()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
        integer(int32), parameter :: n_points = 10, n_dims = 3
        real(real64) :: y0(n_dims, n_points), y(n_dims, n_points), r
        
        call random_number(y0)
        call random_number(y)
        
        call compute_rmse(y0, y, n_points, n_dims, r)
        
        call assert_true(r >= 0.0_real64, "Robustness: RMSE must be non-negative")
        call assert_true(.not. ieee_is_nan(r), "Robustness: RMSE must be finite (not NaN)")
    end subroutine test_rmse_non_negativity

    !> Safety: n_points = 0 should return 0.0
    subroutine test_rmse_zero_points()
        integer(int32), parameter :: n_points = 0
        real(real64) :: y0(1,1), y(1,1), r
        
        ! Note: Arrays are allocated with size 1 to avoid compiler-specific issues
        ! with zero-sized arrays, but we pass n_points = 0.
        call compute_rmse(y0, y, n_points, 1_int32, r)
        call assert_equal_real(r, 0.0_real64, 1d-12, "Safety: Zero points must yield 0.0 RMSE")
    end subroutine test_rmse_zero_points

    ! ==========================================================================
    ! DIAMETER HEURISTIC 
    ! ==========================================================================

    !> Edge case: n < 2 => diameter = 0 (anchors = 1,1)
    subroutine test_diameter_n_lt_2()
        integer(int32), parameter :: dim = 2
        real(real64) :: d
        integer(int32) :: i1, i2
        real(real64) :: field0(dim, 0)
        real(real64) :: field1(dim, 1)

        ! n = 0
        call get_approx_diameter(field0, 0_int32, dim, d, i1, i2)
        call assert_equal_real(d, 0.0_real64, 1d-12, "Diameter n=0 must be 0.0")
        call assert_equal_int(i1, 1_int32, "Diameter n=0: i1 must be 1")
        call assert_equal_int(i2, 1_int32, "Diameter n=0: i2 must be 1")

        ! n = 1
        field1(:,1) = [1.0_real64, 2.0_real64]
        call get_approx_diameter(field1, 1_int32, dim, d, i1, i2)
        call assert_equal_real(d, 0.0_real64, 1d-12, "Diameter n=1 must be 0.0")
        call assert_equal_int(i1, 1_int32, "Diameter n=1: i1 must be 1")
        call assert_equal_int(i2, 1_int32, "Diameter n=1: i2 must be 1")
    end subroutine test_diameter_n_lt_2

    !> Simple 1D line: points at 0 and 10 => diameter 10 (exact)
    subroutine test_diameter_1d()
        integer(int32), parameter :: n = 3, dim = 1
        real(real64) :: field(dim, n), d
        integer(int32) :: i1, i2

        field(1,:) = [0.0_real64, 5.0_real64, 10.0_real64]
        call get_approx_diameter(field, n, dim, d, i1, i2)
        call assert_equal_real(d, 10.0_real64, 1d-12, "Diameter 1D: distance between 0 and 10")

        ! anchors must hit the extremes for this deterministic case
        call assert_equal_int(i1, 3_int32, "Diameter 1D: i1 should be farthest from index 1 (point at 10)")
        call assert_equal_int(i2, 1_int32, "Diameter 1D: i2 should be farthest from i1 (point at 0)")
    end subroutine test_diameter_1d

    !> 2D unit square: diameter is diagonal sqrt(2) (exact)
    subroutine test_diameter_2d_square()
        integer(int32), parameter :: n = 4, dim = 2
        real(real64) :: field(dim, n), d
        integer(int32) :: i1, i2

        field(:,1) = [0.0_real64, 0.0_real64]
        field(:,2) = [1.0_real64, 0.0_real64]
        field(:,3) = [0.0_real64, 1.0_real64]
        field(:,4) = [1.0_real64, 1.0_real64]

        call get_approx_diameter(field, n, dim, d, i1, i2)
        call assert_equal_real(d, sqrt(2.0_real64), 1d-12, "Diameter 2D: diagonal of unit square")

        ! Because of ties, anchors can be any opposite corners. Just validate the distance.
        call assert_true(i1 >= 1_int32 .and. i1 <= n, "Diameter 2D: i1 in range")
        call assert_true(i2 >= 1_int32 .and. i2 <= n, "Diameter 2D: i2 in range")
    end subroutine test_diameter_2d_square

    !> 3D deterministic case (no ties) where 2-sweep from index 1 must find the true diameter.
    !! Construction:
    !!   p1=(0,0,0)
    !!   p2=(10,0,0)
    !!   p3=(0,1,0)
    !! From p1, farthest is p2 (unique). From p2, farthest is p3 (unique).
    !! True diameter between p2 and p3: sqrt((10)^2 + (1)^2) = sqrt(101).
    subroutine test_diameter_3d_deterministic()
        integer(int32), parameter :: n = 3, dim = 3
        real(real64) :: field(dim, n), d
        real(real64) :: expected
        integer(int32) :: i1, i2

        field(:,1) = [0.0_real64, 0.0_real64, 0.0_real64]   ! p1
        field(:,2) = [10.0_real64, 0.0_real64, 0.0_real64]  ! p2
        field(:,3) = [0.0_real64, 1.0_real64, 0.0_real64]   ! p3

        expected = sqrt(101.0_real64)
        call get_approx_diameter(field, n, dim, d, i1, i2)
        call assert_equal_real(d, expected, 1d-12, "Diameter 3D: deterministic 2-sweep must match sqrt(101)")

        call assert_equal_int(i1, 2_int32, "Diameter 3D: i1 must be p2 (farthest from p1)")
        call assert_equal_int(i2, 3_int32, "Diameter 3D: i2 must be p3 (farthest from p2)")
    end subroutine test_diameter_3d_deterministic


    ! ==========================================================================
    ! COVERAGE (FROZEN ANCHORS)
    ! ==========================================================================

    !> Normal ratio: D0=10, Dt=5 => coverage=0.5, penalty=1/(0.5+eps)
    subroutine test_coverage_ratio()
        real(real64) :: d0, y(1,2), cov, pen
        integer(int32), parameter :: n=2, dim=1
        integer(int32) :: i1, i2
        real(real64), parameter :: eps = 1.0e-12_real64

        d0 = 10.0_real64
        i1 = 1_int32
        i2 = 2_int32

        y(:,1) = [0.0_real64]
        y(:,2) = [5.0_real64]   ! dt = 5

        call compute_coverage(d0, i1, i2, y, n, dim, cov, pen)
        call assert_equal_real(cov, 0.5_real64, 1d-12, "Coverage: ratio 5/10 must be 0.5")
        call assert_equal_real(pen, 1.0_real64/(0.5_real64 + eps), 1d-12, "Coverage: penalty must be 1/(cov+eps)")
    end subroutine test_coverage_ratio

    !> Translation invariance: shifting all points does not change dt, hence same coverage
    subroutine test_coverage_translation_invariance()
        real(real64) :: d0, y1(1,2), y2(1,2), cov1, pen1, cov2, pen2
        integer(int32), parameter :: n=2, dim=1
        integer(int32) :: i1, i2

        d0 = 1.0_real64
        i1 = 1_int32
        i2 = 2_int32

        y1(:,1) = [0.0_real64]
        y1(:,2) = [1.0_real64]          ! dt=1
        y2(:,1) = [100.0_real64]
        y2(:,2) = [101.0_real64]        ! dt=1

        call compute_coverage(d0, i1, i2, y1, n, dim, cov1, pen1)
        call compute_coverage(d0, i1, i2, y2, n, dim, cov2, pen2)

        call assert_equal_real(cov1, cov2, 1d-12, "Coverage: translation must not change coverage")
        call assert_equal_real(pen1, pen2, 1d-12, "Coverage: translation must not change penalty")
    end subroutine test_coverage_translation_invariance

    !> Degenerate original: d0=0
    !! - If dt=0 => coverage=1
    !! - If dt>0 => coverage=0
    subroutine test_coverage_degenerate_d0()
        real(real64) :: d0, y(1,2), cov, pen
        integer(int32), parameter :: n=2, dim=1
        integer(int32) :: i1, i2

        d0 = 0.0_real64
        i1 = 1_int32
        i2 = 2_int32

        ! dt = 0
        y(:,1) = [0.0_real64]
        y(:,2) = [0.0_real64]
        call compute_coverage(d0, i1, i2, y, n, dim, cov, pen)
        call assert_equal_real(cov, 1.0_real64, 1d-12, "Degenerate: dt=0 and d0=0 => coverage=1")

        ! dt > 0
        y(:,2) = [10.0_real64]
        call compute_coverage(d0, i1, i2, y, n, dim, cov, pen)
        call assert_equal_real(cov, 0.0_real64, 1d-12, "Degenerate: dt>0 and d0=0 => coverage=0")
    end subroutine test_coverage_degenerate_d0

    !> Collapse case: d0>0 but current field collapsed so dt=0 => coverage=0
    subroutine test_coverage_collapse_dt_zero()
        real(real64) :: d0, y(1,3), cov, pen
        integer(int32), parameter :: n=3, dim=1
        integer(int32) :: i1, i2

        d0 = 10.0_real64
        i1 = 1_int32
        i2 = 2_int32

        y(:,1) = [7.0_real64]
        y(:,2) = [7.0_real64]
        y(:,3) = [7.0_real64]   ! dt = 0 between any anchors

        call compute_coverage(d0, i1, i2, y, n, dim, cov, pen)
        call assert_equal_real(cov, 0.0_real64, 1d-12, "Collapse: dt=0 with d0>0 => coverage=0")
        call assert_true(pen > 1.0e6_real64, "Collapse: penalty must be very large when coverage ~ 0")
    end subroutine test_coverage_collapse_dt_zero

    !> Expansion case: dt>d0 => coverage>1 (no clamping)
    subroutine test_coverage_allows_expansion()
        real(real64) :: d0, y(1,2), cov, pen
        integer(int32), parameter :: n=2, dim=1
        integer(int32) :: i1, i2

        d0 = 2.0_real64
        i1 = 1_int32
        i2 = 2_int32

        y(:,1) = [0.0_real64]
        y(:,2) = [5.0_real64]   ! dt = 5 => coverage = 2.5

        call compute_coverage(d0, i1, i2, y, n, dim, cov, pen)
        call assert_true(cov > 1.0_real64, "Expansion: coverage should be > 1 if dt > d0")
    end subroutine test_coverage_allows_expansion

    !> Scale consistency: if y is scaled by c and d0 is scaled by c, coverage must be unchanged
    subroutine test_coverage_scale_consistency()
        real(real64) :: d0, d0s, c
        real(real64) :: y(1,2), ys(1,2), cov1, pen1, cov2, pen2
        integer(int32), parameter :: n=2, dim=1
        integer(int32) :: i1, i2

        i1 = 1_int32
        i2 = 2_int32

        c  = 3.0_real64
        d0 = 4.0_real64
        d0s = d0 * c

        y(:,1) = [0.0_real64]
        y(:,2) = [2.0_real64]   ! dt = 2
        ys = y * c               ! dt scales to 6

        call compute_coverage(d0,  i1, i2, y,  n, dim, cov1, pen1)
        call compute_coverage(d0s, i1, i2, ys, n, dim, cov2, pen2)

        call assert_equal_real(cov1, cov2, 1d-12, "Scale: coverage must be invariant under consistent scaling")
        call assert_equal_real(pen1, pen2, 1d-12, "Scale: penalty must be invariant under consistent scaling")
    end subroutine test_coverage_scale_consistency

    ! =========================================================================
    ! SMOOTHING SCORE
    ! =========================================================================

    !> t0 sanity: Rt=R0, Et=0, pt=1 => score = 1
    subroutine test_smoothing_score_t0_is_one()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 2.5_real64
        R0  = 2.5_real64
        Et  = 0.0_real64
        D0  = 10.0_real64
        pt  = 1.0_real64
        eps = 1.0e-12_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_equal_real(score, 1.0_real64, 1d-12, "Score t0 must be 1 when Rt=R0, Et=0, pt=1")
    end subroutine test_smoothing_score_t0_is_one

    !> Roughness only: halve Rt with others unchanged => score = (0.5)^(1/3)
    subroutine test_smoothing_score_roughness_halved()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score, expected
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 5.0_real64
        R0  = 10.0_real64       ! r_norm = 0.5
        Et  = 0.0_real64
        D0  = 1.0_real64
        pt  = 1.0_real64
        eps = 1e-30_real64       

        expected = (0.5_real64)**(1.0_real64/3.0_real64)

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_equal_real(score, expected, 1d-12, "Score must scale with cube-root of r_norm")
    end subroutine test_smoothing_score_roughness_halved

    !> Fidelity only: Et/D0 = 0.5 => (1+e_norm)=1.5, others 1 => score = (1.5)^(1/3)
    subroutine test_smoothing_score_fidelity_increase()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score, expected
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 7.0_real64
        R0  = 7.0_real64       
        Et  = 5.0_real64
        D0  = 10.0_real64      
        pt  = 1.0_real64
        eps = 1e-30_real64

        expected = (1.5_real64)**(1.0_real64/3.0_real64)

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_equal_real(score, expected, 1d-12, "Score must increase with cube-root of (1+e_norm)")
    end subroutine test_smoothing_score_fidelity_increase

    !> Coverage penalty only: pt=2, others 1 => score = 2^(1/3)
    subroutine test_smoothing_score_coverage_penalty_increase()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score, expected
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 3.0_real64
        R0  = 3.0_real64
        Et  = 0.0_real64
        D0  = 10.0_real64
        pt  = 2.0_real64
        eps = 1e-30_real64

        expected = (2.0_real64)**(1.0_real64/3.0_real64)

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_equal_real(score, expected, 1d-12, "Score must increase with cube-root of pt")
    end subroutine test_smoothing_score_coverage_penalty_increase

    !> Combined known product: r_norm=0.25, (1+e)=2.0, pt=4.0 => product=2.0 => score = 2^(1/3)
    subroutine test_smoothing_score_known_product()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score, expected
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 1.0_real64
        R0  = 4.0_real64
        Et  = 1.0_real64
        D0  = 1.0_real64
        pt  = 4.0_real64
        eps = 1e-30_real64

        expected = (2.0_real64)**(1.0_real64/3.0_real64)

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_equal_real(score, expected, 1d-12, "Score must match combined product in geometric mode")
    end subroutine test_smoothing_score_known_product

    !> Epsilon stability: with eps>0, t0 should still be ~1
    subroutine test_smoothing_score_t0_with_epsilon_close_to_one()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 1.0_real64
        R0  = 1.0_real64
        Et  = 0.0_real64
        D0  = 1.0_real64
        pt  = 1.0_real64
        eps = 1.0e-12_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_equal_real(score, 1.0_real64, 1d-10, "Score should remain ~1 at t0 even with epsilon")
    end subroutine test_smoothing_score_t0_with_epsilon_close_to_one

    !> Monotonicity in Rt: if Rt decreases, score must decrease
    subroutine test_smoothing_score_decreases_when_roughness_decreases()
        real(real64) :: Rt1, Rt2, R0, Et, D0, pt, eps, s1, s2
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 1 ; wr = 1.0; we = 1.0; wc = 1.0
        R0  = 10.0_real64
        Et  = 2.0_real64
        D0  = 10.0_real64
        pt  = 1.2_real64
        eps = 1e-30_real64

        Rt1 = 8.0_real64
        Rt2 = 4.0_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt1, R0, Et, D0, pt, eps, s1)
        call compute_smoothing_score(m_flag, wr, we, wc, Rt2, R0, Et, D0, pt, eps, s2)

        call assert_true(s2 < s1, "Score must decrease when roughness decreases")
    end subroutine test_smoothing_score_decreases_when_roughness_decreases

    !> Monotonicity in Et: if RMSE increases, score must increase
    subroutine test_smoothing_score_increases_when_rmse_increases()
        real(real64) :: Rt, R0, Et1, Et2, D0, pt, eps, s1, s2
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 1 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 5.0_real64
        R0  = 5.0_real64
        D0  = 10.0_real64
        pt  = 1.0_real64
        eps = 1e-30_real64

        Et1 = 0.0_real64
        Et2 = 5.0_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et1, D0, pt, eps, s1)
        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et2, D0, pt, eps, s2)

        call assert_true(s2 > s1, "Score must increase when RMSE increases")
    end subroutine test_smoothing_score_increases_when_rmse_increases

    !> Monotonicity in pt: if coverage penalty increases, score must increase
    subroutine test_smoothing_score_increases_when_penalty_increases()
        real(real64) :: Rt, R0, Et, D0, pt1, pt2, eps, s1, s2
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 1 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 5.0_real64
        R0  = 5.0_real64
        Et  = 1.0_real64
        D0  = 10.0_real64
        eps = 1e-30_real64

        pt1 = 1.0_real64
        pt2 = 3.0_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt1, eps, s1)
        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt2, eps, s2)

        call assert_true(s2 > s1, "Score must increase when pt increases")
    end subroutine test_smoothing_score_increases_when_penalty_increases

    !> Scale invariance for Rt/R0
    subroutine test_smoothing_score_invariant_to_scaling_of_roughness_pair()
        real(real64) :: Rt, R0, Et, D0, pt, eps, c, s1, s2
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 2.0_real64
        R0  = 8.0_real64
        Et  = 1.0_real64
        D0  = 5.0_real64
        pt  = 1.1_real64
        eps = 1e-30_real64
        c   = 7.0_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, s1)
        call compute_smoothing_score(m_flag, wr, we, wc, Rt*c, R0*c, Et, D0, pt, eps, s2)

        call assert_equal_real(s1, s2, 1d-12, "Scaling Rt and R0 equally must not change score")
    end subroutine test_smoothing_score_invariant_to_scaling_of_roughness_pair

    !> Scale invariance for Et/D0
    subroutine test_smoothing_score_invariant_to_scaling_of_rmse_and_diameter()
        real(real64) :: Rt, R0, Et, D0, pt, eps, c, s1, s2
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 1.0_real64
        R0  = 2.0_real64
        Et  = 3.0_real64
        D0  = 6.0_real64
        pt  = 1.0_real64
        eps = 1e-30_real64
        c   = 4.0_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, s1)
        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et*c, D0*c, pt, eps, s2)

        call assert_equal_real(s1, s2, 1d-12, "Scaling Et and D0 equally must not change score")
    end subroutine test_smoothing_score_invariant_to_scaling_of_rmse_and_diameter

    !> Extreme contraction: pt huge, score should be large but finite
    subroutine test_smoothing_score_extreme_penalty_large_but_finite()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 2 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = 1.0_real64
        R0  = 1.0_real64
        Et  = 0.0_real64
        D0  = 1.0_real64
        pt  = 1.0e12_real64
        eps = 1e-30_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_true(score > 1.0e3_real64, "Extreme pt: score should be very large")
        call assert_true(score < 1.0e5_real64, "Extreme pt: cube-root should still be finite")
    end subroutine test_smoothing_score_extreme_penalty_large_but_finite

    !> Robustness: negative inputs should not crash
    subroutine test_smoothing_score_negative_inputs_do_not_crash()
        real(real64) :: Rt, R0, Et, D0, pt, eps, score
        integer      :: m_flag
        real(real64) :: wr, we, wc
        
        m_flag = 1 ; wr = 1.0; we = 1.0; wc = 1.0
        Rt  = -1.0_real64
        R0  =  1.0_real64
        Et  =  0.0_real64
        D0  =  1.0_real64
        pt  =  1.0_real64
        eps =  1e-30_real64

        call compute_smoothing_score(m_flag, wr, we, wc, Rt, R0, Et, D0, pt, eps, score)
        call assert_true(score >= 0.0_real64, "Negative inputs: score must be non-negative (clamped)")
    end subroutine test_smoothing_score_negative_inputs_do_not_crash



    ! ==========================================================================
    ! RUNNER
    ! ==========================================================================

    subroutine get_all_tests(all_tests)
        type(test_case), intent(out) :: all_tests(40)
        all_tests(1)  = test_case("test_constant_field", test_constant_field)
        all_tests(2)  = test_case("test_translation_invariance", test_translation_invariance)
        all_tests(3)  = test_case("test_quadratic_scaling", test_quadratic_scaling)
        all_tests(4)  = test_case("test_two_points_analytical", test_two_points_analytical)
        all_tests(5)  = test_case("test_zero_distance_mean", test_zero_distance_mean)
        all_tests(6)  = test_case("test_invalid_indices", test_invalid_indices)
        all_tests(7)  = test_case("test_self_loop_exclusion", test_self_loop_exclusion)
        all_tests(8)  = test_case("test_all_invalid_results_zero", test_all_invalid_results_zero)
        all_tests(9)  = test_case("test_nonnegativity_and_finite", test_nonnegativity_and_finite)
        all_tests(10) = test_case("test_mixed_sigma_exclusion", test_mixed_sigma_exclusion)
        all_tests(11) = test_case("test_edge_duplication_invariance", test_edge_duplication_invariance)
        all_tests(12) = test_case("test_permutation_invariance", test_permutation_invariance)
        all_tests(13) = test_case("test_rmse_identity", test_rmse_identity)
        all_tests(14) = test_case("test_rmse_pythagoras", test_rmse_pythagoras)
        all_tests(15) = test_case("test_rmse_joint_translation", test_rmse_joint_translation)
        all_tests(16) = test_case("test_rmse_scaling", test_rmse_scaling)
        all_tests(17) = test_case("test_rmse_non_negativity", test_rmse_non_negativity)
        all_tests(18) = test_case("test_rmse_zero_points", test_rmse_zero_points)
        all_tests(19) = test_case("test_diameter_1d", test_diameter_1d)
        all_tests(20) = test_case("test_diameter_2d_square", test_diameter_2d_square)
        all_tests(21) = test_case("test_diameter_3d_deterministic", test_diameter_3d_deterministic)
        all_tests(22) = test_case("test_coverage_ratio", test_coverage_ratio)
        all_tests(23) = test_case("test_coverage_translation_invariance", test_coverage_translation_invariance)
        all_tests(24) = test_case("test_coverage_degenerate_d0", test_coverage_degenerate_d0)
        all_tests(25) = test_case("test_coverage_collapse_dt_zero", test_coverage_collapse_dt_zero)
        all_tests(26) = test_case("test_coverage_allows_expansion", test_coverage_allows_expansion)
        all_tests(27) = test_case("test_coverage_scale_consistency", test_coverage_scale_consistency)
        all_tests(28) = test_case("test_smoothing_score_t0_is_one", test_smoothing_score_t0_is_one)
        all_tests(29) = test_case("test_smoothing_score_roughness_halved", test_smoothing_score_roughness_halved)
        all_tests(30) = test_case("test_smoothing_score_fidelity_increase", test_smoothing_score_fidelity_increase)
        all_tests(31) = test_case("test_smoothing_score_coverage_penalty_increase", test_smoothing_score_coverage_penalty_increase)
        all_tests(32) = test_case("test_smoothing_score_known_product", test_smoothing_score_known_product)
        all_tests(33) = test_case("test_smoothing_score_t0_with_epsilon_close_to_one", test_smoothing_score_t0_with_epsilon_close_to_one)
        all_tests(34) = test_case("test_smoothing_score_decreases_when_roughness_decreases", test_smoothing_score_decreases_when_roughness_decreases)
        all_tests(35) = test_case("test_smoothing_score_increases_when_rmse_increases", test_smoothing_score_increases_when_rmse_increases)
        all_tests(36) = test_case("test_smoothing_score_increases_when_penalty_increases", test_smoothing_score_increases_when_penalty_increases)
        all_tests(37) = test_case("test_smoothing_score_invariant_to_scaling_of_roughness_pair", test_smoothing_score_invariant_to_scaling_of_roughness_pair)
        all_tests(38) = test_case("test_smoothing_score_invariant_to_scaling_of_rmse_and_diameter", test_smoothing_score_invariant_to_scaling_of_rmse_and_diameter)
        all_tests(39) = test_case("test_smoothing_score_extreme_penalty_large_but_finite", test_smoothing_score_extreme_penalty_large_but_finite)
        all_tests(40) = test_case("test_smoothing_score_negative_inputs_do_not_crash", test_smoothing_score_negative_inputs_do_not_crash)

    end subroutine get_all_tests

  !> Run all tests.
  subroutine run_all_tests_smoothing_metrics()
    type(test_case) :: all_tests(40)
    integer(int32) :: i

    call get_all_tests(all_tests)
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All smoothing metrics tests passed successfully."
  end subroutine run_all_tests_smoothing_metrics

  !> Run specific smoothing metrics tests by name.
  subroutine run_named_tests_smoothing_metrics(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(40)
    integer(int32) :: i, j
    logical :: found

    call get_all_tests(all_tests)

    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_smoothing_metrics

end module mod_test_smoothing_metrics
