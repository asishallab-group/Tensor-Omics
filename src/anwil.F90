module anwil
  use safeguard
  use kd_tree, only: build_kd_index, kd_knn_query
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
  use f42_utils, only: sort_real
  implicit none

  ! === TIMING GLOBALS FOR PERFORMANCE ANALYSIS ===
  real(real64), save :: timing_kd_build = 0.0_real64
  real(real64), save :: timing_knn_queries = 0.0_real64
  real(real64), save :: timing_gaussian_calc = 0.0_real64
  integer(int32), save :: timing_query_count = 0

  ! PERFORMANCE NOTES:
  ! - K-d tree performance degrades significantly above ~10-20 dimensions
  ! - For high-dimensional data (>20D), consider using brute force search or LSH
  ! - Memory usage scales as O(n) for tree construction
  ! - Query time: O(log n) for low dimensions, O(n) for high dimensions
  
  ! EDGE CASE HANDLING:
  ! - First/last points: Treated same as any other point, no special boundary conditions
  ! - Duplicate points: K-d tree handles gracefully, may return same point multiple times
  ! - Single point datasets: Returns the point unchanged (k > num_points triggers error)
  ! - Very sparse data: Adaptive sigma handles large distance variations
  
    interface
    subroutine dpotrf(uplo, n, a, lda, info)
      use iso_fortran_env, only: real64, int32
      character(len=1), intent(in) :: uplo
      integer(int32), intent(in)   :: n, lda
      real(real64), intent(inout)  :: a(lda, *)
      integer(int32), intent(out)  :: info
    end subroutine dpotrf

    subroutine dpotri(uplo, n, a, lda, info)
      use iso_fortran_env, only: real64, int32
      character(len=1), intent(in) :: uplo
      integer(int32), intent(in)   :: n, lda
      real(real64), intent(inout)  :: a(lda, *)
      integer(int32), intent(out)  :: info
    end subroutine dpotri
  end interface

contains

  ! === TIMING UTILITY FUNCTIONS ===
  subroutine reset_timing_stats()
    timing_kd_build = 0.0_real64
    timing_knn_queries = 0.0_real64
    timing_gaussian_calc = 0.0_real64
    timing_query_count = 0
  end subroutine reset_timing_stats

  subroutine print_timing_stats(n_points, k_neighbors)
    integer(int32), intent(in) :: n_points, k_neighbors
    real(real64) :: total_time, queries_per_sec, avg_query_time
    
    total_time = timing_kd_build + timing_knn_queries + timing_gaussian_calc
    if (timing_query_count > 0) then
      avg_query_time = timing_knn_queries / real(timing_query_count, real64)
      queries_per_sec = real(timing_query_count, real64) / timing_knn_queries
    else
      avg_query_time = 0.0_real64
      queries_per_sec = 0.0_real64
    end if
    
    write(*, '(A)') "=== FORTRAN INTERNAL TIMING BREAKDOWN ==="
    write(*, '(A, I0, A, I0)') "Dataset: ", n_points, " points, k=", k_neighbors
    write(*, '(A, F10.6, A, F5.1, A)') "1. K-d tree build:   ", timing_kd_build, " sec (", &
                                      100.0 * timing_kd_build / total_time, "%)"
    write(*, '(A, F10.6, A, F5.1, A)') "2. KNN+Gaussian:     ", timing_knn_queries, " sec (", &
                                      100.0 * timing_knn_queries / total_time, "%)"
    write(*, '(A, F10.6, A, F5.1, A)') "3. (unused):         ", timing_gaussian_calc, " sec (", &
                                      100.0 * timing_gaussian_calc / total_time, "%)"
    write(*, '(A, F10.6, A)') "Total Fortran time:  ", total_time, " sec"
    write(*, '(A, I0)') "Total queries made:  ", timing_query_count
    write(*, '(A, F10.6, A)') "Avg query time:      ", avg_query_time * 1000.0, " ms"
    write(*, '(A, F8.1, A)') "Queries per second:  ", queries_per_sec, " queries/sec"
    write(*, '(A)') "========================================="
  end subroutine ! print_timing_stats

  ! Helper function to get current time
  function get_time() result(time_seconds)
    integer(int32) :: count, count_rate
    real(real64) :: time_seconds
    call system_clock(count, count_rate)
    time_seconds = real(count, real64) / real(count_rate, real64)
  end function get_time


! Mode 0: Isotropic smoothing (default, simple and efficient)
! Mode 1: Diagonal anisotropy (uses anisotropy_factor for scaling)
! Mode 2: Full anisotropy (uses local covariance matrix)

subroutine smooth_vectors_gaussian_adaptive( &
  coords, vectors, smoothed, &
  n_coord_dims, n_vector_dims, n_points, k_neighbors, &
  kd_indices, dimension_order, neighbors, distances, &
  workspace, value_buffer, permutation, permutation_distances, left_stack, right_stack, &
  anisotropy_mode, anisotropy_factor, sd_arr, ierr )

  use iso_fortran_env, only: real64, int32
  use tox_errors, only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
  use f42_utils, only: sort_real
  use kd_tree, only: build_kd_index, kd_knn_query
  implicit none

  ! =====================
  ! INPUT
  ! =====================
  integer(int32), intent(in) :: n_coord_dims
  integer(int32), intent(in) :: n_vector_dims
  integer(int32), intent(in) :: n_points
  integer(int32), intent(in) :: k_neighbors  ! Number of desired real neighbors

  real(real64), intent(inout) :: coords(n_coord_dims, n_points)
  real(real64), intent(inout) :: vectors(n_vector_dims, n_points)

  integer(int32), intent(in) :: anisotropy_mode
  real(real64), intent(in) :: anisotropy_factor

  ! =====================
  ! OUTPUT
  ! =====================
  real(real64), intent(out) :: smoothed(n_vector_dims, n_points)
  integer(int32), intent(out) :: ierr
  real(real64), intent(out) :: sd_arr(n_points)

  ! =====================
  ! WORK BUFFERS (MUST BE SIZE k_neighbors + 1 FOR KNN)
  ! I assume that neighbors/distances have been declared with that size outside.
  ! If not, the code will fail. We will assume that the caller passes the correct size.
  ! =====================
  integer(int32), intent(out) :: kd_indices(n_points)
  integer(int32), intent(out) :: dimension_order(n_coord_dims)
  integer(int32), intent(inout) :: neighbors(k_neighbors) 
  real(real64), intent(inout) :: distances(k_neighbors) 

  integer(int32), intent(inout) :: workspace(n_points)
  real(real64), intent(inout) :: value_buffer(n_points)
  integer(int32), intent(inout) :: permutation(n_points)
  integer(int32), intent(inout) :: permutation_distances(k_neighbors)
  integer(int32), intent(inout) :: left_stack(n_points)
  integer(int32), intent(inout) :: right_stack(n_points)
  

  ! =====================
  ! LOCALS
  ! =====================
  integer(int32) :: i, j, idx, k_eff,k_start,k_real_neighbors, i_perm
  integer(int32) :: recursion_stack(3, n_points)
  integer(int32) :: k_search ! k_neighbors
  integer(int32) :: k_eff_neighbors

  real(real64) :: local_sigma, local_median, c
  real(real64) :: local_w, local_wsum
  real(real64) :: local_work(n_vector_dims)

  real(real64) :: delta(n_coord_dims)
  real(real64) :: sigma_parallel, sigma_orth
  real(real64) :: d2

  ! Covariance (mode 2)
  real(real64) :: cov(n_coord_dims, n_coord_dims)
  real(real64) :: cov_inv(n_coord_dims, n_coord_dims)
  real(real64) :: tmp_vec(n_coord_dims)

  real(real64) :: trace_cov, eps_reg, scale_factor
  real(real64) :: Iden(n_coord_dims, n_coord_dims)
  integer(int32), parameter :: real_neighbor_offset = 1 
  integer(int32) :: median_idx 
  
  ! =====================
  call set_ok(ierr)

  ! Debugging: ! print input parameters
  ! print *, "Debug: n_coord_dims=", n_coord_dims, ", n_vector_dims=", n_vector_dims
  ! print *, "Debug: n_points=", n_points, ", k_neighbors=", k_neighbors

  ! Validate that the buffers have the minimum size (k_neighbors + 1)
  if (k_neighbors < 1 .or. k_neighbors > n_points) then
    ! print *, "Error: Invalid k_neighbors or n_points. k_neighbors=", k_neighbors, ", n_points=", n_points
    call set_err_once(ierr, ERR_INVALID_INPUT)
    return
  end if
  
  k_search = k_neighbors
  k_eff_neighbors = max(k_search - 1, 1)

  ! =====================
  ! Build KD-tree
  ! =====================
  dimension_order = [(i, i=1,n_coord_dims)]

  ! Debugging: ! print KD-tree construction status
  ! print *, "Debug: Building KD-tree with n_coord_dims=", n_coord_dims, ", n_points=", n_points
  call build_kd_index( &
    coords, n_coord_dims, n_points, &
    kd_indices, dimension_order, &
    workspace, value_buffer, permutation, &
    left_stack, right_stack, recursion_stack, ierr )
  if (.not. is_ok(ierr)) then
    ! print *, "Error: KD-tree construction failed. ierr=", ierr
    return
  end if

  ! =====================
  ! Main loop
  ! =====================
  do i = 1, n_points
    ! Debugging: Processing point
    !print *, "Debug: Processing point", i
    !print *, "Debug: Coordinates of point:", coords(:,i)
    !print *, "Debug: Original vector:", vectors(:,i)

    ! 1. Initialize and query K neighbors (including self)
    distances(1:k_neighbors) = huge(1.0_real64)
    call kd_knn_query( &
      coords, kd_indices, n_coord_dims, n_points, dimension_order, &
      coords(:,i), k_neighbors, neighbors, distances, value_buffer, ierr )
    if (.not. is_ok(ierr)) return

    if (i == 250) then 
      ! Debugging: After kd_knn_query
      !print *, "Debug: kd_knn_query succeeded for point", i
      !print *, "Debug: Neighbors indices:", neighbors
      !print *, "Debug: Distances:", distances
    end if 

    ! 2. Sort distances and obtain permutation_distances
    permutation_distances = [(i_perm, i_perm = 1, k_neighbors)]
    call sort_real(distances(1:k_neighbors), permutation_distances, left_stack, right_stack)

    if (i == 250) then
      ! Debugging: After sorting distances
      !print *, "Debug: Distances sorted for point", i, distances(permutation_distances)
      !print *, "Debug: permutation_distances array:", permutation_distances
    end if

    ! 3. Reorder neighbors using the same permutation_distances
    value_buffer(1:k_neighbors) = real(neighbors(1:k_neighbors), real64)
    do j = 1, k_neighbors
        neighbors(j) = int(value_buffer(permutation_distances(j)), int32)
    end do

    if (i == 250) then
      ! Debugging: After reordering neighbors
      !print *, "Debug: Reordered neighbors:", neighbors
    end if

    ! 4. Calculate the Median (Excluding self and HUGE)
    k_real_neighbors = 0
    k_start = 2
    do j = k_start, k_neighbors
        if (distances(permutation_distances(j)) < huge(1.0_real64)) then
            k_real_neighbors = k_real_neighbors + 1
        else
            exit 
        end if
    end do

    if (k_real_neighbors < 1) then
        local_sigma = 1.0e-12_real64
        smoothed(:,i) = vectors(:,i)
        ! print *, "Debug: No valid neighbors found for point", i
        cycle
    end if

    k_eff_neighbors = k_real_neighbors

    if (k_real_neighbors > 0) then
        if (mod(k_real_neighbors, 2) == 0) then
            local_median = 0.5_real64 * ( &
                distances(permutation_distances(real_neighbor_offset + k_real_neighbors / 2)) + &
                distances(permutation_distances(real_neighbor_offset + k_real_neighbors / 2 + 1)) )
        else
            local_median = distances(permutation_distances(real_neighbor_offset + (k_real_neighbors + 1) / 2))
        end if
    else
        local_median = 0.0_real64
    end if

    if (i == 250) then
      !print *, "Debug: real neighbors", k_eff_neighbors
    end if 

    c = 1.0_real64 / sqrt(log(real(k_eff_neighbors, real64)))
    local_sigma = max(c * local_median, 1.0e-12_real64)
    sd_arr(i) = local_sigma

    ! Debugging: Local sigma and median
    if (i == 250) then
      !print *, "Debug: Local sigma:", local_sigma, "Local median:", local_median
    end if
    if (anisotropy_mode == 2) then
    
      ! Dispersion calculation with respect to the central point (x_i)
      cov(:,:) = 0.0_real64
      k_eff = 0

      ! Dispersion loop (starting from the second neighbor)
      do j = 2, k_neighbors
        idx = neighbors(j)
        ! Calculation with respect to the central point (x_i)
        delta(:) = coords(:,idx) - coords(:,i)
        ! cov += delta * delta^T  
        cov = cov + spread(delta, dim=2, ncopies=n_coord_dims) * spread(delta, dim=1, ncopies=n_coord_dims)
        k_eff = k_eff + 1
      end do
      
      ! Normalize to the average dispersion (if k_eff > 0)
      if (k_eff > 0) then
        cov = cov / real(k_eff, real64)
      else
        cov(:,:) = 0.0_real64
      end if

      ! Trace calculation and regularization
      trace_cov = 0.0_real64
      do j = 1, n_coord_dims
        trace_cov = trace_cov + cov(j,j)
      end do

      call identity_matrix(Iden, n_coord_dims)
      
      ! Add a regularization factor to the matrix
      eps_reg = 1.0e-6_real64 * max(trace_cov / real(max(1_int32,n_coord_dims), real64), 1.0e-12_real64)
      cov = cov + eps_reg * Iden

      ! ******************************************************
      ! SCALING BY THE LOCAL BANDWIDTH (local_sigma)
      ! ******************************************************
      ! Recalculate trace after regularization
      trace_cov = 0.0_real64
      do j = 1, n_coord_dims
        trace_cov = trace_cov + cov(j,j)
      end do
      
      if (trace_cov > 1.0e-12_real64) then
          ! Scale so that the "volume" of the kernel matches the local_sigma
          scale_factor = (local_sigma**2 * real(n_coord_dims, real64)) / trace_cov
          cov = cov * scale_factor
      end if
      ! ******************************************************
      
      ! Inversion
      cov_inv(:,:) = cov
      call invert_small_matrix_spd(cov_inv, n_coord_dims, ierr)
      if (.not. is_ok(ierr)) return
    end if

    ! 5. Weight loop
    local_wsum = 0.0_real64
    local_work = 0.0_real64

    do j = 2, k_neighbors
      idx = neighbors(j)
      if (distances(permutation_distances(j)) == huge(1.0_real64)) cycle

      select case (anisotropy_mode)
      case (0)
        d2 = distances(permutation_distances(j))**2
        local_w = exp(-d2 / (2.0_real64 * local_sigma**2))
      case (1)
        delta = coords(:,idx) - coords(:,i)
        sigma_parallel = local_sigma
        sigma_orth  = anisotropy_factor * local_sigma
        d2 = delta(1)**2 / (2.0_real64 * sigma_parallel**2)
        if (n_coord_dims > 1) then
          d2 = d2 + sum(delta(2:n_coord_dims)**2) / (2.0_real64 * sigma_orth**2)
        end if
        local_w = exp(-d2)
      case (2)
        delta = coords(:,idx) - coords(:,i)
        tmp_vec = matmul(cov_inv, delta)
        d2 = dot_product(delta, tmp_vec)
        local_w = exp(-0.5_real64 * d2)
      end select

      local_w = max(local_w, 1.0e-12_real64)
      local_work = local_work + local_w * vectors(:,idx)
      local_wsum = local_wsum + local_w

      ! Debugging: Neighbor analysis
      if (i == 250) then
        if (j < 50 .or. j > 450) then 
          !print *, "Debug: Point", i, "Neighbor", j
          !print *, "Debug: Neighbor index:", idx
          !print *, "Debug: Distance squared (d2):", d2
          !print *, "Debug: Weight (local_w):", local_w
          !print *, "Debug: Accumulated local_work:", local_work
          !print *, "Debug: Accumulated local_wsum:", local_wsum
        end if 
      end if
    end do

    if (local_wsum > 0.0_real64) then
      smoothed(:,i) = local_work / local_wsum
    else
      smoothed(:,i) = vectors(:,i)
    end if

    ! Debugging: Final smoothed vector
    ! print *, "Debug: Final smoothed vector for point", i, smoothed(:,i)
  end do

  ! print *, "Debug: Exiting main loop"

  end subroutine smooth_vectors_gaussian_adaptive

pure subroutine identity_matrix(A, n)
    integer(int32), intent(in) :: n
    real(real64),   intent(out) :: A(n,n)
    integer(int32) :: ii
    A(:,:) = 0.0_real64
    do ii = 1, n
      A(ii,ii) = 1.0_real64
    end do
  end subroutine identity_matrix

  subroutine invert_small_matrix_spd(A, n, ierr_local)
    ! In-place inverse of SPD matrix using Cholesky (dpotrf+dpotri).
    integer(int32), intent(in) :: n
    real(real64),   intent(inout) :: A(n,n)
    integer(int32), intent(out) :: ierr_local
    integer(int32) :: info, r, c

    call set_ok(ierr_local)

    call dpotrf('L', n, A, n, info)
    if (info /= 0) then
      call set_err_once(ierr_local, ERR_INVALID_INPUT)
      return
    end if

    call dpotri('L', n, A, n, info)
    if (info /= 0) then
      call set_err_once(ierr_local, ERR_INVALID_INPUT)
      return
    end if

    ! Symmetrize: copy lower -> upper
    do r = 1, n
      do c = r+1, n
        A(r,c) = A(c,r)
      end do
    end do
  end subroutine invert_small_matrix_spd



end module anwil

