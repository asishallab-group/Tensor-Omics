module knn_smoothing
  use kd_tree, only: build_kd_index, kd_knn_query
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
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
  end subroutine print_timing_stats

  ! Helper function to get current time
  function get_time() result(time_seconds)
    integer(int32) :: count, count_rate
    real(real64) :: time_seconds
    call system_clock(count, count_rate)
    time_seconds = real(count, real64) / real(count_rate, real64)
  end function get_time


  subroutine smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
                                             n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                                             kd_indices, dimension_order, neighbors, distances, &
                                             workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    ! Smooths vectors using adaptive Gaussian KNN based on coordinate space
    ! Searches for neighbors in coordinate space, applies smoothing to vector space
    ! The Gaussian kernel size adapts to local density: sigma = 3 * std_dev(k_neighbor_distances)
    integer(int32), intent(in) :: n_coord_dims   ! number of coordinate dimensions for neighbor search
    integer(int32), intent(in) :: n_vector_dims  ! number of vector dimensions to smooth
    integer(int32), intent(in) :: n_points       ! total number of points
    integer(int32), intent(in) :: k_neighbors    ! number of nearest neighbors to use
    real(real64), intent(in) :: coords(n_coord_dims, n_points)    ! coordinates for neighbor search
    real(real64), intent(in) :: vectors(n_vector_dims, n_points)  ! vectors to be smoothed
    real(real64), intent(out) :: smoothed(n_vector_dims, n_points) ! output: smoothed vectors
    
    ! Buffers for k-d tree construction
    integer(int32), intent(out) :: kd_indices(n_points), dimension_order(n_coord_dims)
    integer(int32), intent(inout) :: workspace(n_points), permutation(n_points)
    integer(int32), intent(inout) :: left_stack(n_points), right_stack(n_points)
    real(real64), intent(inout) :: value_buffer(n_points)
    
    ! Buffers for smoothing
    integer(int32), intent(inout) :: neighbors(k_neighbors)
    real(real64), intent(inout) :: distances(k_neighbors)
    integer(int32), intent(out) :: ierr
    
    ! Internal variables
    real(real64) :: kd_workspace(n_coord_dims)
    integer(int32) :: recursion_stack(3, n_points)
    integer(int32) :: i, j, query_ierr
    real(real64) :: sigma, wsum, w, work(n_vector_dims)
    real(real64) :: mean_dist, variance, std_dev
    real(real64) :: kd_start_time, kd_end_time, queries_start_time, queries_end_time
    ! Local variables for each thread (automatically private in do concurrent)
    integer(int32) :: local_neighbors(k_neighbors), local_query_ierr, local_j
    real(real64) :: local_distances(k_neighbors), local_kd_workspace(n_coord_dims)
    real(real64) :: local_mean_dist, local_variance, local_std_dev, local_sigma
    real(real64) :: local_wsum, local_w, local_work(n_vector_dims)
    
    call set_ok(ierr)
    
    ! Validation
    if (n_coord_dims < 1 .or. n_vector_dims < 1 .or. n_points < 1 .or. &
        k_neighbors < 1 .or. k_neighbors > n_points) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if
    
    ! Simple dimension order (vectorized)
    dimension_order = [(i, i = 1, n_coord_dims)]
    
    ! Step 1: Build k-d tree using coordinates only (with timing)
    kd_start_time = get_time()
    call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                       workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
    kd_end_time = get_time()
    timing_kd_build = timing_kd_build + (kd_end_time - kd_start_time)
    if (.not. is_ok(ierr)) return
    
    ! Step 2: Query k-d tree and apply Gaussian smoothing in one parallel loop!
    queries_start_time = get_time()
    timing_query_count = n_points  ! We'll process all points
    
    do i = 1, n_points
       
       ! Query k-d tree using coordinate i as query point
       call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                        coords(:, i), k_neighbors, local_neighbors, local_distances, local_kd_workspace, local_query_ierr)
       
       ! Calculate adaptive sigma based on 3 times the standard deviation of distances
       ! The Gaussian kernel adjusts to local dispersion of the k neighbors
       local_mean_dist = sum(local_distances(1:k_neighbors)) / real(k_neighbors, real64)
       local_variance = sum((local_distances(1:k_neighbors) - local_mean_dist)**2) / real(k_neighbors, real64)
       local_std_dev = sqrt(local_variance)
       local_sigma = 3.0_real64 * local_std_dev
       
       ! Fallback to avoid sigma too small
       if (local_sigma <= 1.0e-12_real64) local_sigma = 1.0e-12_real64

       ! Apply Gaussian weighting to smooth the vectors (not coordinates!)
       ! IMPORTANT: Skip neighbors with exactly 0.0 values (placeholders for interpolation)
       local_work(:) = 0.0_real64
       local_wsum    = 0.0_real64

       do local_j = 1, k_neighbors
          ! Skip neighbors that have zero values (likely interpolation placeholders)
          if (all(abs(vectors(:, local_neighbors(local_j))) > 1.0e-15_real64)) then
             local_w = exp(-(local_distances(local_j)**2) / (2.0_real64 * local_sigma**2))
             local_work(:) = local_work(:) + local_w * vectors(:, local_neighbors(local_j))
             local_wsum    = local_wsum + local_w
          end if
       end do

       if (local_wsum > 0.0_real64) then
          ! Normal case: we found valid neighbors with non-zero values
          smoothed(:, i) = local_work(:) / local_wsum
       else
          ! Fallback case: no valid neighbors found (all were zero placeholders)
          ! Keep original value if non-zero, otherwise leave as zero for interpolation detection
          if (all(abs(vectors(:, i)) > 1.0e-15_real64)) then
             smoothed(:, i) = vectors(:, i)   ! Keep known value unchanged
          else
             smoothed(:, i) = 0.0_real64     ! Mark as still needing interpolation
          end if
       end if

    end do
    
    queries_end_time = get_time()
    ! Note: Total time includes both queries and Gaussian calc since they're combined
    timing_knn_queries = timing_knn_queries + (queries_end_time - queries_start_time)
    timing_gaussian_calc = 0.0_real64  ! Set to 0 since it's included in queries now

  end subroutine smooth_vectors_gaussian_adaptive

end module knn_smoothing


SUBROUTINE smooth_vectors_c(n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                            coords_ptr, vectors_ptr_in, vectors_ptr_out, &
                            kd_indices_ptr, dimension_order_ptr, neighbors_ptr, distances_ptr, &
                            workspace_int_ptr, value_buffer_ptr, permutation_ptr, &
                            left_stack_ptr, right_stack_ptr, ierr) BIND(C, name='smooth_vectors_c')
  USE, INTRINSIC :: ISO_C_BINDING
  USE knn_smoothing, only: smooth_vectors_gaussian_adaptive, set_ok, set_err_once, ERR_INVALID_INPUT, real64, int32
  IMPLICIT NONE
  
  ! === ARGUMENTOS ESCALARES ===
  INTEGER(C_INT), INTENT(IN), VALUE :: n_coord_dims
  INTEGER(C_INT), INTENT(IN), VALUE :: n_vector_dims
  INTEGER(C_INT), INTENT(IN), VALUE :: n_points
  INTEGER(C_INT), INTENT(IN), VALUE :: k_neighbors
  
  ! === PUNTEROS PRINCIPALES (REAL(real64)) ===
  TYPE(C_PTR), INTENT(IN), VALUE :: coords_ptr         ! coords (n_coord_dims x n_points)
  TYPE(C_PTR), INTENT(IN), VALUE :: vectors_ptr_in     ! vectors (n_vector_dims x n_points)
  TYPE(C_PTR), INTENT(IN), VALUE :: vectors_ptr_out    ! smoothed (n_vector_dims x n_points)

  ! === PUNTEROS DE BUFFERS ENTEROS (INTEGER(int32)) ===
  TYPE(C_PTR), INTENT(IN), VALUE :: kd_indices_ptr      ! (n_points)
  TYPE(C_PTR), INTENT(IN), VALUE :: dimension_order_ptr ! (n_coord_dims)
  TYPE(C_PTR), INTENT(IN), VALUE :: neighbors_ptr       ! (k_neighbors)
  TYPE(C_PTR), INTENT(IN), VALUE :: workspace_int_ptr   ! workspace (n_points)
  TYPE(C_PTR), INTENT(IN), VALUE :: permutation_ptr     ! permutation (n_points)
  TYPE(C_PTR), INTENT(IN), VALUE :: left_stack_ptr      ! left_stack (n_points)
  TYPE(C_PTR), INTENT(IN), VALUE :: right_stack_ptr     ! right_stack (n_points)
  
  ! === PUNTEROS DE BUFFERS REALES (REAL(real64)) ===
  TYPE(C_PTR), INTENT(IN), VALUE :: distances_ptr       ! (k_neighbors)
  TYPE(C_PTR), INTENT(IN), VALUE :: value_buffer_ptr    ! (n_points)
  
  INTEGER(C_INT), INTENT(OUT) :: ierr
  
  ! === PUNTEROS FORTRAN INTERNOS (Descriptores creados) ===
  REAL(real64), POINTER :: coords_matrix(:, :), vectors_matrix(:, :), smoothed_matrix(:, :)
  REAL(real64), POINTER :: distances(:), value_buffer(:)

  INTEGER(int32), POINTER :: kd_indices(:), dimension_order(:), neighbors(:)
  INTEGER(int32), POINTER :: workspace(:), permutation(:), left_stack(:), right_stack(:)
  
  ! === Workspace faltante en los argumentos de C pero necesario para la llamada Fortran ===
  ! Se asigna aquí, ya que no se pasó un puntero C para él.
  REAL(real64) :: kd_workspace(n_coord_dims) 
  
  INTEGER :: recursion_stack(3, n_points) ! Stack de recursión necesario para build_kd_index
  
  INTEGER :: validation_error = 0
  
  CALL set_ok(ierr)

  ! === VALIDACIÓN DE PUNTEROS NULL ===
  IF (.NOT. C_ASSOCIATED(coords_ptr) .OR. .NOT. C_ASSOCIATED(vectors_ptr_in) .OR. &
      .NOT. C_ASSOCIATED(vectors_ptr_out)) validation_error = 1
  
  IF (validation_error == 1) THEN
      CALL set_err_once(ierr, ERR_INVALID_INPUT)
      RETURN
  END IF

  ! === CONVERSIÓN DE PUNTEROS C A DESCRIPTORES FORTRAN ===
  ! Matrices 2D
  CALL C_F_POINTER(coords_ptr, coords_matrix, SHAPE=[n_coord_dims, n_points])
  CALL C_F_POINTER(vectors_ptr_in, vectors_matrix, SHAPE=[n_vector_dims, n_points])
  CALL C_F_POINTER(vectors_ptr_out, smoothed_matrix, SHAPE=[n_vector_dims, n_points])
  
  ! Arrays 1D (REAL(real64))
  CALL C_F_POINTER(distances_ptr, distances, SHAPE=[k_neighbors])
  CALL C_F_POINTER(value_buffer_ptr, value_buffer, SHAPE=[n_points])

  ! Arrays 1D (INTEGER(int32))
  CALL C_F_POINTER(kd_indices_ptr, kd_indices, SHAPE=[n_points])
  CALL C_F_POINTER(dimension_order_ptr, dimension_order, SHAPE=[n_coord_dims])
  CALL C_F_POINTER(neighbors_ptr, neighbors, SHAPE=[k_neighbors])
  CALL C_F_POINTER(workspace_int_ptr, workspace, SHAPE=[n_points])
  CALL C_F_POINTER(permutation_ptr, permutation, SHAPE=[n_points])
  CALL C_F_POINTER(left_stack_ptr, left_stack, SHAPE=[n_points])
  CALL C_F_POINTER(right_stack_ptr, right_stack, SHAPE=[n_points])
  
  ! === LLAMADA A LA FUNCIÓN PRINCIPAL ===
  ! Se usan los descriptors Fortran creados a partir de los punteros C
  CALL smooth_vectors_gaussian_adaptive(coords_matrix, vectors_matrix, smoothed_matrix, &
                                        n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                                        kd_indices, dimension_order, neighbors, distances, &
                                        workspace, value_buffer, permutation, left_stack, right_stack, ierr)


END SUBROUTINE smooth_vectors_c

!> R interface for KNN smoothing
subroutine smooth_vectors_r(coords, vectors, smoothed, n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                          kd_indices, dimension_order, neighbors, distances, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    use knn_smoothing, only: smooth_vectors_gaussian_adaptive
    use iso_fortran_env, only: int32, real64
    implicit none
    
    integer(int32), intent(in) :: n_coord_dims      !! Number of coordinate dimensions
    integer(int32), intent(in) :: n_vector_dims     !! Number of vector dimensions  
    integer(int32), intent(in) :: n_points          !! Number of points
    integer(int32), intent(in) :: k_neighbors       !! Number of neighbors
    real(real64), intent(in) :: coords(n_coord_dims, n_points)     !! Coordinate matrix
    real(real64), intent(in) :: vectors(n_vector_dims, n_points)   !! Input vectors
    real(real64), intent(out) :: smoothed(n_vector_dims, n_points) !! Output smoothed vectors
    
    ! Workspace arrays
    integer(int32), intent(out) :: kd_indices(n_points)            !! K-d tree indices
    integer(int32), intent(out) :: dimension_order(n_coord_dims)   !! Dimension order
    integer(int32), intent(inout) :: neighbors(k_neighbors)        !! Neighbor buffer
    real(real64), intent(inout) :: distances(k_neighbors)          !! Distance buffer
    integer(int32), intent(inout) :: workspace(n_points)           !! Integer workspace
    real(real64), intent(inout) :: value_buffer(n_points)          !! Real workspace
    integer(int32), intent(inout) :: permutation(n_points)         !! Permutation array
    integer(int32), intent(inout) :: left_stack(n_points)          !! Left stack
    integer(int32), intent(inout) :: right_stack(n_points)         !! Right stack
    integer(int32), intent(out) :: ierr                            !! Error code

    call smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
                                        n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                                        kd_indices, dimension_order, neighbors, distances, &
                                        workspace, value_buffer, permutation, left_stack, right_stack, ierr)
end subroutine smooth_vectors_r

!> R interface for timing statistics
subroutine get_timing_stats_r(kd_build_time, knn_queries_time, gaussian_calc_time, query_count)
    use knn_smoothing, only: timing_kd_build, timing_knn_queries, timing_gaussian_calc, timing_query_count
    use iso_fortran_env, only: int32, real64
    implicit none
    
    real(real64), intent(out) :: kd_build_time      !! K-d tree build time
    real(real64), intent(out) :: knn_queries_time   !! KNN queries time
    real(real64), intent(out) :: gaussian_calc_time !! Gaussian calculation time
    integer(int32), intent(out) :: query_count      !! Total number of queries
    
    kd_build_time = timing_kd_build
    knn_queries_time = timing_knn_queries
    gaussian_calc_time = timing_gaussian_calc
    query_count = timing_query_count
end subroutine get_timing_stats_r

!> R interface for resetting timing statistics
subroutine reset_timing_stats_r()
    use knn_smoothing, only: reset_timing_stats
    implicit none
    
    call reset_timing_stats()
end subroutine reset_timing_stats_r

!> R interface for printing timing statistics
subroutine print_timing_stats_r(n_points, k_neighbors)
    use knn_smoothing, only: print_timing_stats
    use iso_fortran_env, only: int32
    implicit none
    
    integer(int32), intent(in) :: n_points     !! Number of points
    integer(int32), intent(in) :: k_neighbors  !! Number of neighbors
    
    call print_timing_stats(n_points, k_neighbors)
end subroutine print_timing_stats_r