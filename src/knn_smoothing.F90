module knn_smoothing
  use kd_tree, only: build_kd_index, kd_knn_query
  use iso_fortran_env, only: int32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
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

    !---------------------------------------------------------------
  ! Suavizado 1D con penalización de diferencias:
  ! minimiza  sum (y - y_in)^2 + lambda * sum (y_{i+1} - y_i)^2
  !---------------------------------------------------------------
  subroutine smooth_1d_diff_penalty(n, y_in, lambda_pen, y_out)
    implicit none
    integer(int32), intent(in)  :: n
    real(real64),   intent(in)  :: y_in(n)
    real(real64),   intent(in)  :: lambda_pen
    real(real64),   intent(out) :: y_out(n)

    real(real64) :: a(n), b(n), c(n), d_rhs(n)
    real(real64) :: c_star(n), d_star(n), denom
    integer(int32) :: i

    if (n <= 1) then
      if (n == 1) y_out(1) = y_in(1)
      return
    end if

    ! Sistema tridiagonal:
    ! b = diagonal, a = subdiagonal, c = superdiagonal, d_rhs = lado derecho
    b(1)    = 1.0_real64 + lambda_pen
    c(1)    = -lambda_pen
    d_rhs(1)= y_in(1)

    do i = 2, n-1
      a(i)    = -lambda_pen
      b(i)    = 1.0_real64 + 2.0_real64*lambda_pen
      c(i)    = -lambda_pen
      d_rhs(i)= y_in(i)
    end do

    a(n)     = -lambda_pen
    b(n)     = 1.0_real64 + lambda_pen
    c(n)     = 0.0_real64
    d_rhs(n) = y_in(n)

    ! Thomas algorithm (forward sweep)
    c_star(1) = c(1) / b(1)
    d_star(1) = d_rhs(1) / b(1)

    do i = 2, n
      denom    = b(i) - a(i) * c_star(i-1)
      if (abs(denom) < 1.0e-18_real64) denom = 1.0e-18_real64
      c_star(i) = c(i) / denom
      d_star(i) = (d_rhs(i) - a(i) * d_star(i-1)) / denom
    end do

    ! Back substitution
    y_out(n) = d_star(n)
    do i = n-1, 1, -1
      y_out(i) = d_star(i) - c_star(i) * y_out(i+1)
    end do
  end subroutine smooth_1d_diff_penalty

    subroutine smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
                                             n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                                             kd_indices, dimension_order, neighbors, distances, &
                                             workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    ! Smooths vectors using adaptive Gaussian KNN based on coordinate space
    ! En 1D: búsqueda directa + log-space + penalización de diferencias
    ! En >1D: k-d tree + Gaussian KNN (como antes)

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

    ! Local variables (también usados en rama 1D)
    integer(int32) :: local_neighbors(k_neighbors), local_query_ierr, local_j
    real(real64) :: local_distances(k_neighbors), local_kd_workspace(n_coord_dims)
    real(real64) :: local_mean_dist, local_variance, local_std_dev, local_sigma
    real(real64) :: local_wsum, local_w, local_work(n_vector_dims)
    integer(int32) :: ni, left, right, count, d

    integer(int32) :: idx_sorted(n_points), idx_inv(n_points)
    real(real64) :: coords_sorted(n_points)
    real(real64) :: coord_log_sorted(n_points)
    real(real64) :: global_range, mean_spacing, sigma_floor
    real(real64) :: y_in(n_points), y_out(n_points)
    real(real64), parameter :: lambda_pen = 5.0_real64  ! control de suavizado de curvatura
    ! Mezcla local-global
     real(real64), parameter :: alpha_mix = 0.07
     real(real64) :: global_mean, global_variance
    real(real64) :: sigma_global, sigma_global2
    real(real64) :: sigma_final2, local_sigma2
    real(real64) :: inv_two_sigma2
    real(real64) :: log_vectors(n_vector_dims, n_points)


    call set_ok(ierr)
    
    ! Validation
    if (n_coord_dims < 1 .or. n_vector_dims < 1 .or. n_points < 1 .or. &
        k_neighbors < 1 .or. k_neighbors > n_points) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if
    
    ! ============================================================
    ! CASO 1D: coords es 1D -> búsqueda directa + log-space
    ! ============================================================
    if (n_coord_dims == 1) then

      ! Inicializar permutación
      idx_sorted    = [(i, i=1, n_points)]
      coords_sorted = coords(1,:)

      ! Ordenar los puntos por la coordenada original
      call sort_real(coords_sorted, idx_sorted, left_stack, right_stack)

      ! Precomputar log(coord+1) en el orden ordenado
      do i = 1, n_points
        coord_log_sorted(i)      = log(coords(1, idx_sorted(i)) + 1.0_real64)
        idx_inv(idx_sorted(i))   = i
      end do

      ! --------------------------------------------
      ! sigma_global = std global de posiciones log
      ! --------------------------------------------
      global_mean = sum(coord_log_sorted) / real(n_points, real64)
      global_variance = sum((coord_log_sorted - global_mean)**2) / real(n_points-1, real64)
      sigma_global = sqrt(global_variance)
      sigma_global2 = sigma_global*sigma_global

      

      ! Precomputar log(vectores)
      do i = 1, n_points
        do d = 1, n_vector_dims
          if (.not. ieee_is_nan(vectors(d, i))) then
            log_vectors(d, i) = log(vectors(d, i) + 1.0_real64)
          else
            log_vectors(d, i) = 0.0_real64
          end if
        end do
      end do

      queries_start_time = get_time()
      timing_query_count = n_points

      do i = 1, n_points
        ni = idx_inv(i)

        ! Ventana local de k vecinos ordenados
        left  = max(1, ni - k_neighbors/2)
        right = min(n_points, left + k_neighbors - 1)
        if (right - left + 1 < k_neighbors) left = max(1, right - k_neighbors + 1)

        count = 0
        do j = left, right
          count = count + 1
          local_neighbors(count) = idx_sorted(j)
          local_distances(count) = abs(coord_log_sorted(ni) - coord_log_sorted(j))
        end do

        ! ======================================================
        ! sigma_local = 3 * std(local_distances(2..k))
        ! (requerimiento explícito de Asis)
        ! ======================================================
        if (k_neighbors > 1) then
          local_mean_dist = sum(local_distances(2:k_neighbors)) / real(k_neighbors-1, real64)
          local_variance  = sum((local_distances(2:k_neighbors) - local_mean_dist)**2) / &
                            real(k_neighbors-1, real64)
          local_std_dev = sqrt(local_variance)
        else
          local_std_dev = 0.0_real64
        end if

        local_sigma = 3.0_real64 * local_std_dev  ! Asis requirement
        local_sigma2 = local_sigma*local_sigma

        ! ======================================================
        ! mezcla local/global:
        ! sigma_final² = (1-a)sigma_local² + a*sigma_global²
        ! ======================================================
        sigma_final2 = (1.0_real64 - alpha_mix) * local_sigma2 + alpha_mix * sigma_global2
        local_sigma = sqrt(sigma_final2)

        ! ======================================================
        ! smoothing EXCLUYENDO SELF
        ! ======================================================
        inv_two_sigma2 = 0.5_real64 / (local_sigma * local_sigma)

        local_work(:) = 0.0_real64
        local_wsum    = 0.0_real64

        do local_j = 2, k_neighbors   ! self excluido
          do d = 1, n_vector_dims
            if (.not. ieee_is_nan(log_vectors(d, local_neighbors(local_j)))) then
              local_w = exp(-(local_distances(local_j)**2) * inv_two_sigma2)
              local_work(d) = local_work(d) + local_w * log_vectors(d, local_neighbors(local_j))
              if (d == 1) local_wsum = local_wsum + local_w
            end if
          end do
        end do

        if (local_wsum > 0.0_real64) then
          smoothed(:, i) = local_work(:) / local_wsum   ! log-space
        else
          smoothed(:, i) = log_vectors(:, i)
        end if

      end do

      ! Volver a escala original
      do i = 1, n_points
        do d = 1, n_vector_dims
          smoothed(d, i) = exp(smoothed(d, i)) - 1.0_real64
        end do
      end do

      queries_end_time    = get_time()
      timing_kd_build     = 0.0_real64
      timing_knn_queries  = timing_knn_queries + (queries_end_time - queries_start_time)
      timing_gaussian_calc= 0.0_real64
      return
    end if


    ! ============================================================
    ! CASO N-DIM (>1): ruta original con k-d tree (sin cambios grandes)
    ! ============================================================
    dimension_order = [(i, i = 1, n_coord_dims)]
    
    kd_start_time = get_time()
    call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                       workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
    kd_end_time = get_time()
    timing_kd_build = timing_kd_build + (kd_end_time - kd_start_time)
    if (.not. is_ok(ierr)) return
    
    queries_start_time = get_time()
    timing_query_count = n_points
    
    do i = 1, n_points
       call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                        coords(:, i), k_neighbors, local_neighbors, local_distances, local_kd_workspace, local_query_ierr)
       
       local_mean_dist = sum(local_distances(1:k_neighbors)) / real(k_neighbors, real64)
       local_variance  = sum((local_distances(1:k_neighbors) - local_mean_dist)**2) / real(k_neighbors, real64)
       local_std_dev   = sqrt(local_variance)
       local_sigma     = 3.0_real64 * local_std_dev
       
       if (local_sigma <= 1.0e-12_real64) local_sigma = 1.0e-12_real64

       local_work(:) = 0.0_real64
       local_wsum    = 0.0_real64

       do local_j = 1, k_neighbors
          if (all(abs(vectors(:, local_neighbors(local_j))) > 1.0e-15_real64)) then
             local_w = exp(-(local_distances(local_j)**2) / (2.0_real64 * local_sigma**2))
             local_work(:) = local_work(:) + local_w * vectors(:, local_neighbors(local_j))
             local_wsum    = local_wsum + local_w
          end if
       end do

       if (local_wsum > 0.0_real64) then
          smoothed(:, i) = local_work(:) / local_wsum
       else
          if (all(abs(vectors(:, i)) > 1.0e-15_real64)) then
             smoothed(:, i) = vectors(:, i)
          else
             smoothed(:, i) = 0.0_real64
          end if
       end if
    end do
    
    queries_end_time = get_time()
    timing_knn_queries = timing_knn_queries + (queries_end_time - queries_start_time)
    timing_gaussian_calc = 0.0_real64

  end subroutine smooth_vectors_gaussian_adaptive


end module knn_smoothing