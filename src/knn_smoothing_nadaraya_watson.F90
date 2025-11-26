module knn_smoothing_nadaraya_watson
  use kd_tree, only: build_kd_index, kd_knn_query
  use iso_fortran_env, only: int32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use tox_errors, only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
  use f42_utils, only: sort_real
  implicit none

  ! === TIMING GLOBALS FOR PERFORMANCE ANALYSIS ===
  real(real64), save :: timing_kd_build      = 0.0_real64
  real(real64), save :: timing_knn_queries   = 0.0_real64
  real(real64), save :: timing_gaussian_calc = 0.0_real64
  integer(int32), save :: timing_query_count = 0

contains

  ! ---------------------------------------------------------
  ! TIMING UTILS
  ! ---------------------------------------------------------
  subroutine reset_timing_stats()
    timing_kd_build      = 0.0_real64
    timing_knn_queries   = 0.0_real64
    timing_gaussian_calc = 0.0_real64
    timing_query_count   = 0
  end subroutine reset_timing_stats

  subroutine print_timing_stats(n_points, k_neighbors)
    integer(int32), intent(in) :: n_points, k_neighbors
    real(real64) :: total_time, queries_per_sec, avg_query_time

    total_time = timing_kd_build + timing_knn_queries + timing_gaussian_calc
    if (timing_query_count > 0) then
      avg_query_time = timing_knn_queries / real(timing_query_count, real64)
      if (timing_knn_queries > 0.0_real64) then
        queries_per_sec = real(timing_query_count, real64) / timing_knn_queries
      else
        queries_per_sec = 0.0_real64
      end if
    else
      avg_query_time = 0.0_real64
      queries_per_sec = 0.0_real64
    end if

    write(*, '(A)') "=== FORTRAN INTERNAL TIMING BREAKDOWN ==="
    write(*, '(A, I0, A, I0)') "Dataset: ", n_points, " points, k=", k_neighbors
    if (total_time > 0.0_real64) then
      write(*,'(A,F10.6,A,F5.1,A)') "1. K-d tree build:   ", timing_kd_build, " sec (", &
                                     100.0*timing_kd_build/total_time, "%)"
      write(*,'(A,F10.6,A,F5.1,A)') "2. KNN+Gaussian:     ", timing_knn_queries, " sec (", &
                                     100.0*timing_knn_queries/total_time, "%)"
      write(*,'(A,F10.6,A,F5.1,A)') "3. (unused):         ", timing_gaussian_calc, " sec (", &
                                     100.0*timing_gaussian_calc/total_time, "%)"
    else
      write(*,'(A,F10.6,A)') "1. K-d tree build:   ", timing_kd_build, " sec"
      write(*,'(A,F10.6,A)') "2. KNN+Gaussian:     ", timing_knn_queries, " sec"
      write(*,'(A,F10.6,A)') "3. (unused):         ", timing_gaussian_calc, " sec"
    end if
    write(*,'(A,F10.6,A)') "Total Fortran time:  ", total_time, " sec"
    write(*,'(A,I0)')       "Total queries made:  ", timing_query_count
    write(*,'(A,F10.6,A)')  "Avg query time:      ", avg_query_time*1000.0, " ms"
    write(*,'(A,F8.1,A)')   "Queries per second:  ", queries_per_sec, " queries/sec"
    write(*,'(A)') "========================================="
  end subroutine print_timing_stats

  function get_time() result(time_seconds)
    integer(int32) :: count, count_rate
    real(real64)   :: time_seconds
    call system_clock(count, count_rate)
    time_seconds = real(count, real64) / real(count_rate, real64)
  end function get_time

  ! (Se queda por si lo quieres usar en otra parte; no lo usaremos en NW full)
  subroutine smooth_1d_diff_penalty(n, y_in, lambda_pen, y_out)
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

    b(1)     = 1.0_real64 + lambda_pen
    c(1)     = -lambda_pen
    d_rhs(1) = y_in(1)

    do i = 2, n-1
      a(i)     = -lambda_pen
      b(i)     = 1.0_real64 + 2.0_real64*lambda_pen
      c(i)     = -lambda_pen
      d_rhs(i) = y_in(i)
    end do

    a(n)     = -lambda_pen
    b(n)     = 1.0_real64 + lambda_pen
    c(n)     = 0.0_real64
    d_rhs(n) = y_in(n)

    c_star(1) = c(1) / b(1)
    d_star(1) = d_rhs(1) / b(1)

    do i = 2, n
      denom = b(i) - a(i)*c_star(i-1)
      if (abs(denom) < 1.0e-18_real64) denom = 1.0e-18_real64
      c_star(i) = c(i) / denom
      d_star(i) = (d_rhs(i) - a(i)*d_star(i-1)) / denom
    end do

    y_out(n) = d_star(n)
    do i = n-1, 1, -1
      y_out(i) = d_star(i) - c_star(i)*y_out(i+1)
    end do
  end subroutine smooth_1d_diff_penalty

  ! ---------------------------------------------------------
  ! ADAPTIVE GAUSSIAN NADARAYA–WATSON SMOOTHING
  ! ---------------------------------------------------------
  subroutine smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
       n_coord_dims, n_vector_dims, n_points, k_neighbors, &
       kd_indices, dimension_order, neighbors, distances, &
       workspace, value_buffer, permutation, left_stack, right_stack, sigma_factor, ierr)

    integer(int32), intent(in)    :: n_coord_dims
    integer(int32), intent(in)    :: n_vector_dims
    integer(int32), intent(in)    :: n_points
    integer(int32), intent(in)    :: k_neighbors
    real(real64),  intent(in)     :: coords(n_coord_dims, n_points)
    real(real64),  intent(in)     :: vectors(n_vector_dims, n_points)
    real(real64),  intent(out)    :: smoothed(n_vector_dims, n_points)

    integer(int32), intent(out)   :: kd_indices(n_points), dimension_order(n_coord_dims)
    integer(int32), intent(inout) :: workspace(n_points), permutation(n_points)
    integer(int32), intent(inout) :: left_stack(n_points), right_stack(n_points)
    real(real64),  intent(inout)  :: value_buffer(n_points)

    integer(int32), intent(inout) :: neighbors(k_neighbors)
    real(real64),  intent(inout)  :: distances(k_neighbors)
    integer(int32), intent(out)   :: ierr

    ! Internos generales
    integer(int32) :: i, j, d
    real(real64)   :: kd_start_time, kd_end_time, queries_start_time, queries_end_time

    ! Para 1D (NW completo)
    real(real64)   :: all_distances(n_points)
    real(real64)   :: mean_dist, variance, std_dev, sigma_factor
    real(real64)   :: min_pos, inv_two_sigma2
    real(real64)   :: wsum, w
    real(real64)   :: work(n_vector_dims)
    logical        :: neighbor_valid, is_zero
    integer(int32) :: count_nonzero
    real(real64)   :: sigma_eff    ! <<--- sigma efectiva local


    ! Para N-D (ruta anterior basada en k-d tree + KNN)
    real(real64)   :: kd_workspace(n_coord_dims)
    integer(int32) :: recursion_stack(3, n_points)
    integer(int32) :: query_ierr
    integer(int32) :: local_neighbors(k_neighbors)
    real(real64)   :: local_distances(k_neighbors)
    real(real64)   :: local_mean_dist, local_variance, local_std_dev, local_sigma

    call set_ok(ierr)

    ! ----------------- validación -----------------
    if (n_coord_dims < 1 .or. n_vector_dims < 1 .or. n_points < 1 .or. &
        k_neighbors < 1 .or. k_neighbors > n_points) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! ======================================================
    ! CASO 1D: Nadaraya–Watson completo (toma TODOS los puntos)
    ! ======================================================
    if (n_coord_dims == 1) then
      queries_start_time = get_time()
      timing_query_count = n_points

      do i = 1, n_points
        ! 1) Distancias a TODOS los puntos
        do j = 1, n_points
          all_distances(j) = abs(coords(1, i) - coords(1, j))
        end do

        ! 2) Sigma local según Asis: 3 * std(dist(≠ self))
        mean_dist     = 0.0_real64
        variance      = 0.0_real64
        count_nonzero = 0

        do j = 1, n_points
          if (j == i) cycle
          mean_dist = mean_dist + all_distances(j)
          count_nonzero = count_nonzero + 1
        end do

        if (count_nonzero > 0) then
          mean_dist = mean_dist / real(count_nonzero, real64)
          do j = 1, n_points
            if (j == i) cycle
            variance = variance + (all_distances(j) - mean_dist)**2
          end do
          variance = variance / real(count_nonzero, real64)
          if (variance < 0.0_real64) variance = 0.0_real64
          std_dev = sqrt(variance)
        else
          std_dev = 0.0_real64
        end if

        if (std_dev > 0.0_real64) then
          sigma_eff = sigma_factor * std_dev    ! p.ej. 3 * std
        else
          sigma_eff = 0.0_real64
        end if

        ! Fallback si sigma es 0 (puntos muy colapsados)
        if (sigma_eff <= 0.0_real64) then
          min_pos = huge(1.0_real64)
          do j = 1, n_points
            if (j == i) cycle
            if (all_distances(j) > 0.0_real64 .and. all_distances(j) < min_pos) then
              min_pos = all_distances(j)
            end if
          end do
          if (min_pos < huge(1.0_real64)) then
            sigma_eff = 3.0_real64 * min_pos
          else
            sigma_eff = 1.0e-12_real64
          end if
        end if

        inv_two_sigma2 = 0.5_real64 / (sigma_eff * sigma_eff)

        ! 3) Nadaraya–Watson: usar TODOS los puntos con pesos gaussianos
        work(:) = 0.0_real64
        wsum    = 0.0_real64

        do j = 1, n_points
          neighbor_valid = .true.
          is_zero        = .true.

          do d = 1, n_vector_dims
            if (ieee_is_nan(vectors(d, j))) then
              neighbor_valid = .false.
              exit
            end if
            if (abs(vectors(d, j)) > 1.0e-15_real64) is_zero = .false.
          end do
          if (.not. neighbor_valid) cycle
          if (is_zero) cycle

          w    = exp( - (all_distances(j)**2) * inv_two_sigma2 )
          wsum = wsum + w
          do d = 1, n_vector_dims
            work(d) = work(d) + w * vectors(d, j)
          end do
        end do

        if (wsum > 0.0_real64) then
          smoothed(:, i) = work(:) / wsum
        else
          ! Sin vecinos válidos: dejar valor original si no es NaN
          do d = 1, n_vector_dims
            if (.not. ieee_is_nan(vectors(d, i))) then
              smoothed(d, i) = vectors(d, i)
            else
              smoothed(d, i) = 0.0_real64
            end if
          end do
        end if
      end do

      queries_end_time    = get_time()
      timing_kd_build     = 0.0_real64
      timing_knn_queries  = timing_knn_queries + (queries_end_time - queries_start_time)
      timing_gaussian_calc= 0.0_real64
      return
    end if  ! fin caso 1D

    ! ======================================================
    ! CASO N-D: usar k-d tree + mismo kernel Gaussiano,
    !           pero limitado a K vecinos (como antes)
    ! ======================================================
    dimension_order = [(i, i = 1, n_coord_dims)]

    kd_start_time = get_time()
    call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                        workspace, value_buffer, permutation, left_stack, right_stack, &
                        recursion_stack, ierr)
    kd_end_time = get_time()
    timing_kd_build = timing_kd_build + (kd_end_time - kd_start_time)
    if (.not. is_ok(ierr)) return

    queries_start_time = get_time()
    timing_query_count = n_points

    do i = 1, n_points
      call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                        coords(:, i), k_neighbors, local_neighbors, local_distances, &
                        kd_workspace, query_ierr)

      ! Sigma local = 3 * std(dist(2..k)) (Asis requirement)
      if (k_neighbors > 1) then
        local_mean_dist = sum(local_distances(2:k_neighbors)) / real(k_neighbors-1, real64)
        local_variance  = sum((local_distances(2:k_neighbors) - local_mean_dist)**2) / &
                          real(k_neighbors-1, real64)
        if (local_variance < 0.0_real64) local_variance = 0.0_real64
        local_std_dev   = sqrt(local_variance)
      else
        local_std_dev   = 0.0_real64
      end if

      local_sigma = 3.0_real64 * local_std_dev

      if (local_sigma <= 0.0_real64) then
        min_pos = huge(1.0_real64)
        do j = 2, k_neighbors
          if (local_distances(j) > 0.0_real64 .and. local_distances(j) < min_pos) then
            min_pos = local_distances(j)
          end if
        end do
        if (min_pos < huge(1.0_real64)) then
          local_sigma = 3.0_real64 * min_pos
        else
          local_sigma = 1.0e-12_real64
        end if
      end if

      inv_two_sigma2 = 0.5_real64 / (local_sigma * local_sigma)

      work(:) = 0.0_real64
      wsum    = 0.0_real64

      do j = 1, k_neighbors
        neighbor_valid = .true.
        is_zero        = .true.

        do d = 1, n_vector_dims
          if (ieee_is_nan(vectors(d, local_neighbors(j)))) then
            neighbor_valid = .false.
            exit
          end if
          if (abs(vectors(d, local_neighbors(j))) > 1.0e-15_real64) is_zero = .false.
        end do
        if (.not. neighbor_valid) cycle
        if (is_zero) cycle

        w    = exp( - (local_distances(j)**2) * inv_two_sigma2 )
        wsum = wsum + w
        do d = 1, n_vector_dims
          work(d) = work(d) + w * vectors(d, local_neighbors(j))
        end do
      end do

      if (wsum > 0.0_real64) then
        smoothed(:, i) = work(:) / wsum
      else
        do d = 1, n_vector_dims
          if (.not. ieee_is_nan(vectors(d, i))) then
            smoothed(d, i) = vectors(d, i)
          else
            smoothed(d, i) = 0.0_real64
          end if
        end do
      end if
    end do

    queries_end_time    = get_time()
    timing_knn_queries  = timing_knn_queries + (queries_end_time - queries_start_time)
    timing_gaussian_calc= 0.0_real64

  end subroutine smooth_vectors_gaussian_adaptive

end module knn_smoothing_nadaraya_watson

