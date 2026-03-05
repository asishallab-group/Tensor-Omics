module knn_smoothing_nadaraya_watson
  use safeguard
  use kd_tree, only: build_kd_index, kd_knn_query
  use iso_fortran_env, only: int32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use tox_errors, only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
  use f42_utils, only: sort_real
  implicit none

contains

  ! ---------------------------------------------------------
  ! ADAPTIVE GAUSSIAN NADARAYA–WATSON SMOOTHING
  ! ---------------------------------------------------------
  subroutine smooth_vectors_gaussian_adaptive_nw(coords, vectors, smoothed, &
       n_coord_dims, n_vector_dims, n_points, k_neighbors, &
       kd_indices, dimension_order, neighbors, distances, &
       workspace, value_buffer, permutation, left_stack, right_stack, sigma_factor, ierr, smooth_type)

    integer(int32), intent(in)    :: n_coord_dims
    integer(int32), intent(in)    :: n_vector_dims
    integer(int32), intent(in)    :: n_points
    integer(int32), intent(in)    :: k_neighbors
    real(real64),  intent(in)     :: coords(n_coord_dims, n_points), sigma_factor
    real(real64),  intent(in)     :: vectors(n_vector_dims, n_points)
    real(real64),  intent(out)    :: smoothed(n_vector_dims, n_points)

    integer(int32), intent(out)   :: kd_indices(n_points), dimension_order(n_coord_dims)
    integer(int32), intent(inout) :: workspace(n_points), permutation(n_points)
    integer(int32), intent(inout) :: left_stack(n_points), right_stack(n_points)
    real(real64),  intent(inout)  :: value_buffer(n_coord_dims)

    integer(int32), intent(inout) :: neighbors(k_neighbors)
    real(real64),  intent(inout)  :: distances(k_neighbors)
    
    integer(int32), intent(out)   :: ierr
    integer(int32), intent(in), optional :: smooth_type ! 1 = Full, 2 = KNN-based

    ! General internal variables
    integer(int32) :: i, j, d
    real(real64)   :: kd_start_time, kd_end_time, queries_start_time, queries_end_time, dist_k
    integer(int32) :: permutation_distances(n_points)

    ! Variables for 1D (Full NW)
    real(real64)   :: all_distances(n_points)
    real(real64)   :: mean_dist, variance, std_dev
    real(real64)   :: min_pos, inv_two_sigma2
    real(real64)   :: wsum, w
    real(real64)   :: work(n_vector_dims)
    logical        :: neighbor_valid, is_zero
    integer(int32) :: count_nonzero
    real(real64)   :: effective_sigma    ! Effective local sigma

    ! Variables for N-D (previous path based on k-d tree + KNN)
    real(real64)   :: kd_workspace(n_coord_dims)
    integer(int32) :: recursion_stack(3, n_points), idx
    integer(int32) :: query_ierr
    integer(int32) :: local_neighbors(k_neighbors)
    real(real64)   :: local_distances(k_neighbors)
    real(real64)   :: local_mean_dist, local_variance, local_std_dev, local_sigma
    integer(int32) :: p_idx

    call set_ok(ierr)

    ! ----------------- Validation -----------------
    if (n_coord_dims < 1 .or. n_vector_dims < 1 .or. n_points < 1 .or. &
        k_neighbors < 1 .or. k_neighbors > n_points) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! ======================================================
    ! CASE 1D: Full Nadaraya–Watson (Processes ALL points)
    ! ======================================================
    if (smooth_type == 1) then

      do i = 1, n_points
        ! 1) Calculate distances to ALL points
        do j = 1, n_points
          all_distances(j) = abs(coords(1, i) - coords(1, j))
        end do

        ! 2) Local Sigma calculation: sigma_factor * std(distances excluding self)
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
          effective_sigma = sigma_factor * std_dev
        else
          effective_sigma = 0.0_real64
        end if

        ! Fallback if sigma is 0 (collapsed points)
        if (effective_sigma <= 0.0_real64) then
          min_pos = huge(1.0_real64)
          do j = 1, n_points
            if (j == i) cycle
            if (all_distances(j) > 0.0_real64 .and. all_distances(j) < min_pos) then
              min_pos = all_distances(j)
            end if
          end do
          if (min_pos < huge(1.0_real64)) then
            effective_sigma = 3.0_real64 * min_pos
          else
            effective_sigma = 1.0e-12_real64
          end if
        end if

        inv_two_sigma2 = 0.5_real64 / (effective_sigma * effective_sigma)

        ! 3) Nadaraya–Watson: Apply Gaussian weights to ALL points
        work(:) = 0.0_real64
        wsum    = 0.0_real64

        do j = 1, n_points
          if (j == i) cycle
          neighbor_valid = .true.
          is_zero        = .true.

          ! Validate vector data
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
          ! No valid neighbors: keep original value if not NaN
          do d = 1, n_vector_dims
            if (.not. ieee_is_nan(vectors(d, i))) then
              smoothed(d, i) = vectors(d, i)
            else
              smoothed(d, i) = 0.0_real64
            end if
          end do
        end if
      end do

      
    else 
      ! ======================================================
      ! CASE N-D: GEOMETRIC SMOOTHING (Adaptive k-NN Kernel)
      ! ======================================================
      dimension_order = [(i, i=1, n_coord_dims)]
      call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, &
                          recursion_stack, ierr )
      if (.not. is_ok(ierr)) return

      do i = 1, n_points
          ! 1. k-Nearest Neighbors Search
          call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                            coords(:,i), k_neighbors, neighbors, distances, &
                            query_ierr)

          ! ------------------------------------------------------
          ! 2. NEIGHBOR SORTING
          ! ------------------------------------------------------
          ! Initialize permutation vector to maintain distance-index mapping
          permutation_distances = [(j, j = 1, k_neighbors)]
          
          ! Sort distances in ascending order
          call sort_real(distances, permutation_distances, left_stack, right_stack)
          
          ! ------------------------------------------------------
          ! 3. ADAPTIVE SIGMA (Distance to the k-th neighbor)
          ! ------------------------------------------------------
          ! distances(permutation_distances(k_neighbors)) is now the maximum distance
          dist_k = distances(permutation_distances(k_neighbors))
          
          ! Define the bandwidth
          local_sigma = sigma_factor * dist_k
          
          ! Safety floor: if all neighbors are at the same location, sigma cannot be 0
          if (local_sigma < 1.0e-8_real64) local_sigma = 1.0e-8_real64
          inv_two_sigma2 = 0.5_real64 / (local_sigma**2)

          ! 4. WEIGHT ACCUMULATION (Nadaraya-Watson)
          work(:) = 0.0_real64
          wsum    = 0.0_real64
          
          do j = 1, k_neighbors
              ! Use permutation to get neighbor index and distance
              p_idx = permutation_distances(j) 
              idx   = neighbors(p_idx)

              if (idx == i) cycle
              if (idx <= 0) cycle
              
              ! Gaussian weight using correctly indexed distance
              w = exp(-(distances(p_idx)**2) * inv_two_sigma2)
              
              wsum = wsum + w
              do d = 1, n_vector_dims
                  work(d) = work(d) + w * vectors(d, idx)
              end do
          end do

          ! 5. FINAL RESULT 
          if (wsum > 1.0e-18_real64) then
              smoothed(:, i) = work(:) / wsum
          else
              smoothed(:, i) = vectors(:, i)
          end if
      end do
    end if

  end subroutine smooth_vectors_gaussian_adaptive_nw

end module knn_smoothing_nadaraya_watson