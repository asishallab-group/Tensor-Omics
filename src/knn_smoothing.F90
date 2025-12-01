!> Module for KNN smoothing operations
!! This module provides functions and subroutines for performing KNN-based smoothing
!! with Gaussian kernels and adaptive bandwidth selection.
module knn_smoothing
  use kd_tree,      only: build_kd_index, kd_knn_query
  use iso_fortran_env, only: int32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use tox_errors,  only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
  implicit none

contains

  !> Validate K parameters for KNN smoothing
  !! This subroutine ensures that the provided K parameters are within valid ranges.
  subroutine validate_k_params(k_min, k_max, k_grid, n_k, n_points, ierr)
    !| Minimum number of neighbors
    integer(int32), intent(in)  :: k_min
    !| Maximum number of neighbors
    integer(int32), intent(in)  :: k_max
    integer(int32) :: n_k,n_points
    !| Array of K values
    integer(int32), intent(in)  :: k_grid(n_k)
    integer(int32), intent(out) :: ierr
    integer(int32) :: kk

    call set_ok(ierr)

    if (k_min < 1 .or. k_max < 1 .or. k_min > k_max .or. k_max > n_points) then
       call set_err_once(ierr, ERR_INVALID_INPUT)
       return
    end if

    if (n_k < 1) then
       call set_err_once(ierr, ERR_INVALID_INPUT)
       return
    end if

    do kk = 1, n_k
       if (k_grid(kk) < k_min .or. k_grid(kk) > k_max) then
          call set_err_once(ierr, ERR_INVALID_INPUT)
          return
       end if
       if (kk > 1) then
          if (k_grid(kk) <= k_grid(kk-1)) then
             call set_err_once(ierr, ERR_INVALID_INPUT)
             return
          end if
       end if
    end do

    if (k_grid(1) /= k_min) then
       call set_err_once(ierr, ERR_INVALID_INPUT)
       return
    end if
  end subroutine validate_k_params

  !> Compute the Gaussian kernel bandwidth (sigma) for a given K
  !! This subroutine calculates the bandwidth based on the distance to the K-th neighbor.
  subroutine compute_sigma_for_k(distances, k, sigma_factor, sigma)
    !| Array of distances to neighbors
    real(real64),  intent(in)  :: distances(:)
    !| Number of neighbors
    integer(int32), intent(in) :: k
    !| Scaling factor for sigma
    real(real64),  intent(in)  :: sigma_factor
    !| Computed bandwidth
    real(real64),  intent(out) :: sigma

    if (k <= 0 .or. k > size(distances)) then
       sigma = 0.0_real64
       return
    end if

    sigma = abs(sigma_factor * distances(k))
    if (sigma <= 1.0e-18_real64) sigma = 1.0e-18_real64
  end subroutine compute_sigma_for_k

  !> Perform Gaussian KNN smoothing for a single K
  !! This subroutine smooths a point using a Gaussian kernel with a fixed number of neighbors.
  subroutine smooth_point_with_k(neighbors, distances, vectors, n_vector_dims, n_points, &
                               k, sigma_factor, result)
    !| Indices of neighbors
    integer(int32), intent(in) :: neighbors(:)
    !| Distances to neighbors
    real(real64),  intent(in) :: distances(:)
    !| Input vectors to smooth
    real(real64),  intent(in) :: vectors(n_vector_dims, n_points)
    !| Number of dimensions in the vectors
    integer(int32), intent(in) :: n_vector_dims
    !| Total number of points
    integer(int32), intent(in) :: n_points
    !| Number of neighbors
    integer(int32), intent(in) :: k
    !| Scaling factor for sigma
    real(real64),  intent(in) :: sigma_factor
    !| Smoothed vector
    real(real64),  intent(out):: result(n_vector_dims)

    real(real64) :: sigma, w, denom
    real(real64) :: sumvec(n_vector_dims)
    integer(int32) :: n, idx

    if (k <= 0) then
       result(:) = 0.0_real64
       return
    end if

    call compute_sigma_for_k(distances, k, sigma_factor, sigma)

    sumvec = 0.0_real64
    denom  = 0.0_real64

    do n = 1, k
       idx = neighbors(n)
       w   = exp( - (distances(n)*distances(n)) / (2.0_real64 * sigma*sigma) )
       sumvec(:) = sumvec(:) + w * vectors(:, idx)
       denom     = denom     + w
    end do

    if (denom > 0.0_real64) then
       result(:) = sumvec(:) / denom
    else
       result(:) = vectors(:, neighbors(1))
    end if
  end subroutine smooth_point_with_k

  !> Compute the second derivative for curvature calculation
  !! This subroutine calculates the second derivative of a sequence of vectors.
  subroutine compute_second_derivative(current, prev, prev2, a)
    !| Current vector
    real(real64), intent(in)  :: current(:)
    !| Previous vector
    real(real64), intent(in)  :: prev(:)
    !| Vector before the previous one
    real(real64), intent(in)  :: prev2(:)
    !| Computed second derivative
    real(real64), intent(out) :: a(:)

    a = current - 2.0_real64*prev + prev2
  end subroutine compute_second_derivative

  !> Compute curvature for a given K
  !! This subroutine calculates the curvature of smoothed values for a specific K.
  subroutine compute_curvature_for_k(t_hat_k, t_prev, t_prev2, C)
    !| Smoothed values for the current K
    real(real64), intent(in)  :: t_hat_k(:)
    !| Smoothed values for the previous K
    real(real64), intent(in)  :: t_prev(:)
    !| Smoothed values for the K before the previous one
    real(real64), intent(in)  :: t_prev2(:)
    !| Computed curvature
    real(real64), intent(out) :: C

    !| Temporary array to store second derivative
    real(real64) :: a(size(t_hat_k))

    if (any(ieee_is_nan(t_prev)) .or. any(ieee_is_nan(t_prev2))) then
       C = 0.0_real64
       return
    end if

    call compute_second_derivative(t_hat_k, t_prev, t_prev2, a)
    C = sum(a*a)
  end subroutine compute_curvature_for_k

  !> Compute normalized roughness values
  !! This subroutine calculates the normalized roughness values for a given set of curvatures.
  subroutine compute_normalized_roughness(C, n_k, Rtilde)
    !| Array of curvature values
    real(real64),  intent(in)  :: C(n_k)
    !| Number of K values
    integer(int32), intent(in) :: n_k
    !| Output array of normalized roughness values
    real(real64),  intent(out) :: Rtilde(n_k)

    real(real64) :: R0
    integer(int32) :: kk

    R0 = C(1)

    if (R0 <= 0.0_real64) then
       do kk = 1, n_k
          Rtilde(kk) = 1.0_real64
       end do
    else
       do kk = 1, n_k
          Rtilde(kk) = C(kk) / R0
       end do
    end if
  end subroutine compute_normalized_roughness

  !> Apply the elbow rule to determine the best K value
  !! This subroutine applies the elbow rule to identify the optimal K value based on normalized roughness.
  subroutine apply_elbow_rule(Rtilde, n_k, epsilon, best_k_idx)
    !| Array of normalized roughness values
    real(real64),  intent(in)  :: Rtilde(n_k)
    !| Number of K values
    integer(int32), intent(in) :: n_k
    !| Threshold for detecting the elbow
    real(real64),  intent(in)  :: epsilon
    !| Index of the best K value
    integer(int32), intent(out):: best_k_idx

    real(real64) :: delta
    integer(int32) :: kk

    best_k_idx = 1

    do kk = 2, n_k
        delta = Rtilde(kk-1) - Rtilde(kk)

        if (delta < epsilon) then
            best_k_idx = kk - 1
            return
        else
            best_k_idx = kk
        end if
    end do
  end subroutine apply_elbow_rule

  subroutine kd_request_neighbors(coords, kd_indices, n_coord_dims, n_points, &
                                  dimension_order, query_point, k_max, neighbors, distances, ierr)
    real(real64),  intent(in)  :: coords(n_coord_dims, n_points)
    integer(int32), intent(in) :: kd_indices(n_points)
    integer(int32), intent(in) :: n_coord_dims, n_points
    integer(int32), intent(in) :: dimension_order(n_coord_dims)
    real(real64),  intent(in)  :: query_point(n_coord_dims)
    integer(int32), intent(in) :: k_max
    integer(int32), intent(out):: neighbors(k_max)
    real(real64),  intent(out) :: distances(k_max)
    integer(int32), intent(out):: ierr

    real(real64)  :: kd_workspace(n_coord_dims)

    call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                      query_point, k_max, neighbors, distances, kd_workspace, ierr)
    ! Si ocurre error, se propaga en ierr
  end subroutine kd_request_neighbors

  !> Smooth vectors for multiple K values
  !! This subroutine performs Gaussian smoothing for multiple K values and computes curvatures.
  subroutine smooth_multiple_k(neighbors, distances, vectors, smoothed, &
                               n_vector_dims, n_points, i, k_grid, n_k, sigma_factor, t_hat_k, C)

    !| Indices of neighbors
    integer(int32), intent(in) :: neighbors(:)
    !| Distances to neighbors
    real(real64),  intent(in) :: distances(:)
    !| Total number of points
    integer(int32), intent(in) :: n_points
    !| Input vectors matrix
    real(real64),  intent(in) :: vectors(n_vector_dims, n_points)
    !| Number of vector dimensions
    integer(int32), intent(in) :: n_vector_dims
    !| Smoothed vectors matrix
    real(real64),  intent(in) :: smoothed(n_vector_dims, n_points)
    !| Current point index
    integer(int32), intent(in) :: i
    !| Array of K values
    integer(int32), intent(in) :: k_grid(n_k)
    !| Number of K values
    integer(int32), intent(in) :: n_k
    !| Scaling factor for Gaussian kernel width
    real(real64),  intent(in) :: sigma_factor
    !| Output smoothed vectors for each K
    real(real64),  intent(out):: t_hat_k(n_vector_dims, n_k)
    !| Output curvatures for each K
    real(real64),  intent(out):: C(n_k)

    integer(int32) :: kk, k

    do kk = 1, n_k
       k = k_grid(kk)
       call smooth_point_with_k(neighbors, distances, vectors, n_vector_dims, n_points, &
                         k, sigma_factor, t_hat_k(:, kk))

       call compute_curvature_for_k(t_hat_k(:, kk), smoothed(:, i-1), smoothed(:, i-2), C(kk))
    end do
  end subroutine smooth_multiple_k

  !===============================================================
  ! Ada-Carl main routine
  !===============================================================
subroutine smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
      n_coord_dims, n_vector_dims, n_points, &
      k_min, k_max, k_increment, sigma_factor, epsilon, &
      kd_indices, dimension_order, neighbors, distances, &
      workspace, value_buffer, permutation, left_stack, right_stack, ierr, &
      use_global_roughness)

    !| Number of coordinate dimensions
    integer(int32), intent(in) :: n_coord_dims
    !| Number of vector dimensions
    integer(int32), intent(in) :: n_vector_dims
    !| Total number of points
    integer(int32), intent(in) :: n_points
    !| Minimum number of neighbors
    integer(int32), intent(in) :: k_min
    !| Maximum number of neighbors
    integer(int32), intent(in) :: k_max
    !| Increment for generating K grid
    integer(int32), intent(in) :: k_increment
    !| Scaling factor for Gaussian kernel width
    real(real64),  intent(in) :: sigma_factor
    !| Threshold for elbow detection
    real(real64),  intent(in) :: epsilon
    !| Optional flag for using global roughness mode
    logical, intent(in), optional :: use_global_roughness
    !| Input coordinates matrix
    real(real64), intent(in)  :: coords(n_coord_dims, n_points)
    !| Input vectors matrix
    real(real64), intent(in)  :: vectors(n_vector_dims, n_points)
    !| Output smoothed vectors matrix
    real(real64), intent(out) :: smoothed(n_vector_dims, n_points)
    !| KD-tree indices
    integer(int32), intent(out)   :: kd_indices(n_points), dimension_order(n_coord_dims)
    !| Workspace arrays for KD-tree operations
    integer(int32), intent(inout) :: workspace(n_points), permutation(n_points)
    integer(int32), intent(inout) :: left_stack(n_points), right_stack(n_points)
    real(real64),  intent(inout)  :: value_buffer(n_points)
    !| Neighbor indices
    integer(int32), intent(inout) :: neighbors(k_max)
    !| Distances to neighbors
    real(real64),  intent(inout) :: distances(k_max)
    !| Error status
    integer(int32), intent(out)   :: ierr

    integer :: kk
    character(len=1), parameter :: TAB = achar(9)

    integer(int32) :: i, best_k_idx, best_k_global_idx, k_global
    real(real64), allocatable :: C(:), Rtilde(:)
    real(real64), allocatable :: C_all(:,:), R(:), Rtilde_global(:), t_hat_k(:, :)
    integer(int32) :: recursion_stack(3, n_points)
    logical :: global_mode
    integer(int32), allocatable :: k_grid(:)
    integer(int32) :: n_k

    call set_ok(ierr)

    ! Validate k_increment
    if (k_increment <= 0 .or. mod(k_max - k_min, k_increment) /= 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    ! Dynamically generate k_grid based on k_min, k_max, and k_increment
    n_k = 1 + (k_max - k_min) / k_increment
    allocate(k_grid(n_k))
    k_grid = [(k_min + (i-1)*k_increment, i=1,n_k)]

    ! Allocate memory for smoothing variables
    allocate(C(n_k), Rtilde(n_k), t_hat_k(n_vector_dims, n_k))

    ! Initialize smoothed array with default values
    smoothed = 0.0_real64

    ! Validate k parameters
    call validate_k_params(k_min, k_max, k_grid, n_k, n_points, ierr)
    if (.not. is_ok(ierr)) return

    dimension_order = [(i, i = 1, n_coord_dims)]
    call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                        workspace, value_buffer, permutation, left_stack, right_stack, &
                        recursion_stack, ierr)
    if (.not. is_ok(ierr)) return

    ! ============================
    ! seleccionar modo
    ! ============================
    global_mode = .false.
    if (present(use_global_roughness)) global_mode = use_global_roughness

    ! ================================================================
    !                 LOCAL MODE
    ! ================================================================
    if (.not. global_mode) then
        do i = 1, n_points
           call kd_request_neighbors(coords, kd_indices, n_coord_dims, n_points, &
                                     dimension_order, coords(:, i), k_max, neighbors, distances, ierr)
           if (.not. is_ok(ierr)) return

           if (i <= 2) then
              call smooth_point_with_k(neighbors, distances, vectors, &
                                       n_vector_dims, n_points, &
                                       k_min, sigma_factor, smoothed(:, i))
              cycle
           end if

           call smooth_multiple_k(neighbors, distances, vectors, smoothed, &
                                  n_vector_dims, n_points, i, k_grid, n_k, sigma_factor, t_hat_k, C)

           call compute_normalized_roughness(C, n_k, Rtilde)
           call apply_elbow_rule(Rtilde, n_k, epsilon, best_k_idx)

           ! Print best k for this point
           ! write(*,*) "Point:", i, "Best k:", k_grid(best_k_idx)

           ! Export values to check elbow method
           if (i == 1) then
               open(unit=13, file="curvature_values.tsv", status="replace")
               write(13, '(A)') "k"//TAB//"C"
           end if
           do kk = 1, n_k
               write(13, '(I10,A,F12.4)') k_grid(kk), TAB, C(kk)
           end do
           if (i == n_points) close(13)

           ! Export values for each k
           if (i == 1) then
               open(unit=14, file="smoothed_values_per_k.tsv", status="replace")
               write(14, '(A)') "Point_Index"//TAB//"k"//TAB//"Smoothed_Value"
           end if
           do kk = 1, n_k
               write(14, '(I10,A,I10,A,F12.4)') i, TAB, k_grid(kk), TAB, t_hat_k(1, kk)
           end do
           if (i == n_points) close(14)

           smoothed(:, i) = t_hat_k(:, best_k_idx)
        end do

        return
    end if


    ! ================================================================
    !                 GLOBAL MODE
    ! ================================================================
    call global_roughness_helper(n_points, n_k, C_all, R, Rtilde_global)

    ! 1) Collect curvatures for the entire series
    if (n_points > 0) then
        open(unit=13, file="curvature_values_global.tsv", status="replace")
        ! write(13, '(A)') "Point_Index"//TAB//"k"//TAB//"C"
        open(unit=14, file="smoothed_values_global.tsv", status="replace")
        ! write(14, '(A)') "Point_Index"//TAB//"k"//TAB//"Smoothed_Value"
    end if

    do i = 1, n_points
       call kd_request_neighbors(coords, kd_indices, n_coord_dims, n_points, &
                                 dimension_order, coords(:, i), k_max, neighbors, distances, ierr)
       if (.not. is_ok(ierr)) return

       if (i <= 2) then
          C_all(i,:) = 0.0_real64
          cycle
       end if

       call smooth_multiple_k(neighbors, distances, vectors, smoothed, &
                              n_vector_dims, n_points, i, k_grid, n_k, sigma_factor, t_hat_k, C)
       C_all(i,:) = C(:)

       ! Save curvature and smoothed values
       do kk = 1, n_k
           ! write(13, '(I10,A,I10,A,F12.4)') i, TAB, k_grid(kk), TAB, C(kk)
           ! write(14, '(I10,A,I10,A,F12.4)') i, TAB, k_grid(kk), TAB, t_hat_k(1, kk)
       end do
    end do

    if (n_points > 0) then
        close(13)
        close(14)
    end if

    ! 2) Apply global elbow rule
    call compute_global_roughness(C_all, n_points, n_k, R)
    call compute_normalized_roughness(R, n_k, Rtilde_global)
    call apply_elbow_rule(Rtilde_global, n_k, epsilon, best_k_global_idx)
    k_global = k_grid(best_k_global_idx)

    ! Print the best global k selected
    write(*,*) "Best global k selected:", k_global

    ! 3) Final smoothing using the single global k
    do i = 1, n_points
        call kd_request_neighbors(coords, kd_indices, n_coord_dims, n_points, &
                                  dimension_order, coords(:, i), k_max, neighbors, distances, ierr)
        call smooth_point_with_k(neighbors, distances, vectors, &
                                   n_vector_dims, n_points, &
                                   k_global, sigma_factor, smoothed(:, i))
    end do

end subroutine


  !> Helper for global roughness calculations
  !! This subroutine allocates memory for global roughness calculations.
  subroutine global_roughness_helper(n_points, n_k, C_all, R, Rtilde_global)
    !| Number of points
    integer(int32), intent(in) :: n_points
    !| Number of K values
    integer(int32), intent(in) :: n_k
    !| Output matrix for curvature values
    real(real64), allocatable, intent(out) :: C_all(:,:)
    !| Output array for roughness values
    real(real64), allocatable, intent(out) :: R(:)
    !| Output array for normalized roughness values
    real(real64), allocatable, intent(out) :: Rtilde_global(:)

    allocate(C_all(n_points, n_k))
    allocate(R(n_k), Rtilde_global(n_k))
  end subroutine global_roughness_helper

  !> Compute global roughness values
  !! This subroutine calculates the global roughness values for all points and K values.
  subroutine compute_global_roughness(C_all, n_points, n_k, R)
    !| Input matrix of curvature values
    real(real64), intent(in)  :: C_all(n_points, n_k)
    !| Number of points
    integer(int32), intent(in) :: n_points
    !| Number of K values
    integer(int32), intent(in) :: n_k
    !| Output array of roughness values
    real(real64), intent(out) :: R(n_k)

    integer(int32) :: kk, i
    R = 0.0_real64

    do kk = 1, n_k
      do i = 2, n_points
        R(kk) = R(kk) + C_all(i, kk)
      end do
    end do
  end subroutine compute_global_roughness

end module knn_smoothing
