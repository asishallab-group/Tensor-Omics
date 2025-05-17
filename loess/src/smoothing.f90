!> @file smoothing.f90
!> @brief Functionality for LOESS smoothing of data.
!> @details This module implements the LOESS smoothing algorithm, which is a non-parametric regression method.


module loess_module
  implicit none
contains

  !> Finds the indices of the true values in a logical mask.
  !! @param mask Logical array of size n.
  !! @param n Size of the mask.
  !! @param idx_out Integer array to store the indices of true values.
  !! @param m_max Maximum size of idx_out.
  !! @param m_out Actual size of idx_out (number of true values found).
  subroutine which(mask, n, idx_out, m_max, m_out)  
    logical, intent(in) :: mask(n)
    integer, intent(out) :: idx_out(m_max)
    integer, intent(out) :: m_out
    integer :: i, count, m_max, n
    count = 0
    idx_out = 0  ! Initialize to avoid garbage values
    do i = 1, n
      if (mask(i)) then
        count = count + 1
        if (count <= m_max) then
          idx_out(count) = i
        end if
      end if
    end do
    m_out = count
  end subroutine which

  

  !> Performs LOESS smoothing on a set of data points.
  !! @param n_total Total number of reference points.
  !! @param n_target Number of target points to smooth.
  !! @param d Dimensionality of the data.
  !! @param x_ref Reference x-coordinates.
  !! @param y_ref Reference y-coordinates (d x n_total).
  !! @param indices_used Indices of reference points used for smoothing.
  !! @param x_query Target x-coordinates to smooth.
  !! @param kernel_sigma Bandwidth parameter for the kernel.
  !! @param kernel_cutoff Cutoff for the kernel.
  !! @param y_out Output smoothed values (d x n_target).
  !! @param workspace_weights Temporary array for weights.
  !! @param workspace_values Temporary array for values.
  !! @return y_out Smoothed values.
  subroutine loess_smooth(n_total, n_target, d, x_ref, y_ref, indices_used, x_query, &
                          kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                          workspace_values)
    integer, intent(in) :: n_total, n_target, d
    real, intent(in) :: x_ref(n_total), y_ref(d, n_total), x_query(n_target)
    integer, intent(in) :: indices_used(n_target)
    real, intent(in) :: kernel_sigma, kernel_cutoff
    real, intent(out) :: y_out(d, n_target)
    real, intent(inout) :: workspace_weights(n_total), workspace_values(d, n_total)

    integer :: q, i, k, idx, m_out, valid_indices(n_target)
    real :: query_x, ref_x, delta, sum_weights, weight
    logical :: mask(n_target)

    do q = 1, n_target
      query_x = x_query(q)
      sum_weights = 0.0
      y_out(:, q) = 0.0
      ! Create mask for valid indices within cutoff
      mask = .false.
      do i = 1, n_target
        idx = indices_used(i)
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        mask(i) = (delta <= kernel_cutoff * kernel_sigma)
      end do
      ! Get valid indices using which
      call which(mask, n_target, valid_indices, n_target, m_out)

      ! Accumulate weights and values
      do i = 1, m_out
        idx = indices_used(valid_indices(i))
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        weight = exp(-(delta / kernel_sigma)**2)
        sum_weights = sum_weights + weight
        do k = 1, d
          y_out(k, q) = y_out(k, q) + weight * y_ref(k, idx)
        end do
      end do

      ! Normalize or fallback
      if (sum_weights > 0.0) then
        y_out(:, q) = y_out(:, q) / sum_weights
      else
        y_out(:, q) = y_ref(:, indices_used(q))
      end if
    end do

  end subroutine loess_smooth

  !> End of the loess_module.
end module loess_module

!> Program to test the LOESS smoothing implementation.
program test_loess
  use loess_module
  implicit none

  integer, parameter :: n_total = 100, n_target = 50, d = 1
  real :: x_ref(n_total), y_ref(d, n_total), x_query(n_target), y_out(d, n_target)
  integer :: indices_used(n_total)  ! Adjust size to match n_total
  real :: workspace_weights(n_total), workspace_values(d, n_total)
  real :: kernel_sigma, kernel_cutoff
  integer :: i, j

  ! Variables for Test Case 5 (vector field smoothing)
  integer, parameter :: d_vec = 3
  real :: y_ref_vec(d_vec, n_total), y_out_vec(d_vec, n_target)
  real :: workspace_values_vec(d_vec, n_total)

  !> Test Case 1: Constant Input
  x_ref = 5.0
  y_ref = 10.0
  x_query = 5.0
  indices_used = (/ (i, i = 1, n_total) /)  ! Include all reference points
  kernel_sigma = 1.0
  kernel_cutoff = 3.0
  call loess_smooth(n_total, n_target, d, x_ref, y_ref, indices_used, x_query, &
                    kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                    workspace_values)
  print *, 'Test Case 1 Error: ', abs(y_out - 10.0)
  print *, 'Test Case 1: ', all(abs(y_out - 10.0) < 1.0e-6)

  !> Test Case 2: Linear Trend Recovery
  x_ref = (/ (real(i), i = 1, n_total) /)
  y_ref(1, :) = 0.5 * x_ref
  x_query = (/ (real(i) + 0.5, i = 1, n_target) /)
  indices_used = (/ (i, i = 1, n_total) /)  ! Include all reference points
  kernel_sigma = 1.0  ! Adjusted for a wider kernel
  kernel_cutoff = 3.0
  call loess_smooth(n_total, n_target, d, x_ref, y_ref, indices_used, x_query, &
                    kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                    workspace_values)
  print *, 'Test Case 2 Error: ', abs(y_out(1, :) - 0.5 * x_query)
  print *, 'Test Case 2: ', all(abs(y_out(1, :) - 0.5 * x_query) < 0.05)

  !> Test Case 3: Outlier Suppression
  x_ref = (/ (10.0, i = 1, n_total - 1), 100.0 /)  ! Adjusted to match n_total
  y_ref(1, :) = (/ (5.0, i = 1, n_total - 1), 99.0 /)  ! Adjusted to match n_total
  x_query = (/ (10.0, i = 1, n_target) /)  ! Adjusted to match n_target
  indices_used = (/ (i, i = 1, n_total) /)  ! Include all reference points
  kernel_sigma = 1.0
  kernel_cutoff = 3.0
  call loess_smooth(n_total, n_target, d, x_ref, y_ref, indices_used, x_query, &
                    kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                    workspace_values)
  print *, 'Test Case 3 Error: ', abs(y_out(1, 1) - 5.0)
  print *, 'Test Case 3: ', abs(y_out(1, 1) - 5.0) < 0.01

  !> Test Case 4: Sparse Fallback Behavior
  x_ref = (/ (real(i) * 100.0, i = 1, n_total) /)
  y_ref(1, :) = (/ (real(i), i = 1, n_total) /)
  x_query = (/ (real(i) * 100.0 + 50.0, i = 1, n_target) /)
  indices_used = (/ (i, i = 1, n_total) /)  ! Include all reference points
  kernel_sigma = 1.0
  kernel_cutoff = 3.0
  call loess_smooth(n_total, n_target, d, x_ref, y_ref, indices_used, x_query, &
                    kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                    workspace_values)
  print *, 'Test Case 4 Error: ', abs(y_out(1, :) - y_ref(1, indices_used(1:n_target)))
  print *, 'Test Case 4: ', all(abs(y_out(1, :) - y_ref(1, indices_used(1:n_target))) < 1.0e-6)

  !> Test Case 5: Vector Field Smoothing
  do j = 1, d_vec
    do i = 1, n_total
      y_ref_vec(j, i) = real(i) + 0.1 * real(j)
    end do
  end do
  x_ref = (/ (real(i), i = 1, n_total) /)
  x_query = (/ (real(i) + 0.5, i = 1, n_target) /)
  indices_used = (/ (i, i = 1, n_total) /)  ! Include all reference points
  kernel_sigma = 1.0
  kernel_cutoff = 3.0
  call loess_smooth(n_total, n_target, d_vec, x_ref, y_ref_vec, indices_used, &
                    x_query, kernel_sigma, kernel_cutoff, y_out_vec, &
                    workspace_weights, workspace_values_vec)
  print *, 'Test Case 5 Error: ', abs(y_out_vec - y_ref_vec(:, indices_used(1:n_target)))
  print *, 'Test Case 5: ', all(abs(y_out_vec - y_ref_vec(:, indices_used(1:n_target))) < 0.1)

!> End of the test_loess program.
end program test_loess

!> Test Case 2 and Test Case 5 are expected to pass with a tolerance of 0.05 and 0.1 respectively but fail at the moment.
