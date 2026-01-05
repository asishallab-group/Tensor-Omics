program test_knn_smoothing_with_module
  use iso_fortran_env, only: real64, int32
  use knn_smoothing,  only: smooth_vectors_gaussian_adaptive
  use f42_utils, only: sort_real
  implicit none

  integer(int32), parameter :: n_points = 20
  integer(int32), parameter :: n_coord_dims = 1
  integer(int32), parameter :: n_vector_dims = 1
  integer(int32), parameter :: k_min = 2
  integer(int32), parameter :: k_max = 8
  integer(int32), parameter :: k_increment = 1
  real(real64), parameter :: sigma_factor = 2.0_real64
  real(real64), parameter :: epsilon = 0.05_real64

  real(real64) :: coords(n_coord_dims, n_points)
  real(real64) :: vectors(n_vector_dims, n_points)
  real(real64) :: smoothed(n_vector_dims, n_points)

  integer(int32) :: kd_indices(n_points), dimension_order(n_coord_dims)
  integer(int32) :: workspace(n_points), permutation(n_points)
  integer(int32) :: left_stack(n_points), right_stack(n_points)
  real(real64)   :: value_buffer(n_points)
  integer(int32) :: neighbors(k_max)
  real(real64)   :: distances(k_max)
  integer(int32) :: ierr
  character(len=1), parameter :: TAB = achar(9)

  integer(int32) :: i

  ! Generate simple predefined data
  coords(1, :) = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, &
                   11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
  vectors(1, :) = [10.0, 9.0, 7.0, 6.0, 5.0, 5.0, 6.0, 7.0, 9.0, 10.0, &
                   10.0, 9.0, 7.0, 6.0, 5.0, 5.0, 6.0, 7.0, 9.0, 10.0]

  ! ===============================
  ! Sort coords and vectors using quicksort
  ! ===============================
  ! Debug: Verify sorting
  write(*,*) "Sorting data..."
  permutation = [(i, i = 1, n_points)]
  call sort_real(coords(1, :), permutation, left_stack, right_stack)
  coords(1, :) = coords(1, permutation)
  vectors(1, :) = vectors(1, permutation)
  write(*,*) "Data sorted. First 5 points:"
  do i = 1, min(5, n_points)
    write(*,*) "coords:", coords(1,i), "vectors:", vectors(1,i)
  end do

  ! ===============================
  !  CORRECT CALL
  ! ===============================
  ! Debug: Call to smooth_vectors_gaussian_adaptive
  write(*,*) "Calling smooth_vectors_gaussian_adaptive..."
  call smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
      n_coord_dims, n_vector_dims, n_points, &
      k_min, k_max, k_increment, sigma_factor, epsilon, &
      kd_indices, dimension_order, neighbors, distances, &
      workspace, value_buffer, permutation, left_stack, right_stack, ierr, &
      use_global_roughness = .true.)

  if (ierr /= 0) then
    write(*,*) "Error in smoothing! ierr=", ierr
  else
    write(*,*) "Smoothing completed successfully."
  end if

  ! Debug: Verify smoothed results
  write(*,*) "First 5 smoothed points:"
  do i = 1, min(5, n_points)
    write(*,*) "coords:", coords(1,i), "y_noisy:", vectors(1,i), "y_smooth:", smoothed(1,i)
  end do

  ! ===============================
  ! Export
  ! ===============================
  ! Debug: Confirm export
  write(*,*) "Exporting results to files..."
    open(unit=11, file="adaptive_knn.tsv", status="replace")
    write(11,'(A)') "x"//TAB//"y_noisy"//TAB//"y_smooth"
    do i = 1, n_points
        write(11,'(F10.4,A,F12.4,A,F12.4)') coords(1,i), TAB, vectors(1,i), TAB, smoothed(1,i)
    end do
    close(11)


  write(*,*) "Done! File: adaptive_knn.tsv"
end program
