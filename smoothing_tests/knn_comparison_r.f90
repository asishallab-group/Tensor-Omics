program test_knn_smoothing_with_module
  use iso_fortran_env, only: real64, int32
  use knn_smoothing,  only: smooth_vectors_gaussian_adaptive, reset_timing_stats, print_timing_stats
  implicit none

  integer(int32), parameter :: n_points = 500
  integer(int32), parameter :: n_coord_dims = 1
  integer(int32), parameter :: n_vector_dims = 1
  integer(int32), parameter :: max_k = 300

  real(real64) :: x(n_points)
  real(real64) :: y_true(n_points)
  real(real64) :: y_noisy(n_points)
  real(real64) :: shift   ! ⚠️ Ahora declarado

  real(real64) :: coords(n_coord_dims, n_points)
  real(real64) :: vectors(n_vector_dims, n_points)
  real(real64) :: smoothed(n_vector_dims, n_points)

  integer(int32) :: kd_indices(n_points), dimension_order(n_coord_dims)
  integer(int32) :: workspace(n_points), permutation(n_points)
  integer(int32) :: left_stack(n_points), right_stack(n_points)
  real(real64)   :: value_buffer(n_points)
  integer(int32) :: neighbors(max_k)
  real(real64)   :: distances(max_k)
  integer(int32) :: ierr

  integer(int32) :: i, k, unit
  real(real64) :: u1, u2, z

  integer(int32), parameter :: nk = 5
  integer(int32), parameter :: k_list(nk) = [25, 50, 100, 200, 300]

  character(len=64) :: fname
  character(len=1), parameter :: tab = char(9)


  call random_seed()

  ! ============================
  ! 1. Generar x
  ! ============================
  do i = 1, n_points
    call random_number(u1)
    x(i) = 10000.0_real64 * u1
  end do

  ! ============================
  ! 2. Curva real
  ! ============================
  do i = 1, n_points
    y_true(i) = 0.5_real64 * sqrt(x(i))
  end do

  ! ============================
  ! 3. Agregar ruido
  ! ============================
  do i = 1, n_points
    call random_number(u1)
    call random_number(u2)
    if (u1 < 1.0e-12_real64) u1 = 1.0e-12_real64
    z = sqrt(-2.0_real64 * log(u1)) * cos(2.0_real64 * acos(-1.0_real64) * u2)
    y_noisy(i) = y_true(i) + z * 200.0_real64
  end do

  ! ============================
  ! 3.5 Hacer todo positivo
  ! ============================
  shift = abs(minval(y_noisy)) + 1.0_real64
  do i = 1, n_points
    y_noisy(i) = y_noisy(i) + shift
  end do

  ! ============================
  ! Exportar los puntos originales
  ! ============================
  open(newunit=unit, file="synthetic_points.tsv", status="replace", action="write")
  write(unit,'(A)') "x" // tab // "y_noisy"
  do i = 1, n_points
    write(unit,'(F12.4,A,F12.4)') x(i), tab, y_noisy(i)
  end do
  close(unit)
  write(*,*) "Escribí archivo base synthetic_points.tsv"


    ! Copiar a coords/vectors en la forma requerida por el módulo
    do i = 1, n_points
        coords(1, i)  = x(i)
        vectors(1, i) = y_noisy(i)
    end do
  ! ============================
  ! 4. Loop KNN smoothing
  ! ============================
  do k = 1, nk
    call reset_timing_stats()

    call smooth_vectors_gaussian_adaptive( &
         coords        , &
         vectors       , &
         smoothed      , &
         n_coord_dims  , &
         n_vector_dims , &
         n_points      , &
         k_list(k)     , &
         kd_indices    , &
         dimension_order, &
         neighbors(1:k_list(k)), &
         distances(1:k_list(k)), &
         workspace     , &
         value_buffer  , &
         permutation   , &
         left_stack    , &
         right_stack   , &
         ierr)

    if (ierr /= 0) then
      write(*,*) "Error en smoothing, k=", k_list(k), " ierr=", ierr
    end if

    write(fname,'("knn_sim_k",I0,".tsv")') k_list(k)
    open(newunit=unit, file=fname, status='replace', action='write')
    write(unit,'(A)') "Mean_x" // tab // "Y_noisy" // tab // "Y_true" // tab // "Y_smoothed"

    do i = 1, n_points
      write(unit,'(A,F12.4,A,F12.4,A,F12.4,A,F12.4)') &
     "", x(i), tab, y_noisy(i), tab, y_true(i), tab, smoothed(1,i)
    end do

    close(unit)

    write(*,*) "Escribí archivo: ", trim(fname)
    call print_timing_stats(n_points, k_list(k))
  end do

end program test_knn_smoothing_with_module
