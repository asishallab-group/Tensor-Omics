program benchmark_smoothing
  use iso_fortran_env, only: real64, int32
  ! Asegúrate de que knn_smoothing tenga los módulos y subrutinas correctas
  use knn_smoothing, only: smooth_vectors_gaussian_adaptive, reset_timing_stats, get_time
  implicit none

  integer(int32), parameter :: n_coord_dims = 1
  integer(int32), parameter :: n_vector_dims = 1
  
  ! ==========================================================
  ! DEFINICIÓN DE LOS ESCENARIOS DE PRUEBA
  ! ==========================================================
  
  ! 1. Tamaños de N a probar
  integer(int32), parameter :: num_sizes = 6
  integer(int32), parameter :: sizes(num_sizes) = &
        [10000, 50000, 100000, 200000, 500000, 1000000]

  ! 2. Regímenes de K a probar para cada N
  integer(int32), parameter :: num_k_regimes = 3
  character(len=8), parameter :: k_regime_names(num_k_regimes) = &
       ["RELATIVE", "FIXED   ", "ROOT    "]
  
  ! Valores específicos para los regímenes
  real(real64), parameter :: relative_frac = 0.001_real64 ! k = 0.1% de N (Tu prueba original)
  integer(int32), parameter :: fixed_k_value = 30      ! k = 30 (Para validar O(N log N))
  
  ! ==========================================================
  
  character(len=128) :: fname
  character(len=1), parameter :: tab = char(9)

  integer(int32) :: s, r, i
  integer(int32) :: N, unit, k_neighbors, ierr
  real(real64) :: t1, t2, elapsed
  real(real64) :: u1, u2

  ! Datos
  real(real64), allocatable :: coords(:,:), vectors(:,:), smoothed(:,:)

  ! Buffers para el módulo KNN (se asignarán con el N actual)
  integer(int32), allocatable :: kd_indices(:), dimension_order(:)
  integer(int32), allocatable :: workspace(:), permutation(:)
  integer(int32), allocatable :: left_stack(:), right_stack(:), neighbors(:)
  real(real64), allocatable :: value_buffer(:), distances(:)

  ! Archivo de resultados (tiempos)
  open(11, file="benchmark_results_full.tsv", status='replace', action='write')
  write(11,'(A)') "N"//tab//"k"//tab//"k_regime"//tab//"seconds"

  call random_seed()

  ! ==========================================================
  ! LOOP EXTERNO: TAMAÑOS DEL DATASET (N)
  ! ==========================================================
  do s = 1, num_sizes
    N = sizes(s)
    write(*,*) ">>>> STARTING N=", N

    ! ----------------------------
    ! 1. PREPARACIÓN DE MEMORIA
    ! ----------------------------
    ! Liberar si ya estaban alloc
    if (allocated(coords))     deallocate(coords)
    if (allocated(vectors))   deallocate(vectors)
    if (allocated(smoothed))    deallocate(smoothed)

    if (allocated(kd_indices))   deallocate(kd_indices)
    if (allocated(dimension_order)) deallocate(dimension_order)
    if (allocated(workspace))   deallocate(workspace)
    if (allocated(permutation)) deallocate(permutation)
    if (allocated(left_stack))   deallocate(left_stack)
    if (allocated(right_stack)) deallocate(right_stack)
    if (allocated(value_buffer))  deallocate(value_buffer)
    if (allocated(neighbors))   deallocate(neighbors)
    if (allocated(distances))   deallocate(distances)

    ! Alloc nuevo para este N
    allocate(coords(n_coord_dims, N))
    allocate(vectors(n_vector_dims, N))
    allocate(smoothed(n_vector_dims, N))

    allocate(kd_indices(N))
    allocate(dimension_order(n_coord_dims))
    allocate(workspace(N))
    allocate(permutation(N))
    allocate(left_stack(N))
    allocate(right_stack(N))
    allocate(value_buffer(N))

    ! Vecinos/distancias: tamaño máximo posible = N (aunque usamos sub-arrays)
    allocate(neighbors(N))
    allocate(distances(N))

    ! ----------------------------
    ! 2. GENERACIÓN DE DATOS
    ! ----------------------------
    do i = 1, N
      ! coords ~ U(0,10000)
      call random_number(u1)
      coords(1,i) = 10000.0_real64 * u1

      ! vectors ~ algo positivo (p.ej. 300..800)
      call random_number(u2)
      vectors(1,i) = 300.0_real64 + 500.0_real64 * u2
    end do

    ! ==========================================================
    ! LOOP INTERNO: REGÍMENES DE K (FIXED   , RELATIVE, ROOT    )
    ! ==========================================================
    do r = 1, num_k_regimes

      ! 3. Determinar k y la fracción
      select case (trim(k_regime_names(r)))
        case ("RELATIVE")
          k_neighbors = int(relative_frac * real(N, real64))
        case ("FIXED   ")
          k_neighbors = fixed_k_value
        case ("ROOT    ")
          k_neighbors = int(sqrt(real(N, real64)))
        case default
          k_neighbors = 1
      end select
      
      ! Asegurar límites de k (mínimo 1, máximo N)
      if (k_neighbors < 1) k_neighbors = 1
      if (k_neighbors > N) k_neighbors = N

      ! 4. EJECUCIÓN DEL SMOOTHING
      call reset_timing_stats()
      t1 = get_time()

      call smooth_vectors_gaussian_adaptive( &
            coords, vectors, smoothed, &
            n_coord_dims, n_vector_dims, N, &
            k_neighbors, &
            kd_indices, dimension_order, &
            neighbors(1:k_neighbors), distances(1:k_neighbors), &
            workspace, value_buffer, permutation, &
            left_stack, right_stack, ierr)

      t2   = get_time()
      elapsed = t2 - t1

      ! 5. LOGGING Y ESCRITURA DE RESULTADOS
      write(*,'(" | K_Regime=",A,T25,"k=",I0," time=",F10.4," sec")') &
            trim(k_regime_names(r)), k_neighbors, elapsed

      write(11,'(I0,A,I0,A,A,A,F10.4)') &
            N, tab, k_neighbors, tab, trim(k_regime_names(r)), tab, elapsed

      ! 6. Volcado de curvas para análisis visual (solo para N pequeños)
      if (N <= 200000 .and. trim(k_regime_names(r)) == "FIXED   ") then
        ! Solo volcamos el caso FIXED    para visualizar el comportamiento O(N log N)
        write(fname,'("smoothing_N",I0,"_k",I0,"_",A,".tsv")') &
              N, k_neighbors, trim(k_regime_names(r))
        open(newunit=unit, file=trim(fname), status='replace', action='write')
        write(unit,'(A)') "x"//tab//"y_original"//tab//"y_smoothed"
        do i = 1, N
          write(unit,'(F12.4,A,F12.4,A,F12.4)') coords(1,i), tab, vectors(1,i), tab, smoothed(1,i)
        end do
        close(unit)
      end if

    end do ! r (regímenes de k)

    ! ----------------------------
    ! 7. DESPUÉS DE TODOS LOS K's
    ! ----------------------------
    ! Deallocate después de probar todos los K's para un N dado.
    deallocate(coords, vectors, smoothed)
    deallocate(kd_indices, dimension_order, workspace, permutation)
    deallocate(left_stack, right_stack, value_buffer, neighbors, distances)

  end do  ! s (tamaños de N)

  close(11)
  write(*,*) "DONE — resultados en benchmark_results_full.tsv"
end program benchmark_smoothing