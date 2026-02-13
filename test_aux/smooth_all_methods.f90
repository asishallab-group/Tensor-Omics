program smooth_all_methods_functional
    use iso_fortran_env, only: int32, real64
    use anwil,        only: smooth_vectors_gaussian_adaptive, anwil_smooth_sigma
    use manle_module, only: manle_pipeline, anwil_iterative, amanle_pipeline
    use tox_loess,   only: loess_alloc
    use knn_smoothing_nadaraya_watson, only: smooth_vectors_gaussian_adaptive_nw

    implicit none

    ! ============================
    ! GLOBAL PARAMETERS
    ! ============================
    real(real64),  parameter :: tol_manle    = 1.0e-6_real64

    ! LOESS
    integer(int32), parameter :: d_loess = 1
    integer(int32), parameter :: degree_loess = 2
    integer(int32), parameter :: liv_loess = 50
    integer(int32), parameter :: lv_loess  = 1000

    ! ============================
    ! VARIABLES
    ! ============================
    character(len=512) :: infile, outfile, line
    integer :: ios, arg_len, i
    integer(int32) :: n_points, k_neighbors, n_iters_max, method_id, kernel_type, k_neighbors_sigma
    real(real64) :: span_loess

    real(real64), allocatable :: x(:), y_orig(:)
    real(real64), allocatable :: y_loess(:), y_anwil(:), y_nw(:),y_anwil_iterative(:),x_anwil_iterative(:)


    ! -------- ANWIL --------
    real(real64), allocatable :: coords_anwil(:,:), vecs_anwil(:,:), smoothed_anwil(:,:), smoothed_anwil_mode1(:,:), smoothed_anwil_mode2(:,:), coords_nw(:,:), vecs_nw(:,:), smoothed_nw(:,:), smoothed_nw_knn(:,:), sd_arr(:,:), sd_arr_anwil(:,:), sigma_raw(:,:),sigma_raw_anwil(:,:)
    integer(int32), allocatable :: kd_indices(:), dimension_order(:)
    integer(int32), allocatable :: neighbors(:), workspace(:), permutation(:), permutation_distances(:)
    integer(int32), allocatable :: left_stack(:), right_stack(:)
    real(real64),    allocatable :: distances(:), value_buffer(:)

    ! -------- ManLe (1D functional) --------
    integer(int32), parameter :: k_manle = 1, p_manle = 5
    integer(int32), parameter :: l_manle = k_manle + p_manle
    real(real64), allocatable :: Omega(:,:), Y(:,:), Q(:,:), B(:,:)
    real(real64), allocatable :: Stmp(:), Utmp(:,:), tau(:), work_manle(:), work_amanle(:)
    real(real64), allocatable :: data2d(:,:), manifold2d(:,:)
    real(real64), allocatable :: y_manle_nocenter(:), y_manle_center(:), x_manle_nocenter(:), x_manle_center(:), x_amanle_center(:), y_amanle_center(:)
    real(real64), allocatable :: tmp_manle(:,:)
    real(real64), allocatable :: svd_line(:,:) ,svd_line_amanle(:,:) 

    ! -------- LOESS --------
    real(real64),  allocatable :: x_loess(:), w_loess(:), v_loess(:)
    integer(int32), allocatable :: iv_loess(:)
    


    integer(int32) :: ierr
    character(len=512) :: sd_arr_file

    ! ============================
    ! ARGUMENTS
    ! ============================
    call get_command_argument(1, infile, length=arg_len)
    if (arg_len == 0) then
        print *, "Usage:"
        print *, "  ./build/smooth_all  results/data/input.csv k_neighbors n_iters_max method k_neighbors_sigma kernel_type:optional span:optional"
        stop 1
    end if

    ! Read mandatory arguments for k_neighbors and n_iters_max
    call get_command_argument(2, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Error: k_neighbors is a mandatory argument."
        print *, "Usage:"
        print *, "  ./build/smooth_all  results/data/input.csv k_neighbors n_iters_max method k_neighbors_sigma kernel_type:optional span:optional"
        stop 1
    end if
    read(line, *) k_neighbors

    call get_command_argument(3, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Error: n_iters_max is a mandatory argument."
        print *, "Usage:"
        print *, "  ./build/smooth_all  results/data/input.csv k_neighbors n_iters_max method k_neighbors_sigma kernel_type:optional span:optional"
        stop 1
    end if
    read(line, *) n_iters_max

    ! Read argument to select method
    call get_command_argument(4, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Error: method_id is mandatory."
        print *, "Usage:"
        print *, "  ./build/smooth_all results/data/input.csv k_neighbors n_iters_max method k_neighbors_sigma kernel_type:optional span:optional"
        stop 1
    end if
    read(line, *) method_id

    ! Validate method_id
    if (method_id < 0 .or. method_id > 8) then
        print *, "Error: method_id should be between 0 and 7."
        stop 1
    end if

    ! Read argument to select k_neighbors_sigma
    call get_command_argument(5, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Error: k_neighbors_sigma is mandatory."
        print *, "Usage:"
        print *, "  ./build/smooth_all results/data/input.csv k_neighbors n_iters_max method k_neighbors_sigma kernel_type:optional span:optional"
        stop 1
    end if
    read(line, *) k_neighbors_sigma

    ! Read argument to select kernel_type
    call get_command_argument(6, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Warning: kernel_type not provided. Using default value: gaussian"
        kernel_type = 1_int32
    else
        read(line, *) kernel_type
    end if

    ! Read argument for span_loess
    call get_command_argument(7, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Warning: span_loess not provided. Using default value of 0.7."
        span_loess = 0.7_real64
    else
        read(line, *) span_loess
    end if

    infile = trim(infile)
    outfile = trim(infile)
    call make_output_name(outfile, k_neighbors, n_iters_max, 1, span_loess, k_neighbors_sigma, kernel_type)

    print *
    print *, "Reading:  ", trim(infile)
    print *, "Writing:", trim(outfile)
    print *, "k_neighbors:", k_neighbors
    print *, "n_iters_max:", n_iters_max
    print *, "span_loess:", span_loess
    print *, "kernel_type:", kernel_type

    ! ============================
    ! COUNT POINTS
    ! ============================
    open(10, file=trim(infile), status="old", action="read", iostat=ios)
    read(10,'(A)') line
    n_points = 0
    do
        read(10,'(A)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        n_points = n_points + 1
    end do
    close(10)

    print *
    print *, "Points:", n_points

    ! ============================
    ! ALLOCATE
    ! ============================
    allocate(x(n_points), y_orig(n_points))
    allocate(y_loess(n_points), y_anwil(n_points), y_nw(n_points),y_anwil_iterative(n_points),x_anwil_iterative(n_points))
    allocate(svd_line(2,n_points))  ! Asignar memoria para svd_line
    allocate(svd_line_amanle(2,n_points))

    ! ANWIL
    allocate(coords_anwil(2, n_points))
    allocate(sd_arr(2, n_points))
    allocate(sd_arr_anwil(2, n_points))
    allocate(sigma_raw(1, n_points))
    allocate(sigma_raw_anwil(1, n_points))
    allocate(vecs_anwil(2, n_points))
    allocate(smoothed_anwil(2, n_points))
    allocate(smoothed_anwil_mode1(2, n_points))
    allocate(smoothed_anwil_mode2(2, n_points))

    allocate(coords_nw(1, n_points))
    allocate(vecs_nw(1, n_points))
    allocate(smoothed_nw(1, n_points))
    allocate(smoothed_nw_knn(1, n_points))

    allocate(kd_indices(n_points), dimension_order(2))
    allocate(neighbors(k_neighbors), workspace(n_points), permutation(n_points), permutation_distances(k_neighbors))
    allocate(left_stack(n_points), right_stack(n_points))
    allocate(distances(k_neighbors), value_buffer(n_points))

    ! ManLe
    allocate(Omega(n_points, l_manle))
    allocate(Y(1, l_manle), Q(1, l_manle))
    allocate(B(l_manle, n_points))
    allocate(Stmp(l_manle), Utmp(l_manle, l_manle))
    allocate(tau(l_manle))
    allocate(work_manle(max(2000, 10*n_points)))
    allocate(work_amanle(max(2000, 10*n_points)))
    allocate(y_manle_nocenter(n_points))
    allocate(y_manle_center(n_points))
    allocate(x_manle_nocenter(n_points))
    allocate(x_manle_center(n_points))
    allocate(x_amanle_center(n_points))
    allocate(y_amanle_center(n_points))

    allocate(data2d(2, n_points))
    allocate(manifold2d(2, n_points))
    allocate(tmp_manle(2, n_points))


    ! LOESS
    allocate(x_loess(n_points), w_loess(n_points))
    allocate(iv_loess(liv_loess), v_loess(lv_loess))

    ! ============================
    ! READ CSV
    ! ============================
    open(11, file=trim(infile), status="old", action="read")
    read(11,'(A)') line
    do i = 1, n_points
        read(11,*) x(i), y_orig(i)
    end do
    close(11)

    print *, "Number of points:", n_points

    ! ============================
    ! Execute methods according to method_id
    ! ============================
    ! ============================
    ! Update LOESS smoothing to use span_loess
    ! ============================
    print * 

    if (method_id == 0 .or. method_id == 1) then
        print *, "Starting LOESS smoothing with span_loess=", span_loess
        x_loess = x
        w_loess(:) = 1.0_real64

        call loess_alloc(x_loess, y_orig, span_loess, degree_loess, y_loess, 1, 3, ierr )

    end if


    if (method_id == 0 .or. method_id == 2) then
        print *, "Starting ANWIL isotropic smoothing..."
        coords_anwil(1,:) = x
        coords_anwil(2,:) = y_orig
        vecs_anwil(1,:)   = x
        vecs_anwil(2,:)   = y_orig

        call anwil_smooth_sigma( &
           coords_anwil, vecs_anwil, smoothed_anwil, &
            2, 2, n_points, k_neighbors, &
            kd_indices, dimension_order, &
            workspace, value_buffer, permutation, left_stack, right_stack, sigma_raw_anwil, sd_arr_anwil, 0, 1.0_real64, k_neighbors_sigma, kernel_type, ierr )
        ! print *, "sd_arr: ", sd_arr
        !print *, "DEBUG: Antes de subrutina - Punto 250 vecs_anwil:", vecs_anwil(1, 250)
        !call anwil_smooth_sigma( &
        !    coords_anwil, vecs_anwil, smoothed_anwil, &
        !   2, 2, n_points, k_neighbors, &
        !    kd_indices, dimension_order, workspace, value_buffer, permutation, left_stack, right_stack, sigma_raw, sd_arr_anwil, ierr )
        !print *, "sd_arr: ", sd_arr_anwil
        

        print *, "DEBUG: Después de subrutina"
        print *, "   - Punto 250 smoothed_anwil:", smoothed_anwil(1, 250)
        print *, "   - Punto 251 smoothed_anwil:", smoothed_anwil(1, 251)
        print *, "   - Diferencia (debe ser pequeña):", smoothed_anwil(1, 251) - smoothed_anwil(1, 250)
    end if

    if (method_id == 0 .or. method_id == 3) then
        print *, "Starting ANWIL anisotropic smoothing (mode 1)..."
        coords_anwil(1,:) = x
        coords_anwil(2,:) = y_orig
        vecs_anwil(1,:)   = x
        vecs_anwil(2,:)   = y_orig

        call anwil_smooth_sigma( &
            coords_anwil, vecs_anwil, smoothed_anwil_mode1, &
            2, 2, n_points, k_neighbors, &
            kd_indices, dimension_order, &
            workspace, value_buffer, permutation, left_stack, right_stack, sigma_raw, sd_arr, 1, 10.0_real64, k_neighbors_sigma, kernel_type, ierr )
    end if

    if (method_id == 0 .or. method_id == 4) then
        print *, "Starting ANWIL anisotropic smoothing (mode 2)..."
        coords_anwil(1,:) = x
        coords_anwil(2,:) = y_orig
        vecs_anwil(1,:)   = x
        vecs_anwil(2,:)   = y_orig

        call anwil_smooth_sigma( &
            coords_anwil, vecs_anwil, smoothed_anwil_mode2, &
            2, 2, n_points, k_neighbors, &
            kd_indices, dimension_order, &
            workspace, value_buffer, permutation, left_stack, right_stack, sigma_raw, sd_arr, 2, 0.0_real64, k_neighbors_sigma, kernel_type, ierr )
    end if

    if (method_id == 0 .or. method_id == 5) then
        print *, "Starting Nadaraya–Watson smoothing..."
        coords_nw(1,:) = x
        vecs_nw(1,:)   = y_orig

        call smooth_vectors_gaussian_adaptive_nw( &
            coords_nw, vecs_nw, smoothed_nw, &
            1, 1, n_points, k_neighbors, &
            kd_indices, dimension_order, neighbors, distances, &
            workspace, value_buffer, permutation, left_stack, right_stack, &
            span_loess, ierr )

        print *, "Ejecutando NW KNN Local (k = 50)..."
        call smooth_vectors_gaussian_adaptive_nw( &
            coords_nw, vecs_nw, smoothed_nw_knn, &
            1, 1, n_points, k_neighbors, & ! <--- n_coord_dims=2 y k_neighbors=50
            kd_indices, dimension_order, neighbors, distances, &
            workspace, value_buffer, permutation, left_stack, right_stack, &
            0.5_real64, ierr ) ! Usamos el sigma_factor de 3.0
    end if

    if (method_id == 0 .or. method_id == 6) then
        print *, "Starting ManLe smoothing..."
        data2d(1,:) = x
        data2d(2,:) = y_orig

        call manle_pipeline( &
            data2d, n_points, 2, &
            k_neighbors, n_iters_max, tol_manle, &
            Omega, Y, Q, B, Stmp, Utmp, tau, work_manle, size(work_manle), &
            manifold2d, svd_line, k_neighbors_sigma, kernel_type, ierr )

        x_manle_center(:) = manifold2d(1,:)
        y_manle_center(:) = manifold2d(2,:)
    end if

    if (method_id == 0 .or. method_id == 7) then
        print *, "Starting AManLe smoothing..."
        data2d(1,:) = x
        data2d(2,:) = y_orig

        call amanle_pipeline( &
            data2d, n_points, 2, &
            1, k_neighbors, n_iters_max, tol_manle, &
            work_amanle, size(work_amanle), &
            manifold2d, svd_line_amanle, k_neighbors_sigma, kernel_type, ierr )

        x_amanle_center(:) = manifold2d(1,:)
        y_amanle_center(:) = manifold2d(2,:)

        ! ============================
        ! Debugging after AManLe smoothing
        ! ============================
        ! print *, "Debugging after AManLe smoothing..."
        ! print *, "First 5 points of x_amanle_center:", x_amanle_center(1:min(5, n_points))
        ! print *, "First 5 points of y_amanle_center:", y_amanle_center(1:min(5, n_points))
        ! print *, "First 5 points of svd_line_amanle (1st dimension):", svd_line_amanle(1, 1:min(5, n_points))
        ! print *, "First 5 points of svd_line_amanle (2nd dimension):", svd_line_amanle(2, 1:min(5, n_points))

        ! ============================
        ! Debugging: Check buffers and variables before and after amanle_pipeline
        ! ============================
        ! print *, "Debug: Before amanle_pipeline"
        ! print *, "First 5 points of data2d (1st dimension):", data2d(1, 1:min(5, n_points))
        ! print *, "First 5 points of data2d (2nd dimension):", data2d(2, 1:min(5, n_points))
    end if


    if (method_id == 0 .or. method_id == 8) then
        print *, "Starting anwil iterative smoothing..."
        data2d(1,:) = x
        data2d(2,:) = y_orig

        call anwil_iterative( &
            data2d, n_points, 2, &
            k_neighbors, n_iters_max, tol_manle, &
            tmp_manle, manifold2d, k_neighbors_sigma, kernel_type, ierr )

        y_anwil_iterative(:) = manifold2d(2,:)
        x_anwil_iterative(:) = manifold2d(1,:)
    end if



    ! ============================
    ! WRITE FINAL CSV
    ! ============================
    open(12, file=trim(outfile), status="replace", action="write")
    write(12,'(A)') 'x_original,y_original,x_loess,y_loess,x_anwil,y_anwil,x_anwil_mode1,y_anwil_mode1,x_anwil_mode2,y_anwil_mode2,x_anwil_iterative,y_anwil_iterative,x_nw,y_nw,y_nw_knn,x_manle,y_manle,x_manle_svd,y_manle_svd,x_amanle,y_amanle'

    do i = 1, n_points
        write(12,'(F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6)') &
            x(i), y_orig(i), x_loess(i), y_loess(i), smoothed_anwil(1,i), smoothed_anwil(2,i), smoothed_anwil_mode1(1,i), smoothed_anwil_mode1(2,i), smoothed_anwil_mode2(1,i), smoothed_anwil_mode2(2,i), x_anwil_iterative(i), y_anwil_iterative(i), coords_nw(1,i), smoothed_nw(1,i), smoothed_nw_knn(1,i),x_manle_center(i), y_manle_center(i), svd_line(1,i), svd_line(2,i), x_amanle_center(i), y_amanle_center(i)
    end do
    close(12)

    print *
    print *, "Done:", trim(outfile)

    ! ============================
    ! WRITE sd_arr TO TSV FILE
    ! ============================
    

    ! Create the filename based on outfile
    call make_output_name(outfile,k_neighbors,n_iters_max,2, span_loess, k_neighbors_sigma, kernel_type)

    open(20, file=outfile, status="replace", action="write")
    write(20, '(A)') 'x,local_sigma_raw,local_sigma_smooth'

    do i = 1, n_points
        write(20, '(F12.6,",",F12.6,",",F12.6)') x(i), sigma_raw(1,i), sd_arr_anwil(1,i)
    end do

    close(20)
    print *, "sd_arr_anwil saved to", trim(outfile)

contains

    subroutine make_output_name(name, k_neighbors, n_iters_max, option, span, k_neighbors_sigma, kernel_type)
        character(len=*), intent(inout) :: name
        integer(int32), intent(in) :: k_neighbors, n_iters_max
        real(real64), intent(in) :: span
        integer :: dotpos, i, option
        character(len=32) :: k_str, iter_str, span_str, kernel_str, sigma_str
        integer(int32), intent(in) :: kernel_type ! 1 gaussian, 2 tricubic
        integer(int32), intent(in) :: k_neighbors_sigma

        dotpos = 0
        do i = len_trim(name), 1, -1
            if (name(i:i) == '.') then
                dotpos = i
                exit
            end if
        end do

        ! Convert k_neighbors, n_iters_max, and span to chars
        write(k_str, '(I0)') k_neighbors
        write(iter_str, '(I0)') n_iters_max
        write(span_str, '(F4.2)') span
        write(kernel_str, '(I0)') kernel_type
        write(sigma_str, '(I0)') k_neighbors_sigma

        if (dotpos > 0) then
            if (option == 1) then 
                name = trim(name(1:dotpos-1)) // '_smoothed_k' // trim(k_str) // '_iter' // trim(iter_str) // '_span' // trim(adjustl(span_str)) // '_ksigma' // trim(sigma_str) // '_kernel' // trim(kernel_str) //'.csv'
            end if
            if (option == 2) then
                name = trim(name(1:dotpos-1)) // '_anwil_std.csv'
            end if 
        else
            if (option == 1) then 
                name = trim(name) // '_smoothed_k' // trim(k_str) // '_iter' // trim(iter_str) // '_span' // trim(adjustl(span_str)) // '.csv'
            end if 
            if (option == 2) then
                name = trim(name) // '_anwil_std.csv'
            end if 

        end if
    end subroutine make_output_name

end program smooth_all_methods_functional
