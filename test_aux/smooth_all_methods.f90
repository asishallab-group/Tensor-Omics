program smooth_all_methods_functional
    use iso_fortran_env, only: int32, real64
    use anwil,           only: anwil_smooth_sigma
    use manle_module,    only: manle_pipeline, anwil_iterative, amanle_pipeline
    use tox_loess,       only: loess_alloc
    use knn_smoothing_nadaraya_watson, only: smooth_vectors_gaussian_adaptive_nw

    implicit none

    ! ============================
    ! GLOBAL PARAMETERS
    ! ============================
    real(real64),  parameter :: tol_manle    = 1.0e-6_real64

    ! LOESS Parameters
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
    integer    :: method_flag      ! 1 = Arithmetic, 2 = Geometric
    real(real64)  :: w_r, w_e, w_c    ! Weights for Roughness (R), Error (E), and Coverage (C)

    real(real64), allocatable :: x(:), y_orig(:)
    real(real64), allocatable :: y_loess(:), y_anwil(:), y_nw(:), y_anwil_iterative(:), x_anwil_iterative(:)

    ! -------- ANWIL Data --------
    real(real64), allocatable :: coords_anwil(:,:), vecs_anwil(:,:), smoothed_anwil(:,:), &
                                 coords_nw(:,:), vecs_nw(:,:), smoothed_nw(:,:), &
                                 smoothed_nw_knn(:,:), sd_arr(:,:), sd_arr_anwil(:,:), &
                                 sigma_raw(:,:), sigma_raw_anwil(:,:)
    integer(int32), allocatable :: kd_indices(:), dimension_order(:)
    integer(int32), allocatable :: neighbors(:), workspace(:), permutation(:), permutation_distances(:)
    integer(int32), allocatable :: left_stack(:), right_stack(:)
    real(real64),   allocatable :: distances(:), value_buffer(:)

    ! Optimization History
    integer(int32) :: stop_iter, stop_reason, best_iter
    real(real64)   :: best_score
    real(real64), allocatable  :: history_scores(:)
    real(real64), allocatable  :: history_roughness(:)
    real(real64), allocatable  :: history_rmse(:)
    real(real64), allocatable  :: history_coverage(:)
    real(real64), allocatable  :: history_penalty(:)
    integer(int32), parameter :: patience_k = 5_int32
    integer(int32), parameter :: min_iters  = 3_int32
    real(real64),   parameter :: tol_rel    = 5.0e-3_real64   ! 0.5%

    ! -------- ManLe / AManLe Data --------
    integer(int32), parameter :: k_manle = 1, p_manle = 5
    integer(int32), parameter :: l_manle = k_manle + p_manle
    real(real64), allocatable :: Omega(:,:), Y(:,:), Q(:,:), B(:,:)
    real(real64), allocatable :: Stmp(:), Utmp(:,:), tau(:), work_manle(:), work_amanle(:)
    real(real64), allocatable :: data2d(:,:), manifold2d(:,:)
    real(real64), allocatable :: y_manle_nocenter(:), y_manle_center(:), x_manle_nocenter(:), &
                                 x_manle_center(:), x_amanle_center(:), y_amanle_center(:)
    real(real64), allocatable :: tmp_manle(:,:)
    real(real64), allocatable :: svd_line(:,:), svd_line_amanle(:,:) 

    ! -------- LOESS Data --------
    real(real64),  allocatable :: x_loess(:), w_loess(:), v_loess(:)
    integer(int32), allocatable :: iv_loess(:)

    integer(int32) :: ierr

    ! ============================
    ! COMMAND LINE ARGUMENTS
    ! ============================
    call get_command_argument(1, infile, length=arg_len)
    if (arg_len == 0) then
        print *, "Usage:"
        print *, "  ./build/smooth_all input.csv k_neighbors n_iters_max method_id k_sigma kernel_type span score_type w_r w_e w_c"
        stop 1
    end if

    ! Argument parsing (2: k_neighbors)
    call get_command_argument(2, line, length=arg_len)
    if (arg_len == 0) stop "Error: k_neighbors is mandatory."
    read(line, *) k_neighbors

    ! Argument parsing (3: n_iters_max)
    call get_command_argument(3, line, length=arg_len)
    if (arg_len == 0) stop "Error: n_iters_max is mandatory."
    read(line, *) n_iters_max

    ! Argument parsing (4: method_id)
    call get_command_argument(4, line, length=arg_len)
    if (arg_len == 0) stop "Error: method_id is mandatory."
    read(line, *) method_id
    if (method_id < 0 .or. method_id > 8) stop "Error: method_id must be between 0 and 8."

    ! Argument parsing (5: k_neighbors_sigma)
    call get_command_argument(5, line, length=arg_len)
    if (arg_len == 0) stop "Error: k_neighbors_sigma is mandatory."
    read(line, *) k_neighbors_sigma

    ! Argument parsing (6: kernel_type)
    call get_command_argument(6, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Warning: kernel_type not provided. Using default: Gaussian (1)"
        kernel_type = 1_int32
    else
        read(line, *) kernel_type
    end if

    ! Argument parsing (7: span_loess)
    call get_command_argument(7, line, length=arg_len)
    if (arg_len == 0) then
        print *, "Warning: span_loess not provided. Using default: 0.7"
        span_loess = 0.7_real64
    else
        read(line, *) span_loess
    end if

    ! Scoring parameters (8: method_flag, 9-11: weights)
    call get_command_argument(8, line, length=arg_len)
    if (arg_len == 0) stop "Error: score_type (method_flag) is mandatory."
    read(line, *) method_flag

    call get_command_argument(9, line, length=arg_len)
    if (arg_len == 0) stop "Error: w_r is mandatory."
    read(line, *) w_r

    call get_command_argument(10, line, length=arg_len)
    if (arg_len == 0) stop "Error: w_e is mandatory."
    read(line, *) w_e

    call get_command_argument(11, line, length=arg_len)
    if (arg_len == 0) stop "Error: w_c is mandatory."
    read(line, *) w_c

    ! Setup output naming
    infile = trim(infile)
    outfile = trim(infile)
    call make_output_name(outfile, k_neighbors, n_iters_max, 1, span_loess, k_neighbors_sigma, kernel_type, method_flag, w_r, w_e, w_c)

    print *
    print *, "Input file:  ", trim(infile)
    print *, "Output file: ", trim(outfile)
    print *, "k_neighbors: ", k_neighbors
    print *, "n_iters_max: ", n_iters_max
    print *, "span_loess:  ", span_loess
    print *, "kernel_type: ", kernel_type

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

    print *, "Total points detected:", n_points

    ! ============================
    ! ALLOCATION BLOCK
    ! ============================
    allocate(x(n_points), y_orig(n_points))
    allocate(y_loess(n_points), y_anwil(n_points), y_nw(n_points), y_anwil_iterative(n_points), x_anwil_iterative(n_points))
    allocate(svd_line(2, n_points), svd_line_amanle(2, n_points))

    ! ANWIL / NW Buffers
    allocate(coords_anwil(2, n_points), sd_arr(2, n_points), sd_arr_anwil(2, n_points))
    allocate(sigma_raw(1, n_points), sigma_raw_anwil(1, n_points), vecs_anwil(2, n_points), smoothed_anwil(2, n_points))
    allocate(coords_nw(1, n_points), vecs_nw(1, n_points), smoothed_nw(1, n_points), smoothed_nw_knn(1, n_points))
    
    ! History buffers
    allocate(history_scores(n_iters_max), history_coverage(n_iters_max), history_penalty(n_iters_max), &
             history_rmse(n_iters_max), history_roughness(n_iters_max))

    ! Tree/Search buffers
    allocate(kd_indices(n_points), dimension_order(2))
    allocate(neighbors(k_neighbors), workspace(n_points), permutation(n_points), permutation_distances(k_neighbors))
    allocate(left_stack(n_points), right_stack(n_points))
    allocate(distances(k_neighbors), value_buffer(n_points))

    ! ManLe / AManLe Buffers
    allocate(Omega(n_points, l_manle), Y(1, l_manle), Q(1, l_manle), B(l_manle, n_points))
    allocate(Stmp(l_manle), Utmp(l_manle, l_manle), tau(l_manle))
    allocate(work_manle(max(2000, 10*n_points)), work_amanle(max(2000, 10*n_points)))
    allocate(y_manle_nocenter(n_points), y_manle_center(n_points), x_manle_nocenter(n_points), &
             x_manle_center(n_points), x_amanle_center(n_points), y_amanle_center(n_points))
    allocate(data2d(2, n_points), manifold2d(2, n_points), tmp_manle(2, n_points))

    ! LOESS Buffers
    allocate(x_loess(n_points), w_loess(n_points), iv_loess(liv_loess), v_loess(lv_loess))

    ! ============================
    ! READ CSV DATA
    ! ============================
    open(11, file=trim(infile), status="old", action="read")
    read(11,'(A)') line
    do i = 1, n_points
        read(11,*) x(i), y_orig(i)
    end do
    close(11)

    ! ============================
    ! EXECUTION PIPELINE
    ! ============================

    ! --- 1. LOESS ---
    if (method_id == 0 .or. method_id == 1) then
        print *, "Executing LOESS: span=", span_loess
        x_loess = x
        w_loess(:) = 1.0_real64
        call loess_alloc(x_loess, y_orig, span_loess, degree_loess, y_loess, 1, 3, ierr)
    end if

    ! --- 2. ANWIL Isotropic ---
    if (method_id == 0 .or. method_id == 2) then
        print *, "Executing ANWIL Isotropic..."
        coords_anwil(1,:) = x
        coords_anwil(2,:) = y_orig
        vecs_anwil(1,:)   = x
        vecs_anwil(2,:)   = y_orig
        call anwil_smooth_sigma(coords_anwil, vecs_anwil, smoothed_anwil, 2, 2, n_points, k_neighbors, &
             kd_indices, dimension_order, workspace, value_buffer, permutation, left_stack, right_stack, &
             sigma_raw_anwil, sd_arr_anwil, 0, 1.0_real64, k_neighbors_sigma, kernel_type, ierr)
    end if

    ! --- 5. Nadaraya-Watson ---
    if (method_id == 0 .or. method_id == 5) then
        print *, "Executing Nadaraya-Watson (Adaptive Gaussian)..."
        coords_nw(1,:) = x
        vecs_nw(1,:)   = y_orig
        call smooth_vectors_gaussian_adaptive_nw(coords_nw, vecs_nw, smoothed_nw, 1, 1, n_points, k_neighbors, &
             kd_indices, dimension_order, neighbors, distances, workspace, value_buffer, permutation, &
             left_stack, right_stack, span_loess, ierr, 1)

        print *, "Executing NW KNN Local (k=50)..."
        call smooth_vectors_gaussian_adaptive_nw(coords_nw, vecs_nw, smoothed_nw_knn, 1, 1, n_points, k_neighbors, &
             kd_indices, dimension_order, neighbors, distances, workspace, value_buffer, permutation, &
             left_stack, right_stack, 0.5_real64, ierr, 2)
    end if

    ! --- 6. ManLe (Manifold Learning) ---
    if (method_id == 0 .or. method_id == 6) then
        print *, "Executing ManLe pipeline..."
        data2d(1,:) = x
        data2d(2,:) = y_orig
        call manle_pipeline(data2d, n_points, 2, k_neighbors, n_iters_max, tol_manle, &
             Omega, Y, Q, B, Stmp, Utmp, tau, work_manle, size(work_manle), manifold2d, svd_line, &
             k_neighbors_sigma, kernel_type, ierr)
        x_manle_center(:) = manifold2d(1,:)
        y_manle_center(:) = manifold2d(2,:)
    end if

    ! --- 7. AManLe (Anisotropic Manifold Learning) ---
    if (method_id == 0 .or. method_id == 7) then
        print *, "Executing AManLe pipeline..."
        data2d(1,:) = x
        data2d(2,:) = y_orig
        call amanle_pipeline(data2d, n_points, 2, 1, k_neighbors, n_iters_max, tol_manle, &
             work_amanle, size(work_amanle), manifold2d, svd_line_amanle, k_neighbors_sigma, kernel_type, ierr)
        x_amanle_center(:) = manifold2d(1,:)
        y_amanle_center(:) = manifold2d(2,:)
    end if

    ! --- 8. ANWIL Iterative (Optimization) ---
    if (method_id == 0 .or. method_id == 8) then
        print *, "Executing ANWIL Iterative Optimization..."
        data2d(1,:) = x
        data2d(2,:) = y_orig
        call anwil_iterative(data2d, n_points, 2, k_neighbors, n_iters_max, patience_k, tol_rel, min_iters, &
             k_neighbors_sigma, kernel_type, 0, 1.0_real64, manifold2d, ierr, method_flag, w_r, w_e, w_c, &
             history_scores, history_coverage, history_penalty, history_rmse, history_roughness, &
             stop_iter, stop_reason, best_iter, best_score)
        
        if (ierr /= 0_int32) then
            print *, "ERROR: anwil_iterative failed. error_code=", ierr
        else
            y_anwil_iterative(:) = manifold2d(2,:)
            x_anwil_iterative(:) = manifold2d(1,:)
        end if
    end if

    ! ============================
    ! EXPORT RESULTS (CSV)
    ! ============================
    open(12, file=trim(outfile), status="replace", action="write")
    write(12,'(A)') 'x_original,y_original,x_loess,y_loess,x_anwil,y_anwil,x_anwil_iterative,y_anwil_iterative,x_nw,y_nw,y_nw_knn,x_manle,y_manle,x_manle_svd,y_manle_svd,x_amanle,y_amanle'

    do i = 1, n_points
        write(12,'(F12.6,16(",",F12.6))') x(i), y_orig(i), x_loess(i), y_loess(i), &
            smoothed_anwil(1,i), smoothed_anwil(2,i), x_anwil_iterative(i), y_anwil_iterative(i), &
            coords_nw(1,i), smoothed_nw(1,i), smoothed_nw_knn(1,i), x_manle_center(i), y_manle_center(i), &
            svd_line(1,i), svd_line(2,i), x_amanle_center(i), y_amanle_center(i)
    end do
    close(12)

    print *, "Processing complete. Main results saved to:", trim(outfile)

    ! Export Score History
    call make_output_name(outfile, k_neighbors, n_iters_max, 3, span_loess, k_neighbors_sigma, kernel_type, method_flag, w_r, w_e, w_c)
    call write_anwil_score_history_csv(trim(outfile), history_scores, history_coverage, history_penalty, history_rmse, history_roughness, &
                                       stop_iter, stop_reason, best_iter, best_score, ierr)

    ! Export Standard Deviation (Sigma) Arrays
    call make_output_name(outfile, k_neighbors, n_iters_max, 2, span_loess, k_neighbors_sigma, kernel_type, method_flag, w_r, w_e, w_c)
    open(20, file=outfile, status="replace", action="write")
    write(20, '(A)') 'x,local_sigma_raw,local_sigma_smooth'
    do i = 1, n_points
        write(20, '(F12.6,2(",",F12.6))') x(i), sigma_raw_anwil(1,i), sd_arr_anwil(1,i)
    end do
    close(20)
    print *, "Sigma (Standard Deviation) data saved to:", trim(outfile)

contains

    subroutine make_output_name(name, k_neighbors, n_iters_max, option, span, k_neighbors_sigma, kernel_type, method_flag, w_r, w_e, w_c)
        character(len=*), intent(inout) :: name
        integer(int32), intent(in) :: k_neighbors, n_iters_max
        real(real64), intent(in) :: span, w_r, w_e, w_c
        integer :: dotpos, i, option
        character(len=32) :: k_s, it_s, sp_s, ker_s, sig_s, met_s, wr_s, we_s, wc_s
        integer(int32), intent(in) :: kernel_type, k_neighbors_sigma, method_flag
        integer :: anwil_pos

        ! Find file extension
        dotpos = 0
        do i = len_trim(name), 1, -1
            if (name(i:i) == '.') then
                dotpos = i
                exit
            end if
        end do

        ! String conversions
        write(k_s, '(I0)') k_neighbors
        write(it_s, '(I0)') n_iters_max
        write(sp_s, '(F4.2)') span
        write(ker_s, '(I0)') kernel_type
        write(sig_s, '(I0)') k_neighbors_sigma
        write(met_s, '(I0)') method_flag
        write(wr_s, '(F4.2)') w_r
        write(we_s, '(F4.2)') w_e
        write(wc_s, '(F4.2)') w_c

        if (dotpos > 0) then
            if (option == 1) then 
                name = trim(name(1:dotpos-1)) // '_smoothed_k' // trim(k_s) // '_iter' // trim(it_s) // &
                       '_span' // trim(adjustl(sp_s)) // '_ksigma' // trim(sig_s) // '_kernel' // trim(ker_s) // &
                       '_method' // trim(met_s) // '_wr' // trim(wr_s) // '_we' // trim(we_s) // '_wc' // trim(wc_s) // '.csv'
            else if (option == 2) then
                anwil_pos = index(name, 'anwil')
                if (anwil_pos > 0) then
                    name = trim(name(1:anwil_pos+4)) // '_std.csv'
                else
                    name = trim(name(1:dotpos-1)) // '_anwil_std.csv'
                end if
            else if (option == 3) then
                name = trim(name(1:dotpos-1)) // '_anwil_score_history.csv'
            end if
        end if
    end subroutine make_output_name

    subroutine write_anwil_score_history_csv(filename, history_scores, history_coverage, history_penalty, &
                                            history_rmse, history_roughness, stop_iter, stop_reason, best_iter, best_score, ierr)
        use iso_fortran_env, only: real64, int32
        use tox_errors,      only: set_ok
        implicit none
        character(len=*), intent(in) :: filename
        real(real64),     intent(in) :: history_scores(:), history_roughness(:), history_rmse(:), &
                                        history_coverage(:), history_penalty(:)
        integer(int32),   intent(in) :: stop_iter, stop_reason, best_iter
        real(real64),     intent(in) :: best_score
        integer(int32),   intent(out):: ierr

        integer :: unit, ios
        integer(int32) :: t

        call set_ok(ierr)
        open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) then
            ierr = 1_int32
            return
        end if

        write(unit,'(a)') "# ANWIL Iterative Score History"
        write(unit,'(a,i0)') "# stop_iter=", stop_iter
        write(unit,'(a,i0)') "# stop_reason=", stop_reason
        write(unit,'(a,i0)') "# best_iter=", best_iter
        write(unit,'(a,es24.16)') "# best_score=", best_score
        write(unit,'(a)') "iter,score,roughness,rmse,coverage,coverage_penalty"

        do t = 1_int32, stop_iter
            write(unit,'(i0,5(a,es24.16))') t, ",", history_scores(t), ",", history_roughness(t), ",", &
                  history_rmse(t), ",", history_coverage(t), ",", history_penalty(t)
        end do
        close(unit)
    end subroutine write_anwil_score_history_csv

end program smooth_all_methods_functional