program test_lomanle
    use iso_fortran_env, only: int32, real64
    use lomanle_mod,     only: lomanle_compute
    implicit none

    ! --- Control Variables ---
    character(len=512) :: infile, line
    integer :: ios, n_points, i, arg_len, j
    integer(int32) :: k_min, ierr, dim, manifold_dim, lwork
    real(real64)   :: g_threshold, o_max, o_min

    ! --- Work Buffers ---
    integer(int32), allocatable :: kd_indices(:), workspace(:), perm(:)
    integer(int32), allocatable :: l_stack(:), r_stack(:), rec_stack(:,:), dim_order(:)
    integer(int32), allocatable :: n_loc(:) ! <--- NEW: Passed as argument
    real(real64),   allocatable :: val_buf(:), work_lapack(:)
    real(real64),   allocatable :: d_loc(:) ! <--- NEW: Passed as argument
    
    ! --- Data and Results ---
    real(real64), allocatable :: coords(:,:), radii(:), densities(:)
    real(real64), allocatable :: gap_values(:)
    real(real64), allocatable :: tangent_bases(:,:,:), tangent_scales(:,:)
    logical,      allocatable :: is_anchor(:)

    ! 1. Argument Handling
    call get_command_argument(1, infile, length=arg_len)
    if (arg_len == 0) stop "Usage: ./test_lomanle input.csv k_min [d] [g] [o_max] [o_min]"
    
    call get_command_argument(2, line) ; read(line, *, iostat=ios) k_min
    call get_command_argument(3, line, length=arg_len) ; manifold_dim = 1
    if (arg_len > 0) read(line, *) manifold_dim

    call get_command_argument(4, line, length=arg_len) ; g_threshold = 3.0_real64
    if (arg_len > 0) read(line, *) g_threshold

    call get_command_argument(5, line, length=arg_len) ; o_max = 0.30_real64
    if (arg_len > 0) read(line, *) o_max

    call get_command_argument(6, line, length=arg_len) ; o_min = 0.05_real64
    if (arg_len > 0) read(line, *) o_min

    ! 2. Count Points and Determine Dim
    ! Note: Defaulting to dim=2. Could be detected from first line if needed.
    dim = 2 
    lwork = max(1, 3*dim - 1)
    open(10, file=trim(infile), status="old", action="read")
    n_points = 0 ; read(10, *) 
    do
        read(10, *, iostat=ios) 
        if (ios /= 0) exit
        n_points = n_points + 1
    end do
    rewind(10)

    ! 3. Memory Allocation
    allocate(coords(dim, n_points), radii(n_points), densities(n_points))
    allocate(gap_values(n_points), is_anchor(n_points))
    allocate(tangent_bases(dim, manifold_dim, n_points)) 
    allocate(tangent_scales(manifold_dim, n_points))    
    allocate(kd_indices(n_points), workspace(n_points), perm(n_points))
    allocate(l_stack(n_points), r_stack(n_points), rec_stack(3, n_points))
    allocate(val_buf(n_points), dim_order(dim), work_lapack(lwork))
    
    ! Buffers for KNN queries (size n_points is safe, k_limit would also suffice)
    allocate(n_loc(n_points), d_loc(n_points))

    ! 4. Load Data
    read(10, *) ! Skip header
    do i = 1, n_points
        read(10, *) coords(1, i), coords(2, i)
    end do
    close(10)

    ! 5. Run LoManLe
    print *, "Running LoManLe Multi-D (d=", manifold_dim, ")..."
    call lomanle_compute(coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                        o_max, o_min, kd_indices, workspace, val_buf, perm, &
                        l_stack, r_stack, rec_stack, dim_order, work_lapack, lwork, &
                        n_loc, d_loc, & ! <--- Pass allocated buffers here
                        radii, densities, gap_values, is_anchor, &
                        tangent_bases, tangent_scales, ierr)

    if (ierr /= 0) stop "Error in lomanle_compute"

    ! 6. Save Results (Dynamic Multi-D Writing)
    open(20, file="lomanle_output.csv", status="replace", action="write")
    
    ! --- Write Dynamic Header ---
    write(20, '(A)', advance='no') "x,y,radius,density,gap,anchor"
    do j = 1, manifold_dim
        write(20, '(A,I0,A,I0,A,I0)', advance='no') ",v", j, "_x,v", j, "_y,s", j
    end do
    write(20, *) 

    ! --- Write Dynamic Data ---
    do i = 1, n_points
        ! Base columns
        write(20, '(5(F16.8,","), I1)', advance='no') &
            coords(1,i), coords(2,i), radii(i), densities(i), gap_values(i), merge(1, 0, is_anchor(i))
        
        ! Tangent basis columns (v_j_x, v_j_y, s_j)
        do j = 1, manifold_dim
            write(20, '(A, F16.8, A, F16.8, A, F16.8)', advance='no') &
                ",", tangent_bases(1, j, i), ",", tangent_bases(2, j, i), ",", tangent_scales(j, i)
        end do
        write(20, *) 
    end do

    close(20)
    print *, "Success! CSV saved for manifold_dim =", manifold_dim

    ! 7. Cleanup
    deallocate(kd_indices, workspace, perm, l_stack, r_stack, rec_stack, &
               val_buf, dim_order, work_lapack, coords, radii, densities, &
               gap_values, is_anchor, tangent_bases, tangent_scales, n_loc, d_loc)

end program test_lomanle