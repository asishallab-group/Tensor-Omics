program test_lomanle
    use iso_fortran_env, only: int32, real64
    use lomanle_mod,     only: lomanle_compute_alloc
    implicit none

    character(len=512) :: infile, line
    integer :: ios, n_points, i, arg_len, j
    integer(int32) :: k_min, ierr, dim, manifold_dim, max_iterations
    real(real64)   :: g_threshold, o_max, o_min, stability_threshold, scale_factor, relative_conv_tol

    real(real64), allocatable :: coords(:,:), radii(:), densities(:)
    real(real64), allocatable :: gap_values(:)
    real(real64), allocatable :: normal_errors(:), stability_values(:)
    integer(int32), allocatable :: k_selected(:)
    logical,      allocatable :: growth_stopped_complex(:)
    real(real64), allocatable :: skeleton_coords(:,:)
    real(real64), allocatable :: skeleton_iter1(:,:)
    real(real64), allocatable :: radii_1(:), densities_1(:), gap_1(:)
    logical,      allocatable :: is_anchor_1(:)
    integer(int32), allocatable :: labels_1(:)
    real(real64), allocatable :: tangent_bases_1(:,:,:), tangent_scales_1(:,:)
    integer(int32), allocatable :: primary_anchor_iter1(:), secondary_anchor_iter1(:)
    integer(int32), allocatable :: primary_anchor_final(:), secondary_anchor_final(:)

    ! Topological skeleton (backbone) graph, built inside lomanle_compute from the
    ! anchor adjacency graph's minimum spanning tree -- see build_skeleton_edges in
    ! lomanle.F90. Edges are point-index pairs (edge_from(k), edge_to(k)); anchor_role
    ! classifies each anchor as endpoint (1), pass-through (2) or branch/junction (3).
    real(real64),   allocatable :: anchor_centers_1(:,:), anchor_centers_final(:,:)
    integer(int32), allocatable :: anchor_role_1(:), anchor_role_final(:)
    integer(int32), allocatable :: edge_from_1(:), edge_to_1(:)
    integer(int32), allocatable :: edge_from_final(:), edge_to_final(:)
    integer(int32) :: n_edges_1, n_edges_final, e_idx

    call get_command_argument(1, infile, length=arg_len)
    if (arg_len == 0) stop "Usage: ./test_lomanle input.csv k_min [d] [g] [o_max] [o_min] [stability] [scale_factor] [max_iterations] [relative_conv_tol] [dim]"

    call get_command_argument(2, line) ; read(line, *, iostat=ios) k_min

    call get_command_argument(3, line, length=arg_len) ; manifold_dim = 1
    if (arg_len > 0) read(line, *) manifold_dim

    call get_command_argument(4, line, length=arg_len) ; g_threshold = 3.0_real64
    if (arg_len > 0) read(line, *) g_threshold

    call get_command_argument(5, line, length=arg_len) ; o_max = 0.30_real64
    if (arg_len > 0) read(line, *) o_max

    call get_command_argument(6, line, length=arg_len) ; o_min = 0.05_real64
    if (arg_len > 0) read(line, *) o_min

    call get_command_argument(7, line, length=arg_len) ; stability_threshold = 0.90_real64
    if (arg_len > 0) read(line, *) stability_threshold

    call get_command_argument(8, line, length=arg_len) ; scale_factor = 2.5_real64
    if (arg_len > 0) read(line, *) scale_factor

    call get_command_argument(9, line, length=arg_len) ; max_iterations = 50
    if (arg_len > 0) read(line, *) max_iterations

    call get_command_argument(10, line, length=arg_len) ; relative_conv_tol = 0.01_real64
    if (arg_len > 0) read(line, *) relative_conv_tol

    ! Ambient dimension of the input coordinates (2 for datasets under
    ! results/data/2d/, 3 for results/data/3d/). Defaults to 2 since that's
    ! the more common case; previously this was a hand-edited literal
    ! (dim = 3) that had to be changed and the binary rebuilt every time to
    ! switch between 2D and 3D datasets.
    call get_command_argument(11, line, length=arg_len) ; dim = 2
    if (arg_len > 0) read(line, *) dim

    open(10, file=trim(infile), status="old", action="read")
    n_points = 0 ; read(10, *)
    do
        read(10, *, iostat=ios)
        if (ios /= 0) exit
        n_points = n_points + 1
    end do
    rewind(10)

    allocate(coords(dim, n_points), radii(n_points), densities(n_points))
    allocate(gap_values(n_points))
    allocate(normal_errors(n_points), stability_values(n_points))
    allocate(k_selected(n_points), growth_stopped_complex(n_points))
    allocate(skeleton_coords(dim+1, n_points))
    allocate(skeleton_iter1(dim+1, n_points))
    allocate(radii_1(n_points), densities_1(n_points), gap_1(n_points))
    allocate(is_anchor_1(n_points), labels_1(n_points))
    allocate(tangent_bases_1(dim, manifold_dim, n_points))
    allocate(tangent_scales_1(manifold_dim, n_points))
    allocate(primary_anchor_iter1(n_points), secondary_anchor_iter1(n_points))
    allocate(primary_anchor_final(n_points), secondary_anchor_final(n_points))
    allocate(anchor_centers_1(dim, n_points), anchor_centers_final(dim, n_points))
    allocate(anchor_role_1(n_points), anchor_role_final(n_points))

    read(10, *)
    do i = 1, n_points
        read(10, *) coords(1:dim, i)
    end do
    close(10)

    print *, "Running LoManLe Multi-D (d=", manifold_dim, ")..."
    call lomanle_compute_alloc(coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                        o_max, o_min, stability_threshold, scale_factor, &
                        max_iterations, relative_conv_tol, &
                        radii, densities, gap_values, &
                        normal_errors, stability_values, k_selected, growth_stopped_complex, &
                        skeleton_coords, skeleton_iter1, &
                        radii_1, densities_1, gap_1, is_anchor_1, labels_1, &
                        tangent_bases_1, tangent_scales_1, &
                        primary_anchor_final, secondary_anchor_final, &
                        primary_anchor_iter1, secondary_anchor_iter1, &
                        anchor_centers_1, anchor_centers_final, &
                        edge_from_1, edge_to_1, n_edges_1, anchor_role_1, &
                        edge_from_final, edge_to_final, n_edges_final, anchor_role_final, &
                        ierr)
    if (ierr /= 0) stop "Error in lomanle_compute_alloc"

    open(20, file="lomanle_output.csv", status="replace", action="write")
    if (dim == 3) then
        write(20, '(A)', advance='no') &
            "x,y,z,n_anchors,sk_x,sk_y,sk_z,sk_x_f,sk_y_f,sk_z_f,radius,density,gap,anchor,label"
    else
        write(20, '(A)', advance='no') "x,y,n_anchors,sk_x,sk_y,sk_x_f,sk_y_f,radius,density,gap,anchor,label"
    end if
    do j = 1, manifold_dim
        if (dim == 3) then
            write(20, '(A,I0,A,I0,A,I0,A,I0)', advance='no') ",v", j, "_x,v", j, "_y,v", j, "_z,s", j
        else
            write(20, '(A,I0,A,I0,A,I0)', advance='no') ",v", j, "_x,v", j, "_y,s", j
        end if
    end do
    write(20, '(A)') ",k_selected,stability,normal_error,stopped_complex,anchor_role_iter1,anchor_role_final"

    do i = 1, n_points
        if (dim == 3) then
            write(20, '(10(F16.8,","), 3(F16.8,","), I1, ",", I0)', advance='no') &
                coords(1,i), coords(2,i), coords(3,i), &
                skeleton_iter1(1,i), skeleton_iter1(2,i), skeleton_iter1(3,i), skeleton_iter1(4,i), &
                skeleton_coords(2,i), skeleton_coords(3,i), skeleton_coords(4,i), &
                radii_1(i), densities_1(i), gap_1(i), &
                merge(1, 0, is_anchor_1(i)), labels_1(i)
        else
            write(20, '(7(F16.8,","), 3(F16.8,","), I1, ",", I0)', advance='no') &
                coords(1,i), coords(2,i), skeleton_iter1(1,i), skeleton_iter1(2,i), skeleton_iter1(3,i), &
                skeleton_coords(2,i), skeleton_coords(3,i), &
                radii_1(i), densities_1(i), gap_1(i), &
                merge(1, 0, is_anchor_1(i)), labels_1(i)
        end if

        do j = 1, manifold_dim
            if (dim == 3) then
                write(20, '(A, F16.8, A, F16.8, A, F16.8, A, F16.8)', advance='no') &
                    ",", tangent_bases_1(1, j, i), ",", tangent_bases_1(2, j, i), &
                    ",", tangent_bases_1(3, j, i), ",", tangent_scales_1(j, i)
            else
                write(20, '(A, F16.8, A, F16.8, A, F16.8)', advance='no') &
                    ",", tangent_bases_1(1, j, i), ",", tangent_bases_1(2, j, i), ",", tangent_scales_1(j, i)
            end if
        end do

        ! Final-iteration neighbor-selection diagnostics (see smoothing_vecinos.md)
        write(20, '(A, I0, A, F16.8, A, F16.8, A, I0)', advance='no') &
            ",", k_selected(i), ",", stability_values(i), ",", normal_errors(i), &
            ",", merge(1, 0, growth_stopped_complex(i))

        ! Backbone anchor role, from the MST over the anchor graph: 0 = not an
        ! anchor, 1 = endpoint, 2 = pass-through, 3 = branch/junction.
        write(20, '(A, I0, A, I0)') ",", anchor_role_1(i), ",", anchor_role_final(i)
    end do
    close(20)

    ! Skeleton edges: the actual backbone line, as point-to-point segments, so
    ! plotting code only needs to draw geom_segment(x,y,xend,yend) -- no joining
    ! logic lives outside lomanle anymore.
    open(21, file="lomanle_edges.csv", status="replace", action="write")
    if (dim == 3) then
        write(21, '(A)') "stage,edge_id,x,y,z,xend,yend,zend"
        do e_idx = 1, n_edges_1
            write(21, '(A, ",", I0, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8)') &
                "iter1", e_idx, &
                skeleton_iter1(2, edge_from_1(e_idx)), skeleton_iter1(3, edge_from_1(e_idx)), &
                skeleton_iter1(4, edge_from_1(e_idx)), &
                skeleton_iter1(2, edge_to_1(e_idx)),   skeleton_iter1(3, edge_to_1(e_idx)), &
                skeleton_iter1(4, edge_to_1(e_idx))
        end do
        do e_idx = 1, n_edges_final
            write(21, '(A, ",", I0, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8)') &
                "final", e_idx, &
                skeleton_coords(2, edge_from_final(e_idx)), skeleton_coords(3, edge_from_final(e_idx)), &
                skeleton_coords(4, edge_from_final(e_idx)), &
                skeleton_coords(2, edge_to_final(e_idx)),   skeleton_coords(3, edge_to_final(e_idx)), &
                skeleton_coords(4, edge_to_final(e_idx))
        end do
    else
        write(21, '(A)') "stage,edge_id,x,y,xend,yend"
        do e_idx = 1, n_edges_1
            write(21, '(A, ",", I0, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8)') &
                "iter1", e_idx, &
                skeleton_iter1(2, edge_from_1(e_idx)), skeleton_iter1(3, edge_from_1(e_idx)), &
                skeleton_iter1(2, edge_to_1(e_idx)),   skeleton_iter1(3, edge_to_1(e_idx))
        end do
        do e_idx = 1, n_edges_final
            write(21, '(A, ",", I0, ",", F16.8, ",", F16.8, ",", F16.8, ",", F16.8)') &
                "final", e_idx, &
                skeleton_coords(2, edge_from_final(e_idx)), skeleton_coords(3, edge_from_final(e_idx)), &
                skeleton_coords(2, edge_to_final(e_idx)),   skeleton_coords(3, edge_to_final(e_idx))
        end do
    end if
    close(21)

    print *, "Success! CSV saved for manifold_dim =", manifold_dim

    deallocate(coords, radii, densities, &
               gap_values, normal_errors, stability_values, k_selected, growth_stopped_complex, &
               skeleton_coords, skeleton_iter1, radii_1, densities_1, gap_1, is_anchor_1, labels_1, &
               tangent_bases_1, tangent_scales_1, &
               primary_anchor_iter1, secondary_anchor_iter1, primary_anchor_final, secondary_anchor_final, &
               anchor_centers_1, anchor_centers_final, anchor_role_1, anchor_role_final, &
               edge_from_1, edge_to_1, edge_from_final, edge_to_final)

end program test_lomanle
