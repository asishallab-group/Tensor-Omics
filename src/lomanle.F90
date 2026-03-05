module lomanle_mod
    use iso_fortran_env, only: int32, real64
    use kd_tree,         only: build_kd_index, kd_knn_query
    use tox_errors,      only: set_ok, is_ok
    use f42_utils,       only: sort_array
    implicit none

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: int32, real64
            character,        intent(in)    :: jobz, uplo
            integer(int32),   intent(in)    :: n, lda, lwork
            real(real64),     intent(inout) :: a(lda, n)
            real(real64),     intent(out)   :: w(n)
            real(real64),     intent(out)   :: work(lwork)
            integer(int32),   intent(out)   :: info
        end subroutine dsyev
    end interface

contains

    subroutine lomanle_compute(coords, n_points, dim, manifold_dim, k_min, g_threshold, &
                               o_max, o_min, &
                               kd_indices, workspace, val_buf, perm, &
                               l_stack, r_stack, rec_stack, &
                               dim_order, work_lapack, lwork, & 
                               n_loc, d_loc, & ! <--- NEW EXTERNAL BUFFERS
                               sphere_radii, densities, gap_values, &
                               is_anchor, tangent_bases, tangent_scales, ierr)
        
        integer(int32), intent(in) :: n_points, dim, manifold_dim, k_min, lwork
        real(real64),   intent(in) :: coords(dim, n_points)
        real(real64),   intent(in) :: g_threshold, o_max, o_min
        
        integer(int32), intent(inout) :: kd_indices(n_points), workspace(n_points), perm(n_points)
        integer(int32), intent(inout) :: l_stack(n_points), r_stack(n_points), rec_stack(3, n_points)
        integer(int32), intent(inout) :: dim_order(dim)
        real(real64),   intent(inout) :: val_buf(n_points), work_lapack(lwork)

        ! EXTERNAL BUFFERS (allocated in caller with size n_points or k_limit)
        integer(int32), intent(inout) :: n_loc(:) 
        real(real64),   intent(inout) :: d_loc(:)
        
        real(real64),   intent(out) :: sphere_radii(n_points)
        real(real64),   intent(out) :: densities(n_points)
        real(real64),   intent(out) :: gap_values(n_points)   
        logical,        intent(out) :: is_anchor(n_points)
        real(real64),   intent(out) :: tangent_bases(dim, manifold_dim, n_points)
        real(real64),   intent(out) :: tangent_scales(manifold_dim, n_points)
        integer(int32), intent(out) :: ierr

        ! Local variables
        integer(int32) :: i, j, k, row, col, info, k_curr, k_limit, d_idx, tmp_idx, base_idx
        integer(int32) :: n_anchor, i_seed, k_idx, m_idx, num_overlap, total_in_sphere
        logical        :: is_covered(n_points)
        real(real64)   :: current_ratio, dist_sq, dist_centers, sigma_i, sigma2_i, rho_i, s_gap
        real(real64)   :: center(dim), cov(dim, dim), w_eig(dim), p_diff(dim)

        call set_ok(ierr)
        k_limit = n_points / 4 
        is_anchor = .false.
        is_covered = .false.
        tangent_bases = 0.0_real64
        tangent_scales = 0.0_real64
        do i = 1, dim ; dim_order(i) = i ; end do

        ! STEP 0: Build KD-Tree
        call build_kd_index(coords, dim, n_points, kd_indices, dim_order, &
                            workspace, val_buf, perm, l_stack, r_stack, rec_stack, ierr)
        if (.not. is_ok(ierr)) return

        ! --- STEPS 1, 2 and 3: Adaptive Radius and SVD ---
        do i = 1, n_points
            k_curr = k_min
            adaptive_k: do
                ! Query KNN using provided external buffers
                call kd_knn_query(coords, kd_indices, dim, n_points, perm, &
                                  coords(:,i), k_curr, n_loc(1:k_curr), d_loc(1:k_curr), ierr)
                
                sigma_i = maxval(d_loc(1:k_curr))
                sphere_radii(i) = sigma_i
                center = sum(coords(:, n_loc(1:k_curr)), dim=2) / real(k_curr, real64)

                cov = 0.0_real64
                do j = 1, k_curr
                    p_diff = coords(:, n_loc(j)) - center
                    do col = 1, dim
                        do row = 1, dim
                            cov(row, col) = cov(row, col) + p_diff(row) * p_diff(col)
                        end do
                    end do
                end do
                cov = cov / max(1.0_real64, real(k_curr - 1, real64))

                call dsyev('V', 'U', dim, cov, dim, w_eig, work_lapack, lwork, info)
                if (info /= 0) exit adaptive_k

                d_idx = dim - manifold_dim + 1
                if (w_eig(d_idx - 1) > 1.0e-12_real64) then
                    s_gap = sqrt(w_eig(d_idx)) / sqrt(w_eig(d_idx - 1))
                else
                    s_gap = g_threshold + 1.0_real64
                end if
                gap_values(i) = s_gap

                ! Density (rho)
                rho_i = 0.0_real64
                sigma2_i = sigma_i**2
                if (sigma_i > 1.0e-12_real64) then
                    do j = 1, k_curr
                        rho_i = rho_i + exp(-(d_loc(j)**2) / (2.0_real64 * sigma2_i))
                    end do
                else
                    rho_i = real(k_curr, real64)
                end if
                densities(i) = rho_i / (max(1.0e-12_real64, sphere_radii(i))**manifold_dim)

                ! EXIT CONDITION: If spectral gap is sufficient, capture SVD basis
                if (s_gap >= g_threshold .or. k_curr >= k_limit) then
                    do base_idx = 1, manifold_dim
                        tangent_bases(:, base_idx, i) = cov(:, dim - manifold_dim + base_idx)
                        tangent_scales(base_idx, i) = sqrt(max(0.0_real64, w_eig(dim - manifold_dim + base_idx)))
                    end do
                    exit adaptive_k
                end if
                
                k_curr = nint(k_curr * 1.25)
                if (k_curr > size(n_loc)) k_curr = size(n_loc) ! Guard against buffer overflow
            end do adaptive_k
        end do

        ! --- STEP 4: Sort by Density ---
        do i = 1, n_points ; perm(i) = i ; end do
        call sort_array(densities, perm, l_stack, r_stack)
        do i = 1, n_points / 2 ! Reverse to descending
            tmp_idx = perm(i) ; perm(i) = perm(n_points - i + 1) ; perm(n_points - i + 1) = tmp_idx
        end do

        ! --- STEP 5: Atlas Selection (Expansion) ---
        i_seed = perm(1)
        is_anchor(i_seed) = .true.
        n_anchor = 1
        do m_idx = 1, n_points
            if (sum((coords(:, i_seed) - coords(:, m_idx))**2) <= sphere_radii(i_seed)**2) is_covered(m_idx) = .true.
        end do

        do k_idx = 2, n_points
            i = perm(k_idx)
            if (is_anchor(i)) cycle
            num_overlap = 0 ; total_in_sphere = 0
            do j = 1, n_points
                dist_sq = sum((coords(:, i) - coords(:, j))**2)
                if (dist_sq <= sphere_radii(i)**2) then
                    total_in_sphere = total_in_sphere + 1
                    if (is_covered(j)) num_overlap = num_overlap + 1
                end if
            end do
            current_ratio = 0.0_real64
            if (total_in_sphere > 0) current_ratio = real(num_overlap, real64) / real(total_in_sphere, real64)

            if (current_ratio >= o_min .and. current_ratio <= o_max) then
                is_anchor(i) = .true. ; n_anchor = n_anchor + 1
                do m_idx = 1, n_points
                    if (sum((coords(:, i) - coords(:, m_idx))**2) <= sphere_radii(i)**2) is_covered(m_idx) = .true.
                end do
            end if
        end do

        ! Orphans pass
        do k_idx = 1, n_points
            i = perm(k_idx)
            if (is_covered(i) .or. is_anchor(i)) cycle
            is_anchor(i) = .true. ; n_anchor = n_anchor + 1
            do m_idx = 1, n_points
                if (sum((coords(:, i) - coords(:, m_idx))**2) <= sphere_radii(i)**2) is_covered(m_idx) = .true.
            end do
        end do

        ! STEP 5.1: Redundancy Cleanup
        do i = 1, n_points
            if (.not. is_anchor(i)) cycle
            do j = 1, n_points
                if (.not. is_anchor(j) .or. i == j) cycle
                dist_centers = sqrt(sum((coords(:, i) - coords(:, j))**2))
                if (dist_centers + sphere_radii(i) <= sphere_radii(j)) then
                    is_anchor(i) = .false. ; n_anchor = n_anchor - 1
                    exit
                end if
            end do
        end do

    end subroutine lomanle_compute

end module lomanle_mod