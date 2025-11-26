!> Benchmark program for KNN smoothing with different k values
!| Uses the same data loading and statistics as comparison_knn_loess.f90
!| Only runs KNN smoothing for a list of k values and outputs results for each

program benchmark_knn_smoothing
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use knn_smoothing
    implicit none

    integer(int32), parameter :: MAX_GENES = 90000
    integer(int32), parameter :: N_TISSUES = 67
    integer(int32), parameter :: MAX_LINE_LENGTH = 10000
    integer(int32), parameter :: N_K = 7
    integer(int32), parameter :: K_LIST(N_K) = (/ 10, 20, 30, 50, 80, 110, 150 /)

    ! Data arrays
    real(real64), allocatable :: input_matrix(:)
    real(real64), allocatable :: x_coords(:,:), y_values(:)
    character(len=50), allocatable :: gene_ids(:)
    real(real64), allocatable :: gene_means(:), gene_stds(:)
    real(real64), allocatable :: temp_gene_values(:)
    real(real64), allocatable :: y_smoothed_knn(:,:)

    ! KNN workspace arrays (max k)
    integer(int32), allocatable :: kd_indices(:), dimension_order(:)
    integer(int32), allocatable :: neighbors(:), workspace_int(:), permutation(:)
    integer(int32), allocatable :: left_stack(:), right_stack(:)
    real(real64), allocatable :: distances(:), value_buffer(:)

    integer(int32) :: file_unit, iostat_val, n_genes, ierr
    character(len=MAX_LINE_LENGTH) :: line
    character(len=50) :: temp_string
    integer(int32) :: i, j, pos, start_pos, k_idx, k_val
    real(real64) :: temp_val, mean_val, std_val
    integer :: clk_start, clk_end, clk_rate
    real(real64) :: elapsed_sec

    write(*,*) "=== KNN Smoothing Benchmark ==="
    write(*,*) "Loading data from material/kallisto_sex_data.tsv..."

    ! Initialize clock rate for timing
    call system_clock(count_rate=clk_rate)

    allocate(gene_ids(MAX_GENES))
    allocate(input_matrix(MAX_GENES * N_TISSUES))

    open(newunit=file_unit, file='material/kallisto_sex_data.tsv', status='old', action='read', iostat=iostat_val)
    if (iostat_val /= 0) then
        write(*,*) "ERROR: Could not open material/kallisto_sex_data.tsv"
        stop 1
    end if
    read(file_unit, '(A)', iostat=iostat_val) line
    if (iostat_val /= 0) then
        write(*,*) "ERROR: Could not read header line"
        stop 1
    end if
    n_genes = 0
    do while (.true.)
        read(file_unit, '(A)', iostat=iostat_val) line
        if (iostat_val /= 0) exit
        n_genes = n_genes + 1
        if (n_genes > MAX_GENES) then
            write(*,*) "ERROR: Too many genes, increase MAX_GENES"
            stop 1
        end if
        pos = 1
        start_pos = pos
        do while (pos <= len_trim(line) .and. line(pos:pos) /= char(9))
            pos = pos + 1
        end do
        gene_ids(n_genes) = line(start_pos:pos-1)
        pos = pos + 1
        do j = 1, N_TISSUES
            start_pos = pos
            do while (pos <= len_trim(line) .and. line(pos:pos) /= char(9))
                pos = pos + 1
            end do
            if (start_pos <= len_trim(line)) then
                temp_string = line(start_pos:min(pos-1, len_trim(line)))
            else
                temp_string = ''
            end if
            if (len_trim(temp_string) == 0) then
                temp_val = 0.0_real64
            else
                read(temp_string, *, iostat=iostat_val) temp_val
                if (iostat_val /= 0) temp_val = 0.0_real64
            end if
            input_matrix((n_genes-1)*N_TISSUES + j) = temp_val
            if (pos <= len_trim(line) .and. line(pos:pos) == char(9)) then
                pos = pos + 1
            end if
        end do
        if (mod(n_genes, 10000) == 0) then
            write(*,*) "Loaded", n_genes, "genes..."
        end if
    end do
    close(file_unit)
    write(*,*) "Data loaded successfully:", n_genes, "genes,", N_TISSUES, "tissues"

    allocate(gene_means(n_genes))
    allocate(gene_stds(n_genes))
    allocate(temp_gene_values(N_TISSUES))
    allocate(x_coords(1, n_genes))
    allocate(y_values(n_genes))
    allocate(y_smoothed_knn(N_K, n_genes))
    allocate(kd_indices(n_genes))
    allocate(dimension_order(1))
    allocate(neighbors(maxval(K_LIST)))
    allocate(distances(maxval(K_LIST)))
    allocate(workspace_int(n_genes))
    allocate(value_buffer(n_genes))
    allocate(permutation(n_genes))
    allocate(left_stack(n_genes))
    allocate(right_stack(n_genes))
    dimension_order(1) = 1

    write(*,*) "Computing mean and standard deviation for each gene..."
    do i = 1, n_genes
        mean_val = 0.0_real64
        do j = 1, N_TISSUES
            mean_val = mean_val + input_matrix((i-1)*N_TISSUES + j)
        end do
        mean_val = mean_val / real(N_TISSUES, real64)
        gene_means(i) = mean_val
        std_val = 0.0_real64
        do j = 1, N_TISSUES
            std_val = std_val + (input_matrix((i-1)*N_TISSUES + j) - mean_val)**2
        end do
        std_val = sqrt(std_val / real(N_TISSUES - 1, real64))
        gene_stds(i) = std_val
        x_coords(1, i) = gene_means(i)
        y_values(i) = gene_stds(i)
    end do
    write(*,*) "Gene statistics computed successfully!"
    write(*,*) "Mean range: [", minval(gene_means), ", ", maxval(gene_means), "]"
    write(*,*) "Std range: [", minval(gene_stds), ", ", maxval(gene_stds), "]"

    do k_idx = 1, N_K
        k_val = K_LIST(k_idx)
        write(*,*) "\nPerforming KNN smoothing for k=", k_val, "..."
        call system_clock(clk_start)
        call smooth_vectors_gaussian_adaptive(x_coords, reshape(y_values, [1, n_genes]), y_smoothed_knn(k_idx,:), &
                                             1, 1, n_genes, k_val, &
                                             kd_indices, dimension_order, neighbors, distances, &
                                             workspace_int, value_buffer, permutation, &
                                             left_stack, right_stack, ierr)
        call system_clock(clk_end)
        elapsed_sec = real(clk_end - clk_start, real64) / real(clk_rate, real64)
        if (ierr /= 0) then
            write(*,*) "ERROR in KNN smoothing for k=", k_val, ", ierr =", ierr
            stop 1
        end if
        write(*,*) "  KNN smoothing for k=", k_val, "completed. Time: ", elapsed_sec, "sec. Range: [", minval(y_smoothed_knn(k_idx,:)), ", ", maxval(y_smoothed_knn(k_idx,:)), "]"
        ! Output results for this k
        write(*,*) "  Writing results to knn_smoothing_results_k"//trim(adjustl(itoa(k_val)))//".tsv ..."
        open(newunit=file_unit, file='knn_smoothing_results_k'//trim(adjustl(itoa(k_val)))//'.tsv', status='replace', action='write')
        write(file_unit, '(A)') 'GeneID' // char(9) // 'Mean' // char(9) // 'OriginalStd' // char(9) // 'SmoothedStd_KNN_k'
        do i = 1, n_genes
            write(file_unit, '(A,A,F12.6,A,F12.6,A,F12.6)') &
                trim(gene_ids(i)), char(9), &
                gene_means(i), char(9), &
                gene_stds(i), char(9), &
                y_smoothed_knn(k_idx,i)
        end do
        close(file_unit)
    end do

    write(*,*) "\n=== KNN SMOOTHING BENCHMARK SUMMARY ==="
    write(*,*) "Dataset:", n_genes, "genes x", N_TISSUES, "tissues"
    write(*,*) "Smoothing input: X = gene means, Y = gene standard deviations"
    write(*,*) "K tested:", K_LIST
    write(*,*) "Input data ranges:"
    write(*,*) "  - Gene means: [", minval(gene_means), ", ", maxval(gene_means), "]"
    write(*,*) "  - Gene std devs: [", minval(gene_stds), ", ", maxval(gene_stds), "]"
    do k_idx = 1, N_K
        write(*,*) "  - KNN smoothed std (k=", K_LIST(k_idx), "): [", minval(y_smoothed_knn(k_idx,:)), ", ", maxval(y_smoothed_knn(k_idx,:)), "]"
    end do
    write(*,*) "\nOutput files:"
    do k_idx = 1, N_K
        write(*,*) "  - knn_smoothing_results_k", K_LIST(k_idx), ".tsv"
    end do
    write(*,*) "\nKNN smoothing benchmark complete!"

    ! Cleanup
    deallocate(input_matrix, gene_means, gene_stds, temp_gene_values)
    deallocate(gene_ids)
    deallocate(x_coords, y_values)
    deallocate(y_smoothed_knn)
    deallocate(kd_indices, dimension_order, neighbors, distances)
    deallocate(workspace_int, value_buffer, permutation, left_stack, right_stack)

contains
    ! Simple integer to string conversion (Fortran 2003+)
    pure function itoa(i) result(str)
        integer, intent(in) :: i
        character(len=12) :: str
        write(str, '(I0)') i
    end function itoa

end program benchmark_knn_smoothing
