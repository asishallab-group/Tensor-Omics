!> Comparison program for Adaptive KNN smoothing
!| Compares Adaptive KNN smoothing using kallisto_sex_data.tsv dataset
!| Outputs:
!| - adaptive_knn_smoothing_results.tsv: original and smoothed values for each gene using Adaptive KNN

program comparison_adaptive_knn
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use knn_smoothing
    use f42_utils, only: sort_real
    implicit none

    ! Dataset dimensions
    integer(int32), parameter :: MAX_GENES = 90000
    integer(int32), parameter :: N_TISSUES = 67  ! Actual number of tissue columns
    integer(int32), parameter :: MAX_LINE_LENGTH = 10000

    ! Parameters for smoothing
    integer(int32), parameter :: K_MIN = 3
    integer(int32), parameter :: K_MAX = 63
    integer(int32), parameter :: K_INCREMENT = 5
    real(real64), parameter :: SIGMA_FACTOR = 1.5_real64
    real(real64), parameter :: EPSILON = 0.01_real64

    ! Data arrays
    real(real64), allocatable :: input_matrix(:)
    real(real64), allocatable :: x_coords(:,:), y_values(:), w_values(:)
    real(real64), allocatable :: y_smoothed_knn(:,:)
    character(len=50), allocatable :: gene_ids(:)

    ! Analysis arrays
    real(real64), allocatable :: gene_means(:), gene_stds(:)
    real(real64), allocatable :: temp_gene_values(:)

    ! KNN workspace arrays
    integer(int32), allocatable :: kd_indices(:), dimension_order(:)
    integer(int32), allocatable :: neighbors(:), workspace_int(:), permutation(:)
    integer(int32), allocatable :: left_stack(:), right_stack(:)
    real(real64), allocatable :: distances(:), value_buffer(:)

    ! File handling
    integer(int32) :: file_unit, iostat_val, n_genes, ierr
    character(len=MAX_LINE_LENGTH) :: line
    character(len=50) :: temp_string
    integer(int32) :: i, j, pos, start_pos
    real(real64) :: temp_val, mean_val, std_val

    write(*,*) "=== Adaptive KNN Smoothing Comparison ==="
    write(*,*) "Loading data from material/kallisto_sex_data.tsv..."

    ! ========== LOAD DATA ==========
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

    ! ========== ALLOCATE WORKSPACE ARRAYS ==========
    allocate(gene_means(n_genes))
    allocate(gene_stds(n_genes))
    allocate(temp_gene_values(N_TISSUES))
    allocate(x_coords(1, n_genes))
    allocate(y_values(n_genes))
    allocate(w_values(n_genes))
    allocate(y_smoothed_knn(1, n_genes))
    allocate(kd_indices(n_genes))
    allocate(dimension_order(1))
    allocate(neighbors(K_MAX))
    allocate(distances(K_MAX))
    allocate(workspace_int(n_genes))
    allocate(value_buffer(n_genes))
    allocate(permutation(n_genes))
    allocate(left_stack(n_genes))
    allocate(right_stack(n_genes))

    dimension_order(1) = 1

    ! ========== COMPUTE GENE STATISTICS ==========
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
        w_values(i) = 1.0_real64
    end do

    write(*,*) "Gene statistics computed successfully!"

    ! ========== SORT DATA BY MEAN USING F42_UTILS ==========
    write(*,*) "Sorting data by gene mean using f42_utils before smoothing..."

    permutation = [(i, i=1, n_genes)]

    call sort_real(gene_means, permutation, left_stack, right_stack)

    ! Reorder all arrays based on the sorted permutation
    gene_means = gene_means(permutation)
    gene_stds = gene_stds(permutation)
    y_values = y_values(permutation)
    x_coords(1, :) = x_coords(1, permutation)
    gene_ids = gene_ids(permutation)

    deallocate(permutation, left_stack, right_stack)

    write(*,*) "Data sorted successfully using f42_utils!"

    ! ========== PERFORM ADAPTIVE KNN SMOOTHING ==========
    write(*,*) "Performing Adaptive KNN smoothing..."

    ! Reallocate permutation array for use in smoothing
    allocate(permutation(n_genes))
    permutation = [(i, i=1, n_genes)]

    ! Reallocate left_stack and right_stack arrays for use in smoothing
    allocate(left_stack(n_genes), right_stack(n_genes))

    call smooth_vectors_gaussian_adaptive(x_coords, reshape(y_values, [1, n_genes]), y_smoothed_knn, &
                                         1, 1, n_genes, K_MIN, K_MAX, K_INCREMENT, SIGMA_FACTOR, EPSILON, &
                                         kd_indices, dimension_order, neighbors, distances, &
                                         workspace_int, value_buffer, permutation, left_stack, right_stack, ierr, use_global_roughness=.true.)


    if (ierr /= 0) then
        write(*,*) "ERROR in Adaptive KNN smoothing, ierr =", ierr
        stop 1
    end if

    write(*,*) "Adaptive KNN smoothing completed successfully!"

    ! ========== WRITE OUTPUT ==========
    write(*,*) "Writing Adaptive KNN smoothing results to adaptive_knn_smoothing_results.tsv..."

    open(newunit=file_unit, file='adaptive_knn_smoothing_results.tsv', status='replace', action='write')
    write(file_unit, '(A)') 'GeneID' // char(9) // 'Mean' // char(9) // 'OriginalStd' // char(9) // 'SmoothedStd_AdaptiveKNN'

    do i = 1, n_genes
        write(file_unit, '(A,A,F12.6,A,F12.6,A,F12.6)') &
            trim(gene_ids(i)), char(9), &
            gene_means(i), char(9), &
            gene_stds(i), char(9), &
            y_smoothed_knn(1,i)
    end do

    close(file_unit)

    write(*,*) "Adaptive KNN smoothing results written successfully!"

    ! Cleanup
    deallocate(input_matrix, gene_means, gene_stds, temp_gene_values)
    deallocate(gene_ids)
    deallocate(x_coords, y_values, w_values)
    deallocate(y_smoothed_knn)
    deallocate(kd_indices, dimension_order, neighbors, distances)
    deallocate(workspace_int, value_buffer, left_stack, right_stack)

end program comparison_adaptive_knn
