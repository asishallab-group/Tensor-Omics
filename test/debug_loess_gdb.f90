!> @file test_loess_visualization.f90
!> Generates LOESS smoothing data for visualization and comparison

program test_loess_gene_expression
  use iso_fortran_env, only: real64, int32
  use tox_loess
  use tox_data_tools
  implicit none
  
  integer(int32) :: i, j, n_genes, n_samples, d, liv, lv, ierr, n, nq
  real(real64), allocatable :: x(:,:), y(:), w(:), xq(:,:), y_out(:)
  real(real64), allocatable :: expression_vectors(:,:)
  real(real64), allocatable :: gene_means(:), gene_stds(:)
  integer(int32), allocatable :: iv(:)
  real(real64), allocatable :: v(:)
  character(len=100) :: filename
  character(len=256), allocatable :: file_list(:)
  character(len=100), allocatable :: gene_ids(:)
  integer(int32), allocatable :: value_cols(:)
    ! Count valid genes (non-zero in at least some samples)
  integer(int32) :: valid_genes, n_valid_samples
  real(real64) :: sum_val, sum_sq, mean_val, std_val
    ! List of TSV files to read
  integer, parameter :: n_files = 1
  character(len=256) :: input_files(n_files) = ['material/kallisto_test.tsv']
  real(real64) :: min_mean, max_mean
  ! File format parameters
  integer, parameter :: n_header_rows = 1      ! Number of header rows to skip
  integer, parameter :: gene_col = 1           ! Column containing gene IDs
  integer, parameter :: value_start_col = 2    ! First column with expression values
    ! Apply LOESS with different spans
  integer, parameter :: n_spans = 4
  real(real64) :: spans(n_spans) = [0.2_real64, 0.4_real64, 0.6_real64, 0.8_real64]
  real(real64), allocatable :: loess_results(:,:)
  
  print *, "=========================================="
  print *, "GENE EXPRESSION MEAN-STD LOESS ANALYSIS"
  print *, "=========================================="
  
  ! ======================================================================
  ! 1. USER CONFIGURATION - EDIT THESE FOR YOUR DATA
  ! ======================================================================
  
  ! ======================================================================
  ! 2. ALLOCATE MEMORY
  ! ======================================================================
  
  ! Determine number of value columns (assume same for all files)
  ! You may need to adjust this based on your file structure
  
  ! Total samples = n_files * n_value_cols_per_file
  n_samples = 67
  n_genes = 2000
  
  allocate(gene_ids(n_genes))
  allocate(expression_vectors(n_samples, n_genes))
  expression_vectors = 0.0_real64  ! Initialize
  
  ! ======================================================================
  ! 3. READ DATA FROM TSV FILES
  ! ======================================================================
  
  print *, ""
  print *, "=== Reading expression data from TSV files ==="
  
  ! Create value columns array (all columns from value_start_col onward)
  allocate(value_cols(67))
  do i = 1, 67
    value_cols(i) = value_start_col + i - 1
  end do

  call read_gene_ids_from_tsv_file(input_files(1), gene_ids, n_header_rows, gene_col, ierr)
  
  ! Read the files
  call read_expression_vectors_tsv(input_files, gene_ids, expression_vectors, &
                                  n_header_rows, gene_col, value_cols, 1, ierr)
  
  if (ierr /= 0) then
    print *, "Error reading TSV files: ierr = ", ierr
    stop
  end if
  
  print *, "Successfully read ", n_samples, " samples × ", n_genes, " genes"
  
  ! ======================================================================
  ! 4. COMPUTE MEAN AND STANDARD DEVIATION PER GENE
  ! ======================================================================
  
  print *, ""
  print *, "=== Computing mean and standard deviation per gene ==="
  
  allocate(gene_means(n_genes), gene_stds(n_genes))
  
  valid_genes = 0
  do i = 1, n_genes
    sum_val = 0.0_real64
    sum_sq = 0.0_real64
    n_valid_samples = 0
    
    do j = 1, n_samples
      if (abs(expression_vectors(j, i)) > 1.0e-10_real64) then
        sum_val = sum_val + expression_vectors(j, i)
        sum_sq = sum_sq + expression_vectors(j, i)**2
        n_valid_samples = n_valid_samples + 1
      end if
    end do
    
    if (n_valid_samples >= 3) then  ! Need at least 3 samples for meaningful std
      mean_val = sum_val / real(n_valid_samples, real64)
      gene_means(i) = mean_val
      
      if (n_valid_samples > 1) then
        std_val = sqrt((sum_sq - n_valid_samples * mean_val**2) / &
                      real(n_valid_samples - 1, real64))
        gene_stds(i) = std_val
        valid_genes = valid_genes + 1
      else
        gene_stds(i) = 0.0_real64
      end if
    else
      gene_means(i) = 0.0_real64
      gene_stds(i) = 0.0_real64
    end if
  end do
  
  print *, "Computed statistics for ", valid_genes, " valid genes"
  
  ! ======================================================================
  ! 5. PREPARE DATA FOR LOESS (mean vs std)
  ! ======================================================================
  
  print *, ""
  print *, "=== Preparing data for LOESS smoothing (mean vs std) ==="
  
  ! Use only valid genes for LOESS
  d = 1
  n = valid_genes
  
  allocate(x(d, n), y(n), w(n))
  
  ! Fill x with gene means, y with gene stds
  j = 0
  do i = 1, n_genes
    if (gene_stds(i) > 0.0_real64 .and. gene_means(i) /= 0.0_real64) then
      j = j + 1
      x(1, j) = gene_means(i)
      y(j) = gene_stds(i)
      w(j) = 1.0_real64  ! Uniform weights
      
      if (j == n) exit
    end if
  end do
  
  ! ======================================================================
  ! 6. APPLY LOESS SMOOTHING
  ! ======================================================================
  
  print *, ""
  print *, "=== Applying LOESS smoothing ==="
  
  ! Create query points (span the range of mean values)
  nq = 200
  allocate(xq(d, nq), y_out(nq))
  
  
  min_mean = minval(x(1, 1:n))
  max_mean = maxval(x(1, 1:n))
  
  do i = 1, nq
    xq(1, i) = min_mean + (max_mean - min_mean) * real(i-1, real64) / real(nq-1, real64)
  end do
  
  ! Allocate workspace
  call tox_loess_required_workspace(d, n, liv, lv)
  allocate(iv(liv), v(lv))
  

  
  allocate(loess_results(nq, n_spans))
  
  do i = 1, n_spans
    print *, "  Applying LOESS with span = ", spans(i)
    call tox_loess_predict(d, n, x, y, w, spans(i), 2, &
                          iv, liv, v, lv, nq, xq, y_out, ierr)
    
    if (ierr /= 0) then
      print *, "    Warning: ierr = ", ierr, " for span = ", spans(i)
    end if
    
    loess_results(:, i) = y_out
  end do
  
  ! ======================================================================
  ! 7. SAVE RESULTS
  ! ======================================================================
  
  print *, ""
  print *, "=== Saving results to files ==="
  
  ! Save raw mean-std data
  open(unit=10, file='gene_mean_std_raw.tsv', status='replace')
  write(10, '(A)') 'mean	std	gene_id'
  do i = 1, n_genes
    if (gene_stds(i) > 0.0_real64 .and. gene_means(i) /= 0.0_real64) then
      write(10, '(F15.5, A, F20.5, A, A)') &
        gene_means(i), char(9), gene_stds(i), char(9), trim(gene_ids(i))
    end if
  end do
  close(10)
  print *, "Saved raw mean-std data to gene_mean_std_raw.tsv"
  
  ! Save LOESS results
  open(unit=11, file='gene_mean_std_loess.tsv', status='replace')
  write(11, '(A)') 'x_query	y_span_0.2	y_span_0.4	y_span_0.6	y_span_0.8'
  do i = 1, nq
    write(11, '(F20.5, 4(A, F20.5))') &
      xq(1, i), &
      (char(9), loess_results(i, j), j = 1, n_spans)
  end do
  close(11)
  print *, "Saved LOESS results to gene_mean_std_loess.tsv"
  
  ! ======================================================================
  ! 8. CREATE R VISUALIZATION SCRIPT
  ! ======================================================================
  print *, ""
  print *, "To analyze results:"
  print *, "1. Run this Fortran program with your TSV files"
  print *, "2. Run in R: source('plot_gene_mean_std.R')"
  print *, "3. Check generated PDF/PNG files and console output"
  
  ! ======================================================================
  ! 9. CLEANUP
  ! ======================================================================
  
  deallocate(x, y, w, xq, y_out, iv, v)
  deallocate(expression_vectors, gene_means, gene_stds)
  deallocate(gene_ids, value_cols, loess_results)
  
  print *, ""
  print *, "=========================================="
  print *, "ANALYSIS COMPLETE"
  print *, "=========================================="
  
end program test_loess_gene_expression