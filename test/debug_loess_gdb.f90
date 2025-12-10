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
  character(len=256) :: input_files(n_files) = ['material/kallisto_sex_data_no_na.tsv']
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
  n_genes = 88327
  
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
  
  ! Generate dummy gene IDs (replace with actual gene IDs from your files)
  do i = 1, n_genes
    write(gene_ids(i), '(A,I6.6)') 'GENE_', i
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
      write(10, '(F15.10, A, F15.10, A, A)') &
        gene_means(i), char(9), gene_stds(i), char(9), trim(gene_ids(i))
    end if
  end do
  close(10)
  print *, "Saved raw mean-std data to gene_mean_std_raw.tsv"
  
  ! Save LOESS results
  open(unit=11, file='gene_mean_std_loess.tsv', status='replace')
  write(11, '(A)') 'x_query	y_span_0.2	y_span_0.4	y_span_0.6	y_span_0.8'
  do i = 1, nq
    write(11, '(F15.10, 4(A, F15.10))') &
      xq(1, i), &
      (char(9), loess_results(i, j), j = 1, n_spans)
  end do
  close(11)
  print *, "Saved LOESS results to gene_mean_std_loess.tsv"
  
  ! ======================================================================
  ! 8. CREATE R VISUALIZATION SCRIPT
  ! ======================================================================
  
  print *, ""
  print *, "=== Creating R visualization script ==="
  
  open(unit=20, file='plot_gene_mean_std.R', status='replace')
  write(20, '(A)') '# R script to visualize gene expression mean-std relationship'
  write(20, '(A)') '# with LOESS smoothing'
  write(20, '(A)') ''
  write(20, '(A)') '# Load libraries'
  write(20, '(A)') 'library(ggplot2)'
  write(20, '(A)') 'library(dplyr)'
  write(20, '(A)') 'library(tidyr)'
  write(20, '(A)') ''
  write(20, '(A)') '# Read data'
  write(20, '(A)') 'raw_data <- read.table("gene_mean_std_raw.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'loess_data <- read.table("gene_mean_std_loess.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') ''
  write(20, '(A)') '# Basic statistics'
  write(20, '(A)') 'cat("Summary of gene expression means:\n")'
  write(20, '(A)') 'print(summary(raw_data$mean))'
  write(20, '(A)') 'cat("\nSummary of gene expression standard deviations:\n")'
  write(20, '(A)') 'print(summary(raw_data$std))'
  write(20, '(A)') 'cat("\nNumber of genes analyzed:", nrow(raw_data), "\n")'
  write(20, '(A)') ''
  write(20, '(A)') '# Plot 1: Mean vs STD with LOESS smoothing'
  write(20, '(A)') 'p1 <- ggplot(raw_data, aes(x = mean, y = std)) +'
  write(20, '(A)') '  geom_point(alpha = 0.3, size = 0.8, color = "gray50") +'
  write(20, '(A)') '  geom_line(data = loess_data, aes(x = x_query, y = y_span_0.4),'
  write(20, '(A)') '            color = "#E41A1C", linewidth = 1.2, linetype = "solid") +'
  write(20, '(A)') '  geom_line(data = loess_data, aes(x = x_query, y = y_span_0.6),'
  write(20, '(A)') '            color = "#377EB8", linewidth = 1.2, linetype = "dashed") +'
  write(20, '(A)') '  labs(title = "Gene Expression: Mean vs Standard Deviation",'
  write(20, '(A)') '       subtitle = "LOESS smoothing reveals heteroscedasticity pattern",'
  write(20, '(A)') '       x = "Mean Expression Level",'
  write(20, '(A)') '       y = "Standard Deviation",'
  write(20, '(A)') '       caption = paste("n =", nrow(raw_data), "genes")) +'
  write(20, '(A)') '  theme_minimal() +'
  write(20, '(A)') '  theme(legend.position = "none")'
  write(20, '(A)') ''
  write(20, '(A)') '# Plot 2: All LOESS spans comparison'
  write(20, '(A)') 'loess_long <- loess_data %>%'
  write(20, '(A)') '  pivot_longer(cols = -x_query, names_to = "span", values_to = "std_pred")'
  write(20, '(A)') ''
  write(20, '(A)') 'p2 <- ggplot() +'
  write(20, '(A)') '  geom_point(data = raw_data, aes(x = mean, y = std),'
  write(20, '(A)') '             alpha = 0.2, size = 0.5, color = "gray70") +'
  write(20, '(A)') '  geom_line(data = loess_long, aes(x = x_query, y = std_pred, color = span),'
  write(20, '(A)') '            linewidth = 1.2) +'
  write(20, '(A)') '  scale_color_manual('
  write(20, '(A)') '    name = "LOESS span",'
  write(20, '(A)') '    values = c("y_span_0.2" = "#4DAF4A",'
  write(20, '(A)') '               "y_span_0.4" = "#E41A1C",'
  write(20, '(A)') '               "y_span_0.6" = "#377EB8",'
  write(20, '(A)') '               "y_span_0.8" = "#984EA3"),'
  write(20, '(A)') '    labels = c("0.2 (local)", "0.4", "0.6", "0.8 (global)")) +'
  write(20, '(A)') '  labs(title = "LOESS Smoothing with Different Spans",'
  write(20, '(A)') '       subtitle = "Gene expression mean-std relationship",'
  write(20, '(A)') '       x = "Mean Expression Level", y = "Standard Deviation") +'
  write(20, '(A)') '  theme_minimal() +'
  write(20, '(A)') '  theme(legend.position = "bottom")'
  write(20, '(A)') ''
  write(20, '(A)') '# Save plots'
  write(20, '(A)') 'pdf("gene_mean_std_analysis.pdf", width = 10, height = 8)'
  write(20, '(A)') 'print(p1)'
  write(20, '(A)') 'print(p2)'
  write(20, '(A)') 'dev.off()'
  write(20, '(A)') ''
  write(20, '(A)') 'ggsave("gene_mean_std_main.png", p1, width = 8, height = 6, dpi = 300)'
  write(20, '(A)') 'ggsave("gene_mean_std_spans.png", p2, width = 8, height = 6, dpi = 300)'
  write(20, '(A)') ''
  write(20, '(A)') '# Calculate correlation'
  write(20, '(A)') 'correlation <- cor(raw_data$mean, raw_data$std, use = "complete.obs")'
  write(20, '(A)') 'cat("\nCorrelation between mean and std:", round(correlation, 3), "\n")'
  write(20, '(A)') ''
  write(20, '(A)') '# Fit linear model for comparison'
  write(20, '(A)') 'lm_fit <- lm(std ~ mean, data = raw_data)'
  write(20, '(A)') 'cat("\nLinear model summary:\n")'
  write(20, '(A)') 'print(summary(lm_fit))'
  write(20, '(A)') ''
  write(20, '(A)') 'cat("\n=== Analysis Complete ===\n")'
  write(20, '(A)') 'cat("Files generated:\n")'
  write(20, '(A)') 'cat("• gene_mean_std_raw.tsv - Raw mean/std values\n")'
  write(20, '(A)') 'cat("• gene_mean_std_loess.tsv - LOESS smoothed curves\n")'
  write(20, '(A)') 'cat("• gene_mean_std_analysis.pdf - PDF with all plots\n")'
  write(20, '(A)') 'cat("• gene_mean_std_main.png - Main plot (PNG)\n")'
  write(20, '(A)') 'cat("• gene_mean_std_spans.png - Span comparison (PNG)\n")'
  
  close(20)
  
  print *, "Created R script: plot_gene_mean_std.R"
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