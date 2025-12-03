!> @file test_loess_visualization.f90
!> Generates LOESS smoothing data for visualization and comparison

program test_loess_visualization
  use iso_fortran_env, only: real64, int32
  use tox_loess
  implicit none
  
  integer(int32) :: i, j, n, nq, d, liv, lv, ierr
  real(real64), allocatable :: x(:,:), y(:), w(:), xq(:,:), y_out(:)
  integer(int32), allocatable :: iv(:)
  real(real64), allocatable :: v(:)
  character(len=100) :: filename
  
  print *, "=========================================="
  print *, "LOESS SMOOTHING VISUALIZATION TEST"
  print *, "=========================================="
  
  ! Parameters
  d = 1              ! 1D data
  n = 100            ! Number of data points
  nq = 200           ! Number of query points (for smooth curve)
  
  ! Allocate arrays
  allocate(x(d, n), y(n), w(n), xq(d, nq), y_out(nq))
  allocate(iv(1000), v(5000))  ! Conservative workspace sizes
  
  !======================================================================
  ! TEST 1: Multi-scale signal with noise
  !======================================================================
  print *, ""
  print *, "=== Test 1: Multi-scale signal with noise ==="
  
  ! Generate data: slow sine wave + fast sine wave + noise
  do i = 1, n
    x(1, i) = real(i-1, real64) * 0.1_real64  ! x from 0 to 9.9
    ! True signal: 2 sine waves at different frequencies
    y(i) = sin(x(1, i)) + 0.5_real64 * sin(3.0_real64 * x(1, i))
    ! Add some random-ish noise
    y(i) = y(i) + 0.2_real64 * sin(real(i, real64) * 0.7_real64)
  end do
  w = 1.0_real64  ! Uniform weights
  
  ! Generate query points (more densely sampled)
  do i = 1, nq
    xq(1, i) = 0.0_real64 + (9.9_real64 - 0.0_real64) * real(i-1, real64) / real(nq-1, real64)
  end do
  
  ! Save raw data
  open(unit=10, file='raw_data.tsv', status='replace')
  write(10, '(A)') 'x_raw	y_raw'
  do i = 1, n
    write(10, '(F10.6, A, F15.10)') x(1, i), char(9), y(i)
  end do
  close(10)
  print *, "Saved raw data to raw_data.tsv"
  
  !======================================================================
  ! Test different spans and save results
  !======================================================================
  
  ! Span 0.1: Very local smoothing
  call tox_loess_predict(d, n, x, y, w, 0.1_real64, 2, &
                        iv, 1000, v, 5000, nq, xq, y_out, ierr)
  if (ierr /= 0) print *, "Warning: ierr = ", ierr, " for span=0.1"
  
  open(unit=11, file='loess_span_0.1.tsv', status='replace')
  write(11, '(A)') 'x	y_span_0.1'
  do i = 1, nq
    write(11, '(F10.6, A, F15.10)') xq(1, i), char(9), y_out(i)
  end do
  close(11)
  
  ! Span 0.3: Moderate smoothing
  call tox_loess_predict(d, n, x, y, w, 0.3_real64, 2, &
                        iv, 1000, v, 5000, nq, xq, y_out, ierr)
  if (ierr /= 0) print *, "Warning: ierr = ", ierr, " for span=0.3"
  
  open(unit=12, file='loess_span_0.3.tsv', status='replace')
  write(12, '(A)') 'x	y_span_0.3'
  do i = 1, nq
    write(12, '(F10.6, A, F15.10)') xq(1, i), char(9), y_out(i)
  end do
  close(12)
  
  ! Span 0.5: Strong smoothing
  call tox_loess_predict(d, n, x, y, w, 0.5_real64, 2, &
                        iv, 1000, v, 5000, nq, xq, y_out, ierr)
  if (ierr /= 0) print *, "Warning: ierr = ", ierr, " for span=0.5"
  
  open(unit=13, file='loess_span_0.5.tsv', status='replace')
  write(13, '(A)') 'x	y_span_0.5'
  do i = 1, nq
    write(13, '(F10.6, A, F15.10)') xq(1, i), char(9), y_out(i)
  end do
  close(13)
  
  ! Span 0.8: Very strong smoothing (nearly global polynomial)
  call tox_loess_predict(d, n, x, y, w, 0.8_real64, 2, &
                        iv, 1000, v, 5000, nq, xq, y_out, ierr)
  if (ierr /= 0) print *, "Warning: ierr = ", ierr, " for span=0.8"
  
  open(unit=14, file='loess_span_0.8.tsv', status='replace')
  write(14, '(A)') 'x	y_span_0.8'
  do i = 1, nq
    write(14, '(F10.6, A, F15.10)') xq(1, i), char(9), y_out(i)
  end do
  close(14)
  
  print *, "Saved LOESS results for spans: 0.1, 0.3, 0.5, 0.8"
  
  !======================================================================
  ! TEST 2: Sharp discontinuity to test local adaptation
  !======================================================================
  print *, ""
  print *, "=== Test 2: Data with sharp discontinuity ==="
  
  ! Create data with sudden jump
  do i = 1, n
    x(1, i) = real(i-1, real64) * 0.1_real64
    if (x(1, i) < 5.0_real64) then
      y(i) = sin(x(1, i)) + 0.1_real64 * sin(real(i, real64) * 0.5_real64)
    else
      y(i) = 2.0_real64 + cos(x(1, i) - 5.0_real64) + 0.1_real64 * sin(real(i, real64) * 0.5_real64)
    end if
  end do
  
  ! Save raw data with discontinuity
  open(unit=15, file='discontinuity_raw.tsv', status='replace')
  write(15, '(A)') 'x_disc	y_disc'
  do i = 1, n
    write(15, '(F10.6, A, F15.10)') x(1, i), char(9), y(i)
  end do
  close(15)
  
  ! Apply LOESS with moderate span
  call tox_loess_predict(d, n, x, y, w, 0.3_real64, 2, &
                        iv, 1000, v, 5000, nq, xq, y_out, ierr)
  
  open(unit=16, file='discontinuity_loess.tsv', status='replace')
  write(16, '(A)') 'x_disc	y_disc_smoothed'
  do i = 1, nq
    write(16, '(F10.6, A, F15.10)') xq(1, i), char(9), y_out(i)
  end do
  close(16)
  
  print *, "Saved discontinuity test data"
  
  !======================================================================
  ! TEST 3: Varying noise level
  !======================================================================
  print *, ""
  print *, "=== Test 3: Heteroscedastic noise (varying amplitude) ==="
  
  ! Generate data with increasing noise amplitude
  do i = 1, n
    x(1, i) = real(i-1, real64) * 0.1_real64
    ! True signal
    y(i) = exp(-0.2_real64 * x(1, i)) * sin(2.0_real64 * x(1, i))
    ! Noise amplitude increases with x
    y(i) = y(i) + 0.1_real64 * (1.0_real64 + x(1, i)/5.0_real64) * &
           sin(real(i, real64) * 1.3_real64)
  end do
  
  ! Save raw data
  open(unit=17, file='varying_noise_raw.tsv', status='replace')
  write(17, '(A)') 'x_var	y_var'
  do i = 1, n
    write(17, '(F10.6, A, F15.10)') x(1, i), char(9), y(i)
  end do
  close(17)
  
  ! Apply LOESS with different spans
  call tox_loess_predict(d, n, x, y, w, 0.2_real64, 2, &
                        iv, 1000, v, 5000, nq, xq, y_out, ierr)
  
  open(unit=18, file='varying_noise_span_0.2.tsv', status='replace')
  write(18, '(A)') 'x_var	y_var_span_0.2'
  do i = 1, nq
    write(18, '(F10.6, A, F15.10)') xq(1, i), char(9), y_out(i)
  end do
  close(18)
  
  call tox_loess_predict(d, n, x, y, w, 0.4_real64, 2, &
                        iv, 1000, v, 5000, nq, xq, y_out, ierr)
  
  open(unit=19, file='varying_noise_span_0.4.tsv', status='replace')
  write(19, '(A)') 'x_var	y_var_span_0.4'
  do i = 1, nq
    write(19, '(F10.6, A, F15.10)') xq(1, i), char(9), y_out(i)
  end do
  close(19)
  
  print *, "Saved varying noise test data"
  
  !======================================================================
  ! Create R script for visualization
  !======================================================================
  print *, ""
  print *, "=== Creating R visualization script ==="
  
  open(unit=20, file='plot_loess_results.R', status='replace')
  write(20, '(A)') '# R script to visualize LOESS smoothing results'
  write(20, '(A)') '# Generated by Fortran test program'
  write(20, '(A)') ''
  write(20, '(A)') '# Load libraries'
  write(20, '(A)') 'library(ggplot2)'
  write(20, '(A)') 'library(dplyr)'
  write(20, '(A)') 'library(tidyr)'
  write(20, '(A)') ''
  write(20, '(A)') '# Read data'
  write(20, '(A)') 'raw_data <- read.table("raw_data.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'loess_01 <- read.table("loess_span_0.1.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'loess_03 <- read.table("loess_span_0.3.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'loess_05 <- read.table("loess_span_0.5.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'loess_08 <- read.table("loess_span_0.8.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') ''
  write(20, '(A)') '# Combine all LOESS results'
  write(20, '(A)') 'all_loess <- loess_01 %>%'
  write(20, '(A)') '  left_join(loess_03, by="x") %>%'
  write(20, '(A)') '  left_join(loess_05, by="x") %>%'
  write(20, '(A)') '  left_join(loess_08, by="x") %>%'
  write(20, '(A)') '  pivot_longer(cols = -x, names_to = "method", values_to = "y")'
  write(20, '(A)') ''
  write(20, '(A)') '# Plot 1: Multi-scale signal with different spans'
  write(20, '(A)') 'p1 <- ggplot() +'
  write(20, '(A)') '  geom_point(data = raw_data, aes(x = x_raw, y = y_raw),'
  write(20, '(A)') '             alpha = 0.5, size = 1.5, color = "gray60") +'
  write(20, '(A)') '  geom_line(data = all_loess, aes(x = x, y = y, color = method),'
  write(20, '(A)') '            linewidth = 1.2) +'
  write(20, '(A)') '  scale_color_manual('
  write(20, '(A)') '    name = "LOESS span",'
  write(20, '(A)') '    values = c("y_span_0.1" = "#E41A1C",'
  write(20, '(A)') '               "y_span_0.3" = "#377EB8",'
  write(20, '(A)') '               "y_span_0.5" = "#4DAF4A",'
  write(20, '(A)') '               "y_span_0.8" = "#984EA3"),'
  write(20, '(A)') '    labels = c("0.1", "0.3", "0.5", "0.8")) +'
  write(20, '(A)') '  labs(title = "LOESS Smoothing with Different Spans",'
  write(20, '(A)') '       subtitle = "Small span = local fitting, Large span = global smoothing",'
  write(20, '(A)') '       x = "X", y = "Y") +'
  write(20, '(A)') '  theme_minimal() +'
  write(20, '(A)') '  theme(legend.position = "bottom")'
  write(20, '(A)') ''
  write(20, '(A)') '# Plot 2: Discontinuity handling'
  write(20, '(A)') 'disc_raw <- read.table("discontinuity_raw.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'disc_loess <- read.table("discontinuity_loess.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') ''
  write(20, '(A)') 'p2 <- ggplot() +'
  write(20, '(A)') '  geom_point(data = disc_raw, aes(x = x_disc, y = y_disc),'
  write(20, '(A)') '             alpha = 0.5, size = 1.5, color = "gray60") +'
  write(20, '(A)') '  geom_line(data = disc_loess, aes(x = x_disc, y = y_disc_smoothed),'
  write(20, '(A)') '            color = "#D95F02", linewidth = 1.5) +'
  write(20, '(A)') '  geom_vline(xintercept = 5.0, linetype = "dashed", color = "red", alpha = 0.5) +'
  write(20, '(A)') '  labs(title = "LOESS Handling of Discontinuity (span=0.3)",'
  write(20, '(A)') '       subtitle = "Vertical line shows actual discontinuity location",'
  write(20, '(A)') '       x = "X", y = "Y") +'
  write(20, '(A)') '  theme_minimal()'
  write(20, '(A)') ''
  write(20, '(A)') '# Plot 3: Varying noise'
  write(20, '(A)') 'var_raw <- read.table("varying_noise_raw.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'var_02 <- read.table("varying_noise_span_0.2.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') 'var_04 <- read.table("varying_noise_span_0.4.tsv", header=TRUE, sep="\t")'
  write(20, '(A)') ''
  write(20, '(A)') 'p3 <- ggplot() +'
  write(20, '(A)') '  geom_point(data = var_raw, aes(x = x_var, y = y_var),'
  write(20, '(A)') '             alpha = 0.5, size = 1.5, color = "gray60") +'
  write(20, '(A)') '  geom_line(data = var_02, aes(x = x_var, y = y_var_span_0.2),'
  write(20, '(A)') '            color = "#1B9E77", linewidth = 1.2) +'
  write(20, '(A)') '  geom_line(data = var_04, aes(x = x_var, y = y_var_span_0.4),'
  write(20, '(A)') '            color = "#7570B3", linewidth = 1.2) +'
  write(20, '(A)') '  labs(title = "LOESS with Heteroscedastic Noise",'
  write(20, '(A)') '       subtitle = "Noise amplitude increases with X",'
  write(20, '(A)') '       x = "X", y = "Y") +'
  write(20, '(A)') '  theme_minimal()'
  write(20, '(A)') ''
  write(20, '(A)') '# Save plots'
  write(20, '(A)') 'pdf("loess_visualization.pdf", width = 10, height = 12)'
  write(20, '(A)') 'print(p1)'
  write(20, '(A)') 'print(p2)'
  write(20, '(A)') 'print(p3)'
  write(20, '(A)') 'dev.off()'
  write(20, '(A)') ''
  write(20, '(A)') '# Also save individual PNGs'
  write(20, '(A)') 'ggsave("loess_spans_comparison.png", p1, width = 8, height = 6, dpi = 300)'
  write(20, '(A)') 'ggsave("loess_discontinuity.png", p2, width = 8, height = 6, dpi = 300)'
  write(20, '(A)') 'ggsave("loess_varying_noise.png", p3, width = 8, height = 6, dpi = 300)'
  write(20, '(A)') ''
  write(20, '(A)') 'cat("Plots saved as:")'
  write(20, '(A)') 'cat("\n• loess_visualization.pdf (all plots)")'
  write(20, '(A)') 'cat("\n• loess_spans_comparison.png")'
  write(20, '(A)') 'cat("\n• loess_discontinuity.png")'
  write(20, '(A)') 'cat("\n• loess_varying_noise.png")'
  write(20, '(A)') 'cat("\n\nData files available for inspection.")'
  
  close(20)
  
  print *, "Created R script: plot_loess_results.R"
  print *, ""
  print *, "To visualize results:"
  print *, "1. Run this Fortran program"
  print *, "2. Run in R: source('plot_loess_results.R')"
  print *, "3. Check generated PDF and PNG files"
  
  ! Cleanup
  deallocate(x, y, w, xq, y_out, iv, v)
  
  print *, ""
  print *, "=========================================="
  print *, "TEST COMPLETE - Check generated .tsv files"
  print *, "=========================================="
  
end program test_loess_visualization