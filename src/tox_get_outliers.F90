!> Module to identify gene outliers based on their distances to family centroids.
module tox_get_outliers
  use safeguard
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use f42_utils, only: sort_array  
  use tox_errors, only: ERR_INVALID_INPUT, set_ok, set_err_once, check_alloc_stat, is_err
  use tox_loess, only: tox_loess_required_workspace, loess_fit_robust, loess_fit_plain
  implicit none

contains

  !> Compute family scaling factors (dscale) to normalize distances.
  !| Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
  subroutine compute_family_scaling( &
    n_genes, n_families, distances, gene_to_fam, dscale, &
    loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, &
    iv, liv, wv, lv, diagl, w_init, z_mat, rw, ww, res, pi, yhat_tmp, &
    span, degree, mode, n_iters, ierr)

      use, intrinsic :: iso_fortran_env, only: real64, int32
      implicit none

      !| Total number of genes
      integer(int32), intent(in) :: n_genes
      !| Total number of gene families
      integer(int32), intent(in) :: n_families
      !| Array of Euclidean distances for each gene
      real(real64),   intent(in) :: distances(n_genes)
      !| Mapping of each gene to its family (1-based)
      integer(int32), intent(in) :: gene_to_fam(n_genes)

      !| Output: array of scaling factors per family
      real(real64),   intent(out) :: dscale(n_families)

      ! Buffers (reused)
      !| Reference x-coordinates for LOESS smoothing
      real(real64),   intent(inout) :: loess_x(n_families)
      !| Reference y-coordinates for LOESS smoothing
      real(real64),   intent(inout) :: loess_y(n_families)
      !| Indices of reference points used for smoothing
      integer(int32), intent(inout) :: indices_used(n_families)
      !| Permutation array for sorting gene distances
      integer(int32), intent(inout) :: perm_tmp(n_genes)
      !| Stack array for left indices during sorting
      integer(int32), intent(inout) :: stack_left_tmp(n_genes)
      !| Stack array for right indices during sorting
      integer(int32), intent(inout) :: stack_right_tmp(n_genes)
      !| Pre-allocated work array for family distances (dimension n_genes)
      real(real64),   intent(out)   :: family_distances(n_genes)

      ! LOESS workspace
      !| Length of integer workspace
      integer(int32), intent(in)    :: liv, lv
      !| Integer workspace array
      integer(int32), intent(inout) :: iv(liv)
      !| Real workspace array
      real(real64),   intent(inout) :: wv(lv)

      !| Diagonal elements of the weight matrix
      real(real64), intent(inout) :: diagl(:)
      !| Initial weights for LOESS
      real(real64), intent(inout) :: w_init(:)
      !| Z matrix for LOESS fitting
      real(real64), intent(inout) :: z_mat(:,:)
      !| Residuals for robust LOESS fitting
      real(real64), intent(inout) :: rw(:), ww(:), res(:)
      !| Permutation indices for robust LOESS fitting
      integer(int32), intent(inout):: pi(:)
      !| Output array for LOESS predictions
      real(real64), intent(out)    :: yhat_tmp(:)

      !| Span parameter for LOESS smoothing
      real(real64), intent(in)     :: span
      !| Degree of the LOESS polynomial
      integer(int32), intent(in)   :: degree
      !| Mode for LOESS fitting (0=plain, 1=robust)
      integer(int32), intent(in)   :: mode
      !| Number of iterations for robust LOESS fitting
      integer(int32), intent(in)   :: n_iters
      !| Error code: 0=ok, non-zero indicates an error
      integer(int32), intent(out)  :: ierr

      ! Local variables
      integer(int32) :: i, j, family_idx, n_in_family, n_valid
      real(real64)   :: median_dist, stddev_dist, mean_dist, sumsq
      real(real64) :: xmin, xmax
      

      ! Initialize error code and output arrays
      call set_ok(ierr)
      dscale  = 0.0_real64
      loess_x = 0.0_real64
      loess_y = 0.0_real64
      n_valid = 0

      ! Validate family indices
      do i = 1, n_genes
          if (gene_to_fam(i) < 1 .or. gene_to_fam(i) > n_families) then
              dscale = -1.0_real64
              call set_err_once(ierr, ERR_INVALID_INPUT)
              return
          end if
      end do

      ! ------------------------------------------------------------
      ! PASS 1: compute (median, stddev) per family; store COMPACT
      ! ------------------------------------------------------------
      do family_idx = 1, n_families

          n_in_family = 0
          do i = 1, n_genes
              if (gene_to_fam(i) == family_idx) then
                  n_in_family = n_in_family + 1
                  family_distances(n_in_family) = abs(distances(i))
              end if
          end do

          if (n_in_family <= 1) cycle

          do j = 1, n_in_family
              perm_tmp(j) = j
              stack_left_tmp(j) = 0
              stack_right_tmp(j) = 0
          end do

          call sort_array( &
              family_distances(1:n_in_family), perm_tmp(1:n_in_family), &
              stack_left_tmp(1:n_in_family), stack_right_tmp(1:n_in_family))

          if (mod(n_in_family, 2) == 0) then
              median_dist = 0.5_real64 * ( &
                  family_distances(perm_tmp(n_in_family/2)) + &
                  family_distances(perm_tmp(n_in_family/2 + 1)) )
          else
              median_dist = family_distances(perm_tmp((n_in_family + 1)/2))
          end if

          mean_dist = sum(family_distances(1:n_in_family)) / real(n_in_family, real64)
          sumsq     = sum((family_distances(1:n_in_family) - mean_dist)**2)
          stddev_dist = sqrt(sumsq / real(n_in_family - 1, real64))

          n_valid = n_valid + 1
          loess_x(n_valid) = median_dist
          loess_y(n_valid) = stddev_dist
          indices_used(n_valid) = family_idx
      end do

      if (n_valid <= 1) then
          ! not enough points; leave dscale zeros
          return
      end if

      xmin = minval(loess_x(1:n_valid))
      xmax = maxval(loess_x(1:n_valid))
      
      if (xmax == xmin) then
        ! Degenerate: all x are identical => LOESS undefined.
        ! Fallback: just map and back directly 
        do j = 1, n_valid
          family_idx = indices_used(j)
          dscale(family_idx) = loess_y(j)
        end do
        return
      end if


      ! ------------------------------------------------------------
      ! LOESS GLOBAL: smooth y_ref as function of x_ref (once)
      ! ------------------------------------------------------------
      w_init(1:n_valid) = 1.0_real64
      z_mat(1:n_valid,1) = loess_x(1:n_valid)

      if (mode == 0) then
          ! If you have a plain routine, call it; otherwise keep robust always.
          call loess_fit_plain( &
              n_valid, loess_x(1:n_valid), loess_y(1:n_valid), w_init(1:n_valid), z_mat(1:n_valid,1:1), &
              span, degree, n_valid, .false., .false., iv, liv, wv, lv, diagl(1:n_valid), yhat_tmp(1:n_valid), ierr)
      else
          call loess_fit_robust( &
              n_valid, loess_x(1:n_valid), loess_y(1:n_valid), w_init(1:n_valid), z_mat(1:n_valid,1:1), &
              span, degree, n_valid, .false., .false., n_iters, iv, liv, wv, lv, diagl(1:n_valid), &
              rw(1:n_valid), ww(1:n_valid), res(1:n_valid), pi(1:n_valid), yhat_tmp(1:n_valid), ierr)
      end if

      if (is_err(ierr)) then
          dscale = -1.0_real64
          return
      end if

      ! ------------------------------------------------------------
      ! Map compact results back to full dscale
      ! ------------------------------------------------------------
      do j = 1, n_valid
          family_idx = indices_used(j)
          dscale(family_idx) = yhat_tmp(j)
      end do

  end subroutine compute_family_scaling


  !> Helper routine that allocates internal arrays and calls compute_family_scaling.
  !| This makes usage easier since users don't need to care about internal array requirements.
  subroutine compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
    loess_x, loess_y, indices_used, ierr)
    !| Total number of genes
    integer(int32), intent(in) :: n_genes
    !| Total number of gene families
    integer(int32), intent(in) :: n_families
    !| Array of Euclidean distances for each gene
    real(real64), intent(in) :: distances(n_genes)
    !| Mapping of each gene to its family (1-based)
    integer(int32), intent(in) :: gene_to_fam(n_genes)
    !| Output: array of scaling factors per family
    real(real64), intent(out) :: dscale(n_families)
    !| Reference x-coordinates.
    real(real64), intent(inout) :: loess_x(n_families)
    !| Reference y-coordinates (length n_total).
    real(real64), intent(inout) :: loess_y(n_families)
    !| Indices of reference points used for smoothing.
    integer(int32), intent(inout) :: indices_used(n_families)
    !| Error code: 0=ok, 201=invalid family indices
    integer(int32), intent(out) :: ierr

    ! Local work arrays
    real(real64) :: family_distances(n_genes)
    integer(int32) :: perm_tmp(n_genes)
    integer(int32) :: stack_left_tmp(n_genes)
    integer(int32) :: stack_right_tmp(n_genes)

    ! LOESS workspace
    integer(int32) :: liv, lv, istat
    integer(int32), allocatable :: iv(:), pi(:)
    real(real64), allocatable :: wv(:), diagl(:), w_init(:), z_mat(:,:), rw(:), ww(:), res(:), yhat_tmp(:)

    ! LOESS params (defaults)
    real(real64), parameter :: span = 0.7_real64
    integer(int32), parameter :: degree = 2_int32
    integer(int32), parameter :: mode = 1_int32
    integer(int32), parameter :: n_iters = 3_int32
    logical, parameter :: setlf = .false.

    call set_ok(ierr)
    call set_ok(istat)

    ! Workspace sizes 
    call tox_loess_required_workspace(1_int32, n_families, liv, lv, setlf)

    allocate(iv(liv), wv(lv), stat=istat)
    call check_alloc_stat(istat, ierr)
    if(is_err(ierr)) return

    ! For robust we also need arrays sized to n_valid (<= n_families).
    allocate(diagl(n_families), w_init(n_families), z_mat(n_families,1), &
             rw(n_families), ww(n_families), res(n_families), pi(n_families), &
             yhat_tmp(n_families), stat=istat)
    call check_alloc_stat(istat, ierr)
    if(is_err(ierr)) return

    ! Initialize (important for netlib)
    iv = 1_int32
    wv = 0.0_real64
    rw = 1.0_real64
    pi = 0_int32

    call compute_family_scaling( &
        n_genes, n_families, distances, gene_to_fam, dscale, &
        loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, &
        iv, liv, wv, lv, diagl, w_init, z_mat, rw, ww, res, pi, yhat_tmp, &
        span, degree, mode, n_iters, ierr)

    deallocate(iv, wv, diagl, w_init, z_mat, rw, ww, res, pi, yhat_tmp)

  end subroutine compute_family_scaling_alloc

  !> Compute the hybrid RDI for each gene.
  !| RDI = Euclidean distance / family scaling factor
  pure subroutine compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, &
                                    stack_left, stack_right)
    implicit none
    !| Total number of genes
    integer(int32), intent(in) :: n_genes
    !| Array of Euclidean distances for each gene to its centroid
    real(real64), intent(in) :: distances(n_genes)
    !| Gene-to-family mapping (1-based indexing)
    integer(int32), intent(in) :: gene_to_fam(n_genes)
    !| Array of scaling factors for each family
    real(real64), intent(in) :: dscale(:)
    !| Output array of RDI values for each gene
    real(real64), intent(out) :: rdi(n_genes)
    !| Work array for sorting (dimension n_genes)
    real(real64), intent(inout) :: sorted_rdi(n_genes)
    !| Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
    integer(int32), intent(inout) :: perm(n_genes)
    !| Stack array for sorting (dimension n_genes)
    integer(int32), intent(inout) :: stack_left(n_genes)
    !| Stack array for sorting (dimension n_genes)
    integer(int32), intent(inout) :: stack_right(n_genes)
    
    integer(int32) :: i, family_idx
    real(real64), parameter :: tol = epsilon(1.0_real64)
    
    ! Calculate RDI for each gene
    do i = 1, n_genes
      family_idx = gene_to_fam(i)
      
      ! Handle invalid family indices
      if (family_idx < 1 .or. family_idx > size(dscale)) then
        rdi(i) = -1.0_real64  ! Error indicator
        cycle
      end if
      
      ! Detect NaN input (portable)
      if (ieee_is_nan(distances(i))) then
        rdi(i) = distances(i)
      else if (abs(dscale(family_idx)) < tol) then
        rdi(i) = 0.0_real64  ! If scaling is zero, set RDI to zero (not outlier)
      else
        ! Calculate RDI
        rdi(i) = abs(distances(i)) / dscale(family_idx)
      end if
    end do

    ! Create a copy of RDI for sorting (excluding error values)
    sorted_rdi = rdi

    ! Filter out error values (negative RDIs)
    where (sorted_rdi < 0.0_real64)
      sorted_rdi = 0.0_real64
    end where

    ! Sort RDI values using the tox_sorting module
    call sort_array(sorted_rdi, perm, stack_left, stack_right)
    
  end subroutine compute_rdi

  !> Identify gene outliers based on the top percentile of RDI values.
  !| Expects sorted_rdi to be filtered (no negative values) and sorted in ascending order before calling.
  !| If sorted_rdi contains negatives or is not sorted, results may be invalid.
  pure subroutine identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, percentile)
    implicit none

    !| Total number of genes
    integer(int32), intent(in) :: n_genes
    !| Array of RDI values for each gene
    real(real64), intent(in) :: rdi(n_genes)
    !| Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
    real(real64), intent(in) :: sorted_rdi(n_genes)
    !| Output boolean array indicating outliers
    logical, intent(out) :: is_outlier(n_genes)
    !| Output threshold value used for detection
    real(real64), intent(out) :: threshold
    !| (optional) Percentile threshold (default: 95 for top 5%)
    real(real64), intent(in), optional :: percentile
       
    integer(int32) :: i, idx
    real(real64) :: perc_pos, percentile_val

    ! Set default percentile if not present
    if (present(percentile)) then
      percentile_val = percentile
    else
      percentile_val = 95.0_real64
    end if

    ! Initialize output
    is_outlier = .false.

    ! Calculate the position corresponding to the desired percentile
    perc_pos = (n_genes * percentile_val) / 100.0_real64
    idx = ceiling(perc_pos)
    ! Clamp idx to valid range
    if (idx < 1) idx = 1
    if (idx > n_genes) idx = n_genes

    ! Get the threshold value from the sorted array (sorted_rdi must be ascending)
    threshold = sorted_rdi(idx)

    ! Mark genes as outliers if their RDI exceeds the threshold (and is positive)
    do i = 1, n_genes
      is_outlier(i) = (rdi(i) >= threshold .and. rdi(i) > 0.0_real64)
    end do
  end subroutine identify_outliers

  !> Main routine to detect outliers using RDI and LOESS-based scaling.
  subroutine detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                            work_array, perm, stack_left, stack_right, &
                            is_outlier, loess_x, loess_y, loess_n, ierr, &
                            percentile)
    implicit none

    !| Total number of genes
    integer(int32), intent(in) :: n_genes
    !| Total number of gene families
    integer(int32), intent(in) :: n_families
    !| Array of Euclidean distances for each gene to its centroid
    real(real64), intent(in) :: distances(n_genes)
    !| Gene-to-family mapping (1-based indexing)
    integer(int32), intent(in) :: gene_to_fam(n_genes)
    !| Work array for sorting (dimension n_genes)
    real(real64), intent(inout) :: work_array(n_genes)
    !| Permutation array for sorting (dimension n_genes)
    integer(int32), intent(inout) :: perm(n_genes)
    !| Stack array for left indices during sorting
    integer(int32), intent(inout) :: stack_left(n_genes)
    !| Stack array for right indices during sorting
    integer(int32), intent(inout) :: stack_right(n_genes)
    !| Output boolean array indicating outliers
    logical, intent(out) :: is_outlier(n_genes)
    !| Reference x-coordinates.
    real(real64), intent(inout) :: loess_x(n_families)
    !| Reference y-coordinates (length n_total).
    real(real64), intent(inout) :: loess_y(n_families)
    !| Indices of reference points used for smoothing.
    integer(int32), intent(inout) :: loess_n(n_families)
    !| Error code: 0=ok, 201=invalid family indices
    integer(int32), intent(out) :: ierr
    !| (optional) Percentile threshold for outlier detection (default: 95)
    real(real64), intent(in), optional :: percentile

    ! Local variables
    real(real64) :: dscale(n_families)
    real(real64) :: rdi(n_genes)
    real(real64) :: threshold
    integer(int32) :: i
    real(real64) :: percentile_val

    ! Set default percentile if not present
    if (present(percentile)) then
      percentile_val = percentile
    else
      percentile_val = 95.0_real64
    end if

    ! Always initialize permutation array
    do i = 1, n_genes
      perm(i) = i
    end do

    call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, loess_n, ierr)
    if (is_err(ierr)) return
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, work_array, perm, stack_left, stack_right)
    call identify_outliers(n_genes, rdi, work_array, is_outlier, threshold, percentile_val)
  end subroutine detect_outliers
end module tox_get_outliers

!> R wrapper for compute_family_scaling (main version with automatic allocation).
!| Calls compute_family_scaling_alloc with standard Fortran types for R interface.
!| This is the recommended version for most users as it handles memory allocation automatically.
subroutine compute_family_scaling_r(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, ierr)
  use tox_get_outliers, only: compute_family_scaling_alloc
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  !| Total number of genes
  integer(int32), intent(in) :: n_genes
  !| Total number of gene families
  integer(int32), intent(in) :: n_families
  !| Array of Euclidean distances for each gene
  real(real64), intent(in) :: distances(n_genes)
  !| Mapping of each gene to its family (1-based)
  integer(int32), intent(in) :: gene_to_fam(n_genes)
  !| Output: array of scaling factors per family
  real(real64), intent(out) :: dscale(n_families)
  !| Reference x-coordinates.
  real(real64), intent(inout) :: loess_x(n_families)
  !| Reference y-coordinates (length n_total).
  real(real64), intent(inout) :: loess_y(n_families)
  !| Indices of reference points used for smoothing.
  integer(int32), intent(inout) :: indices_used(n_families)
  !| Error code: 0=ok, 201=invalid family indices
  integer(int32), intent(out) :: ierr
  call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, ierr)
end subroutine compute_family_scaling_r

!> R wrapper for compute_rdi.
!| Calls compute_rdi with standard Fortran types for R interface.
!| Outputs both unsorted and sorted RDI, permutation, and sorting workspace arrays for downstream use.
subroutine compute_rdi_r(n_genes, n_families, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
  use tox_get_outliers, only: compute_rdi
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  !| Total number of genes
  integer(int32), intent(in) :: n_genes
  !| Total number of families
  integer(int32), intent(in) :: n_families
  !| Array of Euclidean distances for each gene to its centroid
  real(real64), intent(in) :: distances(n_genes)
  !| Gene-to-family mapping (1-based indexing)
  integer(int32), intent(in) :: gene_to_fam(n_genes)
  !| Array of scaling factors for each family
  real(real64), intent(in) :: dscale(n_families)
  !| Output array of RDI values for each gene
  real(real64), intent(out) :: rdi(n_genes)
  !| Work array for sorting (dimension n_genes)
  real(real64), intent(inout) :: sorted_rdi(n_genes)
  !| Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
  integer(int32), intent(inout) :: perm(n_genes)
  !| Stack array for sorting (dimension n_genes)
  integer(int32), intent(inout) :: stack_left(n_genes)
  !| Stack array for sorting (dimension n_genes)
  integer(int32), intent(inout) :: stack_right(n_genes)
  call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
end subroutine compute_rdi_r

!> R wrapper for identify_outliers.
!| Calls identify_outliers with standard Fortran types for R interface.
subroutine identify_outliers_r(n_genes, rdi, sorted_rdi, is_outlier, threshold, percentile)
  use tox_get_outliers, only: identify_outliers
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  !| Total number of genes
  integer(int32), intent(in) :: n_genes
  !| Array of RDI values for each gene
  real(real64), intent(in) :: rdi(n_genes)
  !| Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
  real(real64), intent(in) :: sorted_rdi(n_genes)
  !| Output boolean array indicating outliers
  logical, intent(out) :: is_outlier(n_genes)
  !| Output threshold value used for detection
  real(real64), intent(out) :: threshold
  !| Percentile threshold (default: 95 for top 5%)
  real(real64), intent(in) :: percentile
  call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, percentile)
end subroutine identify_outliers_r

!> R wrapper for detect_outliers.
!| Calls detect_outliers with standard Fortran types for R interface.
subroutine detect_outliers_r(n_genes, n_families, distances, gene_to_fam, &
                            work_array, perm, stack_left, stack_right, &
                            is_outlier, loess_x, loess_y, loess_n, ierr, &
                            percentile)
  use tox_get_outliers, only: detect_outliers
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  !| Total number of genes
  integer(int32), intent(in) :: n_genes
  !| Total number of gene families
  integer(int32), intent(in) :: n_families
  !| Array of Euclidean distances for each gene to its centroid
  real(real64), intent(in) :: distances(n_genes)
  !| Gene-to-family mapping (1-based indexing)
  integer(int32), intent(in) :: gene_to_fam(n_genes)
  !| Work array for sorting (dimension n_genes)
  real(real64), intent(inout) :: work_array(n_genes)
  !| Permutation array for sorting (dimension n_genes)
  integer(int32), intent(inout) :: perm(n_genes)
  !| Stack array for left indices during sorting
  integer(int32), intent(inout) :: stack_left(n_genes)
  !| Stack array for right indices during sorting
  integer(int32), intent(inout) :: stack_right(n_genes)
  !| Output boolean array indicating outliers
  logical, intent(out) :: is_outlier(n_genes)
  !| Reference x-coordinates.
  real(real64), intent(inout) :: loess_x(n_families)
  !| Reference y-coordinates (length n_total).
  real(real64), intent(inout) :: loess_y(n_families)
  !| Indices of reference points used for smoothing.
  integer(int32), intent(inout) :: loess_n(n_families)
  !| Error code: 0=ok, 201=invalid family indices
  integer(int32), intent(out) :: ierr
  !| (optional) Percentile threshold for outlier detection (default: 95)
  real(real64), intent(in), optional :: percentile
  call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                      work_array, perm, stack_left, stack_right, &
                      is_outlier, loess_x, loess_y, loess_n, ierr, &
                      percentile)
end subroutine detect_outliers_r

!> C wrapper for compute_family_scaling (main version with automatic allocation).
!| Calls compute_family_scaling_alloc with C-compatible types for external interface.
!| This is the recommended version for most users as it handles memory allocation automatically.
subroutine compute_family_scaling_c(n_genes, n_families, distances, gene_to_fam, dscale, &
  loess_x, loess_y, indices_used, ierr) bind(C, name="compute_family_scaling_c")
use, intrinsic :: iso_c_binding, only : c_int, c_double
use tox_get_outliers, only: compute_family_scaling_alloc
implicit none
!| Total number of genes
integer(c_int), intent(in), value :: n_genes
!| Total number of families
integer(c_int), intent(in), value :: n_families
!| Array of Euclidean distances for each gene
real(c_double), intent(in), target :: distances(n_genes)
!| Mapping of each gene to its family (1-based)
integer(c_int), intent(in), target :: gene_to_fam(n_genes)
!| Output: array of scaling factors per family
real(c_double), intent(out), target :: dscale(n_families)
!| Reference x-coordinates for LOESS
real(c_double), intent(inout), target :: loess_x(n_families)
!| Reference y-coordinates for LOESS
real(c_double), intent(inout), target :: loess_y(n_families)
!| Indices of reference points used for smoothing
integer(c_int), intent(inout), target :: indices_used(n_families)
!| Error code: 0=ok, 201=invalid family indices
integer(c_int), intent(out) :: ierr
  call compute_family_scaling_alloc(n_genes, n_families, distances, gene_to_fam, dscale, &
    loess_x, loess_y, indices_used, ierr)
end subroutine compute_family_scaling_c

!> C wrapper for compute_rdi.
!| Calls compute_rdi with C-compatible types for external interface.
!| Outputs both unsorted and sorted RDI, permutation, and sorting workspace arrays for downstream use.
subroutine compute_rdi_c(n_genes, n_families, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right) bind(C, name="compute_rdi_c")
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use tox_get_outliers, only: compute_rdi
  implicit none

  !| Total number of genes
  integer(c_int), intent(in), value :: n_genes
  !| Total number of families
  integer(c_int), intent(in), value :: n_families
  !| Array of Euclidean distances for each gene to its centroid
  real(c_double), intent(in), target :: distances(n_genes)
  !| Gene-to-family mapping (1-based indexing)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  !| Array of scaling factors for each family
  real(c_double), intent(in), target :: dscale(n_families)
  !| Output array of RDI values for each gene (unsorted)
  real(c_double), intent(out), target :: rdi(n_genes)
  !| Output array of sorted RDI values (filtered, sorted)
  real(c_double), intent(out), target :: sorted_rdi(n_genes)
  !| Output permutation array for sorting (dimension n_genes)
  integer(c_int), intent(out), target :: perm(n_genes)
  !| Output stack array for left indices during sorting
  integer(c_int), intent(out), target :: stack_left(n_genes)
  !| Output stack array for right indices during sorting
  integer(c_int), intent(out), target :: stack_right(n_genes)
  call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
end subroutine compute_rdi_c

!> C wrapper for identify_outliers.
!| Calls identify_outliers with C-compatible types for external interface.
subroutine identify_outliers_c(n_genes, rdi, sorted_rdi, is_outlier_int, threshold, percentile) &
                              bind(C, name="identify_outliers_c")
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use tox_get_outliers, only: identify_outliers
  implicit none

  !| Total number of genes
  integer(c_int), intent(in), value :: n_genes
  !| Array of RDI values for each gene
  real(c_double), intent(in), target :: rdi(n_genes)
  !| Sorted RDI array (must be sorted in ascending order before calling)
  real(c_double), intent(in), target :: sorted_rdi(n_genes)
  !| Output integer array indicating outliers (1=outlier, 0=not)
  integer(c_int), intent(out), target :: is_outlier_int(n_genes)
  !| Output threshold value used for detection
  real(c_double), intent(out) :: threshold
  !| Percentile threshold for outlier detection
  real(c_double), intent(in), value :: percentile
  logical :: is_outlier(n_genes)
  integer :: i

  ! Convert integer (0/1) to logical (.false./.true.) for is_outlier
  do i = 1, n_genes
    is_outlier(i) = (is_outlier_int(i) /= 0)
  end do

  call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, percentile)
  ! Convert logical (.true./.false.) to integer (1/0)
  do i = 1, n_genes
    if (is_outlier(i)) then
      is_outlier_int(i) = 1
    else
      is_outlier_int(i) = 0
    end if
  end do
end subroutine identify_outliers_c


!> C wrapper for detect_outliers.
!| Calls detect_outliers with C-compatible types for external interface.
subroutine detect_outliers_c(n_genes, n_families, distances, gene_to_fam, &
                          work_array, perm, stack_left, stack_right, &
                          is_outlier_int, loess_x, loess_y, loess_n, ierr, &
                          percentile) bind(C, name="detect_outliers_c")
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use tox_get_outliers, only: detect_outliers
  implicit none

  !| Total number of genes
  integer(c_int), intent(in), value :: n_genes, n_families
  !| Array of Euclidean distances for each gene to its centroid
  real(c_double), intent(in), target :: distances(n_genes)
  !| Gene-to-family mapping (1-based indexing)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  !| Work array for sorting (dimension n_genes)
  real(c_double), intent(inout), target :: work_array(n_genes)
  !| Permutation array for sorting (dimension n_genes)
  integer(c_int), intent(inout), target :: perm(n_genes)
  !| Stack array for left indices during sorting
  integer(c_int), intent(inout), target :: stack_left(n_genes)
  !| Stack array for right indices during sorting
  integer(c_int), intent(inout), target :: stack_right(n_genes)
  !| Output integer array indicating outliers (1=outlier, 0=not)
  integer(c_int), intent(out), target :: is_outlier_int(n_genes)
  !| Reference x-coordinates for LOESS
  real(c_double), intent(inout), target :: loess_x(n_families)
  !| Reference y-coordinates for LOESS
  real(c_double), intent(inout), target :: loess_y(n_families)
  !| Indices of reference points used for smoothing
  integer(c_int), intent(inout), target :: loess_n(n_families)
  !| Error code: 0=ok, 201=invalid family indices
  integer(c_int), intent(out) :: ierr
  !| Percentile threshold for outlier detection
  real(c_double), intent(in), value :: percentile
  logical :: is_outlier(n_genes)
  integer :: i
  call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                    work_array, perm, stack_left, stack_right, &
                    is_outlier, loess_x, loess_y, loess_n, ierr, &
                    percentile)

  ! Convert logical (.true./.false.) to integer (1/0) for is_outlier
  do i = 1, n_genes
    if (is_outlier(i)) then
      is_outlier_int(i) = 1
    else
      is_outlier_int(i) = 0
    end if
  end do
end subroutine detect_outliers_c

!> C wrapper for compute_family_scaling expert version.
!| Calls compute_family_scaling with C-compatible types for external interface.
!| This wrapper is designed for external use, providing additional arguments for advanced configurations.
subroutine compute_family_scaling_expert_c(n_genes, n_families, distances, gene_to_fam, dscale, &
  loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, &
  iv, liv, wv, lv, diagl, w_init, z_mat, rw, ww, res, pi, yhat_tmp, &
  span, degree, mode, n_iters, ierr) bind(C, name="compute_family_scaling_expert_c")

  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use tox_get_outliers, only: compute_family_scaling
  implicit none
  
  !| Total number of genes
  integer(c_int), intent(in) :: n_genes
  !| Total number of families
  integer(c_int), intent(in) :: n_families
  !| Array of Euclidean distances for each gene
  real(c_double), intent(in), target :: distances(n_genes)
  !| Mapping of each gene to its family (1-based)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  !| Output: array of scaling factors per family
  real(c_double), intent(out), target :: dscale(n_families)
  !| Reference x-coordinates for LOESS
  real(c_double), intent(inout), target :: loess_x(n_families)
  !| Reference y-coordinates for LOESS
  real(c_double), intent(inout), target :: loess_y(n_families)
  !| Indices of reference points used for smoothing
  integer(c_int), intent(inout), target :: indices_used(n_families)
  !| Temporary array for permutation
  integer(c_int), intent(inout), target :: perm_tmp(n_genes)
  !| Temporary array for left stack
  integer(c_int), intent(inout), target :: stack_left_tmp(n_genes)
  !| Temporary array for right stack
  integer(c_int), intent(inout), target :: stack_right_tmp(n_genes)
  !| Temporary array for family distances
  real(c_double), intent(inout), target :: family_distances(n_genes)
  !| Integer workspace array for LOESS
  integer(c_int), intent(inout), target :: iv(liv)
  !| Length of integer workspace array
  integer(c_int), intent(in) :: liv
  !| Length of real workspace array
  integer(c_int), intent(in) :: lv
  !| Real workspace array for LOESS
  real(c_double), intent(inout), target :: wv(lv)
  !| Diagonal elements for LOESS
  real(c_double), intent(inout), target :: diagl(n_genes)
  !| Initial weights for LOESS
  real(c_double), intent(inout), target :: w_init(n_genes)
  !| Z matrix for LOESS
  real(c_double), intent(inout), target :: z_mat(n_genes, 1)
  !| Residual weights for LOESS
  real(c_double), intent(inout), target :: rw(n_genes)
  !| Working weights for LOESS
  real(c_double), intent(inout), target :: ww(n_genes)
  !| Residuals for LOESS
  real(c_double), intent(inout), target :: res(n_genes)
  !| Pi values for LOESS
  integer(c_int), intent(inout), target :: pi(n_genes)
  !| Temporary array for predicted values
  real(c_double), intent(inout), target :: yhat_tmp(n_genes)
  !| Span parameter for LOESS
  real(c_double), intent(in) :: span
  !| Degree of polynomial for LOESS
  integer(c_int), intent(in) :: degree
  !| Mode for LOESS
  integer(c_int), intent(in) :: mode
  !| Number of iterations for LOESS
  integer(c_int), intent(in) :: n_iters
  !| Error code: 0=ok, 201=invalid family indices
  integer(c_int), intent(out) :: ierr

  call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
    loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, &
    iv, liv, wv, lv, diagl, w_init, z_mat, rw, ww, res, pi, yhat_tmp, &
    span, degree, mode, n_iters, ierr)
end subroutine compute_family_scaling_expert_c
