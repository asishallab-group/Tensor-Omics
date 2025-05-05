!> @file TensorOmics_mod.f90
!! @brief Core functions for data storage management


!> @brief Calculates memory requirements for the data container
!!
!! Computes total memory needed for both numerical arrays and metadata storage.
!! @param[in] n_tissues Number of tissues
!! @param[in] n_genes Number of genes
!! @param[in] meta_n_rows Number of metadata rows
!! @param[in] meta_max_char Maximum character length in metadata strings
!! @param[in] meta_col_types Array of column type codes (1=int, 2=real, 3=char)
!! @param[in] n_cols Number of metadata columns
!! @param[out] mem_bytes Total required memory in bytes
subroutine calculate_memory_requirements(n_tissues, n_genes, meta_n_rows, meta_max_char, meta_col_types, n_cols, mem_bytes)
  use TensorOmics_Interface
  integer(int32), intent(in) :: n_tissues, n_genes, meta_n_rows, meta_max_char, n_cols
  integer(int32), intent(in) :: meta_col_types(n_cols)
  integer(int32), intent(out) :: mem_bytes
  type(DataEstimatesType) :: est
  integer(int32) :: i, r, c, k
  integer(int32) :: dt_mem

  est%n_tissues = n_tissues
  est%n_genes = n_genes
  est%meta_n_rows = meta_n_rows
  est%meta_max_char = meta_max_char
  allocate(est%meta_col_types(n_cols))
  est%meta_col_types = meta_col_types

  mem_bytes = 2 * storage_size(0.0_real64)/8 * n_tissues * n_genes

  k = size(est%meta_col_types)
  i = count(est%meta_col_types == 1)
  r = count(est%meta_col_types == 2)
  c = count(est%meta_col_types == 3)

  dt_mem = 2 * k * (storage_size(1_int32)/8) + &
           est%meta_n_rows * i * (storage_size(1_int32)/8) + &
           est%meta_n_rows * r * (storage_size(0.0_real64)/8) + &
           est%meta_n_rows * c * est%meta_max_char * (storage_size('a')/8)

  mem_bytes = mem_bytes + dt_mem
end subroutine

!> @brief Initializes the TensorOmics data container
!!
!! Allocates and zeros the primary storage matrices.
!! @param[in] n_tissues Number of tissues
!! @param[in] n_genes Number of genes
!! @param[in] meta_n_rows Number of metadata rows
!! @param[in] meta_max_char Maximum character length in metadata strings
!! @param[in] meta_col_types Array of column type codes
!! @param[in] n_cols Number of metadata columns
!! @param[out] vec_out Allocated vector container (n_tissues × n_genes)
!! @param[out] shift_out Allocated shift vectors (n_tissues × n_genes)
!! @param[out] next_idx Initial insertion index (always 1)
subroutine init(n_tissues, n_genes, meta_n_rows, meta_max_char, meta_col_types, n_cols, vec_out, shift_out, next_idx)
  use TensorOmics_Interface
  integer(int32), intent(in) :: n_tissues, n_genes, meta_n_rows, meta_max_char, n_cols
  integer(int32), intent(in) :: meta_col_types(n_cols)
  real(real64), intent(out) :: vec_out(n_tissues, n_genes)
  real(real64), intent(out) :: shift_out(n_tissues, n_genes)
  integer(int32), intent(out) :: next_idx
  type(TensorOmicsData) :: self

  allocate(self%vec_container(n_tissues, n_genes), source=0.0_real64)
  allocate(self%shift_vecs(n_tissues, n_genes), source=0.0_real64)
  self%next_idx = 1

  vec_out = self%vec_container
  shift_out = self%shift_vecs
  next_idx = self%next_idx
end subroutine

!> @brief Updates the container with new data patches
!!
!! Inserts new columns while maintaining running statistics.
!! @param[in] n_tissues Number of tissues
!! @param[in] n_genes Number of genes (max capacity)
!! @param[in] patch Data patch to insert (n_tissues × n_patch)
!! @param[in] n_patch Number of columns in patch
!! @param[in,out] vec_container Primary data matrix
!! @param[in,out] shift_vecs Shifted data matrix
!! @param[in,out] next_idx Next insertion index (auto-incremented)
!! @param[out] indices Array of inserted column indices
subroutine update(n_tissues, n_genes, patch, n_patch, vec_container, shift_vecs, next_idx, indices)
  use TensorOmics_Interface
  integer(int32), intent(in) :: n_tissues, n_genes, n_patch
  real(real64), intent(in) :: patch(n_tissues, n_patch)
  real(real64), intent(inout) :: vec_container(n_tissues, n_genes)
  real(real64), intent(inout) :: shift_vecs(n_tissues, n_genes)
  integer(int32), intent(inout) :: next_idx
  integer(int32), intent(out) :: indices(n_patch)
  integer(int32) :: i

  do i = 1, n_patch
    if (next_idx > n_genes) exit
    vec_container(:, next_idx) = patch(:,i)
    shift_vecs(:, next_idx) = patch(:,i) - sum(vec_container, dim=2)/n_genes
    indices(i) = next_idx
    next_idx = next_idx + 1
  end do
end subroutine

!> @brief Saves data to binary file (stream format)
!!
!! Uses byte-level I/O for cross-platform compatibility with R.
!! @param[in] n_tissues Number of tissues
!! @param[in] n_genes Number of genes
!! @param[in] vec_container Primary data matrix
!! @param[in] shift_vecs Shifted data matrix
!! @param[in] filename_bytes Filename as ASCII byte array
!! @param[in] filename_len Length of filename
subroutine save(n_tissues, n_genes, vec_container, shift_vecs, filename_bytes, filename_len)
  use TensorOmics_Interface
  integer(int32), intent(in) :: n_tissues, n_genes, filename_len
  real(real64), intent(in) :: vec_container(n_tissues, n_genes)
  real(real64), intent(in) :: shift_vecs(n_tissues, n_genes)
  integer(int32), intent(in) :: filename_bytes(filename_len)
  character(len=filename_len) :: filename
  integer(int32) :: i, u

  ! Convert bytes to filename string
  do i = 1, filename_len
    filename(i:i) = achar(filename_bytes(i))
  end do

  open(newunit=u, file=filename, access='stream', form='unformatted', status='replace')
  write(u) vec_container, shift_vecs
  close(u)
end subroutine
