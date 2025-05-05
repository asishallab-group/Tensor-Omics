! === Memory Calculation ===
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

! === Initialization ===
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

! === Update ===
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

! === Save ===
subroutine save(n_tissues, n_genes, vec_container, shift_vecs, filename, filename_len)
  use TensorOmics_Interface
  integer(int32), intent(in) :: n_tissues, n_genes, filename_len
  real(real64), intent(in) :: vec_container(n_tissues, n_genes)
  real(real64), intent(in) :: shift_vecs(n_tissues, n_genes)
  character(len=filename_len), intent(in) :: filename
  integer(int32) :: u

  open(newunit=u, file=filename, form='unformatted', status='replace')
  write(u) vec_container, shift_vecs
  close(u)
end subroutine
