subroutine calculate_memory_requirements(estimates, mem_bytes)
  use TensorOmics_Interface, only: DataEstimatesType, real64, int32, int64
  implicit none
  type(DataEstimatesType), intent(in) :: estimates
  integer(int64), intent(out) :: mem_bytes
  integer(int32) :: i, r, c, k
  integer(int64) :: dt_mem

  mem_bytes = 2_int64 * (storage_size(0.0_real64)/8) * estimates%n_tissues * estimates%n_genes

  k = size(estimates%meta_col_types)
  i = count(estimates%meta_col_types == 1)
  r = count(estimates%meta_col_types == 2)
  c = count(estimates%meta_col_types == 3)

  dt_mem = 2_int64 * k * (storage_size(1_int32)/8) + &
           estimates%meta_n_rows * i * (storage_size(1_int32)/8) + &
           estimates%meta_n_rows * r * (storage_size(0.0_real64)/8) + &
           estimates%meta_n_rows * c * estimates%meta_max_char * (storage_size('a')/8) + &
           k * 256 * (storage_size('a')/8)

  mem_bytes = mem_bytes + dt_mem
end subroutine

subroutine init(self, estimates)
  use TensorOmics_Interface, only: TensorOmicsData, DataEstimatesType, real64, int32
  implicit none
  type(TensorOmicsData), intent(out) :: self
  type(DataEstimatesType), intent(in) :: estimates
  integer(int32) :: k, i, r, c

  allocate(self%vec_container(estimates%n_tissues, estimates%n_genes), source=0.0_real64)
  allocate(self%shift_vecs(estimates%n_tissues, estimates%n_genes), source=0.0_real64)
  self%next_idx = 1

  k = size(estimates%meta_col_types)
  allocate(self%meta_type_map(2, k))
  self%meta_type_map(1,:) = estimates%meta_col_types
  self%meta_type_map(2,:) = 0

  i = count(estimates%meta_col_types == 1)
  r = count(estimates%meta_col_types == 2)
  c = count(estimates%meta_col_types == 3)

  if (i > 0) allocate(self%meta_int(estimates%meta_n_rows, i))
  if (r > 0) allocate(self%meta_real(estimates%meta_n_rows, r))
  if (c > 0) allocate(self%meta_char(estimates%meta_n_rows), source=repeat(' ', estimates%meta_max_char))
  allocate(self%meta_names(k), source=repeat(' ', 256))
end subroutine

subroutine update(self, patch, indices)
  use TensorOmics_Interface, only: TensorOmicsData, real64, int32
  implicit none
  type(TensorOmicsData), intent(inout) :: self
  real(real64), intent(in) :: patch(:,:)
  integer(int32), intent(out) :: indices(:)
  integer(int32) :: i, n_patch

  n_patch = size(patch, 2)
  do i = 1, n_patch
    if (self%next_idx > size(self%vec_container, 2)) exit
    self%vec_container(:, self%next_idx) = patch(:,i)
    self%shift_vecs(:, self%next_idx) = patch(:,i) - sum(self%vec_container, dim=2)/size(self%vec_container,2)
    indices(i) = self%next_idx
    self%meta_type_map(2, self%next_idx) = self%next_idx
    self%next_idx = self%next_idx + 1
  end do
end subroutine

subroutine save(self, filename)
  use TensorOmics_Interface, only: TensorOmicsData, real64, int32
  implicit none
  type(TensorOmicsData), intent(in) :: self
  character(len=*), intent(in) :: filename
  integer(int32) :: u

  open(newunit=u, file=filename, form='unformatted', status='replace')
  write(u) self%vec_container, self%shift_vecs
  write(u) self%meta_type_map, self%meta_real, self%meta_int, self%meta_char, self%meta_names
  close(u)
end subroutine
