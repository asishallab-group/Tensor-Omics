module TensorOmics_Interface
  use, intrinsic :: iso_fortran_env, only: real64, int32, int64
  implicit none

  public :: DataEstimatesType, TensorOmicsData

  type :: DataEstimatesType
    integer(int32) :: n_genes, n_tissues
    integer(int32) :: meta_n_rows, meta_max_char
    integer(int32), allocatable :: meta_col_types(:)
  end type

  type :: TensorOmicsData
    real(real64), allocatable :: vec_container(:,:)
    real(real64), allocatable :: shift_vecs(:,:)
    integer(int32), allocatable :: meta_type_map(:,:)
    real(real64), allocatable :: meta_real(:,:)
    integer(int32), allocatable :: meta_int(:,:)
    character(len=256), allocatable :: meta_names(:)
    character(len=:), allocatable :: meta_char(:)
    integer(int32) :: next_idx
  end type

end module TensorOmics_Interface
