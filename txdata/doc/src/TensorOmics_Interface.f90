!> @file TensorOmics_Interface.f90
!! @brief Core data types for TensorOmics

module TensorOmics_Interface
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  !> @brief Metadata estimates container
  type :: DataEstimatesType
    integer(int32) :: n_genes      !< Number of genes
    integer(int32) :: n_tissues    !< Number of tissues
    integer(int32) :: meta_n_rows  !< Rows in metadata
    integer(int32) :: meta_max_char !< Max chars in metadata strings
    integer(int32), allocatable :: meta_col_types(:) !< Column types (1=int, 2=real, 3=char)
  end type

  !> @brief Main data container
  type :: TensorOmicsData
    real(real64), allocatable :: vec_container(:,:)   !< Primary data matrix
    real(real64), allocatable :: shift_vecs(:,:)      !< Shifted data matrix
    integer(int32) :: next_idx                        !< Next insertion index
  end type
end module