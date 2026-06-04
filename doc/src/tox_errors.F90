#include "macros.h"

!> error handling module for tensor-omics
module tox_errors
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  implicit none
  public   ! <-- expose all names (constants + procedures)

  !------------------------------
  ! Success
  !------------------------------
  integer(int32), parameter :: ERR_OK                 = 0
    !! no error, operation successful
  !------------------------------
  ! 1xx: I/O / File reading
  !------------------------------
  integer(int32), parameter :: ERR_FILE_OPEN          = 101   
    !! could not open file
  integer(int32), parameter :: ERR_READ_MAGIC         = 102
    !! could not read magic number
  integer(int32), parameter :: ERR_READ_TYPE          = 103   
    !! could not read array type code
  integer(int32), parameter :: ERR_READ_NDIMS         = 104   
    !! could not read number of dimensions
  integer(int32), parameter :: ERR_READ_DIMS          = 105  
    !! could not read array dimensions
  integer(int32), parameter :: ERR_READ_CHARLEN       = 106
    !! could not read character length
  integer(int32), parameter :: ERR_READ_DATA          = 107   
    !! could not read array data
  integer(int32), parameter :: ERR_WRITE_MAGIC        = 112
    !! could not write magic number
  integer(int32), parameter :: ERR_WRITE_TYPE         = 113
    !! could not write array type code
  integer(int32), parameter :: ERR_WRITE_NDIMS        = 114
    !! could not write number of dimensions
  integer(int32), parameter :: ERR_WRITE_DIMS         = 115
    !! could not write array dimensions
  integer(int32), parameter :: ERR_WRITE_CHARLEN      = 116
    !! could not write character length
  integer(int32), parameter :: ERR_WRITE_DATA         = 117
    !! could not write array data
  integer(int32), parameter :: ERR_FILE_ADD           = 121
    !! Could not add file to archive
  integer(int32), parameter :: ERR_FILE_EXTRACT       = 122
    !! Could not extract file from archive
  integer(int32), parameter :: ERR_MISSING_MANIFEST   = 123
    !! Manifest in zip file is missing 
  integer(int32), parameter :: ERR_FILE_CLOSE         = 124
    !! Failed to close the file
  !------------------------------
  ! 2xx: Format / Input validation
  !------------------------------
  integer(int32), parameter :: ERR_INVALID_FORMAT     = 200 
    !! invalid format detected
  integer(int32), parameter :: ERR_INVALID_INPUT      = 201 
    !! invalid input arguments
  integer(int32), parameter :: ERR_EMPTY_INPUT        = 202 
    !! empty input arrays
  integer(int32), parameter :: ERR_DIM_MISMATCH       = 203  
    !! dimensions do not match expected shape
  integer(int32), parameter :: ERR_NAN_INF            = 204  
    !! NaN or Inf found where not allowed
  integer(int32), parameter :: ERR_UNSUPPORTED_TYPE   = 205 
    !! unsupported data type encountered
  integer(int32), parameter :: ERR_SIZE_MISMATCH      = 206
    !! Array size mismatch
  integer(int32), parameter :: ERR_TYPE_MISMATCH      = 207
    !! Array type read does not match expected type
  integer(int32), parameter :: ERR_STRING_TOO_LONG    = 208
    !! String exceeds buffer size 
  integer(int32), parameter :: ERR_IDX_OUT_OF_BOUNDS  = 209
    !! Array index out of bounds
  integer(int32), parameter :: ERR_DIVISION_BY_ZERO   = 210
    !! Division by zero encountered

  !------------------------------
  ! 3xx: Memory
  !------------------------------
  integer(int32), parameter :: ERR_ALLOC_FAIL         = 301 
    !! memory allocation failed
  integer(int32), parameter :: ERR_POINTER_NULL       = 302 
    !! null pointer dereference

  !------------------------------
  ! 5xxx: Fortran runtime / Unit state
  ! (Keep 5002 for compatibility with existing R mapping)
  !------------------------------
  integer(int32), parameter :: ERR_UNIT_NOT_CONNECTED = 5002
    !! Fortran runtime error: unit not connected
    
  !------------------------------
  ! 9xxx: Internal / Unknown
  !------------------------------
  integer(int32), parameter :: ERR_INTERNAL           = 9001
    !! unexpected internal state or logic error
  integer(int32), parameter :: ERR_UNKNOWN            = 9999 
    !! unknown error

contains

  !> set the error code to OK, use at beginning of procedures
  pure subroutine set_ok(ierr)
    integer(int32), intent(out) :: ierr
    ierr = ERR_OK
  end subroutine set_ok

  !> set the error code to specific code
  pure subroutine set_err(ierr, code)
    integer(int32), intent(inout) :: ierr
    integer(int32), intent(in)    :: code
    if (ierr == ERR_OK) ierr = code
  end subroutine set_err

  !> set the error code only if it is currently OK, use to prevent overwriting first error
  pure subroutine set_err_once(ierr, code)
    integer(int32), intent(inout) :: ierr
    integer(int32), intent(in)    :: code
    if (ierr == ERR_OK) call set_err(ierr, code)
  end subroutine set_err_once

  !> check if the error code indicates error
  pure logical function is_err(ierr) result(not_ok)
    integer(int32), intent(in) :: ierr
    not_ok = (ierr /= ERR_OK)
  end function is_err

  !> check if the error code indicates success
  pure logical function is_ok(ierr) result(ok)
    integer(int32), intent(in) :: ierr
    ok = .not. is_err(ierr)
  end function is_ok

  !Checks if allocation is successful
  pure subroutine check_io_stat(ios, ierr)
    integer(int32), intent(in) :: ios
    integer(int32), intent(out) :: ierr
    if(is_err(ios)) call set_err(ierr, ERR_ALLOC_FAIL)
  end subroutine check_io_stat

  pure subroutine check_alloc_stat(stat, ierr)
    integer(int32), intent(in)    :: stat
    integer(int32), intent(inout) :: ierr
    if (is_err(stat)) call set_err_once(ierr, ERR_ALLOC_FAIL)
  end subroutine check_alloc_stat

  
  pure subroutine validate_dimension_size(n, ierr)
    integer(int32), intent(in), optional :: n
    integer(int32), intent(inout) :: ierr

    if (present(n)) then
      if(n < 0) call set_err(ierr, ERR_INVALID_INPUT)  
      if(n == 0) call set_err(ierr, ERR_EMPTY_INPUT)
    end if
  end subroutine

  !> validate that the actual type code matches the expected type code
  !! Closes the connected unit in case of an error
  subroutine validate_type_code(actual, expected, unit, ierr)
    integer(int32), intent(in) :: actual
    !! Actual type code read from file
    integer(int32), intent(in) :: expected
    !! Expected type code
    integer(int32), intent(in) :: unit
    !! Unit connection
    integer(int32), intent(out) :: ierr
    !! Error code
    if(actual /= expected) then
      close(unit)
      call set_err(ierr, ERR_TYPE_MISMATCH)
    end if
  end subroutine

  !> Validates min<=e<=max for all elements e of an array
  pure subroutine validate_in_range_int(val, ierr, min, max)
    integer(int32), intent(in), optional :: val
      !! value to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
    integer(int32), intent(in), optional :: min
      !! lower bound for a value, default is lowest 32-bit integer -> -huge(1_int32)
    integer(int32), intent(in), optional :: max
      !! upper bound for a value, default is largest 32-bit integer -> huge(1_int32)

    integer(int32) :: actual_min, actual_max

    if (present(val)) then
      M_DEFAULT_VAL(min, actual_min, -huge(1_int32))
      M_DEFAULT_VAL(max, actual_max, huge(1_int32))
    
      if ((val < actual_min) .or. (val > actual_max)) then
        call set_err(ierr, ERR_INVALID_INPUT)
      end if
    end if
  end subroutine validate_in_range_int

  !> Validates min<=e<=max for all elements e of an array
  pure subroutine validate_all_in_range_int(array, n_elements, ierr, min, max)
    integer(int32), intent(in) :: n_elements
      !! Size of `array`
    integer(int32), dimension(n_elements), intent(in), optional :: array
      !! Array to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
    integer(int32), intent(in), optional :: min
      !! lower bound for a value, default is lowest 32-bit integer -> -huge(1_int32)
    integer(int32), intent(in), optional :: max
      !! upper bound for a value, default is largest 32-bit integer -> huge(1_int32)

    integer(int32) :: i_element
    
    if (present(array)) then
      do concurrent (i_element = 1:n_elements) shared(ierr, min, max, array)
        call validate_in_range_int(array(i_element), ierr, min, max)
      end do
    end if
  end subroutine validate_all_in_range_int

  !> Validates min<=e<=max AND e/=NaN for all elements e of an array.
  !|
  !| @note
  !| This validation is inclusive: `min<=val<=max`<br>
  !| To achieve exclusive bounds, us above/below from f42_utils,
  !| like: `validate_in_range_real(x, ierr, min=above(0.0_real64), max=below(100.0_real64))` for `0<x<100`
  !| @endnote
  pure subroutine validate_in_range_real(val, ierr, min, max)
    real(real64), intent(in), optional :: val
      !! value to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
    real(real64), intent(in), optional :: min
      !! lower bound for a value, default is lowest 64-bit float -> -huge(1.0_real64)
    real(real64), intent(in), optional :: max
      !! upper bound for a value, default is largest 64-bit float -> huge(1.0_real64)

    real(real64) :: actual_min, actual_max

    if (present(val)) then
      M_DEFAULT_VAL(min, actual_min, -huge(1.0_real64))
      M_DEFAULT_VAL(max, actual_max, huge(1.0_real64))
    
      if (ieee_is_nan(val) .or. .not. ieee_is_finite(val)) then
        call set_err(ierr, ERR_NAN_INF)
      else if ((val < actual_min) .or. (val > actual_max)) then
        call set_err(ierr, ERR_INVALID_INPUT)
      end if
    end if
  end subroutine validate_in_range_real

  !> Validates min<=e<=max AND e/=NaN for all elements e of an array
  pure subroutine validate_all_in_range_real(array, n_elements, ierr, min, max)
    integer(int32), intent(in) :: n_elements
      !! Size of `array`
    real(real64), dimension(n_elements), intent(in), optional :: array
      !! Array to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
    real(real64), intent(in), optional :: min
      !! lower bound for a value, default is lowest 64-bit float -> -huge(1.0_real64)
    real(real64), intent(in), optional :: max
      !! upper bound for a value, default is largest 64-bit float -> huge(1.0_real64)

    integer(int32) :: i_element
  
    if (present(array)) then
      do concurrent (i_element = 1:n_elements) shared(ierr, min, max, array)
        call validate_in_range_real(array(i_element), ierr, min, max)
      end do
    end if
  end subroutine validate_all_in_range_real

  !> Strictly validates a distance matrix of euclidean distances.
  !| This means that distance X->X is exactly zero (no tolerance)
  !| and X->Y is exactly the same as Y->X (no tolerance).
  pure subroutine validate_distance_matrix(distances, n, ierr, min, max)
    integer(int32), intent(in) :: n
      !! Number of columns and rows of `distances`
    real(real64), dimension(n, n), intent(in), optional :: distances
      !! Matrix to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
    real(real64), intent(in), optional :: min
      !! lower bound for a distance value, default is zero to have only positives
    real(real64), intent(in), optional :: max
      !! upper bound for a distance value, default is largest 64-bit float -> huge(1.0_real64)

    integer(int32) :: i_col, i_row
    real(real64) :: actual_min

    if (present(distances)) then
      M_DEFAULT_VAL(min, actual_min, 0.0_real64)

      call validate_all_in_range_real(distances, n * n, ierr, actual_min, max)
      if (is_err(ierr)) return

      do i_col = 1, n
        if (distances(i_col, i_col) /= 0.0_real64) call set_err(ierr, ERR_INVALID_INPUT)
        do i_row = 1, n
          if (distances(i_row, i_col) /= distances(i_col, i_row)) then
            call set_err(ierr, ERR_INVALID_INPUT)
          end if
        end do
      end do
    end if
  end subroutine validate_distance_matrix

end module tox_errors
