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
    !! could not read array type error
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
    !! could not write array type error
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

  !> Creates the error code for an error and a specific argument position (if argument related)
  pure integer(int32) function create_err_code(error, arg_pos) result(ierr)
    integer(int32), intent(in)    :: error
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    if (present(arg_pos)) then
      ierr = 10000 * arg_pos + error
    else
      ierr = error
    end if
  end function create_err_code

  !> Extracts the error from the error code
  pure integer(int32) function get_err_code(ierr) result(error)
    integer(int32), intent(in) :: ierr
      !! Error code

    error = mod(ierr, 10000)
  end function get_err_code

  !> Extracts the argument position from the error code
  pure integer(int32) function get_err_arg_pos(ierr) result(arg_pos)
    integer(int32), intent(in) :: ierr
      !! Error code

    arg_pos = ierr / 10000
  end function get_err_arg_pos

  !> Maps a specific argument position in the encoded Error code to a new one. Helpful for nested subroutine calls.
  pure subroutine map_err_arg_pos(ierr, old, new)
    integer(int32), intent(inout) :: ierr
      !! Error code
    integer(int32), intent(in) :: old
      !! Argument position that should be mapped
    integer(int32), intent(in) :: new
      !! New argument position
  
    if (get_err_arg_pos(ierr) == old) then
      ierr = create_err_code(get_err_code(ierr), new)
    end if
  end subroutine map_err_arg_pos

  !> set the error code to OK, use at beginning of procedures
  pure subroutine set_ok(ierr)
    integer(int32), intent(out) :: ierr
      !! Error code
    ierr = ERR_OK
  end subroutine set_ok

  !> set the error code to specific error
  pure subroutine set_err(ierr, error, arg_pos)
    integer(int32), intent(inout) :: ierr
      !! Error code
    integer(int32), intent(in)    :: error
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    ierr = create_err_code(error, arg_pos)
  end subroutine set_err

  !> set the error code only if it is currently OK, use to prevent overwriting first error
  pure subroutine set_err_once(ierr, error, arg_pos)
    integer(int32), intent(inout) :: ierr
      !! Error code
    integer(int32), intent(in)    :: error
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related
    if (.not. is_err(ierr)) call set_err(ierr, error, arg_pos)
  end subroutine set_err_once

  !> check if the error code indicates error
  pure logical function is_err(ierr) result(not_ok)
    integer(int32), intent(in) :: ierr
      !! Error code
    not_ok = (get_err_code(ierr) /= ERR_OK)
  end function is_err

  !> check if the error code indicates success
  pure logical function is_ok(ierr) result(ok)
    integer(int32), intent(in) :: ierr
      !! Error code
    ok = .not. is_err(ierr)
  end function is_ok

  !Checks if allocation is successful
  pure subroutine check_io_stat(ios, ierr)
    integer(int32), intent(in) :: ios
    integer(int32), intent(inout) :: ierr
      !! Error code
    if(is_err(ios)) call set_err(ierr, ERR_ALLOC_FAIL)
  end subroutine check_io_stat
  
  pure subroutine validate_dimension_size(n, ierr, arg_pos)
    integer(int32), intent(in), optional :: n
    integer(int32), intent(inout) :: ierr
      !! Error code
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    if (present(n)) then
      if(n < 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)  
      else if(n == 0) then
        call set_err_once(ierr, ERR_EMPTY_INPUT, arg_pos)
      end if
    end if
  end subroutine

  !> validate that the actual type error matches the expected type error
  !! Closes the connected unit in case of an error
  subroutine validate_type_code(actual, expected, unit, ierr, arg_pos)
    integer(int32), intent(in) :: actual
    !! Actual type error read from file
    integer(int32), intent(in) :: expected
    !! Expected type error
    integer(int32), intent(in) :: unit
    !! Unit connection
    integer(int32), intent(inout) :: ierr
      !! Error code
    !! Error code
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    if(actual /= expected) then
      close(unit)
      call set_err_once(ierr, ERR_TYPE_MISMATCH, arg_pos)
    end if
  end subroutine

  !> Validates min<=e<=max for all elements e of an array
  pure subroutine validate_in_range_int(val, ierr, arg_pos, min, max)
    integer(int32), intent(in), optional :: val
      !! value to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
      !! Error code
    integer(int32), intent(in), optional :: min
      !! lower bound for a value, default is lowest 32-bit integer -> -huge(1_int32)
    integer(int32), intent(in), optional :: max
      !! upper bound for a value, default is largest 32-bit integer -> huge(1_int32)
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    integer(int32) :: actual_min, actual_max

    if (present(val)) then
      M_DEFAULT_VAL(min, actual_min, -huge(1_int32))
      M_DEFAULT_VAL(max, actual_max, huge(1_int32))
    
      if ((val < actual_min) .or. (val > actual_max)) then
        call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
      end if
    end if
  end subroutine validate_in_range_int

  !> Validates min<=e<=max for all elements e of an array
  pure subroutine validate_all_in_range_int(array, n_elements, ierr, arg_pos, min, max)
    integer(int32), intent(in) :: n_elements
      !! Size of `array`
    integer(int32), dimension(n_elements), intent(in), optional :: array
      !! Array to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
      !! Error code
    integer(int32), intent(in), optional :: min
      !! lower bound for a value, default is lowest 32-bit integer -> -huge(1_int32)
    integer(int32), intent(in), optional :: max
      !! upper bound for a value, default is largest 32-bit integer -> huge(1_int32)
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    integer(int32) :: i_element
    
    if (present(array)) then
      do concurrent (i_element = 1:n_elements) shared(ierr, arg_pos, min, max, array)
        call validate_in_range_int(array(i_element), ierr, arg_pos, min, max)
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
  pure subroutine validate_in_range_real(val, ierr, arg_pos, min, max)
    real(real64), intent(in), optional :: val
      !! value to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
      !! Error code
    real(real64), intent(in), optional :: min
      !! lower bound for a value, default is lowest 64-bit float -> -huge(1.0_real64)
    real(real64), intent(in), optional :: max
      !! upper bound for a value, default is largest 64-bit float -> huge(1.0_real64)
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    real(real64) :: actual_min, actual_max

    if (present(val)) then
      M_DEFAULT_VAL(min, actual_min, -huge(1.0_real64))
      M_DEFAULT_VAL(max, actual_max, huge(1.0_real64))
    
      if (ieee_is_nan(val) .or. .not. ieee_is_finite(val)) then
        call set_err_once(ierr, ERR_NAN_INF, arg_pos)
      else if ((val < actual_min) .or. (val > actual_max)) then
        call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
      end if
    end if
  end subroutine validate_in_range_real

  !> Validates min<=e<=max AND e/=NaN for all elements e of an array
  pure subroutine validate_all_in_range_real(array, n_elements, ierr, arg_pos, min, max)
    integer(int32), intent(in) :: n_elements
      !! Size of `array`
    real(real64), dimension(n_elements), intent(in), optional :: array
      !! Array to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
      !! Error code
    real(real64), intent(in), optional :: min
      !! lower bound for a value, default is lowest 64-bit float -> -huge(1.0_real64)
    real(real64), intent(in), optional :: max
      !! upper bound for a value, default is largest 64-bit float -> huge(1.0_real64)
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    integer(int32) :: i_element
  
    if (present(array)) then
      do concurrent (i_element = 1:n_elements) shared(ierr, arg_pos, min, max, array)
        call validate_in_range_real(array(i_element), ierr, arg_pos, min, max)
      end do
    end if
  end subroutine validate_all_in_range_real

  !> Strictly validates a distance matrix of euclidean distances.
  !| This means that distance X->X is exactly zero (no tolerance)
  !| and X->Y is exactly the same as Y->X (no tolerance).
  pure subroutine validate_distance_matrix(distances, n, ierr, arg_pos, min, max)
    integer(int32), intent(in) :: n
      !! Number of columns and rows of `distances`
    real(real64), dimension(n, n), intent(in), optional :: distances
      !! Matrix to be validated
    integer(int32), intent(inout) :: ierr
      !! Error code
      !! Error code
    real(real64), intent(in), optional :: min
      !! lower bound for a distance value, default is zero to have only positives
    real(real64), intent(in), optional :: max
      !! upper bound for a distance value, default is largest 64-bit float -> huge(1.0_real64)
    integer(int32), intent(in), optional :: arg_pos
      !! Position of the validated argument that triggered the error, default: 0 -> not argument related

    integer(int32) :: i_col, i_row
    real(real64) :: actual_min

    if (present(distances)) then
      M_DEFAULT_VAL(min, actual_min, 0.0_real64)

      call validate_all_in_range_real(distances, n * n, ierr, arg_pos, actual_min, max)
      if (is_err(ierr)) return

      do concurrent (i_col = 1:n) shared(ierr, distances, arg_pos)
        if (distances(i_col, i_col) /= 0.0_real64) then
          call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
        else
          do concurrent (i_row = 1:n) shared(ierr, i_col, distances, arg_pos)
            if (distances(i_row, i_col) /= distances(i_col, i_row)) then
              call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
            end if
          end do
        end if
      end do
    end if
  end subroutine validate_distance_matrix

end module tox_errors
