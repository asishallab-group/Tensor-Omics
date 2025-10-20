!> error handling module for tensor-omics
module tox_errors
  use iso_fortran_env, only: int32
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

  pure subroutine validate_dimension_size(n, ierr)
    integer(int32), intent(in) :: n
    integer(int32), intent(out) :: ierr
    if(n < 0) call set_err(ierr, ERR_INVALID_INPUT)  
    if(n == 0) call set_err(ierr, ERR_EMPTY_INPUT)
  end subroutine

end module tox_errors