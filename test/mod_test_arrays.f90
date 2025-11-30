module mod_test_arrays
    use asserts
    use array_utils, only: get_array_metadata
    use int_deserialize_mod
    use real_deserialize_mod
    use char_deserialize_mod
    use serialize_char
    use serialize_int
    use serialize_real
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use iso_c_binding
    use tox_errors
    implicit none
    PUBLIC

    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    type :: test_case
        character(len=64) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

contains

  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(18)
    all_tests(1) = test_case("test_integer_array_1d", test_integer_array_1d)
    all_tests(2) = test_case("test_integer_array_2d", test_integer_array_2d)
    all_tests(3) = test_case("test_integer_array_3d", test_integer_array_3d)
    all_tests(4) = test_case("test_integer_array_4d", test_integer_array_4d)
    all_tests(5) = test_case("test_integer_array_5d", test_integer_array_5d)
    all_tests(6) = test_case("test_real_array_1d", test_real_array_1d)
    all_tests(7) = test_case("test_real_array_2d", test_real_array_2d)
    all_tests(8) = test_case("test_real_array_3d", test_real_array_3d)
    all_tests(9) = test_case("test_real_array_4d", test_real_array_4d)
    all_tests(10) = test_case("test_real_array_5d", test_real_array_5d)
    all_tests(11) = test_case("test_char_array_1d", test_char_array_1d)
    all_tests(12) = test_case("test_char_array_2d", test_char_array_2d)
    all_tests(13) = test_case("test_char_array_3d", test_char_array_3d)
    all_tests(14) = test_case("test_char_array_4d", test_char_array_4d)
    all_tests(15) = test_case("test_char_array_5d", test_char_array_5d)
    all_tests(16) = test_case("test_integer_array_1x1", test_integer_array_1x1)
    all_tests(17) = test_case("test_integer_array_empty", test_integer_array_empty)
    all_tests(18) = test_case("test_char_array_protein", test_char_array_protein)
  end function get_all_tests

  subroutine run_all_tests_array()
    type(test_case) :: all_tests(18)
    integer :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All array tests passed successfully."
  end subroutine run_all_tests_array

  subroutine run_named_tests_array(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(18)
    integer :: i, j
    logical :: found
    all_tests = get_all_tests()
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) print *, "Unknown test: ", trim(test_names(i))
    end do
  end subroutine run_named_tests_array

  ! ================================================================
  ! Integer tests
  ! ================================================================
  subroutine test_integer_array_1d()
      integer(int32), allocatable :: iarr1d(:), iarr1d2(:)
      character(len=100) :: fname
      integer(int32) :: ierr, ndims, dims(5)
      call set_ok(ierr)
      allocate(iarr1d(5)); iarr1d = [10,20,30,40,50]

      fname = "test_iarr1d.bin"
      call serialize_int_1d(iarr1d, fname, ierr)
      if (.not. is_ok(ierr)) error stop

      ! Metadata auslesen
      call get_array_metadata(fname, dims, 5, ndims, ierr)
      if (.not. is_ok(ierr)) error stop
      
      ! Array basierend auf Metadaten allokieren
      allocate(iarr1d2(dims(1)))
      
      call deserialize_int_1d(iarr1d2, fname, ierr)
      if (.not. is_ok(ierr)) error stop
      call assert_equal_array_int(iarr1d, iarr1d2, size(iarr1d), "Mismatch")
  end subroutine test_integer_array_1d

  subroutine test_integer_array_2d()
    integer(int32), allocatable :: iarr(:,:), iarr2(:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    call set_ok(ierr)
    allocate(iarr(2,3)); iarr = reshape([1,2,3,4,5,6],[2,3])

    fname = "test_iarr2d.bin"
    call serialize_int_2d(iarr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if(.not. is_ok(ierr)) error stop
    allocate(iarr2(dims(1),dims(2)))
    call deserialize_int_2d(iarr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
  end subroutine test_integer_array_2d

  subroutine test_integer_array_3d()
    integer(int32), allocatable :: iarr(:,:,:), iarr2(:,:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    call set_ok(ierr)
    allocate(iarr(2,2,2)); iarr = reshape([1,2,3,4,5,6,7,8],[2,2,2])

    fname = "test_iarr3d.bin"
    call serialize_int_3d(iarr, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    
    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop

    allocate(iarr2(dims(1),dims(2),dims(3)))
    call deserialize_int_3d(iarr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
  end subroutine test_integer_array_3d

  subroutine test_integer_array_4d()
    integer(int32), allocatable :: iarr(:,:,:,:), iarr2(:,:,:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    integer(int32) :: i
    call set_ok(ierr)
    allocate(iarr(2,2,1,2)); iarr = reshape([(i, i=1,8)],[2,2,1,2])
    fname = "test_iarr4d.bin"

    call serialize_int_4d(iarr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(iarr2(dims(1),dims(2),dims(3),dims(4)))

    call deserialize_int_4d(iarr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
  end subroutine test_integer_array_4d

  subroutine test_integer_array_5d()
    integer(int32), allocatable :: iarr(:,:,:,:,:), iarr2(:,:,:,:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    integer(int32) :: i
    call set_ok(ierr)
    allocate(iarr(2,1,2,1,2)); iarr = reshape([(i,i=1,8)],[2,1,2,1,2])
    fname = "test_iarr5d.bin"
    call serialize_int_5d(iarr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(iarr2(dims(1),dims(2),dims(3),dims(4),dims(5)))

    call deserialize_int_5d(iarr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
  end subroutine test_integer_array_5d

  subroutine test_integer_array_1x1()
    integer(int32), allocatable :: iarr(:,:), iarr2(:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    call set_ok(ierr)
    allocate(iarr(1,1)); iarr = 42
    fname = "test_iarr_1x1.bin"

    call serialize_int_2d(iarr, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    
    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop

    allocate(iarr2(dims(1),dims(2)))
    call deserialize_int_2d(iarr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
  end subroutine test_integer_array_1x1

  subroutine test_integer_array_empty()
    integer(int32), allocatable :: iarr(:,:), iarr2(:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    call set_ok(ierr)
    allocate(iarr(0,3))
    fname = "test_iarr_empty.bin"
    call serialize_int_2d(iarr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(iarr2(dims(1),dims(2)))
    call deserialize_int_2d(iarr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
  end subroutine test_integer_array_empty

  ! ================================================================
  ! Real tests
  ! ================================================================
  subroutine test_real_array_1d()
    real(real64), allocatable :: arr(:), arr2(:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    call set_ok(ierr)
    allocate(arr(4)); arr = [1.1d0,2.2d0,3.3d0,4.4d0]
    fname = "test_rarr1d.bin"
    call serialize_real_1d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(arr2(dims(1)))
    call deserialize_real_1d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
  end subroutine test_real_array_1d

  subroutine test_real_array_2d()
    real(real64), allocatable :: arr(:,:), arr2(:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    call set_ok(ierr)
    allocate(arr(2,2)); arr = reshape([1.5d0,2.5d0,3.5d0,4.5d0],[2,2])
    fname = "test_rarr2d.bin"
    call serialize_real_2d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(arr2(dims(1),dims(2)))
    call deserialize_real_2d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
  end subroutine test_real_array_2d

  subroutine test_real_array_3d()
    real(real64), allocatable :: arr(:,:,:), arr2(:,:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    call set_ok(ierr)
    allocate(arr(2,2,2)); arr = reshape([1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0],[2,2,2])
    fname = "test_rarr3d.bin"
    call serialize_real_3d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(arr2(dims(1),dims(2),dims(3)))
    call deserialize_real_3d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
  end subroutine test_real_array_3d

  subroutine test_real_array_4d()
    real(real64), allocatable :: arr(:,:,:,:), arr2(:,:,:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    integer(int32) :: i
    call set_ok(ierr)
    allocate(arr(2,2,1,2)); arr = reshape([(real(i,real64),i=1,8)],[2,2,1,2])
    fname = "test_rarr4d.bin"
    call serialize_real_4d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(arr2(dims(1),dims(2),dims(3),dims(4)))
    call deserialize_real_4d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
  end subroutine test_real_array_4d

  subroutine test_real_array_5d()
    real(real64), allocatable :: arr(:,:,:,:,:), arr2(:,:,:,:,:)
    character(len=100) :: fname
    integer(int32) :: ierr, ndims, dims(5)
    integer(int32) :: i
    call set_ok(ierr)

    allocate(arr(2,1,2,1,2)); arr = reshape([(real(i,real64),i=1,8)],[2,1,2,1,2])
    fname = "test_rarr5d.bin"
    call serialize_real_5d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr)
    if (.not. is_ok(ierr)) error stop
    allocate(arr2(dims(1),dims(2),dims(3),dims(4),dims(5)))
    call deserialize_real_5d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
  end subroutine test_real_array_5d

  ! ================================================================
  ! Char tests
  ! ================================================================
  subroutine test_char_array_1d()
    character(len=:), allocatable :: arr(:), arr2(:)
    character(len=100) :: fname
    integer(int32) :: clen, ierr, ndims, dims(5)
    clen = 3
    call set_ok(ierr)
    allocate(character(len=clen) :: arr(3))
    arr = ['foo','bar','baz']
    fname = "test_carr1d.bin"
    call serialize_char_1d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr, clen)
    if (.not. is_ok(ierr)) error stop
    allocate(character(len=clen) :: arr2(dims(1)))
    call deserialize_char_1d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_char(arr, arr2, clen, size(arr), "Mismatch")
  end subroutine test_char_array_1d

  subroutine test_char_array_2d()
    character(len=:), allocatable :: arr(:,:), arr2(:,:)
    character(len=100) :: fname
    integer(int32) :: clen, ierr, ndims, dims(5)
    clen = 5
    call set_ok(ierr)
    allocate(character(len=clen) :: arr(2,2))
    arr = reshape(['foo  ','bar  ','baz  ','quxxx'],[2,2])
    fname = "test_carr2d.bin"
    call serialize_char_2d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call get_array_metadata(fname, dims, 5, ndims, ierr, clen)
    if (.not. is_ok(ierr)) error stop

    allocate(character(len=clen) :: arr2(dims(1),dims(2)))
    call deserialize_char_2d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_char(arr, arr2, clen, size(arr), "Mismatch")
  end subroutine test_char_array_2d

  subroutine test_char_array_3d()
    character(len=:), allocatable :: arr(:,:,:), arr2(:,:,:)
    character(len=100) :: fname
    integer(int32) :: clen, ierr, ndims, dims(5)
    clen = 5
    call set_ok(ierr)
    allocate(character(len=clen) :: arr(2,2,1))
    arr = reshape(['foo  ','bar  ','baz  ','qux  '],[2,2,1])
    fname = "test_carr3d.bin"
    call serialize_char_3d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr, clen)
    if (.not. is_ok(ierr)) error stop  

    allocate(character(len=clen) :: arr2(dims(1),dims(2),dims(3)))
    call deserialize_char_3d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_char(arr, arr2, clen, size(arr), "Mismatch")
  end subroutine test_char_array_3d

  subroutine test_char_array_4d()
    character(len=:), allocatable :: arr(:,:,:,:), arr2(:,:,:,:)
    character(len=100) :: fname
    integer(int32) :: clen, ierr, ndims, dims(5)
    clen = 5
    call set_ok(ierr)
    allocate(character(len=clen) :: arr(2,1,1,2))
    arr = reshape(['foo  ','bar  ','baz  ','qux  '],[2,1,1,2])
    fname = "test_carr4d.bin"
    call serialize_char_4d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop

    call get_array_metadata(fname, dims, 5, ndims, ierr, clen)
    if (.not. is_ok(ierr)) error stop

    allocate(character(len=clen) :: arr2(dims(1),dims(2),dims(3),dims(4)))
    call deserialize_char_4d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_char(arr, arr2, clen, size(arr), "Mismatch")
  end subroutine test_char_array_4d

  subroutine test_char_array_5d()
    character(len=:), allocatable :: arr(:,:,:,:,:), arr2(:,:,:,:,:)
    character(len=100) :: fname
    integer(int32) :: clen, ierr, ndims, dims(5)
    clen = 5
    call set_ok(ierr)
    allocate(character(len=clen) :: arr(2,1,2,1,2))
    arr = reshape(['foo  ','bar  ','baz  ','qux  ','aaa  ','bbb  ','ccc  ','ddd  '],[2,1,2,1,2])
    fname = "test_carr5d.bin"
    call serialize_char_5d(arr, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call get_array_metadata(fname, dims, 5, ndims, ierr, clen)
    if (.not. is_ok(ierr)) error stop

    allocate(character(len=clen) :: arr2(dims(1),dims(2),dims(3),dims(4),dims(5)))
    call deserialize_char_5d(arr2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_char(arr, arr2, clen, size(arr), "Mismatch")
  end subroutine test_char_array_5d

  subroutine test_char_array_protein()
    character(len=:), allocatable :: protein(:,:,:,:,:), protein2(:,:,:,:,:)
    integer(int32) :: ierr, clen, ndims, dims(5)
    character(len=*), parameter :: fname = "test_proteins_arr.bin"
    call set_ok(ierr)
    clen=10
    allocate(character(len=clen) :: protein(2,1,2,1,2))
    protein = reshape(['METHIONINE','GLYCINE   ','SERINE    ','LYSINE    ', &
                       'VALINE    ','HISTIDINE ','PROLINE   ','LEUCINE   '], [2,1,2,1,2])
    call serialize_char_5d(protein, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call get_array_metadata(fname, dims, 5, ndims, ierr, clen)
    if (.not. is_ok(ierr)) error stop

    allocate(character(len=clen) :: protein2(dims(1),dims(2),dims(3),dims(4),dims(5)))
    call deserialize_char_5d(protein2, fname, ierr)
    if (.not. is_ok(ierr)) error stop
    call assert_equal_array_char(protein, protein2, clen, size(protein), "Mismatch")
  end subroutine test_char_array_protein

end module mod_test_arrays
