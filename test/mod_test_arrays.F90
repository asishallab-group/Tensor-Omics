!> Unit test suite for array utilities.
module mod_test_arrays
    use asserts
    use f42_serde_arrays_utils
    use f42_serde_arrays_deserialize
    use f42_serde_arrays_serialize
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use iso_c_binding
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_arrays() result(all_tests)

        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(28))
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
        all_tests(19) = test_case("test_logical_array_1d", test_logical_array_1d)
        all_tests(20) = test_case("test_logical_array_2d", test_logical_array_2d)
        all_tests(21) = test_case("test_logical_array_3d", test_logical_array_3d)
        all_tests(22) = test_case("test_logical_array_4d", test_logical_array_4d)
        all_tests(23) = test_case("test_logical_array_5d", test_logical_array_5d)
        all_tests(24) = test_case("test_complex_array_1d", test_complex_array_1d)
        all_tests(25) = test_case("test_complex_array_2d", test_complex_array_2d)
        all_tests(26) = test_case("test_complex_array_3d", test_complex_array_3d)
        all_tests(27) = test_case("test_complex_array_4d", test_complex_array_4d)
        all_tests(28) = test_case("test_complex_array_5d", test_complex_array_5d)
    end function get_all_tests_arrays

    ! ================================================================
    ! Integer tests
    ! ================================================================

    !> Test integer array for 1D
    subroutine test_integer_array_1d()
        integer(int32), allocatable :: iarr1d(:), iarr1d2(:)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (iarr1d(5)); iarr1d = [10, 20, 30, 40, 50]

        fname = "iarr1d.test.bin"
        call serialize_int_1d(iarr1d, fname, ierr)
        print *, "Serialized integer 1D array to ", trim(fname)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        ! Metadata auslesen
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, INTEGER_TYPE_CODE, "Type code mismatch")

        ! Array basierend auf Metadaten allokieren
        allocate (iarr1d2(dims(1)))

        call deserialize_int_1d(iarr1d2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_int(iarr1d, iarr1d2, size(iarr1d), "Mismatch")
    end subroutine test_integer_array_1d

    !> Test integer array for 2D
    subroutine test_integer_array_2d()
        integer(int32), allocatable :: iarr(:, :), iarr2(:, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (iarr(2, 3)); iarr = reshape([1, 2, 3, 4, 5, 6], [2, 3])

        fname = "iarr2d.test.bin"
        call serialize_int_2d(iarr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, INTEGER_TYPE_CODE, "Type code mismatch")
        allocate (iarr2(dims(1), dims(2)))
        call deserialize_int_2d(iarr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
    end subroutine test_integer_array_2d

    !> Test integer array for 3D
    subroutine test_integer_array_3d()
        integer(int32), allocatable :: iarr(:, :, :), iarr2(:, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (iarr(2, 2, 2)); iarr = reshape([1, 2, 3, 4, 5, 6, 7, 8], [2, 2, 2])

        fname = "iarr3d.test.bin"
        call serialize_int_3d(iarr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, INTEGER_TYPE_CODE, "Type code mismatch")

        allocate (iarr2(dims(1), dims(2), dims(3)))
        call deserialize_int_3d(iarr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
    end subroutine test_integer_array_3d

    !> Test integer array for 4D
    subroutine test_integer_array_4d()
        integer(int32), allocatable :: iarr(:, :, :, :), iarr2(:, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        integer(int32) :: i
        call set_ok(ierr)
        allocate (iarr(2, 2, 1, 2)); iarr = reshape([(i, i=1, 8)], [2, 2, 1, 2])
        fname = "iarr4d.test.bin"

        call serialize_int_4d(iarr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, INTEGER_TYPE_CODE, "Type code mismatch")
        allocate (iarr2(dims(1), dims(2), dims(3), dims(4)))

        call deserialize_int_4d(iarr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
    end subroutine test_integer_array_4d

    !> Test integer array for 5D
    subroutine test_integer_array_5d()
        integer(int32), allocatable :: iarr(:, :, :, :, :), iarr2(:, :, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        integer(int32) :: i
        call set_ok(ierr)
        allocate (iarr(2, 1, 2, 1, 2)); iarr = reshape([(i, i=1, 8)], [2, 1, 2, 1, 2])
        fname = "iarr5d.test.bin"
        call serialize_int_5d(iarr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, INTEGER_TYPE_CODE, "Type code mismatch")
        allocate (iarr2(dims(1), dims(2), dims(3), dims(4), dims(5)))

        call deserialize_int_5d(iarr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
    end subroutine test_integer_array_5d

    !> Test integer array for 1x1 (edge case)
    subroutine test_integer_array_1x1()
        integer(int32), allocatable :: iarr(:, :), iarr2(:, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (iarr(1, 1)); iarr = 42
        fname = "iarr_1x1.test.bin"

        call serialize_int_2d(iarr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, INTEGER_TYPE_CODE, "Type code mismatch")

        allocate (iarr2(dims(1), dims(2)))
        call deserialize_int_2d(iarr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
    end subroutine test_integer_array_1x1

    !> Test integer array for empty case
    subroutine test_integer_array_empty()
        integer(int32), allocatable :: iarr(:, :), iarr2(:, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (iarr(0, 3))
        fname = "iarr_empty.test.bin"
        call serialize_int_2d(iarr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, INTEGER_TYPE_CODE, "Type code mismatch")
        allocate (iarr2(dims(1), dims(2)))
        call deserialize_int_2d(iarr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_int(iarr, iarr2, size(iarr), "Mismatch")
    end subroutine test_integer_array_empty

    ! ================================================================
    ! Real tests
    ! ================================================================
    !> Test real array for 1D
    subroutine test_real_array_1d()
        real(real64), allocatable :: arr(:), arr2(:)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (arr(4)); arr = [1.1d0, 2.2d0, 3.3d0, 4.4d0]
        fname = "rarr1d.test.bin"
        call serialize_real_1d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, REAL_TYPE_CODE, "Type code mismatch")
        allocate (arr2(dims(1)))
        call deserialize_real_1d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
    end subroutine test_real_array_1d

    !> Test real array for 2D
    subroutine test_real_array_2d()
        real(real64), allocatable :: arr(:, :), arr2(:, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (arr(2, 2)); arr = reshape([1.5d0, 2.5d0, 3.5d0, 4.5d0], [2, 2])
        fname = "rarr2d.test.bin"
        call serialize_real_2d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, REAL_TYPE_CODE, "Type code mismatch")
        allocate (arr2(dims(1), dims(2)))
        call deserialize_real_2d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
    end subroutine test_real_array_2d

    !> Test real array for 3D
    subroutine test_real_array_3d()
        real(real64), allocatable :: arr(:, :, :), arr2(:, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (arr(2, 2, 2)); arr = reshape([1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0], [2, 2, 2])
        fname = "rarr3d.test.bin"
        call serialize_real_3d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, REAL_TYPE_CODE, "Type code mismatch")
        allocate (arr2(dims(1), dims(2), dims(3)))
        call deserialize_real_3d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
    end subroutine test_real_array_3d

    !> Test real array for 4D
    subroutine test_real_array_4d()
        real(real64), allocatable :: arr(:, :, :, :), arr2(:, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        integer(int32) :: i
        call set_ok(ierr)
        allocate (arr(2, 2, 1, 2)); arr = reshape([(real(i, real64), i=1, 8)], [2, 2, 1, 2])
        fname = "rarr4d.test.bin"
        call serialize_real_4d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, REAL_TYPE_CODE, "Type code mismatch")
        allocate (arr2(dims(1), dims(2), dims(3), dims(4)))
        call deserialize_real_4d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
    end subroutine test_real_array_4d

    !> Test real array for 5D
    subroutine test_real_array_5d()
        real(real64), allocatable :: arr(:, :, :, :, :), arr2(:, :, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        integer(int32) :: i
        call set_ok(ierr)

        allocate (arr(2, 1, 2, 1, 2)); arr = reshape([(real(i, real64), i=1, 8)], [2, 1, 2, 1, 2])
        fname = "rarr5d.test.bin"
        call serialize_real_5d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, REAL_TYPE_CODE, "Type code mismatch")
        allocate (arr2(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call deserialize_real_5d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_real(arr, arr2, size(arr), 1d-12, "Mismatch")
    end subroutine test_real_array_5d

    ! ================================================================
    ! Char tests
    ! ================================================================
    !> Test char array for 1D
    subroutine test_char_array_1d()
        character(len=3) :: arr(3), arr2(3)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        arr = ['foo', 'bar', 'baz']
        fname = "carr1d.test.bin"
        call serialize_char_1d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, len(arr), "Type code mismatch")
        call assert_equal_int(ndims, rank(arr), "ndims mismatch")
        call assert_equal_array_int(dims, shape(arr), rank(arr), "shape mismatch")

        call deserialize_char_1d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_char(arr, arr2, len(arr), size(arr), "Mismatch")
    end subroutine test_char_array_1d

    !> Test char array for 2D
    subroutine test_char_array_2d()
        character(len=5) :: arr(2, 2), arr2(2, 2)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)

        arr = reshape(['foo  ', 'bar  ', 'baz  ', 'quxxx'], [2, 2])
        fname = "carr2d.test.bin"
        call serialize_char_2d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, len(arr), "Type code mismatch")
        call assert_equal_int(ndims, rank(arr), "ndims mismatch")
        call assert_equal_array_int(dims, shape(arr), rank(arr), "shape mismatch")

        call deserialize_char_2d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_char(arr, arr2, len(arr), size(arr), "Mismatch")
    end subroutine test_char_array_2d

    !> Test char array for 3D
    subroutine test_char_array_3d()
        character(len=5) :: arr(2, 2, 1), arr2(2, 2, 1)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        arr = reshape(['foo  ', 'bar  ', 'baz  ', 'qux  '], [2, 2, 1])
        fname = "carr3d.test.bin"
        call serialize_char_3d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, len(arr), "Type code mismatch")
        call assert_equal_int(ndims, rank(arr), "ndims mismatch")

        call deserialize_char_3d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_char(arr, arr2, len(arr), size(arr), "Mismatch")
    end subroutine test_char_array_3d

    !> Test char array for 4D
    subroutine test_char_array_4d()
        character(len=5) :: arr(2, 1, 1, 2), arr2(2, 1, 1, 2)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        arr = reshape(['foo  ', 'bar  ', 'baz  ', 'qux  '], [2, 1, 1, 2])
        fname = "carr4d.test.bin"
        call serialize_char_4d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, len(arr), "Type code mismatch")
        call assert_equal_int(ndims, rank(arr), "ndims mismatch")

        call deserialize_char_4d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_char(arr, arr2, len(arr), size(arr), "Mismatch")
    end subroutine test_char_array_4d

    !> Test char array for 5D
    subroutine test_char_array_5d()
        character(len=5) :: arr(2, 1, 2, 1, 2), arr2(2, 1, 2, 1, 2)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        arr = reshape(['foo  ', 'bar  ', 'baz  ', 'qux  ', 'aaa  ', 'bbb  ', 'ccc  ', 'ddd  '], [2, 1, 2, 1, 2])
        fname = "carr5d.test.bin"
        call serialize_char_5d(arr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, len(arr), "Type code mismatch")
        call assert_equal_int(ndims, rank(arr), "ndims mismatch")

        call deserialize_char_5d(arr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_char(arr, arr2, len(arr), size(arr), "Mismatch")
    end subroutine test_char_array_5d

    !> Test char array with protein names (longer strings)
    subroutine test_char_array_protein()
        character(len=10) :: protein(2, 1, 2, 1, 2), protein2(2, 1, 2, 1, 2)
        integer(int32) :: ierr, type_code, ndims, dims(5)
        character(len=*), parameter :: fname = "proteins_arr.test.bin"
        call set_ok(ierr)
        protein = reshape(['METHIONINE', 'GLYCINE   ', 'SERINE    ', 'LYSINE    ', &
                           'VALINE    ', 'HISTIDINE ', 'PROLINE   ', 'LEUCINE   '], [2, 1, 2, 1, 2])
        call serialize_char_5d(protein, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, len(protein), "Type code mismatch")
        call assert_equal_int(ndims, rank(protein), "ndims mismatch")

        call deserialize_char_5d(protein2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_char(protein, protein2, len(protein), size(protein), "Mismatch")
    end subroutine test_char_array_protein

    !> Test logical array for 1D
    subroutine test_logical_array_1d()
        logical, allocatable :: larr1d(:), larr1d2(:)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (larr1d(4)); larr1d = [.true., .false., .true., .false.]

        fname = "larr1d.test.bin"
        call serialize_logical_1d(larr1d, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        ! Metadata auslesen
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, LOGICAL_TYPE_CODE, "Type code mismatch")

        ! Array basierend auf Metadaten allokieren
        allocate (larr1d2(dims(1)))

        call deserialize_logical_1d(larr1d2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_logical(larr1d, larr1d2, size(larr1d), "Logical 1D Mismatch")
    end subroutine test_logical_array_1d

    !> Test logical array for 2D
    subroutine test_logical_array_2d()
        logical, allocatable :: larr(:, :), larr2(:, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (larr(2, 3)); larr = reshape([.true., .false., .true., .false., .true., .false.], [2, 3])

        fname = "larr2d.test.bin"
        call serialize_logical_2d(larr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, LOGICAL_TYPE_CODE, "Type code mismatch")
        allocate (larr2(dims(1), dims(2)))
        call deserialize_logical_2d(larr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_logical(larr, larr2, size(larr), "Logical 2D Mismatch")
    end subroutine test_logical_array_2d

    !> Test logical array for 3D
    subroutine test_logical_array_3d()
        logical, allocatable :: larr(:, :, :), larr2(:, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (larr(2, 2, 2)); larr = reshape([.true., .false., .true., .false., .true., .false., .true., .false.], [2, 2, 2])

        fname = "larr3d.test.bin"
        call serialize_logical_3d(larr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, LOGICAL_TYPE_CODE, "Type code mismatch")

        allocate (larr2(dims(1), dims(2), dims(3)))
        call deserialize_logical_3d(larr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_logical(larr, larr2, size(larr), "Logical 3D Mismatch")
    end subroutine test_logical_array_3d

    !> Test logical array for 4D
    subroutine test_logical_array_4d()
        logical, allocatable :: larr(:, :, :, :), larr2(:, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (larr(2, 2, 1, 2)); larr = reshape([.true., .false., .true., .false., .true., .false., .true., .false.], [2, 2, 1, 2])
        fname = "larr4d.test.bin"

        call serialize_logical_4d(larr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, LOGICAL_TYPE_CODE, "Type code mismatch")
        allocate (larr2(dims(1), dims(2), dims(3), dims(4)))

        call deserialize_logical_4d(larr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_logical(larr, larr2, size(larr), "Logical 4D Mismatch")
    end subroutine test_logical_array_4d

    !> Test logical array for 5D
    subroutine test_logical_array_5d()
        logical, allocatable :: larr(:, :, :, :, :), larr2(:, :, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (larr(2, 1, 2, 1, 2)); larr = reshape([.true., .false., .true., .false., .true., .false., .true., .false.], [2, 1, 2, 1, 2])
        fname = "larr5d.test.bin"
        call serialize_logical_5d(larr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, LOGICAL_TYPE_CODE, "Type code mismatch")
        allocate (larr2(dims(1), dims(2), dims(3), dims(4), dims(5)))

        call deserialize_logical_5d(larr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_logical(larr, larr2, size(larr), "Logical 5D Mismatch")
    end subroutine test_logical_array_5d

    ! ================================================================
    ! Complex tests
    ! ================================================================
    !> Test complex array for 1D
    subroutine test_complex_array_1d()
        complex(real64), allocatable :: carr(:), carr2(:)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (carr(3)); 
        carr = [(1.0d0, 2.0d0), (3.0d0, 4.0d0), (5.0d0, 6.0d0)]
        fname = "carr1d.test.bin"
        call serialize_complex_1d(carr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, COMPLEX_TYPE_CODE, "Type code mismatch")
        allocate (carr2(dims(1)))
        call deserialize_complex_1d(carr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_complex(carr, carr2, size(carr), 1d-12, "Complex 1D Mismatch")
    end subroutine test_complex_array_1d

    !> Test complex array for 2D
    subroutine test_complex_array_2d()
        complex(real64), allocatable :: carr(:, :), carr2(:, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (carr(2, 2)); 
        carr = reshape([(1.5d0, 0.5d0), (2.5d0, 1.5d0), (3.5d0, 2.5d0), (4.5d0, 3.5d0)], [2, 2])
        fname = "carr2d.test.bin"
        call serialize_complex_2d(carr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, COMPLEX_TYPE_CODE, "Type code mismatch")
        allocate (carr2(dims(1), dims(2)))
        call deserialize_complex_2d(carr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_complex(carr, carr2, size(carr), 1d-12, "Complex 2D Mismatch")
    end subroutine test_complex_array_2d

    !> Test complex array for 3D
    subroutine test_complex_array_3d()
        complex(real64), allocatable :: carr(:, :, :), carr2(:, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (carr(2, 2, 2)); 
        carr = reshape([(1.0d0, 0.0d0), (0.0d0, 1.0d0), (1.0d0, 1.0d0), (0.0d0, 0.0d0), &
                        (2.0d0, 0.0d0), (0.0d0, 2.0d0), (2.0d0, 2.0d0), (0.0d0, 0.0d0)], [2, 2, 2])
        fname = "carr3d.test.bin"
        call serialize_complex_3d(carr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, COMPLEX_TYPE_CODE, "Type code mismatch")
        allocate (carr2(dims(1), dims(2), dims(3)))
        call deserialize_complex_3d(carr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_complex(carr, carr2, size(carr), 1d-12, "Complex 3D Mismatch")
    end subroutine test_complex_array_3d

    !> Test complex array for 4D
    subroutine test_complex_array_4d()
        complex(real64), allocatable :: carr(:, :, :, :), carr2(:, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)
        allocate (carr(2, 2, 1, 2))

        ! Manual assignment for complex array
        carr(1, 1, 1, 1) = (1.0d0, 0.5d0)
        carr(2, 1, 1, 1) = (2.0d0, 1.0d0)
        carr(1, 2, 1, 1) = (3.0d0, 1.5d0)
        carr(2, 2, 1, 1) = (4.0d0, 2.0d0)
        carr(1, 1, 1, 2) = (5.0d0, 2.5d0)
        carr(2, 1, 1, 2) = (6.0d0, 3.0d0)
        carr(1, 2, 1, 2) = (7.0d0, 3.5d0)
        carr(2, 2, 1, 2) = (8.0d0, 4.0d0)

        fname = "carr4d.test.bin"
        call serialize_complex_4d(carr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, COMPLEX_TYPE_CODE, "Type code mismatch")
        allocate (carr2(dims(1), dims(2), dims(3), dims(4)))
        call deserialize_complex_4d(carr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_complex(carr, carr2, size(carr), 1d-12, "Complex 4D Mismatch")
    end subroutine test_complex_array_4d

    !> Test complex array for 5D
    subroutine test_complex_array_5d()
        complex(real64), allocatable :: carr(:, :, :, :, :), carr2(:, :, :, :, :)
        character(len=100) :: fname
        integer(int32) :: ierr, type_code, ndims, dims(5)
        call set_ok(ierr)

        allocate (carr(2, 1, 2, 1, 2))

        ! Manual assignment for complex array
        carr(1, 1, 1, 1, 1) = (1.0d0, 0.25d0)
        carr(2, 1, 1, 1, 1) = (2.0d0, 0.5d0)
        carr(1, 1, 2, 1, 1) = (3.0d0, 0.75d0)
        carr(2, 1, 2, 1, 1) = (4.0d0, 1.0d0)
        carr(1, 1, 1, 1, 2) = (5.0d0, 1.25d0)
        carr(2, 1, 1, 1, 2) = (6.0d0, 1.5d0)
        carr(1, 1, 2, 1, 2) = (7.0d0, 1.75d0)
        carr(2, 1, 2, 1, 2) = (8.0d0, 2.0d0)

        fname = "carr5d.test.bin"
        call serialize_complex_5d(carr, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")

        call get_array_metadata(fname, dims, 5, ndims, type_code, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_int(type_code, COMPLEX_TYPE_CODE, "Type code mismatch")
        allocate (carr2(dims(1), dims(2), dims(3), dims(4), dims(5)))
        call deserialize_complex_5d(carr2, fname, ierr)
        call assert_equal_int(ierr, ERR_OK, "Unexpected error")
        call assert_equal_array_complex(carr, carr2, size(carr), 1d-12, "Complex 5D Mismatch")
    end subroutine test_complex_array_5d

end module mod_test_arrays
