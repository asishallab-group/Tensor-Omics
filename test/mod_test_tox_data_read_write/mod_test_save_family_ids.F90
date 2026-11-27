!> The `save_family_ids` cases: family ids, a 1-d array of strings, written in the published layout
!| byte for byte, and a longer file in the way replaced. The checks are the suite fixtures'
!| `check_*` for a 1-d array of strings, shared with every `save_*` of that kind.
module mod_test_save_family_ids
    use tox_data_read_write, only: save_family_ids
    use test_suite, only: test_case
    use mod_test_tox_data_read_write_fixtures, only: check_save_refuses_an_unwritable_path_strings
    use mod_test_tox_data_read_write_fixtures, only: check_save_replaces_strings
    use mod_test_tox_data_read_write_fixtures, only: check_saves_strings
    implicit none
    private
    public :: get_all_tests_save_family_ids

contains

    !> Every case of `save_family_ids`.
    function get_all_tests_save_family_ids() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(3))
        all_tests(1) = test_case("test_save_family_ids_layout", test_save_family_ids_layout)
        all_tests(2) = test_case("test_save_family_ids_replaces_a_file", test_save_family_ids_replaces_a_file)
        all_tests(3) = test_case("test_save_family_ids_refuses_an_unwritable_path", &
                                 test_save_family_ids_refuses_an_unwritable_path)
    end function get_all_tests_save_family_ids

    !> In place of the old `test_read_write_family_ids`, which round-tripped the shipped family ids but
    !| compared only the first half of the array: it passed the character(len=256) ids to
    !| assert_equal_array_char's character(len=128) a(n), which reads the storage as n strings of 128,
    !| so ids n/2+1..n were never compared. Here: family ids are written in the published layout, byte
    !| for byte.
    subroutine test_save_family_ids_layout()
        call check_saves_strings(save_family_ids, "test_save_family_ids_layout")
    end subroutine test_save_family_ids_layout

    !> A longer file in the way is replaced, not overwritten in part.
    subroutine test_save_family_ids_replaces_a_file()
        call check_save_replaces_strings(save_family_ids, "test_save_family_ids_replaces_a_file")
    end subroutine test_save_family_ids_replaces_a_file

    !> A path in a directory that does not exist is refused with ERR_FILE_OPEN, and the error unit
    !| stays connected.
    subroutine test_save_family_ids_refuses_an_unwritable_path()
        call check_save_refuses_an_unwritable_path_strings(save_family_ids, &
            "test_save_family_ids_refuses_an_unwritable_path")
    end subroutine test_save_family_ids_refuses_an_unwritable_path
end module mod_test_save_family_ids
