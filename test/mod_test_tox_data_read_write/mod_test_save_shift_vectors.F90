!> The `save_shift_vectors` cases: shift vectors, a 2-d real64 array, written in the published
!| layout byte for byte, and a longer file in the way replaced. The checks are the suite fixtures'
!| `check_*` for a 2-d real64 array, shared with every `save_*` of that kind.
module mod_test_save_shift_vectors
    use tox_data_read_write, only: save_shift_vectors
    use test_suite, only: test_case
    use mod_test_tox_data_read_write_fixtures, only: check_save_replaces_reals
    use mod_test_tox_data_read_write_fixtures, only: check_saves_non_finite_reals
    use mod_test_tox_data_read_write_fixtures, only: check_saves_reals
    implicit none
    private
    public :: get_all_tests_save_shift_vectors

contains

    !> Every case of `save_shift_vectors`.
    function get_all_tests_save_shift_vectors() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(3))
        all_tests(1) = test_case("test_save_shift_vectors_layout", test_save_shift_vectors_layout)
        all_tests(2) = test_case("test_save_shift_vectors_replaces_a_file", &
                                 test_save_shift_vectors_replaces_a_file)
        all_tests(3) = test_case("test_save_shift_vectors_non_finite", test_save_shift_vectors_non_finite)
    end function get_all_tests_save_shift_vectors

    !> In place of the old `test_write_read_shift_vectors`, which round-tripped the shipped shift
    !| vectors through save_expression_vectors and load_expression_vectors, so that no case called
    !| this procedure: shift vectors are written in the published layout, byte for byte.
    subroutine test_save_shift_vectors_layout()
        call check_saves_reals(save_shift_vectors, "test_save_shift_vectors_layout")
    end subroutine test_save_shift_vectors_layout

    !> A longer file in the way is replaced, not overwritten in part.
    subroutine test_save_shift_vectors_replaces_a_file()
        call check_save_replaces_reals(save_shift_vectors, "test_save_shift_vectors_replaces_a_file")
    end subroutine test_save_shift_vectors_replaces_a_file

    !> Pins today's behaviour: NaN and Inf are written bit for bit, without an error.
    subroutine test_save_shift_vectors_non_finite()
        call check_saves_non_finite_reals(save_shift_vectors, "test_save_shift_vectors_non_finite")
    end subroutine test_save_shift_vectors_non_finite
end module mod_test_save_shift_vectors
