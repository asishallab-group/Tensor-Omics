!> The `load_shift_vectors` cases: shift vectors, a 2-d real64 array, read from a file written byte
!| for byte in the published layout, every broken file refused with its code, and a round trip
!| through `save_shift_vectors`. The checks are the suite fixtures' `check_*` for a 2-d real64
!| array, shared with every `load_*` of that kind.
module mod_test_load_shift_vectors
    use tox_data_read_write, only: load_shift_vectors, save_shift_vectors
    use test_suite, only: test_case
    use mod_test_tox_data_read_write_fixtures, only: check_load_refuses_reals
    use mod_test_tox_data_read_write_fixtures, only: check_loads_non_finite_reals
    use mod_test_tox_data_read_write_fixtures, only: check_loads_reals
    use mod_test_tox_data_read_write_fixtures, only: check_round_trip_reals
    implicit none
    private
    public :: get_all_tests_load_shift_vectors

contains

    !> Every case of `load_shift_vectors`.
    function get_all_tests_load_shift_vectors() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(4))
        all_tests(1) = test_case("test_load_shift_vectors_layout", test_load_shift_vectors_layout)
        all_tests(2) = test_case("test_load_shift_vectors_refuses_broken_files", &
                                 test_load_shift_vectors_refuses_broken_files)
        all_tests(3) = test_case("test_load_shift_vectors_round_trip", test_load_shift_vectors_round_trip)
        all_tests(4) = test_case("test_load_shift_vectors_non_finite", test_load_shift_vectors_non_finite)
    end function get_all_tests_load_shift_vectors

    !> In place of the old `test_write_read_shift_vectors`, which round-tripped the shipped shift
    !| vectors through save_expression_vectors and load_expression_vectors, so that no case called
    !| this procedure: shift vectors are read from a file written byte for byte in the published
    !| layout.
    subroutine test_load_shift_vectors_layout()
        call check_loads_reals(load_shift_vectors, "test_load_shift_vectors_layout")
    end subroutine test_load_shift_vectors_layout

    !> Every broken file is refused with its code: an empty one, a short header at each of its
    !| fields, a wrong magic number, a number of dimensions out of range, another type, rank or
    !| extent, and a short payload.
    subroutine test_load_shift_vectors_refuses_broken_files()
        call check_load_refuses_reals(load_shift_vectors, "test_load_shift_vectors_refuses_broken_files")
    end subroutine test_load_shift_vectors_refuses_broken_files

    !> What `save_shift_vectors` writes is read back unchanged.
    subroutine test_load_shift_vectors_round_trip()
        call check_round_trip_reals(save_shift_vectors, load_shift_vectors, &
                                    "test_load_shift_vectors_round_trip")
    end subroutine test_load_shift_vectors_round_trip

    !> Pins today's behaviour: NaN and Inf are read bit for bit, without an error.
    subroutine test_load_shift_vectors_non_finite()
        call check_loads_non_finite_reals(load_shift_vectors, "test_load_shift_vectors_non_finite")
    end subroutine test_load_shift_vectors_non_finite
end module mod_test_load_shift_vectors
