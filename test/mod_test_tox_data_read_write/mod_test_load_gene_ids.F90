!> The `load_gene_ids` cases: gene ids, a 1-d array of strings, read from a file written byte for
!| byte in the published layout, every broken file refused with its code, and a round trip through
!| `save_gene_ids`. The checks are the suite fixtures' `check_*` for a 1-d array of strings, shared
!| with every `load_*` of that kind.
module mod_test_load_gene_ids
    use tox_data_read_write, only: load_gene_ids, save_gene_ids
    use test_suite, only: test_case
    use mod_test_tox_data_read_write_fixtures, only: check_load_refuses_strings
    use mod_test_tox_data_read_write_fixtures, only: check_loads_strings
    use mod_test_tox_data_read_write_fixtures, only: check_round_trip_strings
    implicit none
    private
    public :: get_all_tests_load_gene_ids

contains

    !> Every case of `load_gene_ids`.
    function get_all_tests_load_gene_ids() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(3))
        all_tests(1) = test_case("test_load_gene_ids_layout", test_load_gene_ids_layout)
        all_tests(2) = test_case("test_load_gene_ids_refuses_broken_files", &
                                 test_load_gene_ids_refuses_broken_files)
        all_tests(3) = test_case("test_load_gene_ids_round_trip", test_load_gene_ids_round_trip)
    end function get_all_tests_load_gene_ids

    !> In place of the old `test_read_write_gene_ids`, which round-tripped the shipped gene ids and
    !| compared only the first 128 characters of each 256: gene ids are read from a file written
    !| byte for byte in the published layout.
    subroutine test_load_gene_ids_layout()
        call check_loads_strings(load_gene_ids, "test_load_gene_ids_layout")
    end subroutine test_load_gene_ids_layout

    !> Every broken file is refused with its code: an empty one, a short header at each of its
    !| fields, a wrong magic number, a number of dimensions out of range, another type, rank or
    !| extent, and a short payload.
    subroutine test_load_gene_ids_refuses_broken_files()
        call check_load_refuses_strings(load_gene_ids, "test_load_gene_ids_refuses_broken_files")
    end subroutine test_load_gene_ids_refuses_broken_files

    !> What `save_gene_ids` writes is read back unchanged.
    subroutine test_load_gene_ids_round_trip()
        call check_round_trip_strings(save_gene_ids, load_gene_ids, "test_load_gene_ids_round_trip")
    end subroutine test_load_gene_ids_round_trip
end module mod_test_load_gene_ids
