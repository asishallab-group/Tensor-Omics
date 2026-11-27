!> The `load_gene_to_family` cases: the gene-to-family mapping, a 1-d int32 array, read from a file
!| written byte for byte in the published layout, every broken file refused with its code, and a
!| round trip through `save_gene_to_family`. The checks are the suite fixtures' `check_*` for a 1-d
!| int32 array, shared with every `load_*` of that kind.
module mod_test_load_gene_to_family
    use tox_data_read_write, only: load_gene_to_family, save_gene_to_family
    use test_suite, only: test_case
    use mod_test_tox_data_read_write_fixtures, only: check_load_refuses_a_missing_file_integers
    use mod_test_tox_data_read_write_fixtures, only: check_load_refuses_integers
    use mod_test_tox_data_read_write_fixtures, only: check_loads_integers
    use mod_test_tox_data_read_write_fixtures, only: check_round_trip_integers
    implicit none
    private
    public :: get_all_tests_load_gene_to_family

contains

    !> Every case of `load_gene_to_family`.
    function get_all_tests_load_gene_to_family() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(4))
        all_tests(1) = test_case("test_load_gene_to_family_layout", test_load_gene_to_family_layout)
        all_tests(2) = test_case("test_load_gene_to_family_refuses_broken_files", &
                                 test_load_gene_to_family_refuses_broken_files)
        all_tests(3) = test_case("test_load_gene_to_family_round_trip", test_load_gene_to_family_round_trip)
        all_tests(4) = test_case("test_load_gene_to_family_refuses_a_missing_file", &
                                 test_load_gene_to_family_refuses_a_missing_file)
    end function get_all_tests_load_gene_to_family

    !> In place of the old `test_read_write_gene_to_fam`, which round-tripped the shipped mapping:
    !| the gene-to-family mapping are read from a file written byte for byte in the published
    !| layout.
    subroutine test_load_gene_to_family_layout()
        call check_loads_integers(load_gene_to_family, "test_load_gene_to_family_layout")
    end subroutine test_load_gene_to_family_layout

    !> Every broken file is refused with its code: an empty one, a short header at each of its
    !| fields, a wrong magic number, a number of dimensions out of range, another type, rank or
    !| extent, and a short payload.
    subroutine test_load_gene_to_family_refuses_broken_files()
        call check_load_refuses_integers(load_gene_to_family, "test_load_gene_to_family_refuses_broken_files")
    end subroutine test_load_gene_to_family_refuses_broken_files

    !> What `save_gene_to_family` writes is read back unchanged.
    subroutine test_load_gene_to_family_round_trip()
        call check_round_trip_integers(save_gene_to_family, load_gene_to_family, &
                                    "test_load_gene_to_family_round_trip")
    end subroutine test_load_gene_to_family_round_trip

    !> A file that does not exist is refused with ERR_FILE_OPEN, and the error unit stays connected.
    subroutine test_load_gene_to_family_refuses_a_missing_file()
        call check_load_refuses_a_missing_file_integers(load_gene_to_family, &
            "test_load_gene_to_family_refuses_a_missing_file")
    end subroutine test_load_gene_to_family_refuses_a_missing_file
end module mod_test_load_gene_to_family
