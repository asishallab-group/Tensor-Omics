!> The `get_family_index` cases: the index of the first id of `family_ids` equal to `family`,
!| blanks around either not counted, or 0 when there is none.
!|
!| Unlike `get_gene_index`, a blank `family` is not refused: it finds a blank id of the list.
!| Whether it should is not documented, so that is not pinned here.
module mod_test_get_family_index
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_data_accessors, only: get_family_index
    use test_suite, only: test_case
    use mod_test_tox_data_accessors_fixtures, only: FIXTURE_FAMILY_IDS
    implicit none
    private
    public :: get_all_tests_get_family_index

contains

    !> Every case of `get_family_index`.
    function get_all_tests_get_family_index() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(3))
        all_tests(1) = test_case("test_get_family_index_found", test_get_family_index_found)
        all_tests(2) = test_case("test_get_family_index_not_found", test_get_family_index_not_found)
        all_tests(3) = test_case("test_get_family_index_first_of_two", test_get_family_index_first_of_two)
    end function get_all_tests_get_family_index

    !> Kept from the old `test_data_accessors` (OG0000001 is family 2 of the shipped data): each
    !| id gives its position, and blanks around the id, in the list or in the query, do not count.
    subroutine test_get_family_index_found()
        character(len=8) :: family_ids(2), padded(2)

        family_ids = FIXTURE_FAMILY_IDS
        call assert_equal_int(get_family_index(family_ids, "OG1"), 1, "test_get_family_index_found: OG1")
        call assert_equal_int(get_family_index(family_ids, "OG2"), 2, "test_get_family_index_found: OG2")
        call assert_equal_int(get_family_index(family_ids, " OG2 "), 2, "test_get_family_index_found: a padded query")
        padded = ["  OG7", " OG8 "]
        call assert_equal_int(get_family_index(padded, "OG8"), 2, "test_get_family_index_found: a padded id")
    end subroutine test_get_family_index_found

    !> An id not in the list is 0: one that is not there, one that differs in case, a longer one
    !| whose start is an id, and any id in a list of none.
    subroutine test_get_family_index_not_found()
        character(len=8) :: family_ids(2), no_ids(0)

        family_ids = FIXTURE_FAMILY_IDS
        call assert_equal_int(get_family_index(family_ids, "OG3"), 0, "test_get_family_index_not_found: OG3")
        call assert_equal_int(get_family_index(family_ids, "og1"), 0, "test_get_family_index_not_found: og1")
        call assert_equal_int(get_family_index(family_ids, "OG10"), 0, "test_get_family_index_not_found: OG10")
        call assert_equal_int(get_family_index(no_ids, "OG1"), 0, "test_get_family_index_not_found: no id")
    end subroutine test_get_family_index_not_found

    !> Of an id listed twice, the first is found.
    subroutine test_get_family_index_first_of_two()
        character(len=8) :: family_ids(3)

        family_ids = ["OG1", "OG2", "OG1"]
        call assert_equal_int(get_family_index(family_ids, "OG1"), 1, "test_get_family_index_first_of_two")
    end subroutine test_get_family_index_first_of_two
end module mod_test_get_family_index
