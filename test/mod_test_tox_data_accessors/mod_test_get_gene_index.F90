!> The `get_gene_index` cases: the index of the first id of `gene_ids` equal to `gene`, blanks
!| around either not counted, or 0 when there is none; a blank `gene` is never found.
module mod_test_get_gene_index
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_data_accessors, only: get_gene_index
    use test_suite, only: test_case
    use mod_test_tox_data_accessors_fixtures, only: FIXTURE_GENE_IDS
    implicit none
    private
    public :: get_all_tests_get_gene_index

contains

    !> Every case of `get_gene_index`.
    function get_all_tests_get_gene_index() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(4))
        all_tests(1) = test_case("test_get_gene_index_found", test_get_gene_index_found)
        all_tests(2) = test_case("test_get_gene_index_not_found", test_get_gene_index_not_found)
        all_tests(3) = test_case("test_get_gene_index_blank", test_get_gene_index_blank)
        all_tests(4) = test_case("test_get_gene_index_first_of_two", test_get_gene_index_first_of_two)
    end function get_all_tests_get_gene_index

    !> Kept from the old `test_data_accessors` (NP_001000001.1 is gene 2 of the shipped data):
    !| each id gives its position, and blanks around the id, in the list or in the query, do not
    !| count.
    subroutine test_get_gene_index_found()
        character(len=8) :: gene_ids(3), padded(2)

        gene_ids = FIXTURE_GENE_IDS
        call assert_equal_int(get_gene_index(gene_ids, "g1"), 1, "test_get_gene_index_found: g1")
        call assert_equal_int(get_gene_index(gene_ids, "g2"), 2, "test_get_gene_index_found: g2")
        call assert_equal_int(get_gene_index(gene_ids, "g3"), 3, "test_get_gene_index_found: g3")
        call assert_equal_int(get_gene_index(gene_ids, "  g2  "), 2, "test_get_gene_index_found: a padded query")
        padded = ["  a", " b "]
        call assert_equal_int(get_gene_index(padded, "b"), 2, "test_get_gene_index_found: a padded id")
    end subroutine test_get_gene_index_found

    !> An id not in the list is 0: one that is not there, one that differs in case, a longer one
    !| whose start is an id (no prefix match), and any id in a list of none.
    subroutine test_get_gene_index_not_found()
        character(len=8) :: gene_ids(3), no_ids(0)

        gene_ids = FIXTURE_GENE_IDS
        call assert_equal_int(get_gene_index(gene_ids, "g4"), 0, "test_get_gene_index_not_found: g4")
        call assert_equal_int(get_gene_index(gene_ids, "G1"), 0, "test_get_gene_index_not_found: G1")
        call assert_equal_int(get_gene_index(gene_ids, "g10"), 0, "test_get_gene_index_not_found: g10")
        call assert_equal_int(get_gene_index(no_ids, "g1"), 0, "test_get_gene_index_not_found: no id")
    end subroutine test_get_gene_index_not_found

    !> A blank query is never found, not even in a list with a blank id.
    subroutine test_get_gene_index_blank()
        character(len=8) :: gene_ids(2)

        gene_ids = ["g1", "  "]
        call assert_equal_int(get_gene_index(gene_ids, ""), 0, "test_get_gene_index_blank: empty")
        call assert_equal_int(get_gene_index(gene_ids, "   "), 0, "test_get_gene_index_blank: blanks")
    end subroutine test_get_gene_index_blank

    !> Of an id listed twice, the first is found.
    subroutine test_get_gene_index_first_of_two()
        character(len=8) :: gene_ids(3)

        gene_ids = ["g1", "g2", "g2"]
        call assert_equal_int(get_gene_index(gene_ids, "g2"), 2, "test_get_gene_index_first_of_two")
    end subroutine test_get_gene_index_first_of_two
end module mod_test_get_gene_index
