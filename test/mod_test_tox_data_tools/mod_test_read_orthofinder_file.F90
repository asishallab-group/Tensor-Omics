!> The `read_orthofinder_file` cases: an OrthoFinder `Orthogroups.tsv` (a header line, then one
!| family per line: its id, then a tab-separated column per species, each a comma-separated list
!| of genes) gives each family its line's index, and each gene of `gene_ids` the index of the
!| family that lists it, or 0 (M_GENE_TO_FAM_SENTINEL) when none does; genes not in `gene_ids`,
!| empty species columns, a file of no family, a missing or empty file, and the shipped files.
!|
!| Each file but the shipped one is written by the case into a scratch directory, and the
!| expected indices are read off its content.
module mod_test_read_orthofinder_file
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_data_tools, only: read_orthofinder_file, read_gene_ids_from_tsv_file
    use tox_errors
    use test_suite, only: test_case
    use scratch_directory
    use mod_test_tox_data_tools_fixtures, only: TAB, LF
    implicit none
    private
    public :: get_all_tests_read_orthofinder_file

    !> The gene_to_fam of a gene in no family.
    integer(int32), parameter :: UNASSIGNED = 0

contains

    !> Every case of `read_orthofinder_file`.
    function get_all_tests_read_orthofinder_file() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(9))
        all_tests(1) = test_case("test_read_orthofinder_file_families", test_read_orthofinder_file_families)
        all_tests(2) = test_case("test_read_orthofinder_file_skips_the_header", test_read_orthofinder_file_skips_the_header)
        all_tests(3) = test_case("test_read_orthofinder_file_empty_species_columns", &
                                 test_read_orthofinder_file_empty_species_columns)
        all_tests(4) = test_case("test_read_orthofinder_file_no_family", test_read_orthofinder_file_no_family)
        all_tests(5) = test_case("test_read_orthofinder_file_missing_or_empty_file", &
                                 test_read_orthofinder_file_missing_or_empty_file)
        all_tests(6) = test_case("test_read_orthofinder_file_shipped_material", test_read_orthofinder_file_shipped_material)
        all_tests(7) = test_case("test_read_orthofinder_file_reads_a_long_line_whole", &
                                 test_read_orthofinder_file_reads_a_long_line_whole)
        all_tests(8) = test_case("test_read_orthofinder_file_skips_a_blank_line", &
                                 test_read_orthofinder_file_skips_a_blank_line)
        all_tests(9) = test_case("test_read_orthofinder_file_checks_the_sizes", &
                                 test_read_orthofinder_file_checks_the_sizes)
    end function get_all_tests_read_orthofinder_file

    !> Three families over two species. A gene's index is its family's line (OG7 is family 1,
    !| whatever its id says); genes are found in either species column and in any position of
    !| a list, blanks after a comma included; a gene in no family gets 0; a listed gene not in
    !| `gene_ids` (h9) is ignored. `gene_ids` is in another order than the file.
    subroutine test_read_orthofinder_file_families()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: orthogroups
        character(len=4) :: gene_ids(6)
        character(len=8) :: family_ids(3), expected_family_ids(3)
        integer(int32) :: gene_to_fam(6), expected_gene_to_fam(6), ierr

        call create_scratch_directory(scratch, "orthofinder_families")
        orthogroups = scratch_file(scratch, "Orthogroups.tsv")
        call write_scratch_file(orthogroups, "Orthogroup"//TAB//"human"//TAB//"mouse"//LF// &
                                "OG7"//TAB//"g3, g1"//TAB//"m2"//LF// &
                                "OG2"//TAB//"g2"//TAB//"h9, m1"//LF// &
                                "OG5"//TAB//"g5"//TAB//"m3"//LF)
        gene_ids = ["m1", "g1", "g2", "g3", "g4", "m2"]

        call read_orthofinder_file(orthogroups, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_families: ierr")
        expected_family_ids = ["OG7", "OG2", "OG5"]
        call assert_equal_array_char(family_ids, expected_family_ids, 8, 3, "test_read_orthofinder_file_families: family_ids")
        expected_gene_to_fam = [2, 1, 2, 1, UNASSIGNED, 1]
        call assert_equal_array_int(gene_to_fam, expected_gene_to_fam, 6, "test_read_orthofinder_file_families: gene_to_fam")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_families

    !> The first line is the header and is skipped, even when it reads like a family.
    subroutine test_read_orthofinder_file_skips_the_header()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: orthogroups
        character(len=4) :: gene_ids(2)
        character(len=8) :: family_ids(1), expected_family_ids(1)
        integer(int32) :: gene_to_fam(2), expected_gene_to_fam(2), ierr

        call create_scratch_directory(scratch, "orthofinder_header")
        orthogroups = scratch_file(scratch, "Orthogroups.tsv")
        call write_scratch_file(orthogroups, "OG1"//TAB//"g1"//LF//"OG2"//TAB//"g2"//LF)
        gene_ids = ["g1", "g2"]

        call read_orthofinder_file(orthogroups, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_skips_the_header: ierr")
        expected_family_ids = ["OG2"]
        call assert_equal_array_char(family_ids, expected_family_ids, 8, 1, &
                                     "test_read_orthofinder_file_skips_the_header: family_ids")
        expected_gene_to_fam = [UNASSIGNED, 1]
        call assert_equal_array_int(gene_to_fam, expected_gene_to_fam, 2, &
                                    "test_read_orthofinder_file_skips_the_header: gene_to_fam")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_skips_the_header

    !> A family with no gene of some species has an empty column there, as OrthoFinder writes
    !| it: two tabs in a row, or a tab at the end of the line. The genes of the other columns
    !| are still found.
    subroutine test_read_orthofinder_file_empty_species_columns()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: orthogroups
        character(len=4) :: gene_ids(3)
        character(len=8) :: family_ids(2)
        integer(int32) :: gene_to_fam(3), expected_gene_to_fam(3), ierr

        call create_scratch_directory(scratch, "orthofinder_empty_columns")
        orthogroups = scratch_file(scratch, "Orthogroups.tsv")
        call write_scratch_file(orthogroups, "Orthogroup"//TAB//"a"//TAB//"b"//TAB//"c"//LF// &
                                "OG1"//TAB//TAB//"g2"//TAB//"g3"//LF// &
                                "OG2"//TAB//"g1"//TAB//TAB//LF)
        gene_ids = ["g1", "g2", "g3"]

        call read_orthofinder_file(orthogroups, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_empty_species_columns: ierr")
        expected_gene_to_fam = [2, 1, 1]
        call assert_equal_array_int(gene_to_fam, expected_gene_to_fam, 3, "test_read_orthofinder_file_empty_species_columns")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_empty_species_columns

    !> A file of nothing but its header lists no family: every gene gets 0, without an error,
    !| also with no gene and no family at all.
    subroutine test_read_orthofinder_file_no_family()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: orthogroups
        character(len=4) :: gene_ids(2), no_gene_ids(0)
        character(len=8) :: family_ids(1), no_family_ids(0)
        integer(int32) :: gene_to_fam(2), no_gene_to_fam(0), expected_gene_to_fam(2), ierr

        call create_scratch_directory(scratch, "orthofinder_no_family")
        orthogroups = scratch_file(scratch, "Orthogroups.tsv")
        call write_scratch_file(orthogroups, "Orthogroup"//TAB//"a"//LF)
        gene_ids = ["g1", "g2"]

        gene_to_fam = 7
        call read_orthofinder_file(orthogroups, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_no_family: ierr")
        expected_gene_to_fam = UNASSIGNED
        call assert_equal_array_int(gene_to_fam, expected_gene_to_fam, 2, "test_read_orthofinder_file_no_family: gene_to_fam")

        call read_orthofinder_file(orthogroups, no_gene_ids, no_family_ids, no_gene_to_fam, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_no_family: no gene, no family")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_no_family

    !> A file that does not exist is refused with ERR_FILE_OPEN, and an empty one, which lacks
    !| the header, with ERR_READ_DATA.
    subroutine test_read_orthofinder_file_missing_or_empty_file()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: missing, empty
        character(len=4) :: gene_ids(1)
        character(len=8) :: family_ids(1)
        integer(int32) :: gene_to_fam(1), ierr

        call create_scratch_directory(scratch, "orthofinder_missing_or_empty")
        missing = scratch_file(scratch, "missing.tsv")
        empty = scratch_file(scratch, "empty.tsv")
        call write_scratch_file(empty, "")
        gene_ids = ["g1"]

        call read_orthofinder_file(missing, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "test_read_orthofinder_file_missing_or_empty_file: missing")
        call read_orthofinder_file(empty, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_READ_DATA, "test_read_orthofinder_file_missing_or_empty_file: empty")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_missing_or_empty_file

    !> Kept from the old `test_data_accessors` and the R suite's `test_2_read_orthofinder_file`
    !| (the families of the first ten genes): the shipped files, read as the pipeline reads them.
    !| Derived with awk from the files: the first ten genes of the expression file are in
    !| families 48 (OG0000047; genes 1-4), 2146 (OG0002145), 231 (OG0000230) and 128 (OG0000127;
    !| genes 7-10); gene 80623, XP_063120092.1, is the last gene of family 1's line (every line
    !| of Orthogroups.tsv ends in CR LF, so this checks the CR is not part of it); the last family
    !| is 15512, OG0015511; and 81960 of the 88327 genes are in a family.
    subroutine test_read_orthofinder_file_shipped_material()
        integer(int32), parameter :: N_GENES = 88327, N_FAMILIES = 15512
        character(len=32), allocatable :: gene_ids(:), family_ids(:)
        integer(int32), allocatable :: gene_to_fam(:)
        integer(int32) :: expected_first_ten(10), ierr

        allocate (gene_ids(N_GENES), family_ids(N_FAMILIES), gene_to_fam(N_GENES))
        call read_gene_ids_from_tsv_file("material/kallisto_sex_data_no_na.tsv", gene_ids, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_shipped_material: the gene ids")

        call read_orthofinder_file("material/Orthogroups.tsv", gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_shipped_material: ierr")
        expected_first_ten = [48, 48, 48, 48, 2146, 231, 128, 128, 128, 128]
        call assert_equal_array_int(gene_to_fam(1:10), expected_first_ten, 10, &
                                    "test_read_orthofinder_file_shipped_material: the first ten genes")
        call assert_string_equal(family_ids(48), "OG0000047", "test_read_orthofinder_file_shipped_material: family 48")
        call assert_string_equal(family_ids(N_FAMILIES), "OG0015511", "test_read_orthofinder_file_shipped_material: last")
        call assert_string_equal(gene_ids(80623), "XP_063120092.1", "test_read_orthofinder_file_shipped_material: gene 80623")
        call assert_equal_int(gene_to_fam(80623), 1, "test_read_orthofinder_file_shipped_material: the last gene of a line")
        call assert_equal_int(count(gene_to_fam /= UNASSIGNED), 81960, &
                              "test_read_orthofinder_file_shipped_material: genes in a family")
    end subroutine test_read_orthofinder_file_shipped_material

    !> BUG (fails today): a line is read into a buffer of 4096 characters, and what is past it is
    !| dropped without a word, so a family listing many genes loses its last ones: of 600 genes
    !| of 8 characters each ("gene0001, ..."), the 191 past character 4096 are left in no family.
    !| OrthoFinder writes such lines for large families over many species. The reader has to
    !| read the whole line or refuse it; leaving genes unassigned with ERR_OK is wrong.
    subroutine test_read_orthofinder_file_reads_a_long_line_whole()
        integer(int32), parameter :: N_GENES = 600
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: orthogroups, family_line
        character(len=8) :: gene_ids(N_GENES), family_ids(1)
        integer(int32) :: gene_to_fam(N_GENES), i_gene, ierr

        call create_scratch_directory(scratch, "orthofinder_long_line")
        family_line = "OG1"//TAB
        do i_gene = 1, N_GENES
            write (gene_ids(i_gene), '("gene", i4.4)') i_gene
            family_line = family_line//gene_ids(i_gene)//", "
        end do
        orthogroups = scratch_file(scratch, "Orthogroups.tsv")
        call write_scratch_file(orthogroups, "Orthogroup"//TAB//"a"//LF//family_line//LF)

        call read_orthofinder_file(orthogroups, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_false(is_ok(ierr) .and. count(gene_to_fam == 1) /= N_GENES, &
                          "test_read_orthofinder_file_reads_a_long_line_whole: genes past character 4096 left in no family")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_reads_a_long_line_whole

    !> BUG (fails today): a line with fewer than two fields is skipped, but counted as a family
    !| first, so a blank line between OG1 and OG2 leaves family_ids(2) unwritten and puts OG2's
    !| genes in family 3, with ERR_OK. (family_ids has room for three here; with exactly two the
    !| same file writes past its end and aborts a bounds-checked run, which is evidence of this
    !| batch, not a case.) The reader has to refuse the line, or skip it without counting it.
    subroutine test_read_orthofinder_file_skips_a_blank_line()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: orthogroups
        character(len=4) :: gene_ids(2)
        character(len=8) :: family_ids(3)
        integer(int32) :: gene_to_fam(2), ierr

        call create_scratch_directory(scratch, "orthofinder_blank_line")
        orthogroups = scratch_file(scratch, "Orthogroups.tsv")
        call write_scratch_file(orthogroups, "Orthogroup"//TAB//"sp"//LF//"OG1"//TAB//"g1"//LF//LF//"OG2"//TAB//"g2"//LF)
        gene_ids = ["g1", "g2"]
        family_ids = "UNSET"

        call read_orthofinder_file(orthogroups, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_true(is_err(ierr) .or. (gene_to_fam(2) == 2 .and. family_ids(2) == "OG2"), &
                         "test_read_orthofinder_file_skips_a_blank_line: g2 in family "//gene_to_fam(2)// &
                         ", family 2 '"//trim(family_ids(2))//"', with ERR_OK")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_skips_a_blank_line

    !> `family_ids` needs room for every family of the file, and `gene_to_fam` an entry per gene:
    !| two families for one family id, and a gene_to_fam one shorter or one longer than gene_ids,
    !| are refused with ERR_SIZE_MISMATCH naming the argument (family_ids is argument 3,
    !| gene_to_fam argument 4). They used to write past the arrays' ends.
    subroutine test_read_orthofinder_file_checks_the_sizes()
        type(scratch_directory_type) :: scratch
        character(len=:), allocatable :: orthogroups
        character(len=4) :: gene_ids(2)
        character(len=8) :: one_family_id(1), family_ids(2)
        integer(int32) :: gene_to_fam(2), one_short(1), one_long(3), ierr

        call create_scratch_directory(scratch, "orthofinder_sizes")
        orthogroups = scratch_file(scratch, "Orthogroups.tsv")
        call write_scratch_file(orthogroups, "Orthogroup"//TAB//"sp"//LF//"OG1"//TAB//"g1"//LF//"OG2"//TAB//"g2"//LF)
        gene_ids = ["g1", "g2"]

        call read_orthofinder_file(orthogroups, gene_ids, one_family_id, gene_to_fam, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_orthofinder_file_checks_the_sizes: two families, one id", 3)
        call read_orthofinder_file(orthogroups, gene_ids, family_ids, one_short, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_orthofinder_file_checks_the_sizes: gene_to_fam one short", 4)
        call read_orthofinder_file(orthogroups, gene_ids, family_ids, one_long, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_read_orthofinder_file_checks_the_sizes: gene_to_fam one long", 4)
        call read_orthofinder_file(orthogroups, gene_ids, family_ids, gene_to_fam, ierr)
        call assert_err(ierr, ERR_OK, "test_read_orthofinder_file_checks_the_sizes: as many as needed")

        call remove_scratch_directory(scratch)
    end subroutine test_read_orthofinder_file_checks_the_sizes
end module mod_test_read_orthofinder_file
