!> Unit tests for save_flyer_json, the flyer JSON export in tox_data_flyer_json.
!|
!| The reference data set (3 axes, 3 families of which one is empty, 5 genes of which one has
!| no family) uses only reals whose 17-digit form is exact, so its file is compared byte for
!| byte with the hand-written test/test_files/flyer_expected.json on every compiler. Every
!| refused input is checked for its error code, its argument position, and that no file was
!| created.
module mod_test_flyer_json
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use asserts
    use test_suite
    use tox_errors, only: ERR_OK, ERR_EMPTY_INPUT, ERR_INVALID_INPUT, ERR_NAN_INF, ERR_INVALID_UTF8, ERR_FILE_OPEN
    use tox_data_flyer_json, only: save_flyer_json
    use mod_test_json_serialize, only: read_whole_file, remove_file, assert_file_equals, assert_no_file
    implicit none
    private
    public :: get_all_tests_flyer_json

    character(len=*), parameter :: LF = achar(10)
    character(len=*), parameter :: EXPECTED_FILE = "test/test_files/flyer_expected.json"
    integer(int32), parameter :: N_AXES = 3, N_GENES = 5, N_FAMILIES = 3

contains

    function get_all_tests_flyer_json() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(19))
        all_tests(1) = test_case("test_flyer_matches_expected_file", test_flyer_matches_expected_file)
        all_tests(2) = test_case("test_flyer_structure", test_flyer_structure)
        all_tests(3) = test_case("test_flyer_zero_families", test_flyer_zero_families)
        all_tests(4) = test_case("test_flyer_zero_genes", test_flyer_zero_genes)
        all_tests(5) = test_case("test_flyer_gene_indices_ascending", test_flyer_gene_indices_ascending)
        all_tests(6) = test_case("test_flyer_empty_species_and_type", test_flyer_empty_species_and_type)
        all_tests(7) = test_case("test_flyer_zero_axes", test_flyer_zero_axes)
        all_tests(8) = test_case("test_flyer_negative_counts", test_flyer_negative_counts)
        all_tests(9) = test_case("test_flyer_non_finite_expression", test_flyer_non_finite_expression)
        all_tests(10) = test_case("test_flyer_non_finite_centroid", test_flyer_non_finite_centroid)
        all_tests(11) = test_case("test_flyer_gene_to_fam_out_of_range", test_flyer_gene_to_fam_out_of_range)
        all_tests(12) = test_case("test_flyer_duplicate_axis_labels", test_flyer_duplicate_axis_labels)
        all_tests(13) = test_case("test_flyer_duplicate_family_ids", test_flyer_duplicate_family_ids)
        all_tests(14) = test_case("test_flyer_duplicate_gene_ids", test_flyer_duplicate_gene_ids)
        all_tests(15) = test_case("test_flyer_empty_ids", test_flyer_empty_ids)
        all_tests(16) = test_case("test_flyer_invalid_utf8", test_flyer_invalid_utf8)
        all_tests(17) = test_case("test_flyer_existing_file", test_flyer_existing_file)
        all_tests(18) = test_case("test_flyer_unwritable_path", test_flyer_unwritable_path)
        all_tests(19) = test_case("test_flyer_first_error_wins", test_flyer_first_error_wins)
    end function get_all_tests_flyer_json

    !> The reference data set, identical to the one in the Python and R suites.
    subroutine reference_data(expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, &
                              family_ids, gene_ids, gene_species, gene_types)
        real(real64), intent(out) :: expression_vectors(N_AXES, N_GENES), family_centroids(N_AXES, N_FAMILIES)
        integer(int32), intent(out) :: gene_to_fam(N_GENES)
        logical(c_bool), intent(out) :: is_outlier(N_GENES)
        character(len=8), intent(out) :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        character(len=12), intent(out) :: gene_species(N_GENES), gene_types(N_GENES)

        expression_vectors(:, 1) = [1.0_real64, 0.5_real64, -2.0_real64]
        expression_vectors(:, 2) = [0.0_real64, -0.0_real64, 0.25_real64]
        expression_vectors(:, 3) = [1024.0_real64, -0.125_real64, 3.0_real64]
        expression_vectors(:, 4) = [1.5_real64, 2.5_real64, -3.5_real64]
        expression_vectors(:, 5) = [0.75_real64, 8.0_real64, -16.0_real64]
        family_centroids(:, 1) = [0.375_real64, 4.25_real64, -8.125_real64]
        family_centroids(:, 2) = [512.5_real64, 0.1875_real64, 0.5_real64]
        family_centroids(:, 3) = [-0.0_real64, 0.0_real64, 0.0625_real64]
        gene_to_fam = [2, 1, 2, 0, 1]
        is_outlier = [.false._c_bool, .true._c_bool, .false._c_bool, .false._c_bool, .true._c_bool]
        axis_labels = [character(len=8) :: "liver", "brain", "heart"]
        family_ids = [character(len=8) :: "F1", "F2", "F3"]
        ! g-grave-e-ne, written as its UTF-8 bytes
        gene_ids = [character(len=8) :: "g1", "g2", "g"//char(195)//char(168)//"ne", "g4", "g5"]
        gene_species = [character(len=12) :: "human", "mouse", "human", "fly", 'say "hi"']
        gene_types = [character(len=12) :: "ortholog", "paralog", "ortholog", "ortholog", "a"//achar(92)//"b"]
    end subroutine reference_data

    !> Write the reference data set, with whatever a test changed in it.
    subroutine write_reference(filename, ierr, expression_vectors, family_centroids, gene_to_fam, axis_labels, &
                               family_ids, gene_ids, gene_species, gene_types)
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        real(real64), intent(in), optional :: expression_vectors(N_AXES, N_GENES), family_centroids(N_AXES, N_FAMILIES)
        integer(int32), intent(in), optional :: gene_to_fam(N_GENES)
        character(len=8), intent(in), optional :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        character(len=12), intent(in), optional :: gene_species(N_GENES), gene_types(N_GENES)
        real(real64) :: ref_expression(N_AXES, N_GENES), ref_centroids(N_AXES, N_FAMILIES)
        integer(int32) :: ref_gene_to_fam(N_GENES)
        logical(c_bool) :: ref_is_outlier(N_GENES)
        character(len=8) :: ref_axis_labels(N_AXES), ref_family_ids(N_FAMILIES), ref_gene_ids(N_GENES)
        character(len=12) :: ref_species(N_GENES), ref_types(N_GENES)

        call reference_data(ref_expression, ref_centroids, ref_gene_to_fam, ref_is_outlier, ref_axis_labels, &
                            ref_family_ids, ref_gene_ids, ref_species, ref_types)
        if (present(expression_vectors)) ref_expression = expression_vectors
        if (present(family_centroids)) ref_centroids = family_centroids
        if (present(gene_to_fam)) ref_gene_to_fam = gene_to_fam
        if (present(axis_labels)) ref_axis_labels = axis_labels
        if (present(family_ids)) ref_family_ids = family_ids
        if (present(gene_ids)) ref_gene_ids = gene_ids
        if (present(gene_species)) ref_species = gene_species
        if (present(gene_types)) ref_types = gene_types
        call remove_file(filename)
        call save_flyer_json(N_AXES, N_GENES, N_FAMILIES, ref_expression, ref_centroids, ref_gene_to_fam, ref_is_outlier, &
                             ref_axis_labels, ref_family_ids, ref_gene_ids, ref_species, ref_types, filename, ierr)
    end subroutine write_reference

    !> The reference data set, byte for byte as the hand-written file.
    subroutine test_flyer_matches_expected_file()
        character(len=*), parameter :: filename = "flyer_reference.test.json"
        character(len=:), allocatable :: expected
        logical :: found
        integer(int32) :: ierr

        call write_reference(filename, ierr)
        call assert_err(ierr, ERR_OK, "reference data set")
        call read_whole_file(EXPECTED_FILE, expected, found)
        call assert_true(found, "the expected file exists")
        call assert_file_equals(filename, expected, "reference data set against "//EXPECTED_FILE)
        call remove_file(filename)
    end subroutine test_flyer_matches_expected_file

    !> The same file, checked for its structure rather than its bytes.
    subroutine test_flyer_structure()
        character(len=*), parameter :: filename = "flyer_structure.test.json"
        character(len=:), allocatable :: content
        logical :: found
        integer(int32) :: ierr

        call write_reference(filename, ierr)
        call assert_err(ierr, ERR_OK, "reference data set")
        call read_whole_file(filename, content, found)
        call assert_true(found, "file written")
        call assert_true(index(content, '{"kind":"tox-flyer","version":1,"tissues":["liver","brain","heart"],') == 1, &
                         "starts with kind, version and tissues")
        call assert_true(content(len(content) - 2:) == "]}"//LF, "ends with the genes and one line feed")
        call assert_equal_int(count_occurrences(content, LF), 1, "a single line")
        call assert_string_contains(content, '{"family":"F1","gene_indices":[2,5],', "F1 holds genes 2 and 5")
        call assert_string_contains(content, '{"family":"F2","gene_indices":[1,3],', "F2 holds genes 1 and 3")
        call assert_string_contains(content, '{"family":"F3","gene_indices":[],', "F3 is empty")
        call assert_equal_int(count_occurrences(content, '"id":'), N_GENES, "one id per gene")
        call assert_equal_int(count_occurrences(content, '"family":null'), 1, "one gene without family")
        call assert_equal_int(count_occurrences(content, '"is_outlier":true'), 2, "two outliers")
        call assert_equal_int(count_occurrences(content, '"is_outlier":false'), 3, "three others")
        call assert_string_contains(content, '"id":"g4","family":null,', "g4 has a null family")
        call remove_file(filename)
    end subroutine test_flyer_structure

    subroutine test_flyer_zero_families()
        character(len=*), parameter :: filename = "flyer_zero_families.test.json"
        real(real64) :: expression_vectors(2, 2), family_centroids(2, 0)
        integer(int32) :: gene_to_fam(2), ierr
        logical(c_bool) :: is_outlier(2)
        character(len=1) :: no_family_ids(0), axis_labels(2), gene_ids(2), gene_species(2), gene_types(2)

        expression_vectors(:, 1) = [1.0_real64, 2.0_real64]
        expression_vectors(:, 2) = [3.0_real64, 4.0_real64]
        gene_to_fam = 0
        is_outlier = [.true._c_bool, .false._c_bool]
        axis_labels = ["x", "y"]
        gene_ids = ["a", "b"]
        gene_species = "s"
        gene_types = "t"
        call remove_file(filename)
        call save_flyer_json(2, 2, 0, expression_vectors, family_centroids, gene_to_fam, is_outlier, &
                             axis_labels, no_family_ids, gene_ids, gene_species, gene_types, filename, ierr)
        call assert_err(ierr, ERR_OK, "zero families")
        call assert_file_equals(filename, '{"kind":"tox-flyer","version":1,"tissues":["x","y"],"families":[],"genes":[' &
                                //'{"coordinates":[1.0000000000000000E+000,2.0000000000000000E+000],"id":"a",' &
                                //'"family":null,"species":"s","is_outlier":true,"type":"t"},' &
                                //'{"coordinates":[3.0000000000000000E+000,4.0000000000000000E+000],"id":"b",' &
                                //'"family":null,"species":"s","is_outlier":false,"type":"t"}]}'//LF, "zero families")
        call remove_file(filename)
    end subroutine test_flyer_zero_families

    subroutine test_flyer_zero_genes()
        character(len=*), parameter :: filename = "flyer_zero_genes.test.json"
        real(real64) :: expression_vectors(1, 0), family_centroids(1, 2)
        integer(int32) :: no_gene_to_fam(0)
        logical(c_bool) :: no_is_outlier(0)
        character(len=4) :: no_strings(0), axis_labels(1), family_ids(2)
        integer(int32) :: ierr

        family_centroids(1, :) = [0.5_real64, -0.0_real64]
        axis_labels = "axis"
        family_ids = ["F1", "F2"]
        call remove_file(filename)
        call save_flyer_json(1, 0, 2, expression_vectors, family_centroids, no_gene_to_fam, no_is_outlier, &
                             axis_labels, family_ids, no_strings, no_strings, no_strings, filename, ierr)
        call assert_err(ierr, ERR_OK, "zero genes")
        call assert_file_equals(filename, '{"kind":"tox-flyer","version":1,"tissues":["axis"],"families":[' &
                                //'{"family":"F1","gene_indices":[],"centroid":[5.0000000000000000E-001]},' &
                                //'{"family":"F2","gene_indices":[],"centroid":[-0.0]}],"genes":[]}'//LF, "zero genes")
        call remove_file(filename)
    end subroutine test_flyer_zero_genes

    !> Families interleaved over many genes: each lists exactly its genes, ascending.
    subroutine test_flyer_gene_indices_ascending()
        character(len=*), parameter :: filename = "flyer_gene_indices.test.json"
        integer(int32), parameter :: n_genes = 50, n_families = 4
        real(real64) :: expression_vectors(1, n_genes), family_centroids(1, n_families)
        integer(int32) :: gene_to_fam(n_genes)
        logical(c_bool) :: is_outlier(n_genes)
        character(len=4) :: family_ids(n_families), gene_ids(n_genes), labels(n_genes), axis_labels(1)
        character(len=:), allocatable :: content, expected_members
        character(len=12) :: number
        logical :: found
        integer(int32) :: ierr, i_gene, i_family

        expression_vectors = 1.0_real64
        family_centroids = 2.0_real64
        is_outlier = .false._c_bool
        labels = "x"
        axis_labels = "axis"
        do i_gene = 1, n_genes
            gene_to_fam(i_gene) = mod(7*i_gene, n_families + 1)
            write (gene_ids(i_gene), "(A,I0)") "g", i_gene
        end do
        do i_family = 1, n_families
            write (family_ids(i_family), "(A,I0)") "F", i_family
        end do
        call remove_file(filename)
        call save_flyer_json(1, n_genes, n_families, expression_vectors, family_centroids, gene_to_fam, is_outlier, &
                             axis_labels, family_ids, gene_ids, labels, labels, filename, ierr)
        call assert_err(ierr, ERR_OK, "interleaved families")
        call read_whole_file(filename, content, found)
        do i_family = 1, n_families
            expected_members = ""
            do i_gene = 1, n_genes
                if (gene_to_fam(i_gene) /= i_family) cycle
                write (number, "(I0)") i_gene
                if (len(expected_members) > 0) expected_members = expected_members//","
                expected_members = expected_members//trim(number)
            end do
            call assert_string_contains(content, '"family":"'//trim(family_ids(i_family))//'","gene_indices":[' &
                                        //expected_members//"]", "members of family "//i_family)
        end do
        call assert_equal_int(count_occurrences(content, '"family":null'), count(gene_to_fam == 0), &
                              "genes without family")
        call remove_file(filename)
    end subroutine test_flyer_gene_indices_ascending

    !> Species and type are not identifiers: empty ones are written as "".
    subroutine test_flyer_empty_species_and_type()
        character(len=*), parameter :: filename = "flyer_empty_species.test.json"
        character(len=12) :: gene_species(N_GENES), gene_types(N_GENES)
        character(len=:), allocatable :: content
        logical :: found
        integer(int32) :: ierr

        gene_species = ""
        gene_types = ""
        call write_reference(filename, ierr, gene_species=gene_species, gene_types=gene_types)
        call assert_err(ierr, ERR_OK, "empty species and types are allowed")
        call read_whole_file(filename, content, found)
        call assert_equal_int(count_occurrences(content, '"species":"","is_outlier"'), N_GENES, "empty species")
        call assert_equal_int(count_occurrences(content, '"type":""}'), N_GENES, "empty types")
        call remove_file(filename)
    end subroutine test_flyer_empty_species_and_type

    ! ============================================================================================
    ! Refused inputs: code, position, and no file
    ! ============================================================================================

    subroutine test_flyer_zero_axes()
        character(len=*), parameter :: filename = "flyer_zero_axes.test.json"
        real(real64) :: expression_vectors(0, 1), family_centroids(0, 1)
        integer(int32) :: gene_to_fam(1), ierr
        logical(c_bool) :: is_outlier(1)
        character(len=1) :: no_labels(0), family_ids(1), gene_ids(1), gene_species(1), gene_types(1)

        gene_to_fam = 1
        is_outlier = .false._c_bool
        family_ids = "F"
        gene_ids = "g"
        gene_species = "s"
        gene_types = "t"
        call remove_file(filename)
        call save_flyer_json(0, 1, 1, expression_vectors, family_centroids, gene_to_fam, is_outlier, no_labels, &
                             family_ids, gene_ids, gene_species, gene_types, filename, ierr)
        call assert_refused(ierr, ERR_EMPTY_INPUT, 1, filename, "zero axes")
    end subroutine test_flyer_zero_axes

    subroutine test_flyer_negative_counts()
        character(len=*), parameter :: filename = "flyer_negative.test.json"
        real(real64) :: expression_vectors(1, 1), family_centroids(1, 1)
        integer(int32) :: gene_to_fam(1), ierr
        logical(c_bool) :: is_outlier(1)
        character(len=1) :: axis_labels(1), family_ids(1), gene_ids(1), gene_species(1), gene_types(1)

        expression_vectors = 1.0_real64
        family_centroids = 1.0_real64
        gene_to_fam = 0
        is_outlier = .false._c_bool
        axis_labels = "a"
        family_ids = "F"
        gene_ids = "g"
        gene_species = "s"
        gene_types = "t"
        call remove_file(filename)
        call save_flyer_json(1, -1, 1, expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, &
                             family_ids, gene_ids, gene_species, gene_types, filename, ierr)
        call assert_refused(ierr, ERR_INVALID_INPUT, 2, filename, "negative gene count")
        call save_flyer_json(1, 1, -1, expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, &
                             family_ids, gene_ids, gene_species, gene_types, filename, ierr)
        call assert_refused(ierr, ERR_INVALID_INPUT, 3, filename, "negative family count")
    end subroutine test_flyer_negative_counts

    subroutine test_flyer_non_finite_expression()
        character(len=*), parameter :: filename = "flyer_nan_expression.test.json"
        real(real64) :: expression_vectors(N_AXES, N_GENES), family_centroids(N_AXES, N_FAMILIES)
        integer(int32) :: gene_to_fam(N_GENES), ierr
        logical(c_bool) :: is_outlier(N_GENES)
        character(len=8) :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        character(len=12) :: gene_species(N_GENES), gene_types(N_GENES)

        call reference_data(expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, family_ids, &
                            gene_ids, gene_species, gene_types)
        expression_vectors(2, 3) = ieee_value(1.0_real64, ieee_quiet_nan)
        call write_reference(filename, ierr, expression_vectors=expression_vectors)
        call assert_refused(ierr, ERR_NAN_INF, 4, filename, "NaN in expression_vectors")
        expression_vectors(2, 3) = ieee_value(1.0_real64, ieee_positive_inf)
        call write_reference(filename, ierr, expression_vectors=expression_vectors)
        call assert_refused(ierr, ERR_NAN_INF, 4, filename, "+Inf in expression_vectors")
    end subroutine test_flyer_non_finite_expression

    subroutine test_flyer_non_finite_centroid()
        character(len=*), parameter :: filename = "flyer_inf_centroid.test.json"
        real(real64) :: expression_vectors(N_AXES, N_GENES), family_centroids(N_AXES, N_FAMILIES)
        integer(int32) :: gene_to_fam(N_GENES), ierr
        logical(c_bool) :: is_outlier(N_GENES)
        character(len=8) :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        character(len=12) :: gene_species(N_GENES), gene_types(N_GENES)

        call reference_data(expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, family_ids, &
                            gene_ids, gene_species, gene_types)
        ! the empty family's centroid is checked like any other
        family_centroids(3, 3) = ieee_value(1.0_real64, ieee_negative_inf)
        call write_reference(filename, ierr, family_centroids=family_centroids)
        call assert_refused(ierr, ERR_NAN_INF, 5, filename, "-Inf in family_centroids")
    end subroutine test_flyer_non_finite_centroid

    subroutine test_flyer_gene_to_fam_out_of_range()
        character(len=*), parameter :: filename = "flyer_gene_to_fam.test.json"
        integer(int32) :: gene_to_fam(N_GENES), ierr

        gene_to_fam = [2, 1, -1, 0, 1]
        call write_reference(filename, ierr, gene_to_fam=gene_to_fam)
        call assert_refused(ierr, ERR_INVALID_INPUT, 6, filename, "gene_to_fam below 0")
        gene_to_fam = [2, 1, 2, 0, N_FAMILIES + 1]
        call write_reference(filename, ierr, gene_to_fam=gene_to_fam)
        call assert_refused(ierr, ERR_INVALID_INPUT, 6, filename, "gene_to_fam above n_families")
    end subroutine test_flyer_gene_to_fam_out_of_range

    subroutine test_flyer_duplicate_axis_labels()
        character(len=*), parameter :: filename = "flyer_dup_axes.test.json"
        character(len=8) :: axis_labels(N_AXES)
        integer(int32) :: ierr

        axis_labels = [character(len=8) :: "liver", "brain", "liver"]
        call write_reference(filename, ierr, axis_labels=axis_labels)
        call assert_refused(ierr, ERR_INVALID_INPUT, 8, filename, "repeated axis label")
    end subroutine test_flyer_duplicate_axis_labels

    subroutine test_flyer_duplicate_family_ids()
        character(len=*), parameter :: filename = "flyer_dup_families.test.json"
        character(len=8) :: family_ids(N_FAMILIES)
        integer(int32) :: ierr

        family_ids = [character(len=8) :: "F2", "F1", "F2"]
        call write_reference(filename, ierr, family_ids=family_ids)
        call assert_refused(ierr, ERR_INVALID_INPUT, 9, filename, "repeated family id")
    end subroutine test_flyer_duplicate_family_ids

    subroutine test_flyer_duplicate_gene_ids()
        character(len=*), parameter :: filename = "flyer_dup_genes.test.json"
        character(len=8) :: gene_ids(N_GENES)
        integer(int32) :: ierr

        gene_ids = [character(len=8) :: "g1", "g2", "g3", "g4", "g2"]
        call write_reference(filename, ierr, gene_ids=gene_ids)
        call assert_refused(ierr, ERR_INVALID_INPUT, 10, filename, "repeated gene id")
    end subroutine test_flyer_duplicate_gene_ids

    subroutine test_flyer_empty_ids()
        character(len=*), parameter :: filename = "flyer_empty_ids.test.json"
        character(len=8) :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        integer(int32) :: ierr

        axis_labels = [character(len=8) :: "liver", "  ", "heart"]
        call write_reference(filename, ierr, axis_labels=axis_labels)
        call assert_refused(ierr, ERR_INVALID_INPUT, 8, filename, "empty axis label")
        family_ids = [character(len=8) :: "F1", "F2", ""]
        call write_reference(filename, ierr, family_ids=family_ids)
        call assert_refused(ierr, ERR_INVALID_INPUT, 9, filename, "empty family id")
        gene_ids = [character(len=8) :: "", "g2", "g3", "g4", "g5"]
        call write_reference(filename, ierr, gene_ids=gene_ids)
        call assert_refused(ierr, ERR_INVALID_INPUT, 10, filename, "empty gene id")
    end subroutine test_flyer_empty_ids

    subroutine test_flyer_invalid_utf8()
        character(len=*), parameter :: filename = "flyer_utf8.test.json"
        character(len=*), parameter :: BAD = "x"//char(237)//char(160)//char(128)   ! a surrogate
        character(len=8) :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        character(len=12) :: gene_species(N_GENES), gene_types(N_GENES)
        integer(int32) :: ierr

        axis_labels = [character(len=8) :: "liver", BAD, "heart"]
        call write_reference(filename, ierr, axis_labels=axis_labels)
        call assert_refused(ierr, ERR_INVALID_UTF8, 8, filename, "invalid UTF-8 in axis_labels")
        family_ids = [character(len=8) :: "F1", "F2", BAD]
        call write_reference(filename, ierr, family_ids=family_ids)
        call assert_refused(ierr, ERR_INVALID_UTF8, 9, filename, "invalid UTF-8 in family_ids")
        gene_ids = [character(len=8) :: "g1", "g2", "g3", BAD, "g5"]
        call write_reference(filename, ierr, gene_ids=gene_ids)
        call assert_refused(ierr, ERR_INVALID_UTF8, 10, filename, "invalid UTF-8 in gene_ids")
        gene_species = [character(len=12) :: "a", "b", "c", "d", "e"//char(200)]
        call write_reference(filename, ierr, gene_species=gene_species)
        call assert_refused(ierr, ERR_INVALID_UTF8, 11, filename, "invalid UTF-8 in gene_species")
        gene_types = [character(len=12) :: char(128), "b", "c", "d", "e"]
        call write_reference(filename, ierr, gene_types=gene_types)
        call assert_refused(ierr, ERR_INVALID_UTF8, 12, filename, "invalid UTF-8 in gene_types")
    end subroutine test_flyer_invalid_utf8

    !> An existing file is never overwritten.
    subroutine test_flyer_existing_file()
        character(len=*), parameter :: filename = "flyer_existing.test.json"
        integer(int32) :: ierr, unit

        call remove_file(filename)
        open (newunit=unit, file=filename, access="stream", form="unformatted", status="new", action="write")
        write (unit) "keep me"
        close (unit)
        call write_existing(filename, ierr)
        call assert_err(ierr, ERR_FILE_OPEN, "existing file", 13)
        call assert_file_equals(filename, "keep me", "the existing file is untouched")
        call remove_file(filename)
    end subroutine test_flyer_existing_file

    subroutine test_flyer_unwritable_path()
        character(len=*), parameter :: filename = "flyer_no_such_directory.test/out.json"
        integer(int32) :: ierr

        call write_reference(filename, ierr)
        call assert_refused(ierr, ERR_FILE_OPEN, 13, filename, "path in a missing directory")
    end subroutine test_flyer_unwritable_path

    !> Checks run in argument order, and the first failure is the one reported.
    subroutine test_flyer_first_error_wins()
        character(len=*), parameter :: filename = "flyer_first_error.test.json"
        real(real64) :: expression_vectors(N_AXES, N_GENES), family_centroids(N_AXES, N_FAMILIES)
        integer(int32) :: gene_to_fam(N_GENES), ierr
        logical(c_bool) :: is_outlier(N_GENES)
        character(len=8) :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        character(len=12) :: gene_species(N_GENES), gene_types(N_GENES)

        call reference_data(expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, family_ids, &
                            gene_ids, gene_species, gene_types)
        family_centroids(1, 1) = ieee_value(1.0_real64, ieee_quiet_nan)
        gene_ids(2) = gene_ids(1)
        call write_reference(filename, ierr, family_centroids=family_centroids, gene_ids=gene_ids)
        call assert_refused(ierr, ERR_NAN_INF, 5, filename, "NaN before a duplicate id")
    end subroutine test_flyer_first_error_wins

    ! ============================================================================================
    ! Helpers
    ! ============================================================================================

    !> Write the reference data set without removing an existing file first.
    subroutine write_existing(filename, ierr)
        character(len=*), intent(in) :: filename
        integer(int32), intent(out) :: ierr
        real(real64) :: expression_vectors(N_AXES, N_GENES), family_centroids(N_AXES, N_FAMILIES)
        integer(int32) :: gene_to_fam(N_GENES)
        logical(c_bool) :: is_outlier(N_GENES)
        character(len=8) :: axis_labels(N_AXES), family_ids(N_FAMILIES), gene_ids(N_GENES)
        character(len=12) :: gene_species(N_GENES), gene_types(N_GENES)

        call reference_data(expression_vectors, family_centroids, gene_to_fam, is_outlier, axis_labels, family_ids, &
                            gene_ids, gene_species, gene_types)
        call save_flyer_json(N_AXES, N_GENES, N_FAMILIES, expression_vectors, family_centroids, gene_to_fam, is_outlier, &
                             axis_labels, family_ids, gene_ids, gene_species, gene_types, filename, ierr)
    end subroutine write_existing

    !> A refused input: its code, its argument position, and no file.
    subroutine assert_refused(ierr, expected_code, arg_pos, filename, msg)
        integer(int32), intent(in) :: ierr, expected_code, arg_pos
        character(len=*), intent(in) :: filename, msg

        call assert_err(ierr, expected_code, msg, arg_pos)
        call assert_no_file(filename, msg)
    end subroutine assert_refused

    !> How often `pattern` occurs in `text`, without overlaps.
    integer(int32) function count_occurrences(text, pattern)
        character(len=*), intent(in) :: text, pattern
        integer :: search_start, found_at

        count_occurrences = 0
        search_start = 1
        do
            found_at = index(text(search_start:), pattern)
            if (found_at == 0) exit
            count_occurrences = count_occurrences + 1
            search_start = search_start + found_at + len(pattern) - 1
            if (search_start > len(text)) exit
        end do
    end function count_occurrences

end module mod_test_flyer_json
