!> Unit test suite for f42_json routine.
module mod_test_json
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use f42_json
    use tox_errors
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = epsilon(1.0_real64)

    character(len=119), parameter :: TEST_STRING = 'we test "quoting", /slashing and \backslashing, 	tabbing, ' //&
                                                    achar(8) //  'backspacing, ' //&
                                                    achar(10) // 'new-lining, ' //&
                                                    achar(13) // 'carriage-returning, ' //&
                                                    achar(12) // 'form-feeding'

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(2)

        all_tests(1) = test_case("test_f42_json_serialization", test_serialization)
        all_tests(2) = test_case("test_f42_flyer_serialization", test_flyer_serialization)
    end function get_all_tests

    subroutine test_serialization
        ! Test idea:
        !   - an array of all possible types
        !   - for simplicity the array type is a self reference and tests recursion simultaneously
        !   - also, the self-containing array is reused as values of the object, and the object is also included in the array
        ! 
        type(json_value), dimension(6), target :: elements
        type(json_array), target :: array
        type(json_object), target :: object
        character(len=:), dimension(:), allocatable, target :: keys
        integer(int32) :: i_element
        integer(int32), target:: test_integer
        real(real64), target:: test_real
        complex(real64), target:: test_complex
        logical, target:: test_logical

        allocate(character(len=7) :: keys(6))

        keys(1) = "integer"
        test_integer = -huge(1_int32)
        elements(1)%value => test_integer

        keys(2) = "real"
        test_real = huge(1.0_real64)
        elements(2)%value => test_real

        keys(3) = "logical"
        test_logical = .true.
        elements(3)%value => test_logical

        keys(4) = "complex"
        test_complex = cmplx(1.0_real64, -1.0_real64)
        elements(4)%value => test_complex

        keys(5) = "array"
        array%elements => elements
        elements(5)%value => array

        keys(6) = "object"
        object%keys => keys
        object%values => elements
        elements(6)%value => object

        ! Case 1: basic object serialization test with recursion
        call helper_test_serialization(&
            object,&
            '{"integer":-2147483647,"real": 1.7976931348623157E+308,"logical":true,"complex":[ 1.0000000000000000E+000,-1.0000000000000000E+000],"array":[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[-2147483647,',&
            "test_serialization: case 1 Object"&
        )

        ! Case 2: basic array serialization test with recursion
        call helper_test_serialization(&
            array,&
            '[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[-2147483647,',&
            "test_serialization: case 2 Array"&
        )

        ! Case 3: basic array serialization test with custom recursion limit
        call helper_test_serialization(&
            array,&
            '[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[null,null,null,null,null,null],{"integer":null,"real":null,"logical":null,"complex":null,"array":null,"object":null}]',&
            "test_serialization: case 3 custom recursion limit",&
            2_int32&
        )

        ! Case 4: NaN/Inf handling
        test_complex = cmplx(ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf))
        test_real = ieee_value(1.0_real64, ieee_quiet_nan)

        call helper_test_serialization(&
            array,&
            '[-2147483647,null,true,[null,null],[-2147483647,',&
            "test_serialization: case 4 NaN/Inf handling"&
        )

        ! Case 5: Different lengths of json object's key-value arrays
        object%keys => keys(1:4)
        call helper_test_serialization(&
            object,&
            '{"integer":-2147483647,"real":null,"logical":true,"complex":[null,null]}',&
            "test_serialization: case 5 different size key-value arrays"&
        )

        ! Case 6: unassigned array values
        do i_element = 1, 4
            nullify(elements(i_element)%value)
        end do
        nullify(object%keys)
        nullify(elements(5)%value)

        call helper_test_serialization(&
            array,&
            '[null,null,null,null,null,{}]',&
            "test_serialization: case 6 unassigned array values"&
        )

        ! Case 7: unassigned object values
        elements(5)%value => array
        object%keys => keys
        nullify(array%elements)
        nullify(elements(6)%value)

        call helper_test_serialization(&
            object,&
            '{"integer":null,"real":null,"logical":null,"complex":null,"array":[],"object":null}',&
            "test_serialization: case 7 unassigned object values"&
        )

        ! Case 8: strings
        deallocate(keys)
        allocate(character(len=len_trim(TEST_STRING)) :: keys(1))
        keys(1) = TEST_STRING
        elements(1)%value => keys(1)
        object%keys => keys(1:1)
        object%values => keys(1:1)

        call helper_test_serialization(&
            object,&
            '{"we test \"quoting\", \/slashing and \\backslashing, \ttabbing, \bbackspacing, \nnew-lining, \rcarriage-returning, \fform-feeding"' //&
            ':"we test \"quoting\", \/slashing and \\backslashing, \ttabbing, \bbackspacing, \nnew-lining, \rcarriage-returning, \fform-feeding"}',&
            "test_serialization: case 8 strings"&
        )
    end subroutine test_serialization

    subroutine test_flyer_serialization
        integer(int32), parameter :: n_tissues = 2, n_families = 4, n_genes = 5
        character(len=29) :: filename
        character(len=7), dimension(n_tissues) :: tissues
        character(len=9), dimension(n_families) :: family_ids
        character(len=14), dimension(n_genes) :: gene_ids
        character(len=8), dimension(n_genes) :: gene_types
        character(len=21), dimension(n_genes) :: gene_species
        real(real64), dimension(n_tissues, n_families) :: centroids
        real(real64), dimension(n_tissues, n_genes) :: genes
        integer(int32), dimension(n_genes) :: gene_to_fam
        integer(int32), dimension(n_genes) :: sorted_gene_to_fam_perm
        logical, dimension(n_genes) :: gene_outliers
        integer(int32) :: ierr, unit
        logical :: exists

        filename = "test_flyer_serialization.json"
        tissues = ["Adipose", "Thyroid"]
        family_ids = ["EMPTYFAM1", "OG0000000", "OG0000001", "EMPTYFAM2"]
        gene_ids = ["NP_001243379.1", "UNASSIGNED.1  ", "XP_038381480.1", "NP_001243386.1", "XP_038421312.1"]
        gene_species = ["Canis_lupus_protein.1", "Whatever_protein     ", "Canis_lupus_protein.2", "Canis_lupus_protein.3", "Canis_lupus_protein.4"]
        gene_types = ["ortholog", "ortholog", "paralog ", "ortholog", "paralog "]
        gene_to_fam = [2, 0, 3, 2, 3]
        sorted_gene_to_fam_perm = [2, 1, 4, 3, 5]
        gene_outliers = [.true., .false., .false., .false., .true.]
        genes(:, 1) = [0.0541530138142222_real64, 0.0041979991981664_real64]
        genes(:, 2) = [3.1415926535897932_real64, 2.7182818284590452_real64]
        genes(:, 3) = [1.115620102135_real64, 0.308001439923219_real64]
        genes(:, 4) = [0.0060563875402208_real64, 0.147377684857407_real64]
        genes(:, 5) = [0.0060563875402208_real64, 0.0041979991981664_real64]
        centroids(:, 1) = 0.0_real64
        centroids(:, 2) = [0.03010470067722_real64, 0.07578784202779_real64]
        centroids(:, 3) = [0.56083824483761_real64, 0.15609971956069_real64]
        centroids(:, 4) = 0.0_real64

        call serialize_tox_data_as_flyer_json(filename, tissues, n_tissues, family_ids, n_families, centroids, gene_ids, n_genes, genes, gene_to_fam, sorted_gene_to_fam_perm, gene_outliers, gene_species, gene_types, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_flyer_serialization: Unexpected error when serializing")
        inquire(file=filename, exist=exists)
        call assert_true(exists, "test_flyer_serialization: '" // filename // "' does not exist")

        open(newunit=unit, file=filename, status="old", action="read")
        call helper_check_fragment(&
            unit,&
            '{' //&
                '"tissues":["Adipose","Thyroid"],' //&
                '"families":[' //&
                    '{"family":"EMPTYFAM1","gene_indices":[],"centroid":[ 0.0000000000000000E+000, 0.0000000000000000E+000]},' //&
                    '{"family":"OG0000000","gene_indices":[1,4],"centroid":[ 3.0104700677220000E-002, 7.5787842027790001E-002]},' //&
                    '{"family":"OG0000001","gene_indices":[3,5],"centroid":[ 5.6083824483761002E-001, 1.5609971956068999E-001]},' //&
                    '{"family":"EMPTYFAM2","gene_indices":[],"centroid":[ 0.0000000000000000E+000, 0.0000000000000000E+000]}' //&
                '],' //&
                '"genes":[' //&
                    '{"coordinates":[ 5.4153013814222203E-002, 4.1979991981664000E-003],"id":"NP_001243379.1","family":"OG0000000","species":"Canis_lupus_protein.1","is_outlier":true,"type":"ortholog"},' //&
                    '{"coordinates":[ 3.1415926535897931E+000, 2.7182818284590451E+000],"id":"UNASSIGNED.1","family":null,"species":"Whatever_protein","is_outlier":false,"type":"ortholog"},' //&
                    '{"coordinates":[ 1.1156201021350001E+000, 3.0800143992321899E-001],"id":"XP_038381480.1","family":"OG0000001","species":"Canis_lupus_protein.2","is_outlier":false,"type":"paralog"},' //&
                    '{"coordinates":[ 6.0563875402208003E-003, 1.4737768485740699E-001],"id":"NP_001243386.1","family":"OG0000000","species":"Canis_lupus_protein.3","is_outlier":false,"type":"ortholog"},' //&
                    '{"coordinates":[ 6.0563875402208003E-003, 4.1979991981664000E-003],"id":"XP_038421312.1","family":"OG0000001","species":"Canis_lupus_protein.4","is_outlier":true,"type":"paralog"}' //&
                ']' //&
            '}',&
            "test_flyer_serialization"&
        )
    end subroutine test_flyer_serialization

    subroutine helper_test_serialization(json_var, expected_fragment, test_case, max_depth)
        class(*), intent(in) :: json_var
        character(len=*), intent(in) :: expected_fragment
        character(len=*), intent(in) :: test_case
        integer(int32), intent(in), optional :: max_depth

        integer(int32) :: ierr, unit

        open(newunit=unit, file="test_json.json", form='formatted', access='stream', status='replace', iostat=ierr)
        call assert_equal_int(ierr, ERR_OK, trim(test_case) // ": could not open file")

        select type (json_var)
            type is (json_array)
                call serialize_json_array(json_var, unit, max_depth)
            type is (json_object)
                call serialize_json_object(json_var, unit, max_depth)
            class default
                error stop trim(test_case) // ": helper_test_serialization: Unsupported type"
        end select

        rewind(unit)

        call helper_check_fragment(unit, expected_fragment, test_case)
    end subroutine helper_test_serialization

    subroutine helper_check_fragment(unit, expected_fragment, test_case)
        character(len=*), intent(in) :: expected_fragment
        integer(int32), intent(in) :: unit
        character(len=*), intent(in) :: test_case

        character(len=:), allocatable :: actual_fragment

        allocate(actual_fragment, mold=expected_fragment)
        read (unit, "(A)") actual_fragment
        call assert_string_equal(actual_fragment, expected_fragment, trim(test_case) // ": fragments differ")
        close(unit, status="delete")
    end subroutine helper_check_fragment

    !> Run all f42_json tests.
    subroutine run_all_tests_f42_json
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All f42_json tests passed successfully."
    end subroutine run_all_tests_f42_json

    !> Run specific f42_json tests by name.
    subroutine run_named_tests_f42_json(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_f42_json
end module mod_test_json
