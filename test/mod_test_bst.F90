!> Unit test suite for binary_search_tree module.
module mod_test_bst
    use f42_binary_search_tree
    use tox_errors
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_bst() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(7))
        all_tests(1) = test_case("test_bst_index_construction", test_bst_index_construction)
        all_tests(2) = test_case("test_bst_sorted_values", test_bst_sorted_values)
        all_tests(3) = test_case("test_bst_range_query", test_bst_range_query)
        all_tests(4) = test_case("test_bst_empty_array", test_bst_empty_array)
        all_tests(5) = test_case("test_bst_single_element", test_bst_single_element)
        all_tests(6) = test_case("test_bst_identical_values", test_bst_identical_values)
        all_tests(7) = test_case("test_bst_large_random", test_bst_large_random)
    end function get_all_tests_bst

    !> Test BST index construction and monotonicity.
    subroutine test_bst_index_construction()
        integer(int32), parameter :: n = 100
        real(real64) :: x(n)
        integer(int32) :: ix(n), expected_ix(n)
        integer(int32) :: i, ierr
        logical :: is_sorted

        call set_ok(ierr)

        ! call random_array(x, n)
        x = [ &
            .7553780843536994_real64, &
            .3521157186151643_real64, &
            .8386715726401098_real64, &
            .4632533650736017_real64, &
            .2106637605030177_real64, &
            .1456505080593249_real64, &
            .6452301046087202_real64, &
            .9192324923978923_real64, &
            .8201528711770617_real64, &
            .1893330822201963_real64, &
            .3848617084851692_real64, &
            .2708002789961279_real64, &
            .8344516576865629_real64, &
            .4391079346296249_real64, &
            .1328552697285532_real64, &
            .2627126874078891_real64, &
            .7119541072827280_real64, &
            .2808800497440635E-01_real64, &
            .1830050086395000_real64, &
            .4650572894745821_real64, &
            .8372880332961116_real64, &
            .3196942043369664_real64, &
            .5400293254770759E-02_real64, &
            .9969762506629233_real64, &
            .2318563543762030E-01_real64, &
            .9735861019020967_real64, &
            .9475508553636347_real64, &
            .2784847480576497_real64, &
            .5557827857516401_real64, &
            .1384686002367318_real64, &
            .9060128503530698_real64, &
            .9320127662369427_real64, &
            .8995943728208159_real64, &
            .7605201819186155E-01_real64, &
            .1059943730987242_real64, &
            .3127855474999036_real64, &
            .8654963791217601_real64, &
            .6831660126657728_real64, &
            .5772796995326747_real64, &
            .2227835925000740_real64, &
            .1131337017828432_real64, &
            .1835807126967049_real64, &
            .2075011206965867_real64, &
            .4369390728603207_real64, &
            .6447316495711857_real64, &
            .9422622644772234_real64, &
            .8429110262763843_real64, &
            .9862119722310516_real64, &
            .9369999853172315_real64, &
            .5891046589584525_real64, &
            .6494055880231188_real64, &
            .9662886085615155_real64, &
            .9711014840489356_real64, &
            .5257859317081990_real64, &
            .6882911596934989_real64, &
            .5927138702853939_real64, &
            .5935388759015137_real64, &
            .2790030691378096_real64, &
            .7553594825815190_real64, &
            .5091780658253159_real64, &
            .4399963740257965_real64, &
            .9209253584400991E-01_real64, &
            .2040092341326103_real64, &
            .6148179509954169_real64, &
            .5084340270687314_real64, &
            .8918018540382170_real64, &
            .8530268340871098_real64, &
            .8124369420377242_real64, &
            .1156591227422585_real64, &
            .9123843589577203_real64, &
            .6584092951215744_real64, &
            .8168219441687045_real64, &
            .7420874098676381_real64, &
            .8228361261734121_real64, &
            .2948578228535660_real64, &
            .3764548185275206_real64, &
            .9267980329589121_real64, &
            .7806132786684322_real64, &
            .9407241982228834_real64, &
            .7712658739451300_real64, &
            .9638930987244980_real64, &
            .5316412538241150_real64, &
            .8397493657556793_real64, &
            .5174919785870313_real64, &
            .8328148768233419_real64, &
            .8417288831188118_real64, &
            .7833064271048843_real64, &
            .9183904365930627_real64, &
            .9040383956596534_real64, &
            .2786629813194053_real64, &
            .9063587859461516_real64, &
            .4205534144989393_real64, &
            .7338328060581281_real64, &
            .3022583251315897_real64, &
            .7130571108357289_real64, &
            .2954993695567572_real64, &
            .7061266927145262_real64, &
            .2331531852567660_real64, &
            .6652184420002460_real64, &
            .5525761190657353_real64 &
        ]

        expected_ix = [ 23, 25, 18, 34, 62, 35, 41, 69, 15, 30, 6, 19, 42, 10, 63, 43, 5, 40, 98, 16, 12, 28, 90, 58, 75, 96, 94, 36, 22, 2, 76, 11, 92, 44, 14, 61, 4, 20, 65, 60, 84, 54, 82, 100, 29, 39, 50, 56, 57, 64, 45, 7, 51, 71, 99, 38, 55, 97, 17, 95, 93, 73, 59, 1, 80, 78, 87, 68, 72, 9, 74, 85, 13, 21, 3, 83, 86, 47, 67, 37, 66, 33, 89, 31, 91, 70, 88, 8, 77, 32, 49, 79, 46, 27, 81, 52, 53, 26, 48, 24 ]
        
        call build_bst_index(x, n, ix, ierr)
        call assert_equal_int(ierr, ERR_OK, "bst construction failed")
        call assert_permutation(ix, n, "ix not permutation")

        ! Check monotonicity of x
        is_sorted = .true.
        do i = 2, n
            if (x(ix(i)) < x(ix(i - 1))) then
                is_sorted = .false.
                exit
            end if
        end do
        call assert_true(is_sorted, "BST index test FAILED: x(ix) is not monotonic.")

        call assert_true(is_sorted, "BST index test FAILED: x(ix) is not monotonic.")
    end subroutine test_bst_index_construction

    !> Test get_sorted_value function.
    subroutine test_bst_sorted_values()
        integer(int32), parameter :: n = 10
        real(real64) :: x(n) = [3.0d0, 1.0d0, 4.0d0, 2.0d0, 5.0d0, 7.0d0, 6.0d0, 9.0d0, 8.0d0, 10.0d0]
        integer(int32) :: ix(n)
        real(real64) :: val
        integer(int32) :: ierr

        call set_ok(ierr)

        call build_bst_index(x, n, ix, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Build bst index failed for sorted values')
        val = get_sorted_value(x, ix, 3, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Get sorted value failed')

        call assert_equal_real(val, 3.0d0, 1d-12, "get_sorted_value returned incorrect value")
    end subroutine test_bst_sorted_values

    !> Test BST range query functionality.
    subroutine test_bst_range_query()
        integer(int32), parameter :: n = 10
        real(real64) :: x(n) = [3.0d0, 1.0d0, 4.0d0, 2.0d0, 5.0d0, 7.0d0, 6.0d0, 9.0d0, 8.0d0, 10.0d0]
        integer(int32) :: ix(n)
        integer(int32) :: res_ix(n), res_n, ierr

        call set_ok(ierr)

        call build_bst_index(x, n, ix, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Build bst index failed')
        call bst_range_query(x, ix, n, 2.5d0, 7.5d0, res_ix, res_n, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Bst range query failed')

        call assert_true(res_n == 5, "BST range query returned incorrect count")
    end subroutine test_bst_range_query

    !> Test BST with empty array.
    subroutine test_bst_empty_array()
        integer(int32), parameter :: n = 0
        real(real64) :: x(n)
        integer(int32) :: ix(n), ierr

        call set_ok(ierr)

        call build_bst_index(x, n, ix, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, arg_pos=2_int32), 'Expected error for emtpy input')
        call assert_true(.true., "BST empty array handling")
    end subroutine test_bst_empty_array

    !> Test BST with single element.
    subroutine test_bst_single_element()
        integer(int32), parameter :: n = 1
        real(real64) :: x(n) = [42.0d0]
        integer(int32) :: ix(n), ierr

        call build_bst_index(x, n, ix, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Build bst index failed for single element')
        call assert_equal_int(ix(1), 1, "BST single element index incorrect")
        call assert_equal_real(x(ix(1)), 42.0d0, 1d-12, "BST single element value incorrect")
    end subroutine test_bst_single_element

    !> Test BST with identical values.
    subroutine test_bst_identical_values()
        integer(int32), parameter :: n = 5
        real(real64) :: x(n) = 7.0d0
        integer(int32) :: ix(n)
        integer(int32) :: i, ierr

        call set_ok(ierr)

        call build_bst_index(x, n, ix, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Build bst index failed for identical values')
        ! Should still be a valid permutation
        do i = 1, n
            call assert_true(ix(i) >= 1 .and. ix(i) <= n, "BST identical values index out of bounds")
        end do
    end subroutine test_bst_identical_values

    !> Test BST with large random array.
    subroutine test_bst_large_random()
        integer(int32), parameter :: n = 1000
        real(real64) :: x(n)
        integer(int32) :: ix(n)
        integer(int32) :: i, ierr
        logical :: is_sorted

        call set_ok(ierr)
        call random_array(x, n)
        call build_bst_index(x, n, ix, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Build bst index failed for large random values')
        ! Check monotonicity of x(ix)
        is_sorted = .true.
        do i = 2, n
            if (x(ix(i)) < x(ix(i - 1))) then
                is_sorted = .false.
                exit
            end if
        end do

        call assert_true(is_sorted, "Large BST index test FAILED: x(ix) is not monotonic.")
    end subroutine test_bst_large_random

    !> Fill an array with random real values in [0,1).
    subroutine random_array(arr, nval)
        real(real64), intent(out) :: arr(:)
        integer(int32), intent(in) :: nval
        integer(int32) :: j
        call random_seed()
        do j = 1, nval
            call random_number(arr(j))
        end do
    end subroutine random_array

end module mod_test_bst
