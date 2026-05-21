! filepath: test/test_suite.f90
!> Common definitions for all test suites.
module test_suite
    implicit none

    !> Common interface for all test procedures.
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    !> Type to hold information about a single test.
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    !> Abstract interface for getting all tests in a suite.
    abstract interface
        function get_all_interface() result(all_tests)
            import :: test_case
            type(test_case), allocatable :: all_tests(:)
        end function get_all_interface
    end interface

    !> Type to hold suite information.
    type :: suite_entry
        character(len=128) :: name
        procedure(get_all_interface), pointer, nopass :: get_all => null()
    end type suite_entry

    type(suite_entry), allocatable :: available_suites(:)

    character(len=11) :: COLOR_GREEN = achar(27) // "[38;5;154m"
    character(len=11) :: COLOR_COPPER = achar(27) // "[38;5;214m"
    character(len=11) :: COLOR_DARK_COPPER = achar(27) // "[38;5;208m"
    character(len=11) :: COLOR_RED = achar(27) // "[38;5;196m"
    character(len=11) :: COLOR_LIGHT_GRAY = achar(27) // "[38;5;252m"
    character(len=11) :: COLOR_YELLOW = achar(27) // "[38;5;226m"
    character(len=11) :: COLOR_CREAM = achar(27) // "[38;5;255m"
    character(len=11) :: COLOR_ERROR = achar(27) // "[38;5;222m"
    character(len=4) :: COLOR_RESET = achar(27) // "[0m"

contains

    !> Register all test suites here.
    subroutine initialize_suites()
        allocate (available_suites(0))
    end subroutine initialize_suites

    !> Add a new test suite to the registry.
    subroutine add_suite(name, get_all_proc)
        character(len=*), intent(in) :: name
        procedure(get_all_interface) :: get_all_proc
        type(suite_entry), allocatable :: temp_suites(:)
        integer :: n

        n = size(available_suites)
        allocate (temp_suites(n + 1))

        if (n > 0) then
            temp_suites(1:n) = available_suites(1:n)
        end if

        temp_suites(n + 1)%name = name
        temp_suites(n + 1)%get_all => get_all_proc

        call move_alloc(temp_suites, available_suites)
    end subroutine add_suite

    !> Run all tests in a given suite.
    subroutine run_all_tests(suite_name, all_tests)
        character(len=*), intent(in) :: suite_name
        type(test_case), intent(in) :: all_tests(:)
        integer :: i

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(A)", COLOR_COPPER // " " // trim(all_tests(i)%name) // COLOR_GREEN // " passed." // COLOR_RESET
        end do

        print "(A)", COLOR_CREAM // "All '" // COLOR_DARK_COPPER // trim(suite_name) // COLOR_CREAM // "' tests passed successfully." // COLOR_RESET
    end subroutine run_all_tests

    !> Run selected tests by name from a given suite.
    subroutine run_named_tests(test_names, all_tests)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), intent(in) :: all_tests(:)
        integer :: i, j
        logical :: found, some_not_found

        some_not_found = .false.
        do i = 1, size(test_names)
            found = .false.

            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(A)", COLOR_COPPER // " " // trim(test_names(i)) // COLOR_GREEN // " passed." // COLOR_RESET
                    found = .true.
                    exit
                end if
            end do

            if (.not. found) then
                print "(A)", COLOR_RED // "Unknown test" // COLOR_CREAM // ": " // COLOR_COPPER // trim(test_names(i)) // COLOR_RESET
                some_not_found = .true.
            end if
        end do

        if (some_not_found) then
            stop 1
        end if
    end subroutine run_named_tests

    !> Run every registered suite.
    subroutine run_all_suites()
        integer :: i
        type(test_case), allocatable :: all_tests(:)

        do i = 1, size(available_suites)
            print "(A)", COLOR_CREAM // "Running suite: '" // COLOR_DARK_COPPER // trim(available_suites(i)%name) // "'" // COLOR_RESET
            all_tests = available_suites(i)%get_all()
            call run_all_tests(trim(available_suites(i)%name), all_tests)
        end do
    end subroutine run_all_suites

    !> Run one whole suite by name.
    subroutine run_suite_all(requested_suite)
        character(len=*), intent(in) :: requested_suite
        integer :: i
        type(test_case), allocatable :: all_tests(:)

        do i = 1, size(available_suites)
            if (trim(available_suites(i)%name) == trim(requested_suite)) then
                print "(A)", COLOR_CREAM // "Running suite: '" // COLOR_DARK_COPPER // trim(available_suites(i)%name) // "'" // COLOR_RESET
                all_tests = available_suites(i)%get_all()
                call run_all_tests(trim(requested_suite), all_tests)
                return
            end if
        end do

        print "(A)", COLOR_RED // "Unknown test suite" // COLOR_CREAM // ": " // COLOR_DARK_COPPER // trim(requested_suite) // COLOR_RESET
        stop 1
    end subroutine run_suite_all

    !> Run selected tests from one suite.
    subroutine run_suite_named(requested_suite, test_list)
        character(len=*), intent(in) :: requested_suite, test_list
        integer :: i
        type(test_case), allocatable :: all_tests(:)
        character(len=128), allocatable :: test_names(:)

        do i = 1, size(available_suites)
            if (trim(available_suites(i)%name) == trim(requested_suite)) then
                all_tests = available_suites(i)%get_all()
                call split_test_list(test_list, test_names)
                call run_named_tests(test_names, all_tests)
                return
            end if
        end do

        print "(A)", COLOR_RED // "Unknown test suite" // COLOR_CREAM // ": " // COLOR_DARK_COPPER // trim(requested_suite) // COLOR_RESET
        stop 1
    end subroutine run_suite_named

    !> Split a comma-separated list of test names into an array.
    subroutine split_test_list(test_list, test_names)
        character(len=*), intent(in) :: test_list
        character(len=128), allocatable, intent(out) :: test_names(:)

        integer :: i, n_items, lenstr
        integer :: start_pos, comma_pos, item_idx

        lenstr = len_trim(test_list)

        if (lenstr == 0) then
            allocate (test_names(0))
            return
        end if

        n_items = 1
        do i = 1, lenstr
            if (test_list(i:i) == ',') n_items = n_items + 1
        end do

        allocate (test_names(n_items))
        test_names = ""

        start_pos = 1
        item_idx = 1

        do
            comma_pos = index(test_list(start_pos:), ',')

            if (comma_pos == 0) then
                test_names(item_idx) = trim(adjustl(test_list(start_pos:lenstr)))
                exit
            else
                test_names(item_idx) = trim(adjustl(test_list(start_pos:start_pos + comma_pos - 2)))
                start_pos = start_pos + comma_pos
                item_idx = item_idx + 1
            end if
        end do
    end subroutine split_test_list

end module test_suite
