# Test Suite Framework

## Overview

This framework provides a robust and scalable system for organizing and executing unit tests in Fortran. It allows running individual tests, complete test suites, or all project tests with simple and clear syntax.

## System Architecture

### Main Components

1. **`run_tests.f90`** - Main program that handles command line arguments
2. **Test Modules** - Each module (suite) contains tests for a specific functionality
3. **`asserts.f90`** - Assertion function library for validating results

### File Structure
```
test/
├── run_tests.f90                     # Main program
├── asserts.f90                       # Assertion library
├── mod_test_normalize_by_std_dev.f90 # Normalization tests
├── mod_test_sorting.f90              # Sorting tests
└── mod_test_quantile_normalization.f90 # Quantile normalization tests
```

## System Usage

### Command Syntax

```bash
# Run all tests from all suites
./test_runner.sh

# Run all tests from a specific suite
./test_runner.sh <suite_name>

# Run specific tests from a suite
./test_runner.sh <suite_name> <test1,test2,test3>
```

### Usage Examples

```bash
# Run all tests
./test_runner.sh

# Run all normalization tests
./test_runner.sh normalization

# Run specific normalization tests
./test_runner.sh normalization test_identity_matrix,test_zero_rows

# Run multiple normalization tests
./test_runner.sh normalization test_normalize_by_std_dev_basic,test_large_random_matrix,test_single_nonzero
```

## How to Add a New Test

### 1. Add Test to an Existing Module

To add a new test to an existing module (e.g., `mod_test_normalize_by_std_dev.f90`):

#### Step 1: Write the test subroutine
```fortran
!> @brief Test description
subroutine test_my_new_test()
  ! Declare variables
  real(8), dimension(2,2) :: mat, result, expected
  
  ! Set up test data
  mat = reshape([1.0d0, 2.0d0, 3.0d0, 4.0d0], [2,2])
  
  ! Execute function under test
  call my_function_to_test(2, 2, mat, result)
  
  ! Define expected result
  expected = reshape([0.5d0, 1.0d0, 1.5d0, 2.0d0], [2,2])
  
  ! Validate result
  call assert_equal_array_real(result, expected, 4, 1d-12, "my_new_test: result failure")
end subroutine test_my_new_test
```

#### Step 2: Register the test in the `get_all_tests()` function
```fortran
function get_all_tests() result(all_tests)
  type(test_case) :: all_tests(14)  ! Increment the number
  
  ! ...existing tests...
  all_tests(13) = test_case("test_symmetric_rows", test_symmetric_rows)
  all_tests(14) = test_case("test_my_new_test", test_my_new_test)  ! Add this line
end function get_all_tests
```

#### Step 3: Update the counter in functions that use the array
```fortran
subroutine run_all_tests_normalize_by_std_dev()
  type(test_case) :: all_tests(14)  ! Change from 13 to 14
  ! ...rest of function...
end subroutine

subroutine run_named_tests_normalize_by_std_dev(test_names)
  type(test_case) :: all_tests(14)  ! Change from 13 to 14
  ! ...rest of function...
end subroutine
```

### 2. Create a New Test Module

#### Step 1: Create the module file
Create `test/mod_test_my_new_functionality.f90`:

```fortran
! filepath: test/mod_test_my_new_functionality.f90
!> @brief Unit test suite for my_new_functionality.
module mod_test_my_new_functionality
  use asserts
  implicit none
  public

  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> @brief Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(2)  ! Number of tests
    
    all_tests(1) = test_case("test_my_first_test", test_my_first_test)
    all_tests(2) = test_case("test_my_second_test", test_my_second_test)
  end function get_all_tests

  !> @brief Run all tests in this module.
  subroutine run_all_tests_my_new_functionality()
    type(test_case) :: all_tests(2)
    integer :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All my_new_functionality tests passed successfully."
  end subroutine run_all_tests_my_new_functionality

  !> @brief Run specific tests by name.
  subroutine run_named_tests_my_new_functionality(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(2)
    integer :: i, j
    logical :: found
    
    all_tests = get_all_tests()
    
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_my_new_functionality

  !> @brief My first test
  subroutine test_my_first_test()
    ! Implement test here
    call assert_true(.true., "test_my_first_test: dummy test")
  end subroutine test_my_first_test

  !> @brief My second test
  subroutine test_my_second_test()
    ! Implement test here
    call assert_true(.true., "test_my_second_test: dummy test")
  end subroutine test_my_second_test

end module mod_test_my_new_functionality
```

#### Step 2: Register the module in `run_tests.f90`

**Add the use statement:**
```fortran
program main
  use mod_test_normalize_by_std_dev
  use mod_test_my_new_functionality  ! Add this line
  implicit none
```

**Register the suite in `initialize_suites()`:**
```fortran
subroutine initialize_suites()
  ! Start with empty registry
  allocate(available_suites(0))
  
  ! Add each suite
  call add_suite("normalization", run_all_tests_normalize_by_std_dev, run_named_tests_normalize_by_std_dev)
  call add_suite("my_new_functionality", run_all_tests_my_new_functionality, run_named_tests_my_new_functionality)  ! Add this line
end subroutine initialize_suites
```

## Available Assertion Functions

The `asserts` module provides the following functions for validating results:

### Basic Assertions
- `assert_true(condition, message)` - Verifies that a condition is true
- `assert_false(condition, message)` - Verifies that a condition is false

### Numerical Assertions
- `assert_equal_real(actual, expected, tolerance, message)` - Compares real numbers
- `assert_equal_int(actual, expected, message)` - Compares integers
- `assert_in_range_real(value, min_val, max_val, message)` - Verifies range

### Array Assertions
- `assert_equal_array_real(actual, expected, size, tolerance, message)` - Compares real arrays
- `assert_no_nan_real(array, size, message)` - Verifies absence of NaN

### Assertion Usage Example
```fortran
subroutine test_assertion_examples()
  real(8) :: result, expected
  integer :: int_result
  real(8), dimension(3) :: array_result, array_expected
  
  ! Basic equality test
  result = 2.0d0
  expected = 2.0d0
  call assert_equal_real(result, expected, 1d-12, "Real equality failed")
  
  ! Integer test
  int_result = 5
  call assert_equal_int(int_result, 5, "Integer equality failed")
  
  ! Range test
  call assert_in_range_real(result, 1.5d0, 2.5d0, "Value out of range")
  
  ! Array test
  array_result = [1.0d0, 2.0d0, 3.0d0]
  array_expected = [1.0d0, 2.0d0, 3.0d0]
  call assert_equal_array_real(array_result, array_expected, 3, 1d-12, "Arrays differ")
  
  ! Condition test
  call assert_true(result > 1.0d0, "Condition false")
end subroutine test_assertion_examples
```

## Requirements

### Naming Convention
- **Modules:** `mod_test_<functionality>.f90` - **REQUIRED:** Must start with `mod_test_` because compilation occurs in alphabetical order and test modules must compile before `run_tests.f90`
- **Test subroutines:** `test_<specific_description>()`
- **Suite names:** use underscores to separate words

### Module Structure
- **Each module represents one main Fortran function. You can use the _r wrapper to call your subroutine.**
- **Each test within that module must be a separate subroutine**
- **All suites must implement two required subroutines:**
  - `run_all_tests_<suite_name>()` - Runs all tests in the suite
  - `run_named_tests_<suite_name>(test_names)` - Runs specific tests by name

### Test Structure
1. **Arrange:** Set up test data
2. **Act:** Execute the function under test
3. **Assert:** Validate the results

## System Advantages

1. **Scalable:** Adding new tests requires minimal code
2. **No merge conflicts:** Each developer can add suites independently
3. **Flexible:** Allows running individual tests or in groups
