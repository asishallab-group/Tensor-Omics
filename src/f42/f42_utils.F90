#include <src/macros.h>

!> Utility module for data analysis.
!| This module provides general-purpose utility functions for data analysis, to be used as needed.
module f42_utils
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: validate_all_in_range_int, ERR_DIVISION_BY_ZERO, set_ok, set_err, validate_in_range_real, is_err, validate_in_range_int, validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int, ERR_ALLOC_FAIL
    use, intrinsic :: ieee_arithmetic, only: ieee_next_after, ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_is_finite, ieee_is_nan
    implicit none
    public :: init_random, rand_range
    public :: sort_real, sort_integer, sort_character
    public :: sort_array
    public :: compute_edf, compute_edf_alloc

    interface sort_array
        module procedure sort_real, sort_integer, sort_character
    end interface sort_array

    interface sort_array_heapsort
        module procedure sort_real_heapsort, sort_integer_heapsort, sort_character_heapsort
        module procedure sort_real_heapsort_expl_size
    end interface sort_array_heapsort

    interface shuffle_vector
        module procedure shuffle_vector_real, shuffle_vector_int
    end interface shuffle_vector

    interface clamp
        module procedure clamp_real, clamp_int
    end interface clamp

    interface operator(.lessthan.)
        module procedure real_less
    end interface operator(.lessthan.)

    interface operator(.greaterthan.)
        module procedure real_greater
    end interface operator(.greaterthan.)

    interface operator(.isclose.)
        module procedure is_close
    end interface operator(.isclose.)

#define CM_EPS epsilon(1.0_real64)

    real(real64), parameter :: PI = 4.0_real64*atan(1.0_real64)
    real(real64), parameter :: EPS = CM_EPS
    real(real64), parameter :: LOG_2 = log(2.0_real64)
contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Initializes a permutation vector with ascending indices
    pure subroutine init_perm(perm)
        integer(int32), dimension(:), intent(out) :: perm
            !! Permutation vector to initialize with indices

        integer(int32) :: i_perm

        do concurrent(i_perm=1:size(perm, kind=int32)) shared(perm)
            perm(i_perm) = i_perm
        end do
    end subroutine init_perm

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Calculates the arithmetic mean of vector
    pure real(real64) function mean(vec)
        real(real64), dimension(:), intent(in) :: vec
            !! Vector to compute the mean value from

        mean = sum(vec)/real(size(vec, kind=int32), real64)
    end function mean

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Calculates the standard deviation of vector, with or without Bessel's correction
    pure real(real64) function std_dev(vec, do_bessel_correction)
        real(real64), dimension(:), intent(in) :: vec
            !! Vector to compute the standard deviation value from
        logical, intent(in), optional :: do_bessel_correction
            !! Tells whether to apply the bessel's correction or not, default: `.false.`
            !!
            !! |    Case     |                                                Formula                                                      |
            !! |-------------|-------------------------------------------------------------------------------------------------------------|
            !! |  `.true.`   | \(\frac{1}{\texttt{size}(vec) - 1} \cdot \sum_{i=1}^{\texttt{size}(vec)} (vec(i) - \texttt{mean}(i))^{2}\)  |
            !! |  `.false.`  |  \(\frac{1}{\texttt{size}(vec)} \cdot \sum_{i=1}^{\texttt{size}(vec)} vec(i)^{2} - \texttt{mean}(i)^{2}\)   |

        logical :: bessel
        integer(int32) :: n_elements, i_element
        real(real64) :: mean_val, squares_sum

        M_DEFAULT_VAL(do_bessel_correction, bessel, .false.)

        mean_val = mean(vec)
        n_elements = size(vec, kind=int32)
        if (bessel) then
            squares_sum = 0.0_real64
            do concurrent(i_element=1:n_elements) shared(vec, mean_val) reduce(+:squares_sum)
                squares_sum = squares_sum + (vec(i_element) - mean_val)**2
            end do
            std_dev = sqrt(squares_sum/real(n_elements - 1, kind=real64))
        else
            squares_sum = 0.0_real64
            do concurrent(i_element=1:n_elements) shared(vec) reduce(+:squares_sum)
                squares_sum = squares_sum + vec(i_element)**2
            end do
            std_dev = sqrt(max(0.0_real64, squares_sum/real(n_elements, kind=real64) - mean_val**2))
        end if
    end function std_dev

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Function to find the position to place a value in a sorted array using binary search
    pure integer(int32) function binary_search_insertion(arr, perm, value, lower_idx, upper_idx) result(idx)
        real(real64), dimension(:), contiguous, intent(in) :: arr
            !! Array of values
        integer(int32), dimension(size(arr, kind=int32)), intent(in) :: perm
            !! Sorting permutation of `arr`, NaNs sorted to the end
        real(real64), intent(in) :: value
            !! Value to find
        integer(int32), intent(in), optional :: lower_idx
            !! The lower index to start searching in the array, default `1` -> searching in `perm(lower_idx:upper_idx)`
        integer(int32), intent(in), optional :: upper_idx
            !! The upper index to stop searching in the array, default `size(arr)` -> searching in `perm(lower_idx:upper_idx)`

        integer(int32) :: actual_lower_idx, actual_upper_idx, mid_idx, n

        n = size(arr, kind=int32)

        ! NaN is sorted to the end -> If value and last element is NaN, return n, else it is not found
        if (ieee_is_nan(value)) then
            if (ieee_is_nan(arr(perm(n)))) then
                idx = n
            else
                idx = n + 1
            end if
            return
        end if

        M_DEFAULT_VAL(lower_idx, actual_lower_idx, 1_int32)
        actual_lower_idx = clamp(actual_lower_idx, min_val=1_int32, max_val=n)
        M_DEFAULT_VAL(upper_idx, actual_upper_idx, n + 1)
        actual_upper_idx = clamp(actual_upper_idx, min_val=actual_lower_idx - 1, max_val=n) + 1

        do while (actual_lower_idx < actual_upper_idx)
            mid_idx = (actual_lower_idx + actual_upper_idx)/2

            ! Update bounds: Note that NaN is sorted to the end -> if condition fails because of NaN, the upper bound will be reduced
            if (arr(perm(mid_idx)) < value) then
                actual_lower_idx = mid_idx + 1
            else
                actual_upper_idx = mid_idx
            end if
        end do

        idx = actual_lower_idx
    end function binary_search_insertion

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Function to find a value in a sorted array using binary search. Returns -1 if not found
    pure integer(int32) function binary_search(arr, perm, value) result(idx)
        real(real64), dimension(:), contiguous, intent(in) :: arr
            !! Array of values
        integer(int32), dimension(size(arr, kind=int32)), intent(in) :: perm
            !! Sorting permutation of `arr`, NaNs sorted to the end
        real(real64), intent(in) :: value
            !! Value to find

        integer(int32) :: n
        real(real64) :: found

        n = size(arr, kind=int32)

        ! Nothing found in empty arrays
        if (n <= 0) then
            idx = -1
        else
            idx = binary_search_insertion(arr, perm, value)

            ! nothing found if inserted as new element after last one
            if (idx == n + 1) then
                idx = -1
            else
                found = arr(perm(idx))

                ! allow NaN search -> if value doesn't match found one and both aren't NaN, it is not found
                if (found /= value .and. .not. (ieee_is_nan(found) .and. ieee_is_nan(value))) then
                    idx = -1
                end if
            end if
        end if
    end function binary_search

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Clamps a value into a range `min_val <= val <= max_val`. If `max_val < min_val`, `min_val` is returned
    pure real(real64) function clamp_real(val, min_val, max_val) result(clamped)
        real(real64), intent(in) :: val
            !! Value to be clamped
        real(real64), intent(in) :: min_val
            !! Lower bound
        real(real64), intent(in) :: max_val
            !! Upper bound

        clamped = max(min_val, min(val, max_val))
    end function clamp_real

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Clamps a value into a range `min_val <= val <= max_val`. If `max_val < min_val`, `min_val` is returned
    pure integer(int32) function clamp_int(val, min_val, max_val) result(clamped)
        integer(int32), intent(in) :: val
            !! Value to be clamped
        integer(int32), intent(in) :: min_val
            !! Lower bound
        integer(int32), intent(in) :: max_val
            !! Upper bound

        clamped = max(min_val, min(val, max_val))
    end function clamp_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Compute logarithm for any base
    pure subroutine logx(val, base, exponent, ierr)
        real(real64), intent(in) :: val
            !! Value (`x` in \( b^y = x \))
        real(real64), intent(in) :: base
            !! Base (`b` in \( b^y = x \))
        real(real64), intent(out) :: exponent
            !! Exponent (`y` in \( b^y = x \))
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_in_range_real(val, ierr, min=above(0.0_real64), arg_pos=1_int32)
        call validate_in_range_real(base, ierr, min=above(0.0_real64), arg_pos=2_int32)
        if (is_close(base, 1.0_real64)) call set_err(ierr, ERR_DIVISION_BY_ZERO, arg_pos=2_int32)

        if (is_err(ierr)) return

        if (base == 2.0_real64) then
            exponent = log(val)/LOG_2
        else
            exponent = log(val)/log(base)
        end if
    end subroutine logx

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Initialize Fortran's random number generator
    subroutine init_random(seed)
        integer(int32), intent(in), optional :: seed
            !! optional random seed, default: 42
            !!
            !! @note
            !! This subroutine uses the intrinsic `random_seed` subroutine that expects default kind integer.
            !! To map `int32` to default kind, it is being clamped by modulo to be in range.
            !! @endnote

        integer(int32) :: actual_seed

        ! IMPORTANT: these locals need to be default kind integer
        integer :: seed_default_kind, size, i
        integer(int32), parameter :: max_default_kind_val = int(huge(1), kind=int32)

        M_DEFAULT_VAL(seed, actual_seed, 42_int32)

        ! clamp to default kind range
        seed_default_kind = int(mod(actual_seed, max_default_kind_val))

        ! determine needed array size to seed random numbers
        call random_seed(size=size)

        ! create the seeding array, has negligible size
        call random_seed(put=[(seed_default_kind, i=1, size)])
    end subroutine init_random

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns a random real number `min <= rand_num < max`. If `min > max`, it will be `max <= rand_num < min`. If `min == max`, it will be `min`.
    real(real64) function rand_range(min, max) result(rand_num)
        real(real64), intent(in) :: min
            !! Lower bound
        real(real64), intent(in) :: max
            !! Upper bound

        call random_number(rand_num)
        rand_num = min + rand_num*(max - min)
    end function rand_range

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Shuffle a vector in-place, using Fisher-Yates shuffle
    subroutine shuffle_vector_real(vec)
        real(real64), dimension(:), intent(inout) :: vec
            !! Output permutation array

        integer(int32) :: i, rand_idx

        ! Fisher-Yates shuffle
        do i = size(vec, kind=int32), 2, -1
            ! Generate random integer in range [1, i]
            rand_idx = int(rand_range(1.0_real64, real(i, real64)), int32)

            call swap_real(vec(i), vec(rand_idx))
        end do
    end subroutine shuffle_vector_real

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Shuffle a vector in-place, using Fisher-Yates shuffle
    subroutine shuffle_vector_int(vec)
        integer(int32), dimension(:), intent(inout) :: vec
            !! Output permutation array

        integer(int32) :: i, rand_idx

        ! Fisher-Yates shuffle
        do i = size(vec, kind=int32), 2, -1
            ! Generate random integer in range [1, i]
            rand_idx = int(rand_range(1.0_real64, real(i, real64)), int32)

            call swap_int(vec(i), vec(rand_idx))
        end do
    end subroutine shuffle_vector_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns the next representable float lower than a value. Helpful for exclusive upper bounds in ranges. Doesn't return denormals, thus `below(0.0_real64)==-tiny(1.0_real64)` and `below(tiny(1.0_real64))==0.0_real64`
    pure real(real64) function below(val)
        real(real64), intent(in) :: val

        if (val == 0.0_real64) then
            below = -tiny(1.0_real64)
        else if (val == tiny(1.0_real64)) then
            below = 0.0_real64
        else
            below = ieee_next_after(val, M_NEG_INF)
        end if
    end function below

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns the next representable float greater than a value. Helpful for exclusive upper bounds in ranges. Doesn't return denormals, thus `above(0.0_real64)==tiny(1.0_real64)` and `above(-tiny(1.0_real64))==0.0_real64`
    pure real(real64) function above(val)
        real(real64), intent(in) :: val

        if (val == 0.0_real64) then
            above = tiny(1.0_real64)
        else if (val == -tiny(1.0_real64)) then
            above = 0.0_real64
        else
            above = ieee_next_after(val, M_POS_INF)
        end if
    end function above

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Checks if two values are close to eachother, using a tolerance of `max(1d-12, EPS*max(abs(a), abs(b)))` with `EPS=CM_EPS`
    pure logical function is_close(a, b)
        real(real64), intent(in) :: a
            !! First variable of comparison `a==b`
        real(real64), intent(in) :: b
            !! Second variable of comparison `a==b`

        real(real64) :: rel_tolerance

        if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
            rel_tolerance = EPS*max(abs(a), abs(b))
            is_close = abs(a - b) <= max(rel_tolerance, 1d-12)
        else
            is_close = a == b
        end if
    end function is_close

    !> AUTHOR_AARON_SCHROEDER
    !| Find the next power of two greater than or equal to n
    function next_power_of_two(n) result(power)
        integer(int32), intent(in) :: n
            !! input value
        integer(int32) :: power
            !! next greater value that is a power of two

        power = 2**(bit_size(n) - leadz(n - 1))
    end function next_power_of_two

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Computes the radian angle between two vectors
    pure subroutine angle_between(v1, v2, n_dims, angle, ierr)
        integer(int32), intent(in) :: n_dims
            !! number of elements in `v1` and `v2`
        real(real64), dimension(n_dims), intent(in) :: v1
            !! first vector for angle calculation
        real(real64), dimension(n_dims), intent(in) :: v2
            !! second vector for angle calculation
        real(real64), intent(out) :: angle
            !! will hold calculated angle
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_dim
        real(real64) :: theta, dot_product, norm1, norm2, norm_product

        call set_ok(ierr)

        dot_product = 0.0_real64
        norm1 = 0.0_real64
        norm2 = 0.0_real64
        do concurrent(i_dim=1:n_dims) shared(v1, v2) reduce(+:dot_product, norm1, norm2)
            dot_product = dot_product + v1(i_dim)*v2(i_dim)
            norm1 = norm1 + v1(i_dim)**2
            norm2 = norm2 + v2(i_dim)**2
        end do

        norm_product = sqrt(norm1)*sqrt(norm2)
        if (is_close(norm_product, 0.0_real64)) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        theta = dot_product/norm_product
        theta = clamp(theta, -1.0_real64, 1.0_real64)
        angle = acos(theta)
    end subroutine angle_between

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns the given degrees in positive radian value \( -90^{\circ} \Rightarrow \frac{3\cdot \pi}{2}, \text{not} -\frac{\pi}{2} \)
    pure real(real64) function radians(degrees)
        real(real64), intent(in) :: degrees
            !! degrees to be converted

        radians = modulo(degrees, 360.0_real64)*PI/180
    end function radians

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns the given radians in positive degree value \( -\frac{\pi}{2} \Rightarrow 270^{\circ}, \text{not} -90^{\circ} \)
    pure real(real64) function degrees(radians)
        real(real64), intent(in) :: radians
            !! radians to be converted

        degrees = modulo(radians, 2*PI)*180/PI
    end function degrees

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Calculates the euclidean norm of a vector
    pure real(real64) function norm(vector)
        real(real64), dimension(:), intent(in) :: vector
            !! Input vector the norm will be calcuated for

        integer(int32) :: i_dim
        real(real64) :: norm_val

        norm_val = 0.0_real64
        do concurrent(i_dim=1:size(vector)) shared(vector) reduce(+:norm_val)
            norm_val = norm_val + vector(i_dim)**2
        end do
        norm = sqrt(norm_val)
    end function norm

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Adds two vectors in-place
    pure subroutine add_vector(vector, to_be_added)
        real(real64), dimension(:), intent(inout) :: vector
            !! First vector, it will be modified in-place
        real(real64), dimension(:), intent(in) :: to_be_added
            !! Vector that should be added to `vector`

        integer(int32) :: i_dim

        do concurrent(i_dim=1:size(vector)) shared(vector, to_be_added)
            vector(i_dim) = vector(i_dim) + to_be_added(i_dim)
        end do
    end subroutine add_vector

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Subtracts two vectors in-place
    pure subroutine subtract_vector(vector, to_be_subtracted)
        real(real64), dimension(:), intent(inout) :: vector
            !! First vector, it will be modified in-place
        real(real64), dimension(:), intent(in) :: to_be_subtracted
            !! Vector that should be subtracted from `vector`

        integer(int32) :: i_dim

        do concurrent(i_dim=1:size(vector)) shared(vector, to_be_subtracted)
            vector(i_dim) = vector(i_dim) - to_be_subtracted(i_dim)
        end do
    end subroutine subtract_vector

    !> AUTHOR_VIVIAN_BASS
    !| Sort a real array indirectly using quicksort.
    !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
    pure subroutine sort_real(array, perm, tmp_stack_left, tmp_stack_right)
        real(real64), intent(in) :: array(:)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(:)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(:)
            !! Manual stack of right indices for quicksort recursion
        call quicksort_real(array, perm, int(size(array), int32), tmp_stack_left, tmp_stack_right)
    end subroutine sort_real

    !> AUTHOR_VIVIAN_BASS
    !| Sort an integer array indirectly using quicksort.
    !| Similar to `sort_real`, but for integer input.
    pure subroutine sort_integer(array, perm, tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: array(:)
            !! Integer input array to sort
        integer(int32), intent(inout) :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(:)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(:)
            !! Manual stack of right indices for quicksort recursion
        call quicksort_int(array, perm, int(size(array), int32), tmp_stack_left, tmp_stack_right)
    end subroutine sort_integer

    !> AUTHOR_VIVIAN_BASS
    !| Sort a character array indirectly using quicksort.
    !| Uses lexicographic ordering and permutation vector sorting.
    pure subroutine sort_character(array, perm, tmp_stack_left, tmp_stack_right)
        character(len=*), intent(in) :: array(:)
            !! Character input array to sort
        integer(int32), intent(inout) :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(:)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(:)
            !! Manual stack of right indices for quicksort recursion
        call quicksort_char(array, perm, int(size(array), int32), tmp_stack_left, tmp_stack_right)
    end subroutine sort_character

    !> AUTHOR_MOHAMED_AKDI
    !| Sort a real array indirectly using heapsort.
    !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
    pure subroutine sort_real_heapsort(array, perm)
        real(real64), intent(in), contiguous :: array(:)
            !! Real input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        call heapsort_real(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_real_heapsort

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Sort a real explicit-size array indirectly using heapsort.
    !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
    pure subroutine sort_real_heapsort_expl_size(array, perm, n)
        integer(int32), intent(in) :: n
            !! Size of `array`
        real(real64), intent(in) :: array(n)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        call heapsort_real(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_real_heapsort_expl_size

    !> AUTHOR_MOHAMED_AKDI
    !| Sort an integer array indirectly using heapsort.
    !| Similar to `sort_real_heapsort`, but for integer input.
    pure subroutine sort_integer_heapsort(array, perm)
        integer(int32), intent(in), contiguous :: array(:)
            !! Integer input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        call heapsort_integer(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_integer_heapsort

    !> AUTHOR_MOHAMED_AKDI
    !| Sort a character array indirectly using heapsort.
    !| Uses lexicographic ordering and permutation vector sorting.
    pure subroutine sort_character_heapsort(array, perm)
        character(len=*), intent(in), contiguous :: array(:)
            !! Character input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        call heapsort_character(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_character_heapsort

    !> AUTHOR_VIVIAN_BASS
    !| Internal quicksort implementation for real arrays.
    !| Sorts indirectly using the permutation vector `perm`. Manual stack replaces recursion.
    pure subroutine quicksort_real(array, perm, n, tmp_stack_left, tmp_stack_right)
        real(real64), intent(in) :: array(n)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: n
            !! Size of the array
        integer(int32), intent(out) :: tmp_stack_left(n)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(n)
            !! Manual stack of right indices for quicksort recursion

        integer(int32) :: left, right, i, j, top, pivot_idx
        real(real64) :: pivot_val

        if (n == 0) return

        top = 1
        tmp_stack_left(top) = 1
        tmp_stack_right(top) = n

        ! Iterative quicksort using explicit stack
        do while (top > 0)
            left = tmp_stack_left(top)
            right = tmp_stack_right(top)
            top = top - 1

            if (left >= right) cycle

            ! Select pivot and initialize pointers
            pivot_idx = (left + right)/2
            pivot_val = array(perm(pivot_idx))
            i = left
            j = right

            ! Partitioning loop
            do
                do while (array(perm(i)) .lessthan.pivot_val)
                    i = i + 1
                end do
                do while (array(perm(j)) .greaterthan.pivot_val)
                    j = j - 1
                end do
                if (i <= j) then
                    call swap_int(perm(i), perm(j))
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do

            ! Push new ranges onto stack
            if (left < j) then
                top = top + 1
                tmp_stack_left(top) = left
                tmp_stack_right(top) = j
            end if
            if (i < right) then
                top = top + 1
                tmp_stack_left(top) = i
                tmp_stack_right(top) = right
            end if
        end do
    end subroutine quicksort_real

    !> AUTHOR_VIVIAN_BASS
    !| Internal quicksort implementation for integer arrays.
    !| Indirectly sorts `array` using `perm`, same algorithm as `quicksort_real`.
    pure subroutine quicksort_int(array, perm, n, tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: array(n)
            !! Integer input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: n
            !! Size of the array
        integer(int32), intent(out) :: tmp_stack_left(n)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(n)
            !! Manual stack of right indices for quicksort recursion

        integer(int32) :: left, right, i, j, top, pivot_idx
        integer(int32) :: pivot_val

        if (n == 0) return

        top = 1
        tmp_stack_left(top) = 1
        tmp_stack_right(top) = n

        do while (top > 0)
            left = tmp_stack_left(top)
            right = tmp_stack_right(top)
            top = top - 1

            if (left >= right) cycle

            pivot_idx = (left + right)/2
            pivot_val = array(perm(pivot_idx))
            i = left
            j = right

            do
                do while (array(perm(i)) < pivot_val)
                    i = i + 1
                end do
                do while (array(perm(j)) > pivot_val)
                    j = j - 1
                end do
                if (i <= j) then
                    call swap_int(perm(i), perm(j))
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do

            if (left < j) then
                top = top + 1
                tmp_stack_left(top) = left
                tmp_stack_right(top) = j
            end if
            if (i < right) then
                top = top + 1
                tmp_stack_left(top) = i
                tmp_stack_right(top) = right
            end if
        end do
    end subroutine quicksort_int

    !> AUTHOR_VIVIAN_BASS
    !| Internal quicksort implementation for character arrays.
    !| Lexicographic quicksort using string comparison, indirect via `perm`.
    pure subroutine quicksort_char(array, perm, n, tmp_stack_left, tmp_stack_right)
        character(len=*), intent(in) :: array(n)
            !! Character input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: n
            !! Size of the array
        integer(int32), intent(out) :: tmp_stack_left(n)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(n)
            !! Manual stack of right indices for quicksort recursion
        integer(int32) :: left, right, i, j, top, pivot_idx
            !! Temporary variables
        character(len=len(array)) :: pivot_val

        if (n == 0) return

        top = 1
        tmp_stack_left(top) = 1
        tmp_stack_right(top) = n

        do while (top > 0)
            left = tmp_stack_left(top)
            right = tmp_stack_right(top)
            top = top - 1

            if (left >= right) cycle

            pivot_idx = (left + right)/2
            pivot_val = array(perm(pivot_idx))
            i = left
            j = right

            do
                do while (array(perm(i)) < pivot_val)
                    i = i + 1
                end do
                do while (array(perm(j)) > pivot_val)
                    j = j - 1
                end do
                if (i <= j) then
                    call swap_int(perm(i), perm(j))
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do

            if (left < j) then
                top = top + 1
                tmp_stack_left(top) = left
                tmp_stack_right(top) = j
            end if
            if (i < right) then
                top = top + 1
                tmp_stack_left(top) = i
                tmp_stack_right(top) = right
            end if
        end do
    end subroutine quicksort_char

    !Internal heapsort implementations for real arrays.
    !> AUTHOR_MOHAMED_AKDI
    !| Sorts indirectly using the permutation vector `perm`. Uses `heapify_real` to maintain heap property.
    pure subroutine heapsort_real(array, perm, n_arr, n_perm)
        integer(int32), intent(in) :: n_arr
            !! Size of `arr`
        integer(int32), intent(in) :: n_perm
            !! Size of `perm`
        real(real64), intent(in) :: array(n_arr)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(n_perm)
            !! Permutation vector that will be sorted
        integer(int32) :: n, i
        n = n_perm

        ! Build max heap
        do i = n/2, 1, -1
            call heapify_real(array, perm, n, i)
        end do
        ! Heap sort
        do i = n, 2, -1
            call swap_int(perm(1), perm(i))
            call heapify_real(array, perm, i - 1, 1)
        end do
    end subroutine heapsort_real

    !> AUTHOR_MOHAMED_AKDI
    !| Iterative heapify (non-recursive, pure)
    !| Restore the max-heap property for the subtree rooted at `root`.
    !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
    !| We reorder the permutation vector `perm` (indices) rather than moving
    !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_real(array, perm, heap_size, root)
        real(real64), intent(in), contiguous :: array(:)
            !! Real input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: heap_size
            !! Size of the heap
        integer(int32), intent(in) :: root
            !! Root index

        integer(int32) :: current, next_idx, largest_idx

        current = root

        do
            ! Use `largest_idx` directly as the left-child index: left = 2*current
            largest_idx = 2*current
            ! If there is no left child the subtree is a leaf; we're done
            if (largest_idx > heap_size) exit

            next_idx = largest_idx + 1

            ! Only compare the right child (next_idx) when it actually exists
            if (next_idx <= heap_size) then
                if (array(perm(next_idx)) .greaterthan.array(perm(largest_idx))) then
                    largest_idx = next_idx
                end if
            end if

            ! If the larger child is greater than current, swap permutation entries
            if (array(perm(largest_idx)) .greaterthan.array(perm(current))) then
                call swap_int(perm(current), perm(largest_idx))
                current = largest_idx
            else
                exit
            end if
        end do
    end subroutine heapify_real

    !> AUTHOR_MOHAMED_AKDI
    !| Internal heapsort implementation for integer arrays.
    !| Indirectly sorts `array` using `perm`, same algorithm as `heapsort_real`.
    pure subroutine heapsort_integer(array, perm, n_arr, n_perm)
        integer(int32), intent(in) :: n_arr
            !! Size of `arr`
        integer(int32), intent(in) :: n_perm
            !! Size of `perm`
        integer(int32), intent(in) :: array(n_arr)
            !! Integer input array to sort
        integer(int32), intent(inout) :: perm(n_perm)
            !! Permutation vector that will be sorted

        integer(int32) :: n, i

        n = n_perm

        ! Build max-heap
        do i = n/2, 1, -1
            call heapify_integer(array, perm, n, i)
        end do

        ! Heap sort
        do i = n, 2, -1
            call swap_int(perm(1), perm(i))
            call heapify_integer(array, perm, i - 1, 1)
        end do
    end subroutine heapsort_integer

    !> AUTHOR_MOHAMED_AKDI
    !| Iterative heapify (non-recursive, pure)
    !| Restore the max-heap property for the subtree rooted at `root`.
    !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
    !| We reorder the permutation vector `perm` (indices) rather than moving
    !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_integer(array, perm, heap_size, root)
        integer(int32), intent(in), contiguous :: array(:)
            !! Integer input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: heap_size
            !! Size of the heap
        integer(int32), intent(in) :: root
            !! Root index

        integer(int32) :: current, next_idx, largest_idx

        current = root

        do
            ! Compute left-child index (use largest_idx as left to avoid extra var)
            largest_idx = 2*current
            ! If there is no left child the subtree is a leaf; nothing to do
            if (largest_idx > heap_size) exit

            ! Compute right-child index (may be out of heap bounds)
            next_idx = largest_idx + 1

            ! Only compare right-child when it actually exists (guarded access)
            if (next_idx <= heap_size) then
                if (array(perm(next_idx)) > array(perm(largest_idx))) then
                    ! Right child is larger than left child
                    largest_idx = next_idx
                end if
            end if

            ! Compare the selected child with the current node; if the child is
            ! greater, swap permutation indices so the larger value moves up the
            ! heap. We swap entries of `perm` (indices), not the array values.
            if (array(perm(largest_idx)) > array(perm(current))) then
                call swap_int(perm(current), perm(largest_idx))
                current = largest_idx
            else
                exit
            end if
        end do
    end subroutine heapify_integer

    !> AUTHOR_MOHAMED_AKDI
    !| Internal heapsort implementation for character arrays.
    !| Lexicographic heapsort using string comparison, indirect via `perm`.
    pure subroutine heapsort_character(array, perm, n_arr, n_perm)
        integer(int32), intent(in) :: n_arr
            !! Size of `arr`
        integer(int32), intent(in) :: n_perm
            !! Size of `perm`
        character(len=*), intent(in)    :: array(n_arr)
            !! Character input array to sort
        integer(int32), intent(inout) :: perm(n_perm)
            !! Permutation vector that will be sorted

        integer(int32) :: n, i

        n = n_perm

        ! Build max-heap
        do i = n/2, 1, -1
            call heapify_character(array, perm, n, i)
        end do

        ! Heap sort
        do i = n, 2, -1
            call swap_int(perm(1), perm(i))
            call heapify_character(array, perm, i - 1, 1)
        end do
    end subroutine heapsort_character

    !> AUTHOR_MOHAMED_AKDI
    !| Iterative heapify (non-recursive, pure)
    !| Restore the max-heap property for the subtree rooted at `root`.
    !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
    !| We reorder the permutation vector `perm` (indices) rather than moving
    !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_character(array, perm, heap_size, root)
        character(len=*), intent(in), contiguous :: array(:)
            !! Character input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: heap_size
            !! Size of the heap
        integer(int32), intent(in) :: root
            !! Root index

        integer(int32) :: current, next_idx, largest_idx

        current = root

        do
            ! Compute left-child index (use largest_idx directly as left = 2*current)
            largest_idx = 2*current
            ! If there is no left child the subtree is a leaf; nothing to do
            if (largest_idx > heap_size) exit

            ! Potential right child index
            next_idx = largest_idx + 1

            ! Only compare right child when it exists to avoid OOB access
            if (next_idx <= heap_size) then
                if (array(perm(next_idx)) > array(perm(largest_idx))) then
                    largest_idx = next_idx
                end if
            end if

            ! If the selected child is larger than current, swap permutation
            ! indices so the larger element moves up. We swap entries in `perm`.
            if (array(perm(largest_idx)) > array(perm(current))) then
                call swap_int(perm(current), perm(largest_idx))
                current = largest_idx
            else
                exit
            end if
        end do
    end subroutine heapify_character

    !> AUTHOR_VIVIAN_BASS
    !| Swap two integer values in-place.
    pure subroutine swap_int(a, b)
        integer(int32), intent(inout) :: a
            !! First integer to swap
        integer(int32), intent(inout) :: b
            !! Second integer to swap

        integer(int32) :: temp
        temp = a; a = b; b = temp
    end subroutine swap_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Swap two real values in-place.
    pure subroutine swap_real(a, b)
        real(real64), intent(inout) :: a
            !! First real to swap
        real(real64), intent(inout) :: b
            !! Second real to swap
        real(real64) :: temp
        temp = a; a = b; b = temp
    end subroutine swap_real

    !> AUTHOR_MOHAMED_AKDI
    !| Helper: NaN-aware comparisons for real(real64). Treats NaN as greater as every real number, while `NaN` equals `NaN`.
    pure logical function real_less(a, b)
        real(real64), intent(in) :: a
            !! left side of the comparison `a < b`
        real(real64), intent(in) :: b
            !! right side of the comparison `a < b`

        if (ieee_is_nan(a) .and. ieee_is_nan(b)) then
            real_less = .false.
        else if (ieee_is_nan(a)) then
            real_less = .false.
        else if (ieee_is_nan(b)) then
            real_less = .true.
        else
            real_less = (a < b)
        end if
    end function real_less

    !> AUTHOR_MOHAMED_AKDI
    !| Helper: NaN-aware comparisons for real(real64). Treats NaN as greater as every real number, while `NaN` equals `NaN`.
    pure logical function real_greater(a, b)
        real(real64), intent(in) :: a
            !! left side of the comparison `a > b`
        real(real64), intent(in) :: b
            !! right side of the comparison `a > b`

        if (ieee_is_nan(a) .and. ieee_is_nan(b)) then
            real_greater = .false.
        else if (ieee_is_nan(a)) then
            real_greater = .true.
        else if (ieee_is_nan(b)) then
            real_greater = .false.
        else
            real_greater = (a > b)
        end if
    end function real_greater

    !> AUTHOR_VIVIAN_BASS
    !| Finds the indices of the true values in a logical mask.
    pure subroutine which(mask, n, idx_out, m_max, m_out, ierr)
        integer(int32), intent(in) :: m_max
            !! Maximum size of `idx_out`.
        integer(int32), intent(in) :: n
            !! Size of the mask.
        logical, intent(in) :: mask(n)
            !! Logical array of size n.
        integer(int32), intent(out) :: idx_out(m_max)
            !! Integer array to store the indices of true values.
        integer(int32), intent(out) :: m_out
            !! Actual size of `idx_out` (number of true values found).
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i, count

        ! Initialize error code
        call set_ok(ierr)

        ! Validate inputs
        call validate_dimension_size(n, ierr, arg_pos=2_int32)
        call validate_in_range_int(m_max, ierr, min=1_int32, arg_pos=4_int32)

        if (is_err(ierr)) return

        count = 0
        idx_out = 0  ! Initialize to avoid garbage values
        do i = 1, n
            if (mask(i)) then
                count = count + 1
                if (count <= m_max) then
                    idx_out(count) = i
                end if
            end if
        end do
        m_out = count
    end subroutine which

    !> AUTHOR_VIVIAN_BASS
    !| Performs LOESS smoothing on a set of data points.
    !| Smooths `y_ref` at `x_query` using reference points `x_ref`, `y_ref`, and kernel parameters.
    !| The user must pre-filter data and provide only valid indices in indices_used.
    pure subroutine loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
                                    kernel_sigma, kernel_cutoff, y_out, ierr)
        integer(int32), intent(in) :: n_total
            !! Total number of reference points.
        integer(int32), intent(in) :: n_target
            !! Number of target points to smooth.
        real(real64), intent(in) :: x_ref(n_total)
            !! Reference x-coordinates.
        real(real64), intent(in) :: y_ref(n_total)
            !! Reference y-coordinates (length n_total).
        integer(int32), intent(in) :: indices_used(n_used)
            !! Indices of reference points used for smoothing (only valid indices).
        integer(int32), intent(in) :: n_used
            !! Number of indices actually used for smoothing.
        real(real64), intent(in) :: x_query(n_target)
            !! Target x-coordinates to smooth.
        real(real64), intent(in) :: kernel_sigma
            !! Bandwidth parameter for the kernel.
        real(real64), intent(in) :: kernel_cutoff
            !! Cutoff for the kernel, not used if zero
        real(real64), intent(out) :: y_out(n_target)
            !! Output smoothed values (length n_target).
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_target, i_used, used_idx
        real(real64) :: query_x, ref_x, delta, sum_weights, weight
        real(real64) :: min_dist
        integer(int32) :: min_idx
        logical :: exact_match_found, use_kernel

        ! Initialize error code
        call set_ok(ierr)

        ! Input validation
        call validate_dimension_size(n_total, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_target, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_used, ierr, arg_pos=6_int32)
        call validate_in_range_real(kernel_sigma, ierr, min=0.0_real64, arg_pos=8_int32)
        call validate_in_range_real(kernel_cutoff, ierr, min=0.0_real64, arg_pos=9_int32)
        call validate_all_in_range_int(indices_used, n_used, ierr, min=1_int32, max=n_total, arg_pos=5_int32)

        if (is_err(ierr)) return

        ! Check if we should use kernel smoothing
        use_kernel = (kernel_sigma > 0.0_real64)

        do concurrent(i_target=1:n_target) &
            local(query_x, sum_weights, min_dist, min_idx, exact_match_found, used_idx, ref_x, delta, i_used, weight) &
            shared(y_out, n_used, indices_used, x_ref, x_query, y_ref, use_kernel, kernel_cutoff, kernel_sigma)

            query_x = x_query(i_target)
            sum_weights = 0.0_real64
            y_out(i_target) = 0.0_real64
            min_dist = huge(1.0_real64)
            min_idx = indices_used(1)
            exact_match_found = .false.

            ! Process all reference points
            do i_used = 1, n_used
                used_idx = indices_used(i_used)
                ref_x = x_ref(used_idx)
                delta = abs(query_x - ref_x)

                ! Check for exact match
                if (delta == 0.0_real64) then
                    y_out(i_target) = y_ref(used_idx)
                    exact_match_found = .true.
                    exit
                end if

                ! Track closest point for potential fallback
                if (delta < min_dist) then
                    min_dist = delta
                    min_idx = used_idx
                end if

                ! Apply kernel smoothing if enabled and within cutoff
                if (use_kernel .and. delta <= kernel_cutoff*kernel_sigma) then
                    weight = exp(-(delta/kernel_sigma)**2)
                    sum_weights = sum_weights + weight
                    y_out(i_target) = y_out(i_target) + weight*y_ref(used_idx)
                end if
            end do

            ! Finalize result if no exact match was found
            if (.not. exact_match_found) then
                if (sum_weights > 0.0_real64) then
                    ! We have weighted average from kernel smoothing
                    y_out(i_target) = y_out(i_target)/sum_weights
                else
                    ! Fallback: use nearest neighbor
                    y_out(i_target) = y_ref(min_idx)
                end if
            end if
        end do
    end subroutine loess_smooth_2d

    !> AUTHOR_JITU_DABA
    !| Compute the Empirical Distribution Function (EDF) from pre-sorted permutation.
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| Assumes `values` is already sorted by `values[perm]`. Caller controls sorting algorithm.
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    pure subroutine compute_edf(values, n_values, perm, unique_values, cdf_values, n_unique, ierr)
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        integer(int32), intent(in) :: perm(n_values)
            !! Pre-sorted permutation indices (must be sorted by values[perm]).
        real(real64), intent(out) :: unique_values(n_values)
            !! Sorted unique data values.
        real(real64), intent(out) :: cdf_values(n_values)
            !! Corresponding cumulative frequencies between 0 and 1.
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Initialize error code and outputs
        call set_ok(ierr)

        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(perm, n_values, ierr, min=1_int32, max=n_values, arg_pos=3_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)

        if (is_err(ierr)) return

        call compute_edf_helper(values, n_values, perm, unique_values, cdf_values, n_unique)
    end subroutine compute_edf

    !> AUTHOR_JITU_DABA
    !| (no input validation) Compute the Empirical Distribution Function (EDF) from pre-sorted permutation.
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| Assumes `values` is already sorted by `values[perm]`. Caller controls sorting algorithm.
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    pure subroutine compute_edf_helper(values, n_values, perm, unique_values, cdf_values, n_unique)
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        integer(int32), intent(in) :: perm(n_values)
            !! Pre-sorted permutation indices (must be sorted by values[perm]).
        real(real64), intent(out) :: unique_values(n_values)
            !! Sorted unique data values.
        real(real64), intent(out) :: cdf_values(n_values)
            !! Corresponding cumulative frequencies between 0 and 1.
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)

        integer(int32) :: i_value
        real(real64) :: current_val, cumulative_count

        unique_values = 0.0_real64
        cdf_values = 0.0_real64

        ! Identify unique values and compute cumulative frequencies
        n_unique = 0
        cumulative_count = 0.0_real64
        current_val = -huge(1.0_real64)  ! Initialize to lowest possible value

        do i_value = 1, n_values
            ! Check if this is a new unique value (exact comparison)
            if (values(perm(i_value)) /= current_val) then
                ! New unique value found
                current_val = values(perm(i_value))
                n_unique = n_unique + 1

                unique_values(n_unique) = current_val
            end if

            ! Update cumulative count and set CDF value
            cumulative_count = cumulative_count + 1.0_real64
            cdf_values(n_unique) = cumulative_count/real(n_values, real64)
        end do
    end subroutine compute_edf_helper

    !> AUTHOR_JITU_DABA
    !| Helper routine that sorts and calls compute_edf.
    !| Allocates workspace internally and performs sorting before computing EDF.
    !| Use this for convenience; use compute_edf directly for custom sorting.
    pure subroutine compute_edf_alloc(values, n_values, unique_values, cdf_values, n_unique, ierr)
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        real(real64), intent(out) :: unique_values(n_values)
            !! Sorted unique data values.
        real(real64), intent(out) :: cdf_values(n_values)
            !! Corresponding cumulative frequencies between 0 and 1.
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Local workspace arrays with explicit size
        integer(int32), dimension(:), allocatable :: perm

        call set_ok(ierr)

        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(perm(n_values))

        call init_perm(perm)

        call sort_array_heapsort(values, perm)

        ! Compute EDF with sorted permutation
        call compute_edf_helper(values, n_values, perm, unique_values, cdf_values, n_unique)
    end subroutine compute_edf_alloc

    !> AUTHOR_AARON_SCHROEDER
    !| Calculate the percentile of an array given a sorted permutation.
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_percentile(array, permutation, percentile, value, ierr)
        real(real64), intent(in) :: array(:)
            !! input array
        integer(int32), intent(in) :: permutation(:)
            !! permutation vector representing sorted order
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)
        real(real64), intent(out) :: value
            !! output percentile value
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: n

        ! Initialize error
        call set_ok(ierr)

        ! Input validation
        n = size(permutation, kind=int32)
        call validate_dimension_size(n, ierr, arg_pos=2_int32)
        call validate_in_range_int(size(array, kind=int32), ierr, min=1_int32, max=n, arg_pos=1_int32)
        call validate_all_in_range_int(permutation, n, ierr, min=1_int32, max=size(array, kind=int32), arg_pos=2_int32)
        call validate_in_range_real(percentile, ierr, min=0.0_real64, max=100.0_real64, arg_pos=3_int32)

        if (is_err(ierr)) return

        call calc_percentile_helper(array, permutation, percentile, value)
    end subroutine calc_percentile

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Calculate the fractional index using linear interpolation method
    pure real(real64) function calc_percentile_rank(percentile, n) result(rank)
        integer(int32), intent(in) :: n
            !! Sample size
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)

        rank = (percentile/100.0_real64)*real(n - 1, real64) + 1.0_real64
    end function calc_percentile_rank

    !> AUTHOR_AARON_SCHROEDER
    !| (no input validation) Calculate the percentile of an array given a sorted permutation.
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_percentile_helper(array, permutation, percentile, value)
        real(real64), intent(in) :: array(:)
            !! input array
        integer(int32), intent(in) :: permutation(:)
            !! permutation vector representing sorted order
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)
        real(real64), intent(out) :: value
            !! output percentile value

        integer(int32) :: n, lower_index
        real(real64) :: index, fraction, lower_value, upper_value

        n = size(permutation, kind=int32)

        ! Handle single element case
        if (size(array, kind=int32) == 1) then
            value = array(1)
            return
        end if

        index = calc_percentile_rank(percentile, n)
        lower_index = floor(index)
        fraction = index - real(lower_index, real64)

        ! Handle edge cases for indices
        if (lower_index < 1) then
            value = array(permutation(1))  ! Smallest value in sorted order
        else if (lower_index >= n) then
            value = array(permutation(n))  ! Largest value in sorted order
        else
            ! Linear interpolation between adjacent values using permuted indices
            lower_value = array(permutation(lower_index))
            upper_value = array(permutation(lower_index + 1))
            value = lower_value + fraction*(upper_value - lower_value)
        end if
    end subroutine calc_percentile_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Calculate the percentile of an array, allocating necessary arrays when no sorting permutation is given
    !| @note This subroutine uses quicksort internally which may cause a spike in memory usage for large arrays.
    pure subroutine calc_percentile_alloc(array, percentile, value, ierr)
        real(real64), intent(in) :: array(:)
            !! Input array
        real(real64), intent(in) :: percentile
            !! Desired percentile (0-100)
        real(real64), intent(out) :: value
            !! Output percentile value
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: n
        integer(int32), allocatable :: perm(:)

        n = size(array, kind=int32)
        ! Initialize error
        call set_ok(ierr)

        call validate_dimension_size(n, ierr, arg_pos=1_int32)
        call validate_in_range_real(percentile, ierr, min=0.0_real64, max=100.0_real64, arg_pos=2_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(perm(n))

        call init_perm(perm)

        ! Sort the array indirectly
        call sort_real_heapsort(array, perm)

        ! Calculate percentile using sorted permutation
        call calc_percentile_helper(array, perm, percentile, value)
    end subroutine calc_percentile_alloc

    !> AUTHOR_VIVIAN_BASS
    !| Calculate empirical p-values for scaled expression distances (RDI).
    !|
    !| Implements:
    !|   P(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    !|
    !| Because distances are non-negative, a one-sided upper-tail empirical p-value is used.
    !|
    !| Assumptions / preconditions:
    !| - sorted_rdi(1:n_genes) contains the empirical distribution D.
    !| - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    pure subroutine compute_empirical_p_values(n_genes, rdi, sorted_rdi, perm, p_values, c_const)
        integer(int32), intent(in) :: n_genes
            !! Number of genes being processed.
        real(real64), intent(in) :: rdi(n_genes)
            !! empirical distribution D
        real(real64), intent(in) :: sorted_rdi(n_genes)
            !! empirical distribution D with non negative values
        real(real64), intent(out) :: p_values(n_genes)
            !! Output array to store the computed p-values for each gene.
        real(real64), intent(in) :: c_const
            !! Constant used in the computation, typically 1
        integer(int32), intent(in) :: perm(n_genes)
            !! Permutation array with sorted indices for sorted_rdi

        integer(int32) :: i, first_ge, count_ge
        real(real64) :: denom, d

        denom = real(n_genes, real64) + c_const
        if (denom <= 0.0_real64) then
            p_values = 1.0_real64
            return
        end if

        do i = 1, n_genes
            d = rdi(i)

            ! Invalid / negative => not an outlier: p = 1
            if (d < 0.0_real64) then
                p_values(i) = 1.0_real64
                cycle
            end if

            first_ge = lower_bound_ge(sorted_rdi, perm, n_genes, d)

            if (first_ge <= n_genes) then
                count_ge = n_genes - first_ge + 1_int32
            else
                count_ge = 0_int32
            end if

            p_values(i) = (real(count_ge, real64) + c_const)/denom
        end do

    end subroutine compute_empirical_p_values

    !> AUTHOR_VIVIAN_BASS
    !| First position pos in [1..n] such that sorted_rdi(perm(pos)) >= x. Returns n+1 if none.
    pure integer(int32) function lower_bound_ge(vals, p, n, x) result(pos)
        real(real64), intent(in) :: vals(n)
            !! Input array of values to be searched
        integer(int32), intent(in) :: p(n)
            !! Permutation array with sorted indices
        integer(int32), intent(in) :: n
            !! Number of elements in the `vals` and `p` arrays
        real(real64), intent(in) :: x
            !! Input value to be searched for or compared against within `vals`

        integer(int32) :: l, h, mid

        l = 1_int32
        h = n
        pos = n + 1_int32

        do while (l <= h)
            mid = l + (h - l)/2_int32
            if (vals(p(mid)) >= x) then
                pos = mid
                h = mid - 1_int32
            else
                l = mid + 1_int32
            end if
        end do
    end function lower_bound_ge

end module f42_utils

! === C WRAPPERS ===

!> C wrapper for which.
!| Converts integer mask to logical and calls which.
pure subroutine which_c(mask, n, idx_out, m_max, m_out, ierr) bind(C, name="which_c")
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_fortran_env, only: int32
    use f42_utils, only: which
    use tox_conversions, only: c_int_as_logical
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n
        !! Size of the mask.
    integer(c_int), intent(in), target :: m_max
        !! Maximum size of idx_out.
    integer(c_int), intent(in), target :: mask(n)
        !! Integer mask array (0/1 values).
    integer(c_int), intent(out), target :: idx_out(m_max)
        !! Output array for indices of true values.
    integer(c_int), intent(out), target :: m_out
        !! Actual size of idx_out (number of true values found).
    integer(c_int), intent(out), target :: ierr
        !! Error code: 0=ok, 201=invalid input, 202=empty input
    logical :: mask_f(n)
    integer(int32) :: ierr_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n)
    M_CHECK_NON_NULL(m_max)
    M_CHECK_NON_NULL(mask)
    M_CHECK_NON_NULL(idx_out)

    ! Use tox_conversions utility for c_int to logical conversion
    call c_int_as_logical(mask, mask_f)
    call which(mask_f, n, idx_out, m_max, m_out, ierr_f)
    ierr = ierr_f
end subroutine which_c

!> C wrapper for loess_smooth_2d.
!| Direct wrapper - user must pre-filter indices in C before calling.
pure subroutine loess_smooth_2d_c(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
                                  kernel_sigma, kernel_cutoff, y_out, ierr) bind(C, name="loess_smooth_2d_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int32
    use f42_utils, only: loess_smooth_2d
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_total
        !! Total number of reference points.
    integer(c_int), intent(in), target :: n_target
        !! Number of target points to smooth.
    real(c_double), intent(in), target :: x_ref(n_total)
        !! Reference x-coordinates.
    real(c_double), intent(in), target :: y_ref(n_total)
        !! Reference y-coordinates (length n_total).
    integer(c_int), intent(in), target :: indices_used(n_used)
        !! Indices of reference points used for smoothing (pre-filtered).
    integer(c_int), intent(in), target :: n_used
        !! Number of indices actually used for smoothing.
    real(c_double), intent(in), target :: x_query(n_target)
        !! Target x-coordinates to smooth.
    real(c_double), intent(in), target :: kernel_sigma
        !! Bandwidth parameter for the kernel.
    real(c_double), intent(in), target :: kernel_cutoff
        !! Cutoff for the kernel.
    real(c_double), intent(out), target :: y_out(n_target)
        !! Output smoothed values (length n_target).
    integer(c_int), intent(out), target :: ierr
        !! Error code: 0=ok, 201=invalid input, 202=empty input

    integer(int32) :: ierr_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_total)
    M_CHECK_NON_NULL(n_target)
    M_CHECK_NON_NULL(x_ref)
    M_CHECK_NON_NULL(y_ref)
    M_CHECK_NON_NULL(indices_used)
    M_CHECK_NON_NULL(n_used)
    M_CHECK_NON_NULL(x_query)
    M_CHECK_NON_NULL(kernel_sigma)
    M_CHECK_NON_NULL(kernel_cutoff)
    M_CHECK_NON_NULL(y_out)

    call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
                         kernel_sigma, kernel_cutoff, y_out, ierr_f)
    ierr = ierr_f

end subroutine loess_smooth_2d_c

!> C wrapper for compute_edf.
!| Allocates workspace internally and exposes a simple interface with C types.
!| The number of unique values is returned via n_unique output parameter.
subroutine compute_edf_c(values, n_values, unique_values, cdf_values, n_unique, ierr) &
    bind(C, name="compute_edf_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_utils, only: compute_edf_alloc
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_values
        !! Number of values in the input array.
    real(c_double), intent(in), target :: values(n_values)
        !! Array of observed data values (e.g., contributions or spikes).
    real(c_double), intent(out), target :: unique_values(n_values)
        !! Sorted unique data values (sized to n_values).
    real(c_double), intent(out), target :: cdf_values(n_values)
        !! Corresponding cumulative frequencies between 0 and 1 (sized to n_values).
    integer(c_int), intent(out), target :: n_unique
        !! Number of unique values found.
    integer(c_int), intent(out), target :: ierr
        !! Error code: 0=ok, 201=invalid input, 202=empty input

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_values)
    M_CHECK_NON_NULL(values)
    M_CHECK_NON_NULL(unique_values)
    M_CHECK_NON_NULL(cdf_values)

    call compute_edf_alloc(values, n_values, unique_values, cdf_values, n_unique, ierr)
end subroutine compute_edf_c

!> Expert C wrapper for compute_edf.
!| Direct interface to compute_edf for users who have already sorted their data
!| or have a custom permutation vector. This skips the internal sorting step
!| for better performance when the caller has full control.
subroutine compute_edf_expert_c(values, n_values, perm, unique_values, cdf_values, n_unique, ierr) &
    bind(C, name="compute_edf_expert_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_utils, only: compute_edf
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_values
        !! Number of values in the input array.
    real(c_double), intent(in), target :: values(n_values)
        !! Array of observed data values (e.g., contributions or spikes).
    integer(c_int), intent(in), target :: perm(n_values)
        !! Pre-sorted permutation indices (must be sorted by values[perm]).
        !! Caller is responsible for sorting this array before calling.
    real(c_double), intent(out), target :: unique_values(n_values)
        !! Sorted unique data values (sized to n_values).
    real(c_double), intent(out), target :: cdf_values(n_values)
        !! Corresponding cumulative frequencies between 0 and 1 (sized to n_values).
    integer(c_int), intent(out), target :: n_unique
        !! Number of unique values found.
    integer(c_int), intent(out), target :: ierr
        !! Error code: 0=ok, 201=invalid input, 202=empty input

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_values)
    M_CHECK_NON_NULL(values)
    M_CHECK_NON_NULL(perm)
    M_CHECK_NON_NULL(unique_values)
    M_CHECK_NON_NULL(cdf_values)

    call compute_edf(values, n_values, perm, unique_values, cdf_values, n_unique, ierr)
end subroutine compute_edf_expert_c

!> C-compatible wrapper for compute_empirical_p_values
subroutine compute_empirical_p_values_c(n_genes, rdi, sorted_rdi, perm, p_values, c_const) bind(C, name="empirical_p_values_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use f42_utils, only: compute_empirical_p_values
    implicit none

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes being processed.
    real(c_double), intent(in), target :: rdi(n_genes)
        !! empirical distribution D
    real(c_double), intent(in), target :: sorted_rdi(n_genes)
        !! empirical distribution D with non negative values
    integer(c_int), intent(in), target :: perm(n_genes)
        !! Permutation array with sorted indices for sorted_rdi
    real(c_double), intent(out), target :: p_values(n_genes)
        !! Output array to store the computed p-values for each gene.
    real(c_double), intent(in), target :: c_const
        !! Constant used in the computation, typically 1

    call compute_empirical_p_values(n_genes, rdi, sorted_rdi, perm, p_values, c_const)

end subroutine compute_empirical_p_values_c
