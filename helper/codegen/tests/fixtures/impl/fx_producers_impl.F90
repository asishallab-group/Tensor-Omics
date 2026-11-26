#include <src/macros.h>

!> summary: Implementation fixture for producers that are generated wrappers
!| author: A Developer
!| Each implementation here sizes something a consumer takes, and the consumer names its
!| generated wrapper -- the plain-named one, in `fx_producers` -- as the producer. The consumers
!| are in [[fx_sized_impl(module)]], and one is here, beside the producer it calls: a generated
!| module never `use`s itself. Compiled and run by `test_end_to_end_impl.py`.
module fx_producers_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    M_IMPLICIT_NONE
    private

    public :: fx_scaled_size_impl, fx_capped_size_impl, fx_extent_size_impl, fx_two_sizes_impl
    public :: fx_vector_count_impl, fx_local_total_impl, fx_flagged_size_impl
    public :: fx_column_norms_impl, fx_large_norms_impl, fx_padded_size_impl, fx_rounded_impl
    public :: fx_strided_size_impl, fx_method_size_impl, fx_any_method_size_impl
    public :: METHOD_ONE, METHOD_TWO, METHOD_THREE

    !> the first method: one slot per item
    integer(int32), parameter :: METHOD_ONE = 1_int32
    !> the second method: two slots per item
    integer(int32), parameter :: METHOD_TWO = 2_int32
    !> the third method: three slots per item
    integer(int32), parameter :: METHOD_THREE = 3_int32

contains

    !> summary: A size for a number of items and a factor
    !| author: A Developer
    !| Their product.
    pure subroutine fx_scaled_size_impl(n_items, factor, scaled_size)
        integer(int32), intent(in) :: n_items
            !! number of items
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: factor
            !! multiplier
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: scaled_size
            !! `n_items` times `factor`
        scaled_size = n_items*factor
    end subroutine fx_scaled_size_impl

    !> summary: A size one larger than a count that may be at most five
    !| author: A Developer
    pure subroutine fx_capped_size_impl(n_count, capped_size)
        integer(int32), intent(in) :: n_count
            !! the count
            !! DM_MAX(5_int32)
        integer(int32), intent(out) :: capped_size
            !! one more than `n_count`
        capped_size = n_count + 1_int32
    end subroutine fx_capped_size_impl

    !> summary: A size one larger than the number of values, which its wrapper reads off them
    !| author: A Developer
    pure subroutine fx_extent_size_impl(n, x, extent_size)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: x
            !! the values
        integer(int32), intent(out) :: extent_size
            !! one more than `n`
        extent_size = n + 1_int32 + 0_int32*int(sum(x), int32)
    end subroutine fx_extent_size_impl

    !> summary: Two sizes for a number of items: the number itself, and twice it
    !| author: A Developer
    pure subroutine fx_two_sizes_impl(n_items, first_size, second_size)
        integer(int32), intent(in) :: n_items
            !! number of items
        integer(int32), intent(out) :: first_size
            !! `n_items`
        integer(int32), intent(out) :: second_size
            !! twice `n_items`
        first_size = n_items
        second_size = 2_int32*n_items
    end subroutine fx_two_sizes_impl

    !> summary: A size, a weight and a flag: two outputs of kinds an integer consumer never names
    !| author: A Developer
    !| Its `second_size` is a real, where `fx_two_sizes`' is an integer.
    pure subroutine fx_flagged_size_impl(n_items, first_size, second_size, flag)
        integer(int32), intent(in) :: n_items
            !! number of items
        integer(int32), intent(out) :: first_size
            !! `n_items` plus one
        real(real64), intent(out) :: second_size
            !! half of `n_items`
        logical(c_bool), intent(out) :: flag
            !! whether there are any
        first_size = n_items + 1_int32
        second_size = 0.5_real64*n_items
        flag = n_items > 0
    end subroutine fx_flagged_size_impl

    !> summary: The Euclidean norm of each column of a matrix of at most three rows
    !| author: A Developer
    !| An array producer. Its bound on `n_dims` is one its consumers leave out, so its own error
    !| reaches their callers.
    pure subroutine fx_column_norms_impl(n_dims, n_items, x, norms)
        integer(int32), intent(in) :: n_dims
            !! number of rows
            !! DM_MIN(1_int32)
            !! DM_MAX(3_int32)
        integer(int32), intent(in) :: n_items
            !! number of columns
            !! DM_MIN(0_int32)
        real(real64), dimension(n_dims, n_items), intent(in) :: x
            !! the matrix
        real(real64), dimension(n_items), intent(out) :: norms
            !! the norm of each column
        norms = sqrt(sum(x**2, dim=1))
    end subroutine fx_column_norms_impl

    !> summary: The norm of each column of a matrix, and how many of them exceed one
    !| author: A Developer
    !| One call gives an array and a size that depends on its values.
    pure subroutine fx_large_norms_impl(n_dims, n_items, x, norms, n_large)
        integer(int32), intent(in) :: n_dims
            !! number of rows
            !! DM_MIN(1_int32)
            !! DM_MAX(3_int32)
        integer(int32), intent(in) :: n_items
            !! number of columns
            !! DM_MIN(0_int32)
        real(real64), dimension(n_dims, n_items), intent(in) :: x
            !! the matrix
        real(real64), dimension(n_items), intent(out) :: norms
            !! the norm of each column
        integer(int32), intent(out) :: n_large
            !! how many of `norms` exceed one
        norms = sqrt(sum(x**2, dim=1))
        n_large = count(norms > 1.0_real64)
    end subroutine fx_large_norms_impl

    !> summary: Each value rounded to the nearest integer
    !| author: A Developer
    pure subroutine fx_rounded_impl(n_items, x, rounded)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
            !! DM_MIN(-1000.0_real64)
            !! DM_MAX(1000.0_real64)
        integer(int32), dimension(n_items), intent(out) :: rounded
            !! `x`, rounded
        rounded = nint(x, int32)
    end subroutine fx_rounded_impl

    !> summary: A size of a number of items, plus a padding when one is given
    !| author: A Developer
    pure subroutine fx_padded_size_impl(n_items, padding, padded_size)
        integer(int32), intent(in) :: n_items
            !! number of items
        integer(int32), intent(in), optional :: padding
            !! extra slots, none when it is left out
            !! DM_MIN(0_int32)
        integer(int32), intent(out) :: padded_size
            !! `n_items`, plus `padding`
        padded_size = n_items
        if (present(padding)) padded_size = padded_size + padding
    end subroutine fx_padded_size_impl

    !> summary: A size of a number of items, times a stride of three unless one is given
    !| author: A Developer
    pure subroutine fx_strided_size_impl(n_items, stride, strided_size)
        integer(int32), intent(in) :: n_items
            !! number of items
        integer(int32), intent(in), optional :: stride
            !! slots per item
            !! DM_DEFAULT(3_int32)
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: strided_size
            !! `n_items` times `stride`
        integer(int32) :: slots
        M_DEFAULT_VAL(stride, slots, 3_int32)
        strided_size = n_items*slots
    end subroutine fx_strided_size_impl

    !> summary: A size of a number of items, times the code of a method
    !| author: A Developer
    pure subroutine fx_method_size_impl(n_items, method, method_size)
        integer(int32), intent(in) :: n_items
            !! number of items
        integer(int32), intent(in) :: method
            !! how many slots per item
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | one slot | [[fx_producers_impl(module):METHOD_ONE(variable)]] |
            !! | two slots | [[fx_producers_impl(module):METHOD_TWO(variable)]] |
            !! | three slots | [[fx_producers_impl(module):METHOD_THREE(variable)]] |
        integer(int32), intent(out) :: method_size
            !! `n_items` times the method's code
        method_size = n_items*method
    end subroutine fx_method_size_impl

    !> summary: A size of a number of items, times the code of a method, the second unless given
    !| author: A Developer
    pure subroutine fx_any_method_size_impl(n_items, method, method_size)
        integer(int32), intent(in) :: n_items
            !! number of items
        integer(int32), intent(in), optional :: method
            !! how many slots per item
            !! DM_DEFAULT(METHOD_TWO)
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | one slot | [[fx_producers_impl(module):METHOD_ONE(variable)]] |
            !! | two slots | [[fx_producers_impl(module):METHOD_TWO(variable)]] |
            !! | three slots | [[fx_producers_impl(module):METHOD_THREE(variable)]] |
        integer(int32), intent(out) :: method_size
            !! `n_items` times the method's code
        integer(int32) :: code
        M_DEFAULT_VAL(method, code, METHOD_TWO)
        method_size = n_items*code
    end subroutine fx_any_method_size_impl

    !> summary: The number of vectors in an array of a given shape: every extent but the first
    !| author: A Developer
    pure subroutine fx_vector_count_impl(n_dims, arr_shape, n_vectors)
        integer(int32), intent(in) :: n_dims
            !! the rank
        integer(int32), dimension(n_dims), intent(in) :: arr_shape
            !! the shape
        integer(int32), intent(out) :: n_vectors
            !! the product of every extent but the first
        n_vectors = product(arr_shape(2:))
    end subroutine fx_vector_count_impl

    !> summary: The sum of the values, with a work array sized by a producer in this module
    !| author: A Developer
    pure subroutine fx_local_total_impl(n_items, x, scale, work_size, tmp_buf, total)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: scale
            !! how many slots of scratch per value
        integer(int32), intent(in) :: work_size
            !! size of the scratch
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | factor         | scale       |
        real(real64), dimension(work_size), intent(out) :: tmp_buf
            !! scratch
        real(real64), intent(out) :: total
            !! the sum of `x`, plus the size of the scratch
        tmp_buf = 1.0_real64
        total = sum(x) + sum(tmp_buf)
    end subroutine fx_local_total_impl
end module fx_producers_impl
