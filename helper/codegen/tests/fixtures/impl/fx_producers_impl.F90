#include <src/macros.h>

!> summary: Implementation fixture for producers that are generated wrappers
!| author: A Developer
!| Each implementation here sizes something a consumer takes, and the consumer names its
!| generated wrapper -- the plain-named one, in `fx_producers` -- as the producer. The consumers
!| are in [[fx_sized_impl(module)]], and one is here, beside the producer it calls: a generated
!| module never `use`s itself. Compiled and run by `test_end_to_end_impl.py`.
module fx_producers_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: fx_scaled_size_impl, fx_capped_size_impl, fx_extent_size_impl, fx_two_sizes_impl
    public :: fx_vector_count_impl, fx_local_total_impl

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
