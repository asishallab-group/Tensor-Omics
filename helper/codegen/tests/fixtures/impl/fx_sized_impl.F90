#include <src/macros.h>

!> summary: Implementation fixture for consumers of producers that are generated wrappers
!| author: A Developer
!| Each argument sized by a producer names a generated wrapper in `fx_producers` (from
!| [[fx_producers_impl(module)]]). Where the size only sizes a work array, the generated Fortran
!| wrapper calls the producer itself; where it sizes a result, Python and R call the producer's
!| binding. Either way a producer's error names the argument the caller passed. Compiled and run
!| by `test_end_to_end_impl.py`.
module fx_sized_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: fx_fill_impl, fx_work_total_impl, fx_two_producers_impl, fx_extent_total_impl
    public :: fx_first_of_two_impl, fx_status_per_vector_impl, fx_info_impl, fx_both_firsts_impl
    public :: fx_two_scales_impl, fx_norm_rank_impl, fx_large_total_impl

contains

    !> summary: A result sized by a producer the bindings call, with a renamed input
    !| author: A Developer
    pure subroutine fx_fill_impl(n_items, x, scale, out_size, filled)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: scale
            !! how many results per value
        integer(int32), intent(in) :: out_size
            !! size of the result
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | factor         | scale       |
        real(real64), dimension(out_size), intent(out) :: filled
            !! `scale`, every element
        filled = real(scale, real64) + 0.0_real64*sum(x)
    end subroutine fx_fill_impl

    !> summary: A work array sized by a producer the generated wrapper calls, with a renamed input
    !| author: A Developer
    pure subroutine fx_work_total_impl(n_items, x, scale, work_size, tmp_buf, total)
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
    end subroutine fx_work_total_impl

    !> summary: Two work arrays, each sized by its own producer; the first producer may refuse
    !| author: A Developer
    pure subroutine fx_two_producers_impl(n_items, x, n_count, scale, first_size, tmp_first, &
                                          second_size, tmp_second, total)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: n_count
            !! what the first scratch is sized by
        integer(int32), intent(in) :: scale
            !! what the second scratch is sized by, per value
        integer(int32), intent(in) :: first_size
            !! size of the first scratch
            !! DM_OUTPUT_FROM(capped_size, fx_capped_size, fx_producers, AUTO)
        real(real64), dimension(first_size), intent(out) :: tmp_first
            !! scratch
        integer(int32), intent(in) :: second_size
            !! size of the second scratch
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | factor         | scale       |
        real(real64), dimension(second_size), intent(out) :: tmp_second
            !! scratch
        real(real64), intent(out) :: total
            !! the sum of `x`, plus the sizes of both scratches
        tmp_first = 1.0_real64
        tmp_second = 1.0_real64
        total = sum(x) + sum(tmp_first) + sum(tmp_second)
    end subroutine fx_two_producers_impl

    !> summary: A work array sized by a producer whose extent the consumer calls by another name
    !| author: A Developer
    !| The producer's binding reads its `n` off `x`, so Python and R pass it only `x`; the
    !| generated Fortran wrapper calls the producer directly and has to pass `n` as well -- the
    !| table says it is this procedure's `m`.
    pure subroutine fx_extent_total_impl(m, x, extent_size, tmp_buf, total)
        integer(int32), intent(in) :: m
            !! number of values
        real(real64), dimension(m), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: extent_size
            !! size of the scratch
            !! DM_OUTPUT_FROM(extent_size, fx_extent_size, fx_producers, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n              | m           |
        real(real64), dimension(extent_size), intent(out) :: tmp_buf
            !! scratch
        real(real64), intent(out) :: total
            !! the sum of `x`, plus the size of the scratch
        tmp_buf = 1.0_real64
        total = sum(x) + sum(tmp_buf)
    end subroutine fx_extent_total_impl

    !> summary: A work array sized by one of a producer's two outputs; the other goes unused
    !| author: A Developer
    pure subroutine fx_first_of_two_impl(n_items, x, first_size, tmp_buf, total)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: first_size
            !! size of the scratch
            !! DM_OUTPUT_FROM(first_size, fx_two_sizes, fx_producers, AUTO)
        real(real64), dimension(first_size), intent(out) :: tmp_buf
            !! scratch
        real(real64), intent(out) :: total
            !! the sum of `x`, plus the size of the scratch
        tmp_buf = 1.0_real64
        total = sum(x) + sum(tmp_buf)
    end subroutine fx_first_of_two_impl

    !> summary: Two work arrays sized by two producers, whose dropped outputs share a name
    !| author: A Developer
    !| Each producer's `second_size` goes unused -- an integer from one, a real from the other --
    !| and so does `fx_flagged_size`'s `logical(c_bool)` flag.
    pure subroutine fx_both_firsts_impl(n_items, x, a_size, tmp_a, b_size, tmp_b, total)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: a_size
            !! size of the first scratch
            !! DM_OUTPUT_FROM(first_size, fx_two_sizes, fx_producers, AUTO)
        real(real64), dimension(a_size), intent(out) :: tmp_a
            !! scratch
        integer(int32), intent(in) :: b_size
            !! size of the second scratch
            !! DM_OUTPUT_FROM(first_size, fx_flagged_size, fx_producers, AUTO)
        real(real64), dimension(b_size), intent(out) :: tmp_b
            !! scratch
        real(real64), intent(out) :: total
            !! the sum of `x`, plus the sizes of both scratches
        tmp_a = 1.0_real64
        tmp_b = 1.0_real64
        total = sum(x) + sum(tmp_a) + sum(tmp_b)
    end subroutine fx_both_firsts_impl

    !> summary: Two work arrays sized by one producer, given a different input for each
    !| author: A Developer
    !| Two calls of `fx_scaled_size`, one per table, each into its own size.
    pure subroutine fx_two_scales_impl(n_items, x, scale_a, scale_b, size_a, tmp_a, size_b, tmp_b, total)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: scale_a
            !! slots of the first scratch per value
        integer(int32), intent(in) :: scale_b
            !! slots of the second scratch per value
        integer(int32), intent(in) :: size_a
            !! size of the first scratch
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | factor         | scale_a     |
        real(real64), dimension(size_a), intent(out) :: tmp_a
            !! scratch
        integer(int32), intent(in) :: size_b
            !! size of the second scratch
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | factor         | scale_b     |
        real(real64), dimension(size_b), intent(out) :: tmp_b
            !! scratch
        real(real64), intent(out) :: total
            !! the sum of `x`, plus the first scratch's size, plus a thousand times the second's
        tmp_a = 1.0_real64
        tmp_b = 1000.0_real64
        total = sum(x) + sum(tmp_a) + sum(tmp_b)
    end subroutine fx_two_scales_impl

    !> summary: The column norms in ascending order, beside scratch sized by a second producer
    !| author: A Developer
    !| `norms` is an array a producer fills, which `norms_perm` orders; `work_size` is a scalar
    !| from another producer.
    pure subroutine fx_norm_rank_impl(n_dims, n_items, x, scale, norms, norms_perm, work_size, &
                                      tmp_buf, ranked, total)
        integer(int32), intent(in) :: n_dims
            !! number of rows
        integer(int32), intent(in) :: n_items
            !! number of columns
        real(real64), dimension(n_dims, n_items), intent(in) :: x
            !! the matrix
        integer(int32), intent(in) :: scale
            !! slots of scratch per column
        real(real64), dimension(n_items), intent(in) :: norms
            !! the norm of each column
            !! DM_OUTPUT_FROM(norms, fx_column_norms, fx_producers, AUTO)
        integer(int32), dimension(n_items), intent(in) :: norms_perm
            !! the order that sorts `norms` ascending
        integer(int32), intent(in) :: work_size
            !! size of the scratch
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | factor         | scale       |
        real(real64), dimension(work_size), intent(out) :: tmp_buf
            !! scratch
        real(real64), dimension(n_items), intent(out) :: ranked
            !! `norms` in ascending order
        real(real64), intent(out) :: total
            !! the size of the scratch
        tmp_buf = 1.0_real64
        ranked = norms(norms_perm)
        total = sum(tmp_buf)
    end subroutine fx_norm_rank_impl

    !> summary: The sum of the column norms, plus a thousand for each norm above one
    !| author: A Developer
    !| One producer call gives both `norms` and the size `n_large` of the scratch, which depends
    !| on their values. Any number of columns, none included.
    pure subroutine fx_large_total_impl(n_dims, n_items, x, norms, n_large, tmp_large, total)
        integer(int32), intent(in) :: n_dims
            !! number of rows
        integer(int32), intent(in) :: n_items
            !! number of columns
            !! DM_MIN(0_int32)
        real(real64), dimension(n_dims, n_items), intent(in) :: x
            !! the matrix
        real(real64), dimension(n_items), intent(in) :: norms
            !! the norm of each column
            !! DM_OUTPUT_FROM(norms, fx_large_norms, fx_producers, AUTO)
        integer(int32), intent(in) :: n_large
            !! how many norms exceed one
            !! DM_OUTPUT_FROM(n_large, fx_large_norms, fx_producers, AUTO)
        real(real64), dimension(n_large), intent(out) :: tmp_large
            !! scratch, one slot per large norm
        real(real64), intent(out) :: total
            !! the sum of `norms`, plus a thousand per large norm
        tmp_large = 1000.0_real64
        total = sum(norms) + sum(tmp_large)
    end subroutine fx_large_total_impl

    !> summary: A result sized by a producer that takes the shape of an input
    !| author: A Developer
    pure subroutine fx_status_per_vector_impl(n_elements, arr, n_dims, arr_shape, n_vectors, status)
        integer(int32), intent(in) :: n_elements
            !! number of elements
        real(real64), dimension(n_elements), intent(in) :: arr
            !! the array, flattened
        integer(int32), intent(in) :: n_dims
            !! its rank
        integer(int32), dimension(n_dims), intent(in) :: arr_shape
            !! its shape
        integer(int32), intent(in) :: n_vectors
            !! how many vectors it holds
            !! DM_OUTPUT_FROM(n_vectors, fx_vector_count, fx_producers, AUTO)
        integer(int32), dimension(n_vectors), intent(out) :: status
            !! the length of each vector
        status = arr_shape(1) + 0_int32*int(sum(arr), int32)
    end subroutine fx_status_per_vector_impl

    !> summary: A result whose size the caller passes, pointed at its producer by `JUST_INFO`
    !| author: A Developer
    pure subroutine fx_info_impl(n_items, x, out_size, ones)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), intent(in) :: out_size
            !! size of the result
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, JUST_INFO)
        real(real64), dimension(out_size), intent(out) :: ones
            !! ones
        ones = 1.0_real64 + 0.0_real64*sum(x)
    end subroutine fx_info_impl
end module fx_sized_impl
