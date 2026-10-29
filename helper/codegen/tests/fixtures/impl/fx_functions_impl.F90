#include <src/macros.h>

!> summary: Implementation fixture for functions, whose result every tier publishes before `ierr`
!| author: A Developer
!| A scalar, an array and a `logical(c_bool)` result, an `elemental` function, one with its own
!| `ierr`, and one with a work array, so that both tiers carry the result; and a `complex(real64)`
!| result, beside subroutines with a `complex(real64)` output, input, optional input and
!| `intent(inout)` argument. Compiled and run by
!| `test_end_to_end_fortran.py`, and through C, Python and R by `test_end_to_end_impl.py`.
module fx_functions_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err, ERR_DIVISION_BY_ZERO
    M_IMPLICIT_NONE
    private

    public :: fx_total_impl, fx_doubled_impl, fx_all_positive_impl, fx_halve_impl
    public :: fx_share_impl, fx_smallest_impl, fx_complex_impl, fx_make_complex_impl
    public :: fx_parts_impl, fx_offset_impl, fx_twice_impl

contains

    !> summary: The sum of the values
    !| author: A Developer
    !| An edge case: a scalar function result, published as `total` right before `ierr`.
    pure function fx_total_impl(n_values, values) result(total)
        integer(int32), intent(in) :: n_values
            !! number of values
        real(real64), intent(in) :: values(n_values)
            !! the values
        real(real64) :: total
            !! their sum
        total = sum(values)
    end function fx_total_impl

    !> summary: Every value doubled
    !| author: A Developer
    !| An edge case: an array result sized by an argument.
    pure function fx_doubled_impl(n_values, values) result(doubled)
        integer(int32), intent(in) :: n_values
            !! number of values
        real(real64), intent(in) :: values(n_values)
            !! the values
        real(real64) :: doubled(n_values)
            !! twice each value
        doubled = 2.0_real64*values
    end function fx_doubled_impl

    !> summary: Whether every value is positive
    !| author: A Developer
    !| An edge case: a `logical(c_bool)` result.
    pure function fx_all_positive_impl(n_values, values) result(all_positive)
        integer(int32), intent(in) :: n_values
            !! number of values
        real(real64), intent(in) :: values(n_values)
            !! the values
        logical(c_bool) :: all_positive
            !! whether every value is greater than zero
        all_positive = logical(all(values > 0.0_real64), c_bool)
    end function fx_all_positive_impl

    !> summary: Half of a value
    !| author: A Developer
    !| An edge case: an `elemental` function, wrapped like an elemental subroutine.
    elemental function fx_halve_impl(value) result(half)
        real(real64), intent(in) :: value
            !! the value
        real(real64) :: half
            !! half of it
        half = value/2.0_real64
    end function fx_halve_impl

    !> summary: The first value's share of the total
    !| author: A Developer
    !| An edge case: a function with its own `ierr`, after which nothing may come -- the
    !| result goes before it.
    function fx_share_impl(n_values, values, ierr) result(share)
        integer(int32), intent(in) :: n_values
            !! number of values
        real(real64), intent(in) :: values(n_values)
            !! the values
        integer(int32), intent(out) :: ierr
            !! error code; `ERR_DIVISION_BY_ZERO` when the values sum to zero
        real(real64) :: share
            !! the first value divided by the sum of all
        call set_ok(ierr)
        share = 0.0_real64
        if (sum(values) == 0.0_real64) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if
        share = values(1)/sum(values)
    end function fx_share_impl

    !> summary: The smallest value, found in a work array
    !| author: A Developer
    !| An edge case: a work array, so both tiers exist and both publish the result. Not `pure`,
    !| which a function with an `intent(out)` argument cannot be.
    function fx_smallest_impl(n_values, values, tmp_copy) result(smallest)
        integer(int32), intent(in) :: n_values
            !! number of values
        real(real64), intent(in) :: values(n_values)
            !! the values
        real(real64), intent(out) :: tmp_copy(n_values)
            !! scratch space
        real(real64) :: smallest
            !! the smallest value
        tmp_copy = values
        smallest = minval(tmp_copy)
    end function fx_smallest_impl

    !> summary: A complex number from its parts
    !| author: A Developer
    !| An edge case: a `complex(real64)` result, which is two doubles to C.
    pure function fx_complex_impl(re, im) result(z)
        real(real64), intent(in) :: re
            !! the real part
        real(real64), intent(in) :: im
            !! the imaginary part
        complex(real64) :: z
            !! `re + im i`
        z = cmplx(re, im, kind=real64)
    end function fx_complex_impl

    !> summary: A complex number from its parts, as an output argument
    !| author: A Developer
    !| An edge case: a subroutine's `complex(real64)` scalar output, the twin of `fx_complex_impl`.
    pure subroutine fx_make_complex_impl(re, im, z)
        real(real64), intent(in) :: re
            !! the real part
        real(real64), intent(in) :: im
            !! the imaginary part
        complex(real64), intent(out) :: z
            !! `re + im i`
        z = cmplx(re, im, kind=real64)
    end subroutine fx_make_complex_impl

    !> summary: The parts of a complex number
    !| author: A Developer
    !| An edge case: a `complex(real64)` scalar input, which is two doubles to C.
    pure subroutine fx_parts_impl(z, re, im)
        complex(real64), intent(in) :: z
            !! the number
        real(real64), intent(out) :: re
            !! its real part
        real(real64), intent(out) :: im
            !! its imaginary part
        re = real(z, kind=real64)
        im = aimag(z)
    end subroutine fx_parts_impl

    !> summary: A complex number from its parts, plus an offset if one is given
    !| author: A Developer
    !| An edge case: an optional `complex(real64)` input, which C may pass as null.
    pure subroutine fx_offset_impl(re, im, z, offset)
        real(real64), intent(in) :: re
            !! the real part
        real(real64), intent(in) :: im
            !! the imaginary part
        complex(real64), intent(out) :: z
            !! `re + im i`, plus `offset` if it is given
        complex(real64), intent(in), optional :: offset
            !! added to the number
        z = cmplx(re, im, kind=real64)
        if (present(offset)) z = z + offset
    end subroutine fx_offset_impl

    !> summary: Twice a complex number, in place
    !| author: A Developer
    !| An edge case: a `complex(real64)` scalar that is `intent(inout)`, returned by Python and R.
    pure subroutine fx_twice_impl(z)
        complex(real64), intent(inout) :: z
            !! the number, doubled on return
        z = 2.0_real64*z
    end subroutine fx_twice_impl

end module fx_functions_impl
