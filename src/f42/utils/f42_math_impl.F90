#include <src/macros.h>

!> Elementary mathematics: the constants, approximate equality, clamping, angles and logarithms.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_math_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_DIVISION_BY_ZERO, set_ok, set_err, validate_in_range_real
    use tox_errors, only: is_err
    use, intrinsic :: ieee_arithmetic, only: ieee_next_after, ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_is_finite
    M_IMPLICIT_NONE

    !> Generic clamp of a scalar into `[min_val, max_val]`, dispatches on the value's type.
    interface clamp
        module procedure clamp_real, clamp_int
    end interface clamp

    !> Tolerance-based approximate equality for `real(real64)`.
    !| Called as `is_close(a, b)` with the machine-epsilon tolerance, or as
    !| `is_close(a, b, eps)` where the caller has its own notion of "the same" -- a domain
    !| epsilon such as a smoothing floor, which is coarser than anything the arithmetic
    !| itself would justify. See [[f42_utils(module):is_close_within(function)]].
    interface is_close
        module procedure is_close_default
        module procedure is_close_within
    end interface is_close

    !> Tolerance-based approximate equality operator for `real(real64)`; see [[f42_utils(module):is_close_default(function)]].
    !| An operator takes exactly its two operands, so the custom-tolerance form is reachable
    !| only by calling [[f42_utils(module):is_close(interface)]] directly.
    interface operator(.isclose.)
        module procedure is_close_default
    end interface operator(.isclose.)

#define CM_EPS epsilon(1.0_real64)

    real(real64), parameter :: PI = 4.0_real64*atan(1.0_real64)
        !! The mathematical constant \( \pi \).
    real(real64), parameter :: EPS = CM_EPS
        !! Machine epsilon for `real64`, the base tolerance used by [[f42_utils(module):is_close(function)]].
    real(real64), parameter :: LOG_2 = log(2.0_real64)
        !! Natural logarithm of 2, used to compute base-2 logarithms in [[f42_utils(module):logx(subroutine)]].

contains

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
    pure logical function is_close_default(a, b)
        real(real64), intent(in) :: a
            !! First variable of comparison `a==b`
        real(real64), intent(in) :: b
            !! Second variable of comparison `a==b`

        is_close_default = is_close_within(a, b, EPS)
    end function is_close_default

    pure logical function is_close_within(a, b, eps)
        real(real64), intent(in) :: a
            !! First variable of comparison `a==b`
        real(real64), intent(in) :: b
            !! Second variable of comparison `a==b`
        real(real64), intent(in) :: eps
            !! Relative tolerance factor, scaled by the larger operand. Pass a domain epsilon
            !! to compare on that domain's terms rather than the arithmetic's.

        real(real64) :: rel_tolerance

        if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
            ! The absolute floor keeps the comparison meaningful where both operands are
            ! near zero and the relative term collapses with them.
            rel_tolerance = eps*max(abs(a), abs(b))
            is_close_within = abs(a - b) <= max(rel_tolerance, 1d-12)
        else
            is_close_within = a == b
        end if
    end function is_close_within

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

    !> AUTHOR_AARON_SCHROEDER
    !| Find the next power of two greater than or equal to n
    function next_power_of_two(n) result(power)
        integer(int32), intent(in) :: n
            !! input value
        integer(int32) :: power
            !! next greater value that is a power of two

        ! Guard against n<=0: for n==0, `n-1==-1` has all bits set, so `leadz(-1)==0` and the
        ! unguarded formula would compute `2**bit_size(n)`, an out-of-range shift. Negative n is
        ! similarly undefined. The smallest power of two is 1, so clamp to that.
        if (n <= 0) then
            power = 1
            return
        end if

        power = 2**(bit_size(n) - leadz(n - 1))
    end function next_power_of_two
end module f42_math_impl
