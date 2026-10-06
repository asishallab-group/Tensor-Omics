#include <src/macros.h>

!> Elementary mathematics: the constants, approximate equality, clamping, angles and logarithms.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_math_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
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
    !| itself would justify. See [[f42_math_impl(module):is_close_within(function)]].
    interface is_close
        module procedure is_close_default
        module procedure is_close_within
    end interface is_close

    !> Tolerance-based approximate equality operator for `real(real64)`; see [[f42_math_impl(module):is_close_default(function)]].
    !| An operator takes exactly its two operands, so the custom-tolerance form is reachable
    !| only by calling [[f42_math_impl(module):is_close(interface)]] directly.
    interface operator(.isclose.)
        module procedure is_close_default
    end interface operator(.isclose.)

#define CM_EPS epsilon(1.0_real64)

    real(real64), parameter :: PI = 4.0_real64*atan(1.0_real64)
        !! The mathematical constant \( \pi \).
    real(real64), parameter :: EPS = CM_EPS
        !! Machine epsilon for `real64`, the base tolerance used by [[f42_math_impl(module):is_close(interface)]].
    real(real64), parameter :: LOG_2 = log(2.0_real64)
        !! Natural logarithm of 2, used to compute base-2 logarithms in [[f42_math_impl(module):logx(subroutine)]].

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
        logical(c_bool), intent(in), optional :: do_bessel_correction
            !! Tells whether to apply the bessel's correction or not, default: `.false.`
            !!
            !! |    Case     |                                                Formula                                                      |
            !! |-------------|-------------------------------------------------------------------------------------------------------------|
            !! |  `.true.`   | \(\frac{1}{\texttt{size}(vec) - 1} \cdot \sum_{i=1}^{\texttt{size}(vec)} (vec(i) - \texttt{mean}(i))^{2}\)  |
            !! |  `.false.`  |  \(\frac{1}{\texttt{size}(vec)} \cdot \sum_{i=1}^{\texttt{size}(vec)} vec(i)^{2} - \texttt{mean}(i)^{2}\)   |

        logical(c_bool) :: bessel
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

        call logx_helper(val, base, exponent)
    end subroutine logx

    !> AUTHOR_FRANZ_ERIC_SILL
    !| (no input validation) Compute logarithm for any base: `log_base(val) = log(val)/log(base)`.
    !| Ensure `val > 0`, `base > 0` and `base /= 1`, yields a NaN/Inf result otherwise.
    !|
    !| For a caller that has already established those three facts about a whole array -- so the
    !| per-element check would only re-derive what is known -- and therefore wants the transform
    !| in a `do concurrent` with no shared `ierr` to write.
    pure subroutine logx_helper(val, base, exponent)
        real(real64), intent(in) :: val
            !! Value (`x` in \( b^y = x \)), must be `> 0`
        real(real64), intent(in) :: base
            !! Base (`b` in \( b^y = x \)), must be `> 0` and `/= 1`
        real(real64), intent(out) :: exponent
            !! Exponent (`y` in \( b^y = x \))

        if (base == 2.0_real64) then
            exponent = log(val)/LOG_2
        else
            exponent = log(val)/log(base)
        end if
    end subroutine logx_helper

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

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns the given angle wrapped into \( (-\pi, \pi] \), the same direction on the circle:
    !| \( \pi \) stays \( \pi \), \( -\pi \Rightarrow \pi \), \( \frac{3\cdot \pi}{2} \Rightarrow -\frac{\pi}{2} \).
    !| An angle already in \( (-\pi, \pi] \) is returned unchanged, to the last bit.
    !| The angle must be finite; a NaN or infinite one gives NaN.
    pure elemental real(real64) function wrap_angle(angle) result(wrapped)
        real(real64), intent(in) :: angle
            !! angle to be wrapped, in radians

        if (angle > -PI .and. angle <= PI) then
            wrapped = angle
            return
        end if

        ! `modulo` lies in [0, 2*PI), so `PI - modulo` lies in (-PI, PI]. Rounding can still land
        ! a value a hair beyond PI on exactly -PI, which is the same direction as PI.
        wrapped = PI - modulo(PI - angle, 2*PI)
        if (wrapped <= -PI) wrapped = PI
    end function wrap_angle

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns the absolute distance around the circle between two angles, in \( [0, \pi] \): an
    !| angle just below \( \pi \) and one just above \( -\pi \) are close. Their difference is
    !| wrapped by [[f42_math_impl(module):wrap_angle(function)]], so it must be finite.
    pure elemental real(real64) function angular_distance(angle, reference_angle) result(distance)
        real(real64), intent(in) :: angle
            !! An angle in radians
        real(real64), intent(in) :: reference_angle
            !! The angle it is measured from, in radians

        distance = abs(wrap_angle(angle - reference_angle))
    end function angular_distance

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns \( 1 - \cos x \), computed as \( 2\sin^2(x/2) \), which keeps its full relative
    !| precision for small \( x \), where \( 1 - \cos x \) cancels.
    pure elemental real(real64) function one_minus_cosine(angle)
        real(real64), intent(in) :: angle
            !! Angle in radians

        one_minus_cosine = 2.0_real64*sin(0.5_real64*angle)**2
    end function one_minus_cosine

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Returns \( \ln(1 - x) \), accurate also where \( x \) is small and \( 1 - x \) would round
    !| away its digits. Below \( x = 10^{-3} \) it sums the series \( -(x + x^2/2 + \dots + x^6/6) \),
    !| whose truncation error for \( |x| < 10^{-3} \) is under \( x^6/7 < 1.5 \times 10^{-19} \)
    !| relative to its value: the result is within about an ulp. From \( 10^{-3} \) on it takes
    !| `log(1 - x)`, whose only extra error is the rounding of \( 1 - x \), at most \( 2^{-54} \):
    !| under \( 5.6 \times 10^{-14} \) relative at \( x = 10^{-3} \), falling as \( x \) grows, and
    !| none from \( x = 1/2 \) on, where \( 1 - x \) is exact.
    !|
    !| (no input validation) Ensure \( -10^{-3} < x < 1 \) for the accuracy above. For
    !| \( -1 < x \le -10^{-3} \) the series still converges, but cut after six terms it is off by
    !| about \( |x|^6/7 \) relative (\( 2 \times 10^{-3} \) at \( x = -1/2 \)); `x >= 1` yields an
    !| infinite or NaN result.
    pure elemental real(real64) function log_one_minus(x)
        real(real64), intent(in) :: x
            !! Argument, in `(-1e-3, 1)`

        if (x < 1.0e-3_real64) then
            log_one_minus = -x*(1.0_real64 + x*(0.5_real64 + x*(1.0_real64/3.0_real64 &
                            + x*(0.25_real64 + x*(0.2_real64 + x/6.0_real64)))))
        else
            log_one_minus = log(1.0_real64 - x)
        end if
    end function log_one_minus

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
