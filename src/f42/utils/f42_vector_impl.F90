#include <src/macros.h>

!> Vector geometry in n dimensions: lengths, angles, and element-wise arithmetic.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_vector_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: ERR_DIVISION_BY_ZERO, set_ok, set_err
    use f42_math_impl, only: clamp, is_close
    M_IMPLICIT_NONE

contains

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
        real(real64) :: theta, dot_product, norm1_sq, norm2_sq, norm_product

        call set_ok(ierr)

        dot_product = 0.0_real64
        norm1_sq = 0.0_real64
        norm2_sq = 0.0_real64
        do concurrent(i_dim=1:n_dims) shared(v1, v2) reduce(+:dot_product, norm1_sq, norm2_sq)
            dot_product = dot_product + v1(i_dim)*v2(i_dim)
            norm1_sq = norm1_sq + v1(i_dim)**2
            norm2_sq = norm2_sq + v2(i_dim)**2
        end do

        norm_product = sqrt(norm1_sq)*sqrt(norm2_sq)
        if (is_close(norm_product, 0.0_real64)) then
            angle = 0.0_real64
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        theta = dot_product/norm_product
        theta = clamp(theta, -1.0_real64, 1.0_real64)
        angle = acos(theta)
    end subroutine angle_between

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

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Length of a vector as `scaled_norm / factor`, where `scaled_norm` is the Euclidean length
    !| of `vector * factor` and `factor` a power of two that brings the largest component into
    !| \([2^{-500}, 2^{500}]\) -- 1 for most data. The plain sum of squares overflows once a
    !| component passes about `1e154` and underflows to zero once all lie below about `1e-162`;
    !| scaled, it does neither for any finite vector, and `scaled_norm` stays within
    !| \([2^{-500}, 2^{500} \sqrt{n}]\). `vector(i) * factor / scaled_norm` is then the unit vector's
    !| component.
    !|
    !| Multiplying by a power of two is exact, so scaling loses nothing: a vector and any power-of-two
    !| multiple of it give the same unit vector to the last bit, and lengths that differ by exactly
    !| that power, as long as no square in either sum falls below `tiny`. That is why the factor is a
    !| constant rather than the reciprocal of the largest component, which would round, and which is
    !| subnormal for components near `huge` -- zero in a floating-point environment that flushes
    !| subnormals.
    !|
    !| A vector whose largest component is subnormal counts as the zero vector, so that the answer
    !| does not depend on that environment: a build with a fast floating-point model, or a host
    !| program that sets flush-to-zero or denormals-are-zero, already reads such components as zero.
    pure subroutine scaled_length(n_dims, vector, factor, scaled_norm)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `vector`
        real(real64), intent(in) :: vector(n_dims)
            !! The vector
        real(real64), intent(out) :: factor
            !! Power of two the vector is multiplied by before its length is taken
        real(real64), intent(out) :: scaled_norm
            !! Euclidean length of `vector * factor`; zero exactly for the zero vector and for one
            !! whose components are all subnormal

        real(real64), parameter :: LARGEST_UNSCALED = 2.0_real64**500, SMALLEST_UNSCALED = 2.0_real64**(-500)
        real(real64), parameter :: SCALE_DOWN = 2.0_real64**(-600), SCALE_UP = 2.0_real64**600
        real(real64) :: largest
        integer(int32) :: i_dim

        largest = 0.0_real64
        do i_dim = 1, n_dims
            largest = max(largest, abs(vector(i_dim)))
        end do

        factor = 1.0_real64
        scaled_norm = 0.0_real64
        if (largest < tiny(1.0_real64)) return
        if (largest > LARGEST_UNSCALED) factor = SCALE_DOWN
        if (largest < SMALLEST_UNSCALED) factor = SCALE_UP

        do concurrent (i_dim = 1:n_dims) shared(vector, factor) reduce(+:scaled_norm)
            scaled_norm = scaled_norm + (vector(i_dim)*factor)**2
        end do
        scaled_norm = sqrt(scaled_norm)
    end subroutine scaled_length

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether a vector has a direction: it is not the zero vector, nor one whose components are
    !| all subnormal -- the rule of [[f42_vector_impl(module):scaled_length(subroutine)]].
    pure logical(c_bool) function has_direction(n_dims, vector)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `vector`
        real(real64), intent(in) :: vector(n_dims)
            !! The vector

        real(real64) :: factor, scaled_norm

        call scaled_length(n_dims, vector, factor, scaled_norm)
        has_direction = scaled_norm > 0.0_real64
    end function has_direction

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The angle in radians between the directions of two vectors, in \([0, \pi]\), or `-1.0`
    !| where either has no direction (see [[f42_vector_impl(module):has_direction(function)]]).
    !|
    !| For the unit vectors \(u\) and \(v\) it is \(2\,\mathrm{atan2}(|u - v|, |u + v|)\): unlike
    !| the acos of their dot product, which [[f42_vector_impl(module):angle_between(subroutine)]]
    !| takes, this stays accurate near 0 and \(\pi\), where acos turns one rounding error of the
    !| cosine into an angle of about `1e-8`. Both vectors are scaled as
    !| [[f42_vector_impl(module):scaled_length(subroutine)]] does, so no finite magnitude overflows
    !| or underflows.
    pure real(real64) function angle_to_direction(n_dims, vector, direction) result(angle)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in both vectors
        real(real64), intent(in) :: vector(n_dims)
            !! The vector, of any length
        real(real64), intent(in) :: direction(n_dims)
            !! The direction, of any length

        real(real64) :: vector_factor, vector_scaled_norm, direction_factor, direction_scaled_norm
        real(real64) :: vector_unit, direction_unit, difference_squared, sum_squared
        integer(int32) :: i_dim

        angle = -1.0_real64
        call scaled_length(n_dims, vector, vector_factor, vector_scaled_norm)
        call scaled_length(n_dims, direction, direction_factor, direction_scaled_norm)
        if (vector_scaled_norm == 0.0_real64 .or. direction_scaled_norm == 0.0_real64) return

        difference_squared = 0.0_real64
        sum_squared = 0.0_real64
        do i_dim = 1, n_dims
            vector_unit = (vector(i_dim)*vector_factor)/vector_scaled_norm
            direction_unit = (direction(i_dim)*direction_factor)/direction_scaled_norm
            difference_squared = difference_squared + (vector_unit - direction_unit)**2
            sum_squared = sum_squared + (vector_unit + direction_unit)**2
        end do
        angle = 2.0_real64*atan2(sqrt(difference_squared), sqrt(sum_squared))
    end function angle_to_direction
end module f42_vector_impl
