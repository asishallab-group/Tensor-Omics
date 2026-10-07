#include <src/macros.h>

!> Vector geometry in n dimensions: lengths, angles, and element-wise arithmetic.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_vector_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_DIVISION_BY_ZERO, set_ok, set_err
    use f42_math_impl, only: clamp, is_close, scaling_exponent
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
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        theta = dot_product/norm_product
        theta = clamp(theta, -1.0_real64, 1.0_real64)
        angle = acos(theta)
    end subroutine angle_between

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Calculates the euclidean norm of a vector, as
    !| [[f42_vector_impl(module):scaled_length(subroutine)]] scaled back: squaring an entry neither
    !| overflows (past about 1e154) nor underflows to zero (below about 1e-162), so a vector of any
    !| finite magnitude gets its norm. The result is Inf only where the true norm exceeds `huge`, as
    !| for `[huge, huge]`; the zero vector has norm 0.
    pure real(real64) function norm(vector)
        real(real64), dimension(:), contiguous, intent(in) :: vector
            !! Input vector the norm will be calculated for

        integer(int32) :: exponent
        real(real64) :: scaled_norm

        call scaled_length(size(vector, kind=int32), vector, exponent, scaled_norm)
        norm = scale(scaled_norm, -exponent)
    end function norm

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Length of a vector in scaled coordinates: `scaled_norm` is the euclidean length of
    !| `scale(vector, exponent)`, and `scale(scaled_norm, -exponent)` the vector's own length. The
    !| exponent comes from [[f42_math_impl(module):scaling_exponent(function)]], so it is 0 for
    !| ordinary data, and the sum of squares overflows or underflows for no finite vector.
    !|
    !| `scale(vector(i), exponent)/scaled_norm` is the unit vector's component. Computed that way,
    !| it never passes through the vector's length itself, which may overflow or be subnormal.
    !| `scaled_norm` is exactly 0 for the zero vector, and positive for every other one -- but see
    !| the next paragraph for a vector made only of subnormal numbers.
    !|
    !| An input made only of subnormal numbers can give different results in builds that flush
    !| subnormals to zero (a fast floating-point model, or a host program that sets flush-to-zero or
    !| denormals-are-zero): there it reads as the zero vector, with a `scaled_norm` of 0.
    pure subroutine scaled_length(n_dims, vector, exponent, scaled_norm)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `vector`
        real(real64), dimension(n_dims), intent(in) :: vector
            !! The vector, of any finite magnitude
        integer(int32), intent(out) :: exponent
            !! Power of two the vector is scaled by before its length is taken; meaningless for
            !! the zero vector
        real(real64), intent(out) :: scaled_norm
            !! Euclidean length of `scale(vector, exponent)`; exactly 0 for the zero vector

        integer(int32) :: i_dim
        real(real64) :: scaled_squares_sum

        exponent = scaling_exponent(vector)

        scaled_squares_sum = 0.0_real64
        do concurrent(i_dim=1:n_dims) shared(vector, exponent) reduce(+:scaled_squares_sum)
            scaled_squares_sum = scaled_squares_sum + scale(vector(i_dim), exponent)**2
        end do
        scaled_norm = sqrt(scaled_squares_sum)
    end subroutine scaled_length

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
end module f42_vector_impl
