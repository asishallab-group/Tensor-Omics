#include <src/macros.h>

!> Vector geometry in n dimensions: lengths, angles, and element-wise arithmetic.
!|
!| One of the modules [[f42_utils(module)]] gathers; `use f42_utils` reaches all of them.
module f42_vector
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_DIVISION_BY_ZERO, set_ok, set_err
    use f42_math, only: clamp, is_close
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
end module f42_vector
