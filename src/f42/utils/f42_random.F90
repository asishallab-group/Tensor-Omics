#include <src/macros.h>

!> Randomness: seeding, uniform draws in a range, and the in-place Fisher-Yates shuffle.
!|
!| One of the modules [[f42_utils(module)]] gathers; `use f42_utils` reaches all of them.
module f42_random
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use f42_sort, only: swap_int, swap_real
    M_IMPLICIT_NONE

    !> Generic in-place Fisher-Yates shuffle, dispatches on the vector's element type.
    interface shuffle_vector
        module procedure shuffle_vector_real, shuffle_vector_int
    end interface shuffle_vector

contains

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
end module f42_random
