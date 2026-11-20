#include <src/macros.h>

!> Randomness: seeding, uniform draws in a range, and the in-place Fisher-Yates shuffle.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_random_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use f42_sort_impl, only: swap_int, swap_real
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
            !! Every seed word is derived from the `int32` seed and only converted at `put=`, so any
            !! `int32` value is accepted.
            !! @endnote

        integer(int32) :: actual_seed

        ! IMPORTANT: these locals need to be default kind integer
        integer :: size, i

        M_DEFAULT_VAL(seed, actual_seed, 42_int32)

        ! determine needed array size to seed random numbers
        call random_seed(size=size)

        ! One distinct word per slot, all derived from the one seed. Repeating the seed in every
        ! word, as this did before, leaves some generators in a degenerate state: nvfortran's
        ! (34 words) then draws 0.00016 again and again, so a permutation test barely permutes.
        ! Still a pure function of `seed`, so a seeded run stays reproducible.
        call random_seed(put=[(int(seed_word(actual_seed, int(i, int32))), i=1, size)])
    end subroutine init_random

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The `word`-th seed word derived from `seed`: `word` steps of the MINSTD generator
    !| (multiplier 48271, modulus 2^31 - 1), started from a nonzero value. The products stay far
    !| below the int64 range, and every word lands in [1, 2^31 - 2], a valid default integer.
    pure integer(int32) function seed_word(seed, word) result(value)
        integer(int32), intent(in) :: seed
            !! the seed every word is derived from
        integer(int32), intent(in) :: word
            !! which word, 1-based
        integer(int64), parameter :: MODULUS = 2147483647_int64
        integer(int64), parameter :: MULTIPLIER = 48271_int64
        integer(int64) :: state
        integer(int32) :: step

        state = modulo(int(seed, int64), MODULUS - 1) + 1
        do step = 1, word
            state = modulo(state*MULTIPLIER, MODULUS)
        end do
        value = int(state, int32)
    end function seed_word

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
end module f42_random_impl
