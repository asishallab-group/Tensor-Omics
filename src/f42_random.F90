!> This module implements a custom random number generator to allow compiler independent reproducibility and usage in pure routines.
!|
!| It uses the squares CBRNG, inspired by reference: https://squaresrng.wixsite.com/rand   (version 8 - Feb 4, 2024)
module f42_random
    use, intrinsic :: iso_fortran_env, only: int32, real64, int64
    implicit none

    ! TODO: if Fortran ever supports unsigned integers, change int64 to uint64.
    ! According to this proposal from 2024, unsigned integers might become part of the next standard (today: 23rd, Feb. 2026): https://j3-fortran.org/doc/year/24/24-116.txt
contains
    !> Produces a `int32` of random bits.
    !| Less efficient than [[f42_random(module):random_int64(function)]].
    pure integer(int32) function random_int32(key, count) result(res)
        integer(int32), intent(in) :: key
        integer(int32), intent(in) :: count

        integer(int64) :: x, y, z

        call init(key, count, x, y, z)
        call round(x, y)
        call round(x, z)
        call round(x, y)
        call square(x, z)

        res = int(ishft(x, -32), int32)
    end function random_int32

    pure integer(int64) function random_int64(key, count) result(res)
        integer(int32), intent(in) :: key
        integer(int32), intent(in) :: count

        integer(int64) :: x, y, z

        call init(key, count, x, y, z)
        call round(x, y)
        call round(x, z)
        call round(x, y)
        res = x
        call square(res, z)
        x = res
        call swap_high_low(x)
        call square(x, y)

        res = ieor(res, ishft(x, -32))
    end function random_int64

    pure real(real64) function random_uniform(key, count) result(res)
        integer(int32), intent(in) :: key
        integer(int32), intent(in) :: count

        res = abs(real(random_int64(key, count), real64) / real(huge(1_int64), real64))
    end function random_uniform

    pure real(real64) function rand_range(min, max, key, count) result(res)
        integer(int32), intent(in) :: key
        integer(int32), intent(in) :: count
        real(real64), intent(in) :: min
        real(real64), intent(in) :: max

        res = min + random_uniform(key, count) * (max - min)
    end function rand_range

    pure subroutine swap_high_low(x)
        integer(int64), intent(inout) :: x

        integer(int64) :: high

        high = ishft(x, -32) ! x >> 32
        x = ior(high, ishft(x, 32)) ! (x >> 32) | (x << 32)
    end subroutine swap_high_low

    pure integer(int64) function square64(x) result(res)
        integer(int64), intent(in) :: x
        integer(int64) :: low, high
        integer(int64) :: low_low, low_high, high_low

        low = iand(x, z'FFFFFFFF')
        high = ishft(x, -32)

        res = low * low + ishft(2 * high * low, 32)
    end function square64

    pure integer(int64) function add64(a, b) result(res)
        integer(int64), intent(in) :: a
        integer(int64), intent(in) :: b
        integer(int64) :: a_low, a_high, b_low, b_high, carry
        integer(int64) :: low, high_low, low_high, high

        a_low = iand(a, z'FFFFFFFF')
        a_high = ishft(a, -32)
        b_low = iand(b, z'FFFFFFFF')
        b_high = ishft(b, -32)

        low = a_low + b_low
        low = iand(low, z'FFFFFFFF')    ! keep low 32 bits

        carry = ishft(low, -32)        ! carry = lo >> 32
      
        high = a_high + b_high + carry
        high = ishft(high, 32)    ! move low bit to high

        ! Combine the pieces, all modulo 2^64
        res = ior(high, low)
    end function add64

    !> perform `x * x + o`
    pure subroutine square(x, o)
        integer(int64), intent(inout) :: x
        integer(int64), intent(in) :: o

        x = add64(square64(x), o)
    end subroutine square

    pure subroutine round(x, o)
        integer(int64), intent(inout) :: x
        integer(int64), intent(in) :: o
    
        call square(x, o)
        call swap_high_low(x)
    end subroutine round

    pure subroutine init(key, counter, x, y, z)
        integer(int32), intent(in) :: key
        integer(int32), intent(in) :: counter
        integer(int64), intent(out) :: x
        integer(int64), intent(out) :: y
        integer(int64), intent(out) :: z
        
        x = int(key, int64) * int(counter, int64)
        y = x
        z = y + int(key, int64)
    end subroutine init
end module f42_random