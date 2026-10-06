!> The three candidate spellings of `subnormal_as_zero`, in their own compilation unit as in the
!| real library, where the helper lives in `f42_math_impl` and its callers in `f42_vector_impl`.
module bench_helpers
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_class, ieee_class_type, ieee_positive_denormal, &
                                             ieee_negative_denormal, operator(==)
    implicit none
    private
    public :: saz_exponent, saz_class, saz_bits

    integer(int64), parameter :: MAGNITUDE_MASK = huge(0_int64)           ! every bit but the sign
    integer(int64), parameter :: SMALLEST_NORMAL_BITS = transfer(tiny(1.0_real64), 0_int64)
contains
    pure real(real64) function saz_exponent(val)
        real(real64), intent(in) :: val
        saz_exponent = merge(0.0_real64, val, exponent(val) < minexponent(val))
    end function saz_exponent

    pure real(real64) function saz_class(val)
        real(real64), intent(in) :: val
        saz_class = merge(0.0_real64, val, ieee_class(val) == ieee_positive_denormal .or. &
                                           ieee_class(val) == ieee_negative_denormal)
    end function saz_class

    !> Positive IEEE doubles order like their bit patterns, so a non-zero magnitude below the
    !| bits of `tiny` is exactly a subnormal. Being integer arithmetic on the representation, it
    !| still classifies correctly with denormals-are-zero on; zero of either sign passes through.
    pure real(real64) function saz_bits(val)
        real(real64), intent(in) :: val
        integer(int64) :: magnitude
        magnitude = iand(transfer(val, 0_int64), MAGNITUDE_MASK)
        saz_bits = merge(0.0_real64, val, magnitude > 0_int64 .and. magnitude < SMALLEST_NORMAL_BITS)
    end function saz_bits
end module bench_helpers
