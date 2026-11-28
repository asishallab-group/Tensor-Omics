!> What the children of the `tox_clustering` suite share: the values an output holds before a
!| call, and the three non-finite values every procedure has to refuse.
module mod_test_tox_clustering_fixtures
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    implicit none
    private
    public :: UNWRITTEN_INT, UNWRITTEN_REAL, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value

    !> What an integer output (a label, a count, a merge, a cluster size) holds before a call, so
    !| that an element the procedure never writes shows. No expected value in this suite equals it.
    integer(int32), parameter :: UNWRITTEN_INT = -777_int32
    !> What a real output (a height) holds before a call; no expected value equals it either.
    real(real64), parameter :: UNWRITTEN_REAL = -777.0_real64

    !> How many non-finite values `non_finite_value` hands out.
    integer(int32), parameter :: N_NON_FINITE = 3
    !> The names of the non-finite values, for the assertion messages.
    character(*), parameter :: NON_FINITE_NAMES(N_NON_FINITE) = ["NaN ", "+Inf", "-Inf"]

contains

    !> The `i_value`-th non-finite value: NaN, +Inf, -Inf.
    function non_finite_value(i_value) result(value)
        integer(int32), intent(in) :: i_value
        real(real64) :: value

        select case (i_value)
        case (1)
            value = ieee_value(1.0_real64, ieee_quiet_nan)
        case (2)
            value = ieee_value(1.0_real64, ieee_positive_inf)
        case default
            value = ieee_value(1.0_real64, ieee_negative_inf)
        end select
    end function non_finite_value
end module mod_test_tox_clustering_fixtures
