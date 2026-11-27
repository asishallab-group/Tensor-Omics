!> What the children of the `tox_data_tools` suite share: the two characters the text files are
!| built from, and the value an output holds before a call.
!|
!| Every file a case reads is written by the case itself, into a scratch directory, with its
!| content spelled out in the case: the expected values are read off that content.
module mod_test_tox_data_tools_fixtures
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    private
    public :: TAB, LF, UNWRITTEN

    !> The TSV field delimiter.
    character(len=1), parameter :: TAB = char(9)
    !> The line ending of every file in this suite.
    character(len=1), parameter :: LF = char(10)

    !> What an output holds before a call, so that an element the procedure never writes shows:
    !| an expected value equal to it says "left alone". No value in a file a case writes equals it.
    real(real64), parameter :: UNWRITTEN = -777.0_real64
end module mod_test_tox_data_tools_fixtures
