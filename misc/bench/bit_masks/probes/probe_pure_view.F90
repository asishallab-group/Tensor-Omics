! Can a pure implementation hold a mask type as a view of its arguments? TensorOmics
! implementations are pure, and most masks reach them as intent(in) arguments. Compile with
! exactly one of:
!   -DVIEW_OF_INOUT   control: a view of words the procedure may modify
!   -DVIEW_OF_INPUT   a view of intent(in) words: every read-only mask argument, and every lazy
!                     expression node, which has to point at its operands
! The expected outcome is in ../RESULTS.md.
module probe_pure_view_mod
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none

    type :: view_mask
        integer(int32), pointer, contiguous :: words(:) => null()
    end type view_mask

contains

#ifdef VIEW_OF_INOUT
    pure subroutine count_bits(n_words, words, n_set)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(inout), target :: words(n_words)
        integer(int32), intent(out) :: n_set
        type(view_mask) :: view

        view%words => words
        n_set = sum(popcnt(view%words))
    end subroutine count_bits
#endif

#ifdef VIEW_OF_INPUT
    pure subroutine count_bits(n_words, words, n_set)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in), target :: words(n_words)
        integer(int32), intent(out) :: n_set
        type(view_mask) :: view

        view%words => words
        n_set = sum(popcnt(view%words))
    end subroutine count_bits
#endif

end module probe_pure_view_mod

program probe_pure_view
    use probe_pure_view_mod
    implicit none
    integer(int32), target :: words(2) = [7_int32, 1_int32]
    integer(int32) :: n_set

    call count_bits(2_int32, words, n_set)
    print '(a,i0)', 'compiled and ran, count = ', n_set
end program probe_pure_view
