! Can a mask type with a pointer component own its words? Operators must return new masks, so
! someone has to free them. The type frees owned words in a final procedure. Compile with one of:
!   -DINTRINSIC_ASSIGNMENT   `=` copies the pointer (the default for a pointer component)
!   -DDEFINED_ASSIGNMENT     `=` is a defined assignment that copies the words
! Then `c = a .and. b` runs N_ROUNDS times, and the program reports how many buffers are still
! live. Correct is 3 (a, b and c). The expected outcome is in ../RESULTS.md.
module probe_owning_mask_mod
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none

    integer(int32) :: live_buffers = 0_int32

    type :: owning_mask
        integer(int32), pointer, contiguous :: words(:) => null()
        logical :: owns = .false.
    contains
        procedure, private :: and_masks
        generic :: operator(.and.) => and_masks
#ifdef DEFINED_ASSIGNMENT
        procedure, private :: copy_words
        generic :: assignment(=) => copy_words
#endif
        final :: release
    end type owning_mask

contains

    function new_mask(n_words) result(mask)
        integer(int32), intent(in) :: n_words
        type(owning_mask) :: mask

        allocate (mask%words(n_words), source=5_int32)
        mask%owns = .true.
        live_buffers = live_buffers + 1_int32
    end function new_mask

    function and_masks(left, right) result(conjunction)
        class(owning_mask), intent(in) :: left, right
        type(owning_mask) :: conjunction

        allocate (conjunction%words(size(left%words)))
        conjunction%words = iand(left%words, right%words)
        conjunction%owns = .true.
        live_buffers = live_buffers + 1_int32
    end function and_masks

#ifdef DEFINED_ASSIGNMENT
    subroutine copy_words(lhs, rhs)
        class(owning_mask), intent(inout) :: lhs
        class(owning_mask), intent(in) :: rhs

        if (.not. associated(lhs%words)) then
            allocate (lhs%words(size(rhs%words)))
            lhs%owns = .true.
            live_buffers = live_buffers + 1_int32
        end if
        lhs%words = rhs%words
    end subroutine copy_words
#endif

    subroutine release(self)
        type(owning_mask), intent(inout) :: self

        if (self%owns .and. associated(self%words)) then
            deallocate (self%words)
            live_buffers = live_buffers - 1_int32
        end if
        nullify (self%words)
        self%owns = .false.
    end subroutine release

end module probe_owning_mask_mod

program probe_owning_mask
    use probe_owning_mask_mod
    implicit none
    integer(int32), parameter :: N_ROUNDS = 100_int32, N_WORDS = 1000000_int32

    call rounds()
    print '(a,i0,a)', 'after leaving the scope: ', live_buffers, ' buffers live (correct: 0)'

contains

    subroutine rounds()
        type(owning_mask) :: a, b, c
        integer(int32) :: i_round

        a = new_mask(N_WORDS)
        b = new_mask(N_WORDS)
        do i_round = 1, N_ROUNDS
            c = a .and. b
        end do
        print '(a,i0,a,i0,a)', 'after ', N_ROUNDS, ' rounds of c = a .and. b: ', live_buffers, &
            ' buffers live (correct: 3)'
        print '(a,l1)', 'c holds a .and. b: ', all(c%words == 5_int32)
    end subroutine rounds

end program probe_owning_mask
