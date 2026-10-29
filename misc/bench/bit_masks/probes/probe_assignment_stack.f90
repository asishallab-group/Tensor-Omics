! `mask = flags` as a defined assignment, at 10^7 bits, under the default 8 MB stack. The
! standard passes the right-hand side of a defined assignment as if it were in parentheses, so
! a compiler may copy it first -- onto the stack. The core's `call core_pack(...)` on the same
! data is the control. Link with mask_candidates.o. The expected outcome is in ../RESULTS.md.
program probe_assignment_stack
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use mask_candidates
    implicit none
    integer(int32), parameter :: N_BITS = 10000000_int32
    logical(c_bool), allocatable :: flags(:)
    integer(int32), allocatable :: words(:)
    type(owned_mask) :: mask

    allocate (flags(N_BITS), source=.true._c_bool)
    allocate (words(words_for_bits(N_BITS)))
    call core_pack(N_BITS, words_for_bits(N_BITS), flags, words)
    print '(a,i0)', 'core_pack:  count = ', core_count(words_for_bits(N_BITS), words)
    flush (6)
    mask = flags
    print '(a,i0)', 'mask = flags: count = ', owned_count(mask)
end program probe_assignment_stack
