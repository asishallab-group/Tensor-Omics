! The consumer side: what an implementation using f42_bit_masks would write, one kernel per
! way of doing each operation. A separate unit from mask_candidates, so every call into it is a
! real call, and from the driver, so nothing here can be elided (see ../README.md).
!
! Each operation comes in the forms the two designs offer, with logical(c_bool) -- what the
! code base uses today -- as the reference:
!   flags   logical(c_bool) arrays
!   core    the chosen design: plain words, in-place or fused procedures
!   eager   owned_mask's operator or assignment, which allocates its result
!   view    view_mask, the words behind a pointer component
!   lazy    mask_expression, evaluated per bit or per word on demand
module bench_kernels
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use, intrinsic :: iso_c_binding, only: c_bool
    use mask_candidates
    implicit none
    private

    public :: and_count_flags, and_count_core, and_count_eager, and_count_lazy_bits, and_count_lazy_words
    public :: and_flags, and_core, and_eager, and_view
    public :: copy_flags, copy_core, copy_owned
    public :: pack_core, pack_assignment
    public :: merge_flags, merge_core

contains

    ! ------------------------------------------------ count(a .and. b): Shatter's merge test

    integer(int64) function and_count_flags(n_bits, left, right)
        integer(int32), intent(in) :: n_bits
        logical(c_bool), intent(in) :: left(n_bits), right(n_bits)

        and_count_flags = count(left .and. right, kind=int64)
    end function and_count_flags

    integer(int64) function and_count_core(n_words, left, right)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in) :: left(n_words), right(n_words)

        and_count_core = core_and_count(n_words, left, right)
    end function and_count_core

    !> `count(a .and. b)` with the type's operator: the conjunction is materialised first.
    integer(int64) function and_count_eager(left, right)
        type(owned_mask), intent(in) :: left, right
        type(owned_mask) :: conjunction

        conjunction = left .and. right
        and_count_eager = owned_count(conjunction)
    end function and_count_eager

    integer(int64) function and_count_lazy_bits(left, right)
        type(owned_mask), target, intent(in) :: left, right
        class(mask_expression), allocatable :: conjunction
        integer(int32) :: i_bit

        allocate (conjunction, source=leaf_of(left) .and. leaf_of(right))
        and_count_lazy_bits = 0_int64
        do i_bit = 1, left%n_bits
            if (conjunction%bit(i_bit)) and_count_lazy_bits = and_count_lazy_bits + 1_int64
        end do
    end function and_count_lazy_bits

    integer(int64) function and_count_lazy_words(left, right)
        type(owned_mask), target, intent(in) :: left, right
        class(mask_expression), allocatable :: conjunction
        integer(int32) :: i_word

        allocate (conjunction, source=leaf_of(left) .and. leaf_of(right))
        and_count_lazy_words = 0_int64
        do i_word = 1, size(left%words, kind=int32)
            and_count_lazy_words = and_count_lazy_words + popcnt(conjunction%word(i_word))
        end do
    end function and_count_lazy_words

    ! ------------------------------------------------ c = a .and. b, into an existing mask

    subroutine and_flags(n_bits, left, right, conjunction)
        integer(int32), intent(in) :: n_bits
        logical(c_bool), intent(in) :: left(n_bits), right(n_bits)
        logical(c_bool), intent(inout) :: conjunction(n_bits)

        conjunction = left .and. right
    end subroutine and_flags

    subroutine and_core(n_words, left, right, conjunction)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in) :: left(n_words), right(n_words)
        integer(int32), intent(inout) :: conjunction(n_words)

        call core_and_into(n_words, left, right, conjunction)
    end subroutine and_core

    subroutine and_eager(left, right, conjunction)
        type(owned_mask), intent(in) :: left, right
        type(owned_mask), intent(inout) :: conjunction

        conjunction = left .and. right
    end subroutine and_eager

    subroutine and_view(left, right, conjunction)
        type(view_mask), intent(in) :: left, right
        type(view_mask), intent(inout) :: conjunction

        call view_and_into(left, right, conjunction)
    end subroutine and_view

    ! ------------------------------------------------ the per-iteration backup copy

    subroutine copy_flags(n_bits, source, destination)
        integer(int32), intent(in) :: n_bits
        logical(c_bool), intent(in) :: source(n_bits)
        logical(c_bool), intent(inout) :: destination(n_bits)

        destination = source
    end subroutine copy_flags

    subroutine copy_core(n_words, source, destination)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in) :: source(n_words)
        integer(int32), intent(inout) :: destination(n_words)

        call core_copy_into(n_words, source, destination)
    end subroutine copy_core

    !> Intrinsic assignment of a type with an allocatable component: the standard deallocates
    !| and reallocates the component whenever the shapes differ, and compilers do it anyway.
    subroutine copy_owned(source, destination)
        type(owned_mask), intent(in) :: source
        type(owned_mask), intent(inout) :: destination

        destination = source
    end subroutine copy_owned

    ! ------------------------------------------------ logical(c_bool) -> mask, at the API edge

    subroutine pack_core(n_bits, n_words, flags, words)
        integer(int32), intent(in) :: n_bits, n_words
        logical(c_bool), intent(in) :: flags(n_bits)
        integer(int32), intent(inout) :: words(n_words)

        call core_pack(n_bits, n_words, flags, words)
    end subroutine pack_core

    !> `mask = flags`: the right-hand side of a defined assignment is passed as if it were in
    !| parentheses, which lets a compiler copy it first.
    subroutine pack_assignment(flags, mask)
        logical(c_bool), intent(in) :: flags(:)
        type(owned_mask), intent(inout) :: mask

        mask = flags
    end subroutine pack_assignment

    ! ------------------------------------------------ Shatter's merge: every pair of columns

    integer(int64) function merge_flags(n_bits, n_masks, masks)
        integer(int32), intent(in) :: n_bits, n_masks
        logical(c_bool), intent(in) :: masks(n_bits, n_masks)
        integer(int32) :: i_mask, j_mask

        merge_flags = 0_int64
        do i_mask = 1, n_masks - 1_int32
            do j_mask = i_mask + 1_int32, n_masks
                merge_flags = merge_flags + count(masks(:, i_mask) .and. masks(:, j_mask), kind=int64)
            end do
        end do
    end function merge_flags

    integer(int64) function merge_core(n_words, n_masks, masks)
        integer(int32), intent(in) :: n_words, n_masks
        integer(int32), intent(in) :: masks(n_words, n_masks)
        integer(int32) :: i_mask, j_mask

        merge_core = 0_int64
        do i_mask = 1, n_masks - 1_int32
            do j_mask = i_mask + 1_int32, n_masks
                merge_core = merge_core + core_and_count(n_words, masks(:, i_mask), masks(:, j_mask))
            end do
        end do
    end function merge_core

end module bench_kernels
