! The two designs weighed for f42_bit_masks (issue #167), in their own compilation unit, as the
! real module would be: a consumer calling into it pays a real call (see ../README.md).
!
! The core, which was chosen: pure procedures on plain int32 words. A mask of n_bits bits is
! words(words_for_bits(n_bits)); bit i is bit mod(i - 1, 32) of word (i - 1)/32 + 1, and the
! bits past n_bits stay zero. A column of an n_bits x n_masks mask is simply words(:, i_mask).
!
! The type, which was not: a derived type `bit_mask` that behaves like a logical array. It is
! here in the three forms it could take:
!   - owned_mask: the words in an allocatable component, with an eager `.and.` returning a new
!     mask, and a defined assignment from a logical(c_bool) array;
!   - view_mask: the words behind a pointer component, as first sketched, so that a mask can
!     be a view of words a caller owns;
!   - mask_expression: lazy evaluation, where `a .and. b` builds a node that evaluates a bit
!     (or a word) only when it is asked for, instead of materialising the result.
module mask_candidates
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_bool
    implicit none
    private

    integer(int32), parameter, public :: BITS_PER_WORD = 32_int32

    public :: words_for_bits, core_pack, core_and_into, core_copy_into, core_count, core_and_count
    public :: owned_mask, owned_count, view_mask, view_and_into
    public :: mask_expression, mask_leaf, mask_and_node, leaf_of
    public :: operator(.and.), assignment(=)

    type :: owned_mask
        integer(int32), allocatable :: words(:)
        integer(int32) :: n_bits = 0
    end type owned_mask

    type :: view_mask
        integer(int32), pointer, contiguous :: words(:) => null()
        integer(int32) :: n_bits = 0
    end type view_mask

    !> A lazily evaluated mask: a node answers one bit, or one whole word, on demand.
    type, abstract :: mask_expression
    contains
        procedure(bit_of_expression), deferred :: bit
        procedure(word_of_expression), deferred :: word
    end type mask_expression

    abstract interface
        pure logical function bit_of_expression(self, i_bit)
            import :: mask_expression, int32
            class(mask_expression), intent(in) :: self
            integer(int32), intent(in) :: i_bit
        end function bit_of_expression

        pure integer(int32) function word_of_expression(self, i_word)
            import :: mask_expression, int32
            class(mask_expression), intent(in) :: self
            integer(int32), intent(in) :: i_word
        end function word_of_expression
    end interface

    !> A leaf of an expression: points at a real mask's words.
    type, extends(mask_expression) :: mask_leaf
        integer(int32), pointer, contiguous :: words(:) => null()
    contains
        procedure :: bit => leaf_bit
        procedure :: word => leaf_word
    end type mask_leaf

    !> `left .and. right`, evaluated on access. The children are allocatable copies of the
    !| sub-expressions: pointers to them would dangle, because an operator's operands are
    !| temporaries.
    type, extends(mask_expression) :: mask_and_node
        class(mask_expression), allocatable :: left, right
    contains
        procedure :: bit => and_node_bit
        procedure :: word => and_node_word
    end type mask_and_node

    interface operator(.and.)
        module procedure owned_and, expression_and
    end interface

    interface assignment(=)
        module procedure owned_from_flags
    end interface

contains

    ! ---------------------------------------------------------------- the core

    pure integer(int32) function words_for_bits(n_bits)
        integer(int32), intent(in) :: n_bits

        words_for_bits = (n_bits + BITS_PER_WORD - 1_int32)/BITS_PER_WORD
    end function words_for_bits

    !> Packs a logical(c_bool) array into words.
    pure subroutine core_pack(n_bits, n_words, flags, words)
        integer(int32), intent(in) :: n_bits, n_words
        logical(c_bool), intent(in) :: flags(n_bits)
        integer(int32), intent(out) :: words(n_words)
        integer(int32) :: i_word, i_bit, first_flag

        do i_word = 1, n_words
            words(i_word) = 0_int32
            first_flag = (i_word - 1_int32)*BITS_PER_WORD
            do i_bit = 0_int32, min(BITS_PER_WORD, n_bits - first_flag) - 1_int32
                if (flags(first_flag + i_bit + 1_int32)) words(i_word) = ibset(words(i_word), i_bit)
            end do
        end do
    end subroutine core_pack

    pure subroutine core_and_into(n_words, left, right, result_words)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in) :: left(n_words), right(n_words)
        integer(int32), intent(inout) :: result_words(n_words)
        integer(int32) :: i_word

        do i_word = 1, n_words
            result_words(i_word) = iand(left(i_word), right(i_word))
        end do
    end subroutine core_and_into

    pure subroutine core_copy_into(n_words, source, destination)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in) :: source(n_words)
        integer(int32), intent(inout) :: destination(n_words)
        integer(int32) :: i_word

        do i_word = 1, n_words
            destination(i_word) = source(i_word)
        end do
    end subroutine core_copy_into

    pure integer(int32) function core_count(n_words, words)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in) :: words(n_words)
        integer(int32) :: i_word

        core_count = 0_int32
        do i_word = 1, n_words
            core_count = core_count + int(popcnt(words(i_word)), int32)
        end do
    end function core_count

    !> count(left .and. right) without materialising `left .and. right`.
    pure integer(int32) function core_and_count(n_words, left, right)
        integer(int32), intent(in) :: n_words
        integer(int32), intent(in) :: left(n_words), right(n_words)
        integer(int32) :: i_word

        core_and_count = 0_int32
        do i_word = 1, n_words
            core_and_count = core_and_count + int(popcnt(iand(left(i_word), right(i_word))), int32)
        end do
    end function core_and_count

    ! ---------------------------------------------------------------- owned_mask

    !> The eager operator: every call allocates a new result.
    pure function owned_and(left, right) result(conjunction)
        type(owned_mask), intent(in) :: left, right
        type(owned_mask) :: conjunction

        allocate (conjunction%words(size(left%words)))
        conjunction%n_bits = left%n_bits
        call core_and_into(size(left%words, kind=int32), left%words, right%words, conjunction%words)
    end function owned_and

    !> `mask = flags`: the same packing loop as core_pack, behind a defined assignment.
    pure subroutine owned_from_flags(mask, flags)
        type(owned_mask), intent(inout) :: mask
        logical(c_bool), intent(in) :: flags(:)
        integer(int32) :: n_bits

        n_bits = size(flags, kind=int32)
        if (allocated(mask%words)) then
            if (size(mask%words, kind=int32) /= words_for_bits(n_bits)) deallocate (mask%words)
        end if
        if (.not. allocated(mask%words)) allocate (mask%words(words_for_bits(n_bits)))
        mask%n_bits = n_bits
        call core_pack(n_bits, size(mask%words, kind=int32), flags, mask%words)
    end subroutine owned_from_flags

    pure integer(int32) function owned_count(mask)
        type(owned_mask), intent(in) :: mask

        owned_count = core_count(size(mask%words, kind=int32), mask%words)
    end function owned_count

    ! ---------------------------------------------------------------- view_mask

    !> The same loop as core_and_into, through pointer components.
    pure subroutine view_and_into(left, right, conjunction)
        type(view_mask), intent(in) :: left, right
        type(view_mask), intent(inout) :: conjunction
        integer(int32) :: i_word

        do i_word = 1, size(left%words, kind=int32)
            conjunction%words(i_word) = iand(left%words(i_word), right%words(i_word))
        end do
    end subroutine view_and_into

    ! ---------------------------------------------------------------- mask_expression

    function leaf_of(mask) result(leaf)
        type(owned_mask), target, intent(in) :: mask
        type(mask_leaf) :: leaf

        leaf%words => mask%words
    end function leaf_of

    function expression_and(left, right) result(node)
        class(mask_expression), intent(in) :: left, right
        type(mask_and_node) :: node

        allocate (node%left, source=left)
        allocate (node%right, source=right)
    end function expression_and

    pure logical function leaf_bit(self, i_bit)
        class(mask_leaf), intent(in) :: self
        integer(int32), intent(in) :: i_bit

        leaf_bit = btest(self%words((i_bit - 1_int32)/BITS_PER_WORD + 1_int32), mod(i_bit - 1_int32, BITS_PER_WORD))
    end function leaf_bit

    pure integer(int32) function leaf_word(self, i_word)
        class(mask_leaf), intent(in) :: self
        integer(int32), intent(in) :: i_word

        leaf_word = self%words(i_word)
    end function leaf_word

    pure logical function and_node_bit(self, i_bit)
        class(mask_and_node), intent(in) :: self
        integer(int32), intent(in) :: i_bit

        and_node_bit = self%left%bit(i_bit)
        if (and_node_bit) and_node_bit = self%right%bit(i_bit)
    end function and_node_bit

    pure integer(int32) function and_node_word(self, i_word)
        class(mask_and_node), intent(in) :: self
        integer(int32), intent(in) :: i_word

        and_node_word = iand(self%left%word(i_word), self%right%word(i_word))
    end function and_node_word

end module mask_candidates
