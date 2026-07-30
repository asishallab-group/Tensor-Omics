!> AUTHOR_AARON_SCHROEDER
!| String-keyed hashmap (`character -> integer(int32)`) and hashset (`character` membership),
!| implemented as power-of-two-sized bucket arrays with separate chaining and the external XXH3
!| algorithm for hashing. Both containers grow automatically via [[f42_xxh3_hashmap(module):resize_hashmap(subroutine)]]
!| / [[f42_xxh3_hashmap(module):resize_hashset(subroutine)]] once the load factor exceeds
!| [[f42_xxh3_hashmap(module):MAX_LOAD_FACTOR(variable)]].
module f42_xxh3_hashmap
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_loc
    use iso_fortran_env, only: int32, int64
    use f42_utils, only: next_power_of_two
    use config, only: DEBUG, debug_hashing
    use tox_errors, only: set_ok, set_err, ERR_INVALID_INPUT
    implicit none
    private

    public :: hashmap_type, hashmap_create, hashmap_destroy, hashmap_get, hashmap_put
    public :: hashset_type, hashset_create, hashset_destroy, hashset_put, is_in_hashset

    !> External XXH3 64-bit hash function (from the bundled xxhash C library), used by
    !| [[f42_xxh3_hashmap(module):xxh3_hash_fortran(function)]] to hash string keys.
    interface
        function xxh3_hash_c(key, length) bind(C, name="XXH3_64bits")
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_int64_t
            implicit none
            type(c_ptr), value :: key
                !! C pointer to the start of the key's character data.
            integer(c_int), value :: length
                !! Number of bytes to hash starting at `key`.
            integer(c_int64_t) :: xxh3_hash_c
                !! The resulting 64-bit hash value.
        end function xxh3_hash_c
    end interface

    !> Singly-linked list node used for separate chaining in [[f42_xxh3_hashmap(module):hashmap_type(type)]].
    type :: hashmap_node_type
        character(len=:), allocatable :: key
            !! The stored (trimmed) string key.
        integer(int32) :: value
            !! The value associated with `key`.
        type(hashmap_node_type), pointer :: next => null()
            !! Next node in this bucket's chain, or null if this is the last node.
    end type hashmap_node_type

    !> Singly-linked list node used for separate chaining in [[f42_xxh3_hashmap(module):hashset_type(type)]].
    type :: hashset_node_type
        character(len=:), allocatable :: key
            !! The stored (trimmed) string key.
        type(hashset_node_type), pointer :: next => null()
            !! Next node in this bucket's chain, or null if this is the last node.
    end type hashset_node_type

    ! Parameters
    integer, parameter :: DEFAULT_KEY_LENGTH = 256
        !! Reserved for a future fixed-length key buffer; currently unused (keys are stored as
        !! allocatable strings, not fixed-length buffers).
    real, parameter :: MAX_LOAD_FACTOR = 0.75
        !! Fraction of buckets that may hold at least one entry (`count/size`) before a `put`
        !! triggers a doubling resize. 0.75 is the conventional trade-off for chained hashing:
        !! low enough to keep chain lengths near O(1), high enough to not waste much bucket space.

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Base type for hashmap and hashset (currently)
    type :: hashmap_base
        integer(int32) :: size = 0
            !! Number of buckets (always a power of two, see [[f42_utils(module):next_power_of_two(function)]]).
        integer(int32) :: count = 0
            !! Number of key-value pairs currently stored.
    end type hashmap_base

    !> Hashmap of `character` keys to `integer(int32)` values, using separate chaining over a
    !| power-of-two-sized bucket array.
    type, extends(hashmap_base) :: hashmap_type
        type(hashmap_node_type), pointer :: buckets(:) => null()
            !! Bucket array; `buckets(i)%next` heads the chain for bucket `i` (the bucket
            !! elements themselves only ever use their `next` pointer, never `key`/`value`).
    end type hashmap_type

    !> Hashset of `character` keys, using separate chaining over a power-of-two-sized bucket array.
    type, extends(hashmap_base) :: hashset_type
        type(hashset_node_type), pointer :: buckets(:) => null()
            !! Bucket array; `buckets(i)%next` heads the chain for bucket `i` (the bucket
            !! elements themselves only ever use their `next` pointer, never `key`/`value`).
    end type hashset_type

contains

    !> AUTHOR_AARON_SCHROEDER
    !| Create the hashmap
    subroutine hashmap_create(map, initial_size)
        type(hashmap_type), intent(out) :: map
            !! Hashmap object to create
        integer(int32), intent(in), optional :: initial_size
            !! Desired initial element capacity; the actual bucket count is derived from this
            !! (scaled by MAX_LOAD_FACTOR and rounded up to a power of two, floored at 128),
            !! default: 1024 buckets

        call hashmap_base_create(map, initial_size)

        ! Allocate buckets
        allocate (map%buckets(map%size))

        if (DEBUG) print *, "Hashmap created with size:", map%size
    end subroutine hashmap_create

    !> AUTHOR_AARON_SCHROEDER
    !| Create the hashset
    subroutine hashset_create(set, initial_size)
        type(hashset_type), intent(out) :: set
            !! Hashset object to create
        integer(int32), intent(in), optional :: initial_size
            !! Desired initial element capacity; the actual bucket count is derived from this
            !! (scaled by MAX_LOAD_FACTOR and rounded up to a power of two, floored at 128),
            !! default: 1024 buckets

        call hashmap_base_create(set, initial_size)

        ! Allocate buckets
        allocate (set%buckets(set%size))

        if (DEBUG) print *, "Hashset created with size:", set%size
    end subroutine hashset_create

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper for creating hashmap-likes
    subroutine hashmap_base_create(map, initial_size)
        class(hashmap_base), intent(out) :: map
            !! Hashmap object to create
        integer(int32), intent(in), optional :: initial_size
            !! Desired initial element capacity; the actual bucket count is derived from this
            !! (scaled by MAX_LOAD_FACTOR and rounded up to a power of two, floored at 128),
            !! default: 1024 buckets

        integer(int32) :: table_size, i

        ! Calculate table size (power of two)
        if (present(initial_size)) then
            table_size = next_power_of_two(int(initial_size/MAX_LOAD_FACTOR, int32))
        else
            table_size = 1024  ! Default size
        end if

        ! Ensure minimum size
        table_size = max(table_size, 128)

        map%size = table_size
        map%count = 0
    end subroutine hashmap_base_create

    !> AUTHOR_AARON_SCHROEDER
    !| Destroy the hashmap
    subroutine hashmap_destroy(map)
        type(hashmap_type), intent(inout) :: map
            !! The hashmap to delete

        integer(int32) :: i
        type(hashmap_node_type), pointer :: current, next

        if (associated(map%buckets)) then
            do i = 1, map%size
                current => map%buckets(i)%next
                do while (associated(current))
                    next => current%next
                    deallocate (current)
                    current => next
                end do
            end do
            deallocate (map%buckets)
            nullify (map%buckets)
        end if

        map%size = 0
        map%count = 0
    end subroutine hashmap_destroy

    !> AUTHOR_AARON_SCHROEDER
    !| Destroy the set
    subroutine hashset_destroy(set)
        type(hashset_type), intent(inout) :: set
            !! The set to delete

        integer(int32) :: i
        type(hashset_node_type), pointer :: current, next

        if (associated(set%buckets)) then
            do i = 1, set%size
                current => set%buckets(i)%next
                do while (associated(current))
                    next => current%next
                    deallocate (current)
                    current => next
                end do
            end do
            deallocate (set%buckets)
            nullify (set%buckets)
        end if

        set%size = 0
        set%count = 0
    end subroutine hashset_destroy

    !> AUTHOR_AARON_SCHROEDER
    !| Compute XXH3 hash of a string
    function xxh3_hash_fortran(key, table_size) result(hash_idx)
        character(len=*), intent(in) :: key
            !! key to hash
        integer(int32), intent(in) :: table_size
            !! table size
        integer(int32) :: hash_idx
            !! resulting hash

        integer(int64) :: hash_val
        integer :: key_len
        character(len=:), allocatable, target :: trimmed_key
        trimmed_key = trim(key)

        key_len = len(trimmed_key)

        ! An empty/all-blank key yields a zero-length `trimmed_key`; taking c_loc of a zero-size
        ! target is not well-defined, so special-case it with a fixed hash instead of calling into C.
        if (key_len == 0) then
            hash_val = 0_int64
        else
            hash_val = xxh3_hash_c(c_loc(trimmed_key), key_len)
        end if
        ! table_size is always a power of two (see next_power_of_two in hashmap_create/hashset_create),
        ! so `hash mod table_size` can be computed as a bitmask AND against (table_size - 1) instead of
        ! the more expensive integer modulo/division.
        hash_idx = int(iand(hash_val, int(table_size - 1, int64)) + 1, int32)
    end function xxh3_hash_fortran

    !> AUTHOR_AARON_SCHROEDER
    !| Insert a key-value pair
    subroutine hashmap_put(map, key, value)
        type(hashmap_type), intent(inout) :: map
            !! hashmap to insert into
        character(len=*), intent(in) :: key
            !! Key to store
        integer(int32), intent(in) :: value
            !! value to store

        integer(int32) :: hash_idx
        type(hashmap_node_type), pointer :: new_node, current
        character(len=:), allocatable :: normalized_key

        ! Normalize key
        normalized_key = trim(key)

        if (debug_hashing) print *, "PUT: ", trim(normalized_key), " -> ", value

        ! Check if we need to resize
        if (real(map%count)/real(map%size) > MAX_LOAD_FACTOR) then
            call resize_hashmap(map)
        end if

        ! Compute hash index
        hash_idx = xxh3_hash_fortran(normalized_key, map%size)

        if (debug_hashing) print *, "  Hash index:", hash_idx

        ! Check if key already exists
        current => map%buckets(hash_idx)%next
        do while (associated(current))
            if (current%key == normalized_key) then
                ! Key exists - update value
                current%value = value
                if (DEBUG) print *, "Warning: Duplicate key updated value: ", normalized_key
                return
            end if
            current => current%next
        end do

        ! Key doesn't exist - create new node
        allocate (new_node)
        new_node%key = normalized_key
        new_node%value = value
        new_node%next => map%buckets(hash_idx)%next
        map%buckets(hash_idx)%next => new_node
        map%count = map%count + 1

        if (debug_hashing) print *, "  Added new key, count:", map%count
    end subroutine hashmap_put

    !> AUTHOR_AARON_SCHROEDER
    !| Insert a key-value pair
    subroutine hashset_put(set, key, ierr)
        type(hashset_type), intent(inout) :: set
            !! hashmap to insert into
        character(len=*), intent(in) :: key
            !! Key to store
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: hash_idx
        type(hashset_node_type), pointer :: new_node, current
        call set_ok(ierr)

        if (debug_hashing) print *, "PUT: ", trim(key)

        ! Check if we need to resize
        if (real(set%count)/real(set%size) > MAX_LOAD_FACTOR) then
            call resize_hashset(set)
        end if

        ! Compute hash index
        hash_idx = xxh3_hash_fortran(key, set%size)

        if (debug_hashing) print *, "  Hash index:", hash_idx

        ! Check if key already exists
        current => set%buckets(hash_idx)%next
        do while (associated(current))
            if (current%key == key) then
                call set_err(ierr, ERR_INVALID_INPUT)
                if (DEBUG) print *, "Warning: Duplicate key"
                return
            end if
            current => current%next
        end do

        ! Key doesn't exist - create new node
        allocate (new_node)
        new_node%key = trim(key)
        new_node%next => set%buckets(hash_idx)%next
        set%buckets(hash_idx)%next => new_node
        set%count = set%count + 1

        if (debug_hashing) print *, "  Added new key, count:", set%count
    end subroutine hashset_put

    !> AUTHOR_AARON_SCHROEDER
    !| Checks whether `key` is a member of `hashset`.
    logical function is_in_hashset(hashset, key) result(res)
        type(hashset_type), intent(in) :: hashset
            !! Hashset to query
        character(len=*), intent(in) :: key
            !! Key to look for
        type(hashset_node_type), pointer :: current
        integer(int32) :: hash_idx

        hash_idx = xxh3_hash_fortran(key, hashset%size)

        res = .false.
        current => hashset%buckets(hash_idx)%next
        do while (associated(current))
            if (current%key == key) then
                res = .true.
                if (debug_hashing) print *, "FOUND"
                return
            end if
            current => current%next
        end do
    end function

    !> AUTHOR_AARON_SCHROEDER
    !| Lookup a key
    function hashmap_get(map, key) result(value)
        type(hashmap_type), intent(in) :: map
            !! hashmap
        character(len=*), intent(in) :: key
            !! key to look for
        integer(int32) :: value
            !! return value

        integer(int32) :: hash_idx
        type(hashmap_node_type), pointer :: current
        character(len=:), allocatable :: normalized_key

        value = -1  ! Default: not found

        ! Normalize key
        normalized_key = trim(key)

        if (debug_hashing) print *, "GET: ", trim(normalized_key)

        ! Compute hash index
        hash_idx = xxh3_hash_fortran(normalized_key, map%size)

        if (debug_hashing) print *, "  Hash index:", hash_idx

        ! Search for the key in the bucket
        current => map%buckets(hash_idx)%next
        do while (associated(current))
            if (current%key == key) then
                value = current%value
                if (debug_hashing) print *, "  FOUND: ", value
                return
            end if
            current => current%next
        end do

        if (DEBUG) print *, "  NOT FOUND"
    end function hashmap_get

    !> AUTHOR_AARON_SCHROEDER
    !| Resize the hashmap when load factor is too high
    subroutine resize_hashmap(map)
        type(hashmap_type), intent(inout) :: map
            !! map to resize

        type(hashmap_type) :: new_map
        integer(int32) :: i, new_size
        type(hashmap_node_type), pointer :: current, next

        if (DEBUG) print *, "Resizing hashmap from ", map%size, " to ", map%size*2

        ! Create a new larger map
        new_size = map%size*2
        call hashmap_create(new_map, new_size)

        ! Reinsert all elements
        do i = 1, map%size
            current => map%buckets(i)%next
            do while (associated(current))
                call hashmap_put(new_map, current%key, current%value)
                next => current%next
                deallocate (current)
                current => next
            end do
        end do

        ! Replace the old map with the new one
        deallocate (map%buckets)
        map%size = new_map%size
        map%count = new_map%count
        map%buckets => new_map%buckets

        ! Nullify the new_map's buckets to prevent deallocation
        nullify (new_map%buckets)
    end subroutine resize_hashmap

    !> AUTHOR_AARON_SCHROEDER
    !| Resize the hashmap when load factor is too high
    subroutine resize_hashset(set)
        type(hashset_type), intent(inout) :: set
            !! set to resize

        type(hashset_type) :: new_set
        integer(int32) :: i, new_size, ierr
        type(hashset_node_type), pointer :: current, next
        call set_ok(ierr)

        if (DEBUG) print *, "Resizing hashmap from ", set%size, " to ", set%size*2

        ! Create a new larger set
        new_size = set%size*2
        call hashset_create(new_set, new_size)

        ! Reinsert all elements
        do i = 1, set%size
            current => set%buckets(i)%next
            do while (associated(current))
                call hashset_put(new_set, current%key, ierr)
                next => current%next
                deallocate (current)
                current => next
            end do
        end do

        ! Replace the old set with the new one
        deallocate (set%buckets)
        set%size = new_set%size
        set%count = new_set%count
        set%buckets => new_set%buckets

        ! Nullify the new_map's buckets to prevent deallocation
        nullify (new_set%buckets)
    end subroutine resize_hashset

end module f42_xxh3_hashmap
