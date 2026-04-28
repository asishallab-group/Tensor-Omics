module config
    implicit none
#ifdef DEFAULT_ALIGNMENT
    integer, parameter :: alignment = DEFAULT_ALIGNMENT
#else
    integer, parameter :: alignment = 32  ! fallback
#endif
        !! Selected alignment
    logical :: DEBUG = .false.
        !! Toggles debug messages for e.g. [[tox_archive(module)]]
    logical :: debug_hashing = .false.
        !! Toggles debug messages for [[f42_xxh3_hashmap(module)]]
end module config
