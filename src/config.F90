module config
    implicit none
    logical :: DEBUG = .false.
        !! Toggles debug messages for e.g. [[tox_data_archive(module)]]
    logical :: debug_hashing = .false.
        !! Toggles debug messages for [[f42_xxh3_hashmap(module)]]
end module config
