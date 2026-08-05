#include <src/macros.h>

!> Global, mutable runtime configuration flags (module variables, not parameters).
!| Intended to be toggled by client code before calling into the library to enable
!| diagnostic `print` output from the affected modules; both flags default to off.
module f42_config
    M_IMPLICIT_NONE
    logical :: DEBUG = .false.
        !! Toggles debug messages for e.g. [[tox_data_archive(module)]]
    logical :: debug_hashing = .false.
        !! Toggles debug messages for [[f42_xxh3_hashmap(module)]]
end module f42_config
