#include <authors.h>

!> AUTHOR_FRANZ_ERIC_SILL
!| Main array serialization module, re-exporting all serialization routines
module f42_serde_arrays_serialize
    use f42_serde_arrays_serialize_int
    use f42_serde_arrays_serialize_char
    use f42_serde_arrays_serialize_real
    use f42_serde_arrays_serialize_logical
    use f42_serde_arrays_serialize_complex
end module f42_serde_arrays_serialize