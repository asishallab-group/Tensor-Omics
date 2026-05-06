#include "authors.h"

!> AUTHOR_FRANZ_ERIC_SILL
!| Main array deserialization module, re-exporting all deserialization routines
module f42_serde_arrays_deserialize
    use f42_serde_arrays_deserialize_int
    use f42_serde_arrays_deserialize_char
    use f42_serde_arrays_deserialize_real
    use f42_serde_arrays_deserialize_logical
    use f42_serde_arrays_deserialize_complex
end module f42_serde_arrays_deserialize