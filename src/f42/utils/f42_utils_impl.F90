#include <src/macros.h>

!> General-purpose utilities for data analysis: elementary mathematics, indirect sorting and
!| searching, randomness, vector geometry, and descriptive statistics.
!|
!| It holds nothing itself and only gathers the modules below, so one `use` reaches every
!| utility and which module a given one is defined in stays out of the way.
module f42_utils_impl
    use f42_sort_impl
    use f42_math_impl
    use f42_random_impl
    use f42_vector_impl
    use f42_stats_impl
end module f42_utils_impl
