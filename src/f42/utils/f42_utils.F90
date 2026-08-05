#include <src/macros.h>

!> Utility module for data analysis.
!| This module provides general-purpose utility functions for data analysis, to be used as needed.
!|
!| It holds nothing itself: it re-exports the modules below, which is the whole of its job.
!| `use f42_utils` therefore still reaches every utility, and which module a given one is
!| defined in stays an implementation detail.
module f42_utils
    use f42_sort
    use f42_math
    use f42_random
    use f42_vector
    use f42_stats
end module f42_utils
