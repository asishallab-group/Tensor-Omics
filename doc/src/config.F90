#include "precompiler_constants.F90"

module config
  implicit none
#ifdef DEFAULT_ALIGNMENT
  integer, parameter :: alignment = DEFAULT_ALIGNMENT
#else
  integer, parameter :: alignment = 32  ! fallback
#endif
end module config
