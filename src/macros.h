#include "../authors.h"

#define M_CHECK_IERR_NON_NULL if (.not. c_associated(c_loc(ierr))) return

#define M_CHECK_NON_NULL(ARG) if (.not. c_associated(c_loc(ARG))) then; call set_err(ierr, ERR_POINTER_NULL); return; endif

#define M_USE_NULL_VALIDATION use, intrinsic :: iso_c_binding, only: c_associated, c_loc; use tox_errors, only: set_err, ERR_POINTER_NULL
#define M_USE_ALLOCATION use tox_errors, only: ERR_ALLOC_FAIL, is_err, set_err

#define M_DEFAULT_VAL(OPT_ARG, LOC_VAR, DEFAULT_VAL) if (present(OPT_ARG)) then; LOC_VAR = OPT_ARG; else; LOC_VAR = DEFAULT_VAL; endif

#define M_ALLOCATE(SINGLE_VAR_DECL) allocate(SINGLE_VAR_DECL, stat=ierr); if (is_err(ierr)) then; call set_err(ierr, ERR_ALLOC_FAIL); return; endif

#define M_NAN ieee_value(1.0_real64, ieee_quiet_nan)
#define M_NEG_INF ieee_value(1.0_real64, ieee_negative_inf)
#define M_POS_INF ieee_value(1.0_real64, ieee_positive_inf)

#define M_DOC_DEFAULT(DEFAULT_VAL) The default value is `DEFAULT_VAL`.
#define M_DOC_NO_DEFAULT This argument will be ignored if not present.
#define M_DOC_REQUIRED_IF_MODE(MODULE, MODE) This optional argument needs to be passed if used mode is [[MODULE(module):MODE(variable)]].
#define M_DOC_REQUIRED_IF_METHOD(MODULE, METHOD) This optional argument needs to be passed if used method is [[MODULE(module):METHOD(variable)]].
