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



#define DM_OPTIONAL_OUTPUT This output will only be present if desired.
#define DM_REQUIRED_IF_MODE(MODE_VAR, MODULE, MODE) This optional argument needs to be passed if used mode (`MODE_VAR`) is [[MODULE(module):MODE(variable)]].
#define DM_DEFAULT(DEFAULT_VAL) The default value is `DEFAULT_VAL`.


#define DM_OUTPUT_FROM(ARGUMENT, PROCEDURE, MODULE, MODE) DM_OUTPUT_FROM_##MODE to compute this argument using [[MODULE(module):PROCEDURE]]'s output `ARGUMENT`.
#define DM_OUTPUT_FROM_AUTO It is *VERY IMPORTANT*
#define DM_OUTPUT_FROM_JUST_INFO It is recommended

#define DM_RESULT_SIZE_IS(ARGUMENT) The first `ARGUMENT` elements will hold the results.
