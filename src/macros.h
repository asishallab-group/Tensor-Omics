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


! Documentation macros

! Use only for optional intent(out) arguments, if it is unavoidable to have optional outputs.
#define DM_OPTIONAL_OUTPUT This output will only be present if desired.
! If an argument must be present, depending on the used mode, use this macro
#define DM_REQUIRED_IF_MODE(MODE_VAR, MODULE, MODE) This optional argument needs to be passed if used mode (`MODE_VAR`) is [[MODULE(module):MODE(variable)]].
! If the presence of the input argument is not dependent on other arguments, it always needs a default. Use Fortran expressions as DEFAULT_VAL
#define DM_DEFAULT(DEFAULT_VAL) The default value is `DEFAULT_VAL`.


! Use DM_FROM(..., AUTO) if ARGUMENT is output of another procedure and should be calculated automatically in Python/R, like functions for work array size calculations or other generic things.
! Use DM_FROM(..., JUST_INFO) if ARGUMENT is output of another procedure, but should not be calculated automatically, like normalized inputs coming from normalization functions. It is up to the user to use them beforehand.
#define M_DM_FROM_ DM_FROM_
#define DM_FROM(ARGUMENT, PROCEDURE, MODULE, MODE) M_DM_FROM_##MODE to compute this argument using [[MODULE(module):PROCEDURE]]'s output `ARGUMENT`.
#define DM_FROM_AUTO It is *VERY IMPORTANT*
#define DM_FROM_JUST_INFO It is recommended

! Use for output array arguments whose result count is specified by another output argument, like `output(n), n_results` where `0<=n_results=n`, so actual output is slice `output(:n_results)`
! Note that `ARGUMENT` should be related to the last extent `output(:, :, ..., :ARGUMENT)`.
#define DM_RESULT_SIZE_IS(ARGUMENT) The first `ARGUMENT` elements will hold the results.
