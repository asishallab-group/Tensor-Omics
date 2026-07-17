#include <authors.h>

! `ierr` is the only channel available to report an error back to the caller, so if its own
! pointer is null there is no way to signal failure at all -- just return silently instead of
! risking a null dereference by calling set_err(ierr, ...).
#define M_CHECK_IERR_NON_NULL if (.not. c_associated(c_loc(ierr))) return

!TODO codegen: c_loc's target must not have zero size per the Fortran standard. C wrapper subroutines call
!              M_CHECK_NON_NULL(array_arg) before the array's declared extent (e.g. n_points, n_values) has been
!              validated to be > 0, so a caller passing a legitimately empty array (size 0) hits c_loc on a zero-size
!              target here, which is not standard-conforming / is processor-dependent.
#define M_CHECK_NON_NULL(ARG) if (.not. c_associated(c_loc(ARG))) then; call set_err(ierr, ERR_POINTER_NULL); return; endif

#define M_USE_NULL_VALIDATION use, intrinsic :: iso_c_binding, only: c_associated, c_loc; use tox_errors, only: set_err, ERR_POINTER_NULL
#define M_USE_ALLOCATION use tox_errors, only: ERR_ALLOC_FAIL, is_err, set_err

#define M_DEFAULT_VAL(OPT_ARG, LOC_VAR, DEFAULT_VAL) if (present(OPT_ARG)) then; LOC_VAR = OPT_ARG; else; LOC_VAR = DEFAULT_VAL; endif

! not using is_err here, as our error code encoding might lead to unexpected behavior. The ISO standard defines success value as zero.
#define M_ALLOCATE(SINGLE_VAR_DECL) allocate(SINGLE_VAR_DECL, stat=ierr); if (ierr /= 0) then; call set_err(ierr, ERR_ALLOC_FAIL); return; endif

#define M_NAN ieee_value(1.0_real64, ieee_quiet_nan)
#define M_NEG_INF ieee_value(1.0_real64, ieee_negative_inf)
#define M_POS_INF ieee_value(1.0_real64, ieee_positive_inf)

#define M_GENE_TO_FAM_SENTINEL 0_int32
#define M_GENE_TO_FAM_DOC(GENES_TARGET_ARG) Index mapping -> each index `i` holds the family index for the corresponding gene in `GENES_TARGET_ARG`, using `M_GENE_TO_FAM_SENTINEL` for unassigned genes

! not using is_err here, as our error code encoding might lead to unexpected behavior. The ISO standard defines success value as zero.
#define M_CHECK_IO_ERR(ERR_CODE) if (ierr /= 0) then; call set_err(ierr, ERR_CODE); close(unit); return; end if

! An error code packs the position of the argument that caused it in as
! `M_ERR_ARG_POS_FACTOR*arg_pos + error`, so no `ERR_*` value may reach the factor.
! Defined here rather than written into tox_errors, because the code generator reads it
! from this file to decode `ierr` for the Python and R error modules -- a literal in the
! Fortran would mean the generator had to hardcode the same number and could drift from it.
#define M_ERR_ARG_POS_FACTOR 10000

! =============================================================================
! Documentation macros (DM_)
!
! These are written inside Ford comments (`!!` / `!|`) and expand to prose, so
! the rendered documentation reads naturally while the code generator can still
! recognise the statement. They carry the information the generator cannot infer
! from a signature alone: defaults for optionals, conditional requirements, and
! where an argument's value is meant to come from.
!
! The generator derives its patterns by expanding these very definitions, so
! rewording one here cannot desynchronise it from the generator. Only the macro
! name and its parameter order are load-bearing.
!
! Prefixes: M_ code macros (this file), CM_ file-local code macros, DM_ doc macros.
! See helper/codegen_reworked/README.md for the full contract.
! =============================================================================

! `DEFAULT_VAL` is the value an optional argument takes when it is omitted. It must
! be a constant expression, as the generator evaluates it at generation time to pass
! it on from the interfacing languages.
#define DM_DEFAULT(DEFAULT_VAL) The default value is `DEFAULT_VAL`.

! For an optional argument that has no default but is required in one specific mode:
! `MODE_ARG` names the mode argument, `MODE_PARAM` the mode parameter in `MODULE`.
#define DM_REQUIRED_IF_MODE(MODE_ARG, MODULE, MODE_PARAM) This optional argument needs to be passed if used mode (`MODE_ARG`) is [[MODULE(module):MODE_PARAM(variable)]].

! An `intent(out)` argument the caller may decline to receive.
#define DM_OPTIONAL_OUTPUT This output will only be present if desired.

! For a result array that is filled only partially: `ARGUMENT` names the scalar
! integer argument holding how many leading elements actually carry results.
#define DM_RESULT_SIZE_IS(ARGUMENT) The first `ARGUMENT` elements will hold the results.

! The value of this argument comes from another procedure, typically a work array
! size that cannot be foreseen. `MODE` is AUTO when the interfacing languages should
! call `PROCEDURE` themselves, or JUST_INFO when the caller has to do it.
#define DM_OUTPUT_FROM(ARGUMENT, PROCEDURE, MODULE, MODE) DM_OUTPUT_FROM_##MODE to compute this argument from the `ARGUMENT` output of [[MODULE(module):PROCEDURE]].
#define DM_OUTPUT_FROM_AUTO It is *VERY IMPORTANT*
#define DM_OUTPUT_FROM_JUST_INFO It is recommended
