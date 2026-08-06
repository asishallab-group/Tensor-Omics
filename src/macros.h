#include <authors.h>

! `ierr` is the only channel available to report an error back to the caller, so if its own
! pointer is null there is no way to signal failure at all -- just return silently instead of
! risking a null dereference by calling set_err(ierr, ...).
#define M_CHECK_IERR_NON_NULL if (.not. c_associated(c_loc(ierr))) return

#define M_CHECK_NON_NULL(ARG) if (.not. c_associated(c_loc(ARG))) then; call set_err(ierr, ERR_POINTER_NULL); return; endif

! `c_loc` may not be given a zero-size target, so an array cannot be null-checked until
! the extents that size it are known -- which is what the TODO above is about. The
! generator emits the checks in an order that makes this safe (ierr, then every scalar,
! then the arrays), and guards each array with the element count it has just made
! readable. `N` is that count, an expression over already-checked scalars.
!
! An empty array is therefore left alone: a caller passing null for a legitimately
! zero-size array gets through, and the callee's own validate_dimension_size decides
! whether empty is an error for that routine -- which is where that policy already lives.
#define M_CHECK_ARRAY_NON_NULL(ARG, N) if ((N) > 0) then; M_CHECK_NON_NULL(ARG); end if

! A plain `implicit none` constrains only *variables*. A call to a procedure that does not
! exist -- a typo, or a helper that was renamed -- is still accepted as an implicit
! external and fails at link time, with a message naming a symbol rather than a line. The
! `(type, external)` form (F2018) makes it a compile error at the call site instead.
!
! Use this instead of a bare `implicit none` in every module.
#define M_IMPLICIT_NONE implicit none (type, external)

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

! Marks a procedure for export to C, Python and R. Written in the procedure's Ford
! pre-comment: `!> M_EXPORT_C`. Expands to a Ford `category` meta tag, so Ford still parses
! it, and the code generator reads the category value from this macro rather than
! hardcoding the string -- change it here and both the sources and the generator follow.
#define M_EXPORT_C category: C-binding

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
! codegen_guide.md (repository root) is the author's guide; helper/codegen/README.md is the
! generator's own reference.
! =============================================================================

! `DEFAULT_VAL` is the value an optional argument takes when it is omitted. It must
! be a constant expression, as the generator evaluates it at generation time to pass
! it on from the binding languages.
#define DM_DEFAULT(DEFAULT_VAL) The default value is `DEFAULT_VAL`.

! For an optional argument that has no default but is required in one specific mode:
! `MODE_ARG` names the mode argument, `MODE_PARAM` the mode parameter in `MODULE`.
#define DM_REQUIRED_IF_MODE(MODE_ARG, MODULE, MODE_PARAM) This optional argument needs to be passed if used mode (`MODE_ARG`) is [[MODULE(module):MODE_PARAM(variable)]].

! For a result array that is filled only partially: `ARGUMENT` names the scalar
! integer argument holding how many leading elements actually carry results.
#define DM_RESULT_SIZE_IS(ARGUMENT) The first `ARGUMENT` elements will hold the results.

! The value of this argument comes from another procedure, typically a work array
! size that cannot be foreseen. `MODE` is AUTO when the binding languages should
! call `PROCEDURE` themselves, or JUST_INFO when the caller has to do it.
#define DM_OUTPUT_FROM(ARGUMENT, PROCEDURE, MODULE, MODE) DM_OUTPUT_FROM_##MODE to compute this argument from the `ARGUMENT` output produced by [[MODULE(module):PROCEDURE]].
#define DM_OUTPUT_FROM_AUTO It is *VERY IMPORTANT*
#define DM_OUTPUT_FROM_JUST_INFO It is recommended

! A prologue is the sugar the *allocating* wrapper adds: it runs after the work arrays are
! prepared and before the kernel, and may handle the call itself -- writing the outputs and
! reporting `handled`, so the kernel is skipped. It derives what the expert tier lets a
! caller pass in, exactly as the `<base>_perm` convention seeds and heapsorts a permutation.
! There is no scope. Work that every wrapper needs is work every *caller of the kernel*
! needs, and both wrappers call the kernel -- so it belongs at the top of the kernel, where
! it needs no directive at all. A kernel with no work arrays generates no allocating wrapper,
! so a prologue on one would never run, and that is an error.
! The prologue's dummies are supplied by name from the kernel's arguments -- the work arrays
! included -- plus `handled` and `ierr`; a name that matches neither is an error, so rename
! the prologue's dummy to whatever the kernel calls the same thing. It must declare
! `logical, intent(out) :: handled` and set it on every path.
! It runs below the allocations, so it may not produce anything they, the permutation sorts
! or the recommend calls above it read.
#define DM_PROLOGUE(PROCEDURE, MODULE) The allocating wrapper runs [[MODULE(module):PROCEDURE]] first, which may handle the call and skip this one.

! `DM_MIN` / `DM_MAX` document the inclusive valid range of a numeric argument; the
! generator turns them into a `validate_in_range_*` call in the generated wrapper. `EXPR`
! is Fortran source and may refer to other arguments or module constants; wrap it in
! `above(...)` / `below(...)` for an exclusive bound.
#define DM_MIN(EXPR) The minimum valid value is `EXPR`.
#define DM_MAX(EXPR) The maximum valid value is `EXPR`.

! A value accepted regardless of the range -- e.g. an "unassigned" marker that is not a
! real datum. Passed through to the validator's `sentinel=` argument.
#define DM_SENTINEL(EXPR) The value `EXPR` is additionally accepted.

! Finiteness is the framework's default contract: the generated validation rejects NaN and
! infinity for every real argument. These opt one argument out of that rejection, each
! failure mode separately, so an argument that legitimately carries NaN (say a masked-out
! mean) says so where the tolerance actually lives.
#define DM_ALLOW_NAN NaN is permitted for this value.
#define DM_ALLOW_INFINITE Infinite values are permitted for this value.
