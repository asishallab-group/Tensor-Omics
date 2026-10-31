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
! Characters are viewed rather than copied; see M_CHECK_CHARACTER_VIEW.
#define M_CHECK_ARRAY_NON_NULL(ARG, N) if ((N) > 0) then; M_CHECK_NON_NULL(ARG); end if

! A character buffer is not waved through as the array check above waves every other kind.
! Every other kind is copied or allocated into a local; a string is a *view*, and the
! `c_f_pointer(c_loc(ARG), ...)` that builds it has to run before the callee is called.
! `c_loc` may not be given a zero-size target, and an unassociated pointer may not be passed
! to a non-pointer dummy -- gfortran calls that one out ("Pointer actual argument is not
! associated") and ifx segfaults.
!
! A single string of zero width has no address to view, so the wrapper reaches the verdict
! validate_dimension_size would for a zero extent, and returns it. `N` is the width. (Python
! and R never send one: they pad an empty string to one blank.)
#define M_CHECK_CHARACTER_VIEW(ARG, N) if ((N) <= 0) then; call set_err(ierr, ERR_EMPTY_INPUT); return; endif

! An array of strings may hold none -- a data set with zero families -- and is then viewed at
! the C wrapper module's `character(1), target` anchor with shape [0], so the callee decides
! about it as about any empty array. What has no address to view is a zero `WIDTH` with
! strings in it, and a negative `COUNT` is no count at all; both are refused.
#define M_CHECK_CHARACTER_ARRAY_VIEW(ARG, WIDTH, COUNT) if ((COUNT) < 0 .or. ((COUNT) > 0 .and. (WIDTH) <= 0)) then; call set_err(ierr, ERR_EMPTY_INPUT); return; endif

! A plain `implicit none` constrains only *variables*. A call to a procedure that does not
! exist -- a typo, or a helper that was renamed -- is still accepted as an implicit
! external and fails at link time, with a message naming a symbol rather than a line. The
! `(type, external)` form (F2018) makes it a compile error at the call site instead.
!
! Use this instead of a bare `implicit none` in every module.
#define M_IMPLICIT_NONE implicit none (type, external)

#define M_DEFAULT_VAL(OPT_ARG, LOC_VAR, DEFAULT_VAL) if (present(OPT_ARG)) then; LOC_VAR = OPT_ARG; else; LOC_VAR = DEFAULT_VAL; endif

! not using is_err here, as our error code encoding might lead to unexpected behavior. The ISO standard defines success value as zero.
#define M_ALLOCATE(SINGLE_VAR_DECL) allocate(SINGLE_VAR_DECL, stat=ierr); if (ierr /= 0) then; call set_err(ierr, ERR_ALLOC_FAIL); return; endif

#define M_NAN ieee_value(1.0_real64, ieee_quiet_nan)
#define M_NEG_INF ieee_value(1.0_real64, ieee_negative_inf)
#define M_POS_INF ieee_value(1.0_real64, ieee_positive_inf)

#define M_GENE_TO_FAM_SENTINEL 0_int32
#define M_GENE_TO_FAM_DOC(GENES_TARGET_ARG) Index mapping -> each index `i` holds the family index for the corresponding gene in `GENES_TARGET_ARG`, using `M_GENE_TO_FAM_SENTINEL` for unassigned genes

! OpenMP (#218). `_OPENMP` is defined by the compiler, not here: `./build.sh --openmp` passes
! the OpenMP flag, which defines it, and that is the one switch. Nothing else enables OpenMP.
!
! The purity of a procedure that holds a parallel loop: an OpenMP directive is not allowed in
! a pure procedure, so it is pure without OpenMP and impure with it. Written where `pure` would
! be, in the procedure and in every wrapper that calls it.
!
! The thread count the allocating wrapper sizes the per-thread work arrays with and passes on:
! the caller's `N` with OpenMP, and one without, where a single thread runs the loop.
#ifdef _OPENMP
#define M_OMP_PURE impure
#define M_OMP_N_THREADS(N) N
#else
#define M_OMP_PURE pure
#define M_OMP_N_THREADS(N) 1_int32
#endif

! The parallel loop of a root procedure, so no clause is forgotten. Written as a comment line
! directly above the `do`, with the `!` in the source and the loop's locality clauses after it
! on the same line, e.g. `!M_OMP_PARALLEL(n_items > 1) shared(values, sums) private(i_thread)`.
! It expands to an OpenMP directive, which a build without OpenMP reads as a comment. Threads
! run only when the caller allows more than one (`n_omp_threads`, an argument of the root) and
! the loop has enough to share (`CONDITION`); `default(none)` makes every variable the loop
! touches need a clause. The `!` cannot move into the body: ifx's preprocessor reads it as the
! start of a comment, and the macro expands to nothing.
#define M_OMP_PARALLEL(CONDITION) $omp parallel do default(none) num_threads(n_omp_threads) if(n_omp_threads > 1 .and. (CONDITION)) schedule(dynamic)

! not using is_err here, as our error code encoding might lead to unexpected behavior. The ISO standard defines success value as zero.
#define M_CHECK_IO_ERR(ERR_CODE) if (ierr /= 0) then; call set_err(ierr, ERR_CODE); close(unit); return; end if

! An error code packs the position of the argument that caused it in as
! `M_ERR_ARG_POS_FACTOR*arg_pos + error`, so no `ERR_*` value may reach the factor.
! Defined here rather than written into tox_errors, because the code generator reads it
! from this file to decode `ierr` for the Python and R error modules -- a literal in the
! Fortran would mean the generator had to hardcode the same number and could drift from it.
#define M_ERR_ARG_POS_FACTOR 10000

! Build switches the generator's output honours, defined by the *build* rather than here:
!   NO_C_BINDING        omit the bind(C) wrappers, and with them the R shims
!   NO_R_BINDING        omit the R shims alone
!   NO_INPUT_VALIDATION omit the generated wrappers' input checks
! The last is for a caller who has already established that the inputs are good -- an inner
! loop over data it produced itself. It gives up every diagnostic the framework offers, so it
! is a whole-build decision rather than a per-call one. `call set_ok(ierr)` survives it: that
! is not a check, it is what leaves `ierr` defined, and an implementation's own runtime errors are
! still reported. Pass it with `./build.sh --directive=NO_INPUT_VALIDATION`.

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

! A prologue is the sugar `foo` adds over `foo_expert`: it runs after the work arrays are
! prepared and before the implementation, and may handle the call itself -- writing the
! outputs and reporting `handled`, so the implementation is skipped. It derives what the
! expert tier lets a caller pass in, exactly as the `<base>_perm` convention seeds and
! heapsorts a permutation.
! There is no scope. Work that every wrapper needs is work every *caller of the
! implementation* needs, and both wrappers call it -- so it belongs at the top of the
! implementation, where it needs no directive at all. An implementation with no work arrays
! generates only `foo`, so a prologue on one would never run, and that is an error.
! The prologue's dummies are supplied by name from the implementation's arguments -- the work
! arrays included -- plus `handled` and `ierr`; a name that matches none of these becomes an
! argument of `foo` itself, which is how a prologue asks for what it derives from. A name
! that is merely *close* to one of them is an error, so rename the prologue's dummy to
! whatever the implementation calls the same thing. It must declare
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
