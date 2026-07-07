#include <authors.h>

! `ierr` is the only channel available to report an error back to the caller, so if its own
! pointer is null there is no way to signal failure at all -- just return silently instead of
! risking a null dereference by calling set_err(ierr, ...).
#define M_CHECK_IERR_NON_NULL if (.not. c_associated(c_loc(ierr))) return

!CLAUDE: c_loc's target must not have zero size per the Fortran standard. C wrapper subroutines call
!CLAUDE: M_CHECK_NON_NULL(array_arg) before the array's declared extent (e.g. n_points, n_values) has been
!CLAUDE: validated to be > 0, so a caller passing a legitimately empty array (size 0) hits c_loc on a zero-size
!CLAUDE: target here, which is not standard-conforming / is processor-dependent.
#define M_CHECK_NON_NULL(ARG) if (.not. c_associated(c_loc(ARG))) then; call set_err(ierr, ERR_POINTER_NULL); return; endif

#define M_USE_NULL_VALIDATION use, intrinsic :: iso_c_binding, only: c_associated, c_loc; use tox_errors, only: set_err, ERR_POINTER_NULL
#define M_USE_ALLOCATION use tox_errors, only: ERR_ALLOC_FAIL, is_err, set_err

#define M_DEFAULT_VAL(OPT_ARG, LOC_VAR, DEFAULT_VAL) if (present(OPT_ARG)) then; LOC_VAR = OPT_ARG; else; LOC_VAR = DEFAULT_VAL; endif

!CLAUDE: `stat=ierr` writes a raw, compiler-defined Fortran allocation-status integer into `ierr`, which
!CLAUDE: `is_err` then decodes via `mod(ierr,10000)` (this module's own ERR_* encoding scheme). This only works
!CLAUDE: by coincidence for small stat values; a compiler returning a stat that happens to be a multiple of
!CLAUDE: 10000 would be misreported as success. Same class of issue as M_CHECK_IO_ERR below.
#define M_ALLOCATE(SINGLE_VAR_DECL) allocate(SINGLE_VAR_DECL, stat=ierr); if (is_err(ierr)) then; call set_err(ierr, ERR_ALLOC_FAIL); return; endif

#define M_NAN ieee_value(1.0_real64, ieee_quiet_nan)
#define M_NEG_INF ieee_value(1.0_real64, ieee_negative_inf)
#define M_POS_INF ieee_value(1.0_real64, ieee_positive_inf)

#define M_GENE_TO_FAM_SENTINEL 0_int32
#define M_GENE_TO_FAM_DOC(GENES_TARGET_ARG) Index mapping -> each index `i` holds the family index for the corresponding gene in `GENES_TARGET_ARG`, using `M_GENE_TO_FAM_SENTINEL` for unassigned genes

!CLAUDE: same fragile coupling as M_ALLOCATE above: this treats a raw `iostat` value already sitting in `ierr`
!CLAUDE: as if it were already encoded via this module's `mod(ierr,10000)` error scheme, rather than mapping it
!CLAUDE: through set_err/create_err_code first.
#define M_CHECK_IO_ERR(ERR_CODE) if (is_err(ierr)) then; call set_err(ierr, ERR_CODE); close(unit); return; end if
