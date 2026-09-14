/* Correctness prototype — pure-C .Call marshalling of the awkward types.
 *
 * Proves that the two conversions the Rcpp path centralises in tox_marshal.h
 * (R int-logicals <-> Fortran c_bool bytes, and R strings <-> Fortran fixed-width
 * column-major char buffers) can be done in plain C, with the transient buffers
 * taken from R_alloc (freed automatically when .Call returns — no manual free,
 * no RAII, nothing to leak on an error longjmp).
 *
 * Entry points mirror the generated Rcpp *exports* (they take R-native args and
 * derive extents / marshal internally), not the flat Fortran C signatures.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <string.h>

/* Flat Fortran C ABI. c_bool is one byte; unsigned char* has the same ABI as
 * _Bool* for a pointer arg, so we avoid <stdbool.h> entirely. */
extern void serialize_logical_helper_c(const unsigned char*, const int*, const int*,
                                       const int*, const char*, const int*, int*);
extern void deserialize_logical_helper_c(unsigned char*, const int*, const int*,
                                         const int*, const char*, const int*, int*);
extern void serialize_char_helper_c(const char*, const int*, const int*, const int*,
                                    const int*, const char*, const int*, int*);
extern void deserialize_char_helper_c(char*, const int*, const int*, const int*,
                                      const int*, const char*, const int*, int*);

/* ---- marshalling helpers (the pure-C tox_marshal.h) ------------------------ */

/* longest element, skipping NA — the fixed width a character(len=n) array needs */
static int tox_max_strlen(SEXP x) {
    int longest = 0, n = (int) LENGTH(x);
    for (int i = 0; i < n; ++i) {
        SEXP e = STRING_ELT(x, i);
        if (e == NA_STRING) continue;
        int m = (int) LENGTH(e);
        if (m > longest) longest = m;
    }
    return longest;
}

/* STRSXP -> char[len*n], column-major, each slot zero-padded / truncated to len */
static char *tox_char_in(SEXP x, int len) {
    int n = (int) LENGTH(x);
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char *buf = (char *) R_alloc(total, 1);
    memset(buf, 0, total);
    for (int i = 0; i < n; ++i) {
        SEXP e = STRING_ELT(x, i);
        if (e == NA_STRING) continue;
        int slen = (int) LENGTH(e);
        int m = slen < len ? slen : len;
        memcpy(buf + (size_t) i * len, CHAR(e), (size_t) m);
    }
    return buf;
}

/* char[len*n] fixed-width -> STRSXP, reading each slot up to the first '\0' */
static SEXP tox_char_out(const char *buf, int len, int n) {
    SEXP out = Rf_allocVector(STRSXP, n);   /* used immediately by caller; no alloc before return */
    for (int i = 0; i < n; ++i) {
        const char *p = buf + (size_t) i * len;
        int m = 0;
        while (m < len && p[m] != '\0') ++m;
        SET_STRING_ELT(out, i, Rf_mkCharLen(p, m));
    }
    return out;
}

/* LGLSXP (int 0/1) -> c_bool byte buffer */
static unsigned char *tox_bool_in(SEXP x) {
    int n = (int) LENGTH(x);
    unsigned char *buf = (unsigned char *) R_alloc(n > 0 ? n : 1, 1);
    const int *px = LOGICAL(x);
    for (int i = 0; i < n; ++i) buf[i] = (px[i] == TRUE) ? 1 : 0;
    return buf;
}

/* c_bool byte buffer -> LGLSXP. Read the raw byte: ifx writes 0xFF for true, so
 * test != 0 rather than widening a _Bool (which would be UB and print 255). */
static SEXP tox_bool_out(const unsigned char *buf, int n) {
    SEXP out = Rf_allocVector(LGLSXP, n);
    int *po = LOGICAL(out);
    for (int i = 0; i < n; ++i) po[i] = buf[i] != 0 ? TRUE : FALSE;
    return out;
}

/* arr_shape from a dim attribute, or [length] when there is none */
static SEXP tox_shape_of(SEXP arr) {
    SEXP dim = Rf_getAttrib(arr, R_DimSymbol);
    if (dim != R_NilValue) return dim;             /* dims are always INTSXP */
    SEXP s = Rf_allocVector(INTSXP, 1);
    INTEGER(s)[0] = (int) LENGTH(arr);
    return s;
}

static int tox_prod(SEXP shape) {
    int p = 1, n = (int) LENGTH(shape);
    const int *ps = INTEGER(shape);
    for (int i = 0; i < n; ++i) p *= ps[i];
    return p;
}

static SEXP tox_named2(const char *n0, SEXP v0, const char *n1, SEXP v1) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, v0);
    SET_VECTOR_ELT(out, 1, v1);
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar(n0));
    SET_STRING_ELT(nms, 1, Rf_mkChar(n1));
    Rf_setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(2);
    return out;
}

/* ---- .Call entry points ---------------------------------------------------- */

SEXP serialize_logical_ccall(SEXP arr, SEXP filename) {
    int nprot = 0;
    SEXP arr_shape = PROTECT(tox_shape_of(arr)); nprot++;
    int n_elements = (int) LENGTH(arr);
    int n_arr_shape_elements = (int) LENGTH(arr_shape);
    int filename_strlen = tox_max_strlen(filename);

    const unsigned char *arr_c = tox_bool_in(arr);
    const char *filename_c = tox_char_in(filename, filename_strlen);
    int ierr = 0;

    serialize_logical_helper_c(arr_c, &n_elements, INTEGER(arr_shape),
                               &n_arr_shape_elements, filename_c, &filename_strlen, &ierr);

    UNPROTECT(nprot);
    return Rf_ScalarInteger(ierr);
}

SEXP deserialize_logical_ccall(SEXP arr_shape_in, SEXP filename) {
    int nprot = 0;
    SEXP arr_shape = PROTECT(Rf_coerceVector(arr_shape_in, INTSXP)); nprot++;
    int n_elements = tox_prod(arr_shape);
    int n_arr_shape_elements = (int) LENGTH(arr_shape);
    int filename_strlen = tox_max_strlen(filename);
    const char *filename_c = tox_char_in(filename, filename_strlen);

    unsigned char *arr_c = (unsigned char *) R_alloc(n_elements > 0 ? n_elements : 1, 1);
    memset(arr_c, 0, (size_t)(n_elements > 0 ? n_elements : 1));
    int ierr = 0;

    deserialize_logical_helper_c(arr_c, &n_elements, INTEGER(arr_shape),
                                 &n_arr_shape_elements, filename_c, &filename_strlen, &ierr);

    SEXP arr = PROTECT(tox_bool_out(arr_c, n_elements)); nprot++;
    Rf_setAttrib(arr, R_DimSymbol, arr_shape);
    SEXP out = PROTECT(tox_named2("arr", arr, "ierr", Rf_ScalarInteger(ierr))); nprot++;
    UNPROTECT(nprot);
    return out;
}

SEXP serialize_char_ccall(SEXP arr, SEXP filename) {
    int nprot = 0;
    SEXP arr_shape = PROTECT(tox_shape_of(arr)); nprot++;
    int arr_strlen = tox_max_strlen(arr);
    int n_strings = (int) LENGTH(arr);
    int n_arr_shape_elements = (int) LENGTH(arr_shape);
    int filename_strlen = tox_max_strlen(filename);

    const char *arr_c = tox_char_in(arr, arr_strlen);
    const char *filename_c = tox_char_in(filename, filename_strlen);
    int ierr = 0;

    serialize_char_helper_c(arr_c, &arr_strlen, &n_strings, INTEGER(arr_shape),
                            &n_arr_shape_elements, filename_c, &filename_strlen, &ierr);

    UNPROTECT(nprot);
    return Rf_ScalarInteger(ierr);
}

SEXP deserialize_char_ccall(SEXP strlen_in, SEXP arr_shape_in, SEXP filename) {
    int nprot = 0;
    int strlen = Rf_asInteger(strlen_in);
    SEXP arr_shape = PROTECT(Rf_coerceVector(arr_shape_in, INTSXP)); nprot++;
    int n_strings = tox_prod(arr_shape);
    int n_arr_shape_elements = (int) LENGTH(arr_shape);
    int filename_strlen = tox_max_strlen(filename);
    const char *filename_c = tox_char_in(filename, filename_strlen);

    size_t total = (size_t) strlen * (n_strings > 0 ? n_strings : 1);
    char *arr_c = (char *) R_alloc(total, 1);
    memset(arr_c, 0, total);
    int ierr = 0;

    deserialize_char_helper_c(arr_c, &n_strings, &strlen, INTEGER(arr_shape),
                              &n_arr_shape_elements, filename_c, &filename_strlen, &ierr);

    SEXP arr = PROTECT(tox_char_out(arr_c, strlen, n_strings)); nprot++;
    Rf_setAttrib(arr, R_DimSymbol, arr_shape);
    SEXP out = PROTECT(tox_named2("arr", arr, "ierr", Rf_ScalarInteger(ierr))); nprot++;
    UNPROTECT(nprot);
    return out;
}

static const R_CallMethodDef CallEntries[] = {
    {"serialize_logical_ccall",   (DL_FUNC) &serialize_logical_ccall,   2},
    {"deserialize_logical_ccall", (DL_FUNC) &deserialize_logical_ccall, 2},
    {"serialize_char_ccall",      (DL_FUNC) &serialize_char_ccall,      2},
    {"deserialize_char_ccall",    (DL_FUNC) &deserialize_char_ccall,    3},
    {NULL, NULL, 0}
};

void R_init_serde_c(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
