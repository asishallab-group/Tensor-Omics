// Generated. Do not edit.
//
// Marshalling helpers for the R binding. Pure C, R's C API only. These convert what
// Fortran cannot take from R directly -- R's int-based logicals and its strings -- plus a
// couple of small shape helpers. Transient buffers come from R_alloc and are freed when the
// .Call returns, so there is nothing to free by hand and nothing to leak on an error
// longjmp. The R layer has already validated and rejected NA, so these are straight copies.
//
// Strings cross fixed-width, column-major and blank-padded, both ways: the Fortran wrapper
// takes a character(len=strlen) pointer view of the buffer rather than converting it, so
// the padding is part of the value and it has to be the padding Fortran trims.
#ifndef TOX_MARSHAL_H
#define TOX_MARSHAL_H

#include <R.h>
#include <Rinternals.h>
#include <string.h>

// weak so the one library still loads into a non-R host (Python/ctypes);
// see helper/codegen/emit/c_call.py
#pragma weak COMPLEX
#pragma weak INTEGER
#pragma weak LENGTH
#pragma weak LOGICAL
#pragma weak REAL
#pragma weak STRING_ELT
#pragma weak TYPEOF
#pragma weak XLENGTH
#pragma weak SET_STRING_ELT
#pragma weak SET_VECTOR_ELT
#pragma weak R_CHAR
#pragma weak R_alloc
#pragma weak R_DimSymbol
#pragma weak R_NamesSymbol
#pragma weak R_NilValue
#pragma weak R_NaString
#pragma weak R_registerRoutines
#pragma weak R_useDynamicSymbols
#pragma weak Rf_allocVector
#pragma weak Rf_asInteger
#pragma weak Rf_asLogical
#pragma weak Rf_asReal
#pragma weak Rf_coerceVector
#pragma weak Rf_duplicate
#pragma weak Rf_getAttrib
#pragma weak Rf_length
#pragma weak Rf_mkChar
#pragma weak Rf_mkCharLen
#pragma weak Rf_protect
#pragma weak Rf_setAttrib
#pragma weak Rf_unprotect
#pragma weak Rf_ScalarComplex
#pragma weak Rf_ScalarInteger
#pragma weak Rf_ScalarLogical
#pragma weak Rf_ScalarReal

static inline int tox_imax(int a, int b) { return a > b ? a : b; }
static inline int tox_imin(int a, int b) { return a < b ? a : b; }

// The width a character(len=n) array needs: the longest element, NA skipped. 0 for a
// non-character or absent (R_NilValue) argument, so an omitted optional reports no width.
static inline int tox_max_strlen(SEXP x) {
    if (x == R_NilValue || TYPEOF(x) != STRSXP) return 0;
    int longest = 0, n = (int) XLENGTH(x);
    for (int i = 0; i < n; ++i) {
        SEXP e = STRING_ELT(x, i);
        if (e == NA_STRING) continue;
        int m = (int) LENGTH(e);
        if (m > longest) longest = m;
    }
    return longest;
}

// Fortran carries a string's length as the leading extent: n strings of length len are
// char[len * n], column-major, each **blank**-padded. The Fortran wrapper does not convert
// this buffer, it takes a `character(len=len)` pointer view of it -- so the padding byte
// becomes part of the value, and blanks are the padding Fortran's own `trim`, `len_trim`
// and string comparison understand. A NUL would survive into the string instead. Anything
// longer than len is truncated.
static inline char* tox_char_in(SEXP x, int len) {
    int n = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char* buf = (char*) R_alloc(total, 1);
    memset(buf, ' ', total);
    for (int i = 0; i < n; ++i) {
        SEXP e = STRING_ELT(x, i);
        if (e == NA_STRING) continue;
        int slen = (int) LENGTH(e);
        int m = slen < len ? slen : len;
        memcpy(buf + (size_t) i * len, CHAR(e), (size_t) m);
    }
    return buf;
}

// A blank-filled c_char output buffer of len * n bytes. Blank and not zero, deliberately:
// the callee may leave an element untouched (an early error return, a partly filled array),
// and tox_char_out reads it back by trimming trailing blanks. A zero fill would leave NULs
// there, which Rf_mkCharLen rejects outright with an embedded-nul error.
static inline char* tox_char_alloc(int len, int n) {
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char* buf = (char*) R_alloc(total, 1);
    memset(buf, ' ', total);
    return buf;
}

// char[len * n] fixed-width -> STRSXP, each slot with its trailing blanks removed: Fortran
// blank-pads whatever it assigns into a character(len=n), and the wrapper hands that buffer
// straight through. Trailing NULs are deliberately not stripped -- nothing writes them any
// more, and Rf_mkCharLen turning a stray one into a loud R error is the right answer.
// Returned unprotected: the caller protects it straight into a result slot.
static inline SEXP tox_char_out(const char* buf, int len, int n) {
    SEXP out = Rf_allocVector(STRSXP, n);
    for (int i = 0; i < n; ++i) {
        const char* p = buf + (size_t) i * len;
        int m = len;
        while (m > 0 && p[m - 1] == ' ') --m;
        SET_STRING_ELT(out, i, Rf_mkCharLen(p, m));
    }
    return out;
}

// R logical is a 4-byte int (0/1); Fortran c_bool is one byte. Copy across.
static inline unsigned char* tox_bool_in(SEXP x) {
    int n = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    unsigned char* buf = (unsigned char*) R_alloc(n > 0 ? n : 1, 1);
    if (n > 0) {
        const int* px = LOGICAL(x);
        for (int i = 0; i < n; ++i) buf[i] = (px[i] == TRUE) ? 1 : 0;
    }
    return buf;
}

// A zero-filled c_bool output buffer of n bytes.
static inline unsigned char* tox_bool_alloc(int n) {
    unsigned char* buf = (unsigned char*) R_alloc(n > 0 ? n : 1, 1);
    memset(buf, 0, (size_t) (n > 0 ? n : 1));
    return buf;
}

// c_bool byte buffer -> LGLSXP. Read the raw byte and test != 0: ifx writes 0xFF for true,
// which must map to R's TRUE (1), not stay 255. Returned unprotected.
static inline SEXP tox_bool_out(const unsigned char* buf, int n) {
    SEXP out = Rf_allocVector(LGLSXP, n);
    int* po = LOGICAL(out);
    for (int i = 0; i < n; ++i) po[i] = buf[i] != 0 ? TRUE : FALSE;
    return out;
}

// The dim attribute of x, or a length-1 vector holding its length when it has none. A plain
// R vector carries no dim, and its shape is its length. Returned unprotected.
static inline SEXP tox_shape_of(SEXP x) {
    SEXP dim = Rf_getAttrib(x, R_DimSymbol);
    if (dim != R_NilValue) return dim;
    SEXP s = Rf_allocVector(INTSXP, 1);
    INTEGER(s)[0] = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    return s;
}

// The product of an integer shape vector: the flat element count of the array it describes.
static inline int tox_prod(SEXP shape) {
    int p = 1, n = (int) XLENGTH(shape);
    const int* ps = INTEGER(shape);
    for (int i = 0; i < n; ++i) p *= ps[i];
    return p;
}

// Count the TRUEs in a logical vector; NA has already been rejected in R.
static inline int tox_sum_true(SEXP x) {
    int n = (int) XLENGTH(x), total = 0;
    const int* px = LOGICAL(x);
    for (int i = 0; i < n; ++i) total += (px[i] == TRUE) ? 1 : 0;
    return total;
}

#endif
