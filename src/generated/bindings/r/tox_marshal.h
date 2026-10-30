// Generated. Do not edit.
//
// Marshalling helpers for the R binding. Pure C, R's C API only. These convert what
// Fortran cannot take from R directly -- R's int-based logicals, its strings, and an integer
// array that arrives as doubles -- plus a couple of small shape helpers. Transient buffers
// come from R_alloc and are freed when the .Call returns, so there is nothing to free by hand
// and nothing to leak on an error longjmp.
//
// The conversions are where the input contract's element checks live (#221): the R wrapper
// has checked each argument's kind and rank, and every element of an integer, logical or
// character array is checked here, in the pass that converts it, which costs no extra pass.
// A refusal is raised as the R wrapper raises its own: a classed condition naming the
// argument (tox_refuse).
//
// Strings cross fixed-width, column-major and blank-padded, both ways: the Fortran wrapper
// takes a character(len=strlen) pointer view of the buffer rather than converting it, so
// the padding is part of the value and it has to be the padding Fortran trims.
#ifndef TOX_MARSHAL_H
#define TOX_MARSHAL_H

#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

// weak so the one library still loads into a non-R host (Python/ctypes);
// see helper/codegen/emit/c_call.py
#pragma weak COMPLEX
#pragma weak INTEGER
#pragma weak INTEGER_NO_NA
#pragma weak LENGTH
#pragma weak LOGICAL
#pragma weak REAL
#pragma weak STRING_ELT
#pragma weak TYPEOF
#pragma weak XLENGTH
#pragma weak SET_STRING_ELT
#pragma weak SET_VECTOR_ELT
#pragma weak R_BaseEnv
#pragma weak R_CHAR
#pragma weak STRING_PTR_RO
#pragma weak R_ClassSymbol
#pragma weak R_IsNA
#pragma weak R_NaInt
#pragma weak R_NegInf
#pragma weak R_PosInf
#pragma weak R_alloc
#pragma weak R_DimSymbol
#pragma weak R_NamesSymbol
#pragma weak R_NilValue
#pragma weak R_NaString
#pragma weak R_registerRoutines
#pragma weak R_useDynamicSymbols
#pragma weak Rf_allocVector
#pragma weak Rf_asInteger
#pragma weak Rf_charIsASCII
#pragma weak Rf_asLogical
#pragma weak Rf_asReal
#pragma weak Rf_coerceVector
#pragma weak Rf_duplicate
#pragma weak Rf_eval
#pragma weak Rf_getAttrib
#pragma weak Rf_getCharCE
#pragma weak Rf_install
#pragma weak Rf_lang2
#pragma weak Rf_length
#pragma weak Rf_mkChar
#pragma weak Rf_mkCharLen
#pragma weak Rf_mkCharLenCE
#pragma weak Rf_mkString
#pragma weak Rf_protect
#pragma weak Rf_setAttrib
#pragma weak Rf_translateCharUTF8
#pragma weak Rf_type2char
#pragma weak Rf_unprotect
#pragma weak Rf_ScalarComplex
#pragma weak Rf_ScalarInteger
#pragma weak Rf_ScalarLogical
#pragma weak Rf_ScalarReal
#pragma weak Rf_ScalarString

static inline int tox_imax(int a, int b) { return a > b ? a : b; }
static inline int tox_imin(int a, int b) { return a < b ? a : b; }

// Refuse an argument: signal the condition the R wrapper's own refusals are
// (.tox_argument_error in error_handling.R) -- class c(<condition>, "tox_input_error", "tox_error", "error", "condition"),
// code NA, the argument named -- through base R's stop(), which needs nothing of the
// environment the wrappers were sourced into. It does not return. On the refusal path only.
static inline void tox_refuse(const char* argument, const char* condition, const char* format, ...) {
    char message[512];
    va_list values;
    va_start(values, format);
    vsnprintf(message, sizeof message, format, values);
    va_end(values);
    const char* fields[] = {"message", "call", "code", "name", "argument"};
    const char* classes[] = {condition, "tox_input_error", "tox_error", "error", "condition"};
    int n_fields = (int) (sizeof fields / sizeof fields[0]);
    int n_classes = (int) (sizeof classes / sizeof classes[0]);
    SEXP refusal = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SET_VECTOR_ELT(refusal, 0, Rf_mkString(message));
    SET_VECTOR_ELT(refusal, 1, R_NilValue);
    SET_VECTOR_ELT(refusal, 2, Rf_ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(refusal, 3, Rf_ScalarString(NA_STRING));
    SET_VECTOR_ELT(refusal, 4, Rf_mkString(argument));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    for (int i = 0; i < n_fields; ++i) SET_STRING_ELT(names, i, Rf_mkChar(fields[i]));
    Rf_setAttrib(refusal, R_NamesSymbol, names);
    SEXP class_ = PROTECT(Rf_allocVector(STRSXP, n_classes));
    for (int i = 0; i < n_classes; ++i) SET_STRING_ELT(class_, i, Rf_mkChar(classes[i]));
    Rf_setAttrib(refusal, R_ClassSymbol, class_);
    SEXP call = PROTECT(Rf_lang2(Rf_install("stop"), refusal));
    Rf_eval(call, R_BaseEnv);
    UNPROTECT(4);
}

static inline void tox_refuse_na(const char* argument) {
    tox_refuse(argument, "tox_na_error", "'%s' contains NA, which tensor-omics cannot represent", argument);
}

// A value as a refusal shows it: R's own spelling of NA, NaN and the infinities
static inline void tox_show(char* text, size_t size, double value) {
    if (ISNA(value)) snprintf(text, size, "NA");
    else if (ISNAN(value)) snprintf(text, size, "NaN");
    else if (value == R_PosInf) snprintf(text, size, "Inf");
    else if (value == R_NegInf) snprintf(text, size, "-Inf");
    else snprintf(text, size, "%.15g", value);
}

// An element of x that is not what `expected` says, 0-based i shown 1-based as R counts
static inline void tox_refuse_element(const char* argument, const char* expected, R_xlen_t i,
                                      double value) {
    char shown[32];
    tox_show(shown, sizeof shown, value);
    tox_refuse(argument, "tox_type_error", "'%s' must hold %s, but element %lld is %s",
               argument, expected, (long long) i + 1, shown);
}

// What the R wrapper guarantees, if something called the .Call entry directly instead
static inline void tox_refuse_type(const char* argument, const char* expected, SEXP x) {
    tox_refuse(argument, "tox_type_error", "'%s' must be %s, not an R %s",
               argument, expected, Rf_type2char(TYPEOF(x)));
}

// An integer array the caller passed, as the int32 Fortran is handed: R's integer vector as
// it stands, NA refused, or a double vector whose every element is a whole number in R's
// integer range -- converted only after the range is checked, so the cast is always defined.
// NULL for an absent (R_NilValue) optional. INTEGER_NO_NA answers at once for a sequence like
// 1:n, which R knows holds no NA.
static inline const int* tox_int_in(SEXP x, const char* argument) {
    if (x == R_NilValue) return NULL;
    R_xlen_t n = XLENGTH(x);
    if (TYPEOF(x) == INTSXP) {
        const int* values = INTEGER(x);
        if (!INTEGER_NO_NA(x))
            for (R_xlen_t i = 0; i < n; ++i)
                if (values[i] == NA_INTEGER) tox_refuse_na(argument);
        return values;
    }
    if (TYPEOF(x) != REALSXP) tox_refuse_type(argument, "an integer vector", x);
    const double* values = REAL(x);
    int* converted = (int*) R_alloc((size_t) (n > 0 ? n : 1), sizeof(int));
    for (R_xlen_t i = 0; i < n; ++i) {
        double value = values[i];
        if (!(value >= -INT_MAX && value <= INT_MAX) || value != (double) (int) value) {
            if (ISNA(value)) tox_refuse_na(argument);
            tox_refuse_element(argument, "whole numbers in the integer range", i, value);
        }
        converted[i] = (int) value;
    }
    return converted;
}

// An integer array the callee modifies, as a fresh integer vector with x's attributes:
// checked as tox_int_in checks it, then copied or converted
static inline SEXP tox_int_copy(SEXP x, const char* argument) {
    tox_int_in(x, argument);
    return TYPEOF(x) == INTSXP ? Rf_duplicate(x) : Rf_coerceVector(x, INTSXP);
}

// Whether s holds n bytes of well-formed UTF-8: no overlong form, no surrogate, nothing
// past U+10FFFF -- what R's validUTF8() accepts
static inline int tox_utf8_valid(const unsigned char* s, size_t n) {
    size_t i = 0;
    while (i < n) {
        unsigned char c = s[i];
        size_t length;
        unsigned char low = 0x80, high = 0xBF;
        if (c < 0x80) { ++i; continue; }
        else if (c >= 0xC2 && c <= 0xDF) length = 2;
        else if (c >= 0xE0 && c <= 0xEF) {
            length = 3;
            if (c == 0xE0) low = 0xA0;
            if (c == 0xED) high = 0x9F;
        } else if (c >= 0xF0 && c <= 0xF4) {
            length = 4;
            if (c == 0xF0) low = 0x90;
            if (c == 0xF4) high = 0x8F;
        } else return 0;
        if (i + length > n) return 0;
        if (s[i + 1] < low || s[i + 1] > high) return 0;
        for (size_t k = 2; k < length; ++k)
            if (s[i + k] < 0x80 || s[i + k] > 0xBF) return 0;
        i += length;
    }
    return 1;
}

// A string as the UTF-8 bytes Fortran is handed, with their count in *length -- or NULL, and
// *problem saying why it cannot cross. By what R marks it as, the same in every locale:
// - ASCII, or marked UTF-8: as it stands, if it is well-formed;
// - marked latin1: translated, as enc2utf8() translates it (R reads latin1 as Windows-1252),
//   unless it holds one of the bytes that code page leaves undefined;
// - unmarked (native): as it stands if its bytes are well-formed UTF-8, and refused otherwise.
//   Never translated from the locale: in the C locale R takes every non-ASCII byte for one it
//   cannot convert and writes it out as escape text, so the UTF-8 bytes c3 a9 of an e-acute
//   reached Fortran as the eight characters "<c3><a9>";
// - marked "bytes": it has no encoding to read it in.
// So no escape text R makes up can reach Fortran in place of the caller's string.
static inline const char* tox_utf8(SEXP e, size_t* length, const char** problem) {
    const char* text = CHAR(e);
    *length = (size_t) LENGTH(e);
    *problem = NULL;
    if (Rf_charIsASCII(e)) return text;  // and without a NUL: R cannot hold one
    cetype_t encoding = Rf_getCharCE(e);
    if (encoding == CE_BYTES) {
        *problem = "is marked as bytes, which has no encoding to read it in";
        return NULL;
    }
    if (encoding == CE_LATIN1) {
        const char* translated = Rf_translateCharUTF8(e);
        // every latin1 byte is one character, so a translation holding more characters than
        // the string has bytes is one where R wrote an undefined byte out as "<81>"
        size_t characters = 0;
        for (const unsigned char* p = (const unsigned char*) translated; *p; ++p)
            characters += (*p & 0xC0) != 0x80;
        if (characters != *length) {
            *problem = "is marked latin1, but holds a byte that is no character in it";
            return NULL;
        }
        *length = strlen(translated);
        return translated;
    }
    if (!tox_utf8_valid((const unsigned char*) text, *length)) {
        *problem = "is not valid UTF-8 (a string in another encoding needs marking as one, "
                   "such as Encoding(x) <- \"latin1\")";
        return NULL;
    }
    return text;
}

// The width a character(len=n) array needs: the longest element, NA skipped. 0 for a
// non-character or absent (R_NilValue) argument, so an omitted optional reports no width.
//
// A *present* argument reports at least 1, matching what the Python layer does with
// `.ljust(1)` and `max(..., default=0) or 1`. Every element being "" would otherwise give
// width 0, and `R_alloc(0, 1)` returns NULL -- a null buffer the Fortran wrapper cannot make
// a pointer view of. The floor is deliberately inside the present branch: raising the absent
// case to 1 would send width 1 beside a null pointer, and M_CHECK_ARRAY_NON_NULL would then
// reject an optional the caller legitimately omitted.
//
// Measured in UTF-8 bytes, as tox_char_in copies it. An NA, and a string that cannot cross
// as UTF-8, is skipped here and refused there.
static inline int tox_max_strlen(SEXP x) {
    if (x == R_NilValue || TYPEOF(x) != STRSXP) return 0;
    int longest = 0, n = (int) XLENGTH(x);
    const SEXP* elements = STRING_PTR_RO(x);
    for (int i = 0; i < n; ++i) {
        SEXP e = elements[i];
        if (e == NA_STRING) continue;
        size_t length;
        const char* problem;
        if (tox_utf8(e, &length, &problem) == NULL) continue;
        if ((int) length > longest) longest = (int) length;
    }
    return longest > 0 ? longest : 1;
}

// Fortran carries a string's length as the leading extent: n strings of length len are
// char[len * n], column-major, each **blank**-padded. The Fortran wrapper does not convert
// this buffer, it takes a `character(len=len)` pointer view of it -- so the padding byte
// becomes part of the value, and blanks are the padding Fortran's own `trim`, `len_trim`
// and string comparison understand. A NUL would survive into the string instead.
//
// Every element is checked as it is copied: an NA, a string that is not well-formed UTF-8,
// and one longer than len are refused, naming the argument. A longer one could only be cut
// short, which turned `delimiter = "::"` into ":" (len is the longest element's length
// wherever the caller's strings decide it, so this refuses only past a declared length).
static inline char* tox_char_in(SEXP x, int len, const char* argument) {
    int n = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char* buf = (char*) R_alloc(total, 1);
    memset(buf, ' ', total);
    if (n == 0) return buf;
    if (TYPEOF(x) != STRSXP) tox_refuse_type(argument, "a character vector", x);
    const SEXP* elements = STRING_PTR_RO(x);
    for (int i = 0; i < n; ++i) {
        SEXP e = elements[i];
        if (e == NA_STRING) tox_refuse_na(argument);
        size_t slen;
        const char* problem;
        const char* text = tox_utf8(e, &slen, &problem);
        if (text == NULL)
            tox_refuse(argument, "tox_type_error", "'%s' element %d %s", argument, i + 1, problem);
        if (slen > (size_t) len)
            tox_refuse(argument, "tox_type_error",
                       "'%s' element %d is %d bytes long in UTF-8, but it is declared %d long, "
                       "and a longer value would be cut short",
                       argument, i + 1, (int) slen, len);
        memcpy(buf + (size_t) i * len, text, slen);
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
//
// `out` is protected for the whole fill, unlike every other helper here: Rf_mkCharLen
// allocates, so it can trigger a GC on any iteration, and until the caller receives `out`
// nothing else references it. Protecting on return -- which is what the caller's
// PROTECT(tox_char_out(...)) does -- is too late; the window is the loop, and it grows with
// n. Balanced before returning, so the caller's own PROTECT is still the one that keeps it.
static inline SEXP tox_char_out(const char* buf, int len, int n) {
    SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; ++i) {
        const char* p = buf + (size_t) i * len;
        int m = len;
        while (m > 0 && p[m - 1] == ' ') --m;
        // UTF-8 is what crossed, both ways, so that is what the string is marked as: in a
        // non-UTF-8 locale an unmarked one would be read in the native encoding instead
        SET_STRING_ELT(out, i, Rf_mkCharLenCE(p, m, CE_UTF8));
    }
    UNPROTECT(1);
    return out;
}

// R logical is a 4-byte int (0/1); Fortran c_bool is one byte. Copied across, and so is a
// number whose every element is exactly 0 or 1 -- c(1, 0, 1) is R's own literal for a 0/1
// mask -- while an NA, and any other number, is refused by name.
static inline unsigned char* tox_bool_in(SEXP x, const char* argument) {
    R_xlen_t n = (x == R_NilValue) ? 0 : XLENGTH(x);
    unsigned char* buf = (unsigned char*) R_alloc((size_t) (n > 0 ? n : 1), 1);
    if (n == 0) return buf;
    static const char expected[] = "TRUE and FALSE, or the numbers 0 and 1";
    switch (TYPEOF(x)) {
    case LGLSXP: {
        const int* values = LOGICAL(x);
        for (R_xlen_t i = 0; i < n; ++i) {
            if (values[i] == NA_LOGICAL) tox_refuse_na(argument);
            buf[i] = values[i] != 0;
        }
        break;
    }
    case INTSXP: {
        const int* values = INTEGER(x);
        for (R_xlen_t i = 0; i < n; ++i) {
            int value = values[i];
            if (value != 0 && value != 1) {
                if (value == NA_INTEGER) tox_refuse_na(argument);
                tox_refuse_element(argument, expected, i, (double) value);
            }
            buf[i] = (unsigned char) value;
        }
        break;
    }
    case REALSXP: {
        const double* values = REAL(x);
        for (R_xlen_t i = 0; i < n; ++i) {
            double value = values[i];
            if (value != 0.0 && value != 1.0) {
                if (ISNA(value)) tox_refuse_na(argument);
                tox_refuse_element(argument, expected, i, value);
            }
            buf[i] = value != 0.0;
        }
        break;
    }
    default:
        tox_refuse_type(argument, "a logical vector", x);
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
// which must map to R's TRUE (1), not stay 255. Returned unprotected, and safely so: the one
// allocation is the first statement and nothing after it allocates, so no GC can run while
// the result is unreferenced.
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

// Count the elements of a mask that are set: the TRUEs, or the 1s of a 0/1 number vector.
// Counted before the mask is converted, so an element that is neither -- which tox_bool_in
// refuses next -- counts as unset, as R's count of it does (r_wrapper._mask_count).
static inline int tox_sum_true(SEXP x) {
    int n = (int) XLENGTH(x), total = 0;
    if (TYPEOF(x) == REALSXP) {
        const double* values = REAL(x);
        for (int i = 0; i < n; ++i) total += values[i] == 1.0;
    } else if (TYPEOF(x) == LGLSXP || TYPEOF(x) == INTSXP) {
        const int* values = TYPEOF(x) == LGLSXP ? LOGICAL(x) : INTEGER(x);
        for (int i = 0; i < n; ++i) total += values[i] == 1;
    }
    return total;
}

#endif
