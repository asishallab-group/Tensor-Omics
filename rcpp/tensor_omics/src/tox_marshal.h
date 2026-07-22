// Generated. Do not edit.
//
// Marshalling helpers for the R interface. C++ converts only what C cannot take from R
// directly: R's int-based logicals, and its strings. The R layer has already validated
// and rejected NA, so these are straight copies.
#ifndef TOX_MARSHAL_H
#define TOX_MARSHAL_H

#include <Rcpp.h>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace tox {

// R logical is an int (0/1/NA); Fortran c_bool is a 1-byte _Bool, which C++ bool matches.
class BoolBuffer {
    std::unique_ptr<bool[]> buf_;
    std::size_t n_;
public:
    explicit BoolBuffer(std::size_t n) : buf_(new bool[n]()), n_(n) {}
    explicit BoolBuffer(Rcpp::LogicalVector x)
        : buf_(new bool[x.size()]), n_(x.size()) {
        for (std::size_t i = 0; i < n_; ++i) buf_[i] = (x[i] == TRUE);
    }
    bool* data() { return buf_.get(); }
    Rcpp::LogicalVector to_r() const {
        // R stores a logical as an int that is exactly 0 or 1, and compares them with
        // identical(). Fortran's logical(c_bool) only promises "non-zero is true": ifx
        // writes 0xFF, which widens to 255 and prints as TRUE while comparing unequal to
        // TRUE. Test the byte rather than widening it.
        // and read the byte as unsigned char, not as bool: a C++ bool holding anything
        // other than 0 or 1 is undefined behaviour, so the compiler is entitled to fold
        // `buf_[i] ? TRUE : FALSE` down to a plain widening -- which is exactly what it
        // does, leaving the 255 in place.
        const unsigned char* raw = reinterpret_cast<const unsigned char*>(buf_.get());
        Rcpp::LogicalVector out(n_);
        for (std::size_t i = 0; i < n_; ++i) out[i] = raw[i] != 0 ? TRUE : FALSE;
        return out;
    }
};

// Fortran carries a string's length as the leading extent: a vector of n strings of
// length len is char[len * n], column-major, each string zero-padded. Reading stops at
// the first null, so an untouched buffer (zero-filled) yields empty strings, never noise.
// The width a character(len=n) array needs: the longest element. Every string is stored
// in the same n bytes, so anything shorter is null-padded and anything longer would be
// cut off.
inline int max_strlen(Rcpp::CharacterVector x) {
    int longest = 0;
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] == NA_STRING) continue;
        int n = Rf_length(STRING_ELT(x, i));
        if (n > longest) longest = n;
    }
    return longest;
}

class CharBuffer {
    std::vector<char> data_;
    int len_;
    int n_;
public:
    CharBuffer(int len, int n) : data_((std::size_t)len * n, '\0'), len_(len), n_(n) {}
    CharBuffer(Rcpp::CharacterVector x, int len)
        : data_((std::size_t)len * x.size(), '\0'), len_(len), n_(x.size()) {
        for (int i = 0; i < n_; ++i) {
            std::string s = Rcpp::as<std::string>(x[i]);
            int m = std::min<int>((int)s.size(), len_);
            std::memcpy(data_.data() + (std::size_t)i * len_, s.data(), m);
        }
    }
    char* data() { return data_.data(); }
    Rcpp::CharacterVector to_r() const {
        Rcpp::CharacterVector out(n_);
        for (int i = 0; i < n_; ++i) {
            const char* p = data_.data() + (std::size_t)i * len_;
            int m = 0;
            while (m < len_ && p[m] != '\0') ++m;
            out[i] = std::string(p, m);
        }
        return out;
    }
};

}  // namespace tox

#endif
