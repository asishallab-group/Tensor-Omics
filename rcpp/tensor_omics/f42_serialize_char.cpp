#include <Rcpp.h>

using namespace Rcpp;

static_assert(sizeof(Rcomplex) == sizeof(double _Complex) && alignof(Rcomplex) == alignof(double _Complex), "Rcomplex layout is incompatible with C's 'double _Complex' and thus Fortran's 'c_double_complex'");

extern "C" {
    void serialize_char_1d_c(
        const char* arr,
        const int* arr_strlen,
        const int* n_arr_elements,
        const char* filename,
        const int* filename_strlen,
        int* ierr
    );
}

List serialize_char_1d_rcpp(
    const StringVector& arr,
    const String& filename
) {


    IntegerVector arr_shape = arr.attr("dim");
    n_arr_elements = arr_shape[0]

    SEXP* arr_SEXP_p = arr.begin();
    int arr_strlen = 0;
    for (int i = 0; i < arr.size(); ++i)
        arr_strlen = std::max(arr_strlen, (int)std::strlen(CHAR(arr_p[i])));

    std::vector<char> arr_c(arr.size());
    const char* arr_p = arr_c.data();
    // Check for NA values and return ERR_NAN_INF code, and convert
    for (int i = 0; i < arr.size(); ++i) {
        if (arr_SEXP_p[i] == NA_STRING) {
            return List::create(Named("ierr") = 204);
        } else {
            const char* str = CHAR(arr_SEXP_p[i]);
            int len = std::min((int)std::strlen(str), arr_strlen);
            std::memcpy(arr_p + i * arr_strlen, str, len);
        }
    }
    if (filename == NA_STRING) {
        return List::create(Named("ierr") = 204);
    }
    const char* filename_p = filename.get_cstring();
    int filename_strlen = std::strlen(filename_p);

    int ierr = 0;

    serialize_char_1d_c(
        arr_p,
        &arr_strlen,
        &n_arr_elements,
        filename_p,
        &filename_strlen,
        &ierr
    );

//{format(self.arguments, "type_conversion_outputs")}

    return List::create(
        Named("ierr") = ierr
    );
}