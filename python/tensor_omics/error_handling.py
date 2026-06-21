CODES = {
    101: "could not open file",
    102: "could not read magic number",
    103: "could not read array type code",
    104: "could not read number of dimensions",
    105: "could not read array dimensions",
    106: "could not read character length",
    107: "could not read array data",
    112: "could not write magic number",
    113: "could not write array type code",
    114: "could not write number of dimensions",
    115: "could not write array dimensions",
    116: "could not write character length",
    117: "could not write array data",
    121: "Could not add file to archive",
    122: "Could not extract file from archive",
    123: "Manifest in zip file is missing",
    124: "Failed to close the file",
    200: "invalid format detected",
    201: "invalid input arguments",
    202: "empty input arrays",
    203: "dimensions do not match expected shape",
    204: "NaN or Inf found where not allowed",
    205: "unsupported data type encountered",
    206: "Array size mismatch",
    207: "Array type read does not match expected type",
    208: "String exceeds buffer size",
    209: "Array index out of bounds",
    210: "Division by zero encountered",
    301: "memory allocation failed",
    302: "null pointer dereference",
    5002: "Fortran runtime error: unit not connected",
    9001: "unexpected internal state or logic error",
    9999: "unknown error"
}

def check_err_code(ierr: int):
    if ierr == 0:
        return

    msg = CODES.get(ierr, f"Unmapped error code: {ierr}")
    raise RuntimeError(msg)
