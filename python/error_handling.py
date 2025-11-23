def check_err_code(ierr: int) -> None:
    if ierr == 0:
        return
    msg = {
        # Success
        0: "No error, operation successful.",
        
        # I/O errors
        101: "Could not open file.",
        102: "Could not read magic number.",
        103: "Could not read type code.",
        104: "Could not read number of dimensions.",
        105: "Could not read array dimensions.",
        106: "Could not read character length.",
        107: "Could not read array data.",
        112: "Could not write magic number",
        113: "Could not write type code",
        114: "Could not write number of dimensions",
        115: "Could not write dimensions",
        116: "Could not write character length",
        117: "Could not write array data",
        121: "Could not add file to archive.",
        122: "Could not extract file from archive.",
        123: "Manifest in zip file is missing.",
        124: "Failed to close the file.",

        # FORMAT ERRORS
        200: "Invalid format detected.",
        201: "Invalid input provided.",
        202: "Empty input arrays provided.",
        203: "Dimension mismatch detected.",
        204: "NaN or Inf found in input data.",
        205: "Unsupported data type encountered.",
        206: "Array size mismatch detected",
        207: "String exceeds buffer size.",
        208: "Array index out of bounds",
        209: "Division by Zero",

        # MEMORY ERRORS
        301: "Memory allocation failed.",
        302: "Null pointer reference encountered.",

        # FORTRAN RUNTIME ERRORS
        5002: "Fortran runtime error: unit not open / not connected.",

        # Internal errors
        9001: "Internal error: unexpected state.",
        9999: "Unknown error.",
    }.get(ierr, f"Unmapped error code: {ierr}")
    raise RuntimeError(msg)
