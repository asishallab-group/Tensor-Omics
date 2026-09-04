"""Errors raised by the tensor_omics binding.

Generated from tox_errors. Do not edit.
"""

#: How tox_errors packs the position of the offending argument into a code
ARG_POS_FACTOR = 10000

class ToxError(RuntimeError):
    """An error reported by tensor-omics.

    Attributes
    ----------
    code : int
        The tox_errors code, without the argument position packed into it.
    name : str or None
        The name of the code, e.g. 'ERR_INVALID_INPUT'.
    argument : str or None
        The argument the error was raised for, if it named one.
    """

    def __init__(self, message, code, name=None, argument=None):
        super().__init__(message)
        self.code = code
        self.name = name
        self.argument = argument


class ToxIOError(ToxError):
    """Reading or writing a file failed."""


class ToxInputError(ToxError):
    """An argument was not acceptable."""


class ToxMemoryError(ToxError):
    """Allocation failed, or a pointer was null."""


class ToxRuntimeError(ToxError):
    """The Fortran runtime reported a problem."""


class ToxInternalError(ToxError):
    """A bug in tensor-omics."""


# The codes themselves, so a caller can compare against a name
ERR_OK = 0
#: no error, operation successful
ERR_FILE_OPEN = 101
#: could not open file
ERR_READ_MAGIC = 102
#: could not read magic number
ERR_READ_TYPE = 103
#: could not read array type error
ERR_READ_NDIMS = 104
#: could not read number of dimensions
ERR_READ_DIMS = 105
#: could not read array dimensions
ERR_READ_CHARLEN = 106
#: could not read character length
ERR_READ_DATA = 107
#: could not read array data
ERR_WRITE_MAGIC = 112
#: could not write magic number
ERR_WRITE_TYPE = 113
#: could not write array type error
ERR_WRITE_NDIMS = 114
#: could not write number of dimensions
ERR_WRITE_DIMS = 115
#: could not write array dimensions
ERR_WRITE_CHARLEN = 116
#: could not write character length
ERR_WRITE_DATA = 117
#: could not write array data
ERR_FILE_ADD = 121
#: Could not add file to archive
ERR_FILE_EXTRACT = 122
#: Could not extract file from archive
ERR_FILE_CLOSE = 123
#: Failed to close the file
ERR_INVALID_FORMAT = 200
#: invalid format detected
ERR_INVALID_INPUT = 201
#: invalid input arguments
ERR_EMPTY_INPUT = 202
#: empty input arrays
ERR_DIM_MISMATCH = 203
#: dimensions do not match expected shape
ERR_NAN_INF = 204
#: NaN or Inf found where not allowed
ERR_UNSUPPORTED_TYPE = 205
#: unsupported data type encountered
ERR_SIZE_MISMATCH = 206
#: Array size mismatch
ERR_TYPE_MISMATCH = 207
#: Array type read does not match expected type
ERR_STRING_TOO_LONG = 208
#: String exceeds buffer size
ERR_IDX_OUT_OF_BOUNDS = 209
#: Array index out of bounds
ERR_DIVISION_BY_ZERO = 210
#: Division by zero encountered
ERR_ALLOC_FAIL = 301
#: memory allocation failed
ERR_POINTER_NULL = 302
#: null pointer dereference
ERR_UNIT_NOT_CONNECTED = 5002
#: Fortran runtime error: unit not connected
ERR_INTERNAL = 9001
#: unexpected internal state or logic error
ERR_UNKNOWN = 9999
#: unknown error

_MESSAGES = {
    101: "could not open file",
    102: "could not read magic number",
    103: "could not read array type error",
    104: "could not read number of dimensions",
    105: "could not read array dimensions",
    106: "could not read character length",
    107: "could not read array data",
    112: "could not write magic number",
    113: "could not write array type error",
    114: "could not write number of dimensions",
    115: "could not write array dimensions",
    116: "could not write character length",
    117: "could not write array data",
    121: "Could not add file to archive",
    122: "Could not extract file from archive",
    123: "Failed to close the file",
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
    9999: "unknown error",
}

_NAMES = {
    101: "ERR_FILE_OPEN",
    102: "ERR_READ_MAGIC",
    103: "ERR_READ_TYPE",
    104: "ERR_READ_NDIMS",
    105: "ERR_READ_DIMS",
    106: "ERR_READ_CHARLEN",
    107: "ERR_READ_DATA",
    112: "ERR_WRITE_MAGIC",
    113: "ERR_WRITE_TYPE",
    114: "ERR_WRITE_NDIMS",
    115: "ERR_WRITE_DIMS",
    116: "ERR_WRITE_CHARLEN",
    117: "ERR_WRITE_DATA",
    121: "ERR_FILE_ADD",
    122: "ERR_FILE_EXTRACT",
    123: "ERR_FILE_CLOSE",
    200: "ERR_INVALID_FORMAT",
    201: "ERR_INVALID_INPUT",
    202: "ERR_EMPTY_INPUT",
    203: "ERR_DIM_MISMATCH",
    204: "ERR_NAN_INF",
    205: "ERR_UNSUPPORTED_TYPE",
    206: "ERR_SIZE_MISMATCH",
    207: "ERR_TYPE_MISMATCH",
    208: "ERR_STRING_TOO_LONG",
    209: "ERR_IDX_OUT_OF_BOUNDS",
    210: "ERR_DIVISION_BY_ZERO",
    301: "ERR_ALLOC_FAIL",
    302: "ERR_POINTER_NULL",
    5002: "ERR_UNIT_NOT_CONNECTED",
    9001: "ERR_INTERNAL",
    9999: "ERR_UNKNOWN",
}

_EXCEPTIONS = {
    101: ToxIOError,
    102: ToxIOError,
    103: ToxIOError,
    104: ToxIOError,
    105: ToxIOError,
    106: ToxIOError,
    107: ToxIOError,
    112: ToxIOError,
    113: ToxIOError,
    114: ToxIOError,
    115: ToxIOError,
    116: ToxIOError,
    117: ToxIOError,
    121: ToxIOError,
    122: ToxIOError,
    123: ToxIOError,
    200: ToxInputError,
    201: ToxInputError,
    202: ToxInputError,
    203: ToxInputError,
    204: ToxInputError,
    205: ToxInputError,
    206: ToxInputError,
    207: ToxInputError,
    208: ToxInputError,
    209: ToxInputError,
    210: ToxInputError,
    301: ToxMemoryError,
    302: ToxMemoryError,
    5002: ToxRuntimeError,
    9001: ToxInternalError,
    9999: ToxInternalError,
}

#: Status codes are outcomes, not failures, and never raise
_STATUSES = {
}

def check_err_code(ierr, arguments=(), sources=()):
    """Raise if `ierr` reports an error.

    Parameters
    ----------
    ierr : int
        The code as tox_errors encoded it, argument position included.
    arguments : sequence of str, optional
        The names of the wrapped procedure's arguments, in declaration
        order, so the message can name the offending one rather than
        give its number.
    sources : sequence of str or None, optional
        Positionally alongside `arguments`: the argument the caller actually
        passed, where the Fortran one was derived from it. An extent read off
        an array names that array here, so the message can blame something the
        caller wrote. None where the Fortran argument is the caller's own.

    Returns
    -------
    str or None
        The name of the status code, if the call reported one. Status
        codes are outcomes rather than failures and never raise.

    Raises
    ------
    ToxError
        If `ierr` reports an error. The subclass follows the kind of
        error; every one of them is a ToxError.
    """
    code = ierr % ARG_POS_FACTOR
    if code == 0:
        return None

    if code in _STATUSES:
        return _STATUSES[code]

    arg_pos = ierr // ARG_POS_FACTOR
    message = _MESSAGES.get(code, f"unmapped error code {code}")
    argument = None
    if arg_pos > 0:
        if arg_pos <= len(arguments):
            argument = arguments[arg_pos - 1]
            source = sources[arg_pos - 1] if arg_pos <= len(sources) else None
            if source:
                # the caller's word first, then the Fortran argument it
                # was derived from, which is what the signature calls it
                message = f"{message} (argument '{source}', via '{argument}')"
                argument = source
            else:
                message = f"{message} (argument '{argument}')"
        else:
            message = f"{message} (argument {arg_pos})"

    raise _EXCEPTIONS.get(code, ToxError)(
        message, code=code, name=_NAMES.get(code), argument=argument
    )
