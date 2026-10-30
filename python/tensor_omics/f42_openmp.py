r"""f42_openmp

What this build of the library does with OpenMP.

OpenMP is a build option (`./build.sh --openmp`), so whether a call runs on several threads
is decided when the library is built, not by the caller. A procedure that can run on several
threads takes their number as an argument; this module answers how many a build can use.

Python binding, generated from f42_openmp. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.max_omp_threads_c.restype = None
_lib.max_omp_threads_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MAX_OMP_THREADS_ARGUMENTS = ()

def max_omp_threads(
):
    r"""The number of threads a parallel loop of this build can use: 1 without OpenMP.

    A build with OpenMP answers the OpenMP runtime's `omp_get_max_threads()`, which follows
    the environment variable `OMP_NUM_THREADS` and is otherwise the runtime's own choice,
    usually the number of processors.
    A serial build answers 1, as one thread runs every loop there.

    Returns
    -------
    max_threads : int
        The number of threads a parallel loop can use; 1 without OpenMP

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_openmp::max_omp_threads`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    max_threads = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.max_omp_threads_c(
        _ctypes.byref(max_threads),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _MAX_OMP_THREADS_ARGUMENTS)

    return max_threads.value
