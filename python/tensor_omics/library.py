"""Locating and loading the tensor-omics shared library.

Generated. Do not edit.
"""

import ctypes
import os

#: Where the library is expected, relative to the repository root
DEFAULT_LIBRARY = 'build/libtensor-omics.so'

#: Overrides the search, for an installed or relocated build
LIBRARY_ENV_VAR = "TENSOR_OMICS_LIBRARY"

_loaded = None


def nullable(argtype):
    """Allow None for an argtype, so an omitted optional can be a null pointer.

    ctypes checks an argument against its argtype and rejects None, so an
    optional needs a type that accepts it and passes NULL through.
    """
    def from_param(cls, value):
        return None if value is None else argtype.from_param(value)

    return type(
        f"nullable_{getattr(argtype, '__name__', 'argtype')}",
        (argtype,),
        {"from_param": classmethod(from_param)},
    )


def load_library():
    """Load the shared library once, and return it.

    Raises
    ------
    OSError
        If the library cannot be found, with the paths that were tried.
    """
    global _loaded
    if _loaded is not None:
        return _loaded

    override = os.environ.get(LIBRARY_ENV_VAR)
    candidates = [override] if override else []
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(here))
    candidates.append(os.path.join(root, DEFAULT_LIBRARY))

    for candidate in candidates:
        if candidate and os.path.exists(candidate):
            _loaded = ctypes.CDLL(candidate)
            return _loaded

    raise OSError(
        "cannot find the tensor-omics shared library; tried:\n  "
        + "\n  ".join(str(c) for c in candidates)
        + f"\nbuild it, or point {LIBRARY_ENV_VAR} at it"
    )
