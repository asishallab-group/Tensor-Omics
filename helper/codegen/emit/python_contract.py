"""The Python half of the input contract: the package module every generated function checks by.

What an argument accepts is one table for both binding languages (`design/language-layers.md`,
"Decision: what an argument accepts"). This is Python's side of it, written once rather than
into every generated function, as R's `tox_validate.R` is: a generated function calls
`_contract.integers(x, "x", ndim=1)` and gets back the array C is handed, or a
`ToxArgumentError` naming the argument.

The engine is the same for every target that is not a real: look at the input's own kind
first -- an allowlist on numpy's dtype kind, before anything converts it, because numpy's own
conversion parses strings, counts days since 1970 and turns None into False -- and then accept
the value only if converting it to the Fortran type and back gives exactly the same value. One
comparison catches a fraction, NaN, an infinity and a value out of range, and it costs nothing
where the input already has the Fortran type: there is no pass at all then.
"""

from __future__ import annotations

from .errors_python import ARGUMENT_EXCEPTION


def contract_module() -> str:
    """The text of `python/tensor_omics/_contract.py`."""
    return _CONTRACT.replace("__ARGUMENT_EXCEPTION__", ARGUMENT_EXCEPTION)


_CONTRACT = r'''"""The input contract of the Python binding: what an argument accepts, and what it refuses.

Generated. Do not edit. The rules are helper/codegen/design/language-layers.md's; R's
binding follows the same table.

Every check looks at the input's own kind first -- an allowlist on numpy's dtype kind,
before anything converts it -- and accepts a value only if converting it to the Fortran type
and back gives exactly the same value. A refusal raises __ARGUMENT_EXCEPTION__ naming the
argument.
"""

import os

import numpy as np

from .error_handling import __ARGUMENT_EXCEPTION__

#: The dtype kinds a number may arrive as: signed, unsigned and floating. A bool is not a
#: number here, and nothing else -- a string, an object, a datetime, a complex -- is turned
#: into one.
_NUMBER_KINDS = "iuf"

#: What each target accepts, as a refusal states it: (for an array, for a single value)
_INTEGERS = ("an array of whole numbers in the range of np.int32",
             "a whole number in the range of np.int32")
_REALS = ("an array of numbers", "a number")
_LOGICALS = ("an array of booleans, or of the numbers 0 and 1",
             "a boolean, or the number 0 or 1")


def _refuse(message, argument):
    raise __ARGUMENT_EXCEPTION__(message, argument=argument) from None


def _as_array(value, argument):
    """`value` as an ndarray, unconverted; a masked array only if nothing in it is masked."""
    # np.asarray reads a masked array's data and drops its mask, so a masked element would
    # reach the library as whatever value happens to sit underneath it
    if np.ma.isMaskedArray(value) and np.ma.is_masked(value):
        _refuse(f"'{argument}' is a masked array with masked elements, which hold no value "
                f"the library could use; fill them first (`.filled(...)`)", argument)
    try:
        return np.asarray(value)
    except (TypeError, ValueError) as error:
        _refuse(f"'{argument}' must be array-like: {error}", argument)


def _rank(array, argument, ndim):
    """Refuse any number of dimensions but `ndim`; None takes any (a serialized array)."""
    if ndim is None or array.ndim == ndim:
        return
    if ndim == 0:
        _refuse(f"'{argument}' must be a single value, not an array of shape {array.shape}",
                argument)
    _refuse(f"'{argument}' must have {ndim} dimension{'' if ndim == 1 else 's'}, "
            f"but has {array.ndim}", argument)


def _described(value, array):
    """What the caller passed, as a refusal names it."""
    if isinstance(value, np.ndarray) or array.ndim > 0:
        return f"an array of dtype {array.dtype}"
    return f"{value!r} ({type(value).__name__})"


def _number_kind(value, array, argument, wanted, logical=False):
    """Refuse an input whose kind is not a number's (or a bool's, for a logical)."""
    kind = array.dtype.kind
    if kind in _NUMBER_KINDS or (logical and kind == "b"):
        return
    hint = ""
    if kind == "b":
        hint = "; a bool is not a number here (`.astype(int)` makes it one)"
    elif kind == "S":
        hint = "; those are bytes (`.astype(str)` decodes them, but no string is a number)"
    _refuse(f"'{argument}' must be {wanted}, not {_described(value, array)}{hint}", argument)


def _round_trip(array, argument, dtype, order, wanted):
    """`array` as `dtype`, accepted only if converting it back gives every value unchanged."""
    if array.dtype == dtype:
        return np.asarray(array, order=order)
    # a NaN, an infinity or a value out of range casts to something arbitrary, which the
    # comparison below then tells apart from what was passed
    with np.errstate(invalid="ignore", over="ignore"):
        converted = array.astype(dtype, order=order)
    if not np.array_equal(converted, array):
        if array.ndim == 0:
            index, where = (), "it"
        else:
            index = tuple(int(i) for i in np.argwhere(converted != array)[0])
            where = f"element {index[0] if array.ndim == 1 else index}"
        item = array[index].item()
        hint = ""
        if array.dtype == np.uint32 and dtype == np.int32:
            hint = "; a bit-mask word of 2**31 or more passes as `.view(np.int32)`"
        _refuse(f"'{argument}' must be {wanted}, but {where} is {item!r}{hint}", argument)
    return converted


def integers(value, argument, ndim=None, order="C"):
    """An integer(int32) input: integers, or floats whose every value is whole and in range."""
    array = _as_array(value, argument)
    wanted = _INTEGERS[ndim == 0]
    _number_kind(value, array, argument, wanted)
    _rank(array, argument, ndim)
    return _round_trip(array, argument, np.int32, order, wanted)


def integer(value, argument):
    """A scalar integer(int32): the 0-d case of `integers`, as a Python int."""
    return int(integers(value, argument, ndim=0))


def reals(value, argument, ndim=None, order="C"):
    """A real(real64) input: any integers or floats, as the nearest float64."""
    array = _as_array(value, argument)
    _number_kind(value, array, argument, _REALS[ndim == 0])
    _rank(array, argument, ndim)
    return np.asarray(array, dtype=np.float64, order=order)


def real(value, argument):
    """A scalar real(real64): the 0-d case of `reals`, as a Python float."""
    return float(reals(value, argument, ndim=0))


def logicals(value, argument, ndim=None, order="C"):
    """A logical input: booleans, or numbers whose every value is exactly 0 or 1."""
    array = _as_array(value, argument)
    wanted = _LOGICALS[ndim == 0]
    _number_kind(value, array, argument, wanted, logical=True)
    _rank(array, argument, ndim)
    # an integer is 0 or 1 exactly when it lies between them, which two reductions say without
    # the temporary a comparison with the bools would make; a refusal goes on to say which one
    if array.dtype.kind in "iu" and (array.size == 0 or (array.min() >= 0 and array.max() <= 1)):
        return array.astype(np.bool_, order=order)
    return _round_trip(array, argument, np.bool_, order, wanted)


def logical(value, argument):
    """A scalar logical: the 0-d case of `logicals`, as a Python bool."""
    return bool(logicals(value, argument, ndim=0))


def _encoded(text, argument, element, width):
    """`text` as the UTF-8 bytes Fortran is handed, at most `width` of them if that is fixed.

    Trailing blanks are passed on as they are. To Fortran they are not part of the value --
    `'a '` and `'a'` compare equal, and the padding adds more -- which the contract documents
    rather than refuses: Fortran semantics, like 1-based indexing.
    """
    where = "" if element is None else f" (element {element})"
    if "\x00" in text:
        _refuse(f"'{argument}' contains a NUL character{where}, which no Fortran string "
                f"can carry", argument)
    try:
        encoded = text.encode("utf-8")
    except UnicodeEncodeError as error:
        _refuse(f"'{argument}' is not valid Unicode{where}: {error.reason} at position "
                f"{error.start}", argument)
    if width is not None and len(encoded) > width:
        _refuse(f"'{argument}' is {len(encoded)} bytes long in UTF-8{where}, but it is "
                f"declared {width} long, and a longer value would be cut short", argument)
    return encoded


def string(value, argument, width=None):
    """A string: a str, or an os.PathLike whose path is a str; as UTF-8 bytes."""
    if isinstance(value, os.PathLike):
        value = os.fspath(value)
    if not isinstance(value, str):
        hint = "; decode bytes first" if isinstance(value, (bytes, bytearray)) else ""
        _refuse(f"'{argument}' must be a str, not {type(value).__name__}{hint}", argument)
    return _encoded(value, argument, None, width)


def strings(value, argument, width=None, any_rank=False):
    """An array of strings: a list, tuple, ndarray or Series of str; as UTF-8 bytes each.

    `any_rank` is an array whose shape travels separately (a serialized one), which is taken
    in Fortran order and may be a single str: an array of rank 0.
    """
    if any_rank:
        items = _as_array(value, argument).ravel(order="F")
    else:
        if isinstance(value, (str, bytes, bytearray)):
            _refuse(f"'{argument}' must be a sequence of str, not a single "
                    f"{type(value).__name__}; wrap it in a list", argument)
        if isinstance(value, (set, frozenset, dict)):
            # iterable, but in no order the caller chose: element i would mean nothing
            _refuse(f"'{argument}' must be a sequence of str, not a {type(value).__name__}, "
                    f"which has no order", argument)
        if isinstance(value, np.ndarray):
            _rank(value, argument, 1)
        try:
            items = list(value)
        except TypeError:
            _refuse(f"'{argument}' must be a sequence of str, not {type(value).__name__}",
                    argument)
    # At C speed while nothing is wrong: str.encode refuses a non-str (np.str_ is a str) and a
    # lone surrogate itself, and one join finds a NUL. Only a refusal walks the elements one by
    # one, to say which is at fault.
    try:
        encoded = list(map(str.encode, items))
    except (TypeError, UnicodeEncodeError):
        encoded = None
    if (encoded is None or b"\x00" in b"".join(encoded)
            or (width is not None and max(map(len, encoded), default=0) > width)):
        for element, item in enumerate(items):
            if not isinstance(item, str):
                hint = ""
                if isinstance(item, (bytes, bytearray)):
                    hint = "; decode them first (`.astype(str)` for a numpy array)"
                _refuse(f"'{argument}' must hold only str, but element {element} is "
                        f"{type(item).__name__}{hint}", argument)
            _encoded(item, argument, element, width)
    return encoded


def mode(value, argument, choices):
    """A mode: a str that is one of `choices`, case aside; as the lower-case UTF-8 bytes.

    Matched whole: a value is never cut to the width C is handed, so a longer one cannot pass
    as the mode it starts with.
    """
    listed = ", ".join(f"'{choice}'" for choice in choices)
    if not isinstance(value, str):
        _refuse(f"'{argument}' must be one of {listed}, not {type(value).__name__}", argument)
    if value.lower() not in choices:
        _refuse(f"'{argument}' must be one of {listed}, not {value!r}", argument)
    return value.lower().encode()
'''
