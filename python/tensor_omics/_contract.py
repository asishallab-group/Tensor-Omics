"""The input contract of the Python binding: what an argument accepts, and what it refuses.

Generated. Do not edit. The rules are helper/codegen/design/language-layers.md's; R's
binding follows the same table.

Every check looks at the input's own kind first -- an allowlist on numpy's dtype kind,
before anything converts it -- and accepts a value only if converting it to the Fortran type
and back gives exactly the same value, a real's too. A refusal raises ToxArgumentError
naming the argument.
"""

import itertools
import os

import numpy as np

from .error_handling import ToxArgumentError

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

#: Every integer of at most this magnitude is a float64 exactly: its 53-bit significand holds it
_EXACT_IN_FLOAT64 = 2 ** 53


def _refuse(message, argument):
    raise ToxArgumentError(message, argument=argument) from None


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


def _where(array, index):
    """The words a refusal names the element at `index` by."""
    if array.ndim == 0:
        return "it"
    return f"element {index[0] if array.ndim == 1 else index}"


def _first(array, differs):
    """The first element `differs` marks, as (its index, the words a refusal names it by)."""
    if array.ndim == 0:
        return (), "it"
    index = tuple(int(i) for i in np.argwhere(differs)[0])
    return index, _where(array, index)


def _listed_integers(value, array, bound):
    """The integers of a list or tuple that numpy made floats of, as (index, int) pairs.

    numpy decides a list's dtype first, so one integer beside a float turns them all into
    floats: `[2**53 + 1, 0.5]` is float64 `[2**53, 0.5]` before any check sees it, and no check
    of the array can tell the rounding from a value given that way (FES 2026-10-29). These are
    what the targets' integer rules check instead, and only when the array holds a value at or
    past `bound` -- an integer within it cannot have been rounded into or out of the rule -- which
    two reductions say: a list of ordinary numbers costs no walk. Nothing for an ndarray, whose
    dtype says what it is, nor for a list numpy kept as integers.
    """
    if not isinstance(value, (list, tuple)) or array.dtype.kind != "f" or array.size == 0:
        return ()
    # strictly within: an integer past 2**53 can round to 2**53 itself. fmin and fmax pass over
    # a NaN, which no integer is
    if -bound < np.fmin.reduce(array, axis=None) and np.fmax.reduce(array, axis=None) < bound:
        return ()
    items = value
    if array.ndim > 1:  # numpy took the nesting as the shape, in C order
        for _ in range(array.ndim - 1):
            items = itertools.chain.from_iterable(items)
        items = list(items)
    if not any(issubclass(kind, (int, np.integer)) and not issubclass(kind, (bool, np.bool_))
               for kind in set(map(type, items))):
        return ()
    return [(tuple(int(i) for i in np.unravel_index(position, array.shape)), int(item))
            for position, item in enumerate(items)
            if isinstance(item, (int, np.integer)) and not isinstance(item, (bool, np.bool_))]


def _round_trip(array, argument, dtype, order, wanted):
    """`array` as `dtype`, accepted only if converting it back gives every value unchanged."""
    if array.dtype == dtype:
        return np.asarray(array, order=order)
    # a NaN, an infinity or a value out of range casts to something arbitrary, which the
    # comparison below then tells apart from what was passed
    with np.errstate(invalid="ignore", over="ignore"):
        converted = array.astype(dtype, order=order)
    if not np.array_equal(converted, array):
        index, where = _first(array, converted != array)
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
    # an integer of a list numpy made floats of, as it was given rather than as numpy rounded it
    low, high = np.iinfo(np.int32).min, np.iinfo(np.int32).max
    for index, item in _listed_integers(value, array, high):
        if not low <= item <= high:
            _refuse(f"'{argument}' must be {wanted}, but {_where(array, index)} is {item!r}",
                    argument)
    return _round_trip(array, argument, np.int32, order, wanted)


def integer(value, argument):
    """A scalar integer(int32): the 0-d case of `integers`, as a Python int."""
    return int(integers(value, argument, ndim=0))


def _refuse_inexact(argument, wanted, where, shown, rounded):
    """Refuse a value float64 would round: `shown` is it, `rounded` what float64 makes of it."""
    _refuse(f"'{argument}' must be {wanted} that a float64 holds exactly, but {where} is "
            f"{shown}, which it would round to {rounded!r}; convert it first "
            f"(`.astype(np.float64)`, `float()`) where that rounding is acceptable", argument)


def _exact_in_float64(array, converted, argument, wanted, part=None):
    """Refuse a value `converted`, its float64, does not hold exactly.

    `part` names which part of a single complex `array` is, for the refusal to say.

    Compared in `array`'s own type, converted back: compared as floats, numpy would round the
    original the same way and find nothing. A NaN or an infinity of a wider float is the same
    in float64; a value past float64's range becomes an infinity, and so differs.
    """
    if array.size == 0:
        return
    if array.dtype.kind in "iu":
        if array.min() >= -_EXACT_IN_FLOAT64 and array.max() <= _EXACT_IN_FLOAT64:
            return  # two reductions say it for the common case, without a temporary
        # A float64 at or past 2**63 (int64) or 2**64 (uint64) is past the integer's range, so
        # it was rounded up, and is decided here: casting it back is undefined, and the
        # platforms disagree (x86 wraps it, aarch64 saturates it back to the very value it came
        # from). Only what lies within the range is cast back. Nothing below 2**63 rounds down
        # out of int64's range: -2**63 is a float64 exactly.
        bits = 8 * array.dtype.itemsize - (array.dtype.kind == "i")
        past = converted >= 2.0 ** bits
        differs = past | (np.where(past, 0.0, converted).astype(array.dtype) != array)
    else:
        differs = (converted.astype(array.dtype) != array) & ~np.isnan(array)
    if differs.any():
        index, where = _first(array, differs)
        if part is not None:
            where = f"its {part} part"
        # a long double has no Python number to show it by, but its str is its decimal digits
        item = array[index].item()
        shown = repr(item) if isinstance(item, int) else str(array[index])
        _refuse_inexact(argument, wanted, where, shown, converted[index].item())


def _shown_int(item):
    """A Python int as a refusal shows it: past a few hundred digits a number says nothing, and
    str() of one past 4300 digits raises."""
    bits = item.bit_length()
    return repr(item) if bits <= 256 else f"an int of {bits} bits"


def _holds_python_integers(array):
    """Whether `array` is objects that are all Python ints: what numpy makes of one too large for
    any NumPy integer (`2**64`, `-2**63 - 1`)."""
    return array.dtype.kind == "O" and array.size > 0 and all(
        isinstance(item, int) and not isinstance(item, bool) for item in array.flat)


def _python_integers(array, argument, wanted, order):
    """Python ints too large for any NumPy integer, as float64.

    The kind would refuse them as not numbers. They are, and the rule is the same: accepted
    where float64 holds each exactly, which the int's own exact arithmetic says,
    `int(float(x)) == x`.
    """
    converted = np.empty(array.shape, dtype=np.float64, order=order)
    for index in np.ndindex(array.shape):
        item = array[index]
        try:
            rounded = float(item)
        except OverflowError:  # past float64's range
            rounded = float("inf") if item > 0 else float("-inf")
        if rounded in (float("inf"), float("-inf")) or int(rounded) != item:
            _refuse_inexact(argument, wanted, _where(array, index), _shown_int(item), rounded)
        converted[index] = rounded
    return converted


def reals(value, argument, ndim=None, order="C"):
    """A real(real64) input: integers and floats that float64 holds exactly.

    Every float up to float64 and every integer up to 32 bits widens exactly, so it is not
    compared at all, and a float64 is not even converted. A 64-bit integer is compared unless it
    lies within 2**53, and a float wider than float64 (np.longdouble) always is, as is a Python
    int too large for any NumPy integer.
    """
    array = _as_array(value, argument)
    wanted = _REALS[ndim == 0]
    if _holds_python_integers(array):
        _rank(array, argument, ndim)
        return _python_integers(array, argument, wanted, order)
    _number_kind(value, array, argument, wanted)
    _rank(array, argument, ndim)
    # an integer of a list numpy made floats of: rounded already, so checked as it was given
    for index, item in _listed_integers(value, array, _EXACT_IN_FLOAT64):
        if not -_EXACT_IN_FLOAT64 <= item <= _EXACT_IN_FLOAT64 and int(array[index]) != item:
            _refuse_inexact(argument, wanted, _where(array, index), _shown_int(item),
                            array[index].item())
    kind, size = array.dtype.kind, array.dtype.itemsize
    if not ((kind in "iu" and size > 4) or (kind == "f" and size > 8)):
        return np.asarray(array, dtype=np.float64, order=order)
    # a long double past float64's range becomes an infinity, which the comparison refuses
    with np.errstate(over="ignore"):
        converted = np.asarray(array, dtype=np.float64, order=order)
    _exact_in_float64(array, converted, argument, wanted)
    return converted


def real(value, argument):
    """A scalar real(real64): the 0-d case of `reals`, as a Python float."""
    return float(reals(value, argument, ndim=0))


def complex_number(value, argument):
    """A scalar complex(real64): a real number by `real`'s rule, or a complex one whose parts
    a float64 holds exactly; as a Python complex. A bool, a string or an object is refused."""
    array = _as_array(value, argument)
    if array.dtype.kind != "c":
        return complex(real(value, argument))
    _rank(array, argument, 0)
    if array.dtype.itemsize <= 16:  # complex64 and complex128 widen exactly
        return complex(array.item())
    # a complex long double: each part by the real rule
    with np.errstate(over="ignore"):
        converted = array.astype(np.complex128)
    _exact_in_float64(array.real, converted.real, argument, "a number", part="real")
    _exact_in_float64(array.imag, converted.imag, argument, "a number", part="imaginary")
    return complex(converted.item())


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
    as the mode it starts with. Trailing blanks are dropped first, as for every string: to
    Fortran `'ward '` is `'ward'`. A leading blank is part of the value, and refused.
    """
    listed = ", ".join(f"'{choice}'" for choice in choices)
    if not isinstance(value, str):
        _refuse(f"'{argument}' must be one of {listed}, not {type(value).__name__}", argument)
    matched = value.rstrip(" ").lower()
    if matched not in choices:
        _refuse(f"'{argument}' must be one of {listed}, not {value!r}", argument)
    return matched.encode()
