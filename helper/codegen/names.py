"""The names the generated code introduces, and the words a published name may not be.

A dummy argument's name is published verbatim: it is a Python parameter, an R parameter, a C
parameter of the R `.Call` shim, and a dummy of every generated Fortran wrapper. A procedure's
name is published the same way, as a Python and an R function. So each has to be a name every
one of those languages accepts there (`language_word`, which `ir.validate` refuses), and an
argument must not collide with a name the generated code introduces itself.

That collision is designed out wherever a language allows it. Every name the R shim's C and a
generated Python module bind beside the arguments starts with `_`, which no Fortran name can,
so nothing an author writes can meet one. A Fortran wrapper cannot do that -- a Fortran name
cannot start with `_` either -- so its locals keep the `<argument><suffix>` form built here, and
`emit.name_clashes` refuses an argument only where the wrapper actually declares, calls, imports
or uses as a kind the same name. The emitters build every name from this module, and
`tests/test_names.py` reads the generated code back to hold them to it.

What is deliberately *not* refused:

- Fortran keywords. Fortran has no reserved words, and main publishes `value`.
- R function names. An argument cannot shadow a function call in R -- R skips non-function
  bindings when it looks a call up -- and main publishes `mode`, `vector`, `array`,
  `expression`, `factor`, `points` and `weights`.
- Python soft keywords (`match`, `case`, `type`, `_`), which are valid parameter names. `type`
  is a builtin, though, so it is refused as one.
"""

from __future__ import annotations

import keyword

from .config import CONVENTIONS

# -- the words of the binding languages -------------------------------------------------------

#: Python's hard keywords. A parameter named like one is a SyntaxError, and since the package
#: `__init__` imports every module, `import tensor_omics` fails as a whole.
PYTHON_KEYWORDS = frozenset(keyword.kwlist)

#: Python's builtins (`len`, `max`, `str`, `input`, ...), as Python 3.14 has them with `site`
#: loaded. The generated functions call several of them, and a parameter of the same name hides
#: the builtin for the whole function body; a function of the same name hides it for the module.
#: Written out rather than read from the running interpreter, so the same source is accepted or
#: refused alike on every machine (`python -S` alone drops six of these).
PYTHON_BUILTINS = frozenset({
    "ArithmeticError", "AssertionError", "AttributeError", "BaseException",
    "BaseExceptionGroup", "BlockingIOError", "BrokenPipeError", "BufferError", "BytesWarning",
    "ChildProcessError", "ConnectionAbortedError", "ConnectionError", "ConnectionRefusedError",
    "ConnectionResetError", "DeprecationWarning", "EOFError", "Ellipsis", "EncodingWarning",
    "EnvironmentError", "Exception", "ExceptionGroup", "False", "FileExistsError",
    "FileNotFoundError", "FloatingPointError", "FutureWarning", "GeneratorExit", "IOError",
    "ImportError", "ImportWarning", "IndentationError", "IndexError", "InterruptedError",
    "IsADirectoryError", "KeyError", "KeyboardInterrupt", "LookupError", "MemoryError",
    "ModuleNotFoundError", "NameError", "None", "NotADirectoryError", "NotImplemented",
    "NotImplementedError", "OSError", "OverflowError", "PendingDeprecationWarning",
    "PermissionError", "ProcessLookupError", "PythonFinalizationError", "RecursionError",
    "ReferenceError", "ResourceWarning", "RuntimeError", "RuntimeWarning", "StopAsyncIteration",
    "StopIteration", "SyntaxError", "SyntaxWarning", "SystemError", "SystemExit", "TabError",
    "TimeoutError", "True", "TypeError", "UnboundLocalError", "UnicodeDecodeError",
    "UnicodeEncodeError", "UnicodeError", "UnicodeTranslateError", "UnicodeWarning",
    "UserWarning", "ValueError", "Warning", "ZeroDivisionError",
    "abs", "aiter", "all", "anext", "any", "ascii", "bin", "bool", "breakpoint", "bytearray",
    "bytes", "callable", "chr", "classmethod", "compile", "complex", "copyright", "credits",
    "delattr", "dict", "dir", "divmod", "enumerate", "eval", "exec", "exit", "filter", "float",
    "format", "frozenset", "getattr", "globals", "hasattr", "hash", "help", "hex", "id", "input",
    "int", "isinstance", "issubclass", "iter", "len", "license", "list", "locals", "map", "max",
    "memoryview", "min", "next", "object", "oct", "open", "ord", "pow", "print", "property",
    "quit", "range", "repr", "reversed", "round", "set", "setattr", "slice", "sorted",
    "staticmethod", "str", "sum", "super", "tuple", "type", "vars", "zip",
})

#: R's reserved words (`?Reserved`). A parameter or a function named like one is a parse error,
#: which makes `r/load_tensor_omics.R` fail for every wrapper, not only the one that has it.
R_RESERVED = frozenset({
    "if", "else", "repeat", "while", "function", "for", "next", "break", "in",
    "TRUE", "FALSE", "NULL", "Inf", "NaN",
    "NA", "NA_integer_", "NA_real_", "NA_complex_", "NA_character_",
})

#: C's keywords as R's compiler takes them (`gcc -std=gnu2x`: C23, plus GNU's `asm`). Every
#: published argument is a parameter of the R `.Call` shim, and the shims are compiled into the
#: one library, so a C keyword fails the library build. The `_Upper` keywords are listed for
#: completeness; no Fortran name can start with `_`.
C_KEYWORDS = frozenset({
    "auto", "break", "case", "char", "const", "continue", "default", "do", "double",
    "else", "enum", "extern", "float", "for", "goto", "if", "inline", "int", "long",
    "register", "restrict", "return", "short", "signed", "sizeof", "static", "struct",
    "switch", "typedef", "union", "unsigned", "void", "volatile", "while",
    "alignas", "alignof", "bool", "constexpr", "false", "nullptr", "static_assert",
    "thread_local", "true", "typeof", "typeof_unqual",
    "asm",
    "_Alignas", "_Alignof", "_Atomic", "_BitInt", "_Bool", "_Complex", "_Decimal128",
    "_Decimal32", "_Decimal64", "_Generic", "_Imaginary", "_Noreturn", "_Static_assert",
    "_Thread_local",
})


def language_word(name: str, languages=("Python", "R", "C")) -> str | None:
    """What makes `name` unusable where it is published, or None.

    Case-sensitive, as Python, R and C are: the published name keeps its declared spelling, so
    `NA` is R's reserved word while `na` is a name. `languages` narrows it to where the name is
    published: an argument reaches all three, a procedure is a function in Python and R only
    (its C symbols carry a suffix).
    """
    if "Python" in languages and name in PYTHON_KEYWORDS:
        return "a Python keyword"
    if "R" in languages and name in R_RESERVED:
        return "an R reserved word"
    if "C" in languages and name in C_KEYWORDS:
        return "a C keyword"
    if "Python" in languages and name in PYTHON_BUILTINS:
        return "a Python builtin"
    return None


# -- what a generated Python module binds: all of it underscored -------------------------------

#: The standard-library modules a generated Python module imports, as (module, bound name)
PYTHON_STDLIB_IMPORTS = (("ctypes", "_ctypes"), ("os", "_os"))
#: The third-party ones, in a group of their own
PYTHON_THIRD_PARTY_IMPORTS = (("numpy", "_np"),)
#: What it imports from the package itself, as (module, ((name, bound name), ...))
PYTHON_PACKAGE_IMPORTS = (
    (".error_handling", (("ToxArgumentError", "_ToxArgumentError"),
                         ("check_err_code", "_check_err_code"))),
    (".library", (("load_library", "_load_library"), ("nullable", "_nullable"))),
)
#: The names those are bound to, as the generated code spells them
PYTHON_NUMPY = "_np"
PYTHON_CTYPES = "_ctypes"
PYTHON_CHECK_ERR_CODE = "_check_err_code"
#: What a generated function raises for an argument it refuses itself (`emit.errors_python`)
PYTHON_ARGUMENT_ERROR = "_ToxArgumentError"
PYTHON_LOAD_LIBRARY = "_load_library"
PYTHON_NULLABLE = "_nullable"
#: The exception a generated Python function catches while converting an input
PYTHON_CAUGHT_ERROR = "_error"

#: Every name a generated Python module binds beside its wrappers and their parameters
PYTHON_BOUND = frozenset(
    [bound for _, bound in PYTHON_STDLIB_IMPORTS + PYTHON_THIRD_PARTY_IMPORTS]
    + [bound for _, names in PYTHON_PACKAGE_IMPORTS for _, bound in names]
    + [PYTHON_CAUGHT_ERROR]
)

# -- what the R shim declares: all of it underscored --------------------------------------------
#
# A shim local is `_<tag>_<argument>`: the leading underscore no Fortran name can have, then a
# lower-case tag, so the identifier is legal C at block scope whatever case the argument has
# (`_` before an upper-case letter is reserved everywhere; before a lower-case one only at file
# scope, where the shim declares nothing).

#: A by-value scalar input, pulled out of its length-1 vector
SCALAR_VALUE = "v"
#: An optional, as the pointer C is handed (null when omitted)
OPTIONAL_POINTER = "p"
#: An optional's length (0 when omitted)
OPTIONAL_SIZE = "size"
#: An optional array's extents (all 0 when omitted); also the `dim` attribute of a matrix output
DIM = "dim"
#: An optional array's `dim` attribute, as read
OPTIONAL_DIM_ATTRIBUTE = "d"
#: The C buffer a logical or a string is converted into
BUFFER = "c"
#: The copy of an array modified in place, so the caller's stays intact
INOUT_COPY = "out"
#: The count of protected objects, a local of every shim function
C_PROTECT_COUNT = "_nprot"
#: The index of the loops that copy an optional's extents
C_LOOP_INDEX = "_i"


def shim_local(argument: str, tag: str) -> str:
    """The R shim's local for `argument`: `_<tag>_<argument>`."""
    return f"_{tag}_{argument}"


# -- what a Fortran wrapper declares: `<argument><suffix>`, refused where emitted -------------

#: Fortran C wrapper: the local a converted argument is passed to the callee as
CONVERTED = "_f"
#: Fortran C wrapper: the integer parameter a mode string names
MODE_VALUE = "_mode_f"
#: Fortran C wrapper: a nullable mode as the callee is handed it -- a pointer to the
#: `MODE_VALUE` local, left disassociated, and so absent, when C passed null
MODE_POINTER = "_mode_p"
#: Fortran C wrapper: the address an array of strings is viewed at -- the caller's buffer, or
#: `EMPTY_VIEW_ANCHOR` when it holds no string
VIEW_ADDRESS = "_address"
#: Fortran wrapper: an omitted optional's default, materialised for a recommend routine
DEFAULT_VALUE = "_value"
#: Fortran wrapper: the flag a prologue sets, which is the prologue's contract and so a
#: convention
FORTRAN_HANDLED = CONVENTIONS.prologue_handled_arg


#: Fortran C wrapper module: a `character(1), target` an array of strings with no element is
#: viewed at, since `c_loc` may not be given the zero-size buffer itself. A module variable,
#: so an argument of the same name would hide it in the wrapper (`emit.name_clashes`)
EMPTY_VIEW_ANCHOR = "empty_view_anchor"


def local_name(argument: str, suffix: str) -> str:
    """A Fortran wrapper's local for `argument`: `<argument><suffix>`."""
    return f"{argument}{suffix}"


# -- the arguments the C ABI adds -------------------------------------------------------------
#
# `abi.c_abi` adds arguments C needs and Fortran does not have, and refuses one whose name an
# argument already takes (`_Builder._unique`). Their names are built here all the same, so the
# whole vocabulary the generator writes is in one module.

#: The length of a `character(len=*)` argument: `<name>_strlen`
STRLEN = "_strlen"


def extent_argument_name(argument: str, axis: int, rank: int) -> str:
    """The extent C is given for an assumed-shape axis: `n_<x>_elements[_dim_<k>]`."""
    suffix = f"_dim_{axis + 1}" if rank > 1 else ""
    return f"n_{argument}_elements{suffix}"


def strlen_argument_name(argument: str) -> str:
    """The length C is given for a `character(len=*)` argument."""
    return local_name(argument, STRLEN)
