"""The names a published argument may not take, in one table.

A dummy argument's name is published verbatim: it is a Python parameter, an R parameter, a C
parameter of the R `.Call` shim, and a dummy of every generated Fortran wrapper. So it has to
be a name each of those languages accepts there, and it must not collide with a name the
generated code introduces itself. Both are checked by `ir.validate`, from this table, and the
emitters build their own local names from the same table -- so a local an emitter starts
generating is a name the check knows about, rather than a list in `validate` that silently
goes stale (`tests/test_names.py` reads the generated code back to hold the emitters to it).

What is deliberately *not* here:

- Fortran keywords. Fortran has no reserved words, and main publishes `value`.
- R function names. An argument cannot shadow a function call in R -- R skips non-function
  bindings when it looks a call up -- and main publishes `mode`, `vector`, `array`,
  `expression`, `factor`, `points` and `weights`.
- Python soft keywords (`match`, `case`, `type`, `_`), which are valid parameter names. `type`
  is a builtin, though, so it is refused as one.

A later change that makes an emitter introduce a new name adds it here, in the same commit.
"""

from __future__ import annotations

import builtins
import keyword

from .config import CONVENTIONS

# -- the words of the binding languages -------------------------------------------------------

#: Python's hard keywords. A parameter named like one is a SyntaxError, and since the package
#: `__init__` imports every module, `import tensor_omics` fails as a whole.
PYTHON_KEYWORDS = frozenset(keyword.kwlist)

#: Python's builtins (`len`, `max`, `str`, `input`, ...). The generated functions call several
#: of them, and a parameter of the same name hides the builtin for the whole function body.
PYTHON_BUILTINS = frozenset(name for name in vars(builtins) if not name.startswith("_"))

#: R's reserved words (`?Reserved`). A parameter named like one is a parse error, which makes
#: `r/load_tensor_omics.R` fail for every wrapper, not only the one that has it.
R_RESERVED = frozenset({
    "if", "else", "repeat", "while", "function", "for", "next", "break", "in",
    "TRUE", "FALSE", "NULL", "Inf", "NaN",
    "NA", "NA_integer_", "NA_real_", "NA_complex_", "NA_character_",
})

#: C's keywords, C23 included. Every published argument is a parameter of the R `.Call` shim,
#: and the shims are compiled into the one library, so a C keyword fails the library build.
#: The `_Upper` keywords are listed for completeness; no Fortran name can start with `_`.
C_KEYWORDS = frozenset({
    "auto", "break", "case", "char", "const", "continue", "default", "do", "double",
    "else", "enum", "extern", "float", "for", "goto", "if", "inline", "int", "long",
    "register", "restrict", "return", "short", "signed", "sizeof", "static", "struct",
    "switch", "typedef", "union", "unsigned", "void", "volatile", "while",
    "alignas", "alignof", "bool", "constexpr", "false", "nullptr", "static_assert",
    "thread_local", "true", "typeof", "typeof_unqual",
    "_Alignas", "_Alignof", "_Atomic", "_BitInt", "_Bool", "_Complex", "_Decimal128",
    "_Decimal32", "_Decimal64", "_Generic", "_Imaginary", "_Noreturn", "_Static_assert",
    "_Thread_local",
})

# -- what the generated code binds itself -----------------------------------------------------

#: The modules a generated Python module imports, as (module, the name it is bound to)
PYTHON_IMPORTS = (("ctypes", "ctypes"), ("os", "os"), ("numpy", "np"))
#: What a generated Python module imports from the package itself, as (module, names)
PYTHON_PACKAGE_IMPORTS = (
    (".error_handling", ("check_err_code",)),
    (".library", ("load_library", "nullable")),
)
#: The exception a generated Python function catches while converting an input
PYTHON_CAUGHT_ERROR = "error"

#: Every name a generated Python module binds that a parameter could hide. Its own private
#: names all start with `_`, which no Fortran name can.
PYTHON_BOUND = frozenset(
    [bound for _, bound in PYTHON_IMPORTS]
    + [name for _, names in PYTHON_PACKAGE_IMPORTS for name in names]
    + [PYTHON_CAUGHT_ERROR]
)

#: The R `.Call` shim's count of protected objects, a local of every shim function
C_PROTECT_COUNT = "nprot"

#: The flag a generated Fortran wrapper declares for its prologue to set, which is the
#: prologue's contract and so a convention
FORTRAN_HANDLED = CONVENTIONS.prologue_handled_arg

#: Where a generated local lives: the R binding's C, or a Fortran wrapper (either tier, or the
#: C wrapper). Fortran compares names ignoring case.
R_SHIM = "R shim"
FORTRAN = "Fortran wrapper"

#: Every whole name a generated wrapper declares beside the arguments, by where it is declared
GENERATED_LOCALS = {
    R_SHIM: frozenset({C_PROTECT_COUNT}),
    FORTRAN: frozenset({FORTRAN_HANDLED}),
}

# -- the locals generated per argument --------------------------------------------------------
#
# A wrapper needs a local beside some arguments: the R shim pulls a scalar out of its SEXP into
# `<name>_v`, the Fortran C wrapper views a string through `<name>_f`. Another dummy called
# `<name>_v` would then be declared twice. Every such local is `<argument><suffix>`, built
# with `local_name`, and `validate` refuses an argument named like another one plus a suffix.

#: R shim: a by-value scalar input, pulled out of its length-1 vector
SCALAR_VALUE = "_v"
#: R shim: an optional, as the pointer C is handed (null when omitted)
OPTIONAL_POINTER = "_p"
#: R shim: an optional's length (0 when omitted)
OPTIONAL_SIZE = "_size"
#: R shim: an optional array's extents (all 0 when omitted); also the `dim` attribute built
#: for a matrix output, which is block-scoped
DIM = "_dim"
#: R shim: an optional array's `dim` attribute, as read
OPTIONAL_DIM_ATTRIBUTE = "_d"
#: R shim: the C buffer a logical or a string is converted into
BUFFER = "_c"
#: R shim: the copy of an array modified in place, so the caller's stays intact
INOUT_COPY = "_out"
#: Fortran C wrapper: the local a converted argument is passed to the callee as
CONVERTED = "_f"
#: Fortran C wrapper: the integer parameter a mode string names
MODE_VALUE = "_mode_f"
#: Fortran wrapper: an omitted optional's default, materialised for a recommend routine
DEFAULT_VALUE = "_value"

#: Every suffix above, with where its local is declared. `_out` + `_dim` is the `dim` attribute
#: of an in-place copy, built as one name.
GENERATED_SUFFIXES = {
    SCALAR_VALUE: R_SHIM, OPTIONAL_POINTER: R_SHIM, OPTIONAL_SIZE: R_SHIM, DIM: R_SHIM,
    OPTIONAL_DIM_ATTRIBUTE: R_SHIM, BUFFER: R_SHIM, INOUT_COPY: R_SHIM,
    INOUT_COPY + DIM: R_SHIM,
    CONVERTED: FORTRAN, MODE_VALUE: FORTRAN, DEFAULT_VALUE: FORTRAN,
}

# -- the arguments the C ABI adds -------------------------------------------------------------
#
# `abi.c_abi` adds arguments C needs and Fortran does not have, and refuses one whose name an
# argument already takes (`_Builder._unique`). Their names are built here all the same, so the
# whole vocabulary the generator writes is in one table.

#: The length of a `character(len=*)` argument: `<name>_strlen`
STRLEN = "_strlen"


def local_name(argument: str, suffix: str) -> str:
    """The local a wrapper declares for `argument`: `<argument><suffix>`."""
    return f"{argument}{suffix}"


def extent_argument_name(argument: str, axis: int, rank: int) -> str:
    """The extent C is given for an assumed-shape axis: `n_<x>_elements[_dim_<k>]`."""
    suffix = f"_dim_{axis + 1}" if rank > 1 else ""
    return f"n_{argument}_elements{suffix}"


def strlen_argument_name(argument: str) -> str:
    """The length C is given for a `character(len=*)` argument."""
    return local_name(argument, STRLEN)


def language_word(name: str) -> str | None:
    """What makes `name` unusable as a published argument, or None.

    Case-sensitive, as Python, R and C are: the published name keeps its declared spelling, so
    `NA` is R's reserved word while `na` is a name.
    """
    if name in PYTHON_KEYWORDS:
        return "a Python keyword"
    if name in R_RESERVED:
        return "an R reserved word"
    if name in C_KEYWORDS:
        return "a C keyword"
    if name in PYTHON_BUILTINS:
        return "a Python builtin"
    if name in PYTHON_BOUND:
        return "a name the generated Python module binds itself"
    return None


def generated_local(name: str) -> str | None:
    """Which generated wrapper declares a local called `name`, or None."""
    for where, locals_ in GENERATED_LOCALS.items():
        if where == FORTRAN:
            if name.lower() in {local.lower() for local in locals_}:
                return where
        elif name in locals_:
            return where
    return None


def clashing_local(name: str, others) -> tuple[str, str, str] | None:
    """The (other argument, suffix, wrapper) whose generated local `name` would be, or None.

    Compared ignoring case: a Fortran local ignores it, and for a C local refusing `K_v` beside
    `k` costs nobody anything.
    """
    lowered = name.lower()
    for other in others:
        for suffix, where in GENERATED_SUFFIXES.items():
            if lowered == local_name(other, suffix).lower():
                return other, suffix, where
    return None
