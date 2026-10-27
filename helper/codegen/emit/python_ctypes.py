"""Emitting the Python binding.

What the generated function does, and why:

- It asks the caller only for what cannot be worked out. Extents, shape arguments, mask
  counts and work arrays all come from the real inputs, so they are not parameters.
- It checks the shapes. Fortran cannot: `expr(n_genes, n_tissues)` and `weights(n_tissues)`
  share an extent, and a mismatch is only discovered as a wrong answer or a segfault.
  Python is where that has to be caught.
- It avoids copying. `np.asarray` and friends copy only when the array is not already what
  C needs. An `intent(inout)` array is never converted at all -- a copy would silently
  discard the modification the caller asked for -- so it is checked and rejected instead.
- It returns results, not error codes. `ierr` is decoded by `check_err_code`, which raises.
"""

from __future__ import annotations

from ..abi import extents
from ..abi.c_abi import stripped_name
from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..ir.types import BaseType, Intent
from ..names import (
    PYTHON_ARGUMENT_ERROR,
    PYTHON_CAUGHT_ERROR,
    PYTHON_CHECK_ERR_CODE,
    PYTHON_CTYPES as _CTYPES,
    PYTHON_LOAD_LIBRARY,
    PYTHON_NULLABLE,
    PYTHON_NUMPY as _NP,
    PYTHON_PACKAGE_IMPORTS,
    PYTHON_STDLIB_IMPORTS,
    PYTHON_THIRD_PARTY_IMPORTS,
)
from ..render import Writer
from .error_arguments import shapes_of, sources_of
from .doc_literals import render as _render
from .doc_numpydoc import _description, render_docstring
from .shape_checks import Position, declared_shape, extent_provider, reason, shape_checks

def _python_extent(extent: str) -> str:
    """A Fortran extent expression as Python, through the one grammar the bindings evaluate.

    `generate` has refused every extent a binding computes that lies outside it
    (`abi.extents.check_binding_extents`), so this cannot fail on one that is emitted.
    """
    return extents.to_python(extent)

#: Fortran base type -> the ctypes scalar C sees, by its name in `ctypes`
CTYPES = {
    ("integer", "c_int"): "c_int",
    ("integer", "c_int8_t"): "c_int8",
    ("integer", "c_int16_t"): "c_int16",
    ("integer", "c_int32_t"): "c_int32",
    ("integer", "c_int64_t"): "c_int64",
    ("integer", "c_size_t"): "c_size_t",
    ("real", "c_float"): "c_float",
    ("real", "c_double"): "c_double",
    ("complex", "c_float_complex"): "c_float",
    ("complex", "c_double_complex"): "c_double",
    ("logical", "c_bool"): "c_bool",
    ("character", "c_char"): "c_char",
}

#: Fortran base type -> the numpy dtype of an array of it, by its name in `numpy`
DTYPES = {
    ("integer", "c_int"): "int32",
    ("integer", "c_int8_t"): "int8",
    ("integer", "c_int16_t"): "int16",
    ("integer", "c_int32_t"): "int32",
    ("integer", "c_int64_t"): "int64",
    ("real", "c_float"): "float32",
    ("real", "c_double"): "float64",
    ("complex", "c_float_complex"): "complex64",
    ("complex", "c_double_complex"): "complex128",
    ("logical", "c_bool"): "bool_",
}


def ctype_of(argument: CArgument) -> str:
    return f"{_CTYPES}.{CTYPES[(argument.type.base.value, argument.type.kind)]}"


def dtype_of(argument: CArgument) -> str:
    """The numpy dtype of an array of this type.

    A character is `S<len>`: numpy lays `S<n>` items out as n contiguous bytes each, which
    is exactly Fortran's `character(len=1) :: buf(len, n)` in column order. So the strings
    need no repacking, only encoding -- and the blank padding this layer applies itself,
    because the wrapper reads the buffer back through a `character(len=n)` pointer view.
    """
    if argument.type.is_character:
        return _character_dtype(argument)
    return f"{_NP}.{DTYPES[(argument.type.base.value, argument.type.kind)]}"


def dtype_label(argument: CArgument) -> str:
    """The dtype as a reader writes it (`np.float64`), for documentation and messages.

    The generated code binds numpy as `_np` (`names`), which a docstring has no business
    showing.
    """
    if argument.type.is_character:
        return _character_dtype(argument)
    return f"np.{DTYPES[(argument.type.base.value, argument.type.kind)]}"


def _character_dtype(argument: CArgument) -> str:
    """The `S<len>` dtype, as a Python expression.

    Three cases, because the length is known at three different times:

    - `len=8`  ->  `"S8"`, a literal
    - `len=n`  ->  `f"S{n}"`, an f-string; `n` is another argument and is in scope
    - `len=*`  ->  `"S"`, and numpy sizes the item from the data. The wrapper then reads
      the length back off the array, which is the only order that works: for an assumed
      length the string *is* where the length comes from.
    """
    length = argument.source.type.length if argument.source else None
    if length is not None and length.is_assumed and argument.intent.is_input:
        return '"S"'

    extent = character_length(argument)
    if extent.isdigit():
        return f'"S{extent}"'
    return 'f"S{' + _python_extent(extent) + '}"'


def character_length(argument: CArgument) -> str:
    """A character's length, which C carries as its leading extent."""
    return argument.dimension.extents[0]


def _pad_width(argument: CArgument) -> str | None:
    """The byte width an encoded string is blank-padded to, or None where numpy sizes it.

    Deliberately the same text `_character_dtype` interpolates into `S<n>`, so the padding
    and the item size cannot drift apart. `None` is the assumed-length *input* case, where
    the strings themselves are where the length comes from and there is no width to pad to
    until they have been encoded.
    """
    length = argument.source.type.length if argument.source else None
    if length is not None and length.is_assumed and argument.intent.is_input:
        return None
    return _python_extent(character_length(argument))


def character_rank(argument: CArgument) -> int:
    """How many dimensions the *strings* have, the length not being one of them."""
    return argument.rank - 1


class PythonEmitter:
    def __init__(self, library: str = "build/libtensor_omics.so", links=None, tiers=None):
        self.library = library
        #: resolves a Ford link to what it is called here; None in a unit test, where a link
        #: then renders as plain code rather than a cross-reference
        self.links = links
        #: what each half of an expert/allocating pair says about the other, by published
        #: name (`emit.doc_tiers`); empty where a procedure reaches this language only once
        self.tiers = tiers or {}

    # -- package ----------------------------------------------------------------

    def library_module(self) -> str:
        """Loading the shared library, in one place rather than in every module."""
        writer = Writer()
        writer.block(
            '"""Locating and loading the tensor-omics shared library.\n'
            "\n"
            "Generated. Do not edit.\n"
            '"""'
        )
        writer.blank()
        writer.line("import ctypes")
        writer.line("import os")
        writer.blank()
        writer.line("#: Where the library is expected, relative to the repository root")
        writer.line(f"DEFAULT_LIBRARY = {self.library!r}")
        writer.blank()
        writer.line("#: Overrides the search, for an installed or relocated build")
        writer.line('LIBRARY_ENV_VAR = "TENSOR_OMICS_LIBRARY"')
        writer.blank()
        writer.line("_loaded = None")
        writer.blank(collapse=False)
        writer.blank(collapse=False)
        writer.line("def nullable(argtype):")
        with writer.indent():
            writer.block(
                '"""Allow None for an argtype, so an omitted optional can be a null pointer.\n'
                "\n"
                "ctypes checks an argument against its argtype and rejects None, so an\n"
                "optional needs a type that accepts it and passes NULL through.\n"
                '"""'
            )
            writer.line("def from_param(cls, value):")
            with writer.indent():
                writer.line("return None if value is None else argtype.from_param(value)")
            writer.blank()
            writer.line("return type(")
            with writer.indent():
                writer.line('f"nullable_{getattr(argtype, \'__name__\', \'argtype\')}",')
                writer.line("(argtype,),")
                writer.line('{"from_param": classmethod(from_param)},')
            writer.line(")")
        writer.blank(collapse=False)
        writer.blank(collapse=False)
        writer.line("def load_library():")
        with writer.indent():
            writer.block(
                '"""Load the shared library once, and return it.\n'
                "\n"
                "Raises\n"
                "------\n"
                "OSError\n"
                "    If the library cannot be found, with the paths that were tried.\n"
                '"""'
            )
            writer.line("global _loaded")
            writer.line("if _loaded is not None:")
            with writer.indent():
                writer.line("return _loaded")
            writer.blank()
            writer.line("override = os.environ.get(LIBRARY_ENV_VAR)")
            writer.line("candidates = [override] if override else []")
            writer.line("here = os.path.dirname(os.path.abspath(__file__))")
            writer.line("root = os.path.dirname(os.path.dirname(here))")
            writer.line("candidates.append(os.path.join(root, DEFAULT_LIBRARY))")
            writer.blank()
            writer.line("for candidate in candidates:")
            with writer.indent():
                writer.line("if candidate and os.path.exists(candidate):")
                with writer.indent():
                    # The one library also carries the R `.Call` shims, whose R API symbols
                    # are undefined until R loads it; they are marked weak in the shims so
                    # this eager load resolves them to null (Python never calls them) while a
                    # genuinely missing symbol still fails loudly here.
                    writer.line("_loaded = ctypes.CDLL(candidate)")
                    writer.line("return _loaded")
            writer.blank()
            writer.line("raise OSError(")
            with writer.indent():
                writer.line('"cannot find the tensor-omics shared library; tried:\\n  "')
                writer.line('+ "\\n  ".join(str(c) for c in candidates)')
                writer.line(
                    '+ f"\\nbuild it, or point {LIBRARY_ENV_VAR} at it"'
                )
            writer.line(")")
        return writer.render(trailing_newline=True)

    def package_init(self, modules) -> str:
        """The package's `__init__`, re-exporting every generated function."""
        writer = Writer()
        writer.block(
            '"""Python binding to tensor-omics.\n'
            "\n"
            "Generated. Do not edit.\n"
            '"""'
        )
        writer.blank()
        writer.line("from .error_handling import (")
        with writer.indent():
            writer.line("ToxError,")
            writer.line("check_err_code,")
        writer.line(")")

        exported = []
        for module in modules:
            names = [w.stripped_name for w in module]
            exported.extend(names)
            writer.line(f"from .{module.stripped_name} import (")
            with writer.indent():
                for name in names:
                    writer.line(f"{name},")
            writer.line(")")

        writer.blank()
        writer.line("__all__ = [")
        with writer.indent():
            for name in ["ToxError", "check_err_code", *sorted(exported)]:
                writer.line(f'"{name}",')
        writer.line("]")
        return writer.render(trailing_newline=True)

    # -- module -----------------------------------------------------------------

    def module(self, module: CWrapperModule) -> str:
        # a module constant in an extent is its value here: this language cannot name it
        module = module.with_constant_sizes_evaluated()
        writer = Writer()
        # The module's whole documentation, not only its first line. A generated module is
        # the published API, so what the author wrote about the family is what a Python
        # reader should get; taking the summary alone dropped every module doc to one
        # sentence. Rendered through `_description` like any other doc, so Ford links
        # resolve to the bindings a Python caller can actually reach.
        body = Writer()
        _description(body, module.doc, self)
        rendered = body.render()
        # raw, for the same reason a procedure docstring is: an author's prose carries LaTeX
        # (`\(`, `\frac`), which is an invalid escape sequence in a plain string
        writer.line(f'r"""{module.stripped_name}')
        writer.blank()
        if rendered:
            writer.block(rendered)
            writer.blank()
        writer.line(
            f"Python binding, generated from {module.stripped_name}. Do not edit."
        )
        writer.line('"""')
        writer.blank()
        # from the name table, every one bound to a name starting with `_`, which no argument
        # can have (a Fortran name cannot start with one)
        for group in (PYTHON_STDLIB_IMPORTS, PYTHON_THIRD_PARTY_IMPORTS):
            for module_name, bound in group:
                writer.line(f"import {module_name} as {bound}")
            writer.blank()
        for module_name, names in PYTHON_PACKAGE_IMPORTS:
            imported = ", ".join(f"{name} as {bound}" for name, bound in names)
            writer.line(f"from {module_name} import {imported}")
        writer.blank()
        writer.line(f"_lib = {PYTHON_LOAD_LIBRARY}()")
        writer.blank()

        for wrapper in module:
            writer.block(self._signatures(wrapper))
            writer.blank()

        for wrapper in module:
            writer.block(self.function(wrapper))
            writer.blank(collapse=False)

        return writer.render(trailing_newline=True)

    def _signatures(self, wrapper: CWrapper) -> str:
        """The ctypes signature, declared once at import rather than on every call."""
        writer = Writer()
        writer.line(f"_lib.{wrapper.name}.restype = None")
        writer.line(f"_lib.{wrapper.name}.argtypes = (")
        with writer.indent():
            for argument in wrapper:
                writer.line(f"{self._argtype(argument)},")
        writer.line(")")
        writer.blank()
        writer.line(f"#: The wrapped procedure's arguments, so an error can name one")
        names = ", ".join(f'"{a.name}"' for a in wrapper.procedure.arguments)
        writer.line(f"_{wrapper.stripped_name.upper()}_ARGUMENTS = ({names}{',' if names else ''})")
        # and, positionally, the argument the caller passed that each was derived from -- an
        # extent names its array, so a message can blame something the caller actually wrote
        sources = sources_of(wrapper.procedure)
        if any(sources):
            rendered = ", ".join(f'"{s}"' if s else "None" for s in sources)
            writer.line("#: For a derived argument, the one the caller passed it in")
            writer.line(
                f"_{wrapper.stripped_name.upper()}_ARGUMENT_SOURCES = "
                f"({rendered}{',' if rendered else ''})"
            )
        # and the declared shape of an array the C wrapper compares, which its error states
        shapes = shapes_of(wrapper, "python")
        if any(shapes):
            rendered = ", ".join(f'"{s}"' if s else "None" for s in shapes)
            writer.line("#: For an array the C wrapper compares with its declaration, that shape")
            writer.line(
                f"_{wrapper.stripped_name.upper()}_ARGUMENT_SHAPES = "
                f"({rendered}{',' if rendered else ''})"
            )
        return writer.render()

    def _argtype(self, argument: CArgument) -> str:
        if argument.is_scalar:
            argtype = f"{_CTYPES}.POINTER({ctype_of(argument)})"
        elif argument.type.is_character:
            # c_char_p is read-only bytes, so it cannot receive an output. Both
            # directions use a buffer instead. No dtype is stated: an S<len> item size
            # may only be known at call time, and the wrapper builds the buffer itself,
            # so there is nothing here for ctypes to catch.
            argtype = f"{_NP}.ctypeslib.ndpointer(ndim={max(character_rank(argument), 1)})"
        else:
            rank = argument.rank
            # Fortran is column-major; for rank 1 the two orders coincide
            flags = "'F_CONTIGUOUS'" if rank > 1 else "'C_CONTIGUOUS'"
            argtype = (
                f"{_NP}.ctypeslib.ndpointer(dtype={dtype_of(argument)}, ndim={rank}, "
                f"flags={flags})"
            )
        if argument.optional:
            # ctypes rejects None for a checked argtype, so absence needs to be allowed
            argtype = f"{PYTHON_NULLABLE}({argtype})"
        return argtype

    # -- function ---------------------------------------------------------------

    def function(self, wrapper: CWrapper) -> str:
        writer = Writer()
        parameters = self._parameters(wrapper)

        writer.line(f"def {wrapper.stripped_name}(")
        with writer.indent(2):
            for name, default in parameters:
                writer.line(f"{name}={default}," if default else f"{name},")
        writer.line("):")

        with writer.indent():
            writer.block(render_docstring(wrapper, self))
            self._stash_raw(writer, wrapper)
            self._prepare_inputs(writer, wrapper)
            # extents a producer reads must be settled before it is called; extents read
            # off what a producer *makes* can only be settled after. So the derivation
            # brackets the producers: the plain ones first, the produced-from ones after.
            self._derive_extents(writer, wrapper, from_computed=False)
            self._compute_auto(writer, wrapper)
            self._derive_extents(writer, wrapper, from_computed=True)
            self._check_shapes(writer, wrapper)
            self._allocate(writer, wrapper)
            self._call(writer, wrapper)
            self._check_error(writer, wrapper)
            self._freeze_outputs(writer, wrapper)
            self._return(writer, wrapper)

        return writer.render()

    def _parameters(self, wrapper: CWrapper) -> list[tuple[str, str | None]]:
        """What the caller is asked for, required first as Python requires."""
        required, optional = [], []
        for argument in self._inputs(wrapper):
            if argument.optional:
                optional.append((argument.name, "None"))
            elif self._default_of(argument) is not None:
                optional.append((argument.name, self._default_of(argument)))
            else:
                required.append((argument.name, None))
        return required + optional

    def _inputs(self, wrapper: CWrapper) -> list[CArgument]:
        """The arguments a caller actually supplies."""
        return [
            argument
            for argument in wrapper
            if argument.intent.is_input
            and not argument.is_temporary
            and not self._is_derived(argument)
            and (not argument.is_synthesised or self._must_be_supplied(argument, wrapper))
        ]

    def _must_be_supplied(self, argument: CArgument, wrapper: CWrapper) -> bool:
        """Whether a synthesised extent or strlen has to come from the caller after all.

        One that sizes an input is read off that input. One that sizes an *output* --
        the length of a result only the caller can size -- has no such source, so
        leaving it out would drop it from the signature with nothing to compute it from.
        """
        if argument.is_error:
            return False
        return not self._extent_expression(argument, wrapper)

    @staticmethod
    def _is_derived(argument: CArgument) -> bool:
        roles = argument.source.roles if argument.source else None
        return bool(roles and roles.is_derived)

    @staticmethod
    def _default_of(argument: CArgument) -> str | None:
        if not argument.has_default:
            return None
        return python_literal(argument.default)

    # -- body -------------------------------------------------------------------

    def _prepare_inputs(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in self._inputs(wrapper):
            block = self._prepare(argument)
            if block:
                lines.block(block)
        if lines:
            writer.line("# accept anything array-like, converting only when C needs it")
            writer.extend(lines)
            writer.blank()

    def _prepare(self, argument: CArgument) -> str:
        """Prepare one input, reading its shape off first if that travels separately.

        The shape has to come first whichever branch below runs: preparing a character
        rebuilds the buffer to encode it, and every kind is flattened at the call, so by
        either point the caller's shape is gone.
        """
        body = self._prepare_value(argument)
        if argument.shape_arg is None:
            return body
        writer = Writer()
        writer.line(
            f"{argument.shape_arg} = {_NP}.ascontiguousarray("
            f"{_NP}.shape({argument.name}), dtype={_NP}.int32)"
        )
        if body:
            writer.block(body)
        return writer.render()

    def _prepare_value(self, argument: CArgument) -> str:
        name = argument.name
        writer = Writer()

        if argument.conversion is Conversion.MODE:
            # a mode is a character to C, so it goes through the same buffer as any
            # string; only lower-cased first, since the C wrapper matches on the
            # lower-case spelling of the parameter name. Blank-padded like every other
            # string: the C wrapper matches the buffer with `select case`, and a Fortran
            # character comparison blank-pads the shorter operand, so "iqr   " still
            # selects `case ("iqr")`.
            body = Writer()
            body.line(
                f"{name} = {_NP}.array([str({name}).lower().encode()"
                f".ljust({character_length(argument)})], "
                f"dtype={dtype_of(argument)})"
            )
            return self._guarded(argument, body)

        if argument.type.is_character:
            return self._guarded(argument, self._prepare_character(argument))

        if argument.is_scalar:
            if argument.intent is Intent.INOUT and not argument.optional:
                # Fortran writes back through the pointer, so C needs a box to write
                # into rather than a temporary made at the call. The new value is
                # returned, Python having no way to rebind the caller's name.
                body = Writer()
                body.line(f"{name} = {ctype_of(argument)}({name})")
                return body.render()
            return ""

        body = Writer()
        if argument.intent is Intent.INOUT:
            # Converting would copy, and the caller would never see the modification
            body.line(
                f"if not isinstance({name}, {_NP}.ndarray) or {name}.dtype != {dtype_of(argument)}:"
            )
            with body.indent():
                body.line(_refusal(
                    f'"\'{name}\' is modified in place, so it must '
                    f'already be a numpy array of {{}}".format({dtype_of(argument)})', name))
            self._check_rank(argument, body)
            if argument.rank > 1:
                body.line(f"if not {name}.flags.f_contiguous:")
                with body.indent():
                    body.line(_refusal(
                        f'"\'{name}\' is modified in place, so it must '
                        f"already be column-major (order='F')\"", name))
        else:
            # numpy raises without saying which argument was wrong -- it cannot know --
            # so the name is put back on the way out
            converter = (f"{_NP}.asfortranarray"
                         if argument.rank > 1 or argument.shape_arg is not None
                         else f"{_NP}.ascontiguousarray")
            body.line("try:")
            with body.indent():
                body.line(f"{name} = {converter}({name}, dtype={dtype_of(argument)})")
            body.line(f"except (TypeError, ValueError) as {PYTHON_CAUGHT_ERROR}:")
            with body.indent():
                body.line(_refusal(
                    f'f"\'{name}\' must be an array of {dtype_label(argument)}: '
                    f'{{{PYTHON_CAUGHT_ERROR}}}"', name, from_none=True))
            self._check_rank(argument, body)
        return self._guarded(argument, body)

    @staticmethod
    def _numpy_rank(argument: CArgument) -> int:
        """How many dimensions the numpy array for this argument has.

        A character's length is a C extent but not a numpy axis, and a scalar string is
        still carried as a one-element buffer.
        """
        if argument.type.is_character:
            return max(character_rank(argument), 1)
        return argument.rank

    def _check_rank(self, argument: CArgument, body: Writer) -> None:
        """Reject a wrong-rank array by name.

        The `ndpointer` argtype checks this too, but only at the call, and it raises a
        `ctypes.ArgumentError` naming an argument *position*. Checking here says which
        argument, and raises the ValueError the rest of the wrapper raises.
        """
        if argument.shape_arg is not None:
            # its shape travels in a separate argument, so any rank is accepted and the
            # array is flattened at the call
            return
        rank = self._numpy_rank(argument)
        name = argument.name
        body.line(f"if {name}.ndim != {rank}:")
        with body.indent():
            body.line(_refusal(
                f'f"\'{name}\' must have {rank} dimension{"" if rank == 1 else "s"}, '
                f'but has {{{name}.ndim}}"', name))

    def _prepare_character(self, argument: CArgument) -> Writer:
        """Encode strings into the blank-padded buffer C expects.

        The Fortran wrapper takes a `character(len=strlen)` pointer view of this buffer
        rather than converting it, so whatever pads a short string *is* part of the value.
        Blanks are the padding Fortran's own `trim`, `len_trim` and string comparison are
        defined on, and no UTF-8 continuation byte (0x80-0xBF) can be mistaken for one.

        The padding goes on the **encoded bytes**, never on the `str`: the width is a byte
        count everywhere else in this layer, and `str.ljust` would measure it in characters
        -- after which numpy would re-pad a multibyte string with the nulls this replaced.
        """
        name = argument.name
        width = _pad_width(argument)
        body = Writer()
        if character_rank(argument) == 0:
            # Where numpy sizes the item, ljust(1) is the whole job: it is a no-op for
            # every non-empty string, and it stops an empty one becoming numpy's forced
            # single null byte.
            body.line(
                f"{name} = {_NP}.array([str({name}).encode().ljust({width or 1})], "
                f"dtype={dtype_of(argument)})"
            )
        else:
            body.line("try:")
            with body.indent():
                items = (f"{_NP}.asarray({name}).ravel(order='F')"
                         if argument.shape_arg is not None else name)
                if width is None:
                    # numpy sizes the item from the longest element, so the width is not
                    # known until every string is encoded. `default=0` covers an empty
                    # sequence and `or 1` an array of nothing but empty strings, either of
                    # which would otherwise ask numpy for a width of zero.
                    body.line(f"_{name}_bytes = [str(_s).encode() for _s in {items}]")
                    body.line(f"_{name}_width = max(map(len, _{name}_bytes), default=0) or 1")
                    body.line(
                        f"{name} = {_NP}.asarray("
                        f"[_b.ljust(_{name}_width) for _b in _{name}_bytes], "
                        f"dtype={dtype_of(argument)})"
                    )
                else:
                    body.line(
                        f"{name} = {_NP}.asarray("
                        f"[str(_s).encode().ljust({width}) for _s in {items}], "
                        f"dtype={dtype_of(argument)})"
                    )
            body.line(f"except TypeError as {PYTHON_CAUGHT_ERROR}:")
            with body.indent():
                body.line(_refusal(
                    f'f"\'{name}\' must be a sequence of strings: {{{PYTHON_CAUGHT_ERROR}}}"',
                    name, from_none=True))
            self._check_rank(argument, body)
        return body

    @staticmethod
    def _guarded(argument: CArgument, body: Writer) -> str:
        """Wrap in a presence check when the argument may be omitted."""
        if not argument.optional:
            return body.render()
        writer = Writer()
        writer.line(f"if {argument.name} is not None:")
        with writer.indent():
            writer.extend(body)
        return writer.render()

    @staticmethod
    def _producer_value(supply, stashes: dict[str, str]) -> str:
        """What to pass for one producer input: an argument, or a constant outright."""
        if supply.is_constant:
            return python_literal(supply.constant)
        return stashes.get(supply.argument, supply.argument)

    @staticmethod
    def _rebinds_on_prepare(argument: CArgument) -> bool:
        """Whether preparing this argument replaces the caller's value with a buffer."""
        return argument.conversion is Conversion.MODE or argument.type.is_character

    def _raw_stashes(self, wrapper: CWrapper) -> dict[str, str]:
        """Producer inputs that preparation would destroy, and the name to keep them under.

        A producer runs after the extents are derived -- it may need one of them -- but
        by then a character argument has been replaced by a numpy bytes buffer, and the
        producer's own `str()` would turn that into a nonsense value. So the caller's
        original is kept aside first.
        """
        stashes: dict[str, str] = {}
        for argument in wrapper:
            plan = argument.source.roles.computed_from if argument.source and argument.source.roles else None
            if plan is None or not plan.is_automatic:
                continue
            for supply in plan.inputs:
                if supply.is_constant:
                    continue
                name = supply.argument
                supplied = wrapper.argument(name)
                if supplied is not None and self._rebinds_on_prepare(supplied):
                    stashes[name] = f"_{name}_raw"
        return stashes

    def _stash_raw(self, writer: Writer, wrapper: CWrapper) -> None:
        stashes = self._raw_stashes(wrapper)
        if not stashes:
            return
        writer.line("# kept before conversion, for the producers called below")
        for name, alias in stashes.items():
            writer.line(f"{alias} = {name}")
        writer.blank()

    def _compute_auto(self, writer: Writer, wrapper: CWrapper) -> None:
        """Fill each DM_OUTPUT_FROM(AUTO) argument by calling the producer's wrapper.

        The producer's own generated function does its own validation and error
        checking, so this just calls it. Its inputs -- the extents it reads -- are already
        derived (those that do not themselves come from a producer are), so it just passes
        them, taking the stashed pre-conversion value of any input preparation replaced.

        One call serves every argument it supplies. Producers are not cheap -- the one
        behind `read_tox_data_into` unpacks a zip archive -- so calling it once per
        extent would re-do that work for each.
        """
        lines = Writer()
        calls: dict[tuple, str] = {}
        stashes = self._raw_stashes(wrapper)
        consumer_module = wrapper.procedure.module
        for argument in wrapper:
            plan = argument.source.roles.computed_from if argument.source and argument.source.roles else None
            if plan is None or not plan.is_automatic:
                continue
            producer = stripped_name(plan.producer)
            if plan.producer.module is not consumer_module:
                # Imported here rather than at the top of the file: two modules may size
                # each other's outputs, and a module-level import would then be circular.
                lines.line(f"from .{plan.producer.module.name} import {producer}")
            # the keyword is the *producer's* parameter name, the value the consumer's
            # argument; they coincide under name-matching but not when a table renames one
            call_args = ", ".join(
                f"{supply.name}={self._producer_value(supply, stashes)}"
                for supply in plan.inputs)
            single = self._producer_returns_single(plan.producer)
            if single:
                lines.line(f"{argument.name} = {producer}({call_args})")
                continue

            key = (producer, call_args)
            if key not in calls:
                calls[key] = f"_{producer}_result" if len(calls) == 0 else \
                    f"_{producer}_result_{len(calls)}"
                lines.line(f"{calls[key]} = {producer}({call_args})")
                for supply in plan.inout_feedback:
                    # the producer refined a value we handed it (e.g. capped it); adopt it
                    lines.line(f'{supply.argument} = {calls[key]}["{supply.name}"]')
            lines.line(f'{argument.name} = {calls[key]}["{plan.output.name}"]')
        if lines:
            writer.line("# work out what other procedures must supply, per DM_OUTPUT_FROM")
            writer.extend(lines)
            writer.blank()

    @staticmethod
    def _producer_returns_single(producer) -> bool:
        """Whether the producer's wrapper returns one value (not a dict of outputs)."""
        outputs = [
            a for a in producer.arguments
            if a.intent and a.intent.is_output
            and a.name.lower() != "ierr"
            and not (a.roles and a.roles.is_temporary)
        ]
        if producer.result is not None:
            outputs.append(producer.result)
        return len(outputs) == 1

    def _derive_extents(self, writer: Writer, wrapper: CWrapper,
                        from_computed: bool) -> None:
        lines = Writer()
        for argument in wrapper:
            if self._derives_from_computed(argument, wrapper) != from_computed:
                continue
            expression = self._extent_expression(argument, wrapper)
            if expression:
                lines.line(f"{argument.name} = {expression}")
        if lines:
            writer.line("# what the inputs already say, rather than asking for it again")
            writer.extend(lines)
            writer.blank()

    def _derives_from_computed(self, argument: CArgument, wrapper: CWrapper) -> bool:
        """Whether this argument's derived value is read off something a producer makes.

        `n_elements = prod(arr_shape)` is, when `arr_shape` is filled by a DM_OUTPUT_FROM
        producer: it cannot be settled until that producer has run. A plain extent read
        off a caller's array is not, and is settled before the producers so they can use
        it. The two are derived on opposite sides of the producer calls.
        """
        if argument.origin in (Origin.EXTENT, Origin.STRLEN):
            owner = wrapper.argument(argument.sizes)
            return owner is not None and self._is_computed(owner)
        roles = argument.source.roles if argument.source else None
        if roles is None or not roles.is_extent:
            return False
        for owner in roles.extent_of:
            c_owner = wrapper.argument(owner.name)
            if c_owner is None:
                continue
            if self._is_computed(c_owner):
                return True
            shape = wrapper.argument(c_owner.shape_arg) if c_owner.shape_arg else None
            if shape is not None and self._is_computed(shape):
                return True
        return False

    @staticmethod
    def _is_computed(argument: CArgument) -> bool:
        roles = argument.source.roles if argument.source else None
        return bool(roles and roles.is_computed)

    def _extent_expression(self, argument: CArgument, wrapper: CWrapper) -> str:
        """Where a derived argument's value comes from."""
        source = argument.source
        roles = source.roles if source else None

        if argument.origin is Origin.STRLEN:
            # the buffer's item size, not the number of items in it. For len=* numpy
            # sized the item from the string itself, so this is where the length is.
            # Only for an input: an output buffer does not exist yet, so nothing can be
            # read off it and the caller has to say how long its strings are.
            owner = wrapper.argument(argument.sizes)
            if owner is not None and owner.intent.is_input:
                return self._maybe_absent(owner, f"{argument.sizes}.itemsize")
            return ""

        if argument.origin is Origin.EXTENT:
            owner = wrapper.argument(argument.sizes)
            # not from a work array: that is allocated from this extent, not the reverse
            if owner is not None and owner.intent.is_input and not owner.is_temporary:
                return self._maybe_absent(
                    owner, f"{argument.sizes}.shape[{argument.axis}]")
            return ""

        if roles is None:
            return ""

        if roles.is_mask_count:
            return f"int({roles.mask_count_of.name}.sum())"

        if roles.is_shape_arg:
            # settled while its array was prepared, or asked of the caller when the array
            # is an output and the shape is what they want produced
            return ""

        if roles.is_extent:
            # the same predicate the signature uses, so an argument is never both asked
            # for and derived
            if not roles.is_inferable_extent:
                return ""
            provider = self._extent_provider(argument, wrapper)
            return provider or ""

        return ""

    @staticmethod
    def _maybe_absent(owner: CArgument, expression: str) -> str:
        """Guard an extent read off an argument the caller may have left out.

        An omitted optional array is `None`, and asking a `None` for its shape is an
        AttributeError rather than the absent argument the Fortran expects. Zero is
        what the C wrapper passes on for an absent array anyway.
        """
        if not owner.optional:
            return expression
        return f"0 if {owner.name} is None else {expression}"

    def _extent_provider(self, argument: CArgument, wrapper: CWrapper) -> str:
        """The first input array that already knows this extent.

        Only an input's shape may be read: an output array is allocated *from* this
        extent, so it does not exist yet when the extent is settled. When an extent is
        named by both (e.g. `reps_per_tissue(n_tissues)` and the output
        `log_transformed_expr(n_tissues, n_genes)`) the input is the one to read. Which input
        is `shape_checks.extent_provider`'s to say, so that the axis it is read from is the one
        `_check_shapes` leaves out.
        """
        position = extent_provider(argument, wrapper)
        if position is not None:
            return self._position(position)
        # an output sized by a shape argument: the count is the product of that shape
        for owner in argument.source.roles.extent_of:
            c_owner = wrapper.argument(owner.name)
            if c_owner is not None and c_owner.shape_arg is not None:
                return f"int({_NP}.prod({c_owner.shape_arg}))"
        return ""

    @staticmethod
    def _position(position: Position) -> str:
        """An axis of an input as Python reads it, or its element count."""
        name = position.array.name
        if position.axis is None:
            # its extents travel separately, so the only thing an extent of it can mean is
            # how many elements there are altogether
            return f"{name}.size"
        return f"{name}.shape[{position.axis}]"

    def _check_shapes(self, writer: Writer, wrapper: CWrapper) -> None:
        """Compare every input axis with what it must be (`emit.shape_checks`).

        This is the check Fortran cannot make. `expr(n_genes, n_tissues)` and
        `weights(n_tissues)` agree by declaration, and nothing verifies it: a caller
        passing mismatched arrays gets a wrong answer, or Fortran reads past the buffer.
        """
        lines = Writer()
        for check in shape_checks(wrapper):
            array = check.position.array
            actual = self._position(check.position)
            expected = check.expected
            value = str(expected.value) if expected.origin == "constant" else expected.extent
            condition = f"{actual} != {value}"
            if array.optional:
                condition = f"{array.name} is not None and {condition}"
            where = ("elements" if check.position.axis is None
                     else f"along axis {check.position.axis}")
            shape = _render(declared_shape(array), "python")
            lines.line(f"if {condition}:")
            with lines.indent():
                lines.line(_refusal(
                    f'f"\'{array.name}\' has {{{actual}}} {where}, but {reason(expected)}'
                    f'{{{value}}}; \'{array.name}\' is declared with shape ({shape})"', array.name))
        if lines:
            writer.line("# Fortran cannot check that an array has the shape it is declared with; this can")
            writer.extend(lines)
            writer.blank()

    def _allocate(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            if argument.intent.is_input and not argument.is_temporary:
                continue
            if self._is_derived(argument) and argument.intent.is_input:
                continue
            lines.line(f"{argument.name} = {self._new_value(argument)}")
        if lines:
            writer.line("# outputs and work arrays, which the caller never sees")
            writer.extend(lines)
            writer.blank()

    def _new_value(self, argument: CArgument) -> str:
        """A fresh output or work array.

        `np.empty`, not `np.zeros`: Fortran fills these, so zeroing them is an O(n) write
        of data that is immediately overwritten. It is also more honest -- where Fortran
        fills only part of an array, zeros look like results, while uninitialised memory
        looks like the garbage it is. On the error path nothing is returned at all.
        """
        if argument.is_scalar:
            return f"{ctype_of(argument)}(0)"

        if argument.shape_arg is not None:
            # C sees one flat block, so the element count is the product of the extents
            # that travel beside it
            count = f"int({_NP}.prod({argument.shape_arg}))"
            if argument.type.is_character:
                return f"{_NP}.zeros(({count},), dtype={dtype_of(argument)})"
            return f"{_NP}.empty(({count},), dtype={dtype_of(argument)}, order='C')"

        if argument.type.is_character:
            # zeros, not empty, and this is the one case where it matters. Fortran writes
            # through a pointer view, so it blank-fills the whole width of every element it
            # assigns -- but it may assign none at all (an early error return, an array it
            # fills only partly). The zeros are what makes such an element read back as ""
            # rather than as garbage: numpy drops the trailing nulls of an `S` item.
            # Not `np.full(..., b" ")`, which fills each item with one blank and nulls.
            shape = ", ".join(_python_extent(e) for e in argument.dimension.extents[1:]) or "1"
            return f"{_NP}.zeros(({shape},), dtype={dtype_of(argument)})"

        shape = ", ".join(_python_extent(e) for e in argument.dimension.extents)
        order = "'F'" if argument.rank > 1 else "'C'"
        return f"{_NP}.empty(({shape},), dtype={dtype_of(argument)}, order={order})"

    def _call(self, writer: Writer, wrapper: CWrapper) -> None:
        writer.line(f"_lib.{wrapper.name}(")
        with writer.indent():
            for argument in wrapper:
                writer.line(f"{self._actual(argument)},")
        writer.line(")")
        writer.blank()

    def _actual(self, argument: CArgument) -> str:
        name = argument.name
        if argument.shape_arg is not None and not argument.type.is_character:
            # Fortran has it as one flat block whose extents arrive separately; ravel in
            # Fortran order, which is a view of an already column-major array
            return f"{name}.ravel(order='F')"
        if argument.type.is_character:
            return name
        if argument.is_scalar:
            if argument.optional:
                # None must survive to the nullable argtype
                return (
                    f"None if {name} is None "
                    f"else {_CTYPES}.byref({ctype_of(argument)}({name}))"
                )
            if argument.intent.is_output:
                return f"{_CTYPES}.byref({name})"
            return f"{_CTYPES}.byref({ctype_of(argument)}({name}))"
        return name

    def _check_error(self, writer: Writer, wrapper: CWrapper) -> None:
        error = wrapper.error_argument
        prefix = f"_{wrapper.stripped_name.upper()}"
        tables = [f"{prefix}_ARGUMENTS"]
        if any(sources_of(wrapper.procedure)):
            tables.append(f"{prefix}_ARGUMENT_SOURCES")
        if any(shapes_of(wrapper, "python")):
            tables.append(f"shapes={prefix}_ARGUMENT_SHAPES")
        writer.line(f"{PYTHON_CHECK_ERR_CODE}({error.name}.value, {', '.join(tables)})")
        writer.blank()

    def _freeze_outputs(self, writer: Writer, wrapper: CWrapper) -> None:
        """Hand back results the caller cannot modify by accident.

        A result is a value, not a scratch buffer. Setting the flag is O(1) -- no copy --
        and it propagates into the views `_result_expression` builds, so freezing the
        buffer freezes the reshape of a serialized output and the slice of a partially
        filled one along with it.

        Only what this wrapper allocated. An `intent(inout)` array belongs to the caller
        and is modified in place, so freezing it would break the very thing it exists for;
        scalars and decoded strings are immutable in Python already, and a character buffer
        is never returned -- only the strings decoded out of it, which are fresh objects
        per call and so have nothing to protect.

        After the call rather than before: `ctypes` is handed the raw data pointer, so a
        flag set earlier would not stop Fortran writing -- it would only misdescribe it.

        R needs no counterpart. It is copy-on-modify, so a returned vector cannot be
        aliased by anything else, and there the risk this guards against does not exist.
        """
        frozen = [
            argument for argument in self._outputs(wrapper)
            if argument.is_array and not argument.type.is_character
        ]
        if not frozen:
            return
        writer.line("# a result is a value: modify a copy, not this")
        for argument in frozen:
            writer.line(f"{argument.name}.flags.writeable = False")
        writer.blank()

    def _return(self, writer: Writer, wrapper: CWrapper) -> None:
        outputs = self._outputs(wrapper)
        if not outputs:
            writer.line("return None")
            return

        expressions = {a.name: self._result_expression(a, wrapper) for a in outputs}
        if len(outputs) == 1:
            writer.line(f"return {expressions[outputs[0].name]}")
            return

        writer.line("return {")
        with writer.indent():
            for argument in outputs:
                writer.line(f'"{argument.name}": {expressions[argument.name]},')
        writer.line("}")

    def _outputs(self, wrapper: CWrapper) -> list[CArgument]:
        """What the caller gets back.

        Not `ierr`, which became an exception. Not work arrays, which are an
        implementation detail. Not a count that only exists to size another result --
        the returned array is already that long.

        An `intent(inout)` array is modified in place and so is already visible to the
        caller, but an `intent(inout)` scalar cannot be: it is returned instead.
        """
        consumed = {
            a.source.roles.result_size_arg.name
            for a in wrapper
            if a.source and a.source.roles and a.source.roles.result_size_arg
        }
        return [
            argument
            for argument in wrapper
            if (argument.intent is Intent.OUT
                or (argument.intent is Intent.INOUT and argument.is_scalar
                    and not argument.optional))
            and not argument.is_error
            and not argument.is_temporary
            and argument.name not in consumed
        ]

    def _result_expression(self, argument: CArgument, wrapper: CWrapper) -> str:
        roles = argument.source.roles if argument.source else None
        if argument.shape_arg is not None and not argument.is_scalar:
            # C had it as one flat block; the extents that travelled beside it are the
            # shape, so the caller gets the n-d array back rather than the flat buffer
            # and a shape they have to reapply themselves. Column-major, as Fortran wrote.
            if argument.type.is_character:
                # `dtype=str`, or an array of no strings comes back as float64 -- what numpy
                # makes of an empty list -- where one with strings is `<U...` (#221 B3)
                decoded = f"[_s.decode().rstrip(' ') for _s in {argument.name}]"
                return (f"{_NP}.asarray({decoded}, dtype=str)"
                        f".reshape(tuple({argument.shape_arg}), order='F')")
            return f"{argument.name}.reshape(tuple({argument.shape_arg}), order='F')"
        if argument.type.is_character:
            # numpy hands back the bytes with any trailing nulls already dropped; Fortran
            # blank-padded them to the buffer width, so the caller's str is what is left
            # after the blanks. `rstrip(' ')` and not a bare `rstrip()`: that would also eat
            # tabs and newlines, which travel through this API as data.
            if character_rank(argument) == 0:
                return f"{argument.name}[0].decode().rstrip(' ')"
            return f"[_s.decode().rstrip(' ') for _s in {argument.name}]"
        if argument.is_scalar:
            return f"{argument.name}.value"
        if roles is not None and roles.result_size_arg is not None:
            # only the leading elements were filled
            return f"{argument.name}[..., :{roles.result_size_arg.name}.value]"
        return argument.name


def _refusal(message: str, argument: str, from_none: bool = False) -> str:
    """The statement refusing an argument: `message` is a Python expression giving the text.

    Every refusal the binding makes itself raises `ToxArgumentError`, which carries no tox code
    and is a `TypeError` and a `ValueError` as well (`emit.errors_python`), naming the argument
    so `error.argument` says which one, as it does for an error the Fortran reports.
    """
    statement = f'raise {PYTHON_ARGUMENT_ERROR}({message}, argument="{argument}")'
    return f"{statement} from None" if from_none else statement


def python_literal(value) -> str:
    """A Fortran constant as Python source."""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, str):
        return repr(value)
    return repr(value)
