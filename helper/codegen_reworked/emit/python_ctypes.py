"""Emitting the Python interface.

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

from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..ir.types import BaseType, Intent
from ..render import Writer
from .doc_numpydoc import render_docstring

#: Fortran base type -> the ctypes scalar C sees
CTYPES = {
    ("integer", "c_int"): "ctypes.c_int",
    ("integer", "c_int8_t"): "ctypes.c_int8",
    ("integer", "c_int16_t"): "ctypes.c_int16",
    ("integer", "c_int32_t"): "ctypes.c_int32",
    ("integer", "c_int64_t"): "ctypes.c_int64",
    ("integer", "c_size_t"): "ctypes.c_size_t",
    ("real", "c_float"): "ctypes.c_float",
    ("real", "c_double"): "ctypes.c_double",
    ("complex", "c_float_complex"): "ctypes.c_float",
    ("complex", "c_double_complex"): "ctypes.c_double",
    ("logical", "c_bool"): "ctypes.c_bool",
    ("character", "c_char"): "ctypes.c_char",
}

#: Fortran base type -> the numpy dtype of an array of it
DTYPES = {
    ("integer", "c_int"): "np.int32",
    ("integer", "c_int8_t"): "np.int8",
    ("integer", "c_int16_t"): "np.int16",
    ("integer", "c_int32_t"): "np.int32",
    ("integer", "c_int64_t"): "np.int64",
    ("real", "c_float"): "np.float32",
    ("real", "c_double"): "np.float64",
    ("complex", "c_float_complex"): "np.complex64",
    ("complex", "c_double_complex"): "np.complex128",
    ("logical", "c_bool"): "np.bool_",
}


def ctype_of(argument: CArgument) -> str:
    return CTYPES[(argument.type.base.value, argument.type.kind)]


def dtype_of(argument: CArgument) -> str:
    return DTYPES[(argument.type.base.value, argument.type.kind)]


class PythonEmitter:
    def __init__(self, library: str = "build/libtensor-omics.so"):
        self.library = library

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
            '"""Python interface to tensor-omics.\n'
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
        writer = Writer()
        summary = module.doc.summary or f"the {module.stripped_name} module"
        writer.block(
            f'"""Python interface to {summary}\n'
            f"\n"
            f"Generated from {module.stripped_name}. Do not edit.\n"
            f'"""'
        )
        writer.blank()
        writer.line("import ctypes")
        writer.line("import os")
        writer.blank()
        writer.line("import numpy as np")
        writer.blank()
        writer.line("from .error_handling import check_err_code")
        writer.line("from .library import load_library, nullable")
        writer.blank()
        writer.line("_lib = load_library()")
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
        return writer.render()

    def _argtype(self, argument: CArgument) -> str:
        if argument.is_scalar:
            argtype = f"ctypes.POINTER({ctype_of(argument)})"
        elif argument.type.is_character:
            argtype = "ctypes.c_char_p"
        else:
            rank = argument.rank
            # Fortran is column-major; for rank 1 the two orders coincide
            flags = "'F_CONTIGUOUS'" if rank > 1 else "'C_CONTIGUOUS'"
            argtype = (
                f"np.ctypeslib.ndpointer(dtype={dtype_of(argument)}, ndim={rank}, "
                f"flags={flags})"
            )
        if argument.optional:
            # ctypes rejects None for a checked argtype, so absence needs to be allowed
            argtype = f"nullable({argtype})"
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
            writer.block(render_docstring(wrapper))
            self._prepare_inputs(writer, wrapper)
            self._derive_extents(writer, wrapper)
            self._check_shapes(writer, wrapper)
            self._allocate(writer, wrapper)
            self._call(writer, wrapper)
            self._check_error(writer, wrapper)
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
            and not argument.is_synthesised
            and not argument.is_temporary
            and not self._is_derived(argument)
        ]

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
        name = argument.name
        writer = Writer()

        if argument.conversion is Conversion.MODE:
            body = Writer()
            body.line(f"{name} = str({name}).lower().encode()")
            return self._guarded(argument, body)

        if argument.type.is_character:
            body = Writer()
            body.line(f"{name} = str({name}).encode()")
            return self._guarded(argument, body)

        if argument.is_scalar:
            return ""

        body = Writer()
        if argument.intent is Intent.INOUT:
            # Converting would copy, and the caller would never see the modification
            body.line(
                f"if not isinstance({name}, np.ndarray) or {name}.dtype != {dtype_of(argument)}:"
            )
            with body.indent():
                body.line(
                    f'raise TypeError("\'{name}\' is modified in place, so it must '
                    f'already be a numpy array of {{}}".format({dtype_of(argument)}))'
                )
            if argument.rank > 1:
                body.line(f"if not {name}.flags.f_contiguous:")
                with body.indent():
                    body.line(
                        f'raise ValueError("\'{name}\' is modified in place, so it must '
                        f"already be column-major (order='F')\")"
                    )
        else:
            converter = "np.asfortranarray" if argument.rank > 1 else "np.ascontiguousarray"
            body.line(f"{name} = {converter}({name}, dtype={dtype_of(argument)})")
        return self._guarded(argument, body)

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

    def _derive_extents(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            expression = self._extent_expression(argument, wrapper)
            if expression:
                lines.line(f"{argument.name} = {expression}")
        if lines:
            writer.line("# what the inputs already say, rather than asking for it again")
            writer.extend(lines)
            writer.blank()

    def _extent_expression(self, argument: CArgument, wrapper: CWrapper) -> str:
        """Where a derived argument's value comes from."""
        source = argument.source
        roles = source.roles if source else None

        if argument.origin is Origin.STRLEN:
            return f"len({argument.sizes})"

        if argument.origin is Origin.EXTENT:
            owner = wrapper.argument(argument.sizes)
            if owner is not None and owner.intent.is_input:
                return f"{argument.sizes}.shape[{argument.axis}]"
            return ""

        if roles is None:
            return ""

        if roles.is_mask_count:
            return f"int({roles.mask_count_of.name}.sum())"

        if roles.is_shape_arg:
            return f"np.ascontiguousarray({roles.shape_of.name}.shape, dtype=np.int32)"

        if roles.is_extent:
            provider = self._extent_provider(argument, wrapper)
            return provider or ""

        return ""

    def _extent_provider(self, argument: CArgument, wrapper: CWrapper) -> str:
        """The first input array that already knows this extent."""
        for owner in argument.source.roles.extent_of:
            c_owner = wrapper.argument(owner.name)
            if c_owner is None or not c_owner.intent.is_input or c_owner.optional:
                continue
            axis = self._axis_of(argument.name, c_owner)
            if axis is not None:
                return f"{owner.name}.shape[{axis}]"
        return ""

    @staticmethod
    def _axis_of(extent: str, owner: CArgument) -> int | None:
        """Which numpy axis of `owner` carries `extent`.

        A character's length is a C extent but not a numpy axis: the strings are a numpy
        array of one dimension fewer.
        """
        extents = list(owner.dimension.extents)
        if owner.type.is_character and extents:
            extents = extents[1:]
        try:
            return extents.index(extent)
        except ValueError:
            return None

    def _check_shapes(self, writer: Writer, wrapper: CWrapper) -> None:
        """Cross-check every extent that more than one input claims to know.

        This is the check Fortran cannot make. `expr(n_genes, n_tissues)` and
        `weights(n_tissues)` agree by declaration, and nothing verifies it: a caller
        passing mismatched arrays gets a wrong answer or a segfault.
        """
        lines = Writer()
        for argument in wrapper:
            roles = argument.source.roles if argument.source else None
            if roles is None or not roles.is_extent:
                continue

            owners = [
                (owner, self._axis_of(argument.name, wrapper.argument(owner.name)))
                for owner in roles.extent_of
                if wrapper.argument(owner.name) is not None
                and wrapper.argument(owner.name).intent.is_input
                and not wrapper.argument(owner.name).optional
            ]
            owners = [(owner, axis) for owner, axis in owners if axis is not None]
            if len(owners) < 2:
                continue

            first, _ = owners[0]
            for owner, axis in owners[1:]:
                actual = f"{owner.name}.shape[{axis}]"
                lines.line(f"if {actual} != {argument.name}:")
                with lines.indent():
                    lines.line(
                        f'raise ValueError('
                        f'f"\'{owner.name}\' has {{{actual}}} along axis {axis}, but "'
                    )
                    with lines.indent():
                        lines.line(
                            f'f"\'{first.name}\' implies {argument.name} == '
                            f'{{{argument.name}}}"'
                        )
                    lines.line(")")
        if lines:
            writer.line("# Fortran cannot check that shared extents agree; this can")
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
        if argument.is_scalar:
            return f"{ctype_of(argument)}(0)"
        shape = ", ".join(argument.dimension.extents)
        order = "'F'" if argument.rank > 1 else "'C'"
        return f"np.zeros(({shape},), dtype={dtype_of(argument)}, order={order})"

    def _call(self, writer: Writer, wrapper: CWrapper) -> None:
        writer.line(f"_lib.{wrapper.name}(")
        with writer.indent():
            for argument in wrapper:
                writer.line(f"{self._actual(argument)},")
        writer.line(")")
        writer.blank()

    def _actual(self, argument: CArgument) -> str:
        name = argument.name
        if argument.type.is_character:
            return name
        if argument.is_scalar:
            if argument.optional:
                # None must survive to the nullable argtype
                return (
                    f"None if {name} is None "
                    f"else ctypes.byref({ctype_of(argument)}({name}))"
                )
            if argument.intent.is_output:
                return f"ctypes.byref({name})"
            return f"ctypes.byref({ctype_of(argument)}({name}))"
        return name

    def _check_error(self, writer: Writer, wrapper: CWrapper) -> None:
        error = wrapper.error_argument
        arguments = f"_{wrapper.stripped_name.upper()}_ARGUMENTS"
        writer.line(f"check_err_code({error.name}.value, {arguments})")
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
        """
        consumed = {
            a.source.roles.result_size_arg.name
            for a in wrapper
            if a.source and a.source.roles and a.source.roles.result_size_arg
        }
        return [
            argument
            for argument in wrapper
            if argument.intent is Intent.OUT
            and not argument.is_error
            and not argument.is_temporary
            and argument.name not in consumed
        ]

    def _result_expression(self, argument: CArgument, wrapper: CWrapper) -> str:
        roles = argument.source.roles if argument.source else None
        if argument.is_scalar:
            return f"{argument.name}.value"
        if roles is not None and roles.result_size_arg is not None:
            # only the leading elements were filled
            return f"{argument.name}[..., :{roles.result_size_arg.name}.value]"
        return argument.name


def python_literal(value) -> str:
    """A Fortran constant as Python source."""
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, str):
        return repr(value)
    return repr(value)
