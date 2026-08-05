"""Emitting the C half of the R binding.

The split with R is the same one the former Rcpp layer used and is documented in
`design/language-layers.md`: **R decides and raises, C marshals and calls.** By the time a
value reaches C it has been checked and coerced in R, so nothing here validates -- it
converts what Fortran cannot take from R directly (R's int-based logicals, its strings) and
hands everything to the `extern "C"` Fortran wrapper.

Each generated function is a `.Call` entry `<name>_call`, registered in `init.c`. It is
internal: the user-facing `<name>` is the R wrapper on top of it, which calls
`.Call("<name>_call", ...)` and decodes the returned list of every output plus `ierr`.

Only R's C API is used -- no Rcpp, no cpp11. The two conversions C cannot avoid (logical
arrays, strings) go through `tox_marshal.h`, whose buffers come from `R_alloc` and so are
freed automatically when the call returns; there is nothing to leak on an error longjmp.
"""

from __future__ import annotations

import re

from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..ir.types import BaseType, Intent
from ..render import Writer

#: The R shared object this binding is registered into (drives `R_init_<name>`).
R_DLL_NAME = "tensoromics"

#: The R shims compile into the one `libtensor-omics.so`. This guard drops them (and their
#: R.h include) when the binding is built without R (`NO_R_BINDING`) or without the C ABI
#: they call (`NO_C_BINDING`), leaving empty objects that need no R headers.
_GUARD_OPEN = "#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)"
_GUARD_CLOSE = "#endif  // R binding"

#: Every R API symbol the shims reference. R provides them at dyn.load, but the one
#: libtensor-omics.so also loads into a non-R host (Python via ctypes). Marking them all weak
#: lets that load succeed under eager binding -- the undefined R symbols resolve to null, and
#: the R-only code that would use them never runs from Python -- while a *genuinely* missing
#: symbol (a build regression) still fails loudly at load. R binds them to its real
#: definitions at dyn.load. The `test_ctypes_loads_the_r_bundled_library` end-to-end test
#: guards this list: if the emitter starts using an R symbol not here, that load fails.
#: (NA_STRING is R_NaString; PROTECT/UNPROTECT are Rf_protect/Rf_unprotect; CHAR is R_CHAR.)
_WEAK_R_SYMBOLS = (
    "COMPLEX", "INTEGER", "LENGTH", "LOGICAL", "REAL", "STRING_ELT", "TYPEOF", "XLENGTH",
    "SET_STRING_ELT", "SET_VECTOR_ELT",
    "R_CHAR", "R_alloc", "R_DimSymbol", "R_NamesSymbol", "R_NilValue", "R_NaString",
    "R_registerRoutines", "R_useDynamicSymbols",
    "Rf_allocVector", "Rf_asInteger", "Rf_asLogical", "Rf_asReal", "Rf_coerceVector",
    "Rf_duplicate", "Rf_getAttrib", "Rf_length", "Rf_mkChar", "Rf_mkCharLen",
    "Rf_protect", "Rf_setAttrib", "Rf_unprotect",
    "Rf_ScalarComplex", "Rf_ScalarInteger", "Rf_ScalarLogical", "Rf_ScalarReal",
)


def _weak_pragmas() -> list[str]:
    """`#pragma weak` for every R symbol -- see `_WEAK_R_SYMBOLS`. Emitted at the top of a
    translation unit (before any use) so the references resolve weakly."""
    lines = ["// weak so the one library still loads into a non-R host (Python/ctypes);",
             "// see helper/codegen/emit/c_call.py"]
    lines += [f"#pragma weak {symbol}" for symbol in _WEAK_R_SYMBOLS]
    return lines

#: A Fortran numeric literal's kind suffix -- `0_int32`. C has no such thing.
_FORTRAN_KIND_SUFFIX = re.compile(r"\b(\d+(?:\.\d+)?)_[A-Za-z]\w*")
#: `size(arr)` / `size(arr, dim)` -- an extent stated in terms of another argument.
_FORTRAN_SIZE = re.compile(r"\bsize\s*\(\s*([A-Za-z_]\w*)\s*(?:,\s*(\d+)\s*)?\)")
#: `max`/`min` become the integer helpers in `tox_marshal.h`.
_FORTRAN_MINMAX = re.compile(r"\b(max|min)\s*\(")


def _size_as_c(match: re.Match) -> str:
    array, dim = match.group(1), match.group(2)
    if dim is None:
        return f"(int) Rf_length({array})"
    # Fortran dimensions are 1-based, and a matrix carries its extents in `dim`
    return f"INTEGER(Rf_getAttrib({array}, R_DimSymbol))[{int(dim) - 1}]"


def c_extent(extent: str) -> str:
    """A Fortran extent expression as C (kind suffixes, `size(...)`, and `max`/`min`)."""
    extent = _FORTRAN_KIND_SUFFIX.sub(r"\1", extent)
    extent = _FORTRAN_SIZE.sub(_size_as_c, extent)
    return _FORTRAN_MINMAX.sub(r"tox_i\1(", extent)


#: iso_c_binding kind -> the C type it interoperates as, for the extern declaration.
#: c_bool maps to unsigned char (a 1-byte pointer, ABI-identical to _Bool*) so nothing here
#: needs <stdbool.h>; the marshalling reads the byte directly.
C_CTYPE = {
    "c_int": "int",
    "c_int8_t": "int8_t",
    "c_int16_t": "int16_t",
    "c_int32_t": "int32_t",
    "c_int64_t": "int64_t",
    "c_size_t": "size_t",
    "c_float": "float",
    "c_double": "double",
    "c_float_complex": "float _Complex",
    "c_double_complex": "double _Complex",
    "c_bool": "unsigned char",
    "c_char": "char",
}

#: Fortran base type -> the R SEXP allocator type, for a plain (unconverted) array output.
#: LOGICAL and CHARACTER never appear here: they reach R through a tox_ buffer.
R_SEXPTYPE = {
    BaseType.INTEGER: "INTSXP",
    BaseType.REAL: "REALSXP",
    BaseType.COMPLEX: "CPLXSXP",
}

#: Fortran base type -> the `Rf_asXxx` coercion for a by-value scalar input.
_R_AS_SCALAR = {
    BaseType.INTEGER: "Rf_asInteger",
    BaseType.REAL: "Rf_asReal",
}


def c_ctype(argument: CArgument) -> str:
    return C_CTYPE[argument.type.kind]


def _product(extents) -> str:
    """The product of some extents, each parenthesised.

    `*` binds tighter than `-`, so an extent expression like `n - 1` must be parenthesised
    before joining; a raw join would compute a different size and overrun the buffer.
    """
    parts = [e if e.isidentifier() or e.isdigit() else f"({e})" for e in extents]
    return " * ".join(parts) if parts else "1"


class CCallEmitter:
    def __init__(self, marshal_header: str = "tox_marshal.h", dll_name: str = R_DLL_NAME):
        self.marshal_header = marshal_header
        self.dll_name = dll_name

    def marshal_header_content(self) -> str:
        # the pragmas go here (before the helpers, which use some of these symbols); the .c
        # that includes this header sees them before its own uses too
        return _MARSHAL_HEADER.replace("// __WEAK_PRAGMAS__", "\n".join(_weak_pragmas()))

    # -- module -----------------------------------------------------------------

    def module(self, module: CWrapperModule) -> str:
        writer = Writer()
        writer.line("// Generated. Do not edit.")
        writer.line(_GUARD_OPEN)
        writer.line("#include <R.h>")
        writer.line("#include <Rinternals.h>")
        writer.line(f'#include "{self.marshal_header}"')
        writer.blank()

        writer.line("// the Fortran C-ABI symbols this module calls")
        for wrapper in module:
            writer.line(self._extern_declaration(wrapper))
        writer.blank()

        for wrapper in module:
            writer.block(self.function(wrapper))
            writer.blank(collapse=False)

        writer.line(_GUARD_CLOSE)
        return writer.render(trailing_newline=True)

    def _extern_declaration(self, wrapper: CWrapper) -> str:
        params = ", ".join(self._extern_param(a) for a in wrapper)
        return f"void {wrapper.name}({params});"

    def _extern_param(self, argument: CArgument) -> str:
        const = "const " if argument.intent is Intent.IN else ""
        return f"{const}{c_ctype(argument)}*"

    def registration(self, modules) -> str:
        """`init.c`: register every `.Call` entry across all modules, once per library."""
        entries = [w for module in modules for w in module]
        writer = Writer()
        writer.line("// Generated. Do not edit.")
        writer.line(_GUARD_OPEN)
        writer.line("#include <R.h>")
        writer.line("#include <Rinternals.h>")
        writer.line("#include <R_ext/Rdynload.h>")
        for line in _weak_pragmas():
            writer.line(line)
        writer.blank()
        writer.line("// forward declarations of the .Call entry points")
        for wrapper in entries:
            n = len(list(self._inputs(wrapper)))
            names = ", ".join(["SEXP"] * n) if n else "void"
            writer.line(f"SEXP {wrapper.stripped_name}_call({names});")
        writer.blank()
        writer.line("static const R_CallMethodDef CallEntries[] = {")
        with writer.indent():
            for wrapper in entries:
                n = len(list(self._inputs(wrapper)))
                name = f"{wrapper.stripped_name}_call"
                writer.line(f'{{"{name}", (DL_FUNC) &{name}, {n}}},')
            writer.line("{NULL, NULL, 0}")
        writer.line("};")
        writer.blank()
        writer.line(f"void R_init_{self.dll_name}(DllInfo *dll) {{")
        with writer.indent():
            writer.line("R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);")
            writer.line("R_useDynamicSymbols(dll, FALSE);")
        writer.line("}")
        writer.line(_GUARD_CLOSE)
        return writer.render(trailing_newline=True)

    # -- function ---------------------------------------------------------------

    def function(self, wrapper: CWrapper) -> str:
        writer = Writer()
        # scalars pulled into <name>_v locals; extent expressions naming one must read the
        # local, not the SEXP parameter. Held for the duration of emitting this function.
        self._scalar_names = self._scalar_local_names(wrapper)
        params = ", ".join(f"SEXP {a.name}" for a in self._inputs(wrapper))
        if not params:
            params = "void"
        writer.line(f"SEXP {wrapper.stripped_name}_call({params}) {{")

        with writer.indent():
            writer.line("int nprot = 0;")
            self._materialize_optionals(writer, wrapper)
            self._derive(writer, wrapper)
            self._extract_scalars(writer, wrapper)
            self._copy_inout(writer, wrapper)
            self._marshal_inputs(writer, wrapper)
            self._allocate(writer, wrapper)
            self._call(writer, wrapper)
            self._marshal_outputs(writer, wrapper)
            self._return(writer, wrapper)

        writer.line("}")
        return writer.render()

    # -- argument classification (mirrors r_wrapper, so R passes what C declares) ----------

    def _inputs(self, wrapper: CWrapper) -> list[CArgument]:
        """The arguments the R wrapper passes in -- the same set the Python emitter asks for."""
        return [
            argument
            for argument in wrapper
            if argument.intent.is_input
            and not argument.is_temporary
            and not self._is_derived(argument)
            and (not argument.is_synthesised or self._must_be_supplied(argument, wrapper))
        ]

    def _must_be_supplied(self, argument: CArgument, wrapper: CWrapper) -> bool:
        """Whether a synthesised extent or strlen has to come from the caller after all."""
        if argument.is_error:
            return False
        return self._derived_value(argument, wrapper) is None

    @staticmethod
    def _is_derived(argument: CArgument) -> bool:
        """Whether C works this out itself, so it is not a parameter of the .Call entry."""
        roles = argument.source.roles if argument.source else None
        if roles is None:
            return False
        if roles.is_computed:
            return False
        return (roles.is_inferable_extent or roles.is_inferable_shape_arg
                or roles.is_mask_count)

    # -- deriving extents, shapes, counts ---------------------------------------

    def _derive(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()

        def is_shape(argument):
            roles = argument.source.roles if argument.source else None
            return bool(roles and roles.is_shape_arg)

        ordered = sorted(wrapper, key=lambda argument: not is_shape(argument))
        for argument in ordered:
            expression = self._derived_value(argument, wrapper)
            if expression is None:
                continue
            if is_shape(argument):
                # a shape is a fresh SEXP (the dim attribute or a length-1 vector); it must
                # be protected while it sizes outputs and is set back as their dim
                lines.line(f"SEXP {argument.name} = PROTECT({expression}); nprot++;")
            else:
                lines.line(f"int {argument.name} = {expression};")
        if lines:
            writer.line("// derived from the inputs, not asked of the caller")
            writer.extend(lines)
            writer.blank()

    def _derived_value(self, argument: CArgument, wrapper: CWrapper) -> str | None:
        source = argument.source
        roles = source.roles if source else None

        if argument.origin is Origin.STRLEN:
            owner = wrapper.argument(argument.sizes)
            if owner is None or not owner.intent.is_input:
                return None
            # the longest string: character(len=n) is one width for the whole array
            return f"tox_max_strlen({argument.sizes})"

        if argument.origin is Origin.EXTENT:
            owner = wrapper.argument(argument.sizes)
            if owner is None or not owner.intent.is_input or owner.is_temporary:
                return None
            return self._extent_of(owner, argument.axis)

        if roles is None:
            return None

        if roles.is_mask_count:
            return f"tox_sum_true({roles.mask_count_of.name})"

        if roles.is_shape_arg:
            if not roles.shape_of.intent.is_input:
                return None
            # tox_shape_of yields the dim attribute, or [length] for a plain vector
            return f"tox_shape_of({roles.shape_of.name})"

        if roles.is_extent:
            for owner in roles.extent_of:
                c_owner = wrapper.argument(owner.name)
                if c_owner is None or c_owner.shape_arg is None:
                    continue
                # its extents travel separately, so an extent of it is the total count
                if c_owner.intent.is_input:
                    return f"(int) Rf_length({owner.name})"
                return f"tox_prod({c_owner.shape_arg})"
            best = None
            for owner in roles.extent_of:
                c_owner = wrapper.argument(owner.name)
                if c_owner is None or not c_owner.intent.is_input:
                    continue
                if c_owner.is_temporary:
                    continue
                axis = self._axis_of(argument.name, c_owner)
                if axis is None:
                    continue
                if not c_owner.optional:
                    return self._extent_of(c_owner, axis)
                best = best or self._extent_of(c_owner, axis)
            return best

        return None

    def _extent_of(self, owner: CArgument, axis: int) -> str:
        """The size of `owner` along `axis`, as a C expression."""
        if owner.optional:
            if owner.rank > 1 and not owner.type.is_character:
                return f"{owner.name}_dim[{axis}]"
            return f"{owner.name}_size"
        if owner.rank <= 1 or (owner.type.is_character and owner.rank == 2):
            return f"(int) Rf_length({owner.name})"
        return f"INTEGER(Rf_getAttrib({owner.name}, R_DimSymbol))[{axis}]"

    @staticmethod
    def _axis_of(extent: str, owner: CArgument) -> int | None:
        extents = list(owner.dimension.extents)
        if owner.type.is_character and extents:
            extents = extents[1:]
        try:
            return extents.index(extent)
        except ValueError:
            return None

    # -- scalars ----------------------------------------------------------------

    def _scalar_local_names(self, wrapper: CWrapper) -> set[str]:
        """Scalar inputs pulled into a `<name>_v` local before the call.

        Exactly the by-value scalars among the `.Call` parameters (`_inputs`): a derived
        extent/strlen/shape-count is an `int <name>` local already, a scalar output or
        `ierr` is a plain C local, an optional goes through `<name>_p`, a mode/character
        goes through a buffer. Membership in `_inputs` is the reliable discriminator --
        `_is_derived` alone misses STRLEN and shape-count derivations.
        """
        return {
            a.name for a in self._inputs(wrapper)
            if a.is_scalar and not a.optional and not self._converts_via_buffer(a)
            and a.conversion is not Conversion.MODE
        }

    def _scalar_inputs(self, wrapper: CWrapper) -> list[CArgument]:
        names = self._scalar_local_names(wrapper)
        return [a for a in self._inputs(wrapper) if a.name in names]

    def _extent(self, expr: str) -> str:
        """A Fortran extent as C, with any by-value scalar name read from its `_v` local."""
        out = c_extent(expr)
        for name in self._scalar_names:
            out = re.sub(rf"\b{re.escape(name)}\b", f"{name}_v", out)
        return out

    def _extract_scalars(self, writer: Writer, wrapper: CWrapper) -> None:
        """Pull a by-value scalar out of its length-1 SEXP into a C local to point at."""
        lines = Writer()
        for argument in self._scalar_inputs(wrapper):
            lines.line(self._scalar_local(argument))
        if lines:
            writer.line("// scalar inputs, pulled from their length-1 vectors")
            writer.extend(lines)
            writer.blank()

    def _scalar_local(self, argument: CArgument) -> str:
        return (
            f"{self._scalar_local_type(argument)} {argument.name}_v"
            f" = {self._scalar_value(argument)};"
        )

    @staticmethod
    def _scalar_local_type(argument: CArgument) -> str:
        """The C type a by-value scalar is pulled into, which conversion may narrow."""
        base = argument.type.base
        if base is BaseType.LOGICAL:
            return "unsigned char"
        if base is BaseType.COMPLEX:
            return "Rcomplex"
        return c_ctype(argument)

    @staticmethod
    def _scalar_value(argument: CArgument) -> str:
        """The expression reading a by-value scalar out of its length-1 SEXP."""
        name = argument.name
        base = argument.type.base
        if base is BaseType.LOGICAL:
            return f"(Rf_asLogical({name}) == TRUE) ? 1 : 0"
        if base is BaseType.COMPLEX:
            return f"COMPLEX({name})[0]"
        cast = "" if base is BaseType.REAL or base is BaseType.INTEGER else f"({c_ctype(argument)}) "
        return f"{cast}{_R_AS_SCALAR[base]}({name})"

    # -- optionals --------------------------------------------------------------

    def _materialize_optionals(self, writer: Writer, wrapper: CWrapper) -> None:
        """Absent optional -> NULL pointer (Fortran sees it absent) and size 0."""
        lines = Writer()
        for argument in wrapper:
            if not argument.optional:
                continue
            name = argument.name
            ctype = c_ctype(argument)
            rank = argument.rank
            has_axes = rank > 1 and not argument.type.is_character
            # A scalar that needs converting -- a logical, which R stores four bytes wide --
            # has no storage of its own to point at, and the buffer path only covers arrays
            # and characters. It gets a local to convert into and a pointer to that local,
            # left null while the argument is absent.
            converts_in_place = (
                argument.is_scalar
                and argument.needs_conversion
                and not self._converts_via_buffer(argument)
            )
            if converts_in_place:
                local_type = self._scalar_local_type(argument)
                lines.line(f"{local_type} {name}_v = 0;")
                lines.line(f"const {local_type}* {name}_p = NULL;")
            elif not argument.needs_conversion:
                lines.line(f"const {ctype}* {name}_p = NULL;")
            lines.line(f"int {name}_size = 0;")
            if has_axes:
                lines.line(f"int {name}_dim[{rank}]; "
                           f"for (int i = 0; i < {rank}; ++i) {name}_dim[i] = 0;")
            lines.line(f"if ({name} != R_NilValue) {{")
            with lines.indent():
                lines.line(f"{name}_size = (int) Rf_length({name});")
                if converts_in_place:
                    lines.line(f"{name}_v = {self._scalar_value(argument)};")
                    lines.line(f"{name}_p = &{name}_v;")
                elif not argument.needs_conversion:
                    lines.line(f"{name}_p = {self._c_pointer(argument, name)};")
                if has_axes:
                    lines.line(f"SEXP {name}_d = Rf_getAttrib({name}, R_DimSymbol);")
                    lines.line(f"if ({name}_d != R_NilValue) {{")
                    with lines.indent():
                        lines.line(
                            f"for (int i = 0; i < {rank} && i < (int) Rf_length({name}_d); ++i)"
                            f" {name}_dim[i] = INTEGER({name}_d)[i];"
                        )
                    lines.line("}")
            lines.line("}")
        if lines:
            writer.line("// optionals: a null pointer and size 0 when the caller omits them")
            writer.extend(lines)
            writer.blank()

    # -- shared classification helpers -------------------------------------------

    @staticmethod
    def _converts_via_buffer(argument: CArgument) -> bool:
        if not argument.needs_conversion:
            return False
        return argument.is_array or argument.type.is_character

    @staticmethod
    def _is_inout_copy(argument: CArgument) -> bool:
        return (
            argument.intent is Intent.INOUT
            and argument.is_array
            and not argument.needs_conversion
            and not argument.is_temporary
        )

    def _working_name(self, argument: CArgument) -> str:
        return f"{argument.name}_out" if self._is_inout_copy(argument) else argument.name

    def _copy_inout(self, writer: Writer, wrapper: CWrapper) -> None:
        """A plain array modified in place is duplicated, so the caller's stays intact."""
        lines = Writer()
        for argument in wrapper:
            if self._is_inout_copy(argument):
                lines.line(f"SEXP {argument.name}_out = PROTECT(Rf_duplicate({argument.name})); nprot++;")
        if lines:
            writer.line("// copy what is modified in place, so the caller's stays intact")
            writer.extend(lines)
            writer.blank()

    def _marshal_inputs(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            if not argument.intent.is_input or not self._converts_via_buffer(argument):
                continue
            lines.line(self._input_buffer(argument))
        if lines:
            writer.line("// convert what Fortran cannot take from R directly")
            writer.extend(lines)
            writer.blank()

    def _input_buffer(self, argument: CArgument) -> str:
        name = argument.name
        if argument.conversion is Conversion.LOGICAL:
            return f"unsigned char* {name}_c = tox_bool_in({name});"
        length = self._extent(argument.dimension.extents[0])
        return f"char* {name}_c = tox_char_in({name}, {length});"

    # -- allocating outputs and work space --------------------------------------

    def _allocate(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            if argument.intent.is_input and not argument.is_temporary:
                continue
            if self._is_derived(argument) and argument.intent.is_input:
                continue
            for line in self._new_value(argument):
                lines.line(line)
        if lines:
            writer.line("// outputs and work space")
            writer.extend(lines)
            writer.blank()

    def _new_value(self, argument: CArgument) -> list[str]:
        name = argument.name
        if argument.origin is Origin.ERROR:
            return [f"int {name} = 0;"]
        if argument.is_scalar:
            if argument.type.base is BaseType.LOGICAL:
                return [f"unsigned char {name} = 0;"]
            if argument.type.base is BaseType.COMPLEX:
                return [f"Rcomplex {name} = {{0.0, 0.0}};"]
            return [f"{c_ctype(argument)} {name} = 0;"]

        extents = [self._extent(e) for e in argument.dimension.extents]
        if argument.shape_arg is not None:
            size = f"tox_prod({argument.shape_arg})"
        else:
            size = _product(extents)

        if argument.conversion is Conversion.LOGICAL:
            return [f"unsigned char* {name}_c = tox_bool_alloc({size});"]
        if argument.type.is_character:
            length = extents[0]
            count = size if argument.shape_arg is not None else _product(extents[1:])
            return [f"char* {name}_c = tox_char_alloc({length}, {count});"]
        if argument.is_temporary:
            ctype = c_ctype(argument)
            return [f"{ctype}* {name} = ({ctype}*) R_alloc({size}, sizeof({ctype}));"]
        # a plain output: an R vector, protected, with a dim attribute if rank >= 2
        out = [f"SEXP {name} = PROTECT(Rf_allocVector({R_SEXPTYPE[argument.type.base]}, {size})); nprot++;"]
        dim = self._dim_lines(argument)
        return out + dim

    def _dim_lines(self, argument: CArgument) -> list[str]:
        """Give a rank-2-or-more output the `dim` that makes it a matrix in R."""
        target = self._working_name(argument)
        if argument.shape_arg is not None:
            return [f"Rf_setAttrib({target}, R_DimSymbol, {argument.shape_arg});"]
        extents = self._r_extents(argument)
        if len(extents) < 2:
            return []
        rank = len(extents)
        assigns = " ".join(f"INTEGER({target}_dim)[{i}] = {e};" for i, e in enumerate(extents))
        return [
            f"{{ SEXP {target}_dim = PROTECT(Rf_allocVector(INTSXP, {rank})); "
            f"{assigns} Rf_setAttrib({target}, R_DimSymbol, {target}_dim); UNPROTECT(1); }}"
        ]

    def _r_extents(self, argument: CArgument) -> list[str]:
        extents = [self._extent(e) for e in argument.dimension.extents]
        if argument.type.is_character and extents:
            extents = extents[1:]
        return extents

    # -- the call ---------------------------------------------------------------

    def _call(self, writer: Writer, wrapper: CWrapper) -> None:
        writer.line(f"{wrapper.name}(")
        with writer.indent():
            actuals = [self._actual(a) for a in wrapper]
            for index, actual in enumerate(actuals):
                comma = "" if index == len(actuals) - 1 else ","
                writer.line(f"{actual}{comma}")
        writer.line(");")
        writer.blank()

    def _actual(self, argument: CArgument) -> str:
        name = argument.name
        if argument.optional:
            if self._converts_via_buffer(argument):
                return f"{name} != R_NilValue ? {name}_c : NULL"
            return f"{name}_p"
        if self._converts_via_buffer(argument):
            return f"{name}_c"
        if argument.is_scalar:
            # a by-value scalar input lives in <name>_v; a derived extent, scalar output or
            # ierr is a plain <name> local
            ref = f"{name}_v" if name in self._scalar_names else name
            if argument.type.base is BaseType.COMPLEX:
                const = "const " if argument.intent is Intent.IN else ""
                return f"({const}double _Complex*) &{ref}"
            return f"&{ref}"
        if argument.is_temporary:
            return name
        return self._c_pointer(argument, self._working_name(argument))

    def _c_pointer(self, argument: CArgument, sexp: str) -> str:
        """The R vector `sexp` as the pointer the extern declaration asks for."""
        base = argument.type.base
        if base is BaseType.REAL:
            return f"REAL({sexp})"
        if base is BaseType.INTEGER:
            return f"INTEGER({sexp})"
        if base is BaseType.COMPLEX:
            const = "const " if argument.intent is Intent.IN else ""
            return f"({const}double _Complex*) COMPLEX({sexp})"
        # logical/character reach C through a buffer, never this path
        return f"REAL({sexp})"

    # -- marshalling outputs back and returning ---------------------------------

    def _marshal_outputs(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            if argument.intent is Intent.IN or not self._converts_via_buffer(argument):
                continue
            name = argument.name
            count = self._converted_count(argument)
            if argument.conversion is Conversion.LOGICAL:
                lines.line(f"SEXP {name} = PROTECT(tox_bool_out({name}_c, {count})); nprot++;")
            else:
                length = self._extent(argument.dimension.extents[0])
                lines.line(f"SEXP {name} = PROTECT(tox_char_out({name}_c, {length}, {count})); nprot++;")
            for line in self._dim_lines(argument):
                lines.line(line)
        if lines:
            writer.line("// convert the outputs back")
            writer.extend(lines)
            writer.blank()

    def _converted_count(self, argument: CArgument) -> str:
        extents = [self._extent(e) for e in argument.dimension.extents]
        if argument.shape_arg is not None:
            return f"tox_prod({argument.shape_arg})"
        if argument.type.is_character:
            return _product(extents[1:])
        return _product(extents)

    def _return(self, writer: Writer, wrapper: CWrapper) -> None:
        outputs = [a for a in wrapper if self._is_returned(a)]
        error = wrapper.error_argument
        entries = [(a.name, self._result_sexp(a)) for a in outputs]
        entries.append((error.name, f"Rf_ScalarInteger({error.name})"))
        n = len(entries)

        writer.line(f"SEXP _out = PROTECT(Rf_allocVector(VECSXP, {n})); nprot++;")
        for index, (_, value) in enumerate(entries):
            writer.line(f"SET_VECTOR_ELT(_out, {index}, {value});")
        writer.line(f"SEXP _nms = PROTECT(Rf_allocVector(STRSXP, {n})); nprot++;")
        for index, (nm, _) in enumerate(entries):
            writer.line(f'SET_STRING_ELT(_nms, {index}, Rf_mkChar("{nm}"));')
        writer.line("Rf_setAttrib(_out, R_NamesSymbol, _nms);")
        writer.line("UNPROTECT(nprot);")
        writer.line("return _out;")

    def _result_sexp(self, argument: CArgument) -> str:
        """The SEXP for one returned output, wrapping a scalar or naming an allocated one."""
        if argument.is_scalar:
            base = argument.type.base
            # an inout scalar was pulled into <name>_v and written there; a pure output
            # scalar is a plain <name> local
            name = f"{argument.name}_v" if argument.name in self._scalar_names else argument.name
            if base is BaseType.LOGICAL:
                return f"Rf_ScalarLogical({name} != 0)"
            if base is BaseType.REAL:
                return f"Rf_ScalarReal({name})"
            if base is BaseType.COMPLEX:
                return f"Rf_ScalarComplex({name})"
            return f"Rf_ScalarInteger({name})"
        return self._working_name(argument)

    @staticmethod
    def _is_returned(argument: CArgument) -> bool:
        """Every output the R wrapper might use, ierr aside (added separately)."""
        return (
            argument.intent.is_output
            and not argument.is_error
            and not argument.is_temporary
        )


_MARSHAL_HEADER = r'''// Generated. Do not edit.
//
// Marshalling helpers for the R binding. Pure C, R's C API only. These convert what
// Fortran cannot take from R directly -- R's int-based logicals and its strings -- plus a
// couple of small shape helpers. Transient buffers come from R_alloc and are freed when the
// .Call returns, so there is nothing to free by hand and nothing to leak on an error
// longjmp. The R layer has already validated and rejected NA, so these are straight copies.
#ifndef TOX_MARSHAL_H
#define TOX_MARSHAL_H

#include <R.h>
#include <Rinternals.h>
#include <string.h>

// __WEAK_PRAGMAS__

static inline int tox_imax(int a, int b) { return a > b ? a : b; }
static inline int tox_imin(int a, int b) { return a < b ? a : b; }

// The width a character(len=n) array needs: the longest element, NA skipped. 0 for a
// non-character or absent (R_NilValue) argument, so an omitted optional reports no width.
static inline int tox_max_strlen(SEXP x) {
    if (x == R_NilValue || TYPEOF(x) != STRSXP) return 0;
    int longest = 0, n = (int) XLENGTH(x);
    for (int i = 0; i < n; ++i) {
        SEXP e = STRING_ELT(x, i);
        if (e == NA_STRING) continue;
        int m = (int) LENGTH(e);
        if (m > longest) longest = m;
    }
    return longest;
}

// Fortran carries a string's length as the leading extent: n strings of length len are
// char[len * n], column-major, each zero-padded. Reading stops at the first null, so a
// zero-filled buffer yields empty strings. Anything longer than len is truncated.
static inline char* tox_char_in(SEXP x, int len) {
    int n = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char* buf = (char*) R_alloc(total, 1);
    memset(buf, 0, total);
    for (int i = 0; i < n; ++i) {
        SEXP e = STRING_ELT(x, i);
        if (e == NA_STRING) continue;
        int slen = (int) LENGTH(e);
        int m = slen < len ? slen : len;
        memcpy(buf + (size_t) i * len, CHAR(e), (size_t) m);
    }
    return buf;
}

// A zero-filled c_char output buffer of len * n bytes.
static inline char* tox_char_alloc(int len, int n) {
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char* buf = (char*) R_alloc(total, 1);
    memset(buf, 0, total);
    return buf;
}

// char[len * n] fixed-width -> STRSXP, each slot read up to the first null. Returned
// unprotected: the caller protects it straight into a result slot.
static inline SEXP tox_char_out(const char* buf, int len, int n) {
    SEXP out = Rf_allocVector(STRSXP, n);
    for (int i = 0; i < n; ++i) {
        const char* p = buf + (size_t) i * len;
        int m = 0;
        while (m < len && p[m] != '\0') ++m;
        SET_STRING_ELT(out, i, Rf_mkCharLen(p, m));
    }
    return out;
}

// R logical is a 4-byte int (0/1); Fortran c_bool is one byte. Copy across.
static inline unsigned char* tox_bool_in(SEXP x) {
    int n = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    unsigned char* buf = (unsigned char*) R_alloc(n > 0 ? n : 1, 1);
    if (n > 0) {
        const int* px = LOGICAL(x);
        for (int i = 0; i < n; ++i) buf[i] = (px[i] == TRUE) ? 1 : 0;
    }
    return buf;
}

// A zero-filled c_bool output buffer of n bytes.
static inline unsigned char* tox_bool_alloc(int n) {
    unsigned char* buf = (unsigned char*) R_alloc(n > 0 ? n : 1, 1);
    memset(buf, 0, (size_t) (n > 0 ? n : 1));
    return buf;
}

// c_bool byte buffer -> LGLSXP. Read the raw byte and test != 0: ifx writes 0xFF for true,
// which must map to R's TRUE (1), not stay 255. Returned unprotected.
static inline SEXP tox_bool_out(const unsigned char* buf, int n) {
    SEXP out = Rf_allocVector(LGLSXP, n);
    int* po = LOGICAL(out);
    for (int i = 0; i < n; ++i) po[i] = buf[i] != 0 ? TRUE : FALSE;
    return out;
}

// The dim attribute of x, or a length-1 vector holding its length when it has none. A plain
// R vector carries no dim, and its shape is its length. Returned unprotected.
static inline SEXP tox_shape_of(SEXP x) {
    SEXP dim = Rf_getAttrib(x, R_DimSymbol);
    if (dim != R_NilValue) return dim;
    SEXP s = Rf_allocVector(INTSXP, 1);
    INTEGER(s)[0] = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    return s;
}

// The product of an integer shape vector: the flat element count of the array it describes.
static inline int tox_prod(SEXP shape) {
    int p = 1, n = (int) XLENGTH(shape);
    const int* ps = INTEGER(shape);
    for (int i = 0; i < n; ++i) p *= ps[i];
    return p;
}

// Count the TRUEs in a logical vector; NA has already been rejected in R.
static inline int tox_sum_true(SEXP x) {
    int n = (int) XLENGTH(x), total = 0;
    const int* px = LOGICAL(x);
    for (int i = 0; i < n; ++i) total += (px[i] == TRUE) ? 1 : 0;
    return total;
}

#endif
'''
