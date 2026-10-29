"""Emitting the C half of the R binding.

The split with R is documented in `design/language-layers.md`: **R checks what an argument
is, C checks every element while it converts it, and both raise the same conditions.** By
the time a value reaches C, R has checked its kind and rank, coerced a double, and checked a
single value in full; C converts what Fortran cannot take from R directly (R's int-based
logicals, its strings, an integer array that arrives as doubles), checking each element in
the pass it makes anyway, and hands everything to the `extern "C"` Fortran wrapper. A
refusal is raised from C as the classed condition R raises itself (`tox_refuse`).

Each generated function is a `.Call` entry `<name>_call`, registered in `init.c`. It is
internal: the user-facing `<name>` is the R wrapper on top of it, which calls
`.Call("<name>_call", ...)` and decodes the returned list of every output plus `ierr`.

Only R's C API is used -- no Rcpp, no cpp11. The conversions C cannot avoid (logical and
integer arrays, strings) go through `tox_marshal.h`, whose buffers come from `R_alloc` and so
are freed automatically when the call returns; there is nothing to leak on an error longjmp,
which a refusal is.

Strings cross fixed-width, column-major and **blank**-padded, in both directions. The
Fortran wrapper does not convert that buffer -- it takes a `character(len=strlen)` pointer
view of it -- so the padding byte is part of the value, and blanks are what Fortran's
`trim`, `len_trim` and string comparison are defined on. See `design/c-layer.md`.
"""

from __future__ import annotations

import hashlib

from ..abi import extents
from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..ir.types import BaseType, Intent
from ..names import (
    BUFFER,
    C_LOOP_INDEX,
    C_PROTECT_COUNT,
    DIM,
    INOUT_COPY,
    OPTIONAL_DIM_ATTRIBUTE,
    OPTIONAL_POINTER,
    OPTIONAL_SIZE,
    SCALAR_VALUE,
    shim_local,
)
from ..render import Writer
from .errors_r import ARGUMENT_CONDITION_CLASSES, CONDITION_FIELDS, NA_MESSAGE

#: The R shared object this binding is registered into.
#:
#: R derives the initialisation routine it calls on `dyn.load` from the file's own basename,
#: **including the `lib` prefix and with no substitution at all**: `libtensor_omics.so` is
#: entered through `R_init_libtensor_omics`, and nothing else. That is why the library is not
#: called `libtensor-omics.so` any more -- `R_init_libtensor-omics` is not a C identifier, and
#: the one way to emit it (a quoted `__asm__` label) is GNU-as-only: clang and icx put the
#: quote characters *in the symbol name*, silently, and fpm compiles the shims with icx
#: whenever the Fortran compiler is ifx. A name R can reach is the portable fix.
#:
#: Keep this equal to the basename of the artifact `build.sh` produces, which is `lib` plus
#: fpm's project name. `test_end_to_end_r` checks the routines really do get registered.
R_DLL_NAME = "libtensor_omics"

#: The R shims compile into the one `libtensor_omics.so`. This guard drops them (and their
#: R.h include) when the binding is built without R (`NO_R_BINDING`) or without the C ABI
#: they call (`NO_C_BINDING`), leaving empty objects that need no R headers.
_GUARD_OPEN = "#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)"
_GUARD_CLOSE = "#endif  // R binding"

#: Every R API symbol the shims reference. R provides them at dyn.load, but the one
#: libtensor_omics.so also loads into a non-R host (Python via ctypes). Marking them all weak
#: lets that load succeed under eager binding -- the undefined R symbols resolve to null, and
#: the R-only code that would use them never runs from Python -- while a *genuinely* missing
#: symbol (a build regression) still fails loudly at load. R binds them to its real
#: definitions at dyn.load. The `test_ctypes_loads_the_r_bundled_library` end-to-end test
#: guards this list: if the emitter starts using an R symbol not here, that load fails.
#: (NA_STRING is R_NaString; NA_INTEGER and NA_LOGICAL are R_NaInt; ISNA is R_IsNA;
#: PROTECT/UNPROTECT are Rf_protect/Rf_unprotect; CHAR is R_CHAR.)
_WEAK_R_SYMBOLS = (
    "COMPLEX", "INTEGER", "INTEGER_NO_NA", "LENGTH", "LOGICAL", "REAL", "STRING_ELT", "TYPEOF",
    "XLENGTH", "SET_STRING_ELT", "SET_VECTOR_ELT",
    "R_BaseEnv", "R_CHAR", "STRING_PTR_RO", "R_ClassSymbol", "R_IsNA", "R_NaInt", "R_NegInf", "R_PosInf",
    "R_alloc", "R_DimSymbol", "R_NamesSymbol", "R_NilValue", "R_NaString",
    "R_registerRoutines", "R_useDynamicSymbols",
    "Rf_allocVector", "Rf_asInteger", "Rf_charIsASCII", "Rf_asLogical", "Rf_asReal", "Rf_coerceVector",
    "Rf_duplicate", "Rf_eval", "Rf_getAttrib", "Rf_getCharCE", "Rf_install", "Rf_lang2",
    "Rf_length", "Rf_mkChar", "Rf_mkCharLen", "Rf_mkCharLenCE", "Rf_mkString",
    "Rf_protect", "Rf_setAttrib", "Rf_translateCharUTF8", "Rf_type2char", "Rf_unprotect",
    "Rf_ScalarComplex", "Rf_ScalarInteger", "Rf_ScalarLogical", "Rf_ScalarReal",
    "Rf_ScalarString",
)


def _weak_pragmas() -> list[str]:
    """`#pragma weak` for every R symbol -- see `_WEAK_R_SYMBOLS`. Emitted at the top of a
    translation unit (before any use) so the references resolve weakly."""
    lines = ["// weak so the one library still loads into a non-R host (Python/ctypes);",
             "// see helper/codegen/emit/c_call.py"]
    lines += [f"#pragma weak {symbol}" for symbol in _WEAK_R_SYMBOLS]
    return lines

def c_extent(extent: str, name=lambda n: n, owner_rank=lambda array: 2) -> str:
    """A Fortran extent expression as C, through the one grammar the bindings evaluate.

    `name` writes an identifier (a by-value scalar is read from its `_v_` local), and
    `owner_rank` gives an array's rank as R sees it, which decides how `size(a, k)` is read:
    a vector carries no `dim`, so its only extent is its length. `generate` has refused every
    extent a binding computes that lies outside the grammar (`abi.extents`).
    """
    def size(array: str, dimension: int | None) -> str:
        if dimension is None or owner_rank(array) <= 1:
            return f"(int) Rf_length({array})"
        # Fortran dimensions are 1-based, and a matrix carries its extents in `dim`
        return f"INTEGER(Rf_getAttrib({array}, R_DimSymbol))[{dimension - 1}]"

    return extents.to_c(extent, name, size)


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


def _r_marshals(argument: CArgument) -> bool:
    """Whether R must build the bytes rather than point into the SEXP it was handed.

    `CArgument.needs_conversion` answers this for the *Fortran* side, and the two questions
    are not the same. It is false for a `logical(c_bool)` dummy, because Fortran takes the C
    buffer as it stands -- but R never can: its logical vector is `LGLSXP`, an array of
    4-byte `int`, whatever kind the Fortran declares. The answers coincided only for as long
    as every logical was default-kind, and three sites had quietly keyed on the wrong one.
    And an integer array the caller passed is handed over as R holds it, an integer or a
    double vector (`_checks_integer_array`).
    """
    return argument.type.is_logical or argument.needs_conversion or _checks_integer_array(argument)


def _checks_integer_array(argument: CArgument) -> bool:
    """Whether this is an integer array the caller passed, which the shim checks as it converts.

    The R wrapper hands it over unconverted (`tox_validate.R`): a double vector of whole
    numbers is an integer array to the contract, and its elements are checked here, in the
    one pass that converts them, rather than in R (`tox_int_in`). Not a single value, which R
    checks and converts itself, since it uses it; not an output or a work array; and not what
    the binding works out or a producer computes, which are integers already -- nor a shape,
    which the shim uses as an R vector, as the dim of what it allocates.
    """
    if argument.type.base is not BaseType.INTEGER or argument.intent is not Intent.IN:
        return False
    if not argument.is_array or argument.is_temporary:
        return False
    roles = argument.source.roles if argument.source else None
    return roles is None or not (roles.is_derived or roles.is_computed or roles.is_shape_arg)



class CCallEmitter:
    def __init__(self, marshal_header: str = "tox_marshal.h", dll_name: str = R_DLL_NAME):
        self.marshal_header = marshal_header
        self.dll_name = dll_name

    def marshal_header_content(self) -> str:
        # the pragmas go here (before the helpers, which use some of these symbols); the .c
        # that includes this header sees them before its own uses too. The refusal's classes,
        # fields and NA message are R's own (`emit.errors_r`), so the shim's refusals and the
        # R wrapper's cannot drift apart.
        quoted = lambda names: ", ".join(f'"{name}"' for name in names)
        return (_MARSHAL_HEADER
                .replace("// __WEAK_PRAGMAS__", "\n".join(_weak_pragmas()))
                .replace("__ARGUMENT_CLASSES__", quoted(ARGUMENT_CONDITION_CLASSES))
                .replace("__CONDITION_FIELDS__", quoted(CONDITION_FIELDS))
                .replace("__NA_MESSAGE__", NA_MESSAGE))

    def header_stamp(self) -> str:
        """A comment naming the content hash of the marshal header this file includes.

        fpm decides whether to recompile a source from a hash of *that file's own text* --
        never its includes, and never timestamps (fpm_backend.F90; fortran-lang/fpm#358,
        open since 2021). So a change to the header alone recompiled nothing, and the old
        helpers stayed linked into the library. Every file that includes the header carries
        its hash, so any change to the header is a change to each of them, and fpm rebuilds
        exactly those -- with nothing that depends on fpm's build-directory layout.

        The hash is of the content this emitter writes, so the stamp and the header on disk
        cannot disagree after a generator run; `--check` covers the rest.
        """
        digest = hashlib.sha256(self.marshal_header_content().encode()).hexdigest()[:16]
        return (f"// {self.marshal_header} {digest} -- its hash, so that fpm, which only "
                "hashes this file, recompiles it when the header changes")

    # -- module -----------------------------------------------------------------

    def module(self, module: CWrapperModule) -> str:
        # a module constant in an extent is its value here: the shim's C cannot name it
        module = module.with_constant_sizes_evaluated()
        writer = Writer()
        writer.line("// Generated. Do not edit.")
        writer.line(_GUARD_OPEN)
        writer.line("#include <R.h>")
        writer.line("#include <Rinternals.h>")
        writer.line(f'#include "{self.marshal_header}"')
        writer.line(self.header_stamp())
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
        # scalars pulled into `_v_<name>` locals; extent expressions naming one must read the
        # local, not the SEXP parameter. Held for the duration of emitting this function.
        self._scalar_names = self._scalar_local_names(wrapper)
        self._wrapper = wrapper
        params = ", ".join(f"SEXP {a.name}" for a in self._inputs(wrapper))
        if not params:
            params = "void"
        writer.line(f"SEXP {wrapper.stripped_name}_call({params}) {{")

        with writer.indent():
            writer.line(f"int {C_PROTECT_COUNT} = 0;")
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
                lines.line(f"SEXP {argument.name} = PROTECT({expression}); {C_PROTECT_COUNT}++;")
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
                return f"{shim_local(owner.name, DIM)}[{axis}]"
            return shim_local(owner.name, OPTIONAL_SIZE)
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
        """Scalar inputs pulled into a `_v_<name>` local before the call.

        Exactly the by-value scalars among the `.Call` parameters (`_inputs`): a derived
        extent/strlen/shape-count is an `int <name>` local already, a scalar output or
        `ierr` is a plain C local, an optional goes through `_p_<name>`, a mode/character
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
        """A Fortran extent as C, with any by-value scalar name read from its `_v_` local.

        An assumed size (`*`) has no value to write; an array declared so is sized by the shape
        that travels beside it, which is what every caller of this reads instead.
        """
        if expr == "*":
            return expr

        def name(identifier: str) -> str:
            if identifier in self._scalar_names:
                return shim_local(identifier, SCALAR_VALUE)
            return identifier

        def owner_rank(array: str) -> int:
            owner = self._wrapper.argument(array)
            if owner is None:
                return 2
            return owner.rank - 1 if owner.type.is_character else owner.rank

        return c_extent(expr, name, owner_rank)

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
            f"{self._scalar_local_type(argument)} {shim_local(argument.name, SCALAR_VALUE)}"
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
                and _r_marshals(argument)
                and not self._converts_via_buffer(argument)
            )
            value = shim_local(name, SCALAR_VALUE)
            pointer = shim_local(name, OPTIONAL_POINTER)
            size = shim_local(name, OPTIONAL_SIZE)
            dim = shim_local(name, DIM)
            dim_attribute = shim_local(name, OPTIONAL_DIM_ATTRIBUTE)
            i = C_LOOP_INDEX
            if converts_in_place:
                local_type = self._scalar_local_type(argument)
                lines.line(f"{local_type} {value} = 0;")
                lines.line(f"const {local_type}* {pointer} = NULL;")
            elif not _r_marshals(argument):
                lines.line(f"const {ctype}* {pointer} = NULL;")
            lines.line(f"int {size} = 0;")
            if has_axes:
                lines.line(f"int {dim}[{rank}]; "
                           f"for (int {i} = 0; {i} < {rank}; ++{i}) {dim}[{i}] = 0;")
            lines.line(f"if ({name} != R_NilValue) {{")
            with lines.indent():
                lines.line(f"{size} = (int) Rf_length({name});")
                if converts_in_place:
                    lines.line(f"{value} = {self._scalar_value(argument)};")
                    lines.line(f"{pointer} = &{value};")
                elif not _r_marshals(argument):
                    lines.line(f"{pointer} = {self._c_pointer(argument, name)};")
                if has_axes:
                    lines.line(f"SEXP {dim_attribute} = Rf_getAttrib({name}, R_DimSymbol);")
                    lines.line(f"if ({dim_attribute} != R_NilValue) {{")
                    with lines.indent():
                        lines.line(
                            f"for (int {i} = 0; {i} < {rank} && {i} < (int) Rf_length({dim_attribute}); ++{i})"
                            f" {dim}[{i}] = INTEGER({dim_attribute})[{i}];"
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
        """Whether R marshals this argument through a `tox_` buffer, rather than in place."""
        if not _r_marshals(argument):
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
        if self._is_inout_copy(argument):
            return shim_local(argument.name, INOUT_COPY)
        return argument.name

    def _copy_inout(self, writer: Writer, wrapper: CWrapper) -> None:
        """A plain array modified in place is duplicated, so the caller's stays intact."""
        lines = Writer()
        for argument in wrapper:
            if self._is_inout_copy(argument):
                copy = shim_local(argument.name, INOUT_COPY)
                # an integer array arrives as R holds it, and is checked and made one here
                duplicate = (f'tox_int_copy({argument.name}, "{argument.name}")'
                             if argument.type.base is BaseType.INTEGER
                             else f"Rf_duplicate({argument.name})")
                lines.line(f"SEXP {copy} = PROTECT({duplicate}); {C_PROTECT_COUNT}++;")
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
        """Convert one input for Fortran, checking each element as it goes (`tox_marshal.h`).

        The argument's name travels along, so a refusal can say which one it was.
        """
        name = argument.name
        buffer = shim_local(name, BUFFER)
        if argument.type.is_logical:
            return f'unsigned char* {buffer} = tox_bool_in({name}, "{name}");'
        if _checks_integer_array(argument):
            return f'const int* {buffer} = tox_int_in({name}, "{name}");'
        length = self._extent(argument.dimension.extents[0])
        return f'char* {buffer} = tox_char_in({name}, {length}, "{name}");'

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

        if argument.type.is_logical:
            return [f"unsigned char* {shim_local(name, BUFFER)} = tox_bool_alloc({size});"]
        if argument.type.is_character:
            length = extents[0]
            count = size if argument.shape_arg is not None else _product(extents[1:])
            return [f"char* {shim_local(name, BUFFER)} = tox_char_alloc({length}, {count});"]
        if argument.is_temporary:
            ctype = c_ctype(argument)
            return [f"{ctype}* {name} = ({ctype}*) R_alloc({size}, sizeof({ctype}));"]
        # a plain output: an R vector, protected, with a dim attribute if rank >= 2
        sexptype = R_SEXPTYPE[argument.type.base]
        out = [f"SEXP {name} = PROTECT(Rf_allocVector({sexptype}, {size})); {C_PROTECT_COUNT}++;"]
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
        dim = shim_local(argument.name, DIM)
        assigns = " ".join(f"INTEGER({dim})[{i}] = {e};" for i, e in enumerate(extents))
        return [
            f"{{ SEXP {dim} = PROTECT(Rf_allocVector(INTSXP, {rank})); "
            f"{assigns} Rf_setAttrib({target}, R_DimSymbol, {dim}); UNPROTECT(1); }}"
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
                return f"{name} != R_NilValue ? {shim_local(name, BUFFER)} : NULL"
            return shim_local(name, OPTIONAL_POINTER)
        if self._converts_via_buffer(argument):
            return shim_local(name, BUFFER)
        if argument.is_scalar:
            # a by-value scalar input lives in _v_<name>; a derived extent, scalar output or
            # ierr is a plain <name> local
            ref = shim_local(name, SCALAR_VALUE) if name in self._scalar_names else name
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
            buffer = shim_local(name, BUFFER)
            count = self._converted_count(argument)
            if argument.type.is_logical:
                lines.line(
                    f"SEXP {name} = PROTECT(tox_bool_out({buffer}, {count})); {C_PROTECT_COUNT}++;"
                )
            else:
                length = self._extent(argument.dimension.extents[0])
                lines.line(
                    f"SEXP {name} = PROTECT(tox_char_out({buffer}, {length}, {count})); "
                    f"{C_PROTECT_COUNT}++;"
                )
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

        writer.line(f"SEXP _out = PROTECT(Rf_allocVector(VECSXP, {n})); {C_PROTECT_COUNT}++;")
        for index, (_, value) in enumerate(entries):
            writer.line(f"SET_VECTOR_ELT(_out, {index}, {value});")
        writer.line(f"SEXP _nms = PROTECT(Rf_allocVector(STRSXP, {n})); {C_PROTECT_COUNT}++;")
        for index, (nm, _) in enumerate(entries):
            writer.line(f'SET_STRING_ELT(_nms, {index}, Rf_mkChar("{nm}"));')
        writer.line("Rf_setAttrib(_out, R_NamesSymbol, _nms);")
        writer.line(f"UNPROTECT({C_PROTECT_COUNT});")
        writer.line("return _out;")

    def _result_sexp(self, argument: CArgument) -> str:
        """The SEXP for one returned output, wrapping a scalar or naming an allocated one."""
        if argument.is_scalar:
            base = argument.type.base
            # an inout scalar was pulled into _v_<name> and written there; a pure output
            # scalar is a plain <name> local
            name = argument.name
            if name in self._scalar_names:
                name = shim_local(name, SCALAR_VALUE)
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
// Fortran cannot take from R directly -- R's int-based logicals, its strings, and an integer
// array that arrives as doubles -- plus a couple of small shape helpers. Transient buffers
// come from R_alloc and are freed when the .Call returns, so there is nothing to free by hand
// and nothing to leak on an error longjmp.
//
// The conversions are where the input contract's element checks live (#221): the R wrapper
// has checked each argument's kind and rank, and every element of an integer, logical or
// character array is checked here, in the pass that converts it, which costs no extra pass.
// A refusal is raised as the R wrapper raises its own: a classed condition naming the
// argument (tox_refuse).
//
// Strings cross fixed-width, column-major and blank-padded, both ways: the Fortran wrapper
// takes a character(len=strlen) pointer view of the buffer rather than converting it, so
// the padding is part of the value and it has to be the padding Fortran trims.
#ifndef TOX_MARSHAL_H
#define TOX_MARSHAL_H

#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

// __WEAK_PRAGMAS__

static inline int tox_imax(int a, int b) { return a > b ? a : b; }
static inline int tox_imin(int a, int b) { return a < b ? a : b; }

// Refuse an argument: signal the condition the R wrapper's own refusals are
// (.tox_argument_error in error_handling.R) -- class c(<condition>, __ARGUMENT_CLASSES__),
// code NA, the argument named -- through base R's stop(), which needs nothing of the
// environment the wrappers were sourced into. It does not return. On the refusal path only.
static inline void tox_refuse(const char* argument, const char* condition, const char* format, ...) {
    char message[512];
    va_list values;
    va_start(values, format);
    vsnprintf(message, sizeof message, format, values);
    va_end(values);
    const char* fields[] = {__CONDITION_FIELDS__};
    const char* classes[] = {condition, __ARGUMENT_CLASSES__};
    int n_fields = (int) (sizeof fields / sizeof fields[0]);
    int n_classes = (int) (sizeof classes / sizeof classes[0]);
    SEXP refusal = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SET_VECTOR_ELT(refusal, 0, Rf_mkString(message));
    SET_VECTOR_ELT(refusal, 1, R_NilValue);
    SET_VECTOR_ELT(refusal, 2, Rf_ScalarInteger(NA_INTEGER));
    SET_VECTOR_ELT(refusal, 3, Rf_ScalarString(NA_STRING));
    SET_VECTOR_ELT(refusal, 4, Rf_mkString(argument));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    for (int i = 0; i < n_fields; ++i) SET_STRING_ELT(names, i, Rf_mkChar(fields[i]));
    Rf_setAttrib(refusal, R_NamesSymbol, names);
    SEXP class_ = PROTECT(Rf_allocVector(STRSXP, n_classes));
    for (int i = 0; i < n_classes; ++i) SET_STRING_ELT(class_, i, Rf_mkChar(classes[i]));
    Rf_setAttrib(refusal, R_ClassSymbol, class_);
    SEXP call = PROTECT(Rf_lang2(Rf_install("stop"), refusal));
    Rf_eval(call, R_BaseEnv);
    UNPROTECT(4);
}

static inline void tox_refuse_na(const char* argument) {
    tox_refuse(argument, "tox_na_error", "__NA_MESSAGE__", argument);
}

// A value as a refusal shows it: R's own spelling of NA, NaN and the infinities
static inline void tox_show(char* text, size_t size, double value) {
    if (ISNA(value)) snprintf(text, size, "NA");
    else if (ISNAN(value)) snprintf(text, size, "NaN");
    else if (value == R_PosInf) snprintf(text, size, "Inf");
    else if (value == R_NegInf) snprintf(text, size, "-Inf");
    else snprintf(text, size, "%.15g", value);
}

// An element of x that is not what `expected` says, 0-based i shown 1-based as R counts
static inline void tox_refuse_element(const char* argument, const char* expected, R_xlen_t i,
                                      double value) {
    char shown[32];
    tox_show(shown, sizeof shown, value);
    tox_refuse(argument, "tox_type_error", "'%s' must hold %s, but element %lld is %s",
               argument, expected, (long long) i + 1, shown);
}

// What the R wrapper guarantees, if something called the .Call entry directly instead
static inline void tox_refuse_type(const char* argument, const char* expected, SEXP x) {
    tox_refuse(argument, "tox_type_error", "'%s' must be %s, not an R %s",
               argument, expected, Rf_type2char(TYPEOF(x)));
}

// An integer array the caller passed, as the int32 Fortran is handed: R's integer vector as
// it stands, NA refused, or a double vector whose every element is a whole number in R's
// integer range -- converted only after the range is checked, so the cast is always defined.
// NULL for an absent (R_NilValue) optional. INTEGER_NO_NA answers at once for a sequence like
// 1:n, which R knows holds no NA.
static inline const int* tox_int_in(SEXP x, const char* argument) {
    if (x == R_NilValue) return NULL;
    R_xlen_t n = XLENGTH(x);
    if (TYPEOF(x) == INTSXP) {
        const int* values = INTEGER(x);
        if (!INTEGER_NO_NA(x))
            for (R_xlen_t i = 0; i < n; ++i)
                if (values[i] == NA_INTEGER) tox_refuse_na(argument);
        return values;
    }
    if (TYPEOF(x) != REALSXP) tox_refuse_type(argument, "an integer vector", x);
    const double* values = REAL(x);
    int* converted = (int*) R_alloc((size_t) (n > 0 ? n : 1), sizeof(int));
    for (R_xlen_t i = 0; i < n; ++i) {
        double value = values[i];
        if (!(value >= -INT_MAX && value <= INT_MAX) || value != (double) (int) value) {
            if (ISNA(value)) tox_refuse_na(argument);
            tox_refuse_element(argument, "whole numbers in the integer range", i, value);
        }
        converted[i] = (int) value;
    }
    return converted;
}

// An integer array the callee modifies, as a fresh integer vector with x's attributes:
// checked as tox_int_in checks it, then copied or converted
static inline SEXP tox_int_copy(SEXP x, const char* argument) {
    tox_int_in(x, argument);
    return TYPEOF(x) == INTSXP ? Rf_duplicate(x) : Rf_coerceVector(x, INTSXP);
}

// Whether s holds n bytes of well-formed UTF-8: no overlong form, no surrogate, nothing
// past U+10FFFF -- what R's validUTF8() accepts
static inline int tox_utf8_valid(const unsigned char* s, size_t n) {
    size_t i = 0;
    while (i < n) {
        unsigned char c = s[i];
        size_t length;
        unsigned char low = 0x80, high = 0xBF;
        if (c < 0x80) { ++i; continue; }
        else if (c >= 0xC2 && c <= 0xDF) length = 2;
        else if (c >= 0xE0 && c <= 0xEF) {
            length = 3;
            if (c == 0xE0) low = 0xA0;
            if (c == 0xED) high = 0x9F;
        } else if (c >= 0xF0 && c <= 0xF4) {
            length = 4;
            if (c == 0xF0) low = 0x90;
            if (c == 0xF4) high = 0x8F;
        } else return 0;
        if (i + length > n) return 0;
        if (s[i + 1] < low || s[i + 1] > high) return 0;
        for (size_t k = 2; k < length; ++k)
            if (s[i + k] < 0x80 || s[i + k] > 0xBF) return 0;
        i += length;
    }
    return 1;
}

// A string as the UTF-8 bytes Fortran is handed, with their count in *length -- or NULL, and
// *problem saying why it cannot cross. By what R marks it as, the same in every locale:
// - ASCII, or marked UTF-8: as it stands, if it is well-formed;
// - marked latin1: translated, as enc2utf8() translates it (R reads latin1 as Windows-1252),
//   unless it holds one of the bytes that code page leaves undefined;
// - unmarked (native): as it stands if its bytes are well-formed UTF-8, and refused otherwise.
//   Never translated from the locale: in the C locale R takes every non-ASCII byte for one it
//   cannot convert and writes it out as escape text, so the UTF-8 bytes c3 a9 of an e-acute
//   reached Fortran as the eight characters "<c3><a9>";
// - marked "bytes": it has no encoding to read it in.
// So no escape text R makes up can reach Fortran in place of the caller's string.
static inline const char* tox_utf8(SEXP e, size_t* length, const char** problem) {
    const char* text = CHAR(e);
    *length = (size_t) LENGTH(e);
    *problem = NULL;
    if (Rf_charIsASCII(e)) return text;  // and without a NUL: R cannot hold one
    cetype_t encoding = Rf_getCharCE(e);
    if (encoding == CE_BYTES) {
        *problem = "is marked as bytes, which has no encoding to read it in";
        return NULL;
    }
    if (encoding == CE_LATIN1) {
        const char* translated = Rf_translateCharUTF8(e);
        // every latin1 byte is one character, so a translation holding more characters than
        // the string has bytes is one where R wrote an undefined byte out as "<81>"
        size_t characters = 0;
        for (const unsigned char* p = (const unsigned char*) translated; *p; ++p)
            characters += (*p & 0xC0) != 0x80;
        if (characters != *length) {
            *problem = "is marked latin1, but holds a byte that is no character in it";
            return NULL;
        }
        *length = strlen(translated);
        return translated;
    }
    if (!tox_utf8_valid((const unsigned char*) text, *length)) {
        *problem = "is not valid UTF-8 (a string in another encoding needs marking as one, "
                   "such as Encoding(x) <- \"latin1\")";
        return NULL;
    }
    return text;
}

// The width a character(len=n) array needs: the longest element, NA skipped. 0 for a
// non-character or absent (R_NilValue) argument, so an omitted optional reports no width.
//
// A *present* argument reports at least 1, matching what the Python layer does with
// `.ljust(1)` and `max(..., default=0) or 1`. Every element being "" would otherwise give
// width 0, and `R_alloc(0, 1)` returns NULL -- a null buffer the Fortran wrapper cannot make
// a pointer view of. The floor is deliberately inside the present branch: raising the absent
// case to 1 would send width 1 beside a null pointer, and M_CHECK_ARRAY_NON_NULL would then
// reject an optional the caller legitimately omitted.
//
// Measured in UTF-8 bytes, as tox_char_in copies it. An NA, and a string that cannot cross
// as UTF-8, is skipped here and refused there.
static inline int tox_max_strlen(SEXP x) {
    if (x == R_NilValue || TYPEOF(x) != STRSXP) return 0;
    int longest = 0, n = (int) XLENGTH(x);
    const SEXP* elements = STRING_PTR_RO(x);
    for (int i = 0; i < n; ++i) {
        SEXP e = elements[i];
        if (e == NA_STRING) continue;
        size_t length;
        const char* problem;
        if (tox_utf8(e, &length, &problem) == NULL) continue;
        if ((int) length > longest) longest = (int) length;
    }
    return longest > 0 ? longest : 1;
}

// Fortran carries a string's length as the leading extent: n strings of length len are
// char[len * n], column-major, each **blank**-padded. The Fortran wrapper does not convert
// this buffer, it takes a `character(len=len)` pointer view of it -- so the padding byte
// becomes part of the value, and blanks are the padding Fortran's own `trim`, `len_trim`
// and string comparison understand. A NUL would survive into the string instead.
//
// Every element is checked as it is copied: an NA, a string that is not well-formed UTF-8,
// and one longer than len are refused, naming the argument. A longer one could only be cut
// short, which turned `delimiter = "::"` into ":" (len is the longest element's length
// wherever the caller's strings decide it, so this refuses only past a declared length).
static inline char* tox_char_in(SEXP x, int len, const char* argument) {
    int n = (x == R_NilValue) ? 0 : (int) XLENGTH(x);
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char* buf = (char*) R_alloc(total, 1);
    memset(buf, ' ', total);
    if (n == 0) return buf;
    if (TYPEOF(x) != STRSXP) tox_refuse_type(argument, "a character vector", x);
    const SEXP* elements = STRING_PTR_RO(x);
    for (int i = 0; i < n; ++i) {
        SEXP e = elements[i];
        if (e == NA_STRING) tox_refuse_na(argument);
        size_t slen;
        const char* problem;
        const char* text = tox_utf8(e, &slen, &problem);
        if (text == NULL)
            tox_refuse(argument, "tox_type_error", "'%s' element %d %s", argument, i + 1, problem);
        if (slen > (size_t) len)
            tox_refuse(argument, "tox_type_error",
                       "'%s' element %d is %d bytes long in UTF-8, but it is declared %d long, "
                       "and a longer value would be cut short",
                       argument, i + 1, (int) slen, len);
        memcpy(buf + (size_t) i * len, text, slen);
    }
    return buf;
}

// A blank-filled c_char output buffer of len * n bytes. Blank and not zero, deliberately:
// the callee may leave an element untouched (an early error return, a partly filled array),
// and tox_char_out reads it back by trimming trailing blanks. A zero fill would leave NULs
// there, which Rf_mkCharLen rejects outright with an embedded-nul error.
static inline char* tox_char_alloc(int len, int n) {
    size_t total = (size_t) len * (n > 0 ? n : 1);
    char* buf = (char*) R_alloc(total, 1);
    memset(buf, ' ', total);
    return buf;
}

// char[len * n] fixed-width -> STRSXP, each slot with its trailing blanks removed: Fortran
// blank-pads whatever it assigns into a character(len=n), and the wrapper hands that buffer
// straight through. Trailing NULs are deliberately not stripped -- nothing writes them any
// more, and Rf_mkCharLen turning a stray one into a loud R error is the right answer.
//
// `out` is protected for the whole fill, unlike every other helper here: Rf_mkCharLen
// allocates, so it can trigger a GC on any iteration, and until the caller receives `out`
// nothing else references it. Protecting on return -- which is what the caller's
// PROTECT(tox_char_out(...)) does -- is too late; the window is the loop, and it grows with
// n. Balanced before returning, so the caller's own PROTECT is still the one that keeps it.
static inline SEXP tox_char_out(const char* buf, int len, int n) {
    SEXP out = PROTECT(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; ++i) {
        const char* p = buf + (size_t) i * len;
        int m = len;
        while (m > 0 && p[m - 1] == ' ') --m;
        // UTF-8 is what crossed, both ways, so that is what the string is marked as: in a
        // non-UTF-8 locale an unmarked one would be read in the native encoding instead
        SET_STRING_ELT(out, i, Rf_mkCharLenCE(p, m, CE_UTF8));
    }
    UNPROTECT(1);
    return out;
}

// R logical is a 4-byte int (0/1); Fortran c_bool is one byte. Copied across, and so is a
// number whose every element is exactly 0 or 1 -- c(1, 0, 1) is R's own literal for a 0/1
// mask -- while an NA, and any other number, is refused by name.
static inline unsigned char* tox_bool_in(SEXP x, const char* argument) {
    R_xlen_t n = (x == R_NilValue) ? 0 : XLENGTH(x);
    unsigned char* buf = (unsigned char*) R_alloc((size_t) (n > 0 ? n : 1), 1);
    if (n == 0) return buf;
    static const char expected[] = "TRUE and FALSE, or the numbers 0 and 1";
    switch (TYPEOF(x)) {
    case LGLSXP: {
        const int* values = LOGICAL(x);
        for (R_xlen_t i = 0; i < n; ++i) {
            if (values[i] == NA_LOGICAL) tox_refuse_na(argument);
            buf[i] = values[i] != 0;
        }
        break;
    }
    case INTSXP: {
        const int* values = INTEGER(x);
        for (R_xlen_t i = 0; i < n; ++i) {
            int value = values[i];
            if (value != 0 && value != 1) {
                if (value == NA_INTEGER) tox_refuse_na(argument);
                tox_refuse_element(argument, expected, i, (double) value);
            }
            buf[i] = (unsigned char) value;
        }
        break;
    }
    case REALSXP: {
        const double* values = REAL(x);
        for (R_xlen_t i = 0; i < n; ++i) {
            double value = values[i];
            if (value != 0.0 && value != 1.0) {
                if (ISNA(value)) tox_refuse_na(argument);
                tox_refuse_element(argument, expected, i, value);
            }
            buf[i] = value != 0.0;
        }
        break;
    }
    default:
        tox_refuse_type(argument, "a logical vector", x);
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
// which must map to R's TRUE (1), not stay 255. Returned unprotected, and safely so: the one
// allocation is the first statement and nothing after it allocates, so no GC can run while
// the result is unreferenced.
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

// Count the elements of a mask that are set: the TRUEs, or the 1s of a 0/1 number vector.
// Counted before the mask is converted, so an element that is neither -- which tox_bool_in
// refuses next -- counts as unset, as R's count of it does (r_wrapper._mask_count).
static inline int tox_sum_true(SEXP x) {
    int n = (int) XLENGTH(x), total = 0;
    if (TYPEOF(x) == REALSXP) {
        const double* values = REAL(x);
        for (int i = 0; i < n; ++i) total += values[i] == 1.0;
    } else if (TYPEOF(x) == LGLSXP || TYPEOF(x) == INTSXP) {
        const int* values = TYPEOF(x) == LGLSXP ? LOGICAL(x) : INTEGER(x);
        for (int i = 0; i < n; ++i) total += values[i] == 1;
    }
    return total;
}

#endif
'''
