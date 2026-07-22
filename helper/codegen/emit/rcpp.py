"""Emitting the C++ half of the R interface.

The split with R is deliberate and documented in `design/language-layers.md`: **R decides
and raises, C++ marshals and calls.** By the time a value reaches C++ it has been checked
and coerced in R, so nothing here validates -- it converts what C cannot take directly
(R's int-based logicals, its strings) and hands everything to the `extern "C"` wrapper.

Each generated function is `.<name>_rcpp`, dot-prefixed because it is internal: the
user-facing `<name>` is the R wrapper on top of it. It returns a list of every output
plus `ierr`, which the R wrapper decodes.
"""

from __future__ import annotations

import re

from ..abi.model import CArgument, Conversion, CWrapper, CWrapperModule, Origin
from ..ir.types import BaseType, Intent
from ..render import Writer

#: A Fortran numeric literal's kind suffix -- `0_int32`. C++ has no such thing.
_FORTRAN_KIND_SUFFIX = re.compile(r"\b(\d+(?:\.\d+)?)_[A-Za-z]\w*")
#: `size(arr)` / `size(arr, dim)` -- an extent stated in terms of another argument.
_FORTRAN_SIZE = re.compile(r"\bsize\s*\(\s*([A-Za-z_]\w*)\s*(?:,\s*(\d+)\s*)?\)")
#: `max`/`min` resolve to the C++ standard versions, which need qualifying.
_FORTRAN_MINMAX = re.compile(r"\b(max|min)\s*\(")


def _size_as_cpp(match: re.Match) -> str:
    array, dim = match.group(1), match.group(2)
    if dim is None:
        return f"(int) {array}.size()"
    # Fortran dimensions are 1-based, and a matrix carries its extents in `dim`
    return f"(int) IntegerVector({array}.attr(\"dim\"))[{int(dim) - 1}]"


def cpp_extent(extent: str) -> str:
    """A Fortran extent expression as C++.

    The same three incompatibilities the Python emitter handles -- kind suffixes, which
    C++ has no notion of; `size(...)`, which is a member function here; and `max`/`min`,
    which must name the std versions rather than rely on lookup finding something.
    """
    extent = _FORTRAN_KIND_SUFFIX.sub(r"\1", extent)
    extent = _FORTRAN_SIZE.sub(_size_as_cpp, extent)
    return _FORTRAN_MINMAX.sub(r"std::\1(", extent)

#: iso_c_binding kind -> the C type it interoperates as, for the extern "C" declaration
CPP_CTYPE = {
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
    "c_bool": "bool",
    "c_char": "char",
}

#: Fortran base type -> the Rcpp vector class for an array of it
RCPP_VECTOR = {
    BaseType.INTEGER: "IntegerVector",
    BaseType.REAL: "NumericVector",
    BaseType.COMPLEX: "ComplexVector",
    BaseType.LOGICAL: "LogicalVector",
    BaseType.CHARACTER: "CharacterVector",
}

#: Fortran base type -> the C++ scalar for a by-value input
RCPP_SCALAR = {
    BaseType.INTEGER: "int",
    BaseType.REAL: "double",
    BaseType.COMPLEX: "std::complex<double>",
    BaseType.LOGICAL: "bool",
}


def _product(extents) -> str:
    """The product of some extents, each parenthesised.

    An extent may be an expression (`n_timepoints - 1`), and `*` binds tighter than `-`
    in both C++ and Fortran: joining raw gives `n_timepoints - 1 * n_factors`, which is
    a different number and, for a buffer size, a heap overrun waiting to happen.
    """
    parts = [e if e.isidentifier() or e.isdigit() else f"({e})" for e in extents]
    return " * ".join(parts) if parts else "1"


def cpp_ctype(argument: CArgument) -> str:
    return CPP_CTYPE[argument.type.kind]


class RcppEmitter:
    def __init__(self, marshal_header: str = "tox_marshal.h"):
        self.marshal_header = marshal_header

    def marshal_header_content(self) -> str:
        """The buffer helpers, in one header rather than repeated per module.

        These are the only conversions C++ does: R's logicals are `int`-based and its
        strings are `SEXP`s, neither of which Fortran can take. R has already rejected
        `NA`, so the conversions are straight copies.
        """
        return _MARSHAL_HEADER

    # -- module -----------------------------------------------------------------

    def module(self, module: CWrapperModule) -> str:
        writer = Writer()
        writer.line("// Generated. Do not edit.")
        writer.line("#include <Rcpp.h>")
        writer.line("#include <numeric>")
        writer.line("#include <functional>")
        writer.line(f'#include "{self.marshal_header}"')
        writer.blank()
        writer.line("using namespace Rcpp;")
        writer.blank()

        writer.line('extern "C" {')
        with writer.indent():
            for wrapper in module:
                writer.line(self._extern_declaration(wrapper))
        writer.line("}")
        writer.blank()

        for wrapper in module:
            writer.block(self.function(wrapper))
            writer.blank(collapse=False)

        return writer.render(trailing_newline=True)

    def _extern_declaration(self, wrapper: CWrapper) -> str:
        params = ", ".join(self._extern_param(a) for a in wrapper)
        return f"void {wrapper.name}({params});"

    def _extern_param(self, argument: CArgument) -> str:
        const = "const " if argument.intent is Intent.IN else ""
        return f"{const}{cpp_ctype(argument)}*"

    # -- function ---------------------------------------------------------------

    def function(self, wrapper: CWrapper) -> str:
        writer = Writer()
        writer.line(f"// [[Rcpp::export(.{wrapper.stripped_name}_rcpp)]]")
        params = ", ".join(self._param(a) for a in self._inputs(wrapper))
        writer.line(f"List {wrapper.stripped_name}_rcpp({params}) {{")

        with writer.indent():
            self._materialize_optionals(writer, wrapper)
            self._derive(writer, wrapper)
            self._copy_inout(writer, wrapper)
            self._marshal_inputs(writer, wrapper)
            self._allocate(writer, wrapper)
            self._call(writer, wrapper)
            self._marshal_outputs(writer, wrapper)
            self._return(writer, wrapper)

        writer.line("}")
        return writer.render()

    def _inputs(self, wrapper: CWrapper) -> list[CArgument]:
        """The arguments the R wrapper passes in -- the same set Python asks for."""
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

        One that sizes an input is read off that input; one that sizes an *output* has no
        such source, so leaving it out would drop it from the signature with nothing to
        compute it from. Mirrors the Python emitter.
        """
        if argument.is_error:
            return False
        return self._derived_value(argument, wrapper) is None

    @staticmethod
    def _is_derived(argument: CArgument) -> bool:
        """Whether C++ works this out itself, so it is not a parameter of .rcpp.

        A computed (AUTO) argument is *not* derived here: R computes it and passes it in,
        because the producer's wrapper is an R function C++ cannot call. So from C++'s
        side a computed argument is an ordinary input. Extents, shapes and mask counts C++
        does derive, from the input arrays it is handed.
        """
        roles = argument.source.roles if argument.source else None
        if roles is None:
            return False
        if roles.is_computed:
            # computed wins: R passes it in, even though it may also size a work array
            return False
        return (roles.is_inferable_extent or roles.is_inferable_shape_arg
                or roles.is_mask_count)

    def _param(self, argument: CArgument) -> str:
        if argument.optional:
            # a nullable optional arrives as R's NULL when absent, which Rcpp's Nullable
            # carries; the C wrapper then receives a null pointer and Fortran sees it
            # absent (the OPTIONAL bind(C) rule)
            return f"Nullable<{self._rcpp_class(argument)}> {argument.name} = R_NilValue"
        if argument.conversion is Conversion.MODE:
            return f"CharacterVector {argument.name}"
        if argument.is_scalar:
            return f"{RCPP_SCALAR[argument.type.base]} {argument.name}"
        return f"{RCPP_VECTOR[argument.type.base]} {argument.name}"

    def _rcpp_class(self, argument: CArgument) -> str:
        if argument.conversion is Conversion.MODE or argument.type.is_character:
            return "CharacterVector"
        return RCPP_VECTOR[argument.type.base]

    # -- body -------------------------------------------------------------------

    def _derive(self, writer: Writer, wrapper: CWrapper) -> None:
        """Compute the extents, shapes and counts the call needs, from the inputs.

        A shape argument comes first: C++ needs it declared before the extent that is the
        product of it, and the wrapper's argument order is the Fortran one, which says
        nothing about that.
        """
        lines = Writer()

        def is_shape(argument):
            roles = argument.source.roles if argument.source else None
            return bool(roles and roles.is_shape_arg)

        ordered = sorted(wrapper, key=lambda argument: not is_shape(argument))
        for argument in ordered:
            expression = self._derived_value(argument, wrapper)
            if expression is not None:
                lines.line(f"{self._cpp_local_type(argument)} {argument.name} = {expression};")
        if lines:
            writer.line("// derived from the inputs, not asked of the caller")
            writer.extend(lines)
            writer.blank()

    def _derived_value(self, argument: CArgument, wrapper: CWrapper) -> str | None:
        source = argument.source
        roles = source.roles if source else None

        if argument.origin is Origin.STRLEN:
            owner = wrapper.argument(argument.sizes)
            # only off an input: an output buffer is not built yet, so its item length
            # cannot be read back and the caller has to state it
            if owner is None or not owner.intent.is_input:
                return None
            # a Nullable has no .length(): read the value materialized above, which is an
            # empty vector when the caller omitted the argument
            buffer = f"{argument.sizes}_val" if owner.optional else argument.sizes
            # the longest string, not the first: Fortran's character(len=n) is one width
            # for the whole array, and sizing it from element 0 silently truncates every
            # longer one
            return f"tox::max_strlen({buffer})"

        if argument.origin is Origin.EXTENT:
            owner = wrapper.argument(argument.sizes)
            # not off a work array either: that is sized *from* this extent
            if owner is None or not owner.intent.is_input or owner.is_temporary:
                return None
            return self._extent_of(owner, argument.axis)

        if roles is None:
            return None

        if roles.is_mask_count:
            # R has already rejected NA, so a plain sum counts the TRUEs
            return f"(int) sum({roles.mask_count_of.name})"

        if roles.is_shape_arg:
            if not roles.shape_of.intent.is_input:
                # it says what shape the output should have, so the caller states it
                return None
            # a plain R vector carries no dim attribute, and its shape is its length
            owner = roles.shape_of.name
            return (f"Rf_isNull({owner}.attr(\"dim\")) "
                    f"? IntegerVector::create({owner}.size()) "
                    f": IntegerVector({owner}.attr(\"dim\"))")

        if roles.is_extent:
            for owner in roles.extent_of:
                c_owner = wrapper.argument(owner.name)
                if c_owner is None or c_owner.shape_arg is None:
                    continue
                # its extents travel separately, so an extent of it is the total count:
                # read off the array itself when the caller supplies it, otherwise the
                # product of the shape they stated
                if c_owner.intent.is_input:
                    return f"(int) {owner.name}.size()"
                return (f"(int) std::accumulate({c_owner.shape_arg}.begin(), "
                        f"{c_owner.shape_arg}.end(), 1, std::multiplies<int>())")
            # a required owner first, since its size is always readable; fall back to an
            # optional one, whose materialized size is 0 when the caller omits it
            best = None
            for owner in roles.extent_of:
                c_owner = wrapper.argument(owner.name)
                if c_owner is None or not c_owner.intent.is_input:
                    continue
                if c_owner.is_temporary:
                    # sized from this extent, so it cannot also be its source
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
        """The size of `owner` along `axis`, as an Rcpp expression."""
        if owner.optional:
            # materialized to 0 when absent, so nothing reads a buffer that is not there
            if owner.rank > 1 and not owner.type.is_character:
                return f"{owner.name}_dim[{axis}]"
            return f"{owner.name}_size"
        if owner.rank <= 1 or (owner.type.is_character and owner.rank == 2):
            return f"(int) {owner.name}.size()"
        # a matrix or higher: read the dim attribute
        return f"(int) IntegerVector({owner.name}.attr(\"dim\"))[{axis}]"

    @staticmethod
    def _axis_of(extent: str, owner: CArgument) -> int | None:
        extents = list(owner.dimension.extents)
        if owner.type.is_character and extents:
            extents = extents[1:]
        try:
            return extents.index(extent)
        except ValueError:
            return None

    def _cpp_local_type(self, argument: CArgument) -> str:
        if argument.source and argument.source.roles and argument.source.roles.is_shape_arg:
            return "IntegerVector"
        return "int"

    def _materialize_optionals(self, writer: Writer, wrapper: CWrapper) -> None:
        """Turn a nullable optional into a pointer and a size, guarded on presence.

        An absent optional is a null pointer to the C wrapper (Fortran then sees it
        absent), and a size of 0, which is what the extent derived from it reports so no
        one reads a buffer that is not there.
        """
        lines = Writer()
        for argument in wrapper:
            if not argument.optional:
                continue
            cls = self._rcpp_class(argument)
            ctype = cpp_ctype(argument)
            name = argument.name
            rank = argument.rank
            # a character's leading extent is its item length, not an axis, so its
            # strings are a vector however many C extents it has
            has_axes = rank > 1 and not argument.type.is_character
            # A converted argument (character, mode, logical) reaches C through its own
            # buffer, built from `_val` below, so no raw pointer into R's memory is taken.
            if not argument.needs_conversion:
                lines.line(f"const {ctype}* {name}_p = nullptr;")
            lines.line(f"int {name}_size = 0;")
            if has_axes:
                # a matrix needs its extents per axis, not just its total length
                lines.line(f"std::vector<int> {name}_dim({rank}, 0);")
            lines.line(f"{cls} {name}_val;")
            lines.line(f"if ({name}.isNotNull()) {{")
            with lines.indent():
                lines.line(f"{name}_val = {name}.get();")
                if not argument.needs_conversion:
                    lines.line(
                        f"{name}_p = {self._as_c_pointer(argument, f'{name}_val.begin()')};"
                    )
                lines.line(f"{name}_size = {name}_val.size();")
                if has_axes:
                    # R only carries a dim attribute on a matrix; the R wrapper coerces
                    # to one, but a direct .rcpp call need not have
                    lines.line(f"if (!Rf_isNull({name}_val.attr(\"dim\"))) {{")
                    with lines.indent():
                        lines.line(f"IntegerVector {name}_d({name}_val.attr(\"dim\"));")
                        lines.line(
                            f"for (int i = 0; i < {rank} && i < {name}_d.size(); ++i)"
                            f" {name}_dim[i] = {name}_d[i];"
                        )
                    lines.line("}")
            lines.line("}")
        if lines:
            writer.line("// optionals: a null pointer and size 0 when the caller omits them")
            writer.extend(lines)
            writer.blank()

    @staticmethod
    def _converts_via_buffer(argument: CArgument) -> bool:
        """Whether this argument reaches C through a tox:: buffer rather than directly.

        Conversion is an array concern. A logical *scalar* is already C++'s `bool`, which
        is what the C wrapper takes, so it goes by address like any other scalar; only a
        logical array needs unpacking from R's int-based representation. A character is
        always a buffer, scalar or not, because R hands over a SEXP.
        """
        if not argument.needs_conversion:
            return False
        return argument.is_array or argument.type.is_character

    @staticmethod
    def _is_inout_copy(argument: CArgument) -> bool:
        """A plain array modified in place, which must be cloned before Fortran touches it.

        Rcpp shares R's buffer, so writing through it would modify the caller's object --
        which R's copy-on-modify contract forbids. `clone` copies reliably, where R's own
        coercion copies only sometimes (`as.double` on an already-double vector does not).
        Conversion cases build their own buffer, so they need no clone. Neither do work
        arrays: the wrapper allocates those itself, so there is no caller object to
        protect and nothing to clone from.
        """
        return (
            argument.intent is Intent.INOUT
            and argument.is_array
            and not argument.needs_conversion
            and not argument.is_temporary
        )

    def _working_name(self, argument: CArgument) -> str:
        return f"{argument.name}_out" if self._is_inout_copy(argument) else argument.name

    def _copy_inout(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            if self._is_inout_copy(argument):
                rcpp = RCPP_VECTOR[argument.type.base]
                lines.line(f"{rcpp} {argument.name}_out = clone({argument.name});")
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
            writer.line("// convert what C cannot take directly")
            writer.extend(lines)
            writer.blank()

    def _input_buffer(self, argument: CArgument) -> str:
        name = argument.name
        # an optional arrives as a Nullable; the buffer is built from the value
        # materialized above, which is empty when the caller omitted the argument
        value = f"{name}_val" if argument.optional else name
        if argument.conversion is Conversion.LOGICAL:
            return f"tox::BoolBuffer {name}_c({value});"
        # character or mode: a padded c_char buffer
        length = argument.dimension.extents[0]
        return f"tox::CharBuffer {name}_c({value}, {length});"

    def _allocate(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            if argument.intent.is_input and not argument.is_temporary:
                continue
            if self._is_derived(argument) and argument.intent.is_input:
                continue
            lines.line(self._new_value(argument))
            if not argument.is_temporary and not self._converts_via_buffer(argument):
                dim = self._dim_attribute(argument)
                if dim:
                    lines.line(dim)
        if lines:
            writer.line("// outputs and work space")
            writer.extend(lines)
            writer.blank()

    def _r_extents(self, argument: CArgument) -> list[str]:
        """The extents R sees, as C++ expressions.

        A character's leading extent is its item length, which is not an R dimension --
        its strings are an array of one dimension fewer.
        """
        extents = [cpp_extent(e) for e in argument.dimension.extents]
        if argument.type.is_character and extents:
            extents = extents[1:]
        return extents

    def _dim_attribute(self, argument: CArgument) -> str | None:
        """Give a rank-2-or-more output the `dim` that makes it a matrix in R.

        Without it a freshly allocated output comes back as a flat vector, and every
        `dim(x)` or `x[, 1]` on the R side fails. (An `intent(inout)` array escapes this
        by accident: it is cloned from the caller's object, attributes and all.)
        """
        if argument.shape_arg is not None:
            # its extents travelled separately; apply them so R gets the n-d array back
            # rather than a flat vector and a shape to reapply itself. R is column-major,
            # as Fortran wrote the block, so the dim is all that is needed.
            return f'{self._working_name(argument)}.attr("dim") = {argument.shape_arg};'
        extents = self._r_extents(argument)
        if len(extents) < 2:
            return None
        joined = ", ".join(extents)
        return f'{self._working_name(argument)}.attr("dim") = IntegerVector::create({joined});'

    def _new_value(self, argument: CArgument) -> str:
        name = argument.name
        if argument.origin is Origin.ERROR:
            return f"int {name} = 0;"
        if argument.is_scalar:
            return f"{cpp_ctype(argument)} {name} = 0;"

        extents = [cpp_extent(e) for e in argument.dimension.extents]
        if argument.shape_arg is not None:
            # C sees one flat block, so the count is the product of the extents that
            # travel beside it rather than anything the declaration says
            size = (f"(int) std::accumulate({argument.shape_arg}.begin(), "
                    f"{argument.shape_arg}.end(), 1, std::multiplies<int>())")
        else:
            size = _product(extents)
        if argument.conversion is Conversion.LOGICAL:
            return f"tox::BoolBuffer {name}_c({size});"
        if argument.type.is_character:
            length = extents[0]
            # the items are the product of the shape when that travels separately: the
            # declaration only carries `*` in that case
            count = size if argument.shape_arg is not None else _product(extents[1:])
            return f"tox::CharBuffer {name}_c({length}, {count});"
        # a temporary is scratch C++ never shows anyone; an output is an Rcpp vector
        if argument.is_temporary:
            return f"std::vector<{cpp_ctype(argument)}> {name}({size});"
        return f"{RCPP_VECTOR[argument.type.base]} {name}({size});"

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
            # null pointer when absent, so Fortran's OPTIONAL dummy is absent too
            if self._converts_via_buffer(argument):
                return f"{name}.isNotNull() ? {name}_c.data() : nullptr"
            return f"{name}_p"
        if self._converts_via_buffer(argument):
            return f"{name}_c.data()"
        if argument.is_scalar:
            return f"&{name}"
        if argument.is_temporary:
            return f"{name}.data()"
        # an Rcpp vector/matrix: its buffer is already the C layout
        return self._as_c_pointer(argument, f"{self._working_name(argument)}.begin()")

    @staticmethod
    def _as_c_pointer(argument: CArgument, buffer: str) -> str:
        """The Rcpp buffer as the pointer the `extern "C"` declaration asks for.

        Only complex needs saying: a `ComplexVector` iterates as `Rcomplex*`, the pair of
        doubles R stores, which is layout-compatible with C's `double _Complex` but not
        convertible to it. Every other type's iterator already *is* the C pointer.
        """
        if argument.type.base is not BaseType.COMPLEX:
            return buffer
        const = "const " if argument.intent is Intent.IN else ""
        return f"reinterpret_cast<{const}{cpp_ctype(argument)}*>({buffer})"

    def _marshal_outputs(self, writer: Writer, wrapper: CWrapper) -> None:
        lines = Writer()
        for argument in wrapper:
            if argument.intent is Intent.IN or not self._converts_via_buffer(argument):
                continue
            name = argument.name
            rcpp = RCPP_VECTOR[argument.type.base]
            lines.line(f"{rcpp} {name} = {name}_c.to_r();")
            dim = self._dim_attribute(argument)
            if dim:
                lines.line(dim)
        if lines:
            writer.line("// convert the outputs back")
            writer.extend(lines)
            writer.blank()

    def _return(self, writer: Writer, wrapper: CWrapper) -> None:
        outputs = [a for a in wrapper if self._is_returned(a)]
        error = wrapper.error_argument

        writer.line("return List::create(")
        with writer.indent():
            entries = [f'_["{a.name}"] = {self._working_name(a)}' for a in outputs]
            entries.append(f'_["{error.name}"] = {error.name}')
            for index, entry in enumerate(entries):
                comma = "" if index == len(entries) - 1 else ","
                writer.line(f"{entry}{comma}")
        writer.line(");")

    @staticmethod
    def _is_returned(argument: CArgument) -> bool:
        """Every output the R wrapper might use, ierr aside (added separately).

        `intent(inout)` is returned too, unlike in Python. R is copy-on-modify, so it
        cannot hand the modification back through the argument the way Python does -- the
        R wrapper duplicates the input, C++ modifies the copy, and it comes back in the
        list. See `design/language-layers.md`.
        """
        return (
            argument.intent.is_output
            and not argument.is_error
            and not argument.is_temporary
        )


_MARSHAL_HEADER = r'''// Generated. Do not edit.
//
// Marshalling helpers for the R interface. C++ converts only what C cannot take from R
// directly: R's int-based logicals, and its strings. The R layer has already validated
// and rejected NA, so these are straight copies.
#ifndef TOX_MARSHAL_H
#define TOX_MARSHAL_H

#include <Rcpp.h>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace tox {

// R logical is an int (0/1/NA); Fortran c_bool is a 1-byte _Bool, which C++ bool matches.
class BoolBuffer {
    std::unique_ptr<bool[]> buf_;
    std::size_t n_;
public:
    explicit BoolBuffer(std::size_t n) : buf_(new bool[n]()), n_(n) {}
    explicit BoolBuffer(Rcpp::LogicalVector x)
        : buf_(new bool[x.size()]), n_(x.size()) {
        for (std::size_t i = 0; i < n_; ++i) buf_[i] = (x[i] == TRUE);
    }
    bool* data() { return buf_.get(); }
    Rcpp::LogicalVector to_r() const {
        // R stores a logical as an int that is exactly 0 or 1, and compares them with
        // identical(). Fortran's logical(c_bool) only promises "non-zero is true": ifx
        // writes 0xFF, which widens to 255 and prints as TRUE while comparing unequal to
        // TRUE. Test the byte rather than widening it.
        // and read the byte as unsigned char, not as bool: a C++ bool holding anything
        // other than 0 or 1 is undefined behaviour, so the compiler is entitled to fold
        // `buf_[i] ? TRUE : FALSE` down to a plain widening -- which is exactly what it
        // does, leaving the 255 in place.
        const unsigned char* raw = reinterpret_cast<const unsigned char*>(buf_.get());
        Rcpp::LogicalVector out(n_);
        for (std::size_t i = 0; i < n_; ++i) out[i] = raw[i] != 0 ? TRUE : FALSE;
        return out;
    }
};

// Fortran carries a string's length as the leading extent: a vector of n strings of
// length len is char[len * n], column-major, each string zero-padded. Reading stops at
// the first null, so an untouched buffer (zero-filled) yields empty strings, never noise.
// The width a character(len=n) array needs: the longest element. Every string is stored
// in the same n bytes, so anything shorter is null-padded and anything longer would be
// cut off.
inline int max_strlen(Rcpp::CharacterVector x) {
    int longest = 0;
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] == NA_STRING) continue;
        int n = Rf_length(STRING_ELT(x, i));
        if (n > longest) longest = n;
    }
    return longest;
}

class CharBuffer {
    std::vector<char> data_;
    int len_;
    int n_;
public:
    CharBuffer(int len, int n) : data_((std::size_t)len * n, '\0'), len_(len), n_(n) {}
    CharBuffer(Rcpp::CharacterVector x, int len)
        : data_((std::size_t)len * x.size(), '\0'), len_(len), n_(x.size()) {
        for (int i = 0; i < n_; ++i) {
            std::string s = Rcpp::as<std::string>(x[i]);
            int m = std::min<int>((int)s.size(), len_);
            std::memcpy(data_.data() + (std::size_t)i * len_, s.data(), m);
        }
    }
    char* data() { return data_.data(); }
    Rcpp::CharacterVector to_r() const {
        Rcpp::CharacterVector out(n_);
        for (int i = 0; i < n_; ++i) {
            const char* p = data_.data() + (std::size_t)i * len_;
            int m = 0;
            while (m < len_ && p[m] != '\0') ++m;
            out[i] = std::string(p, m);
        }
        return out;
    }
};

}  // namespace tox

#endif
'''
