"""The input contract: what each binding emits to check an argument, and where (#221 WP6).

The behaviour -- every case of the table, accepted and refused -- is tested on the library in
`python/test/mod_test_input_contract.py` and `r/test/mod_test_input_contract.R`. What is here
is what those cannot see: that the check sits where the design puts it (R's element checks in
the C shim, none in the generated R), that a refusal names its argument in every layer, that
the shim's refusals are the R wrapper's own conditions, and the Python engine's edges, run
directly on the generated module.
"""

import importlib
import re
import sys

import pytest

from codegen import names
from codegen.abi.c_abi import build_project, build_wrapper
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.c_call import CCallEmitter, _checks_integer_array
from codegen.emit.doc_contract import TRAILING_BLANKS
from codegen.emit.errors_python import PythonErrorEmitter
from codegen.emit.errors_r import ARGUMENT_CONDITION_CLASSES, CONDITION_FIELDS, NA_MESSAGE, RErrorEmitter
from codegen.emit.python_contract import contract_module
from codegen.emit.python_ctypes import PythonEmitter
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.errors import ErrorCatalogue
from codegen.ir.roles import analyse, analyse_project
from codegen.ir.types import BaseType, Intent
from codegen.ir.validate import validate_project
from codegen.synthesize import synthesize_wrappers
from conftest import REPO_ROOT

import builders as b


def wrapper_of(procedure):
    bag = DiagnosticBag()
    analyse(procedure, bag)
    wrapper = build_wrapper(procedure, bag)
    assert bag.errors == (), bag.render()
    return wrapper


def numbers():
    return wrapper_of(b.procedure(
        "p",
        b.integer("n", Intent.IN, doc="how many"),
        b.integer("k", Intent.IN, doc="a single integer"),
        b.integer("labels", Intent.IN, "(n)", doc="a label per element"),
        b.integer("counts", Intent.INOUT, "(n)", doc="counted up"),
        b.logical("keep", Intent.IN, "(n)", doc="which to keep"),
        b.integer("picked", Intent.IN, "(n)", optional=True, doc="the ones picked"),
        b.character("tag", Intent.IN, length="8", doc="a tag"),
        b.character("names_", Intent.IN, "(n)", doc="a name per element"),
        b.integer("found", Intent.OUT, "(n)", doc="what was found"),
        b.ierr(),
    ))


@pytest.fixture(scope="module")
def real_binding():
    """The binding of the real sources, as `generate` builds it."""
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT), bag).parse()
    synthesis = synthesize_wrappers(parsed.project, CONVENTIONS, bag)
    analyse_project(synthesis.project, bag)
    validate_project(synthesis.project, bag)
    binding = build_project(synthesis.project, bag)
    assert bag.errors == (), bag.render()
    return binding


# -- R: the kind and rank in R, every element in the C shim -------------------------------------

class TestTheShimChecksTheElements:
    def text(self):
        return CCallEmitter().function(numbers())

    def test_an_integer_array_is_converted_and_checked_by_name(self):
        assert 'const int* _c_labels = tox_int_in(labels, "labels");' in self.text()
        assert "_c_labels," in self.text()
        assert "INTEGER(labels)" not in self.text()

    def test_an_optional_integer_array_is_too_and_passed_as_null_when_absent(self):
        text = self.text()
        assert 'const int* _c_picked = tox_int_in(picked, "picked");' in text
        assert "picked != R_NilValue ? _c_picked : NULL" in text

    def test_an_integer_array_modified_in_place_is_checked_into_its_copy(self):
        assert 'SEXP _out_counts = PROTECT(tox_int_copy(counts, "counts"));' in self.text()

    def test_logicals_and_strings_name_their_argument(self):
        text = self.text()
        assert 'unsigned char* _c_keep = tox_bool_in(keep, "keep");' in text
        assert 'char* _c_tag = tox_char_in(tag, 8, "tag");' in text
        assert 'tox_char_in(names_, names__strlen, "names_")' in text

    def test_an_output_and_a_single_value_are_not_converted(self):
        # an output is allocated as an integer vector; a single value is R's to check, as R
        # uses it -- it reaches the shim an integer already
        text = self.text()
        assert "INTEGER(found)" in text
        assert "_v_k = Rf_asInteger(k);" in text
        assert "tox_int_in(k," not in text and "tox_int_in(found" not in text

    def test_what_the_binding_derives_is_never_converted(self, real_binding):
        # a shape is used as an R vector by the shim, an extent or a count is its own local
        for module in real_binding:
            for wrapper in module:
                for argument in wrapper:
                    roles = argument.source.roles if argument.source else None
                    if roles and (roles.is_derived or roles.is_computed or roles.is_shape_arg):
                        assert not _checks_integer_array(argument), (wrapper.name, argument.name)

    def test_every_integer_logical_and_string_array_the_caller_passes_is_checked(self, real_binding):
        emitter = CCallEmitter()
        r = RWrapperEmitter()
        for module in real_binding:
            for wrapper in module:
                text = emitter.function(wrapper.with_constant_sizes_evaluated())
                for argument in r._inputs(wrapper):
                    roles = argument.source.roles if argument.source else None
                    if roles and roles.is_shape_arg:
                        continue
                    base = argument.type.base
                    if base is BaseType.LOGICAL and argument.is_array:
                        helper = "tox_bool_in"
                    elif base is BaseType.INTEGER and argument.is_array:
                        helper = "tox_int_in"
                    elif base is BaseType.CHARACTER:
                        helper = "tox_char_in"
                    else:
                        continue
                    assert re.search(rf'{helper}\({argument.name}\b.*"{argument.name}"\);', text), \
                        (wrapper.name, argument.name)


#: An R expression that reads every element: what the generated R must not do with an integer,
#: logical or character array (A6.5), the C shim doing it in its conversion pass instead
_ELEMENT_PASSES = ("anyNA(", "as.integer(", "as.logical(", 'storage.mode(x) <- "integer"',
                   "enc2utf8(", "validUTF8(", "is.na(x)")

#: The R helpers that take an array of integers, logicals or strings
_ARRAY_HELPERS = tuple(
    f".tox_as_{base}_{shape}"
    for base in ("integer", "logical") for shape in ("vector", "matrix", "array", "shaped")
) + (".tox_as_character", ".tox_as_character_shaped", ".tox_as_string")


def _r_functions(text):
    """Each top-level `name <- function(...) {...}` of an R file, as name -> body."""
    text = re.sub(r"^#.*\n", "", text, flags=re.M)
    starts = [(m.start(), m.group(1)) for m in re.finditer(r"^(\.\w+) <- function", text, re.M)]
    bodies = {}
    for index, (start, name) in enumerate(starts):
        end = starts[index + 1][0] if index + 1 < len(starts) else len(text)
        bodies[name] = text[start:end]
    return bodies


class TestTheGeneratedRScansNoArray:
    def test_no_helper_for_an_array_reads_its_elements(self):
        functions = _r_functions(RWrapperEmitter().validators())
        for helper in _ARRAY_HELPERS:
            body = functions[helper]
            for scan in _ELEMENT_PASSES:
                assert scan not in body, (helper, scan)

    def test_no_wrapper_reads_an_element_either(self, real_binding):
        emitter = RWrapperEmitter()
        for module in real_binding:
            text = emitter.module(module)
            assert "anyNA(" not in text and "as.integer(" not in text, module.name

    def test_a_single_value_is_checked_in_r_which_uses_it(self):
        functions = _r_functions(RWrapperEmitter().validators())
        assert ".tox_check_whole(x, name)" in functions[".tox_as_integer_scalar"]
        assert "as.integer(x)" in functions[".tox_as_integer_scalar"]
        assert "as.logical(x)" in functions[".tox_as_logical_scalar"]

    def test_a_kind_is_checked_by_predicate_never_by_typeof(self):
        validators = RWrapperEmitter().validators()
        assert "typeof(" not in validators and "storage.mode(x) ==" not in validators
        assert 'inherits(x, "integer64")' in validators

    def test_a_single_string_goes_through_the_string_helper(self):
        # it went through .tox_as_character, which took a vector and passed its first element
        text = RWrapperEmitter().function(numbers())
        assert 'tag <- .tox_as_string(tag, "tag")' in text
        assert 'names_ <- .tox_as_character(names_, "names_")' in text

    def test_a_shape_the_caller_passes_is_checked_and_made_integer_in_r(self):
        # the shim uses it as an R vector -- the dim of what it allocates -- so it must reach
        # the shim an integer vector, not a double one for the shim to convert
        wrapper = wrapper_of(b.procedure(
            "p",
            b.real("arr", Intent.OUT, "(*)", doc="the values"),
            b.integer("arr_shape", Intent.IN, "(:)", doc="the extents of `arr`"),
            b.ierr(),
        ))
        assert 'arr_shape <- .tox_as_integer_shape(arr_shape, "arr_shape")' in \
            RWrapperEmitter().function(wrapper)
        shim = CCallEmitter().function(wrapper)
        assert "tox_int_in(arr_shape" not in shim and "tox_prod(arr_shape)" in shim
        body = _r_functions(RWrapperEmitter().validators())[".tox_as_integer_shape"]
        assert ".tox_check_whole(value, name)" in body and "as.integer(x)" in body

    def test_a_mask_is_counted_as_the_shim_counts_it(self):
        # an unchecked mask can hold NA or a 2; the count stays a whole number regardless,
        # which the shape error prints with %d
        text = RWrapperEmitter().function(wrapper_of(b.procedure(
            "p",
            b.integer("n", Intent.IN, doc="how many"),
            b.logical("genes_selection_mask", Intent.IN, "(n)", doc="which genes"),
            b.integer("n_selected_genes", Intent.IN, doc="how many the mask selects"),
            b.real("weights", Intent.IN, "(n_selected_genes)", doc="one per selected gene"),
            b.real("total", Intent.OUT, doc="the sum"),
            b.ierr(),
        )))
        assert "sum(genes_selection_mask == 1, na.rm = TRUE)" in text
        assert "sum(genes_selection_mask)" not in text


class TestTheShimRefusesAsRDoes:
    def header(self):
        return CCallEmitter().marshal_header_content()

    def test_its_classes_fields_and_na_message_are_the_r_wrappers(self):
        header = self.header()
        assert ", ".join(f'"{name}"' for name in ARGUMENT_CONDITION_CLASSES) in header
        assert ", ".join(f'"{name}"' for name in CONDITION_FIELDS) in header
        assert f'"{NA_MESSAGE}"' in header
        assert not re.search(r"__[A-Z_]+__", header), "a placeholder left unfilled"

    def test_r_raises_with_the_same_fields_in_the_same_order(self):
        text = _r_module_text()
        raised = re.search(r"list\((.*?)\)\n", text, re.S).group(1)
        fields = re.findall(r"(\w+) = ", raised)
        assert tuple(fields) == CONDITION_FIELDS

    def test_every_function_the_header_defines_is_in_the_names_table(self):
        # the table is what validate refuses an argument by, so it must be exactly the header's
        defined = set(re.findall(r"^static inline [\w ]+?\**\s*\b(tox_\w+)\(", self.header(), re.M))
        assert defined and defined == set(names.R_SHIM_HELPERS)

    def test_an_argument_named_like_one_would_hide_it_in_the_shim(self):
        # the premise of refusing it: the shim's parameter list carries the argument's name
        text = CCallEmitter().function(wrapper_of(b.procedure(
            "p", b.integer("tox_int_in", Intent.IN, "(n)", doc="x"), b.integer("n", Intent.IN, doc="n"),
            b.ierr())))
        assert "SEXP tox_int_in" in text and "tox_int_in(tox_int_in," in text
        assert names.language_word("tox_int_in") == "a function of the R binding's C shim"

    def test_a_string_returned_to_r_is_marked_utf8(self):
        assert "Rf_mkCharLenCE(p, m, CE_UTF8)" in self.header()

    def test_a_string_is_translated_and_validated_before_it_is_copied(self):
        header = self.header()
        body = header[header.index("static inline char* tox_char_in"):]
        body = body[:body.index("\n}\n")]
        assert body.index("tox_utf8(") < body.index("tox_refuse(") < body.index("memcpy(")
        assert "a longer value would be cut short" in body

    def test_only_a_latin1_string_is_translated_and_a_native_one_never(self):
        # translating a native string in a non-UTF-8 locale wrote its bytes out as "<c3><a9>"
        # escape text, which then passed as valid ASCII (#221 WP8); the run-time cases, in both
        # locales, are r/test/mod_test_input_contract.R's
        header = self.header()
        body = header[header.index("static inline const char* tox_utf8("):]
        body = body[:body.index("\n}\n")]
        assert body.count("Rf_translateCharUTF8(") == 1
        latin1 = body[body.index("CE_LATIN1"):]
        assert latin1.index("Rf_translateCharUTF8(") < latin1.index("return translated")
        assert body.rindex("tox_utf8_valid(") > body.index("return translated")
        assert "Rf_translateCharUTF8(" not in header.replace(body, "")


def _r_module_text():
    """The R error module, built from the real tox_errors."""
    bag = DiagnosticBag()
    real = FordFrontend(Paths(root=REPO_ROOT, src_dir=Paths().src_dir / "f42"), bag).parse()
    catalogue = ErrorCatalogue.from_module(real.project.module("tox_errors"), bag,
                                           real.project.constant_values(),
                                           arg_pos_factor=real.arg_pos_factor)
    return RErrorEmitter(catalogue).module()


# -- Python: every check through the one module -------------------------------------------------

class TestPythonChecksThroughTheContract:
    def text(self):
        return PythonEmitter().function(numbers())

    def test_each_input_is_checked_by_its_kind(self):
        text = self.text()
        assert f'k = {names.PYTHON_CONTRACT}.integer(k, "k")' in text
        assert f'labels = {names.PYTHON_CONTRACT}.integers(labels, "labels", ndim=1)' in text
        assert f'keep = {names.PYTHON_CONTRACT}.logicals(keep, "keep", ndim=1)' in text
        assert f'picked = {names.PYTHON_CONTRACT}.integers(picked, "picked", ndim=1)' in text

    def test_a_string_of_a_declared_length_passes_the_length_to_check(self):
        assert f'{names.PYTHON_CONTRACT}.string(tag, "tag", 8).ljust(8)' in self.text()
        assert f'{names.PYTHON_CONTRACT}.strings(names_, "names_")' in self.text()

    def test_an_inout_array_is_still_refused_rather_than_converted(self):
        assert "'counts' is modified in place" in self.text()

    def test_the_contract_is_bound_under_an_underscored_name(self):
        assert names.PYTHON_CONTRACT.startswith("_") and names.PYTHON_CONTRACT in names.PYTHON_BOUND


class TestTheDocsStateTheStringRules:
    def test_a_string_input_says_its_trailing_blanks_are_dropped(self):
        python = PythonEmitter().function(numbers())
        r = RWrapperEmitter().function(numbers())
        assert python.count(TRAILING_BLANKS) == 2 and r.count(TRAILING_BLANKS) == 2

    def test_a_mode_says_its_trailing_blanks_are_dropped_too(self, real_binding):
        wrapper = next(w for m in real_binding for w in m if w.stripped_name == "linkage_clustering")
        assert TRAILING_BLANKS in PythonEmitter().function(wrapper)
        assert TRAILING_BLANKS in RWrapperEmitter().function(wrapper)

    @pytest.mark.parametrize("procedure, line_before", [
        ("linkage_clustering", "The maximum valid value is `2`."),
        ("compute_family_scaling", "The default value is `'robust'`."),
    ])
    def test_a_mode_s_sentence_follows_its_text_as_a_string_s_does(self, real_binding,
                                                                    procedure, line_before):
        # the skipped mode table left a blank line, which made the sentence a paragraph
        wrapper = next(w for m in real_binding for w in m if w.stripped_name == procedure)
        lines = [line.strip() for line in PythonEmitter().function(wrapper).splitlines()]
        assert lines[lines.index(TRAILING_BLANKS) - 1] == line_before

    def test_a_declared_length_is_stated(self):
        python = PythonEmitter().function(numbers())
        assert "At most 8 bytes long in UTF-8: a longer value is refused, not cut short." in python

    def test_an_r_array_of_strings_says_a_factor_is_taken_by_its_labels(self):
        assert "A factor is taken by its labels." in RWrapperEmitter().function(numbers())


# -- Python: a complex output is two doubles ----------------------------------------------------

class TestAComplexOutputIsTwoDoubles:
    """#221 WP8: a `c_double` box held half of the 16 bytes Fortran writes into it."""

    def text(self):
        return PythonEmitter().function(wrapper_of(b.procedure(
            "p",
            b.integer("n", Intent.IN, doc="how many"),
            b.complex_("z", Intent.OUT, doc="one complex number"),
            b.complex_("zs", Intent.OUT, "(n)", doc="n complex numbers"),
            b.ierr(),
        )))

    def test_a_scalar_is_boxed_in_two_doubles_passed_as_a_pointer_and_returned_as_complex(self):
        text = self.text()
        assert f"z = ({names.PYTHON_CTYPES}.c_double * 2)()" in text
        assert "        z,\n" in text, "the box is passed as it is, which ctypes reads as c_double*"
        assert '"z": complex(z[0], z[1]),' in text

    def test_an_array_is_numpy_s_complex128_as_before(self):
        text = self.text()
        assert f"zs = {names.PYTHON_NUMPY}.empty((n,), dtype={names.PYTHON_NUMPY}.complex128" in text
        assert '"zs": zs,' in text


# -- the Python engine, run on the generated module ----------------------------------------------

@pytest.fixture(scope="module")
def contract(tmp_path_factory):
    np = pytest.importorskip("numpy")
    package = tmp_path_factory.mktemp("contract") / "contract_package"
    package.mkdir()
    (package / "__init__.py").write_text("")
    bag = DiagnosticBag()
    real = FordFrontend(Paths(root=REPO_ROOT, src_dir=Paths().src_dir / "f42"), bag).parse()
    catalogue = ErrorCatalogue.from_module(real.project.module("tox_errors"), bag,
                                           real.project.constant_values(),
                                           arg_pos_factor=real.arg_pos_factor)
    (package / "error_handling.py").write_text(PythonErrorEmitter(catalogue).module())
    (package / f"{names.PYTHON_CONTRACT}.py").write_text(contract_module())
    sys.path.insert(0, str(package.parent))
    try:
        module = importlib.import_module(f"contract_package.{names.PYTHON_CONTRACT}")
    finally:
        sys.path.remove(str(package.parent))
    return module, np


def refuses(function, pattern):
    with pytest.raises(TypeError, match=pattern) as raised:
        function()
    assert raised.value.code is None and raised.value.argument == "x"
    assert isinstance(raised.value, ValueError)


class TestTheEngine:
    def test_an_input_of_the_fortran_type_is_not_copied(self, contract):
        c, np = contract
        for dtype, check in ((np.int32, c.integers), (np.float64, c.reals), (np.bool_, c.logicals)):
            array = np.zeros(5, dtype=dtype)
            assert check(array, "x", ndim=1) is array
            matrix = np.zeros((2, 3), dtype=dtype, order="F")
            assert check(matrix, "x", ndim=2, order="F") is matrix

    def test_the_int32_range_is_exactly_numpys(self, contract):
        c, np = contract
        assert c.integer(-2**31, "x") == -2**31 and c.integer(2**31 - 1, "x") == 2**31 - 1
        refuses(lambda: c.integer(2**31, "x"), "2147483648")
        refuses(lambda: c.integer(-2**31 - 1, "x"), "-2147483649")
        refuses(lambda: c.integers(np.float64([2.0**31]), "x", ndim=1), "element 0 is 2147483648.0")
        refuses(lambda: c.integers(np.uint64([2**63]), "x", ndim=1), "9223372036854775808")

    def test_a_float_is_an_integer_only_if_whole_finite_and_in_range(self, contract):
        c, np = contract
        assert c.integers([1.0, -0.0, 3e9 - 3e9], "x", ndim=1).tolist() == [1, 0, 0]
        for bad in (0.5, np.nan, np.inf, -np.inf, 1e10):
            refuses(lambda: c.integers([1.0, bad], "x", ndim=1), "element 1 is")

    def test_the_first_refused_element_is_named_by_its_index(self, contract):
        c, np = contract
        refuses(lambda: c.integers(np.array([[1.0, 2.0], [3.0, 3.5]]), "x", ndim=2), r"element \(1, 1\) is 3.5")

    def test_a_logical_is_a_bool_or_exactly_0_or_1(self, contract):
        c, np = contract
        assert c.logicals([1, 0, 1.0, -0.0, True], "x", ndim=1).tolist() == [True, False, True, False, True]
        assert c.logical(np.int8(1), "x") is True
        for bad in (2, -1, 0.5, np.nan):
            refuses(lambda: c.logicals([0, bad], "x", ndim=1), "element 1 is")

    def test_kinds_outside_the_allowlist_are_refused_before_converting(self, contract):
        c, np = contract
        for value in (["1"], [b"1"], [None], [1j], np.array(["2020-01-01"], dtype="datetime64[D]"),
                      np.array([1], dtype="timedelta64[s]"), [True]):
            refuses(lambda: c.integers(value, "x", ndim=1), "must be an array of whole numbers")
            refuses(lambda: c.reals(value, "x", ndim=1), "must be an array of numbers")

    def test_a_scalar_is_the_0d_case(self, contract):
        c, np = contract
        assert c.real(np.array(2.5), "x") == 2.5 and isinstance(c.real(1, "x"), float)
        refuses(lambda: c.real([2.5], "x"), r"single value, not an array of shape \(1,\)")

    def test_a_string_is_checked_whole(self, contract):
        c, np = contract
        assert c.string("ab ", "x", 3) == b"ab "
        assert c.string("é", "x") == "é".encode()
        refuses(lambda: c.string("é", "x", 1), "2 bytes long in UTF-8, but it is declared 1 long")
        refuses(lambda: c.strings(["a", "b\x00"], "x"), r"NUL character \(element 1\)")

    def test_a_mode_lists_the_modes(self, contract):
        c, np = contract
        assert c.mode("Ward", "x", ("ward", "single")) == b"ward"
        assert c.mode("ward  ", "x", ("ward", "single")) == b"ward"
        refuses(lambda: c.mode(" ward", "x", ("ward", "single")), "not ' ward'")
        refuses(lambda: c.mode("wardx", "x", ("ward", "single")), "one of 'ward', 'single', not 'wardx'")
