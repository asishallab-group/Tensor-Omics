"""`DM_EQUALS(EXPR)`: an integer that must equal what other arguments compute (#203, rev 4).

The rule: only on a required `intent(in)` scalar integer; the check is emitted like `DM_MIN`/`DM_MAX`,
after the other scalar checks and before the array checks, under the one return; an extent keeps
its dimension check; the bindings document it on the array a hidden extent is read off (axis 0-
based in Python, 1-based in R) and state `EXPR` in the error. The parsing is
`test_ir_directives.TestEquals`; this is everything after it.
"""

import importlib.util

import pytest

from codegen.abi.c_abi import build_project
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.doc_contract import equality_notes
from codegen.emit.error_arguments import equalities_of
from codegen.emit.errors_python import PythonErrorEmitter
from codegen.emit.fortran_wrapper import FortranWrapperEmitter
from codegen.generate import generate
from codegen.ir.directives import Directives, Equals, Maximum, Minimum
from codegen.ir.entities import Meta
from codegen.ir.errors import ErrorCatalogue
from codegen.ir.roles import analyse_project
from codegen.ir.types import Intent
from codegen.ir.validate import validate_project
from codegen.synthesize import synthesize_wrappers

import builders as b
from conftest import REPO_ROOT

IMPL_META = Meta(summary="Sum each group", author="AUTHOR")


def equals(expression, **others):
    return Directives(equals=Equals(expression), **others)


def grouped(n_rows=None, *extra):
    """`grouped_rows_impl` of the fixture: `n_rows` must equal `sum(group_sizes)`."""
    n_rows = n_rows or b.integer("n_rows", directives=equals("sum(group_sizes)"), doc="rows")
    return (
        n_rows,
        b.integer("n_groups", doc="groups"),
        b.integer("group_sizes", Intent.IN, "(n_groups)", doc="sizes",
                  directives=Directives(minimum=Minimum("1_int32"))),
        b.real("values", Intent.IN, "(n_rows)", doc="values"),
        b.real("sums", Intent.OUT, "(n_groups)", doc="sums"),
        *extra,
    )


def impl_project(*arguments):
    impl = b.procedure("grouped_rows_impl", *arguments, meta=IMPL_META)
    return b.project(b.module("tox_demo_impl", impl, doc="demo", path="src/tox/tox_demo_impl.F90"))


def synthesized(project):
    bag = DiagnosticBag()
    synthesis = synthesize_wrappers(project, CONVENTIONS, bag)
    analyse_project(synthesis.project, bag, CONVENTIONS)
    return bag, synthesis


def c_wrapper(project, name):
    bag, synthesis = synthesized(project)
    binding = build_project(synthesis.project, bag, CONVENTIONS)
    return next(w for module in binding for w in module if w.stripped_name == name)


def refusals(*arguments):
    bag, synthesis = synthesized(impl_project(*arguments))
    validate_project(synthesis.project, bag, CONVENTIONS)
    return [d.message for d in bag.errors if "DM_EQUALS" in d.message]


# -- rule 1: where it may go -------------------------------------------------------------------


REFUSED = "argument '{}' has a DM_EQUALS but is not a required intent(in) scalar integer"


class TestWhereItGoes:
    def test_on_a_required_input_integer_it_is_accepted(self):
        bag, synthesis = synthesized(impl_project(*grouped()))
        validate_project(synthesis.project, bag, CONVENTIONS)

        assert bag.errors == (), bag.render()

    def test_on_an_integer_that_sizes_nothing_it_is_accepted(self):
        assert refusals(*grouped(None, b.integer("version", directives=equals("3_int32"),
                                                  doc="version"))) == []

    @pytest.mark.parametrize("argument", [
        b.real("scale", directives=equals("sum(group_sizes)"), doc="scale"),
        b.logical("flag", kind="c_bool", directives=equals("n_groups > 1"), doc="flag"),
        b.integer("count", Intent.OUT, directives=equals("n_groups"), doc="count"),
        b.integer("count", Intent.INOUT, directives=equals("n_groups"), doc="count"),
        b.integer("count", directives=equals("n_groups"), doc="count", optional=True),
        b.integer("labels", Intent.IN, "(n_groups)", directives=equals("n_groups"), doc="labels"),
    ], ids=["real", "logical", "intent(out)", "intent(inout)", "optional", "array"])
    def test_anything_else_is_refused_with_one_message(self, argument):
        assert refusals(*grouped(None, argument)) == [REFUSED.format(argument.name)]

    def test_the_real_pipeline_breaks_when_it_is_misplaced(self, tmp_path):
        # `verify-by-breaking-it`: the check is reached from real source, not only hand-built IR
        source = (REPO_ROOT / "src/tox/tox_normalization_impl.F90").read_text()
        broken = source.replace("        integer(int32), intent(in) :: n_replicates\n",
                                "        integer(int32), intent(in), optional :: n_replicates\n")
        assert broken != source  # the premise
        src = tmp_path / "src" / "tox"
        src.mkdir(parents=True)
        (src / "tox_normalization_impl.F90").write_text(broken)
        (src / "tox_loess_impl.F90").write_text((REPO_ROOT / "src/tox/tox_loess_impl.F90").read_text())
        result = generate(paths=Paths(root=REPO_ROOT, src_dir=tmp_path / "src"),
                          targets=("fortran",))

        assert any(d.message == REFUSED.format("n_replicates")
                   for d in result.diagnostics.errors), result.diagnostics.render()


def test_the_expression_names_a_dummy_as_it_is_declared():
    # the error text and the doc quote it, and Python and R compare names exactly
    from codegen.frontend.ford_frontend import _use_declared_spellings

    arguments = [b.integer("n_rows", directives=equals("sum(GROUP_SIZES)")),
                 b.integer("Group_Sizes", Intent.IN, "(n_groups)")]
    _use_declared_spellings(arguments)

    assert arguments[0].directives.equals.expression == "sum(Group_Sizes)"


# -- rules 2 and 3: what the wrapper checks, and where -----------------------------------------


def validation_block(module_project, name="grouped_rows"):
    _, synthesis = synthesized(module_project)
    text = FortranWrapperEmitter(project=synthesis.project).module(
        synthesis.project.module("tox_demo"))
    body = text[text.index(f"subroutine {name}("):text.index(f"end subroutine {name}")]
    block = body[body.index("#ifndef NO_INPUT_VALIDATION"):body.index("#endif")]
    return [line.strip() for line in block.splitlines()[1:]], text


class TestTheCheck:
    def test_after_the_scalar_checks_before_the_arrays_under_the_one_return(self):
        lines, _ = validation_block(impl_project(
            *grouped(),
            b.integer("degree", doc="degree", directives=Directives(maximum=Maximum("2_int32"))),
        ))

        assert lines == [
            # n_rows keeps its dimension check
            "call validate_dimension_size(n_rows, ierr, arg_pos=1_int32)",
            "call validate_dimension_size(n_groups, ierr, arg_pos=2_int32)",
            "call validate_in_range_int(degree, ierr, arg_pos=6_int32, max=2_int32)",
            # after every scalar check, so an empty n_groups is not reported as a mismatch
            "if (n_rows /= (sum(group_sizes))) call set_err_once(ierr, ERR_SIZE_MISMATCH, arg_pos=1_int32)",
            "call validate_all_in_range_int(group_sizes, n_groups, ierr, arg_pos=3_int32, min=1_int32)",
            "call validate_all_in_range_real(values, n_rows, ierr, arg_pos=4_int32)",
            "if (is_err(ierr)) return",
        ]

    def test_a_wrapper_without_one_is_unchanged(self):
        lines, _ = validation_block(impl_project(
            b.integer("n_rows", doc="rows"), *grouped()[1:]))

        assert not any("ERR_SIZE_MISMATCH" in line for line in lines)
        assert lines.count("if (is_err(ierr)) return") == 1

    def test_the_module_imports_what_the_comparison_names(self):
        _, text = validation_block(impl_project(*grouped()))
        header = text.split("contains", 1)[0]

        assert "ERR_SIZE_MISMATCH" in header
        assert "set_err_once" in header

    def test_a_module_constant_in_the_expression_is_imported(self):
        project = impl_project(
            b.integer("n", directives=equals("FX_PER_GROUP * n_groups"), doc="length"),
            b.integer("n_groups", doc="groups"),
            b.integer("labels", Intent.IN, "(n_groups)", doc="labels"),
            b.real("values", Intent.IN, "(n)", doc="values"),
        )
        project.modules[0].parameters = (b.parameter("FX_PER_GROUP", "3_int32"),)
        project.modules[0]._adopt(project.modules[0].parameters)
        _, text = validation_block(project)

        assert "use tox_demo_impl, only: FX_PER_GROUP, grouped_rows_impl" in text


# -- rule 4: what the bindings say --------------------------------------------------------


@pytest.fixture(scope="module")
def fixture_files():
    """Every target generated from the implementation fixtures, by path name."""
    paths = Paths(root=REPO_ROOT, src_dir=REPO_ROOT / "helper/codegen/tests/fixtures/impl")
    result = generate(paths=paths, targets=("fortran", "c", "python", "r"))
    return {file.path.name: file.content for file in result.files}


class TestTheDocumentation:
    def test_fortran_says_it_on_the_argument(self, fixture_files):
        wrapper = fixture_files["fx_extents.F90"]
        start = wrapper.index("subroutine grouped_rows(")
        declaration = wrapper[start:wrapper.index(":: n_groups", start)]

        assert "!! Must equal `sum(group_sizes)`." in declaration

    def test_python_says_it_on_the_array_the_extent_is_read_off(self, fixture_files):
        python = fixture_files["fx_extents.py"]
        docstring = python[python.index("def grouped_rows("):]
        docstring = docstring[:docstring.index("Returns")]

        assert ("    values : np.ndarray[np.float64] of shape (n_rows,)\n"
                "        the values, group after group\n"
                "        `n_rows`, its axis 0, must equal `sum(group_sizes)`.\n") in docstring
        assert "n_rows :" not in docstring  # the premise: the extent is hidden

    def test_r_says_it_on_the_array_the_extent_is_read_off(self, fixture_files):
        r = fixture_files["fx_extents.R"]
        block = r[r.index("@param values", r.index("The sum of each group of rows")):]

        assert block.startswith("@param values a numeric vector. the values, group after group\n"
                                "#'   `n_rows`, its axis 1, must equal `sum(group_sizes)`.\n")

    def test_the_axis_is_counted_as_the_language_counts_it(self):
        # numpy's shape[0] is R's dim(x)[1]; a square array names both of its axes
        wrapper = c_wrapper(impl_project(
            b.integer("n", directives=equals("2_int32 * n_half"), doc="size"),
            b.integer("n_half", doc="half"),
            b.integer("halves", Intent.IN, "(n_half)", doc="halves"),
            b.real("square", Intent.IN, "(n, n)", doc="square"),
        ), "grouped_rows")
        square = next(a for a in wrapper if a.name == "square")

        assert equality_notes(square, "python") == ["`n`, its axes 0 and 1, must equal `2 * n_half`."]
        assert equality_notes(square, "r") == ["`n`, its axes 1 and 2, must equal `2 * n_half`."]

    def test_an_array_with_no_marked_extent_gets_nothing(self, fixture_files):
        python = fixture_files["fx_extents.py"]

        # one sentence per marked extent, on its array: none on `group_sizes` or `sums`
        assert python.count(", its axis ") == 1


@pytest.fixture(scope="module")
def errors_module(tmp_path_factory):
    """The Python error module, with ERR_SIZE_MISMATCH, imported."""
    bag = DiagnosticBag()
    source = b.module(
        "tox_errors",
        parameters=(
            b.parameter("ERR_OK", "0", doc="no error, operation successful"),
            b.parameter("ERR_INVALID_INPUT", "201", doc="invalid input arguments"),
            b.parameter("ERR_SIZE_MISMATCH", "206", doc="Array size mismatch"),
        ),
    )
    catalogue = ErrorCatalogue.from_module(source, bag)
    assert bag.errors == (), bag.render()
    path = tmp_path_factory.mktemp("errors") / "error_handling.py"
    path.write_text(PythonErrorEmitter(catalogue).module())
    spec = importlib.util.spec_from_file_location("error_handling", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestTheErrorText:
    ARGUMENTS = ("n_rows", "n_groups", "group_sizes", "values", "sums", "ierr")
    SOURCES = ("values", "group_sizes", None, None, None, None)
    EQUALS = (None, None, None, None, None, None)

    def raised(self, module, ierr, equals):
        with pytest.raises(module.ToxError) as caught:
            module.check_err_code(ierr, self.ARGUMENTS, self.SOURCES, equals=equals)
        return caught.value

    def test_python_states_what_the_argument_had_to_equal(self, errors_module):
        equals = ("sum(group_sizes)",) + (None,) * 5
        error = self.raised(errors_module, 206 + 1 * errors_module.ARG_POS_FACTOR, equals)

        assert str(error) == ("Array size mismatch (argument 'values', via 'n_rows'); "
                              "'n_rows' must equal sum(group_sizes)")
        assert error.argument == "values"

    def test_only_for_a_size_mismatch(self, errors_module):
        equals = ("sum(group_sizes)",) + (None,) * 5
        error = self.raised(errors_module, 201 + 1 * errors_module.ARG_POS_FACTOR, equals)

        assert "must equal" not in str(error)

    def test_only_at_the_argument_that_has_one(self, errors_module):
        equals = ("sum(group_sizes)",) + (None,) * 5
        error = self.raised(errors_module, 206 + 3 * errors_module.ARG_POS_FACTOR, equals)

        assert str(error) == "Array size mismatch (argument 'group_sizes')"

    def test_a_kind_suffix_is_written_as_the_language_writes_it(self):
        procedure = c_wrapper(impl_project(
            b.integer("n", directives=equals("2_int32 * n_half"), doc="size"),
            b.integer("n_half", doc="half"),
            b.integer("halves", Intent.IN, "(n_half)", doc="halves"),
            b.real("values", Intent.IN, "(n)", doc="values"),
        ), "grouped_rows").procedure

        assert equalities_of(procedure, "python")[0] == "2 * n_half"
        assert equalities_of(procedure, "r")[0] == "2 * n_half"

    def test_the_generated_tables(self, fixture_files):
        python = fixture_files["fx_extents.py"]
        r = fixture_files["fx_extents.R"]

        assert ('_GROUPED_ROWS_ARGUMENT_EQUALS = ("sum(group_sizes)", None, None, None, None, None,)'
                in python)
        assert "equals=_GROUPED_ROWS_ARGUMENT_EQUALS)" in python
        assert ('.equals <- c("sum(group_sizes)", NA_character_, NA_character_, NA_character_, '
                'NA_character_, NA_character_)') in r
        assert "check_err_code(.result$ierr, .arguments, .sources, equals = .equals)" in r
