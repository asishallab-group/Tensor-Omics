"""The Fortran C wrapper emitter.

Most tests assert on the emitted text from hand-built IR, which is fast. The last class
compiles the wrappers generated from the fixture modules with gfortran, which is the only
thing that can actually prove the output is valid Fortran.
"""

import re
import shutil
import subprocess

import pytest

from codegen.abi.c_abi import build_module, build_wrapper
from codegen.diagnostics import DiagnosticBag
from codegen.emit.fortran_c import FortranCEmitter
from codegen.ir.roles import analyse
from codegen.ir.types import Intent

import builders as b
from conftest import REPO_ROOT


@pytest.fixture
def bag():
    return DiagnosticBag()


@pytest.fixture
def emitter():
    return FortranCEmitter()


def emit(procedure, bag, emitter):
    analyse(procedure, bag)
    return emitter.wrapper(build_wrapper(procedure, bag))


def emit_module(module, bag, emitter):
    for procedure in module.exported_procedures:
        analyse(procedure, bag)
    return emitter.module(build_module(module, bag))


def normalize(procedure_name="fx_normalize"):
    return b.procedure(
        procedure_name,
        b.integer("n_dims", Intent.IN, doc="number of elements in `vector`"),
        b.real("vector", Intent.INOUT, "(n_dims)", doc="the vector"),
        b.ierr(),
    )


class TestSignature:
    def test_the_symbol_is_bound_with_its_c_name(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        assert "subroutine fx_normalize_c(&" in text
        assert 'bind(C, name="fx_normalize_c")' in text

    def test_the_wrapped_procedure_is_imported_by_name(self, bag, emitter):
        assert "use fx_basics, only: fx_normalize" in emit_module(
            b.module("fx_basics", normalize()), bag, emitter
        )

    def test_every_argument_appears_in_the_list(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        arguments = re.search(r"subroutine fx_normalize_c\(&\n(.*?)\n\s*\) bind", text, re.S)
        assert [a.strip().rstrip(",&").rstrip("&") for a in arguments.group(1).split("\n")] == [
            "n_dims",
            "vector",
            "ierr",
        ]


class TestDeclarations:
    def test_arguments_are_declared_with_target_for_c_loc(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        assert "integer(c_int), intent(in), target :: n_dims" in text
        assert "real(c_double), dimension(n_dims), intent(inout), target :: vector" in text

    def test_an_extent_is_declared_before_the_array_using_it(self, bag, emitter):
        # referring to a symbol typed further down is a GNU extension, not standard
        # Fortran: gfortran -std=f2018 rejects it outright
        procedure = b.procedure(
            "p", b.real("vector", Intent.IN, "(n_dims)"), b.integer("n_dims", Intent.IN), b.ierr()
        )

        text = emit(procedure, bag, emitter)

        assert text.index(":: n_dims") < text.index(":: vector")

    def test_documentation_is_inherited_from_the_wrapped_argument(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        assert "!! number of elements in `vector`" in text

    def test_the_wrapper_says_what_it_wraps(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        assert (
            "!> summary: C-wrapper for [[fx_basics(module):fx_normalize(subroutine)]]"
            in emit_module(b.module("fx_basics", normalize()), bag, emitter)
        )


class TestValidation:
    def test_ierr_is_checked_first_and_set_ok(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        body = [line.strip() for line in text.split("\n") if line.strip()]
        start = body.index("M_CHECK_IERR_NON_NULL")

        assert body[start + 1] == "call set_ok(ierr)"

    def test_a_procedure_without_ierr_still_gets_it_set(self, bag, emitter):
        # nothing else would ever set it, so it would return whatever was on the stack
        text = emit(b.procedure("p", b.real("x", Intent.IN)), bag, emitter)

        assert "call set_ok(ierr)" in text

    def test_scalars_are_checked_before_arrays(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        assert text.index("M_CHECK_NON_NULL(n_dims)") < text.index("M_CHECK_ARRAY_NON_NULL(vector")

    def test_an_array_check_is_guarded_by_its_element_count(self, bag, emitter):
        # c_loc may not be given a zero-size target
        text = emit(normalize(), bag, emitter)

        assert "M_CHECK_ARRAY_NON_NULL(vector, n_dims)" in text

    def test_a_matrix_multiplies_its_extents(self, bag, emitter):
        procedure = b.procedure(
            "p", b.integer("n"), b.integer("m"), b.real("mat", Intent.IN, "(n, m)"), b.ierr()
        )

        assert "M_CHECK_ARRAY_NON_NULL(mat, n * m)" in emit(procedure, bag, emitter)

    def test_an_array_with_a_separate_shape_uses_its_product(self, bag, emitter):
        procedure = b.procedure(
            "p", b.real("data", Intent.IN, "(:)"), b.integer("data_shape", Intent.IN, "(:)"), b.ierr()
        )

        text = emit(procedure, bag, emitter)

        assert "M_CHECK_ARRAY_NON_NULL(data, product(data_shape))" in text
        # and the shape itself is checked first, since reading it must be safe
        assert text.index("M_CHECK_ARRAY_NON_NULL(data_shape") < text.index(
            "M_CHECK_ARRAY_NON_NULL(data,"
        )

    def test_optionals_are_not_checked(self, bag, emitter):
        procedure = b.procedure("p", b.real("span", Intent.IN, optional=True), b.ierr())

        assert "M_CHECK_NON_NULL(span)" not in emit(procedure, bag, emitter)


class TestCall:
    def test_arguments_are_passed_by_keyword(self, bag, emitter):
        text = emit(normalize(), bag, emitter)

        assert "call fx_normalize(&" in text
        assert "n_dims = n_dims,&" in text
        assert "vector = vector,&" in text

    def test_a_synthesised_error_is_not_passed_on(self, bag, emitter):
        # the procedure has no ierr to receive it
        text = emit(b.procedure("p", b.real("x", Intent.IN)), bag, emitter)

        assert "ierr = ierr" not in text

    def test_a_declared_error_is_passed_on(self, bag, emitter):
        assert "ierr = ierr" in emit(normalize(), bag, emitter)

    def test_a_function_result_is_assigned(self, bag, emitter):
        procedure = b.procedure(
            "f", b.real("x", Intent.IN), result=b.integer("out", Intent.OUT, is_result=True)
        )

        assert "out = f(&" in emit(procedure, bag, emitter)

    def test_an_assumed_size_array_is_sliced_to_its_real_extent(self, bag, emitter):
        # an assumed-size actual cannot be passed to an assumed-shape dummy, and the
        # count is the product of the shape argument, not its size
        procedure = b.procedure(
            "p", b.real("data", Intent.IN, "(:)"), b.integer("data_shape", Intent.IN, "(:)"), b.ierr()
        )

        assert "data = data(1:product(data_shape)),&" in emit(procedure, bag, emitter)


class TestLogicals:
    def test_a_default_logical_is_converted_by_assignment(self, bag, emitter):
        procedure = b.procedure("p", b.logical("flag", Intent.IN), b.ierr())

        text = emit(procedure, bag, emitter)

        assert "logical(c_bool), intent(in), target :: flag" in text
        assert "logical :: flag_f" in text
        assert "flag_f = flag" in text
        assert "flag = flag_f" in text.replace("flag_f = flag", "")  # only the input way

    def test_a_c_bool_logical_is_passed_straight_through(self, bag, emitter):
        procedure = b.procedure("p", b.logical("flag", Intent.IN, kind="c_bool"), b.ierr())

        text = emit(procedure, bag, emitter)

        assert "logical :: flag_f" not in text
        assert "flag = flag,&" in text

    def test_an_output_logical_is_converted_back(self, bag, emitter):
        procedure = b.procedure("p", b.logical("flag", Intent.OUT), b.ierr())

        text = emit(procedure, bag, emitter)

        assert text.index("flag = flag_f") > text.index("call p(&")


class TestOptionals:
    def test_a_nullable_optional_is_declared_optional_not_target(self, bag, emitter):
        # an OPTIONAL dummy of a bind(C) procedure is absent exactly when C passes a
        # null pointer, and the C prototype stays a plain pointer
        procedure = b.procedure("p", b.real("span", Intent.IN, optional=True), b.ierr())

        text = emit(procedure, bag, emitter)

        assert "real(c_double), intent(in), optional :: span" in text
        assert "target :: span" not in text

    def test_it_is_handed_straight_to_the_callee(self, bag, emitter):
        # an optional associated with an absent optional is itself absent, so presence
        # propagates with no branch in the wrapper
        procedure = b.procedure("p", b.real("span", Intent.IN, optional=True), b.ierr())

        text = emit(procedure, bag, emitter)

        assert "span = span,&" in text
        assert "present(" not in text
        assert "if (c_associated" not in text

    def test_an_optional_with_a_default_is_required_in_c(self, bag, emitter):
        # issue #131: the interfacing languages know the default and pass it, which is
        # what keeps the wrapper flat
        from codegen.ir.directives import Default, Directives

        span = b.real("span", Intent.IN, optional=True)
        span.directives = Directives(default=Default("0.1_real64"))
        text = emit(b.procedure("p", span, b.ierr()), bag, emitter)

        assert "real(c_double), intent(in), target :: span" in text
        assert "optional" not in text
        assert "M_CHECK_NON_NULL(span)" in text


class TestModes:
    def _mode_procedure(self):
        doc = [
            "| Mode | Value |",
            "|------|-------|",
            "| mean it | [[fx_basics(module):MODE_MEAN(variable)]] |",
            "| median it | [[fx_basics(module):MODE_MEDIAN(variable)]] |",
        ]
        return b.procedure("p", b.integer("mode", Intent.IN, doc=doc), b.ierr())

    def test_a_mode_arrives_as_a_string(self, bag, emitter):
        text = emit(self._mode_procedure(), bag, emitter)

        assert "character(len=1, kind=c_char), dimension(6), intent(in), target :: mode" in text

    def test_the_string_is_mapped_to_its_parameter(self, bag, emitter):
        text = emit(self._mode_procedure(), bag, emitter)

        assert 'case ("mean")' in text
        assert "mode_mode_f = MODE_MEAN" in text
        assert 'case ("median")' in text
        assert "mode_mode_f = MODE_MEDIAN" in text

    def test_an_unknown_string_is_an_invalid_input(self, bag, emitter):
        text = emit(self._mode_procedure(), bag, emitter)

        assert "case default" in text
        assert "call set_err(ierr, ERR_INVALID_INPUT)" in text

    def test_the_parameters_are_imported(self, bag, emitter):
        text = emit(self._mode_procedure(), bag, emitter)

        assert "use fx_basics, only: MODE_MEAN, MODE_MEDIAN" in text

    def test_the_integer_is_what_gets_passed(self, bag, emitter):
        assert "mode = mode_mode_f,&" in emit(self._mode_procedure(), bag, emitter)


class TestModule:
    def test_the_module_is_named_and_guarded(self, bag, emitter):
        text = emit_module(b.module("fx_basics", normalize()), bag, emitter)

        assert text.startswith("#ifndef NO_C_INTERFACE")
        assert "#include <src/macros.h>" in text
        assert "module fx_basics_c" in text
        assert text.rstrip().endswith("#endif")

    def test_safeguard_is_used_only_here(self, bag, emitter):
        # the reason the wrappers get their own modules
        assert "use safeguard" in emit_module(b.module("fx_basics", normalize()), bag, emitter)

    def test_implicit_none_forbids_implicit_externals(self, bag, emitter):
        # a bare 'implicit none' constrains only variables, so a call to a procedure
        # that does not exist compiles as an implicit external and fails at link time
        text = emit_module(b.module("fx_basics", normalize()), bag, emitter)

        assert "M_IMPLICIT_NONE" in text
        assert "implicit none\n" not in text

    def test_only_the_wrappers_are_public(self, bag, emitter):
        text = emit_module(b.module("fx_basics", normalize()), bag, emitter)

        assert "private" in text
        assert "public :: fx_normalize_c" in text

    def test_the_module_says_what_it_wraps(self, bag, emitter):
        text = emit_module(b.module("fx_basics", normalize()), bag, emitter)

        assert "!> summary: C-wrappers for [[fx_basics(module)]]" in text


GFORTRAN = shutil.which("gfortran")

#: Flags matching how the project builds, plus -std=f2018 to check conformance. The
#: line-length flag is needed because expanded macros already exceed 132 columns in the
#: existing sources.
COMPILE_FLAGS = ["-cpp", "-I.", "-std=f2018", "-ffree-line-length-none", "-fsyntax-only"]


def compile_fortran(source, module_dir, extra=()):
    return subprocess.run(
        [GFORTRAN, *COMPILE_FLAGS, *extra, f"-J{module_dir}", str(source)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


def must_compile(source, module_dir, extra=()):
    result = compile_fortran(source, module_dir, extra)
    assert result.returncode == 0, f"{source}:\n{result.stderr}"
    return result


@pytest.fixture(scope="session")
def built(tmp_path_factory):
    """Generate the wrappers from the fixtures and compile everything they depend on."""
    from pathlib import Path

    from codegen.abi.c_abi import build_project
    from codegen.config import Paths
    from codegen.frontend.ford_frontend import FordFrontend
    from codegen.ir.roles import analyse_project
    from codegen.ir.validate import validate_project

    out = tmp_path_factory.mktemp("fortran_c")
    bag = DiagnosticBag()
    fixtures = Path("helper/codegen/tests/fixtures/src")
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=fixtures), bag).parse()
    analyse_project(parsed.project, bag)
    validate_project(parsed.project, bag)
    interface = build_project(parsed.project, bag)

    assert bag.errors == (), bag.render()

    emitter = FortranCEmitter()
    for module in interface:
        (out / f"{module.name}.F90").write_text(emitter.module(module))
    for fixture in (REPO_ROOT / fixtures).glob("*.F90"):
        shutil.copy(fixture, out / fixture.name)

    # the wrappers depend on these, so they have to be compiled first
    for source in ("src/tox/tox_errors.F90", "src/tox/tox_conversions.F90", "src/safeguard.F90"):
        must_compile(REPO_ROOT / source, out)
    for name in ("fx_basics", "fx_edges"):
        must_compile(out / f"{name}.F90", out)
    return out


@pytest.mark.skipif(GFORTRAN is None, reason="gfortran is not installed")
class TestItCompiles:
    """Generate from the fixtures and compile. The only proof the output is Fortran."""

    @pytest.mark.parametrize("name", ["fx_basics_c", "fx_edges_c"])
    def test_the_generated_module_compiles(self, built, name):
        must_compile(built / f"{name}.F90", built)

    def test_no_c_interface_compiles_the_whole_thing_out(self, built):
        # the point of the #ifndef: one directive removes the C interface, and with it
        # the need for safeguard
        must_compile(built / "fx_basics_c.F90", built, extra=["-DNO_C_INTERFACE"])

    def test_an_undeclared_procedure_would_not_slip_through(self, built, tmp_path):
        # what implicit none (type, external) is for: a bare implicit none lets a call
        # to a non-existent procedure compile as an implicit external
        source = (built / "fx_basics_c.F90").read_text()
        probe = tmp_path / "probe.F90"
        probe.write_text(source.replace("c_char_1d_as_string(", "c_char_1d_as_string_typo(", 1))

        result = compile_fortran(probe, built)

        assert result.returncode != 0
        assert "not explicitly declared" in result.stderr
