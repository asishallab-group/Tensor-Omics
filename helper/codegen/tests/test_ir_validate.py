import pytest

from codegen.diagnostics import DiagnosticBag
from codegen.ir.entities import Meta
from codegen.ir.roles import analyse
from codegen.ir.types import Intent
from codegen.ir.validate import validate_module, validate_procedure, validate_project

import builders as b


@pytest.fixture
def bag():
    return DiagnosticBag()


def checked(procedure, bag):
    """Analyse then validate, as a real run does."""
    analyse(procedure, bag)
    validate_procedure(procedure, bag)
    return bag


def messages(bag):
    return [d.message for d in bag.errors]


def only_error(bag):
    assert len(bag.errors) == 1, messages(bag)
    return bag.errors[0]


class TestScope:
    def test_only_exported_procedures_are_validated(self, bag):
        # an internal routine may use anything Fortran allows
        module = b.module(
            "m",
            b.procedure(
                "internal",
                b.character("s", Intent.IN, length=":"),
                meta=Meta(summary="s", author="a"),
            ),
        )

        validate_module(module, bag)

        assert bag.errors == ()

    def test_an_exported_procedure_is_validated(self, bag):
        module = b.module("m", b.procedure("exported", b.character("s", Intent.IN, length=":")))

        validate_module(module, bag)

        assert "deferred length" in only_error(bag).message

    def test_validate_project_walks_every_module(self, bag):
        project = b.project(
            b.module("a", b.procedure("p", b.character("s", Intent.IN, length=":"))),
            b.module("b", b.procedure("q", b.character("t", Intent.IN, length=":"))),
        )

        validate_project(project, bag)

        assert len(bag.errors) == 2


class TestCharacters:
    def test_a_deferred_length_is_rejected(self, bag):
        checked(b.procedure("p", b.character("s", Intent.IN, length=":")), bag)

        error = only_error(bag)
        assert "deferred length (len=:)" in error.message
        assert "len=*" in error.note

    @pytest.mark.parametrize("length", ["*", "n", "42"])
    def test_other_lengths_are_accepted(self, bag, length):
        checked(b.procedure("p", b.character("s", Intent.IN, length=length), b.integer("n")), bag)

        assert bag.errors == ()

    def test_a_character_scalar_is_accepted(self, bag):
        checked(b.procedure("p", b.character("s", Intent.IN)), bag)

        assert bag.errors == ()

    def test_a_character_vector_is_accepted(self, bag):
        checked(b.procedure("p", b.character("s", Intent.IN, "(n)"), b.integer("n")), bag)

        assert bag.errors == ()

    def test_a_character_matrix_is_rejected(self, bag):
        # it would need a rank-3 c_char conversion, which tox_conversions has not got
        checked(
            b.procedure("p", b.character("s", Intent.IN, "(n, m)"), b.integer("n"), b.integer("m")),
            bag,
        )

        assert "rank-2 character array" in only_error(bag).message


class TestTypes:
    def test_a_derived_type_is_rejected(self, bag):
        from codegen.ir.entities import Argument
        from codegen.ir.types import BaseType, FortranType

        argument = Argument("p", FortranType(BaseType.DERIVED, derived_name="point"), intent=Intent.IN)
        checked(b.procedure("p", argument), bag)

        assert "derived type" in only_error(bag).message

    @pytest.mark.parametrize("attribute", ["allocatable", "pointer"])
    def test_non_interoperable_attributes_are_rejected(self, bag, attribute):
        checked(
            b.procedure("p", b.real("v", Intent.OUT, "(:)", attributes=(attribute,))), bag
        )

        assert f"is '{attribute}'" in only_error(bag).message

    def test_target_is_fine(self, bag):
        checked(b.procedure("p", b.real("v", Intent.OUT, "(n)", attributes=("target",)), b.integer("n")), bag)

        assert bag.errors == ()


class TestIntent:
    def test_a_missing_intent_is_rejected(self, bag):
        argument = b.real("v", Intent.IN)
        argument.intent = None

        checked(b.procedure("p", argument), bag)

        error = only_error(bag)
        assert "has no intent" in error.message
        assert "constness" in error.note

    def test_a_result_needs_no_declared_intent(self, bag):
        result = b.real("out", Intent.OUT, is_result=True)
        result.intent = None

        checked(b.procedure("f", result=result), bag)

        assert bag.errors == ()


class TestTemporaries:
    def test_an_intent_in_temporary_is_rejected(self, bag):
        checked(b.procedure("p", b.real("tmp_work", Intent.IN, "(n)"), b.integer("n")), bag)

        error = only_error(bag)
        assert "temporary 'tmp_work' is intent(in)" in error.message
        assert "tmp_" in error.note

    @pytest.mark.parametrize("intent", [Intent.OUT, Intent.INOUT])
    def test_out_and_inout_temporaries_are_accepted(self, bag, intent):
        checked(b.procedure("p", b.real("tmp_work", intent, "(n)"), b.integer("n")), bag)

        assert bag.errors == ()


class TestModeArguments:
    def _mode_doc(self):
        return [
            "| Mode | Value |",
            "|------|-------|",
            "| x | [[m(module):MODE_X(variable)]] |",
        ]

    def test_an_integer_scalar_mode_is_accepted(self, bag):
        checked(b.procedure("p", b.integer("mode", Intent.IN, doc=self._mode_doc())), bag)

        assert bag.errors == ()

    def test_a_non_integer_mode_is_rejected(self, bag):
        # the Fortran mode is an integer; only the C wrapper receives a string
        checked(b.procedure("p", b.character("mode", Intent.IN, doc=self._mode_doc())), bag)

        assert "is not a scalar integer" in only_error(bag).message

    def test_an_array_mode_is_rejected(self, bag):
        checked(
            b.procedure("p", b.integer("mode", Intent.IN, "(n)", doc=self._mode_doc()), b.integer("n")),
            bag,
        )

        assert "is not a scalar integer" in only_error(bag).message

    def _split_mode_doc(self):
        return [
            "| Mode | Value | Procedure |",
            "|------|-------|-----------|",
            "| x | [[m(module):MODE_X(variable)]] | p_x |",
        ]

    def _scoped_optional(self, mode_doc):
        from codegen.ir.directives import Default, Directives, RequiredIfMode

        return b.procedure(
            "p",
            b.integer("mode", Intent.IN, doc=mode_doc),
            b.real(
                "tuning",
                Intent.IN,
                optional=True,
                doc="a tuning knob",
                directives=Directives(
                    default=Default("1.0_real64"),
                    required_if_mode=RequiredIfMode("mode", "m", "MODE_X"),
                ),
            ),
        )

    def test_a_default_with_a_required_if_mode_is_rejected_for_a_runtime_mode(self, bag):
        # the argument is always passed on, so "required in that mode" says nothing
        checked(self._scoped_optional(self._mode_doc()), bag)

        assert "both a default and a mode it is required in" in only_error(bag).message

    def test_a_default_with_a_required_if_mode_is_accepted_when_the_mode_splits(self, bag):
        # there the directive scopes the argument to its mode, and the default applies within it
        checked(self._scoped_optional(self._split_mode_doc()), bag)

        assert bag.errors == ()


class TestShapeArguments:
    def _shape_procedure(self, shape_argument, data=None):
        return b.procedure("p", data or b.real("data", Intent.IN, "(:)"), shape_argument)

    def test_a_well_formed_shape_argument_is_accepted(self, bag):
        checked(self._shape_procedure(b.integer("data_shape", Intent.IN, "(:)")), bag)

        assert bag.errors == ()

    def test_a_shape_argument_must_be_intent_in(self, bag):
        checked(self._shape_procedure(b.integer("data_shape", Intent.OUT, "(:)")), bag)

        assert "is not intent(in)" in only_error(bag).message

    def test_a_shape_argument_must_be_a_rank_1_integer(self, bag):
        checked(self._shape_procedure(b.integer("data_shape", Intent.IN)), bag)

        assert "is not a rank-1 integer array" in only_error(bag).message

    def test_a_real_shape_argument_is_rejected(self, bag):
        checked(self._shape_procedure(b.real("data_shape", Intent.IN, "(:)")), bag)

        assert "is not a rank-1 integer array" in only_error(bag).message

    def test_a_shape_argument_must_not_be_optional(self, bag):
        # the c_loc ordering depends on it being there
        checked(self._shape_procedure(b.integer("data_shape", Intent.IN, "(:)", optional=True)), bag)

        error = only_error(bag)
        assert "shape argument 'data_shape' is optional" in error.message
        assert "c_loc" in error.note

    def test_the_described_argument_must_be_flat(self, bag):
        checked(
            self._shape_procedure(
                b.integer("data_shape", Intent.IN, "(:)"),
                data=b.real("data", Intent.IN, "(:, :)"),
            ),
            bag,
        )

        assert "has a shape argument but is rank 2" in only_error(bag).message


class TestExtentArguments:
    def test_an_extent_must_not_be_optional(self, bag):
        checked(
            b.procedure("p", b.integer("n", Intent.IN, optional=True), b.real("v", Intent.IN, "(n)")),
            bag,
        )

        error = only_error(bag)
        assert "extent argument 'n' is optional, but 'v' needs it" in error.message
        assert "c_loc" in error.note

    def test_the_error_names_every_dependent_array(self, bag):
        checked(
            b.procedure(
                "p",
                b.integer("n", Intent.IN, optional=True),
                b.real("a", Intent.IN, "(n)"),
                b.real("c", Intent.IN, "(n)"),
            ),
            bag,
        )

        assert "'a', 'c'" in only_error(bag).message

    def test_a_non_extent_scalar_may_be_optional(self, bag):
        checked(b.procedure("p", b.integer("seed", Intent.IN, optional=True)), bag)

        assert bag.errors == ()


class TestErrorArgument:
    def test_a_well_formed_ierr_is_accepted(self, bag):
        checked(b.procedure("p", b.ierr()), bag)

        assert bag.errors == ()

    def test_an_intent_in_ierr_is_rejected(self, bag):
        checked(b.procedure("p", b.integer("ierr", Intent.IN, doc="Error code")), bag)

        assert "no error can be reported through it" in only_error(bag).message

    def test_a_non_integer_ierr_is_rejected(self, bag):
        checked(b.procedure("p", b.real("ierr", Intent.OUT, doc="Error code")), bag)

        assert "is not a scalar integer" in only_error(bag).message

    def test_an_array_ierr_is_rejected(self, bag):
        checked(b.procedure("p", b.integer("ierr", Intent.OUT, "(2)", doc="Error code")), bag)

        assert "is not a scalar integer" in only_error(bag).message


class TestWarnings:
    def test_a_procedure_without_a_summary_warns(self, bag):
        checked(b.procedure("p", meta=Meta(author="a", category="C-binding")), bag)

        assert bag.errors == ()
        assert "no summary meta tag" in bag.warnings[0].message

    def test_a_procedure_without_an_author_warns(self, bag):
        checked(b.procedure("p", meta=Meta(summary="s", category="C-binding")), bag)

        assert "no author meta tag" in bag.warnings[0].message

    def test_an_undocumented_argument_warns(self, bag):
        checked(b.procedure("p", b.integer("n", Intent.IN)), bag)

        assert bag.errors == ()
        assert "argument 'n' has no documentation" in bag.warnings[0].message

    def test_an_undocumented_result_does_not_warn(self, bag):
        checked(b.procedure("f", result=b.real("out", Intent.OUT, is_result=True)), bag)

        assert bag.warnings == () or all("result" not in w.message for w in bag.warnings)

    def test_an_undocumented_module_warns(self, bag):
        validate_module(b.module("m"), bag)

        assert "module has no documentation" in bag.warnings[0].message

    def test_warnings_never_stop_generation(self, bag):
        checked(b.procedure("p", b.integer("n", Intent.IN), meta=Meta(category="C-binding")), bag)

        bag.raise_if_errors()  # must not raise


class TestEverythingTogether:
    def test_a_clean_procedure_produces_nothing(self, bag):
        procedure = b.procedure(
            "normalize_unit_length",
            b.integer("n_dims", Intent.IN, doc="number of elements in `vector`"),
            b.real("vector", Intent.INOUT, "(n_dims)", doc="Vector to normalise"),
            b.ierr(),
        )

        checked(procedure, bag)

        assert bag.errors == ()
        assert bag.warnings == ()

    def test_every_problem_is_reported_in_one_run(self, bag):
        procedure = b.procedure(
            "p",
            b.character("s", Intent.IN, length=":", doc="a string"),
            b.real("tmp_work", Intent.IN, "(n)", doc="work"),
            b.integer("n", Intent.IN, optional=True, doc="extent"),
        )

        checked(procedure, bag)

        assert len(bag.errors) == 3
        assert any("deferred length" in m for m in messages(bag))
        assert any("is intent(in)" in m for m in messages(bag))
        assert any("is optional" in m for m in messages(bag))
