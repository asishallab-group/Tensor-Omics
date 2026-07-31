import math

import pytest

from codegen.diagnostics import Context, SourceLocation
from codegen.ir.constants import ConstantEvaluator
from codegen.ir.entities import (
    Argument,
    Meta,
    Module,
    Parameter,
    Procedure,
    ProcedureKind,
    Project,
)
from codegen.ir.types import BaseType, Dimension, FortranType, Intent

import builders as b


class TestConstructibleWithoutFord:
    """The property the whole rework rests on."""

    def test_an_argument_stands_alone(self):
        argument = b.real("vector", Intent.INOUT, "(n_dims)")

        assert argument.name == "vector"
        assert argument.rank == 1
        assert argument.intent is Intent.INOUT

    def test_a_procedure_stands_alone(self):
        procedure = b.procedure("normalize", b.real("vector", Intent.INOUT, "(n)"), b.ierr())

        assert procedure.name == "normalize"
        assert len(procedure.arguments) == 2

    def test_a_project_stands_alone(self):
        project = b.project(b.module("tox_x", b.procedure("p", b.ierr())))

        assert project.procedure("tox_x", "p").name == "p"


class TestArgument:
    def test_scalar_and_array(self):
        assert b.integer("n").is_scalar
        assert not b.integer("n").is_array
        assert b.real("v", dimension="(n)").is_array
        assert b.real("m", dimension="(n, k)").rank == 2

    def test_intent_defaults_to_inout(self):
        # matching Fortran, where an omitted intent means the argument may be both
        assert Argument("x", FortranType(BaseType.INTEGER, kind="int32")).intent is Intent.INOUT

    def test_attributes_are_normalised_to_lower_case(self):
        argument = b.real("v", dimension="(n)", attributes=("TARGET", "Allocatable"))

        assert argument.attributes == ("target", "allocatable")
        assert argument.has_attribute("target")
        assert not argument.has_attribute("pointer")

    def test_optional_defaults_to_false(self):
        assert not b.integer("n").optional
        assert b.integer("n", optional=True).optional

    def test_roles_are_unset_until_analysed(self):
        assert b.integer("n").roles is None

    def test_with_name_copies_everything_else(self):
        original = b.real("v", Intent.INOUT, "(n)", doc="a vector")

        renamed = original.with_name("w")

        assert renamed.name == "w"
        assert renamed.type == original.type
        assert renamed.dimension == original.dimension
        assert renamed.intent is original.intent
        assert renamed.doc is original.doc
        assert original.name == "v"

    def test_a_result_is_marked(self):
        assert not b.integer("n").is_result
        assert b.integer("n", is_result=True).is_result


class TestParenting:
    def test_a_procedure_adopts_its_arguments(self):
        argument = b.real("v", dimension="(n)")
        procedure = b.procedure("p", argument)

        assert argument.parent is procedure
        assert argument.procedure is procedure

    def test_a_procedure_adopts_its_result(self):
        result = b.real("out", Intent.OUT, is_result=True)
        procedure = b.procedure("f", result=result)

        assert result.parent is procedure

    def test_a_module_adopts_its_procedures_and_parameters(self):
        procedure = b.procedure("p")
        parameter = b.parameter("MODE_X", "1")
        module = b.module("m", procedure, parameters=(parameter,))

        assert procedure.parent is module
        assert parameter.parent is module

    def test_a_project_adopts_its_modules(self):
        module = b.module("m")
        project = b.project(module)

        assert module.parent is project

    def test_the_chain_reaches_the_project(self):
        argument = b.real("v")
        b.project(b.module("tox_x", b.procedure("p", argument)))

        assert argument.parent.parent.parent.entity_kind == "project"

    def test_diagnostics_can_render_the_chain(self):
        # the reason parents exist at all
        argument = b.real("v")
        b.project(b.module("tox_x", b.procedure("normalize", argument)))

        assert str(Context.of(argument)).startswith(
            "argument 'v' in procedure 'normalize' in module 'tox_x'"
        )

    def test_qualified_name(self):
        argument = b.real("v")
        b.project(b.module("tox_x", b.procedure("normalize", argument)))

        assert argument.qualified_name == "tox_x.normalize.v"

    def test_qualified_name_of_an_orphan_is_just_its_name(self):
        assert b.real("v").qualified_name == "v"


class TestProcedure:
    def test_a_subroutine_has_no_result(self):
        procedure = b.procedure("p", b.ierr())

        assert not procedure.is_function
        assert procedure.kind == ProcedureKind.SUBROUTINE

    def test_a_function_has_a_result(self):
        procedure = b.procedure("f", result=b.real("out", Intent.OUT, is_result=True))

        assert procedure.is_function
        assert procedure.kind == ProcedureKind.FUNCTION

    def test_argument_lookup_is_case_insensitive(self):
        procedure = b.procedure("p", b.integer("n_dims"))

        assert procedure.argument("N_DIMS").name == "n_dims"
        assert procedure.argument("absent") is None

    def test_has_error_argument(self):
        assert b.procedure("p", b.ierr()).has_error_argument
        assert not b.procedure("p", b.integer("n")).has_error_argument

    def test_iteration_walks_arguments(self):
        procedure = b.procedure("p", b.integer("a"), b.integer("b"))

        assert [a.name for a in procedure] == ["a", "b"]


class TestExportSelection:
    def test_the_c_binding_category_marks_a_procedure(self):
        assert b.procedure("p", meta=Meta(category="C-binding")).is_exported

    def test_the_category_is_matched_case_insensitively_and_trimmed(self):
        assert b.procedure("p", meta=Meta(category="  c-binding  ")).is_exported

    def test_another_category_does_not_export(self):
        assert not b.procedure("p", meta=Meta(category="internal")).is_exported

    def test_no_category_does_not_export(self):
        assert not b.procedure("p", meta=Meta()).is_exported

    def test_a_module_lists_only_exported_procedures(self):
        module = b.module(
            "m",
            b.procedure("exported", meta=Meta(category="C-binding")),
            b.procedure("internal", meta=Meta()),
        )

        assert [p.name for p in module.exported_procedures] == ["exported"]
        assert module.has_exports

    def test_a_module_without_exports_reports_none(self):
        assert not b.module("m", b.procedure("internal", meta=Meta())).has_exports


class TestAllocSibling:
    def test_the_expert_half_finds_its_alloc_sibling(self):
        # the old generator searched the procedure's children instead of the module's,
        # so it never found this and the _expert_c naming rule never fired
        module = b.module("m", b.procedure("cluster"), b.procedure("cluster_alloc"))

        assert module.procedure("cluster").has_alloc_sibling
        assert module.procedure("cluster").alloc_sibling.name == "cluster_alloc"

    def test_an_alloc_procedure_is_not_its_own_sibling(self):
        module = b.module("m", b.procedure("cluster"), b.procedure("cluster_alloc"))

        assert module.procedure("cluster_alloc").is_alloc_variant
        assert not module.procedure("cluster_alloc").has_alloc_sibling

    def test_a_lone_procedure_has_no_sibling(self):
        module = b.module("m", b.procedure("cluster"))

        assert not module.procedure("cluster").has_alloc_sibling

    def test_a_sibling_in_another_module_does_not_count(self):
        # the naming rule is explicitly scoped to the same module
        project = b.project(
            b.module("a", b.procedure("cluster")),
            b.module("b", b.procedure("cluster_alloc")),
        )

        assert not project.procedure("a", "cluster").has_alloc_sibling

    def test_an_unparented_procedure_has_no_sibling(self):
        assert not b.procedure("cluster").has_alloc_sibling

    def test_sibling_lookup_is_case_insensitive(self):
        module = b.module("m", b.procedure("Cluster"), b.procedure("CLUSTER_ALLOC"))

        assert module.procedure("cluster").has_alloc_sibling


class TestModule:
    def test_lookup_by_name(self):
        module = b.module("m", b.procedure("p"), parameters=(b.parameter("MODE_X", "1"),))

        assert module.procedure("p").name == "p"
        assert module.parameter("mode_x").name == "MODE_X"
        assert module.procedure("absent") is None
        assert module.parameter("absent") is None


class TestProject:
    def test_modules_are_sorted_for_deterministic_output(self):
        project = b.project(b.module("zeta"), b.module("alpha"), b.module("Mid"))

        assert [m.name for m in project] == ["alpha", "Mid", "zeta"]

    def test_lookup_by_name(self):
        project = b.project(b.module("tox_x", b.procedure("p")))

        assert project.module("TOX_X").name == "tox_x"
        assert project.procedure("tox_x", "p").name == "p"
        assert project.procedure("absent", "p") is None
        assert project.module("absent") is None

    def test_modules_with_exports(self):
        project = b.project(
            b.module("a", b.procedure("p", meta=Meta(category="C-binding"))),
            b.module("b", b.procedure("q", meta=Meta())),
        )

        assert [m.name for m in project.modules_with_exports] == ["a"]

    def test_parameters_are_collected_across_modules(self):
        project = b.project(
            b.module("a", parameters=(b.parameter("X", "1"),)),
            b.module("b", parameters=(b.parameter("Y", "2"),)),
        )

        assert {p.name for p in project.parameters} == {"X", "Y"}

    def test_len_counts_modules(self):
        assert len(b.project(b.module("a"), b.module("b"))) == 2


class TestConstantValues:
    def test_literal_parameters_resolve(self):
        project = b.project(
            b.module(
                "tox_errors",
                parameters=(b.parameter("ERR_OK", "0"), b.parameter("ERR_UNKNOWN", "9999")),
            )
        )

        assert project.constant_values() == {"ERR_OK": 0, "ERR_UNKNOWN": 9999}

    def test_a_parameter_defined_from_another_resolves(self):
        project = b.project(
            b.module(
                "m",
                parameters=(b.parameter("BASE", "100"), b.parameter("DERIVED", "BASE + 1")),
            )
        )

        assert project.constant_values()["DERIVED"] == 101

    def test_resolution_does_not_depend_on_declaration_order(self):
        # the dependent parameter is declared before the one it needs
        project = b.project(
            b.module(
                "m",
                parameters=(b.parameter("DERIVED", "BASE + 1"), b.parameter("BASE", "100")),
            )
        )

        assert project.constant_values()["DERIVED"] == 101

    def test_a_parameter_defined_across_modules_resolves(self):
        project = b.project(
            b.module("b", parameters=(b.parameter("DERIVED", "BASE * 2"),)),
            b.module("a", parameters=(b.parameter("BASE", "21"),)),
        )

        assert project.constant_values()["DERIVED"] == 42

    def test_non_constant_parameters_are_skipped_not_reported(self):
        project = b.project(
            b.module(
                "m",
                parameters=(
                    b.parameter("GOOD", "1"),
                    b.parameter("BAD", "some_function(x)"),
                ),
            )
        )

        values = project.constant_values()

        assert values == {"GOOD": 1}

    def test_a_cycle_terminates_rather_than_hanging(self):
        project = b.project(
            b.module("m", parameters=(b.parameter("A", "B + 1"), b.parameter("B", "A + 1")))
        )

        assert project.constant_values() == {}

    def test_parameters_without_an_expression_are_skipped(self):
        project = b.project(b.module("m", parameters=(b.parameter("X", ""),)))

        assert project.constant_values() == {}

    def test_the_values_feed_the_default_evaluator(self):
        # the point of collecting them: DM_DEFAULT(PI) has to resolve
        project = b.project(
            b.module("config", parameters=(b.parameter("PI", "3.141592653589793"),))
        )

        evaluate = ConstantEvaluator(project.constant_values()).evaluate

        assert evaluate("PI") == pytest.approx(math.pi)
        assert evaluate("2.0_real64 * PI") == pytest.approx(2 * math.pi)
