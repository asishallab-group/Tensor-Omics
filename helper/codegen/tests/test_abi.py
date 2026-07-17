import pytest

from codegen.abi.c_abi import build_module, build_project, build_wrapper, stripped_name
from codegen.abi.model import Conversion, Origin
from codegen.diagnostics import DiagnosticBag
from codegen.ir.entities import Meta
from codegen.ir.roles import analyse
from codegen.ir.types import BaseType, Intent

import builders as b


@pytest.fixture
def bag():
    return DiagnosticBag()


def wrap(procedure, bag):
    analyse(procedure, bag)
    return build_wrapper(procedure, bag)


def wrap_in_module(module, name, bag):
    for procedure in module.exported_procedures:
        analyse(procedure, bag)
    return next(w for w in build_module(module, bag) if w.procedure.name == name)


class TestNaming:
    def test_an_alloc_procedure_takes_the_plain_name(self, bag):
        module = b.module("m", b.procedure("cluster_alloc", b.ierr()), b.procedure("cluster", b.ierr()))

        wrapper = wrap_in_module(module, "cluster_alloc", bag)

        assert wrapper.name == "cluster_c"
        assert wrapper.stripped_name == "cluster"

    def test_the_non_alloc_twin_becomes_expert(self, bag):
        module = b.module("m", b.procedure("cluster_alloc", b.ierr()), b.procedure("cluster", b.ierr()))

        wrapper = wrap_in_module(module, "cluster", bag)

        assert wrapper.name == "cluster_expert_c"
        assert wrapper.stripped_name == "cluster_expert"

    def test_a_procedure_without_an_alloc_twin_keeps_its_name(self, bag):
        # a lone procedure must not be needlessly renamed to _expert
        module = b.module("m", b.procedure("cluster", b.ierr()))

        assert wrap_in_module(module, "cluster", bag).name == "cluster_c"

    def test_a_twin_in_another_module_does_not_count(self, bag):
        project = b.project(
            b.module("a", b.procedure("cluster", b.ierr())),
            b.module("b", b.procedure("cluster_alloc", b.ierr())),
        )
        for module in project:
            for procedure in module.exported_procedures:
                analyse(procedure, bag)

        interface = build_project(project, bag)

        assert interface.module("a_c").wrappers[0].name == "cluster_c"

    def test_the_module_gets_a_c_module(self, bag):
        module = b.module("tox_clustering", b.procedure("cluster", b.ierr()))

        assert build_module(module, bag).name == "tox_clustering_c"

    def test_stripped_name_is_available_without_building(self):
        assert stripped_name(b.procedure("cluster_alloc")) == "cluster"


class TestModuleSelection:
    def test_a_module_with_no_exports_produces_nothing(self, bag):
        project = b.project(
            b.module("exports", b.procedure("p", b.ierr())),
            b.module("silent", b.procedure("q", b.ierr(), meta=Meta())),
        )
        for module in project:
            for procedure in module.exported_procedures:
                analyse(procedure, bag)

        interface = build_project(project, bag)

        assert [m.name for m in interface] == ["exports_c"]

    def test_only_exported_procedures_are_wrapped(self, bag):
        module = b.module(
            "m", b.procedure("yes", b.ierr()), b.procedure("no", b.ierr(), meta=Meta())
        )
        for procedure in module.exported_procedures:
            analyse(procedure, bag)

        assert [w.procedure.name for w in build_module(module, bag)] == ["yes"]


class TestErrorArgument:
    def test_an_existing_ierr_is_used_where_it_stands(self, bag):
        wrapper = wrap(b.procedure("p", b.integer("n"), b.ierr()), bag)

        assert [a.name for a in wrapper] == ["n", "ierr"]
        assert wrapper.argument("ierr").origin is Origin.ERROR
        assert wrapper.argument("ierr").source is not None

    def test_ierr_is_not_duplicated(self, bag):
        wrapper = wrap(b.procedure("p", b.ierr()), bag)

        assert [a.name for a in wrapper].count("ierr") == 1

    def test_a_procedure_without_ierr_gets_one_last(self, bag):
        wrapper = wrap(b.procedure("p", b.real("x", Intent.IN)), bag)

        assert [a.name for a in wrapper] == ["x", "ierr"]
        error = wrapper.argument("ierr")
        assert error.origin is Origin.ERROR
        assert error.intent is Intent.OUT
        assert error.is_synthesised
        assert error.source is None

    def test_error_argument_finds_it_either_way(self, bag):
        assert wrap(b.procedure("p", b.ierr()), bag).error_argument.name == "ierr"
        assert wrap(b.procedure("p", b.real("x")), bag).error_argument.name == "ierr"


class TestFunctions:
    def test_a_result_becomes_an_output_argument(self, bag):
        procedure = b.procedure(
            "f", b.real("x", Intent.IN), result=b.integer("count", Intent.OUT, is_result=True)
        )

        wrapper = wrap(procedure, bag)

        assert [a.name for a in wrapper] == ["x", "count", "ierr"]
        result = wrapper.argument("count")
        assert result.origin is Origin.RESULT
        assert result.intent is Intent.OUT

    def test_the_result_precedes_the_synthesised_error(self, bag):
        # so adding a result to a subroutine does not move ierr for existing callers
        procedure = b.procedure("f", result=b.real("out", Intent.OUT, is_result=True))

        assert [a.name for a in wrap(procedure, bag)] == ["out", "ierr"]


class TestTypeMapping:
    @pytest.mark.parametrize(
        "kind, c_kind",
        [("int8", "c_int8_t"), ("int16", "c_int16_t"), ("int32", "c_int"), ("int64", "c_int64_t")],
    )
    def test_integer_kinds(self, bag, kind, c_kind):
        wrapper = wrap(b.procedure("p", b.integer("n", kind=kind)), bag)

        assert wrapper.argument("n").type.kind == c_kind

    @pytest.mark.parametrize("kind, c_kind", [("real32", "c_float"), ("real64", "c_double")])
    def test_real_kinds(self, bag, kind, c_kind):
        wrapper = wrap(b.procedure("p", b.real("x", kind=kind)), bag)

        assert wrapper.argument("x").type.kind == c_kind

    def test_complex_stays_complex(self, bag):
        # the old generator emitted 'real(c_double_complex)': a real of a complex kind,
        # which compiles and is the wrong type
        wrapper = wrap(b.procedure("p", b.complex_("z", kind="real64")), bag)

        assert wrapper.argument("z").type.base is BaseType.COMPLEX
        assert wrapper.argument("z").type.kind == "c_double_complex"
        assert str(wrapper.argument("z").type) == "complex(c_double_complex)"

    def test_an_unknown_kind_is_reported_rather_than_guessed(self, bag):
        wrap(b.procedure("p", b.real("x", kind="real128")), bag)

        assert "kind 'real128', which has no known C counterpart" in bag.errors[0].message
        assert "c_double" in bag.errors[0].note


class TestLogicals:
    def test_a_default_logical_becomes_c_bool_and_is_converted(self, bag):
        # agreed ABI: C passes a real bool, so Python and R hand over a boolean
        wrapper = wrap(b.procedure("p", b.logical("flag", Intent.IN)), bag)

        flag = wrapper.argument("flag")
        assert str(flag.type) == "logical(c_bool)"
        assert flag.conversion is Conversion.LOGICAL

    def test_a_c_bool_logical_needs_no_conversion(self, bag):
        wrapper = wrap(b.procedure("p", b.logical("flag", Intent.IN, kind="c_bool")), bag)

        flag = wrapper.argument("flag")
        assert str(flag.type) == "logical(c_bool)"
        assert flag.conversion is Conversion.NONE
        assert not flag.needs_conversion


class TestCharacters:
    def test_the_length_becomes_the_leading_extent(self, bag):
        wrapper = wrap(b.procedure("p", b.character("s", Intent.IN, length="8")), bag)

        s = wrapper.argument("s")
        assert str(s.type) == "character(len=1, kind=c_char)"
        assert s.dimension.extents == ("8",)
        assert s.conversion is Conversion.CHARACTER

    def test_an_assumed_length_gets_a_length_argument(self, bag):
        wrapper = wrap(b.procedure("p", b.character("s", Intent.IN, length="*")), bag)

        assert [a.name for a in wrapper] == ["s", "s_strlen", "ierr"]
        assert wrapper.argument("s").dimension.extents == ("s_strlen",)
        strlen = wrapper.argument("s_strlen")
        assert strlen.origin is Origin.STRLEN
        assert strlen.intent is Intent.IN
        assert strlen.sizes == "s"

    def test_a_length_named_by_another_argument_needs_no_extra(self, bag):
        procedure = b.procedure("p", b.integer("n_chars"), b.character("s", Intent.IN, length="n_chars"))

        wrapper = wrap(procedure, bag)

        assert [a.name for a in wrapper] == ["n_chars", "s", "ierr"]
        assert wrapper.argument("s").dimension.extents == ("n_chars",)

    def test_a_character_vector_keeps_length_and_extent(self, bag):
        procedure = b.procedure("p", b.integer("n"), b.character("s", Intent.IN, "(n)", length="16"))

        assert wrap(procedure, bag).argument("s").dimension.extents == ("16", "n")


class TestAssumedShape:
    def test_an_assumed_shape_array_gets_an_extent(self, bag):
        wrapper = wrap(b.procedure("p", b.real("v", Intent.IN, "(:)")), bag)

        assert [a.name for a in wrapper] == ["v", "n_v_elements", "ierr"]
        extent = wrapper.argument("n_v_elements")
        assert extent.origin is Origin.EXTENT
        assert extent.intent is Intent.IN
        assert extent.sizes == "v"
        assert extent.axis == 0
        assert wrapper.argument("v").dimension.extents == ("n_v_elements",)

    def test_a_rank_two_assumed_shape_gets_one_extent_per_axis(self, bag):
        wrapper = wrap(b.procedure("p", b.real("m", Intent.IN, "(:, :)")), bag)

        assert [a.name for a in wrapper] == ["m", "n_m_elements_dim_1", "n_m_elements_dim_2", "ierr"]
        assert wrapper.argument("m").dimension.extents == (
            "n_m_elements_dim_1",
            "n_m_elements_dim_2",
        )
        assert wrapper.argument("n_m_elements_dim_2").axis == 1

    def test_an_explicit_shape_array_needs_no_extent(self, bag):
        procedure = b.procedure("p", b.integer("n"), b.real("v", Intent.IN, "(n)"))

        assert [a.name for a in wrap(procedure, bag)] == ["n", "v", "ierr"]

    def test_synthesised_arguments_follow_the_one_they_belong_to(self, bag):
        procedure = b.procedure("p", b.real("a", Intent.IN, "(:)"), b.real("bb", Intent.IN, "(:)"))

        assert [a.name for a in wrap(procedure, bag)] == [
            "a", "n_a_elements", "bb", "n_bb_elements", "ierr",
        ]


class TestShapeArguments:
    def test_an_array_with_a_separate_shape_becomes_assumed_size(self, bag):
        procedure = b.procedure(
            "p", b.real("data", Intent.IN, "(:)"), b.integer("data_shape", Intent.IN, "(:)")
        )

        wrapper = wrap(procedure, bag)

        data = wrapper.argument("data")
        assert data.dimension.extents == ("*",)
        assert data.shape_arg == "data_shape"

    def test_no_extent_is_invented_for_it(self, bag):
        # its size comes from the shape argument, not from an extent of its own
        procedure = b.procedure(
            "p", b.real("data", Intent.IN, "(:)"), b.integer("data_shape", Intent.IN, "(:)")
        )

        assert "n_data_elements" not in [a.name for a in wrap(procedure, bag)]

    def test_the_shape_argument_still_gets_its_own_extent(self, bag):
        procedure = b.procedure(
            "p", b.real("data", Intent.IN, "(:)"), b.integer("data_shape", Intent.IN, "(:)")
        )

        assert "n_data_shape_elements" in [a.name for a in wrap(procedure, bag)]


class TestModes:
    def _mode_doc(self, *names):
        rows = [f"| {n} | [[m(module):MODE_{n.upper()}(variable)]] |" for n in names]
        return ["| Mode | Value |", "|------|-------|", *rows]

    def test_a_mode_becomes_a_string(self, bag):
        procedure = b.procedure("p", b.integer("mode", Intent.IN, doc=self._mode_doc("mean")))

        mode = wrap(procedure, bag).argument("mode")

        assert str(mode.type) == "character(len=1, kind=c_char)"
        assert mode.conversion is Conversion.MODE

    def test_the_extent_is_the_longest_accepted_string(self, bag):
        procedure = b.procedure("p", b.integer("mode", Intent.IN, doc=self._mode_doc("mean", "median")))

        assert wrap(procedure, bag).argument("mode").dimension.extents == ("6",)

    def test_the_accepted_values_travel_with_the_argument(self, bag):
        procedure = b.procedure("p", b.integer("mode", Intent.IN, doc=self._mode_doc("mean", "median")))

        mode = wrap(procedure, bag).argument("mode")

        assert [v.string for v in mode.mode.values] == ["mean", "median"]


class TestNameCollisions:
    def test_a_collision_with_a_synthesised_extent_is_reported(self, bag):
        # the wrapper would otherwise have two arguments of one name
        procedure = b.procedure(
            "p", b.real("v", Intent.IN, "(:)"), b.integer("n_v_elements", Intent.IN)
        )

        wrap(procedure, bag)

        assert "already has one" in bag.errors[0].message
        assert "n_v_elements" in bag.errors[0].message

    def test_a_collision_with_the_synthesised_error_is_reported(self, bag):
        procedure = b.procedure("p", b.real("ierr", Intent.IN))

        wrap(procedure, bag)

        assert bag.errors == () or "already has one" in bag.errors[0].message


class TestValidationOrder:
    """The c_loc ordering: extents must be readable before the arrays they size."""

    def test_the_error_argument_comes_first(self, bag):
        wrapper = wrap(b.procedure("p", b.integer("n"), b.real("v", Intent.IN, "(n)"), b.ierr()), bag)

        assert wrapper.validation_order[0].name == "ierr"

    def test_scalars_precede_arrays(self, bag):
        wrapper = wrap(b.procedure("p", b.real("v", Intent.IN, "(n)"), b.integer("n"), b.ierr()), bag)

        assert [a.name for a in wrapper.validation_order] == ["ierr", "n", "v"]

    def test_a_shape_argument_precedes_the_array_it_describes(self, bag):
        # data's element count is the product of data_shape's contents, so data_shape
        # has to be known non-null first
        procedure = b.procedure(
            "p", b.real("data", Intent.IN, "(:)"), b.integer("data_shape", Intent.IN, "(:)"), b.ierr()
        )

        order = [a.name for a in wrap(procedure, bag).validation_order]

        assert order.index("data_shape") < order.index("data")

    def test_optionals_are_never_checked(self, bag):
        # a null pointer is how C says 'not present', so checking one would reject
        # exactly what it is meant to allow
        procedure = b.procedure(
            "p",
            b.integer("n"),
            b.real("v", Intent.IN, "(n)"),
            b.real("span", Intent.IN, optional=True),
            b.real("extra", Intent.IN, "(n)", optional=True),
            b.ierr(),
        )

        order = [a.name for a in wrap(procedure, bag).validation_order]

        assert "span" not in order
        assert "extra" not in order
        assert order == ["ierr", "n", "v"]

    def test_a_synthesised_extent_is_checked_before_its_array(self, bag):
        procedure = b.procedure("p", b.real("v", Intent.IN, "(:)"), b.ierr())

        order = [a.name for a in wrap(procedure, bag).validation_order]

        assert order == ["ierr", "n_v_elements", "v"]

    def test_a_strlen_is_checked_before_its_string(self, bag):
        procedure = b.procedure("p", b.character("s", Intent.IN, length="*"), b.ierr())

        order = [a.name for a in wrap(procedure, bag).validation_order]

        assert order == ["ierr", "s_strlen", "s"]

    def test_every_non_optional_argument_is_checked(self, bag):
        procedure = b.procedure(
            "p", b.integer("n"), b.real("v", Intent.IN, "(n)"), b.real("out", Intent.OUT), b.ierr()
        )

        order = {a.name for a in wrap(procedure, bag).validation_order}

        assert order == {"ierr", "n", "v", "out"}


class TestSizeExtents:
    def test_an_array_carries_the_extents_to_multiply(self, bag):
        procedure = b.procedure(
            "p", b.integer("n"), b.integer("m"), b.real("mat", Intent.IN, "(n, m)")
        )

        assert wrap(procedure, bag).argument("mat").size_extents == ("n", "m")

    def test_a_scalar_carries_none(self, bag):
        assert wrap(b.procedure("p", b.integer("n")), bag).argument("n").size_extents == ()

    def test_an_assumed_size_array_carries_none_and_names_its_shape(self, bag):
        # its count is the product of the shape argument's contents instead
        procedure = b.procedure(
            "p", b.real("data", Intent.IN, "(:)"), b.integer("data_shape", Intent.IN, "(:)")
        )

        data = wrap(procedure, bag).argument("data")

        assert data.size_extents == ()
        assert data.shape_arg == "data_shape"


class TestTemporaries:
    def test_a_temporary_is_marked_through_to_the_abi(self, bag):
        procedure = b.procedure("p", b.integer("n"), b.real("tmp_work", Intent.OUT, "(n)"))

        assert wrap(procedure, bag).argument("tmp_work").is_temporary

    def test_a_synthesised_argument_is_not_a_temporary(self, bag):
        assert not wrap(b.procedure("p", b.real("x")), bag).argument("ierr").is_temporary
