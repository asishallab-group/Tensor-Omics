import pytest

from codegen.diagnostics import DiagnosticBag
from codegen.ir.directives import Directives, ResultSizeIs
from codegen.ir.doc import Doc
from codegen.ir.roles import ModeValue, analyse, analyse_project
from codegen.ir.types import Intent

import builders as b


@pytest.fixture
def bag():
    return DiagnosticBag()


def analysed(procedure, bag):
    analyse(procedure, bag)
    return procedure


def mode_doc(*rows, header=("Mode", "Value"), preamble=("How to link.",)):
    lines = [*preamble, f"| {header[0]} | {header[1]} |", "|------|-------|"]
    lines.extend(f"| {description} | {value} |" for description, value in rows)
    return lines


def producer_input_doc(*rows):
    lines = ["Where the size comes from.",
             "| Producer input | Supplied by |",
             "|----------------|-------------|"]
    lines.extend(f"| {wanted} | {given} |" for wanted, given in rows)
    return lines


class TestTemporaries:
    def test_the_tmp_prefix_marks_a_work_array(self):
        procedure = analysed(
            b.procedure("p", b.real("tmp_work", Intent.OUT, "(n)"), b.real("v", Intent.IN, "(n)")),
            DiagnosticBag(),
        )

        assert procedure.argument("tmp_work").roles.is_temporary
        assert not procedure.argument("v").roles.is_temporary

    def test_a_near_miss_prefix_is_not_a_temporary(self):
        # 'temp_' is called out in issue #131 as the mistake to expect
        procedure = analysed(b.procedure("p", b.real("temp_work", Intent.OUT, "(n)")), DiagnosticBag())

        assert not procedure.argument("temp_work").roles.is_temporary


class TestExtents:
    def test_a_scalar_naming_an_extent_is_linked_to_its_array(self):
        procedure = analysed(
            b.procedure("p", b.integer("n_dims"), b.real("vector", Intent.INOUT, "(n_dims)")),
            DiagnosticBag(),
        )

        roles = procedure.argument("n_dims").roles
        assert roles.is_extent
        assert [a.name for a in roles.extent_of] == ["vector"]
        assert roles.is_derived

    def test_one_extent_can_serve_several_arrays(self):
        # this is the cross-check the binding languages must make: Fortran cannot
        procedure = analysed(
            b.procedure(
                "p",
                b.integer("n"),
                b.real("a", Intent.IN, "(n)"),
                b.real("b", Intent.IN, "(n)"),
            ),
            DiagnosticBag(),
        )

        assert [a.name for a in procedure.argument("n").roles.extent_of] == ["a", "b"]

    def test_an_extent_of_a_matrix_records_the_owner_once(self):
        procedure = analysed(
            b.procedure("p", b.integer("n"), b.real("m", Intent.IN, "(n, n)")),
            DiagnosticBag(),
        )

        assert [a.name for a in procedure.argument("n").roles.extent_of] == ["m"]

    def test_matching_is_case_insensitive(self):
        procedure = analysed(
            b.procedure("p", b.integer("N_Dims"), b.real("v", Intent.IN, "(n_dims)")),
            DiagnosticBag(),
        )

        assert procedure.argument("n_dims").roles.is_extent

    def test_a_scalar_naming_nothing_is_not_an_extent(self):
        procedure = analysed(b.procedure("p", b.integer("seed"), b.real("v", Intent.IN, "(n)")), DiagnosticBag())

        assert not procedure.argument("seed").roles.is_extent

    def test_a_non_integer_scalar_is_not_an_extent(self):
        procedure = analysed(
            b.procedure("p", b.real("n"), b.real("v", Intent.IN, "(n)")), DiagnosticBag()
        )

        assert not procedure.argument("n").roles.is_extent

    def test_an_array_is_not_an_extent(self):
        procedure = analysed(
            b.procedure("p", b.integer("n", dimension="(2)"), b.real("v", Intent.IN, "(n)")),
            DiagnosticBag(),
        )

        assert not procedure.argument("n").roles.is_extent

    def test_a_constant_extent_links_to_nothing(self):
        procedure = analysed(b.procedure("p", b.real("v", Intent.IN, "(3)")), DiagnosticBag())

        assert procedure.argument("v").roles.extent_of == ()


class TestShapeArguments:
    def test_a_shape_argument_is_linked_to_its_array(self):
        procedure = analysed(
            b.procedure(
                "p",
                b.real("data", Intent.IN, "(:)"),
                b.integer("data_shape", Intent.IN, "(:)"),
            ),
            DiagnosticBag(),
        )

        shape = procedure.argument("data_shape")
        data = procedure.argument("data")
        assert shape.roles.is_shape_arg
        assert shape.roles.shape_of is data
        assert data.roles.has_shape_arg
        assert data.roles.shape_arg is shape
        assert shape.roles.is_derived

    def test_a_shape_suffix_without_an_owner_is_not_a_shape_argument(self):
        procedure = analysed(b.procedure("p", b.integer("data_shape", Intent.IN, "(:)")), DiagnosticBag())

        assert not procedure.argument("data_shape").roles.is_shape_arg

    def test_an_array_without_a_shape_argument_has_none(self):
        procedure = analysed(b.procedure("p", b.real("data", Intent.IN, "(:)")), DiagnosticBag())

        assert not procedure.argument("data").roles.has_shape_arg


class TestMaskCounts:
    @pytest.mark.parametrize("mask_name", ["genes_mask", "genes_selection_mask"])
    def test_a_count_is_linked_to_its_mask(self, mask_name):
        procedure = analysed(
            b.procedure(
                "p",
                b.logical(mask_name, Intent.IN, "(n)"),
                b.integer("n_selected_genes", Intent.IN),
            ),
            DiagnosticBag(),
        )

        count = procedure.argument("n_selected_genes")
        mask = procedure.argument(mask_name)
        assert count.roles.is_mask_count
        assert count.roles.mask_count_of is mask
        assert mask.roles.count_arg is count
        assert count.roles.is_derived

    def test_a_count_without_a_mask_is_not_linked(self):
        procedure = analysed(b.procedure("p", b.integer("n_selected_genes")), DiagnosticBag())

        assert not procedure.argument("n_selected_genes").roles.is_mask_count

    def test_a_count_for_a_differently_named_mask_is_not_linked(self):
        procedure = analysed(
            b.procedure(
                "p",
                b.logical("families_mask", Intent.IN, "(n)"),
                b.integer("n_selected_genes"),
            ),
            DiagnosticBag(),
        )

        assert not procedure.argument("n_selected_genes").roles.is_mask_count

    def test_a_mask_without_a_count_has_none(self):
        procedure = analysed(b.procedure("p", b.logical("genes_mask", Intent.IN, "(n)")), DiagnosticBag())

        assert procedure.argument("genes_mask").roles.count_arg is None


class TestModeTable:
    def test_the_values_are_read_from_the_table(self, bag):
        procedure = analysed(
            b.procedure(
                "p",
                b.character(
                    "link_method",
                    Intent.IN,
                    doc=mode_doc(
                        ("minimises variance", "[[tox_clustering(module):METHOD_WARD(variable)]]"),
                        ("nearest neighbour", "[[tox_clustering(module):METHOD_SINGLE(variable)]]"),
                        header=("Method", "Value"),
                    ),
                ),
            ),
            bag,
        )

        mode = procedure.argument("link_method").roles.mode
        assert bag.errors == ()
        assert procedure.argument("link_method").roles.is_mode
        assert mode.alias == "method"
        assert mode.values == (
            ModeValue("METHOD_WARD", "tox_clustering", "ward", "minimises variance"),
            ModeValue("METHOD_SINGLE", "tox_clustering", "single", "nearest neighbour"),
        )

    def test_the_mode_string_is_the_parameter_without_its_prefix_in_lower_case(self, bag):
        procedure = analysed(
            b.procedure(
                "p",
                b.character(
                    "mode",
                    Intent.IN,
                    doc=mode_doc(("x", "[[m(module):MODE_GROUP_ORTHOLOGS(variable)]]")),
                ),
            ),
            bag,
        )

        assert procedure.argument("mode").roles.mode.values[0].string == "group_orthologs"

    def test_a_method_argument_expects_method_prefixed_parameters(self, bag):
        procedure = analysed(
            b.procedure(
                "p",
                b.character(
                    "link_method",
                    Intent.IN,
                    doc=mode_doc(
                        ("x", "[[m(module):METHOD_WARD(variable)]]"), header=("Method", "Value")
                    ),
                ),
            ),
            bag,
        )

        assert bag.errors == ()
        assert procedure.argument("link_method").roles.mode.alias == "method"

    def test_max_string_length_sets_the_character_length(self, bag):
        procedure = analysed(
            b.procedure(
                "p",
                b.character(
                    "mode",
                    Intent.IN,
                    doc=mode_doc(
                        ("a", "[[m(module):MODE_WARD(variable)]]"),
                        ("b", "[[m(module):MODE_COMPLETE(variable)]]"),
                    ),
                ),
            ),
            bag,
        )

        assert procedure.argument("mode").roles.mode.max_string_length == len("complete")

    def test_value_for_looks_a_string_up(self, bag):
        procedure = analysed(
            b.procedure(
                "p",
                b.character("mode", Intent.IN, doc=mode_doc(("a", "[[m(module):MODE_WARD(variable)]]"))),
            ),
            bag,
        )

        mode = procedure.argument("mode").roles.mode
        assert mode.value_for("ward").parameter == "MODE_WARD"
        assert mode.value_for("absent") is None

    def test_an_empty_value_cell_documents_a_mode_without_a_parameter(self, bag):
        procedure = analysed(
            b.procedure(
                "p",
                b.character(
                    "mode",
                    Intent.IN,
                    doc=mode_doc(
                        ("implemented", "[[m(module):MODE_WARD(variable)]]"),
                        ("not yet", ""),
                    ),
                ),
            ),
            bag,
        )

        assert bag.errors == ()
        assert len(procedure.argument("mode").roles.mode.values) == 1

    def test_a_table_elsewhere_in_the_doc_is_ignored(self, bag):
        doc = ["Notes.", "| a | b |", "|---|---|", "| 1 | 2 |", ""]
        doc += mode_doc(("x", "[[m(module):MODE_WARD(variable)]]"), preamble=())

        procedure = analysed(b.procedure("p", b.character("mode", Intent.IN, doc=doc)), bag)

        assert bag.errors == ()
        assert len(procedure.argument("mode").roles.mode.values) == 1

    def test_a_non_mode_argument_gets_no_mode_table(self, bag):
        procedure = analysed(b.procedure("p", b.integer("n_dims")), bag)

        assert procedure.argument("n_dims").roles.mode is None
        assert bag.errors == ()


class TestModeTableErrors:
    def test_a_mode_argument_without_a_table_is_an_error(self, bag):
        analyse(b.procedure("p", b.character("link_method", Intent.IN, doc="Just prose.")), bag)

        assert len(bag.errors) == 1
        assert "has no table of accepted values" in bag.errors[0].message
        assert "[[a_module(module):METHOD_BLA(variable)]]" in bag.errors[0].note

    def test_a_malformed_row_names_the_row_rather_than_disqualifying_the_table(self, bag):
        # the old parser bailed out of mode detection on a bad row, and the author was
        # told 'no mode table' -- pointing away from the actual mistake
        analyse(
            b.procedure(
                "p",
                b.character("mode", Intent.IN, doc=mode_doc(("x", "ward"))),
            ),
            bag,
        )

        assert len(bag.errors) == 1
        assert "row 1" in bag.errors[0].message
        assert "does not name a parameter" in bag.errors[0].message

    def test_a_value_cell_with_trailing_prose_is_rejected(self, bag):
        analyse(
            b.procedure(
                "p",
                b.character(
                    "mode", Intent.IN, doc=mode_doc(("x", "[[m(module):MODE_WARD(variable)]] maybe"))
                ),
            ),
            bag,
        )

        assert "does not name a parameter" in bag.errors[0].message

    def test_a_link_to_a_module_rather_than_a_parameter_is_rejected(self, bag):
        analyse(
            b.procedure(
                "p", b.character("mode", Intent.IN, doc=mode_doc(("x", "[[m(module)]]")))
            ),
            bag,
        )

        assert "links to module 'm' rather than to a parameter" in bag.errors[0].message

    @pytest.mark.parametrize("parameter", ["WARD", "MODE_ward", "METHOD_WARD", "MODEWARD"])
    def test_a_parameter_without_the_expected_prefix_is_rejected(self, bag, parameter):
        analyse(
            b.procedure(
                "p",
                b.character(
                    "mode", Intent.IN, doc=mode_doc(("x", f"[[m(module):{parameter}(variable)]]"))
                ),
            ),
            bag,
        )

        assert len(bag.errors) == 1
        assert f"'{parameter}' is not a valid mode parameter" in bag.errors[0].message

    def test_two_parameters_mapping_to_one_string_are_rejected(self, bag):
        analyse(
            b.procedure(
                "p",
                b.character(
                    "mode",
                    Intent.IN,
                    doc=mode_doc(
                        ("x", "[[a(module):MODE_WARD(variable)]]"),
                        ("y", "[[b(module):MODE_WARD(variable)]]"),
                    ),
                ),
            ),
            bag,
        )

        assert "maps ward more than once" in bag.errors[0].message

    def test_a_table_with_no_values_is_rejected(self, bag):
        analyse(b.procedure("p", b.character("mode", Intent.IN, doc=mode_doc())), bag)

        assert "has a table with no values" in bag.errors[0].message

    def test_two_mode_tables_are_rejected(self, bag):
        doc = mode_doc(("x", "[[m(module):MODE_A(variable)]]"))
        doc += ["", *mode_doc(("y", "[[m(module):MODE_B(variable)]]"), preamble=())]

        analyse(b.procedure("p", b.character("mode", Intent.IN, doc=doc)), bag)

        assert "has 2 tables of accepted values" in bag.errors[0].message

    def test_a_header_that_disagrees_with_the_argument_is_rejected(self, bag):
        # the old generator took the expected parameter prefix from the header, so a
        # 'link_method' argument headed 'Mode' with MODE_ parameters passed silently
        analyse(
            b.procedure(
                "p",
                b.character(
                    "link_method", Intent.IN, doc=mode_doc(("x", "[[m(module):MODE_WARD(variable)]]"))
                ),
            ),
            bag,
        )

        assert len(bag.errors) == 1
        assert "is headed 'Mode' but 'link_method' is a method argument" in bag.errors[0].message
        assert "METHOD_<NAME>" in bag.errors[0].note

    def test_an_error_points_at_the_offending_row(self, bag):
        doc = Doc.parse(mode_doc(("x", "ward")), first_line_number=100)
        argument = b.character("mode", Intent.IN)
        argument.doc = doc

        analyse(b.procedure("p", argument), bag)

        assert bag.errors[0].location.line is not None


class TestResultSize:
    def test_the_directive_resolves_to_the_counting_argument(self, bag):
        results = b.real("results", Intent.OUT, "(n)", directives=Directives(result_size_is=ResultSizeIs("n_results")))
        procedure = analysed(b.procedure("p", results, b.integer("n_results", Intent.OUT)), bag)

        assert bag.errors == ()
        assert procedure.argument("results").roles.result_size_arg is procedure.argument("n_results")

    def test_an_unknown_argument_is_reported(self, bag):
        results = b.real("results", Intent.OUT, "(n)", directives=Directives(result_size_is=ResultSizeIs("nope")))

        analyse(b.procedure("p", results), bag)

        assert "which is not an argument of 'p'" in bag.errors[0].message

    def test_a_non_scalar_counter_is_reported(self, bag):
        results = b.real("results", Intent.OUT, "(n)", directives=Directives(result_size_is=ResultSizeIs("counts")))

        analyse(b.procedure("p", results, b.integer("counts", Intent.OUT, "(2)")), bag)

        assert "a result size is a scalar integer" in bag.errors[0].message

    def test_a_non_integer_counter_is_reported(self, bag):
        results = b.real("results", Intent.OUT, "(n)", directives=Directives(result_size_is=ResultSizeIs("count")))

        analyse(b.procedure("p", results, b.real("count", Intent.OUT)), bag)

        assert "a result size is a scalar integer" in bag.errors[0].message

    def test_an_argument_without_the_directive_has_no_counter(self, bag):
        procedure = analysed(b.procedure("p", b.real("results", Intent.OUT, "(n)")), bag)

        assert procedure.argument("results").roles.result_size_arg is None


class TestOutputFrom:
    """DM_OUTPUT_FROM(AUTO): an argument obtained by calling another procedure."""

    def _consumer(self, directive_mode="AUTO", producer_inputs=("n",), doc=None):
        from codegen.ir.directives import Directives, OutputFrom, OutputFromMode

        mode = OutputFromMode.AUTO if directive_mode == "AUTO" else OutputFromMode.JUST_INFO
        work = b.integer(
            "n_work", Intent.IN,
            directives=Directives(output_from=OutputFrom("n_work", "sizer", "m", mode)),
            **({"doc": doc} if doc is not None else {}),
        )
        consumer = b.procedure("uses_work", b.integer("n"), work, b.ierr())
        producer = b.procedure(
            "sizer",
            *[b.integer(name) for name in producer_inputs],
            b.integer("n_work", Intent.OUT),
            b.ierr(),
        )
        return b.project(b.module("m", consumer, producer))

    def test_a_producer_is_resolved_and_inputs_name_matched(self, bag):
        project = self._consumer(producer_inputs=("n",))

        analyse_project(project, bag)

        plan = project.procedure("m", "uses_work").argument("n_work").roles.computed_from
        assert bag.errors == ()
        assert plan is not None
        assert plan.producer.name == "sizer"
        assert plan.output.name == "n_work"
        assert [(s.name, s.argument) for s in plan.inputs] == [("n", "n")]
        assert plan.is_automatic

    def test_an_auto_argument_is_derived_and_not_asked_of_the_caller(self, bag):
        project = self._consumer()

        analyse_project(project, bag)

        roles = project.procedure("m", "uses_work").argument("n_work").roles
        assert roles.is_computed
        assert roles.is_derived

    def test_just_info_is_not_computed_so_the_caller_still_supplies_it(self, bag):
        project = self._consumer(directive_mode="JUST_INFO")

        analyse_project(project, bag)

        roles = project.procedure("m", "uses_work").argument("n_work").roles
        assert roles.computed_from is not None
        assert not roles.is_computed
        assert not roles.is_derived

    def test_a_producer_input_with_no_matching_name_is_reported(self, bag):
        project = self._consumer(producer_inputs=("n_absent",))

        analyse_project(project, bag)

        assert any("n_absent" in e.message and "no argument of the same name" in e.message
                   for e in bag.errors)

    def test_just_info_does_not_require_the_consumer_to_feed_the_producer(self, bag):
        # a JUST_INFO pointer is documentation only: nothing is called, so the consumer
        # need not be able to supply the producer's inputs (the filter-mask case, where the
        # producer needs data the consumer never sees). AUTO in the same shape errors.
        project = self._consumer(directive_mode="JUST_INFO", producer_inputs=("n_absent",))

        analyse_project(project, bag)

        assert bag.errors == ()
        roles = project.procedure("m", "uses_work").argument("n_work").roles
        assert roles.computed_from is not None
        assert roles.computed_from.inputs == ()
        assert not roles.is_computed

    def test_a_producer_in_another_module_is_resolved(self, bag):
        from codegen.ir.directives import Directives, OutputFrom, OutputFromMode

        work = b.integer(
            "n_work", Intent.IN,
            directives=Directives(
                output_from=OutputFrom("n_work", "sizer", "other", OutputFromMode.AUTO)
            ),
        )
        consumer = b.procedure("uses_work", b.integer("n"), work, b.ierr())
        producer = b.procedure(
            "sizer", b.integer("n"), b.integer("n_work", Intent.OUT), b.ierr()
        )
        project = b.project(b.module("m", consumer), b.module("other", producer))

        analyse_project(project, bag)

        plan = project.procedure("m", "uses_work").argument("n_work").roles.computed_from
        assert bag.errors == ()
        assert plan is not None
        assert plan.producer.module.name == "other"

    def test_an_inout_producer_input_is_fed_back(self, bag):
        # the producer refines an input in place (caps it), so after the AUTO call the
        # consumer must adopt the returned value, not the raw one it passed in
        from codegen.ir.directives import Directives, OutputFrom, OutputFromMode

        work = b.integer(
            "n_work", Intent.IN,
            directives=Directives(output_from=OutputFrom("n_work", "sizer", "m",
                                                         OutputFromMode.AUTO)),
        )
        consumer = b.procedure("uses_work", b.integer("cap"), work, b.ierr())
        producer = b.procedure(
            "sizer", b.integer("cap", Intent.INOUT), b.integer("n_work", Intent.OUT), b.ierr()
        )
        project = b.project(b.module("m", consumer, producer))

        analyse_project(project, bag)

        plan = project.procedure("m", "uses_work").argument("n_work").roles.computed_from
        assert bag.errors == ()
        # cap is passed in AND fed back; a plain input is not
        assert [s.name for s in plan.inputs] == ["cap"]
        assert [(s.name, s.argument) for s in plan.inout_feedback] == [("cap", "cap")]

    def test_a_plain_input_is_not_fed_back(self, bag):
        project = self._consumer(producer_inputs=("n",))

        analyse_project(project, bag)

        plan = project.procedure("m", "uses_work").argument("n_work").roles.computed_from
        assert plan.inout_feedback == ()

    def test_a_table_supplies_a_renamed_producer_input(self, bag):
        project = self._consumer(
            producer_inputs=("n_values",),
            doc=producer_input_doc(("n_values", "n")),
        )

        analyse_project(project, bag)

        plan = project.procedure("m", "uses_work").argument("n_work").roles.computed_from
        assert bag.errors == ()
        assert plan is not None
        # the producer's own name on the left, the consumer's argument on the right
        assert [(s.name, s.argument) for s in plan.inputs] == [("n_values", "n")]

    def test_a_table_supplies_a_constant_the_consumer_has_no_argument_for(self, bag):
        project = self._consumer(
            producer_inputs=("n", "factor"),
            doc=producer_input_doc(("factor", "3_int32")),
        )

        analyse_project(project, bag)

        plan = project.procedure("m", "uses_work").argument("n_work").roles.computed_from
        assert bag.errors == ()
        assert plan is not None
        by_name = {supply.name: supply for supply in plan.inputs}
        assert by_name["n"].argument == "n"
        assert by_name["factor"].is_constant
        assert by_name["factor"].constant == 3

    def test_a_cell_that_is_neither_argument_nor_constant_is_reported(self, bag):
        project = self._consumer(
            producer_inputs=("n", "factor"),
            doc=producer_input_doc(("factor", "not_a_thing")),
        )

        analyse_project(project, bag)

        assert any("factor" in e.message and "no argument of the same name" in e.message
                   for e in bag.errors)

    def test_a_table_naming_an_input_the_producer_lacks_is_reported(self, bag):
        project = self._consumer(
            producer_inputs=("n",), doc=producer_input_doc(("n_absent", "n")))

        analyse_project(project, bag)

        assert any("not an argument of 'sizer'" in e.message for e in bag.errors)

    def test_an_empty_cell_in_the_table_is_reported(self, bag):
        project = self._consumer(
            producer_inputs=("n_values",), doc=producer_input_doc(("n_values", "")))

        analyse_project(project, bag)

        assert any("empty cell" in e.message for e in bag.errors)

    def test_a_missing_producer_is_reported(self, bag):
        from codegen.ir.directives import Directives, OutputFrom, OutputFromMode

        work = b.integer(
            "n_work", Intent.IN,
            directives=Directives(
                output_from=OutputFrom("n_work", "nope", "m", OutputFromMode.AUTO)
            ),
        )
        project = b.project(b.module("m", b.procedure("p", work, b.ierr())))

        analyse_project(project, bag)

        assert any("does not exist" in e.message for e in bag.errors)

    def test_a_non_exported_producer_is_reported(self, bag):
        from codegen.ir.directives import Directives, OutputFrom, OutputFromMode
        from codegen.ir.entities import Meta

        work = b.integer(
            "n_work", Intent.IN,
            directives=Directives(
                output_from=OutputFrom("n_work", "sizer", "m", OutputFromMode.AUTO)
            ),
        )
        consumer = b.procedure("uses_work", b.integer("n"), work, b.ierr())
        producer = b.procedure(
            "sizer", b.integer("n"), b.integer("n_work", Intent.OUT), b.ierr(), meta=Meta()
        )
        project = b.project(b.module("m", consumer, producer))

        analyse_project(project, bag)

        assert any("not exported" in e.message for e in bag.errors)


class TestResultsAndProjects:
    def test_a_function_result_is_analysed_too(self, bag):
        result = b.real("tmp_out", Intent.OUT, is_result=True)

        analyse(b.procedure("f", result=result), bag)

        assert result.roles is not None
        assert result.roles.is_temporary

    def test_an_extent_can_belong_to_the_result(self, bag):
        result = b.real("out", Intent.OUT, "(n)", is_result=True)
        procedure = b.procedure("f", b.integer("n"), result=result)

        analyse(procedure, bag)

        assert [a.name for a in procedure.argument("n").roles.extent_of] == ["out"]

    def test_analyse_project_walks_every_procedure(self, bag):
        project = b.project(
            b.module("a", b.procedure("p", b.integer("n"), b.real("v", Intent.IN, "(n)"))),
            b.module("b", b.procedure("q", b.real("tmp_w", Intent.OUT, "(2)"))),
        )

        analyse_project(project, bag)

        assert project.procedure("a", "p").argument("n").roles.is_extent
        assert project.procedure("b", "q").argument("tmp_w").roles.is_temporary

    def test_every_argument_gets_roles(self, bag):
        procedure = analysed(b.procedure("p", b.integer("n"), b.real("v", Intent.IN, "(n)"), b.ierr()), bag)

        assert all(argument.roles is not None for argument in procedure)
