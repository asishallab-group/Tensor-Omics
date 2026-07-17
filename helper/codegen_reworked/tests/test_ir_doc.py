import pytest

from codegen_reworked.ir.doc import (
    Alignment,
    Doc,
    DocLine,
    DocParseError,
    DocTable,
    FordLink,
    Text,
)


class TestFordLink:
    @pytest.mark.parametrize(
        "text, component, component_type, item, item_type",
        [
            ("[[mod]]", "mod", None, None, None),
            ("[[mod(module)]]", "mod", "module", None, None),
            ("[[mod:sub]]", "mod", None, "sub", None),
            ("[[mod(module):sub(subroutine)]]", "mod", "module", "sub", "subroutine"),
            ("[[mod(module):MODE_X(variable)]]", "mod", "module", "MODE_X", "variable"),
            ("[[ mod (module) : sub (function) ]]", "mod", "module", "sub", "function"),
        ],
    )
    def test_parse(self, text, component, component_type, item, item_type):
        link = DocLine.parse(text).links[0]

        assert link.component == component
        assert link.component_type == component_type
        assert link.item == item
        assert link.item_type == item_type

    @pytest.mark.parametrize(
        "text",
        [
            "[[mod(module):sub(subroutine)]]",
            "[[mod(module)]]",
            "[[mod:sub]]",
            "[[mod]]",
        ],
    )
    def test_str_round_trips(self, text):
        assert str(DocLine.parse(text).links[0]) == text

    def test_target_prefers_the_item(self):
        assert FordLink("mod", "module", "sub", "subroutine").target == "sub"
        assert FordLink("mod", "module").target == "mod"

    def test_an_unknown_entity_type_is_not_a_link(self):
        assert DocLine.parse("[[mod(banana)]]").links == ()

    def test_links_compare_by_value(self):
        assert FordLink("mod", "module") == FordLink("mod", "module")


class TestDocLine:
    def test_plain_prose_is_a_single_text_span(self):
        line = DocLine.parse("just prose")

        assert line.spans == (Text("just prose"),)
        assert line.text == "just prose"

    def test_text_around_a_link_is_preserved(self):
        line = DocLine.parse("see [[mod(module):sub]] for details")

        assert [type(span).__name__ for span in line.spans] == ["Text", "FordLink", "Text"]
        assert line.text == "see [[mod(module):sub]] for details"

    def test_several_links_in_one_line(self):
        line = DocLine.parse("[[a(module)]] and [[b(module)]]")

        assert len(line.links) == 2
        assert line.links[0].component == "a"
        assert line.links[1].component == "b"

    def test_a_line_that_is_only_a_link(self):
        line = DocLine.parse("[[mod(module):MODE_X(variable)]]")

        assert len(line.spans) == 1
        assert isinstance(line.spans[0], FordLink)

    def test_parse_strips_surrounding_whitespace(self):
        assert DocLine.parse("   prose   ").text == "prose"

    def test_empty_line_has_no_spans_and_is_blank(self):
        line = DocLine.parse("")

        assert line.spans == ()
        assert line.is_blank
        assert line.text == ""

    def test_whitespace_only_line_is_blank(self):
        assert DocLine.parse("    ").is_blank

    def test_line_numbers_are_recorded(self):
        assert DocLine.parse("x", line_number=42).line_number == 42

    def test_with_spans_returns_a_new_line_and_keeps_the_number(self):
        line = DocLine.parse("old", line_number=7)

        replaced = line.with_spans([Text("new")])

        assert replaced.text == "new"
        assert replaced.line_number == 7
        assert line.text == "old"


class TestAlignment:
    @pytest.mark.parametrize(
        "text, expected",
        [
            ("---", Alignment.DEFAULT),
            (":---", Alignment.LEFT),
            ("---:", Alignment.RIGHT),
            (":---:", Alignment.CENTER),
            ("-", Alignment.DEFAULT),
            (":-:", Alignment.CENTER),
            ("  ---  ", Alignment.DEFAULT),
        ],
    )
    def test_parse(self, text, expected):
        assert Alignment.parse(text) is expected

    def test_default_and_left_stay_distinct(self):
        # so a table survives a parse/render round trip unchanged
        assert Alignment.parse("---") is not Alignment.parse(":---")

    @pytest.mark.parametrize("text", ["---", ":---", "---:", ":---:"])
    def test_marker_round_trips(self, text):
        assert Alignment.parse(text).marker(width=len(text)) == text

    @pytest.mark.parametrize(
        "alignment, expected",
        [
            (Alignment.DEFAULT, "---"),
            (Alignment.LEFT, ":--"),
            (Alignment.RIGHT, "--:"),
            (Alignment.CENTER, ":-:"),
        ],
    )
    def test_marker_defaults_to_width_three(self, alignment, expected):
        assert alignment.marker() == expected

    def test_marker_never_collapses_to_nothing(self):
        assert Alignment.CENTER.marker(width=1) == ":-:"

    @pytest.mark.parametrize("text", ["", "abc", "::", "-a-"])
    def test_parse_rejects_non_alignments(self, text):
        with pytest.raises(DocParseError, match="not a markdown table alignment"):
            Alignment.parse(text)


class TestDocTableParsing:
    def test_a_table_is_recognised(self):
        doc = Doc.parse(
            [
                "| Mode | Value |",
                "|------|-------|",
                "| ward | [[mod(module):METHOD_WARD(variable)]] |",
            ]
        )

        assert len(doc.tables) == 1
        table = doc.tables[0]
        assert table.header_text == ("Mode", "Value")
        assert table.n_columns == 2
        assert len(table.rows) == 1
        assert table.rows[0][0].text == "ward"
        assert table.rows[0][1].links[0].item == "METHOD_WARD"

    def test_unpadded_pipes_are_accepted(self):
        # the old parser split on ' | ' and required padding, erroring out without it
        doc = Doc.parse(["|Mode|Value|", "|----|-----|", "|ward|1|"])

        assert doc.tables[0].header_text == ("Mode", "Value")
        assert doc.tables[0].rows[0][1].text == "1"

    def test_a_table_may_have_no_rows(self):
        doc = Doc.parse(["| Mode | Value |", "|------|-------|"])

        assert doc.tables[0].rows == ()

    def test_alignments_are_captured_per_column(self):
        doc = Doc.parse(["| a | b | c | d |", "|---|:--|--:|:-:|", "| 1 | 2 | 3 | 4 |"])

        assert doc.tables[0].alignment == (
            Alignment.DEFAULT,
            Alignment.LEFT,
            Alignment.RIGHT,
            Alignment.CENTER,
        )

    def test_prose_before_and_after_a_table_is_kept(self):
        doc = Doc.parse(["intro", "| a |", "|---|", "| 1 |", "outro"])

        assert len(doc.blocks) == 3
        assert doc.blocks[0].text == "intro"
        assert isinstance(doc.blocks[1], DocTable)
        assert doc.blocks[2].text == "outro"

    def test_two_tables_separated_by_prose(self):
        doc = Doc.parse(["| a |", "|---|", "| 1 |", "gap", "| b |", "|---|", "| 2 |"])

        assert len(doc.tables) == 2

    def test_a_row_shaped_line_without_an_alignment_row_stays_prose(self):
        # a lone '| note |' in a comment is not a mistake and must not fail generation
        doc = Doc.parse(["| just prose |"])

        assert doc.tables == ()
        assert doc.blocks[0].text == "| just prose |"

    def test_a_row_shaped_last_line_stays_prose(self):
        doc = Doc.parse(["text", "| trailing |"])

        assert doc.tables == ()

    def test_a_table_ends_at_the_first_non_row_line(self):
        doc = Doc.parse(["| a |", "|---|", "| 1 |", "after", "| 2 |"])

        assert len(doc.tables[0].rows) == 1
        assert doc.blocks[1].text == "after"

    def test_table_line_numbers_are_recorded(self):
        doc = Doc.parse(["| a |", "|---|", "| 1 |"], first_line_number=10)

        assert doc.tables[0].line_number == 10
        assert doc.tables[0].rows[0][0].line_number == 12


class TestDocTableValidation:
    def test_a_row_with_too_few_columns_is_rejected(self):
        with pytest.raises(DocParseError, match="row 1 has 1 columns, expected 2"):
            Doc.parse(["| a | b |", "|---|---|", "| 1 |"])

    def test_a_row_with_too_many_columns_is_rejected(self):
        with pytest.raises(DocParseError, match="row 2 has 3 columns, expected 2"):
            Doc.parse(["| a | b |", "|---|---|", "| 1 | 2 |", "| 1 | 2 | 3 |"])

    def test_a_mismatched_alignment_row_is_rejected(self):
        with pytest.raises(DocParseError, match="2 header columns but 3 alignment columns"):
            Doc.parse(["| a | b |", "|---|---|---|", "| 1 | 2 |"])

    def test_the_error_carries_a_line_number_for_placement(self):
        with pytest.raises(DocParseError) as excinfo:
            Doc.parse(["| a | b |", "|---|---|", "| 1 |"], first_line_number=5)

        assert excinfo.value.line_number == 5

    def test_a_table_needs_columns(self):
        with pytest.raises(DocParseError, match="table has no columns"):
            DocTable(header=(), alignment=())


class TestDocTableAccess:
    @pytest.fixture
    def table(self):
        return Doc.parse(
            ["| Mode | Value |", "|------|-------|", "| a | 1 |", "| b | 2 |"]
        ).tables[0]

    def test_column_extracts_cells(self, table):
        assert [cell.text for cell in table.column(0)] == ["a", "b"]
        assert [cell.text for cell in table.column(1)] == ["1", "2"]

    def test_with_rows_returns_a_new_table(self, table):
        replaced = table.with_rows([(DocLine.parse("x"), DocLine.parse("y"))])

        assert len(replaced.rows) == 1
        assert len(table.rows) == 2
        assert replaced.header == table.header

    def test_with_rows_still_validates_column_counts(self, table):
        with pytest.raises(DocParseError, match="expected 2"):
            table.with_rows([(DocLine.parse("x"),)])

    def test_a_table_is_never_blank(self, table):
        assert not table.is_blank


class TestDoc:
    def test_trailing_blank_lines_are_dropped(self):
        doc = Doc.parse(["text", "", "  ", ""])

        assert len(doc.blocks) == 1

    def test_interior_blank_lines_are_kept(self):
        doc = Doc.parse(["a", "", "b"])

        assert len(doc.blocks) == 3
        assert doc.blocks[1].is_blank

    def test_empty_doc_is_falsy(self):
        assert not Doc.parse([])
        assert not Doc.parse(["", "  "])

    def test_doc_with_content_is_truthy(self):
        assert Doc.parse(["text"])

    def test_summary_is_the_first_non_blank_line(self):
        assert Doc.parse(["", "first", "second"]).summary == "first"

    def test_summary_of_an_empty_doc_is_empty(self):
        assert Doc.parse([]).summary == ""

    def test_text_joins_blocks_by_line(self):
        assert Doc.parse(["a", "b"]).text == "a\nb"

    def test_text_renders_a_table_with_its_alignment_row(self):
        text = Doc.parse(["| a | b |", "|---|---|", "| 1 | 2 |"]).text

        assert text == "| a | b |\n|---|---|\n| 1 | 2 |"

    def test_a_rendered_table_parses_back_into_the_same_table(self):
        # an emitter re-emitting docs must not turn a table into something that no
        # longer parses as one
        source = ["| Mode | Value |", "|:-----|------:|", "| ward | 1 |"]
        original = Doc.parse(source)

        reparsed = Doc.parse(original.text.split("\n"))

        assert len(reparsed.tables) == 1
        assert reparsed.tables[0].header_text == original.tables[0].header_text
        assert reparsed.tables[0].alignment == original.tables[0].alignment
        assert reparsed.text == original.text

    def test_a_rendered_doc_round_trips_prose_and_tables_together(self):
        source = ["intro", "", "| a | b |", "|---|---|", "| 1 | 2 |", "", "outro"]

        once = Doc.parse(source).text

        assert Doc.parse(once.split("\n")).text == once

    def test_links_are_collected_from_prose_and_tables(self):
        doc = Doc.parse(
            ["see [[a(module)]]", "| x | [[b(module)]] |", "|---|---|", "| 1 | [[c(module)]] |"]
        )

        assert {link.component for link in doc.links} == {"a", "b", "c"}

    def test_parse_does_not_mutate_the_input(self):
        # the old DocList popped trailing entries off the caller's list
        lines = ["text", "", ""]

        Doc.parse(lines)

        assert lines == ["text", "", ""]

    def test_replace_swaps_a_block_and_leaves_the_original_alone(self):
        doc = Doc.parse(["a", "b"])
        original_block = doc.blocks[0]

        replaced = doc.replace(original_block, DocLine.parse("c"))

        assert replaced.text == "c\nb"
        assert doc.text == "a\nb"

    def test_replace_of_an_absent_block_changes_nothing(self):
        doc = Doc.parse(["a"])

        assert doc.replace(DocLine.parse("zzz"), DocLine.parse("c")).text == "a"

    def test_len_and_iteration_walk_blocks(self):
        doc = Doc.parse(["a", "b"])

        assert len(doc) == 2
        assert [block.text for block in doc] == ["a", "b"]

    def test_lines_and_tables_partition_the_blocks(self):
        doc = Doc.parse(["intro", "| a |", "|---|", "| 1 |"])

        assert len(doc.lines) == 1
        assert len(doc.tables) == 1
        assert len(doc.blocks) == 2
