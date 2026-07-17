from codegen_reworked.render import Writer


def test_lines_are_indented_by_the_current_level():
    w = Writer()
    w.line("module foo")
    with w.indent():
        w.line("implicit none")
    w.line("end module foo")

    assert w.render() == "module foo\n    implicit none\nend module foo"


def test_indent_nests_and_unwinds():
    w = Writer()
    with w.indent():
        w.line("a")
        with w.indent(2):
            w.line("b")
        w.line("c")
    w.line("d")

    assert w.render() == "    a\n            b\n    c\nd"


def test_indent_unwinds_on_exception():
    w = Writer()
    try:
        with w.indent():
            raise RuntimeError("boom")
    except RuntimeError:
        pass
    w.line("a")

    assert w.render() == "a"


def test_blank_lines_carry_no_trailing_whitespace():
    w = Writer()
    with w.indent():
        w.line("a")
        w.blank()
        w.line("b")

    assert w.render() == "    a\n\n    b"


def test_consecutive_blanks_collapse():
    w = Writer()
    w.line("a")
    w.blank()
    w.blank()
    w.blank()
    w.line("b")

    assert w.render() == "a\n\nb"


def test_leading_blanks_are_dropped():
    w = Writer()
    w.blank()
    w.line("a")

    assert w.render() == "a"


def test_collapse_can_be_opted_out_of():
    w = Writer()
    w.line("a")
    w.blank()
    w.blank(collapse=False)
    w.line("b")

    assert w.render() == "a\n\n\nb"


def test_trailing_blanks_are_stripped_by_render():
    w = Writer()
    w.line("a")
    w.blank(collapse=False)
    w.blank(collapse=False)

    assert w.render() == "a"


def test_block_preserves_relative_indentation():
    w = Writer()
    with w.indent():
        w.block("if (x) then\n    y = 1\nend if")

    assert w.render() == "    if (x) then\n        y = 1\n    end if"


def test_block_does_not_indent_its_blank_lines():
    w = Writer()
    with w.indent():
        w.block("a\n\nb")

    assert w.render() == "    a\n\n    b"


def test_line_accepts_multiline_text():
    w = Writer()
    with w.indent():
        w.line("a\nb")

    assert w.render() == "    a\n    b"


def test_extend_splices_at_the_current_indent():
    inner = Writer()
    inner.line("a")
    with inner.indent():
        inner.line("b")

    outer = Writer()
    with outer.indent():
        outer.extend(inner)

    assert outer.render() == "    a\n        b"


def test_render_strips_per_line_trailing_whitespace():
    w = Writer()
    w.line("a   ")

    assert w.render() == "a"


def test_trailing_newline_is_opt_in():
    w = Writer()
    w.line("a")

    assert w.render() == "a"
    assert w.render(trailing_newline=True) == "a\n"


def test_trailing_newline_on_empty_output_stays_empty():
    assert Writer().render(trailing_newline=True) == ""


def test_empty_writer_is_falsy_and_renders_empty():
    w = Writer()

    assert not w
    assert w.render() == ""


def test_writer_with_only_blanks_is_falsy():
    w = Writer()
    w.blank(collapse=False)

    assert not w


def test_custom_indent_string():
    w = Writer(indent="\t")
    with w.indent():
        w.line("a")

    assert w.render() == "\ta"


def test_lines_writes_each_item():
    w = Writer()
    w.lines(["a", "b"])

    assert w.render() == "a\nb"
