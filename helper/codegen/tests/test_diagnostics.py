from pathlib import Path

import pytest

from codegen.diagnostics import (
    CodegenError,
    Context,
    Diagnostic,
    DiagnosticBag,
    Severity,
    SourceLocation,
)
from codegen.ir.entities import Module, Project


class Entity:
    """Stand-in for an IR entity: the chain walk only needs kind/name/parent."""

    def __init__(self, kind, name, parent=None, location=None):
        self.entity_kind = kind
        self.name = name
        self.parent = parent
        self.location = location


@pytest.fixture
def argument():
    module = Entity("module", "tox_clustering")
    procedure = Entity("procedure", "cluster_alloc", parent=module)
    return Entity(
        "argument",
        "link_method",
        parent=procedure,
        location=SourceLocation(Path("src/tox_clustering.F90"), 42),
    )


def test_source_location_renders_as_clickable_file_line():
    assert str(SourceLocation(Path("src/a.F90"), 7)) == "src/a.F90:7"


def test_source_location_without_line_renders_file_only():
    assert str(SourceLocation(Path("src/a.F90"))) == "src/a.F90"


def test_source_location_without_file_is_unknown():
    assert str(SourceLocation()) == "<unknown>"


def test_context_walks_parents_innermost_first(argument):
    assert Context.of(argument).entities == (
        ("argument", "link_method"),
        ("procedure", "cluster_alloc"),
        ("module", "tox_clustering"),
    )


def test_context_renders_as_an_in_chain(argument):
    assert str(Context.of(argument)) == (
        "argument 'link_method' in procedure 'cluster_alloc' in module 'tox_clustering'"
    )


def test_context_skips_entities_without_a_kind():
    anonymous = Entity("module", "m")
    del anonymous.entity_kind
    child = Entity("procedure", "p", parent=anonymous)

    assert Context.of(child).entities == (("procedure", "p"),)


def test_context_names_the_project_it_reaches():
    """A real diagnostic is raised on a module the project has adopted, not a loose one.

    The name is `[extra.ford] project` in fpm.toml, which the Ford frontend passes on.
    """
    module = Module(name="tox_thing")
    Project([module], name="TensorOmics")

    assert str(Context.of(module)) == "module 'tox_thing' in project 'TensorOmics'"


def test_context_skips_an_unnamed_project():
    """Hand-built IR has no project name, and must not be reported `in project ''`.

    Which is what every module-level diagnostic said until the frontend was taught to pass
    the name on -- `Project.name` existed and no caller had ever filled it.
    """
    module = Module(name="tox_thing")
    Project([module])

    assert str(Context.of(module)) == "module 'tox_thing'"


def test_empty_context_is_falsy():
    assert not Context()


def test_error_takes_location_and_context_from_the_entity(argument):
    bag = DiagnosticBag()
    bag.error("no mode table", entity=argument)

    diagnostic = bag.errors[0]
    assert diagnostic.location == SourceLocation(Path("src/tox_clustering.F90"), 42)
    assert diagnostic.context.entities[0] == ("argument", "link_method")


def test_explicit_location_wins_over_the_entity(argument):
    bag = DiagnosticBag()
    bag.error("boom", entity=argument, location=SourceLocation(Path("other.F90"), 1))

    assert bag.errors[0].location.file == Path("other.F90")


def test_entity_without_location_yields_unknown_location():
    bag = DiagnosticBag()
    bag.error("boom", entity=Entity("module", "m"))

    assert bag.errors[0].location == SourceLocation()


def test_bag_separates_errors_from_warnings(argument):
    bag = DiagnosticBag()
    bag.warn("no summary", entity=argument)
    bag.error("no mode table", entity=argument)

    assert len(bag) == 2
    assert len(bag.errors) == 1
    assert len(bag.warnings) == 1
    assert bag.errors[0].severity is Severity.ERROR
    assert bag.warnings[0].severity is Severity.WARNING


def test_diagnostics_keep_insertion_order():
    bag = DiagnosticBag()
    bag.error("first")
    bag.warn("second")
    bag.error("third")

    assert [d.message for d in bag] == ["first", "second", "third"]


def test_extend_merges_another_bag():
    bag = DiagnosticBag()
    bag.error("a")
    other = DiagnosticBag()
    other.error("b")
    bag.extend(other)

    assert [d.message for d in bag] == ["a", "b"]


def test_render_shows_message_location_context_and_note(argument):
    diagnostic = Diagnostic(
        severity=Severity.ERROR,
        message="mode argument has no mode table",
        location=SourceLocation(Path("src/tox_clustering.F90"), 42),
        context=Context.of(argument),
        note="expected a markdown table",
    )

    assert diagnostic.render() == (
        "error: mode argument has no mode table\n"
        "  --> src/tox_clustering.F90:42\n"
        "  argument 'link_method' in procedure 'cluster_alloc' in module 'tox_clustering'\n"
        "  note: expected a markdown table"
    )


def test_render_omits_absent_parts():
    assert Diagnostic(Severity.WARNING, "no author").render() == "warning: no author"


def test_render_indents_multiline_notes():
    diagnostic = Diagnostic(Severity.ERROR, "bad table", note="line one\nline two")

    assert diagnostic.render() == "error: bad table\n  note: line one\n  line two"


def test_render_is_uncoloured_by_default_and_coloured_on_request():
    diagnostic = Diagnostic(Severity.ERROR, "boom")

    assert "\033[" not in diagnostic.render()
    assert "\033[" in diagnostic.render(color=True)


def test_raise_if_errors_is_silent_when_there_are_only_warnings():
    bag = DiagnosticBag()
    bag.warn("no summary")

    bag.raise_if_errors()  # must not raise: a missing docstring is not a build failure


def test_raise_if_errors_reports_every_error_at_once():
    bag = DiagnosticBag()
    bag.error("first problem")
    bag.error("second problem")

    with pytest.raises(CodegenError) as excinfo:
        bag.raise_if_errors()

    message = str(excinfo.value)
    assert "2 errors" in message
    assert "first problem" in message
    assert "second problem" in message
    assert len(excinfo.value.diagnostics) == 2


def test_raise_if_errors_uses_singular_for_one_error():
    bag = DiagnosticBag()
    bag.error("only problem")

    with pytest.raises(CodegenError, match="1 error in Fortran sources"):
        bag.raise_if_errors()
