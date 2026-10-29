"""An argument is refused exactly where its generated Fortran wrapper uses the same name.

`emit.name_clashes` reads the wrappers as emitted, so these drive the whole pipeline from
hand-built projects (`generate(parsed=...)`) with the real macro header, and check both
directions: the name is refused where the wrapper calls, imports, uses as a kind or declares
it, and accepted where it does not -- `count` stays a good argument of a procedure that counts
nothing, as main's `mask_chunk_count` has it.
"""

import pytest

from codegen.frontend.ford_frontend import ParsedProject
from codegen.frontend.macros import MacroTable
from codegen.generate import generate
from codegen.ir.directives import Default, Directives, Prologue
from codegen.ir.entities import Meta, Project
from codegen.ir.types import Intent

from builders import C_BINDING, character, ierr, integer, logical, module, procedure, real
from conftest import REPO_ROOT

MACROS = MacroTable(REPO_ROOT / "src/macros.h", include_paths=(REPO_ROOT,))
IMPL = Meta(summary="s", author="AUTHOR")


def errors(*modules):
    parsed = ParsedProject(project=Project(list(modules)), macros=MACROS, arg_pos_factor=10000)
    return [d.message for d in generate(targets=("c",), parsed=parsed).diagnostics.errors]


def impl(*arguments, name="pick_impl", **kwargs):
    return module("tox_demo_impl", procedure(name, *arguments, meta=IMPL, **kwargs))


def masked(extra):
    return impl(
        integer("n_genes"), real("values", Intent.IN, "(n_genes)"),
        logical("genes_mask", Intent.IN, "(n_genes)", kind="c_bool"),
        integer("n_selected_genes"), extra,
    )


class TestAnIntrinsicTheWrapperCalls:
    def test_count_beside_a_mask_count_is_refused(self):
        assert errors(masked(integer("count", Intent.OUT))) == [
            "argument 'count' is a name the generated Fortran wrapper 'pick' uses: it calls "
            "'count(...)'"
        ]

    def test_count_where_nothing_is_counted_is_accepted(self):
        # main publishes `mask_chunk_count(n_genes, count)`
        assert errors(impl(integer("n"), real("values", Intent.IN, "(n)"),
                           integer("count", Intent.OUT))) == []

    def test_size_beside_an_assumed_shape_is_refused(self):
        assert errors(impl(real("values", Intent.IN, "(:)"), integer("size", Intent.OUT))) == [
            "argument 'size' is a name the generated Fortran wrapper 'pick' uses: it calls "
            "'size(...)'"
        ]


class TestAnImportOrAKind:
    def test_an_imported_helper_is_refused(self):
        assert errors(masked(integer("set_err_once", Intent.OUT))) == [
            "argument 'set_err_once' is a name the generated Fortran wrapper 'pick' uses: it "
            "calls 'set_err_once'"
        ]

    def test_a_kind_the_wrapper_writes_is_refused(self):
        # every validator call carries `arg_pos=<k>_int32`, and the module imports the kind
        found = errors(impl(integer("n"), real("values", Intent.IN, "(n)"),
                            integer("int32", Intent.OUT)))
        assert "argument 'int32' is a name the generated Fortran wrapper 'pick' uses: it " \
               "imports 'int32'" in found

    def test_a_kind_used_only_in_the_body_is_refused(self):
        # `c_bool` reaches the wrapper's body as `.true._c_bool`, whatever it imports
        from codegen.emit.name_clashes import _role

        body = "call validate(x, allow_nan=.true._c_bool)"
        assert _role("c_bool", True, body, set()) == "uses 'c_bool' as a kind"
        assert _role("c_bool", True, "integer(c_bool) :: y", set()) == "uses 'c_bool' as a kind"
        assert _role("c_bool", True, "count(m, kind=c_bool)", set()) == "uses 'c_bool' as a kind"
        assert _role("int32", True, "x1_int32 = 1", set()) is None


class TestALocalTheWrapperDeclares:
    def test_handled_beside_a_prologue_is_refused(self):
        guard = procedure(
            "guard", real("values", Intent.IN, "(n)"), integer("n"),
            logical("handled", Intent.OUT), ierr(), meta=IMPL,
        )
        crunch = procedure(
            "crunch_impl", real("values", Intent.IN, "(n)"), integer("n"),
            real("tmp_scratch", Intent.OUT, "(n)"), real("handled", Intent.OUT), meta=IMPL,
            directives=Directives(prologue=Prologue("guard", "tox_demo_impl")),
        )
        found = errors(module("tox_demo_impl", guard, crunch))
        assert "argument 'handled' is a name the generated Fortran wrapper 'crunch' uses: " \
               "it declares a local called 'handled'" in found

    def test_handled_without_a_prologue_is_accepted(self):
        assert errors(impl(integer("n"), real("values", Intent.IN, "(n)"),
                           real("handled", Intent.OUT))) == []

    def test_a_materialised_default_s_local_is_refused(self):
        from test_synthesize import optional_producer_input_impl_module

        source = optional_producer_input_impl_module()
        crunch = source.procedure("crunch_impl")
        crunch.arguments = (*crunch.arguments, real("exact_value", Intent.OUT))
        found = errors(source)
        assert "argument 'exact_value' is a name the generated Fortran wrapper 'crunch' uses: " \
               "it declares a local called 'exact_value'" in found

    def test_a_converted_argument_s_local_is_refused(self):
        found = errors(module("tox_demo", procedure(
            "p", character("word", Intent.IN, length="8"), integer("word_f", Intent.OUT), ierr(),
        )))
        assert found == [
            "argument 'word_f' is a name the generated C wrapper 'p_c' uses: it declares a local "
            "called 'word_f'"
        ]

    def test_an_array_of_strings_view_address_is_refused(self):
        # #221 B3: an array of strings is viewed at an address chosen before the remap
        found = errors(module("tox_demo", procedure(
            "p", integer("n"), character("words", Intent.IN, "(n)", length="8"),
            integer("words_address", Intent.OUT), ierr(),
        )))
        assert found == [
            "argument 'words_address' is a name the generated C wrapper 'p_c' uses: it declares "
            "a local called 'words_address'"
        ]

    def test_the_empty_view_anchor_is_refused(self):
        # a module variable, so a dummy of that name would hide it in the wrapper
        found = errors(module("tox_demo", procedure(
            "p", integer("n"), character("words", Intent.IN, "(n)", length="8"),
            integer("empty_view_anchor", Intent.OUT), ierr(),
        )))
        assert found == [
            "argument 'empty_view_anchor' is a name the generated C wrapper 'p_c' uses: it "
            "declares a module variable called 'empty_view_anchor'"
        ]

    def test_the_anchor_s_name_is_free_where_no_array_of_strings_is_viewed(self):
        assert errors(module("tox_demo", procedure(
            "p", character("word", Intent.IN, length="8"),
            integer("empty_view_anchor", Intent.OUT), integer("word_address", Intent.OUT),
            ierr(),
        ))) == []

    @pytest.mark.parametrize(
        "base, name",
        [("y", "y_f"), ("p", "p_value"), ("window", "window_size"), ("x", "x_out"),
         ("k", "k_v"), ("x", "x_dim")],
    )
    def test_a_name_whose_local_is_not_emitted_is_accepted(self, base, name):
        assert errors(module("tox_demo", procedure(
            "probe", real(base, Intent.IN), real(name, Intent.IN), real("s", Intent.OUT), ierr(),
            meta=C_BINDING,
        ))) == []
