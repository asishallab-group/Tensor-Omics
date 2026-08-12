"""Resolving a Ford cross-reference for the language that will read it.

A link is markup no target language renders, pointing at a Fortran symbol its callers cannot
reach. These check that it resolves to something they *can* use, and that where it cannot, no
link is claimed -- a dangling `\\link` is an R CMD check failure, not a cosmetic blemish.
"""

import pytest

from codegen.emit.doc_links import LinkResolver, render_link, render_spans
from codegen.ir.doc import DocLine, FordLink


def link(item, item_type="subroutine"):
    return FordLink(component="m", component_type="module", item=item, item_type=item_type)


class TestAProcedureWithABinding:
    resolver = LinkResolver(bindings={"build_kd_index": "build_kd_index"})

    def test_python_gets_a_sphinx_reference(self):
        assert (render_link(link("build_kd_index"), self.resolver, "python")
                == ":func:`tensor_omics.build_kd_index`")

    def test_r_gets_a_help_link(self):
        assert (render_link(link("build_kd_index"), self.resolver, "r")
                == r"\code{\link{build_kd_index}}")


class TestAnImplementation:
    """An implementation is unreachable, but what it names is not: send the reader to the wrapper."""

    resolver = LinkResolver(implementations={"crunch_impl": ["crunch"]})

    def test_it_resolves_to_the_wrapper(self):
        assert (render_link(link("crunch_impl"), self.resolver, "python")
                == ":func:`tensor_omics.crunch`")

    def test_a_mode_split_impl_names_every_variant(self):
        resolver = LinkResolver(
            implementations={"detect_patterns_impl": ["detect_dosage_effect", "detect_subfunctionalization"]}
        )
        rendered = render_link(link("detect_patterns_impl"), resolver, "r")
        assert rendered == r"\code{\link{detect_dosage_effect}}, \code{\link{detect_subfunctionalization}}"

    def test_the_expert_variant_is_not_the_one_offered(self):
        # the plain name is the entry point; the _expert one takes the work arrays
        assert LinkResolver.prefer_plain(["crunch", "crunch_expert"]) == ["crunch"]

    def test_but_it_is_offered_when_it_is_all_there_is(self):
        assert LinkResolver.prefer_plain(["crunch_expert"]) == ["crunch_expert"]


class TestAModeParameter:
    def test_in_a_split_family_it_names_the_procedure(self):
        # there is no mode argument left to pass a value to, so the parameter names a function
        resolver = LinkResolver(bindings={"detect_dosage_effect": "detect_dosage_effect"},
                                modes={"mode_dosage": "detect_dosage_effect"})
        assert (render_link(link("MODE_DOSAGE", "variable"), resolver, "python")
                == ":func:`tensor_omics.detect_dosage_effect`")

    def test_otherwise_it_names_the_value_the_caller_writes(self):
        resolver = LinkResolver(modes={"method_ward": "'ward'"})
        assert render_link(link("METHOD_WARD", "variable"), resolver, "python") == "``'ward'``"

    def test_a_value_is_never_linked(self):
        resolver = LinkResolver(modes={"method_ward": "'ward'"})
        assert "link" not in render_link(link("METHOD_WARD", "variable"), resolver, "r")


class TestNothingToLinkTo:
    @pytest.mark.parametrize("language, expected", [
        ("python", "``loess_degenerate_fit``"),
        ("r", r"\code{loess_degenerate_fit}"),
    ])
    def test_the_name_is_shown_as_code_and_no_link_is_claimed(self, language, expected):
        assert render_link(link("loess_degenerate_fit"), LinkResolver(), language) == expected

    def test_and_the_same_without_a_resolver_at_all(self):
        assert render_link(link("whatever"), None, "python") == "``whatever``"


class TestProseIsUntouched:
    def test_only_the_link_span_changes(self):
        line = DocLine.parse("use [[m(module):crunch_impl(subroutine)]] for this")
        resolver = LinkResolver(implementations={"crunch_impl": ["crunch"]})

        assert (render_spans(line, resolver, "python")
                == "use :func:`tensor_omics.crunch` for this")

    def test_a_line_with_no_link_is_returned_as_it_was(self):
        line = DocLine.parse("the number of genes in the study")

        assert render_spans(line, LinkResolver(), "python") == "the number of genes in the study"
