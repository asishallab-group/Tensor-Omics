"""Resolving a Ford cross-reference for the language that will read it.

An author writes `[[tox_get_outliers_kernel(module):compute_family_scaling_kernel(subroutine)]]`
and Ford turns it into a link. Carried into a Python docstring or an R help page it is markup
nothing renders, pointing at a symbol -- a kernel, a private helper -- that no caller of those
languages can reach.

So resolve rather than strip. What the reader wants is the thing they *can* call: the binding
generated from that kernel. Where there is no such thing, the name is rendered as plain code
and no link is claimed -- an `\\link` to a topic that does not exist is an R CMD check failure,
so the fallback carries weight.

Only the parsed link is touched. Ford links reach here as `FordLink` objects, never as text,
which is what keeps a text pass from mangling R's own `x[[i]]` indexing in generated code.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..ir.doc import FordLink

#: The binding-layer suffix for the variant that takes the work arrays. Where a kernel has
#: both, the plain name is the one to send a reader to.
_EXPERT_SUFFIX = "_expert"


@dataclass
class LinkResolver:
    """What each Ford link target is called in one target language, if anything."""

    #: Fortran procedure name (lower case) -> the binding's name in this language
    bindings: dict[str, str] = field(default_factory=dict)
    #: kernel procedure name (lower case) -> the bindings generated from it, non-expert first
    kernels: dict[str, list[str]] = field(default_factory=dict)
    #: mode parameter (lower case) -> the string a caller passes instead
    modes: dict[str, str] = field(default_factory=dict)

    def targets(self, link: FordLink) -> list[str]:
        """Every name this link resolves to in this language; empty if it resolves to none."""
        item = (link.item or link.component or "").lower()
        if not item:
            return []
        if item in self.bindings:
            return [self.bindings[item]]
        if item in self.kernels:
            return list(self.kernels[item])
        if item in self.modes:
            named = self.modes[item]
            # a per-mode procedure resolves on through to its own binding; a mode string is
            # already what the caller writes
            return [self.bindings.get(named.lower(), named)]
        return []

    @staticmethod
    def prefer_plain(names: list[str]) -> list[str]:
        """Drop the `_expert` variant where the plain one is also present."""
        plain = [n for n in names if not n.endswith(_EXPERT_SUFFIX)]
        return plain or names


def build(binding, specs) -> tuple[dict[str, str], dict[str, list[str]], dict[str, str]]:
    """The three lookups, from the C binding set and the synthesis specs.

    `specs` carries kernel -> wrapper, which no name convention can recover for a mode-split
    family: `detect_patterns_kernel` becomes `detect_dosage_effect` and
    `detect_subfunctionalization`, names that come from the mode table.
    """
    bindings: dict[str, str] = {}
    modes: dict[str, str] = {}
    for module in binding:
        for wrapper in module:
            bindings.setdefault(wrapper.procedure.name.lower(), wrapper.stripped_name)
            for argument in wrapper:
                if argument.mode is None:
                    continue
                for value in argument.mode.values:
                    # In a split family the parameter names no value a caller ever passes --
                    # there is no mode argument left. What it names is a *procedure*, which
                    # the mode table's third column already gives, so link to that. Elsewhere
                    # the mode really is a value, and the string is what the caller writes.
                    if value.procedure_name:
                        modes.setdefault(value.parameter.lower(), value.procedure_name)
                    else:
                        modes.setdefault(value.parameter.lower(), f"'{value.string}'")

    kernels: dict[str, list[str]] = {}
    for spec in specs or ():
        made = [w.name for w in (spec.validating, spec.allocating) if w is not None]
        named = [bindings[n.lower()] for n in made if n.lower() in bindings]
        if named:
            kernels.setdefault(spec.kernel.name.lower(), []).extend(named)
    for name, found in kernels.items():
        # a reader wants the entry point, not the variant that takes the work arrays
        kernels[name] = LinkResolver.prefer_plain(sorted(set(found)))
    return bindings, kernels, modes


def _one(name: str, language: str, linkable: bool) -> str:
    """One resolved (or unresolved) target, in the language's own markup."""
    if language == "r":
        # \link only where the topic is known to exist: a dangling one fails R CMD check
        return f"\\code{{\\link{{{name}}}}}" if linkable else f"\\code{{{name}}}"
    # Sphinx resolves :func: against the package; everything is re-exported from its root
    return f":func:`tensor_omics.{name}`" if linkable else f"``{name}``"


def render_link(link: FordLink, resolver: LinkResolver | None, language: str) -> str:
    """A Ford link as the target language should show it."""
    targets = resolver.targets(link) if resolver is not None else []
    if targets:
        # a mode parameter resolves to the string a caller passes, which is a value and not a
        # topic to link to
        linkable = not all(t.startswith("'") for t in targets)
        return ", ".join(_one(t, language, linkable) for t in targets)
    fallback = link.item or link.component
    return _one(fallback, language, linkable=False)


def render_spans(block, resolver: LinkResolver | None, language: str) -> str:
    """A doc line with its links resolved and its prose untouched."""
    spans = getattr(block, "spans", None)
    if not spans:
        return block.text
    return "".join(
        render_link(span, resolver, language) if isinstance(span, FordLink) else str(span)
        for span in spans
    )
