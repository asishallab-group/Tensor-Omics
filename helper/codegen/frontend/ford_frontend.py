"""Ford's parse tree, converted into the IR.

The only module that imports Ford. Everything downstream sees `ir.Project` and cannot
tell where it came from, which is what lets the test suite build one by hand.

Ford is used as a Fortran parser and nothing more. Its own model is shaped for generating
documentation pages -- `meta.author` arrives as rendered HTML, `dimension` hides in the
attribute list, entities carry no source position -- so this module's job is to translate
all of that into terms the generator was designed around, once, here.
"""

from __future__ import annotations

import contextlib
import io
import os
import re
import warnings
from dataclasses import dataclass, replace
from pathlib import Path

from ..config import CONVENTIONS, Conventions, Paths
from ..diagnostics import DiagnosticBag, SourceLocation
from ..ir.directives import DirectiveError, DirectiveParser, Directives
from ..ir.doc import Doc, DocParseError
from ..ir.entities import Argument, Declaration, Meta, Module, Parameter, Procedure, Project
from ..ir.types import (
    BaseType,
    CharacterLength,
    Dimension,
    FortranType,
    Intent,
    UnsupportedTypeError,
)
from .macros import (
    MacroTable,
    build_directive_patterns,
    error_arg_pos_factor,
    export_category,
)
from .source_index import SourceIndex

_DIMENSION_ATTRIB_RE = re.compile(r"\bdimension\s*\(", re.IGNORECASE)
_LEN_RE = re.compile(r"\blen\s*=\s*(\*|:|[^,)]+)", re.IGNORECASE)
_BASE_TYPE_RE = re.compile(r"\s*([A-Za-z_]+)", re.IGNORECASE)
#: The name inside a derived type spec: `type(hashmap_type)` -> hashmap_type
_DERIVED_NAME_RE = re.compile(r"\s*type\s*\(\s*([A-Za-z]\w*)", re.IGNORECASE)
#: Ford renders meta tags as markdown, so `author` arrives wrapped in an anchor
_HTML_TAG_RE = re.compile(r"<[^>]+>")
#: A macro that survived preprocessing, i.e. one whose definition was never seen.
#:
#: Any of the three prefixes is a finding, and there are no false positives to worry
#: about: the preprocessor does not respect backticks or comments, so a macro that *is*
#: defined always expands, even where it was only meant to be named. A macro-shaped token
#: still standing after preprocessing is therefore undefined -- a typo, or a missing
#: include -- and never a deliberate mention.
_UNEXPANDED_MACRO_RE = re.compile(r"\b(?P<prefix>DM|CM|M)_[A-Z][A-Z_0-9]*")

#: What each prefix means, for the diagnostic
_MACRO_PREFIXES = {
    "M": "a macro from a macro header",
    "CM": "a custom macro defined in the source file itself",
    "DM": "a documentation macro",
}


@dataclass(frozen=True)
class ParsedProject:
    """What a frontend run produces."""

    project: Project
    macros: MacroTable
    #: How tox_errors packs argument positions, read from M_ERR_ARG_POS_FACTOR
    arg_pos_factor: int


class FordFrontend:
    def __init__(
        self,
        paths: Paths = Paths(),
        diagnostics: DiagnosticBag | None = None,
        conventions: Conventions = CONVENTIONS,
    ):
        self.paths = paths
        self.diagnostics = diagnostics if diagnostics is not None else DiagnosticBag()
        self.index = SourceIndex(paths.root)
        self.macros = MacroTable(
            paths.resolve(paths.macros_header), include_paths=(paths.root,)
        )
        self.directives = DirectiveParser(build_directive_patterns(self.macros))
        # the export category comes from M_EXPORT_C, so the marker and what the generator
        # recognises cannot drift; the passed conventions supply everything else
        self.conventions = replace(
            conventions, c_binding_category=export_category(self.macros)
        )

    def parse(self) -> ParsedProject:
        ford_project = self._parse_with_ford()
        modules = [self._module(module) for module in ford_project.modules]
        return ParsedProject(
            # `[extra.ford] project` in fpm.toml -- "TensorOmics". Ford has already read it,
            # so taking it from there keeps one spelling of the name in the repository. Left
            # unset, every diagnostic reported against a module ended with `in project ''`.
            project=Project(modules, name=ford_project.name),
            macros=self.macros,
            arg_pos_factor=error_arg_pos_factor(self.macros),
        )

    # -- Ford -------------------------------------------------------------------

    def _parse_with_ford(self):
        # Imported here so the rest of the package can be used, and tested, without Ford
        from ford.fortran_project import Project as FordProject
        from ford.settings import load_markdown_settings, load_toml_settings

        # Ford animates a rich progress bar from a background refresh thread that writes to
        # whatever `sys.stdout` currently is. Its reader redirects `sys.stdout` to a StringIO
        # around the preprocessor to capture the preprocessed source -- and the spinner
        # thread races into that same buffer, interleaving progress glyphs into the Fortran
        # the parser then reads. That intermittently mangles a declaration (`real(real64)`
        # loses its kind, "a kind is required"), with nothing wrong in the sources. Ford
        # disables the bar (and its thread) entirely when this is set, which is the cure.
        os.environ["FORD_DEBUGGING"] = "1"

        root = self.paths.root or Path()
        settings = load_toml_settings(root)
        if settings is None:
            ford_yml = root / "ford.yml"
            if not ford_yml.is_file():
                raise FileNotFoundError(
                    f"no Ford settings: neither {root / 'fpm.toml'} nor {ford_yml} exists"
                )
            settings, _ = load_markdown_settings(root, ford_yml.read_text(), str(ford_yml))

        # Absolute, so they do not depend on the working directory below
        settings.src_dir = [self.paths.resolve(self.paths.src_dir).resolve()]
        # The whole generated tree, by directory. Ford must not re-read the generator's own
        # output and define `tox_<name>` twice (once parsed, once synthesised), and the one
        # rule "everything generated lives under `generated_dir`" covers every target at
        # once -- including a stale wrapper left behind by a deleted implementation, which
        # a list derived from the implementations that still exist would not name.
        #
        # By directory rather than by file name, which is what this used to be: Ford matches
        # a bare exclude as `**/<name>`, so once an implementation under `src/f42` generates
        # `f42_stats.F90`, excluding that name would also drop the hand-written
        # `src/f42/utils/f42_stats.F90` from the parse -- and every binding generated from
        # it would vanish with no error at all.
        settings.exclude_dir = list(settings.exclude_dir) + [
            str(self.paths.resolve(self.paths.generated_dir).resolve())
        ]

        # Ford's default drops a procedure's own variables, which is right for published
        # documentation and wrong here: a local declared `allocatable` is how the generator
        # sees that a procedure allocates (`validate._check_impl_allocates`), and it never
        # reads a body. Nothing is published from this parse, so keeping them costs a list.
        settings.proc_internals = True

        # fpm.toml configures `pcpp -D__GFORTRAN__ -I.`, whose include path is relative
        # to the working directory. Run from anywhere else, `#include <src/macros.h>`
        # fails, every DM_ directive silently vanishes with it, and the wrappers come out
        # quietly wrong. The include path cannot be corrected from here either: Ford
        # splits the setting on whitespace (`settings.preprocessor.split()`), so an
        # absolute path containing a space -- as this repository's does -- would arrive
        # as several broken arguments. So honour the project's own `-I.` by giving it the
        # directory it is written against.
        with contextlib.chdir(root.resolve()):
            # Ford narrates its progress to stdout and warns freely; neither is this
            # generator's output, and both drown its diagnostics.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with contextlib.redirect_stdout(io.StringIO()):
                    project = FordProject(settings)
                    project.correlate()
        return project

    # -- entities ---------------------------------------------------------------

    def _module(self, ford_module) -> Module:
        path = self._path_of(ford_module)
        name = ford_module.name

        location = self.index.module(path, name)
        doc, _ = self._documentation(ford_module.doc_list, location)

        return Module(
            name=name,
            procedures=[
                self._procedure(routine, path)
                for routine in self._module_routines(ford_module)
            ],
            parameters=[
                self._parameter(variable, path)
                for variable in ford_module.variables
                if getattr(variable, "parameter", False)
            ],
            declarations=list(self._declarations(ford_module, path)),
            doc=doc,
            meta=self._meta(ford_module),
            location=location,
            uses=self._uses_of(ford_module),
        )

    def _declarations(self, ford_module, path: Path):
        """The public generic interfaces and derived types, by name and documentation.

        Nothing is generated from either, so nothing more of them is modelled. They are read
        because documentation names them -- `[[f42_math_impl(module):is_close(interface)]]` --
        and because their own doc comments hold links, which live nowhere else in the IR:
        an interface block's `!>` comment is attached to no procedure and no module.

        Only *generic* interfaces. An interface block also declares the module procedures of a
        submodule, which `_module_routines` already takes as procedures, and the `bind(C)`
        externals of libzip, which are nobody's public API.
        """
        for interface in getattr(ford_module, "interfaces", ()) or ():
            name = getattr(interface, "name", None)
            if not name or not getattr(interface, "generic", False):
                continue
            location = self.index.module(path, ford_module.name)
            doc, _ = self._documentation(getattr(interface, "doc_list", []) or [], location)
            yield Declaration(name=name, kind="interface", doc=doc, location=location)
        for derived in getattr(ford_module, "types", ()) or ():
            name = getattr(derived, "name", None)
            if not name:
                continue
            location = self.index.module(path, ford_module.name)
            doc, _ = self._documentation(getattr(derived, "doc_list", []) or [], location)
            yield Declaration(name=name, kind="type", doc=doc, location=location)

    @staticmethod
    def _uses_of(entity) -> list[str]:
        """The names of the modules `entity` -- a module or a procedure -- `use`s.

        Ford turns `uses` into a set once it has correlated the project, holding the module
        object where it resolved one and the bare name where it did not (the intrinsic
        modules, anything outside the source tree). Sorted, so a re-export module generated
        from these comes out the same on every run.

        A procedure carries its own list: Fortran lets a `use` sit inside a procedure as well
        as at module level, and a rule about what a module may reach has to see both.
        """
        return sorted(
            getattr(used, "name", used) for used in getattr(entity, "uses", ()) or ()
        )

    @staticmethod
    def _module_routines(ford_module):
        """Every procedure the module defines, including ones split into a submodule.

        A procedure whose body lives in a submodule is declared here as a `module
        subroutine`/`module function` inside an `interface` block, carrying its full
        signature and documentation. Ford lists those under `.interfaces`, not
        `.routines`, and reports the implementing submodule with no routines of its own --
        so the interface declaration is the only place the generator can read them from.

        Only interface bodies marked `module` are taken: an `interface` block also holds
        `bind(C)` externals (e.g. the libzip bindings) and abstract interfaces, which are
        not this module's procedures and carry types the wrapper has no business mapping.
        """
        yield from ford_module.routines
        for interface in getattr(ford_module, "interfaces", ()) or ():
            for routine in getattr(interface, "routines", ()) or ():
                if "module" in (getattr(routine, "attribs", ()) or ()):
                    yield routine

    def _procedure(self, ford_procedure, path: Path) -> Procedure:
        name = ford_procedure.name
        location = self.index.procedure(path, name)

        arguments = [
            self._argument(argument, path, name) for argument in ford_procedure.args
        ]

        result = None
        retvar = getattr(ford_procedure, "retvar", None)
        if retvar is not None:
            result = self._argument(retvar, path, name, is_result=True)
            # A function result is an output by definition; Fortran does not say so
            result.intent = Intent.OUT

        doc, directives = self._documentation(ford_procedure.doc_list, location)

        return Procedure(
            name=name,
            arguments=arguments,
            result=result,
            doc=doc,
            directives=directives,
            meta=self._meta(ford_procedure),
            location=location,
            conventions=self.conventions,
            allocatable_locals=self._allocatable_locals(ford_procedure),
            uses=self._uses_of(ford_procedure),
            is_pure=self._is_pure(ford_procedure),
        )

    @staticmethod
    def _is_pure(ford_procedure) -> bool:
        """Whether the procedure is `pure`. `elemental` implies it, and Ford lists both."""
        attributes = {
            attribute.strip().lower()
            for attribute in (getattr(ford_procedure, "attribs", ()) or ())
        }
        return bool(attributes & {"pure", "elemental"})

    @staticmethod
    def _allocatable_locals(ford_procedure) -> tuple[str, ...]:
        """The local variables declared `allocatable`.

        Ford splits a procedure's declarations into `args` and `variables` once the project
        is correlated, so what is left here is local. Only the `allocatable` attribute is
        collected: `pointer` locals are how an implementation aliases a buffer it was handed, which
        allocates nothing and stays allowed.
        """
        return tuple(
            variable.name
            for variable in getattr(ford_procedure, "variables", ()) or ()
            if "allocatable" in {
                attribute.strip().lower()
                for attribute in (getattr(variable, "attribs", ()) or ())
            }
        )

    def _argument(self, variable, path: Path, procedure: str, is_result: bool = False) -> Argument:
        name = variable.name
        location = self.index.argument(path, procedure, name)
        doc_line = self.index.argument_doc(path, procedure, name)

        doc, directives = self._documentation(variable.doc_list, location, doc_line)

        return Argument(
            name=name,
            type=self._type(variable, location),
            dimension=self._dimension(variable),
            intent=self._intent(variable),
            optional=bool(getattr(variable, "optional", False)),
            doc=doc,
            directives=directives,
            attributes=self._attributes(variable),
            location=location,
            is_result=is_result,
        )

    def _parameter(self, variable, path: Path) -> Parameter:
        name = variable.name
        location = self.index.variable(path, name)
        doc, _ = self._documentation(
            variable.doc_list, location, self.index.variable_doc(path, name)
        )
        return Parameter(
            name=name,
            type=self._type(variable, location, report=False),
            expression=(getattr(variable, "initial", None) or "").strip(),
            doc=doc,
            location=location,
        )

    # -- pieces -----------------------------------------------------------------

    def _type(self, variable, location: SourceLocation, report: bool = True):
        # Ford renders a derived type's name as a link to its documentation page, so
        # full_type arrives as `type(<a href='type/hashmap_type.html'>hashmap_type</a>)`
        full_type = _plain_text(getattr(variable, "full_type", None))
        match = _BASE_TYPE_RE.match(full_type)
        if match is None:
            if report:
                self.diagnostics.error(
                    f"cannot read the type of '{variable.name}' from '{full_type}'",
                    location=location,
                )
            return None

        try:
            base = BaseType.parse(match.group(1))
        except UnsupportedTypeError as error:
            if report:
                self.diagnostics.error(str(error), location=location)
            return None

        kind = getattr(variable, "kind", None) or None
        length = None
        if base is BaseType.CHARACTER:
            length_match = _LEN_RE.search(full_type)
            # `character :: x` is len=1 in Fortran, and Ford states no len for it
            length = CharacterLength.parse(
                length_match.group(1) if length_match else "1"
            )
            if kind is not None and kind.lower() != "c_char":
                kind = None

        derived_name = None
        if base is BaseType.DERIVED:
            # The name is inside the type spec, not in .kind
            name_match = _DERIVED_NAME_RE.match(full_type)
            derived_name = name_match.group(1) if name_match else (kind or "unknown")
            kind = None

        try:
            return FortranType(base=base, kind=kind, length=length, derived_name=derived_name)
        except ValueError as error:
            if report:
                self.diagnostics.error(
                    f"'{variable.name}' has an unusable type '{full_type}': {error}",
                    location=location,
                )
            return None

    @staticmethod
    def _dimension(variable) -> Dimension:
        """Read the extents.

        Ford puts `dimension(n)` in the attribute list and leaves `.dimension` empty when
        the attribute form is used, so the attribute has to win.
        """
        for attribute in getattr(variable, "attribs", ()) or ():
            extents = _dimension_attribute(attribute)
            if extents is not None:
                return Dimension.parse(extents)
        return Dimension.parse(getattr(variable, "dimension", "") or "")

    @staticmethod
    def _attributes(variable) -> tuple[str, ...]:
        return tuple(getattr(variable, "attribs", ()) or ())

    @staticmethod
    def _intent(variable) -> Intent | None:
        """The declared intent, or None when there is none.

        None rather than a default: `validate` reports a missing intent, and guessing
        `inout` here would silence it.
        """
        intent = (getattr(variable, "intent", "") or "").strip().lower()
        if not intent:
            return None
        try:
            return Intent(intent)
        except ValueError:
            return None

    def _meta(self, entity) -> Meta:
        meta = getattr(entity, "meta", None)
        return Meta(
            summary=_plain_text(getattr(meta, "summary", None)),
            author=_plain_text(getattr(meta, "author", None)),
            category=_plain_text(getattr(meta, "category", None)),
        )

    def _documentation(self, doc_list, location: SourceLocation,
                       first_line_number: int | None = None) -> tuple[Doc, Directives]:
        """Parse one doc comment and read its directives.

        Both come from the same lines, so they are produced together: parsing twice would
        mean reporting a malformed table twice as well.
        """
        self._check_macros_expanded(doc_list, location, first_line_number)

        try:
            doc = Doc.parse(_clean_doc(doc_list), first_line_number)
        except DocParseError as error:
            self._report(error, location)
            return Doc(), Directives()

        try:
            return doc, self.directives.parse(doc)
        except DirectiveError as error:
            self._report(error, location)
            return doc, Directives()

    def _report(self, error, location: SourceLocation) -> None:
        self.diagnostics.error(
            str(error),
            location=_at_line(location, getattr(error, "line_number", None)),
            note=getattr(error, "note", None),
        )

    def _check_macros_expanded(self, doc_list, location: SourceLocation,
                               first_line_number: int | None = None) -> None:
        """Report documentation that still contains an unexpanded macro.

        A macro-shaped token surviving preprocessing means its definition was never seen:
        a misspelt name, or a file that forgot `#include <src/macros.h>`.

        A `DM_` is an error, because it is silent: the directive is simply not found, the
        default or the output_from goes missing, and the generated wrapper comes out
        quietly wrong. The others are warnings -- a misspelt `M_`/`CM_` leaves the macro
        name standing in the rendered documentation, which is wrong but visible, and
        nothing is mis-generated from it.
        """
        for offset, line in enumerate(_clean_doc(doc_list)):
            for match in _UNEXPANDED_MACRO_RE.finditer(line):
                name = match.group(0)
                prefix = match.group("prefix")
                line_number = (
                    None if first_line_number is None else first_line_number + offset
                )
                report = (
                    self.diagnostics.error if prefix == "DM"
                    else self.diagnostics.warn
                )
                consequence = (
                    "the directive is silently ignored" if prefix == "DM"
                    else "the macro name is left standing in the documentation"
                )
                report(
                    f"'{name}' looks like {_MACRO_PREFIXES[prefix]} but never expanded, "
                    f"so it is not defined",
                    location=_at_line(location, line_number),
                    note=(
                        f"check the spelling, and that the file has "
                        f"'#include <{self.paths.macros_header}>'; otherwise "
                        f"{consequence}"
                    ),
                )

    def _path_of(self, entity) -> Path:
        source_file = getattr(entity, "source_file", None)
        path = getattr(source_file, "path", None)
        return Path(path) if path else Path("<unknown>")


def _dimension_attribute(attribute: str) -> str | None:
    """The extent list of a `dimension(...)` attribute, or None if it is not one.

    Bracket matching rather than a regex, because the extents are expressions and nest:
    `dimension(max(0_int32, n_timepoints - 1), n_factors)` and `dimension(clen,
    product(orig_shape))` both appear in the sources. A `\\(([^)]*)\\)` pattern stops at
    the first closing bracket and returns something unbalanced.
    """
    match = _DIMENSION_ATTRIB_RE.search(attribute)
    if match is None:
        return None

    depth = 0
    start = match.end() - 1
    for index in range(start, len(attribute)):
        character = attribute[index]
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            if depth == 0:
                return attribute[start : index + 1]
    return None


def _clean_doc(doc_list) -> list[str]:
    """Ford indents every documentation line by one space; drop it, keep the rest."""
    return [line[1:] if line.startswith(" ") else line for line in (doc_list or [])]


def _plain_text(value) -> str:
    """Strip the markup Ford renders meta tags into.

    `author` comes back as an anchor element, because Ford resolves the AUTHOR_* macro
    and then renders it as markdown for an HTML page. The generator wants the name.
    """
    if not value:
        return ""
    return _HTML_TAG_RE.sub("", str(value)).strip()


def _at_line(location: SourceLocation, line: int | None) -> SourceLocation:
    return location if line is None else SourceLocation(location.file, line)
