"""The VS Code snippets: one template per procedure, in three languages.

The snippets are data, not code that runs, so these assert on structure -- that each file
parses, that every prefix is namespaced, that a snippet's arguments match the wrapper the
other emitters build. The point is the same as `test_emitter_parity`: a snippet a
developer expands must not drift from the function that actually exists.
"""

import json
import re
from pathlib import Path

import pytest

from codegen.abi.c_abi import build_project
from codegen.config import Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.python_ctypes import PythonEmitter
from codegen.emit.vscode_snippets import SnippetEmitter
from codegen.frontend.ford_frontend import FordFrontend
from codegen.ir.errors import ErrorCatalogue
from codegen.ir.roles import analyse_project
from codegen.ir.validate import validate_project
from conftest import REPO_ROOT

#: `${1:name}` or `${1|a,b|}` -- a tabstop, with its default text or its choices
_TABSTOP = re.compile(r"\$\{(\d+)(?::([^}]*)|\|([^}]*)\|)?\}")


@pytest.fixture(scope="module")
def built():
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=Path("src")), bag).parse()
    analyse_project(parsed.project, bag)
    validate_project(parsed.project, bag)
    interface = build_project(parsed.project, bag)
    assert bag.errors == (), bag.render()
    catalogue = ErrorCatalogue.from_module(
        parsed.project.module("tox_errors"), bag, parsed.project.constant_values()
    )
    return interface, catalogue


@pytest.fixture(scope="module")
def files(built) -> dict:
    """Every generated file, `{file_name: parsed_json}`."""
    interface, catalogue = built
    return {name: json.loads(content)
            for name, content in SnippetEmitter().snippets_files(interface, catalogue).items()}


@pytest.fixture(scope="module")
def snippets(files) -> dict:
    """Every snippet, merged across files (the keys are globally unique)."""
    merged: dict = {}
    for bucket in files.values():
        merged.update(bucket)
    return merged


def _root(module) -> str:
    return module.stripped_name.partition("_")[0]


def _bucket(files: dict, language: str, root: str) -> dict:
    return files.get(f"{language}_{root}_snippets.json", {})


class TestStructure:
    def test_the_six_expected_files_are_produced(self, files):
        assert set(files) == {
            f"{language}_{root}_snippets.json"
            for language in ("Fortran", "Python", "R")
            for root in ("f42", "tox")
        }

    def test_every_snippet_is_well_formed(self, snippets):
        assert snippets
        for name, snippet in snippets.items():
            assert snippet["prefix"], name
            assert "scope" not in snippet, name  # the file name carries the language now
            assert isinstance(snippet["body"], list) and snippet["body"], name

    def test_a_file_only_holds_prefixes_of_its_own_root(self, files):
        # the whole reason for the split: no snippet collides, each file is one namespace
        for name, bucket in files.items():
            root = name.split("_")[1]
            for snippet in bucket.values():
                assert snippet["prefix"].startswith(f"{root}:"), (name, snippet["prefix"])

    def test_error_module_uses_the_tox_prefix(self, files):
        # tox_errors yields root `tox` with no special-casing
        assert _bucket(files, "Fortran", "tox")["tox_errors code"]["prefix"] == "tox:err"


class TestCalls:
    def test_python_and_r_calls_ask_for_the_same_arguments(self, built, files):
        """The parity the other suites guard, seen through the snippets."""
        interface, _ = built
        for module in interface:
            root = _root(module)
            for wrapper in module:
                key = f"{wrapper.stripped_name} (call)"
                python = _argument_names(_bucket(files, "Python", root)[key])
                r = _argument_names(_bucket(files, "R", root)[key])
                assert python == r, wrapper.stripped_name

    def test_a_call_asks_for_exactly_the_wrapper_inputs(self, built, files):
        interface, _ = built
        emitter = PythonEmitter()
        for module in interface:
            root = _root(module)
            for wrapper in module:
                expected = [a.name for a in emitter._inputs(wrapper)]
                snippet = _bucket(files, "Python", root)[f"{wrapper.stripped_name} (call)"]
                assert _argument_names(snippet) == expected, wrapper.stripped_name

    def test_a_mode_argument_becomes_a_choice(self, built, files):
        interface, _ = built
        seen = False
        for module in interface:
            root = _root(module)
            for wrapper in module:
                for argument in wrapper:
                    if argument.mode is None:
                        continue
                    seen = True
                    strings = ",".join(v.string for v in argument.mode.values)
                    snippet = _bucket(files, "Python", root)[f"{wrapper.stripped_name} (call)"]
                    assert f"|{strings}|" in "\n".join(snippet["body"]), wrapper.stripped_name
        assert seen, "no mode argument in src to exercise the choice path"

    def test_a_checked_variant_guards_ierr(self, built, files):
        interface, _ = built
        for module in interface:
            bucket = _bucket(files, "Fortran", _root(module))
            for wrapper in module:
                key = f"{wrapper.procedure.name} (call, checked)"
                has_ierr = wrapper.procedure.argument("ierr") is not None
                if not has_ierr:
                    assert key not in bucket
                    continue
                assert bucket[key]["body"][-1] == "if (is_err(ierr)) return"


class TestModuleAndGeneric:
    def test_each_module_has_use_and_import_snippets(self, built, files):
        interface, _ = built
        for module in interface:
            root = _root(module)
            assert f"{module.stripped_name} (use)" in _bucket(files, "Fortran", root)
            assert f"{module.stripped_name} (import)" in _bucket(files, "Python", root)

    def test_the_generic_aids_are_present(self, files):
        assert "tensor_omics (R setup)" in _bucket(files, "R", "tox")
        assert "tox_error (handler)" in _bucket(files, "Python", "tox")
        assert "tox_error (handler)" in _bucket(files, "R", "tox")
        assert "tox_errors code" in _bucket(files, "Fortran", "tox")

    def test_the_error_picker_is_dropped_without_a_catalogue(self, built):
        interface, _ = built
        produced = SnippetEmitter().snippets_files(interface, ErrorCatalogue(()))
        fortran_tox = json.loads(produced["Fortran_tox_snippets.json"])
        assert "tox_errors code" not in fortran_tox


def _argument_names(snippet: dict) -> list[str]:
    """The keyword-argument names a call snippet asks for, `name=...` before each tabstop."""
    names = []
    for line in snippet["body"]:
        match = re.match(r"\s*([A-Za-z_]\w*)\s*=", line)
        if match and _TABSTOP.search(line):
            names.append(match.group(1))
    return names
