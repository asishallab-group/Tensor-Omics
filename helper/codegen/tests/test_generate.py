"""The pipeline and CLI: parse to written files.

Checks the orchestration -- that every target is produced, that a clean run writes and a
dirty one does not, that `--check` reports without writing. The emitters themselves are
tested elsewhere; here the question is only that they are wired together correctly.
"""

import re
from pathlib import Path

import pytest

from codegen import cli
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.frontend.ford_frontend import FordFrontend
from codegen.generate import generate, generate_and_write

from conftest import REPO_ROOT

FIXTURE_SRC = Path("helper/codegen/tests/fixtures/src")
#: the real tox_errors lives here, so a full generate can build the error module
# The whole source tree: a DM_DEFAULT may reference a parameter from another package
# (max_angle's default is f42_math_impl's PI), and constant resolution only sees the modules
# that are parsed. The generator is always run on the full tree, so these end-to-end tests
# are too. Parsing it per test is slow, so the `real_project` fixture parses it once and the
# tests reuse it through generate(parsed=...); only the CLI tests, which go through argv,
# parse for themselves.
REAL_SRC = Path("src")


def paths(root: Path, src: Path = FIXTURE_SRC) -> Paths:
    return Paths(root=root, src_dir=src)


@pytest.fixture(scope="module")
def real_project():
    """Parse the whole real source once, for the tests that generate from it.

    Ford is the one slow stage, and the rest of the pipeline is idempotent on the parsed
    project (roles are re-derived, not accumulated), so every generate() test can pass this
    in as `parsed=` instead of re-running Ford. Only the CLI tests, which drive argv, parse
    for themselves.
    """
    bag = DiagnosticBag()
    parsed = FordFrontend(paths(REPO_ROOT, REAL_SRC), bag, CONVENTIONS).parse()
    assert bag.errors == (), bag.render()
    return parsed


class TestGenerate:
    def test_the_fixture_run_lacks_tox_errors_and_says_so(self):
        # the fixtures have no tox_errors module, so the error module cannot be built
        result = generate(paths(REPO_ROOT))

        assert not result.ok
        assert any("tox_errors" in e.message for e in result.diagnostics.errors)

    def test_a_source_tree_with_tox_errors_generates_cleanly(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        assert result.ok, result.diagnostics.render()

    def test_every_target_is_produced(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        suffixes = {file.path.suffix for file in result.files}
        assert {".py", ".h", ".R"} <= suffixes

    def test_the_error_module_is_generated_for_each_language(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        names = {file.path.name for file in result.files}
        assert "error_handling.py" in names
        assert "error_handling.R" in names

    def test_a_single_target_produces_only_that_target(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("python",), parsed=real_project)

        suffixes = {file.path.suffix for file in result.files}
        assert suffixes == {".py"}

    def test_auto_output_from_is_generated_not_skipped(self, real_project):
        # every DM_OUTPUT_FROM(AUTO) in the tree resolves, so none is skipped with a warning
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("python",), parsed=real_project)

        assert not any("DM_OUTPUT_FROM" in w.message for w in result.diagnostics.warnings)

    def test_the_python_module_docstring_carries_the_whole_module_doc(self, real_project):
        """Not just its first line. A generated module is the published API, and taking the
        summary alone reduced every module in the package to one sentence."""
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("python",), parsed=real_project)
        written = {f.path.name: f.content for f in result.files}

        docstring = written["tox_loess.py"].split('"""')[1]
        impl = (REPO_ROOT / REAL_SRC / "tox" / "tox_loess_impl.F90").read_text()
        for line in impl.splitlines():
            if line.startswith("!| ") and len(line) > 20:
                assert line[3:] in docstring

    def test_the_python_module_docstring_says_do_not_edit_once(self, real_project):
        """The Fortran wrapper's own "generated from" note is added by its emitter, not
        carried in the IR -- carried, every language would print it beside their own."""
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("python",), parsed=real_project)
        written = {f.path.name: f.content for f in result.files}

        docstring = written["tox_loess.py"].split('"""')[1]
        assert docstring.lower().count("do not edit") == 1

    def test_files_are_not_written_by_generate(self, isolated_repo, real_project):
        result = generate(paths(isolated_repo, REAL_SRC), targets=("python",), parsed=real_project)

        # generate() builds content only; nothing reaches disk
        assert not (isolated_repo / "python").exists()
        assert result.files

    def test_nothing_is_written_when_there_are_errors(self, isolated_repo):
        # a fixture run errors on the missing tox_errors, so no file should appear
        result = generate_and_write(paths(isolated_repo, FIXTURE_SRC), targets=("python",))

        assert not result.ok
        assert not (isolated_repo / "python").exists()


class TestGenerateAndWrite:
    def test_it_writes_the_files(self, isolated_repo, real_project):
        result = generate_and_write(
            paths(isolated_repo, REAL_SRC), targets=("python",), parsed=real_project
        )

        assert result.ok
        assert (isolated_repo / "python" / "tensor_omics" / "error_handling.py").is_file()

    def test_clean_removes_a_stale_file(self, isolated_repo, real_project):
        out = isolated_repo / "python" / "tensor_omics"
        out.mkdir(parents=True)
        stale = out / "was_exported_once.py"
        stale.write_text("# a wrapper for a procedure no longer exported")

        generate_and_write(
            paths(isolated_repo, REAL_SRC), targets=("python",), clean=True, parsed=real_project
        )

        assert not stale.exists()

    def test_no_clean_keeps_a_stale_file(self, isolated_repo, real_project):
        out = isolated_repo / "python" / "tensor_omics"
        out.mkdir(parents=True)
        stale = out / "leftover.py"
        stale.write_text("# left alone")

        generate_and_write(
            paths(isolated_repo, REAL_SRC), targets=("python",), clean=False, parsed=real_project
        )

        assert stale.exists()


class TestCli:
    def test_check_writes_nothing_and_reports_changes(self, isolated_repo, capsys):
        code = cli.main(["--root", str(isolated_repo), "--src", str(REAL_SRC),
                         "--target", "python", "--check", "--no-color"])

        assert code == 1  # nothing generated yet, so files would change
        assert not (isolated_repo / "python").exists()
        assert "would change" in capsys.readouterr().err

    def test_check_is_up_to_date_after_a_write(self, isolated_repo):
        common = ["--root", str(isolated_repo), "--src", str(REAL_SRC), "--target", "python"]
        assert cli.main(common) == 0

        assert cli.main(common + ["--check"]) == 0

    def test_a_run_that_errors_exits_non_zero(self, isolated_repo, capsys):
        # the fixtures lack tox_errors
        code = cli.main(["--root", str(isolated_repo), "--src", str(FIXTURE_SRC),
                         "--no-color"])

        assert code == 1
        assert "nothing was written" in capsys.readouterr().err

    def test_a_clean_run_reports_what_it_wrote(self, isolated_repo, capsys):
        code = cli.main(["--root", str(isolated_repo), "--src", str(REAL_SRC),
                         "--target", "python", "--no-color"])

        assert code == 0
        assert "generated" in capsys.readouterr().err


@pytest.fixture
def isolated_repo(tmp_path):
    """A tmp dir that looks enough like the repo for a generate run.

    The frontend reads fpm.toml, src and macros.h from the root, so those are linked in;
    output goes to the tmp dir, leaving the real tree untouched.
    """
    for name in ("fpm.toml", "src", "helper"):
        (tmp_path / name).symlink_to(REPO_ROOT / name)
    return tmp_path


class TestAParseThatFailed:
    """A source the frontend could not read stops the run, rather than half-building an IR.

    An argument whose type has no mapping is left with no type at all. Carrying that into the
    semantic pass raised an AttributeError over the top of the diagnostic the author was
    about to read -- a traceback instead of a report, for an ordinary authoring mistake.
    """

    def source(self, tmp_path):
        src = tmp_path / "src"
        src.mkdir()
        (src / "fx_untyped.F90").write_text(
            "#include <src/macros.h>\n"
            "!> summary: a module\n"
            "module fx_untyped\n"
            "    use, intrinsic :: iso_fortran_env, only: int32\n"
            "    M_IMPLICIT_NONE\n"
            "contains\n"
            "    !> M_EXPORT_C\n"
            "    !| summary: p\n"
            "    !| author: A\n"
            "    subroutine fx_p(factor, ierr)\n"
            "        double precision, intent(in) :: factor\n"
            "            !! a type the generator has no C mapping for\n"
            "        integer(int32), intent(out) :: ierr\n"
            "            !! Error code\n"
            "    end subroutine fx_p\n"
            "end module fx_untyped\n"
        )
        return src

    def test_it_reports_instead_of_raising(self, tmp_path):
        result = generate(Paths(root=REPO_ROOT, src_dir=self.source(tmp_path)))

        assert not result.ok
        assert any("unsupported type" in d.message for d in result.diagnostics.errors)

    def test_nothing_is_generated_from_it(self, tmp_path):
        result = generate(Paths(root=REPO_ROOT, src_dir=self.source(tmp_path)))

        assert result.files == []


class TestDeterminism:
    """The same sources must produce the same bytes, every run.

    `--check` is what keeps the committed bindings honest, and it can only do that if
    generation is deterministic. Much of the pipeline walks `set`s -- the prologue and
    producer imports, the bound helpers -- and a set that reached the output unsorted would
    make `--check` flap rather than fail, which is worse than either.
    """

    def test_two_runs_of_the_whole_tree_agree_byte_for_byte(self, real_project):
        first = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)
        second = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        assert first.ok and second.ok
        assert {f.path: f.content for f in first.files} == {
            f.path: f.content for f in second.files
        }

    def test_the_file_list_is_ordered_the_same_way(self, real_project):
        first = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)
        second = generate(paths(REPO_ROOT, REAL_SRC), parsed=real_project)

        assert [f.path for f in first.files] == [f.path for f in second.files]


class TestTheAllocatingSignatureAndItsBodyAgree:
    """What the allocating wrapper does not ask for, it must declare -- and vice versa.

    The drop set is computed twice: in `synthesize`, which shapes the signature, and in the
    emitter, which writes the body. They are the same call, but they are separate call sites
    with separate inputs, and every time those two have disagreed it has produced a wrapper
    that either loses an argument or declares one twice. This asserts the invariant over the
    real tree rather than over one fixture.
    """

    def wrappers(self, real_project):
        from codegen.emit.fortran_wrapper import FortranWrapperEmitter, WrapperInfo
        from codegen.ir.roles import analyse_project
        from codegen.synthesize import synthesize_wrappers

        bag = DiagnosticBag()
        synthesis = synthesize_wrappers(real_project.project, CONVENTIONS, bag)
        analyse_project(synthesis.project, bag)
        info = {}
        for spec in synthesis.specs:
            for wrapper in (spec.validating, spec.allocating):
                if wrapper is not None:
                    info[wrapper.name.lower()] = WrapperInfo(spec.impl.name, spec.mode_fix)
        emitter = FortranWrapperEmitter(project=synthesis.project)
        for spec in synthesis.specs:
            if spec.allocating is None:
                continue
            module = synthesis.project.module(spec.module_name)
            emitter._wrapper_info = info
            text = emitter.subroutine(spec.allocating, module)
            yield spec, text

    def exposed(self, spec):
        """What this wrapper pair exposes of the implementation.

        Read off the validating wrapper rather than off the implementation, because a mode-split
        variant exposes less: the mode argument is fixed in the implementation call and dropped, and
        an argument scoped to another mode is absent entirely. `foo` is exactly that set,
        plus the `ierr` the wrapper may have synthesised.
        """
        return [
            argument
            for argument in spec.validating.arguments
            if spec.impl.argument(argument.name) is not None
        ]

    def test_every_exposed_argument_is_a_dummy_or_a_local_and_never_both(self, real_project):
        for spec, text in self.wrappers(real_project):
            signature = text[text.index("(") : text.index(")")]
            body = text[text.index(")") :]
            dummies = {a.name.lower() for a in spec.allocating.arguments}
            for argument in self.exposed(spec):
                if argument.name.lower() in dummies:
                    continue
                assert f":: {argument.name}" in body, (
                    f"{spec.allocating.name}: '{argument.name}' is neither a dummy nor a "
                    f"local -- the signature and the body disagree about the drop set"
                )
                assert argument.name not in signature, spec.allocating.name

    def test_every_exposed_argument_still_reaches_the_impl(self, real_project):
        for spec, text in self.wrappers(real_project):
            call = text[text.index(f"call {spec.impl.name}(") :]
            for argument in self.exposed(spec):
                assert f"{argument.name} = " in call, (
                    f"{spec.allocating.name}: '{argument.name}' is not passed to the implementation"
                )

    def test_every_allocating_dummy_is_accounted_for(self, real_project):
        # the other direction: nothing reaches that signature except the implementation's arguments,
        # the prologue's own, and ierr. Asserting only over the *exposed* set would miss a
        # prologue argument entirely, since it is on the allocating wrapper alone.
        from codegen.config import CONVENTIONS
        from codegen.synthesize import prologue_only_arguments

        for spec, _ in self.wrappers(real_project):
            allowed = {a.name.lower() for a in spec.impl.arguments}
            allowed |= {
                a.name.lower()
                for a in prologue_only_arguments(
                    spec.prologue, spec.impl.arguments, CONVENTIONS
                )
            }
            allowed.add(CONVENTIONS.error_arg.lower())
            for argument in spec.allocating.arguments:
                assert argument.name.lower() in allowed, (
                    f"{spec.allocating.name}: '{argument.name}' is neither the implementation's, "
                    f"the prologue's, nor ierr"
                )

    def test_the_invariant_covers_the_real_families(self, real_project):
        # a property test that silently matched nothing would assert nothing
        specs = [spec for spec, _ in self.wrappers(real_project)]
        assert len(specs) >= 10
        assert any(spec.mode_fix is not None for spec in specs), "no mode-split wrapper seen"


class TestAModeDefaultSpeaksEachLayersLanguage:
    """`DM_DEFAULT` on a mode argument quotes back the author's Fortran, and that is an
    integer. Every layer that types the mode as a string has to say the string instead --
    and the one layer that keeps the integer has to go on saying the integer.

    `fx_modes`' `link_method` is what exercises it: optional, carrying a mode table, and
    `DM_DEFAULT(METHOD_WARD)`. Nothing in the tree had that combination before, which is how
    a Python docstring came to read `default 'ward'` on one line and ``The default value is
    `METHOD_WARD`.`` on the next.
    """

    @pytest.fixture(scope="class")
    def written(self):
        # the fixture tree has no tox_errors, so the run is not ok; the files are still built
        result = generate(paths(REPO_ROOT), targets=("c", "python", "r"))
        return {file.path.name: file.content for file in result.files}

    def test_python_says_the_mode_string(self, written):
        assert "The default value is `'ward'`." in written["fx_basics.py"]
        assert "METHOD_WARD" not in written["fx_basics.py"]

    def test_r_quotes_it_the_way_r_does(self, written):
        # matching its own type line, which lists the modes as "ward", "single"
        assert 'The default value is `"ward"`.' in written["fx_basics.R"]
        assert "METHOD_WARD" not in written["fx_basics.R"]

    def test_the_c_wrapper_says_it_too(self, written):
        # its dummy is character(len=1, kind=c_char), so the parameter names a value no
        # caller of this layer can pass. Its mode table still cites the parameter by name.
        assert "The default value is `'ward'`." in written["fx_basics_c.F90"]
        assert "The default value is `METHOD_WARD`." not in written["fx_basics_c.F90"]

    def test_the_fortran_wrapper_keeps_the_integer(self, real_project):
        """The negative control. There the dummy really is an integer, so the author's own
        literal is the right thing to print and rewriting it would be the same bug pointing
        the other way."""
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("fortran",), parsed=real_project)
        written = {file.path.name: file.content for file in result.files}

        wrapper = written["tox_get_outliers.F90"]
        assert "integer(int32), intent(in), optional :: mode" in wrapper
        assert "The default value is `1_int32`." in wrapper
        assert "The default value is `'robust'`." not in wrapper


class TestTheCLayerAllocatesNothingOnTheStack:
    """A converted local sized by something C supplied is an automatic array unless the
    emitter is careful, and an automatic array is stack storage. At 10 million elements a
    logical mask is 40 MB: ifx segfaults on it against a default 8 MB stack, and gfortran
    escapes only by quietly rehousing large automatics. It is reachable from the published
    API -- `serialize_logical` takes `arr(n_elements)` of the caller's choosing.

    So the rule is a property of the whole emitted layer, not of one wrapper: every such
    local is `allocatable` and reaches the heap through `M_ALLOCATE`, whose failure path
    returns `ERR_ALLOC_FAIL` to the caller instead of crashing inside the wrapper.

    Since strings became pointer views the real project reaches something stronger still:
    no wrapper allocates anything at all, on either the stack or the heap. The two rules
    below are what keeps that from silently becoming an automatic array again if a new
    conversion ever needs a local.
    """

    #: `<type> ... dimension(<something that is not just ':'>) :: name` on a local, i.e. an
    #: automatic array. Dummies are excluded by the `intent(` they all carry.
    AUTOMATIC = re.compile(
        r"^\s+(?:logical|character|integer|real|complex)[^:\n]*dimension\([^:)]", re.M
    )

    @pytest.fixture(scope="class")
    def c_wrappers(self, real_project):
        result = generate(paths(REPO_ROOT, REAL_SRC), targets=("c",), parsed=real_project)
        return {f.path.name: f.content for f in result.files if f.path.suffix == ".F90"}

    def test_no_wrapper_allocates_at_all(self, c_wrappers):
        # It was 16 wrappers, then 3 when logicals moved to c_bool, then none when strings
        # became pointer views of the caller's buffer. Nothing the emitter can express now
        # needs storage of its own: this is the strongest form of the rule, and the two
        # tests below are what will catch a new conversion reaching for the stack instead.
        offenders = [name for name, text in c_wrappers.items() if "M_ALLOCATE" in text]
        assert offenders == []

    def test_no_wrapper_declares_a_runtime_sized_automatic_local(self, c_wrappers):
        offenders = []
        for name, text in c_wrappers.items():
            for line in text.splitlines():
                if "intent(" in line or "allocatable" in line:
                    continue
                if self.AUTOMATIC.match(line):
                    offenders.append(f"{name}: {line.strip()}")
        assert offenders == []

    def test_every_allocation_goes_through_the_checked_macro(self, c_wrappers):
        # a bare `allocate(` has no stat=, so a failure aborts instead of returning ierr
        bare = [f"{name}: {line.strip()}"
                for name, text in c_wrappers.items()
                for line in text.splitlines()
                if line.strip().startswith("allocate(")]
        assert bare == []
