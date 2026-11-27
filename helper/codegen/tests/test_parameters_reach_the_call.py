"""Every parameter of every published Python and R binding reaches its C call unmodified.

`test_emitter_parity` checks the same property by reading the emitters' decisions; this checks
it by running the generated code. Each published function of every fixture library -- the
hand-written exports (`fixtures/src`), the implementations (`fixtures/impl`) and the carrier case
(`fixtures/carrier`) -- is called with a sentinel for
every parameter, its own C call is intercepted, and each parameter is compared with what
arrived. Producers the binding calls on the way run for real. An argument the binding computes
or derives over a parameter (the E2 attack's prototype had R overwrite the caller's `norms` with
the producer's), a parameter that never reaches the call, and a call that is never made all fail.

Python's call is intercepted by giving the module a `_lib` that answers the function's own C
symbol with a recorder; R's by a global `.Call` that records its own entry point and passes every
other through to `base::.Call`. Skipped without gfortran, gcc and R.
"""

import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from codegen.abi.c_abi import build_project
from codegen.abi import extents
from codegen.abi.model import Conversion
from codegen.config import CONVENTIONS, Paths
from codegen.diagnostics import DiagnosticBag
from codegen.emit.python_ctypes import PythonEmitter
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.emit.shape_checks import binding_axes
from codegen.frontend.ford_frontend import FordFrontend
from codegen.generate import _published_to_the_languages
from codegen.ir.roles import analyse_project
from codegen.ir.types import BaseType
from codegen.ir.validate import validate_project
from codegen.synthesize import synthesize_wrappers

import test_end_to_end as e2e_src_python
import test_end_to_end_r as e2e_src_r
from conftest import REPO_ROOT
from test_end_to_end_impl import (
    FIXTURE_CARRIER, GCC, GFORTRAN, RBIN, RSCRIPT, carrier_built, built,  # noqa: F401
)
from test_end_to_end_fortran import FIXTURE_IMPL

pytestmark = pytest.mark.skipif(
    GFORTRAN is None or GCC is None or RSCRIPT is None or RBIN is None,
    reason="needs gfortran, gcc and R",
)

#: The size every extent the caller's arrays carry is given, so arrays that share one agree
SIZE = 3


# -- the binding each library publishes ---------------------------------------------------------


def _published(src: Path, synthesise: bool):
    bag = DiagnosticBag()
    parsed = FordFrontend(Paths(root=REPO_ROOT, src_dir=src), bag).parse()
    if synthesise:
        synthesis = synthesize_wrappers(parsed.project, CONVENTIONS, bag)
        project = synthesis.project
    else:
        synthesis, project = None, parsed.project
    analyse_project(project, bag)
    validate_project(project, bag)
    binding = build_project(project, bag)
    assert bag.errors == (), bag.render()
    if synthesis is not None:
        binding = _published_to_the_languages(binding, synthesis)
    return binding


# -- the sentinels ------------------------------------------------------------------------------


def _scalar_extents(wrapper) -> dict[str, int]:
    """The integer scalars the caller passes that size an array, each given `SIZE`."""
    sized = {}
    for argument in wrapper:
        roles = argument.source.roles if argument.source is not None else None
        if argument.is_scalar and argument.type.base is BaseType.INTEGER and (
                (roles is not None and roles.is_extent) or argument.is_synthesised):
            sized[argument.name.lower()] = SIZE
    return sized


def _shape(argument, sized: dict[str, int]) -> list[int]:
    if argument.shape_arg is not None:
        return [SIZE, SIZE]
    shape = []
    for axis in binding_axes(argument):
        try:
            shape.append(max(extents.evaluate(extents.parse(axis), sized), 1))
        except Exception:  # an axis read off a size() of another array, or the like
            shape.append(SIZE)
    return shape or [1]


def _value(argument, index: int, sized: dict[str, int]) -> dict:
    """A sentinel for one parameter: its kind, shape and elements, in Fortran order."""
    if argument.conversion is Conversion.MODE:
        return {"kind": "string", "shape": [], "flat": [argument.mode.values[-1].string]}
    base = argument.type.base
    character = argument.type.is_character
    rank = argument.rank - (1 if character else 0)
    shape = _shape(argument, sized) if rank > 0 else []
    count = 1
    for extent in shape:
        count *= extent
    if character:
        length = argument.source.type.length if argument.source is not None else None
        width = 1 if (length is not None and str(length) == "1") else 2
        flat = [chr(ord("a") + (index + i) % 26) * width for i in range(count)]
        return {"kind": "string", "shape": shape, "flat": flat}
    if base is BaseType.LOGICAL:
        return {"kind": "logical", "shape": shape, "flat": [True] * count}
    if base is BaseType.INTEGER:
        if argument.is_scalar and argument.name.lower() in sized:
            return {"kind": "integer", "shape": [], "flat": [sized[argument.name.lower()]]}
        # small, so a size computed from it stays small, and never zero
        return {"kind": "integer", "shape": shape, "flat": [1 + (index + i) % 2 for i in range(count)]}
    if base is BaseType.COMPLEX:
        return {"kind": "complex", "shape": shape,
                "flat": [[0.25 + index + i, 0.5 + i] for i in range(count)]}
    return {"kind": "real", "shape": shape, "flat": [0.25 + 10 * index + i for i in range(count)]}


def _specs(binding) -> list[dict]:
    """What to call and with what, per published function, and where each parameter goes."""
    python, r = PythonEmitter(), RWrapperEmitter()
    specs = []
    for module in binding:
        for wrapper in module.with_constant_sizes_evaluated():
            sized = _scalar_extents(wrapper)
            positions = [argument.name for argument in wrapper]
            call_inputs = [argument.name for argument in r._call_inputs(wrapper)]
            parameters = []
            for index, argument in enumerate(python._inputs(wrapper)):
                parameters.append({
                    "name": argument.name,
                    "index": index,
                    "value": _value(argument, index, sized),
                    "c": positions.index(argument.name),
                    "r": call_inputs.index(argument.name) if argument.name in call_inputs else -1,
                })
            specs.append({"module": module.stripped_name, "function": wrapper.stripped_name,
                          "symbol": wrapper.name, "parameters": parameters})
    return specs


# -- Python -------------------------------------------------------------------------------------

_PYTHON = r'''
import ctypes, importlib, json, sys
import numpy as np

package, specs = sys.argv[1], json.loads(sys.stdin.read())

class Captured(Exception):
    pass

class Recorder:
    def __init__(self, library, symbol):
        self.library, self.symbol, self.args = library, symbol, None
    def __getattr__(self, name):
        if name != self.symbol:
            return getattr(self.library, name)
        def record(*args):
            self.args = args
            raise Captured()
        return record

def build(value):
    flat, shape = value["flat"], value["shape"]
    if value["kind"] == "complex":
        flat = [complex(re, im) for re, im in flat]
    if not shape:
        return flat[0]
    array = np.array(flat, dtype=object if value["kind"] == "string" else None)
    array = array.reshape(shape, order="F")
    return array.tolist() if value["kind"] == "string" else np.asfortranarray(array)

def arrived(arg):
    obj = getattr(arg, "_obj", arg)
    if isinstance(obj, ctypes.Array):
        return list(obj)
    if isinstance(obj, np.ndarray):
        return obj
    return getattr(obj, "value", obj)

def same(value, got):
    flat = value["flat"]
    if value["kind"] == "string":
        texts = [bytes(item).decode().rstrip(" \x00") for item in np.ravel(np.asarray(got))]
        return sorted(texts) == sorted(flat)
    if value["kind"] == "complex":
        if not value["shape"]:
            return list(got) == flat[0]
        got = np.ravel(np.asarray(got))
        return sorted((c.real, c.imag) for c in got) == sorted(tuple(c) for c in flat)
    got = np.ravel(np.asarray(got, dtype=float))
    return sorted(got.tolist()) == sorted(float(x) for x in flat)

def resized(parameter, axis, length):
    """The sentinel rebuilt with `length` along `axis`, as `_value` would have built it."""
    value = dict(parameter["value"]); shape = list(value["shape"]); shape[axis] = length
    count = 1
    for extent in shape:
        count *= extent
    k, index = value["kind"], parameter["index"]
    flat = {"real": [0.25 + 10 * index + i for i in range(count)],
            "integer": [1 + (index + i) % 2 for i in range(count)],
            "logical": [True] * count}[k]
    value.update(shape=shape, flat=flat)
    return dict(parameter, value=value)

import re
SHAPE = re.compile(r"'(\w+)' has \d+ along axis (\d+), but .*== (\d+);")

problems = {}
for spec in specs:
    module = importlib.import_module(f"{package}.{spec['module']}")
    recorder = Recorder(module._lib, spec["symbol"])
    module._lib = recorder
    try:
        function = getattr(module, spec["function"])
        for attempt in range(4):
            # an axis a producer sizes here can only be known by calling it: where the binding
            # says what it must be, the sentinel is rebuilt with that length
            try:
                function(**{p["name"]: build(p["value"]) for p in spec["parameters"]})
                break
            except Exception as error:
                match = SHAPE.search(str(error)) if type(error).__name__ == "ToxArgumentError" else None
                if match is None or attempt == 3:
                    raise
                name, axis, length = match.group(1), int(match.group(2)), int(match.group(3))
                spec["parameters"] = [resized(p, axis, length) if p["name"] == name else p
                                      for p in spec["parameters"]]
        problems[spec["function"]] = "returned without calling C"
    except Captured:
        wrong = [p["name"] for p in spec["parameters"]
                 if not same(p["value"], arrived(recorder.args[p["c"]]))]
        if wrong:
            problems[spec["function"]] = f"changed on the way: {wrong}"
    except Exception as error:
        problems[spec["function"]] = f"{type(error).__name__}: {error}"
    finally:
        module._lib = recorder.library
print(json.dumps(problems))
'''


def _run_python(path: Path, package: str, library: Path, specs) -> dict:
    env = dict(os.environ, TENSOR_OMICS_LIBRARY=str(library), PYTHONPATH=str(path))
    result = subprocess.run([sys.executable, "-c", _PYTHON, package], input=json.dumps(specs),
                            cwd=path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


# -- R ------------------------------------------------------------------------------------------


def _r_literal(value: dict) -> str:
    kind, flat, shape = value["kind"], value["flat"], value["shape"]
    if kind == "string":
        items = [json.dumps(s) for s in flat]
    elif kind == "logical":
        items = ["TRUE" if x else "FALSE" for x in flat]
    elif kind == "integer":
        items = [f"{x}L" for x in flat]
    elif kind == "complex":
        items = [f"complex(real = {re!r}, imaginary = {im!r})" for re, im in flat]
    else:
        items = [repr(x) for x in flat]
    vector = f"c({', '.join(items)})"
    if len(shape) >= 2:
        return f"array({vector}, dim = c({', '.join(map(str, shape))}))"
    return vector


def _r_check(value: dict, index: int) -> str:
    """An R expression that is TRUE where the `index`-th `.Call` argument holds the sentinel."""
    got, kind = f".captured[[{index + 1}]]", value["kind"]
    want = _r_literal(value)
    if kind == "string":
        return f"identical(sort(trimws({got}, 'right')), sort({want}))"
    if kind == "logical":
        return f"identical(sort(as.logical({got})), sort({want}))"
    if kind == "complex":
        return f"identical(sort(Re(as.vector({got}))), sort(Re({want})))"
    return f"isTRUE(all.equal(sort(as.numeric({got})), sort(as.numeric({want}))))"


def _r_body(specs) -> str:
    lines = [
        ".target <- ''; .captured <- NULL",
        ".Call <- function(.NAME, ...) {",
        "    if (identical(.NAME, .target)) {",
        "        .captured <<- list(...)",
        "        stop(structure(class = c('captured', 'condition'), list(message = 'captured', call = NULL)))",
        "    }",
        "    base::.Call(.NAME, ...)",
        "}",
        ".problems <- list()",
    ]
    for spec in specs:
        name = spec["function"]
        arguments = ", ".join(f"{p['name']} = {_r_literal(p['value'])}" for p in spec["parameters"])
        missing = [p["name"] for p in spec["parameters"] if p["r"] < 0]
        checks = " ".join(
            f"if (!({_r_check(p['value'], p['r'])})) .wrong <- c(.wrong, '{p['name']}');"
            for p in spec["parameters"] if p["r"] >= 0)
        lines.append(textwrap.dedent(f"""
            .target <- "{name}_call"; .captured <- NULL
            .outcome <- tryCatch({{ {name}({arguments}); "returned without calling C" }},
                captured = function(c) "captured", error = function(e) conditionMessage(e))
            if (identical(.outcome, "captured")) {{
                .wrong <- c({', '.join(repr(m) for m in missing)})
                {checks}
                if (length(.wrong)) .problems[["{name}"]] <- paste("changed on the way:", paste(.wrong, collapse = " "))
            }} else .problems[["{name}"]] <- .outcome
        """))
    lines.append("for (.n in names(.problems)) cat(.n, ':', .problems[[.n]], '\\n')")
    lines.append("cat('checked', " + str(len(specs)) + ", '\\n')")
    return "\n".join(lines)


# -- the libraries ------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def src_python(tmp_path_factory):
    out = tmp_path_factory.mktemp("sentinel_src_python")
    e2e_src_python._generate(out)
    return out, e2e_src_python._build(out)


@pytest.fixture(scope="module")
def src_r(tmp_path_factory):
    out = tmp_path_factory.mktemp("sentinel_src_r")
    e2e_src_r._generate(out)
    library = e2e_src_r._build_lib(out, with_r=True)
    (out / "libfixtures.so").write_bytes(library.read_bytes())
    return out


def _run_r_file(path: Path, header: str, body: str) -> str:
    """Run `header` + `body` as an R script file in `path` (too long for `Rscript -e`)."""
    script = path / "sentinels.R"
    script.write_text(textwrap.dedent(header) + body)
    result = subprocess.run([RSCRIPT, str(script)], cwd=path, capture_output=True, text=True)
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    return result.stdout


#: What R says where an axis a producer sizes was given another length, as `_r_body` prints it
_R_SHAPE = re.compile(r"^(\w+) : '(\w+)' has \d+ along dimension (\d+), but .*== (\d+);", re.M)


def _resized(parameter: dict, axis: int, length: int) -> dict:
    """The sentinel rebuilt with `length` along `axis`, as `_value` would have built it."""
    value = dict(parameter["value"])
    shape = list(value["shape"])
    shape[axis] = length
    count = 1
    for extent in shape:
        count *= extent
    kind, index = value["kind"], parameter["index"]
    flat = {"real": [0.25 + 10 * index + i for i in range(count)],
            "integer": [1 + (index + i) % 2 for i in range(count)],
            "logical": [True] * count}[kind]
    value.update(shape=shape, flat=flat)
    return dict(parameter, value=value)


def _run_r_resizing(path: Path, header: str, specs) -> str:
    """`_run_r_file` over `specs`, rebuilding a sentinel whose axis a producer sizes and rerunning.

    An axis the binding computes through a producer can only be known by calling it; where R
    says what it must be, the sentinel is rebuilt with that length, as the Python run does.
    """
    for _ in range(4):
        out = _run_r_file(path, header, _r_body(specs))
        found = _R_SHAPE.findall(out)
        if not found:
            return out
        for function, name, dimension, length in found:
            for spec in specs:
                if spec["function"] == function:
                    spec["parameters"] = [
                        _resized(p, int(dimension) - 1, int(length)) if p["name"] == name else p
                        for p in spec["parameters"]]
    return out


class TestEveryParameterReachesTheCall:

    def test_the_hand_written_exports_in_python(self, src_python):
        out, library = src_python
        specs = _specs(_published(e2e_src_python.FIXTURE_SRC, synthesise=False))

        assert len(specs) > 10, "the premise: the fixtures publish something"
        assert _run_python(out, "tensor_omics", library, specs) == {}

    def test_the_hand_written_exports_in_r(self, src_r):
        specs = _specs(_published(e2e_src_r.FIXTURE_SRC, synthesise=False))

        out = _run_r_file(src_r, """
            dyn.load("libfixtures.so")
            source("error_handling.R"); source("tox_validate.R")
            source("fx_basics.R"); source("fx_edges.R")
        """, _r_body(specs))

        assert out.strip() == f"checked {len(specs)}", out

    def test_the_implementations_in_python(self, built):
        specs = _specs(_published(FIXTURE_IMPL, synthesise=True))

        assert "fx_norm_rank_expert" in {s["function"] for s in specs}
        assert _run_python(built / "py", "fx_impl", built / "libtensor_omics.so", specs) == {}

    def test_the_implementations_in_r(self, built):
        specs = _specs(_published(FIXTURE_IMPL, synthesise=True))

        modules = " ".join(f'source("r/{p.name}");' for p in sorted((built / "r").glob("fx_*.R")))
        out = _run_r_file(built, f"""
            dyn.load("libtensor_omics.so")
            source("r/error_handling.R"); source("r/tox_validate.R")
            {modules}
        """, _r_body(specs))

        assert out.strip() == f"checked {len(specs)}", out

    def test_the_carrier_case_in_python(self, carrier_built):
        specs = _specs(_published(FIXTURE_CARRIER, synthesise=True))

        assert {s["function"] for s in specs} >= {"fx_capacity_total", "fx_capacity_total_expert"}
        assert _run_python(carrier_built / "py", "fx_impl", carrier_built / "libtensor_omics.so",
                           specs) == {}

    def test_the_carrier_case_in_r(self, carrier_built):
        specs = _specs(_published(FIXTURE_CARRIER, synthesise=True))

        modules = " ".join(f'source("r/{p.name}");'
                           for p in sorted((carrier_built / "r").glob("fx_*.R")))
        out = _run_r_resizing(carrier_built, f"""
            dyn.load("libtensor_omics.so")
            source("r/error_handling.R"); source("r/tox_validate.R")
            {modules}
        """, specs)

        assert out.strip() == f"checked {len(specs)}", out
