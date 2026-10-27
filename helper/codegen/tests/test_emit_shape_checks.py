"""Which input axes the bindings compare, and that Python and R compare the same ones (#221).

The rule (`emit.shape_checks`): every input axis that names an extent is compared with that
extent's value, except the axis the value was read from, and an axis the declaration fixes at a
constant with the constant. The end-to-end tests run it; these pin what is compared.
"""

import pytest

from codegen.abi.c_abi import build_module
from codegen.diagnostics import DiagnosticBag
from codegen.emit.python_ctypes import PythonEmitter
from codegen.emit.r_wrapper import RWrapperEmitter
from codegen.emit.shape_checks import shape_checks
from codegen.ir.roles import analyse
from codegen.ir.types import Intent

import builders as b


def wrapper_of(*arguments, parameters=()):
    module = b.module("m", b.procedure("p", *arguments, b.ierr()), parameters=parameters)
    b.project(module)
    bag = DiagnosticBag()
    analyse(module.procedure("p"), bag)
    built = build_module(module, bag)
    assert bag.errors == (), bag.render()
    return built.with_constant_sizes_evaluated().wrappers[0]


def compared(wrapper):
    """(array, axis, origin, extent or value) for every comparison, in order."""
    return [
        (check.position.array.name, check.position.axis, check.expected.origin,
         check.expected.value if check.expected.origin == "constant" else check.expected.extent)
        for check in shape_checks(wrapper)
    ]


class TestTheRule:
    def test_the_axis_an_extent_is_read_off_is_not_compared_with_itself(self):
        wrapper = wrapper_of(b.integer("n"), b.real("x", Intent.IN, "(n)"),
                             b.real("y", Intent.IN, "(n)"))

        assert compared(wrapper) == [("y", 0, "axis", "n")]

    def test_every_axis_naming_it_is_compared_a_square_one_too(self):
        wrapper = wrapper_of(b.integer("n"), b.real("d", Intent.IN, "(n, n)"))

        assert compared(wrapper) == [("d", 1, "axis", "n")]

    def test_a_lone_array_has_nothing_to_compare(self):
        assert compared(wrapper_of(b.integer("n"), b.real("x", Intent.IN, "(n)"))) == []

    def test_an_optional_array_is_compared_and_never_read_from(self):
        # main skipped optional owners altogether: build_residual_histograms' 1x1 mask over-read
        wrapper = wrapper_of(b.integer("n"), b.logical("mask", Intent.IN, "(n)", kind="c_bool",
                                                       optional=True),
                             b.real("x", Intent.IN, "(n)"))

        assert compared(wrapper) == [("mask", 0, "axis", "n")]
        assert shape_checks(wrapper)[0].expected.position.array.name == "x"

    def test_an_input_sized_by_a_mask_count_is_compared_with_the_count(self):
        wrapper = wrapper_of(
            b.integer("n_genes"), b.logical("genes_mask", Intent.IN, "(n_genes)", kind="c_bool"),
            b.integer("n_selected_genes"), b.real("weights", Intent.IN, "(n_selected_genes)"),
        )

        assert compared(wrapper) == [("weights", 0, "mask", "n_selected_genes")]

    def test_a_constant_axis_is_compared_with_its_value(self):
        wrapper = wrapper_of(b.integer("n"), b.real("fields", Intent.IN, "(n, 2, N_PAIR)"),
                             parameters=(b.parameter("N_PAIR", "3"),))

        assert compared(wrapper) == [("fields", 1, "constant", 2), ("fields", 2, "constant", 3)]

    def test_an_axis_c_is_given_is_the_c_wrapper_s_to_compare(self):
        wrapper = wrapper_of(b.integer("n"), b.real("x", Intent.IN, "(n + 1)"),
                             b.real("y", Intent.IN, "(n)"))

        assert compared(wrapper) == []

    def test_a_string_array_compares_its_count_and_not_its_length(self):
        wrapper = wrapper_of(b.integer("n"), b.real("x", Intent.IN, "(n)"),
                             b.character("names", Intent.IN, "(n)", length="*"))

        assert compared(wrapper) == [("names", 0, "axis", "n")]

    def test_outputs_and_work_arrays_are_allocated_not_compared(self):
        wrapper = wrapper_of(b.integer("n"), b.real("x", Intent.IN, "(n)"),
                             b.real("y", Intent.OUT, "(n)"), b.real("tmp_w", Intent.INOUT, "(n)"))

        assert compared(wrapper) == []


class TestBothBindingsWriteIt:
    @pytest.fixture
    def wrapper(self):
        return wrapper_of(
            b.integer("n"), b.real("x", Intent.IN, "(n, 2)"),
            b.real("w", Intent.IN, "(n)", optional=True),
        )

    def test_python(self, wrapper):
        text = PythonEmitter().function(wrapper)

        assert "if x.shape[1] != 2:" in text
        assert "if w is not None and w.shape[0] != n:" in text
        assert "but its declared shape fixes it at {2}; 'x' is declared with shape (n, 2)" in text
        assert "but 'x' implies n == {n}; 'w' is declared with shape (n)" in text

    def test_r(self, wrapper):
        text = RWrapperEmitter().function(wrapper)

        assert "if (dim(x)[2] != 2L)" in text
        assert "if (!is.null(w) && length(w) != dim(x)[1])" in text
        assert ('.tox_shape_error("x", 2L, dim(x)[2], "its declared shape fixes it at ", 2L, '
                '"n, 2")') in text
