"""The one grammar Python and the R shim compute extents in (#221 B11).

They allocate outputs and work arrays, and pad strings, before the call, so whatever an extent
means in Fortran it has to mean in Python and in C. The parser accepts exactly what does, and
`check_binding_extents` refuses the rest.
"""

import itertools

import pytest

from codegen.abi import extents
from codegen.abi.c_abi import build_module
from codegen.abi.model import CBinding
from codegen.diagnostics import DiagnosticBag
from codegen.ir.roles import analyse
from codegen.ir.types import Intent

import builders as b


def c_size(array, dimension):
    return f"SIZE({array}, {dimension})"


class TestTheGrammar:
    @pytest.mark.parametrize("text", [
        "n", "n_points - 1", "n_samples*n_timepoints", "max(0_int32, n_timepoints - 1)",
        "min(a, b, c)", "size(gene_to_fam)", "size(x, 2)", "(n_reps_S1 + n_reps_S2)*n",
        "-n + 2", "n/2", "+n", "2_c_int",
    ])
    def test_what_means_the_same_everywhere_is_accepted(self, text):
        extents.parse(text)

    @pytest.mark.parametrize("text, reason", [
        ("mod(n, 2)", "a call to 'mod'"),
        ("nint(x)", "a call to 'nint'"),
        ("sum(reps)", "a call to 'sum'"),
        ("n**2", "the operator '**'"),
        ("2.0*n", "the real literal '2.0'"),
        ("n*1.5_real64", "the real literal '1.5'"),
        ("a * -b", "'-' where an operand belongs"),
        ("max(n)", "'max' of a single operand"),
        ("size(x(1))", "'size' of something other than an argument"),
        ("size(x, k)", "a 'size' dimension that is not a positive literal"),
        ("n .and. m", "'.'"),
        ("n +", "an expression that ends early"),
        ("(n", "an expression that ends early"),
        ("n m", "'m' after a complete expression"),
    ])
    def test_anything_else_is_refused_saying_what(self, text, reason):
        with pytest.raises(extents.ExtentError, match=reason.replace("*", r"\*").replace("(", r"\(")):
            extents.parse(text)

    def test_a_sign_binds_below_a_product_as_in_fortran(self):
        assert extents.evaluate(extents.parse("-2*3 + 1")) == -5
        assert extents.parse("-a*b") == extents.Negation(
            extents.Binary("*", extents.Name("a"), extents.Name("b")))

    def test_the_names_read_include_the_arrays_asked_for_their_size(self):
        assert extents.names(extents.parse("max(n, size(x)) - m")) == {"n", "x", "m"}


class TestDivisionTruncatesTowardZero:
    """Fortran's integer `/` truncates; Python's `//` floors. `-7/2` is -3, not -4."""

    @pytest.mark.parametrize("a, b", list(itertools.product(range(-9, 10), (-4, -3, -2, -1, 1, 2, 3, 4))))
    def test_python_and_the_evaluator_agree_with_fortran(self, a, b):
        expected = int(a / b)   # exact here: small operands
        assert extents.evaluate(extents.parse("a/b"), {"a": a, "b": b}) == expected
        assert eval(extents.to_python("a/b")) == expected

    def test_python_is_exact_beyond_a_float(self):
        a, b = -(2**60 + 1), 3
        assert eval(extents.to_python("a/b")) == -((2**60 + 1) // 3)

    def test_c_keeps_its_own_division_which_truncates(self):
        assert extents.to_c("(n - 10)/3 + 5", str, c_size) == "(n - 10) / 3 + 5"


class TestWritingItOut:
    @pytest.mark.parametrize("text, python, c", [
        ("max(0_int32, n_timepoints - 1)", "max(0, n_timepoints - 1)",
         "tox_imax(0, n_timepoints - 1)"),
        ("n_samples*n_timepoints", "n_samples * n_timepoints", "n_samples * n_timepoints"),
        ("size(gene_to_fam)", "gene_to_fam.size", "SIZE(gene_to_fam, None)"),
        ("size(x, 2)", "x.shape[1]", "SIZE(x, 2)"),
        ("min(a, b, c)", "min(a, b, c)", "tox_imin(a, tox_imin(b, c))"),
        ("n - (a - b)", "n - (a - b)", "n - (a - b)"),
        ("(a + b)*c", "(a + b) * c", "(a + b) * c"),
        ("a*(b/c)", "a * (b // c if (b >= 0) == (c > 0) else -(-b // c))", "a * (b / c)"),
    ])
    def test_each_language_gets_the_same_meaning(self, text, python, c):
        assert extents.to_python(text) == python
        assert extents.to_c(text, str, c_size) == c

    def test_the_names_are_the_callers_to_write(self):
        assert extents.to_c("n + 1", lambda name: f"_v_{name}", c_size) == "_v_n + 1"


def _published(*arguments, parameters=()):
    module = b.module("m", b.procedure("p", *arguments, b.ierr()), parameters=parameters)
    b.project(module)
    bag = DiagnosticBag()
    for procedure in module.exported_procedures:
        analyse(procedure, bag)
    binding = CBinding((build_module(module, bag),))
    assert bag.errors == (), bag.render()
    extents.check_binding_extents(binding, bag)
    return bag


class TestWhatTheBindingsComputeIsChecked:
    def test_an_output_in_the_grammar_passes(self):
        bag = _published(b.integer("n"), b.real("x", Intent.OUT, "((n - 10)/3 + 5)"),
                         b.real("y", Intent.OUT, "(N_AXES, max(1, n))"),
                         parameters=(b.parameter("N_AXES", "2"),))
        assert bag.errors == ()

    def test_an_output_outside_it_is_refused_pointing_at_a_producer(self):
        bag = _published(b.integer("n"), b.real("x", Intent.OUT, "(mod(n, 2))"))

        [error] = bag.errors
        assert "'x' of 'p' is sized by 'mod(n, 2)', which Python and R cannot compute" in (
            error.message)
        assert "a call to 'mod'" in error.message
        assert "DM_OUTPUT_FROM" in error.note

    def test_a_work_array_is_held_to_it_too(self):
        bag = _published(b.integer("n"), b.real("tmp_x", Intent.INOUT, "(n**2)"))
        assert "the operator '**'" in bag.errors[0].message

    def test_a_string_length_is_held_to_it_on_an_input(self):
        bag = _published(b.integer("n"), b.character("label", Intent.IN, length="nint(1.5)"))
        assert "'label' of 'p' is sized by 'nint(1.5)'" in bag.errors[0].message

    def test_an_input_array_is_read_never_computed(self):
        # its extents are read off what the caller passes (or passed to C, WP5's A-shape)
        bag = _published(b.integer("n"), b.real("x", Intent.IN, "(n)"))
        assert bag.errors == ()

    def test_an_output_read_by_an_extent_is_refused(self):
        bag = _published(b.integer("n"), b.real("y", Intent.OUT, "(n)"),
                         b.real("x", Intent.OUT, "(size(y))"))
        assert "it reads 'y', which is an output" in bag.errors[0].message

    def test_an_optional_read_by_an_extent_is_refused(self):
        bag = _published(b.integer("n"), b.real("w", Intent.IN, "(n)", optional=True),
                         b.real("x", Intent.OUT, "(size(w))"), b.real("v", Intent.IN, "(n)"))
        assert "it reads the optional 'w', which may be absent" in bag.errors[0].message

    def test_a_name_that_is_neither_argument_nor_constant_is_refused(self):
        bag = _published(b.integer("n"), b.real("x", Intent.OUT, "(n + n_cached)"))
        assert "'n_cached' is neither an argument nor a module constant" in bag.errors[0].message
