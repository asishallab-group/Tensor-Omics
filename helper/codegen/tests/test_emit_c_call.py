"""The C shim R calls into.

Focused checks on how the shim materialises its arguments. Whole-module output is covered
by the end-to-end suites; what is here is the marshalling that a compiler, not a reader,
would otherwise be the first to catch.
"""

import pytest

from codegen.abi.c_abi import build_wrapper
from codegen.diagnostics import DiagnosticBag
from codegen.emit.c_call import CCallEmitter
from codegen.ir.directives import Default, Directives
from codegen.ir.roles import analyse
from codegen.ir.types import Intent

import builders as b


@pytest.fixture
def bag():
    return DiagnosticBag()


def emitted(procedure, bag):
    analyse(procedure, bag)
    wrapper = build_wrapper(procedure, bag)
    assert bag.errors == (), bag.render()
    return CCallEmitter().function(wrapper)


class TestOptionalScalarThatNeedsConverting:
    """A logical is four bytes to R and one to Fortran, so it is converted, not pointed at.

    An optional one has to be both: converted into a local, and passed as a pointer that is
    null when the caller omits it. Getting only the pointer half emits `<name>_p` at the call
    with nothing declaring it, which no Python or Fortran test can see -- it is a C compile
    error in the generated shim.
    """

    def wrapper(self, bag):
        return emitted(
            b.procedure(
                "p",
                b.real("values", Intent.IN, "(n)", doc="the data"),
                b.integer("n", Intent.IN, doc="length"),
                b.logical("refine", Intent.IN, optional=True, doc="whether to refine"),
                b.real("result", Intent.OUT, doc="the answer"),
                b.ierr(),
            ),
            bag,
        )

    def test_it_declares_the_local_it_converts_into(self, bag):
        assert "unsigned char refine_v = 0;" in self.wrapper(bag)

    def test_it_declares_the_pointer_it_passes(self, bag):
        assert "const unsigned char* refine_p = NULL;" in self.wrapper(bag)

    def test_it_converts_and_points_only_when_the_argument_is_present(self, bag):
        text = self.wrapper(bag)
        assert "refine_v = (Rf_asLogical(refine) == TRUE) ? 1 : 0;" in text
        assert "refine_p = &refine_v;" in text
        # the conversion sits inside the presence test, not before it
        assert text.index("if (refine != R_NilValue)") < text.index("refine_p = &refine_v;")

    def test_the_call_passes_the_pointer(self, bag):
        assert "refine_p," in self.wrapper(bag)


class TestScalarThatNeedsConvertingButIsMandatory:
    def test_a_mandatory_logical_is_still_passed_by_address(self, bag):
        text = emitted(
            b.procedure(
                "p",
                b.logical("refine", Intent.IN, doc="whether to refine"),
                b.real("result", Intent.OUT, doc="the answer"),
                b.ierr(),
            ),
            bag,
        )

        assert "unsigned char refine_v = (Rf_asLogical(refine) == TRUE) ? 1 : 0;" in text
        assert "&refine_v," in text


class TestDefaultedOptionalIsNotNullable:
    def test_an_optional_with_a_default_is_passed_by_address(self, bag):
        # a default is supplied by the binding, so C always receives a value -- no pointer
        text = emitted(
            b.procedure(
                "p",
                b.logical(
                    "refine",
                    Intent.IN,
                    optional=True,
                    doc="whether to refine",
                    directives=Directives(default=Default(".false.")),
                ),
                b.real("result", Intent.OUT, doc="the answer"),
                b.ierr(),
            ),
            bag,
        )

        assert "refine_p" not in text
        assert "&refine_v," in text
