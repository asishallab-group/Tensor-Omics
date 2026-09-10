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


class TestMarshalHeaderProtection:
    """R's rule: anything allocated must be protected before the next allocation.

    Not reachable from a runtime test. `gctorture(TRUE)` does collect an unprotected
    `out`, but the freed node is only reused by a later allocation of its own size class,
    and the fill loop allocates CHARSXPs -- so a 400-string round trip comes back intact
    either way. The invariant is therefore asserted on the emitted source, which is also
    the level a static checker like rchk works at.
    """

    def test_tox_char_out_protects_across_the_fill(self):
        header = CCallEmitter().marshal_header_content()
        body = header[header.index("tox_char_out(const char*"):]
        body = body[: body.index("\n}")]

        # Rf_mkCharLen allocates, so `out` must already be protected when the loop runs
        assert "PROTECT(Rf_allocVector(STRSXP, n))" in body
        assert "Rf_mkCharLen" in body
        assert body.index("PROTECT(") < body.index("Rf_mkCharLen")
        assert "UNPROTECT(1)" in body

    def test_the_helpers_that_stay_unprotected_never_allocate_twice(self):
        # tox_bool_out and tox_shape_of hand back an unprotected SEXP too, which is safe
        # only as long as nothing after their single allocation can trigger a GC
        header = CCallEmitter().marshal_header_content()
        for helper in ("tox_bool_out", "tox_shape_of"):
            body = header[header.index(f"{helper}("):]
            body = body[: body.index("\n}")]
            allocators = [call for call in ("Rf_allocVector", "Rf_mkChar", "Rf_mkCharLen",
                                           "Rf_ScalarInteger", "Rf_duplicate", "Rf_coerceVector")
                          if call in body]
            assert len(allocators) <= 1, f"{helper} allocates more than once: {allocators}"
