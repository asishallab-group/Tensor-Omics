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


class TestStringWidthNeverReachesFortranAsZero:
    """`R_alloc(0, 1)` returns NULL, so a width of 0 hands the Fortran wrapper a null buffer
    it cannot make a pointer view of -- gfortran aborts, ifx segfaults. R reaches that with
    nothing exotic: `tox_max_strlen("")` is 0. The Python layer already floors the width
    (`.ljust(1)`, `max(..., default=0) or 1`); this is the R equivalent."""

    def widths(self):
        header = CCallEmitter().marshal_header_content()
        body = header[header.index("tox_max_strlen(SEXP x)"):]
        return body[: body.index("\n}")]

    def test_a_present_argument_reports_at_least_one(self):
        assert "return longest > 0 ? longest : 1;" in self.widths()

    def test_an_absent_argument_still_reports_zero(self):
        """The floor must NOT cover R_NilValue. An omitted optional arrives as a null buffer,
        and width 0 is how the wrapper is told so -- raising it to 1 would put width 1 beside
        a null pointer and make M_CHECK_CHARACTER_VIEW reject an optional the caller left
        out on purpose."""
        body = self.widths()

        early = body[: body.index("int longest")]
        assert "if (x == R_NilValue || TYPEOF(x) != STRSXP) return 0;" in early


class TestHeaderStamp:
    """fpm recompiles a source only when that file's own text changes -- it never hashes
    what the file includes (fortran-lang/fpm#358) -- so a change to `tox_marshal.h` alone
    used to rebuild nothing and leave the old helpers linked in. Every shim that includes the
    header therefore carries its hash."""

    def test_the_stamp_is_the_hash_of_the_header_as_emitted(self):
        import hashlib

        emitter = CCallEmitter()
        digest = hashlib.sha256(emitter.marshal_header_content().encode()).hexdigest()[:16]

        assert digest in emitter.header_stamp()
        assert emitter.header_stamp().startswith("// tox_marshal.h ")

    def test_any_change_to_the_header_changes_the_stamp(self):
        class Edited(CCallEmitter):
            def marshal_header_content(self):
                return super().marshal_header_content() + "\n// one more comment\n"

        assert Edited().header_stamp() != CCallEmitter().header_stamp()

    def test_every_shim_on_disk_carries_the_hash_of_the_header_on_disk(self):
        """The invariant fpm needs, checked on the committed tree: if a shim's stamp and the
        header disagree, the header changed without the shim changing with it."""
        import hashlib

        from conftest import REPO_ROOT

        binding = REPO_ROOT / "src/generated/bindings/r"
        header = (binding / "tox_marshal.h").read_text()
        digest = hashlib.sha256(header.encode()).hexdigest()[:16]
        includers = [c for c in sorted(binding.glob("*.c"))
                     if '#include "tox_marshal.h"' in c.read_text()]

        assert includers, "no shim includes the header -- the premise is gone"
        stale = [c.name for c in includers if digest not in c.read_text()]
        assert not stale, f"stamp does not match tox_marshal.h: {stale}"
