from dataclasses import replace

import pytest

from codegen.diagnostics import DiagnosticBag
from codegen.ir.entities import Meta
from codegen.ir.roles import analyse
from codegen.ir.types import Intent
from codegen.ir.validate import validate_module, validate_procedure, validate_project

import builders as b


@pytest.fixture
def bag():
    return DiagnosticBag()


def checked(procedure, bag):
    """Analyse then validate, as a real run does."""
    analyse(procedure, bag)
    validate_procedure(procedure, bag)
    return bag


def messages(bag):
    return [d.message for d in bag.errors]


def only_error(bag):
    assert len(bag.errors) == 1, messages(bag)
    return bag.errors[0]


class TestScope:
    def test_only_exported_procedures_are_validated(self, bag):
        # an internal routine may use anything Fortran allows
        module = b.module(
            "m",
            b.procedure(
                "internal",
                b.character("s", Intent.IN, length=":"),
                meta=Meta(summary="s", author="a"),
            ),
        )

        validate_module(module, bag)

        assert bag.errors == ()

    def test_an_exported_procedure_is_validated(self, bag):
        module = b.module("m", b.procedure("exported", b.character("s", Intent.IN, length=":")))

        validate_module(module, bag)

        assert "deferred length" in only_error(bag).message

    def test_validate_project_walks_every_module(self, bag):
        project = b.project(
            b.module("a", b.procedure("p", b.character("s", Intent.IN, length=":"))),
            b.module("b", b.procedure("q", b.character("t", Intent.IN, length=":"))),
        )

        validate_project(project, bag)

        assert len(bag.errors) == 2


class TestCharacters:
    def test_a_deferred_length_is_rejected(self, bag):
        checked(b.procedure("p", b.character("s", Intent.IN, length=":")), bag)

        error = only_error(bag)
        assert "deferred length (len=:)" in error.message
        assert "len=*" in error.note

    @pytest.mark.parametrize("length", ["*", "n", "42"])
    def test_other_lengths_are_accepted(self, bag, length):
        checked(b.procedure("p", b.character("s", Intent.IN, length=length), b.integer("n")), bag)

        assert bag.errors == ()

    def test_a_character_scalar_is_accepted(self, bag):
        checked(b.procedure("p", b.character("s", Intent.IN)), bag)

        assert bag.errors == ()

    def test_a_character_vector_is_accepted(self, bag):
        checked(b.procedure("p", b.character("s", Intent.IN, "(n)"), b.integer("n")), bag)

        assert bag.errors == ()

    def test_a_character_matrix_is_rejected(self, bag):
        # it would need a rank-3 c_char conversion, which tox_conversions has not got
        checked(
            b.procedure("p", b.character("s", Intent.IN, "(n, m)"), b.integer("n"), b.integer("m")),
            bag,
        )

        assert "rank-2 character array" in only_error(bag).message


class TestTypes:
    def test_a_derived_type_is_rejected(self, bag):
        from codegen.ir.entities import Argument
        from codegen.ir.types import BaseType, FortranType

        argument = Argument("p", FortranType(BaseType.DERIVED, derived_name="point"), intent=Intent.IN)
        checked(b.procedure("p", argument), bag)

        assert "derived type" in only_error(bag).message

    @pytest.mark.parametrize("attribute", ["allocatable", "pointer"])
    def test_non_interoperable_attributes_are_rejected(self, bag, attribute):
        checked(
            b.procedure("p", b.real("v", Intent.OUT, "(:)", attributes=(attribute,))), bag
        )

        assert f"is '{attribute}'" in only_error(bag).message

    def test_target_is_fine(self, bag):
        checked(b.procedure("p", b.real("v", Intent.OUT, "(n)", attributes=("target",)), b.integer("n")), bag)

        assert bag.errors == ()


class TestIntent:
    def test_a_missing_intent_is_rejected(self, bag):
        argument = b.real("v", Intent.IN)
        argument.intent = None

        checked(b.procedure("p", argument), bag)

        error = only_error(bag)
        assert "has no intent" in error.message
        assert "constness" in error.note

    def test_a_result_needs_no_declared_intent(self, bag):
        result = b.real("out", Intent.OUT, is_result=True)
        result.intent = None

        checked(b.procedure("f", result=result), bag)

        assert bag.errors == ()


class TestTemporaries:
    def test_an_intent_in_temporary_is_rejected(self, bag):
        checked(b.procedure("p", b.real("tmp_work", Intent.IN, "(n)"), b.integer("n")), bag)

        error = only_error(bag)
        assert "temporary 'tmp_work' is intent(in)" in error.message
        assert "tmp_" in error.note

    @pytest.mark.parametrize("intent", [Intent.OUT, Intent.INOUT])
    def test_out_and_inout_temporaries_are_accepted(self, bag, intent):
        checked(b.procedure("p", b.real("tmp_work", intent, "(n)"), b.integer("n")), bag)

        assert bag.errors == ()


class TestModeArguments:
    def _mode_doc(self):
        return [
            "| Mode | Value |",
            "|------|-------|",
            "| x | [[m(module):MODE_X(variable)]] |",
        ]

    def test_an_integer_scalar_mode_is_accepted(self, bag):
        checked(b.procedure("p", b.integer("mode", Intent.IN, doc=self._mode_doc())), bag)

        assert bag.errors == ()

    def test_a_non_integer_mode_is_rejected(self, bag):
        # the Fortran mode is an integer; only the C wrapper receives a string
        checked(b.procedure("p", b.character("mode", Intent.IN, doc=self._mode_doc())), bag)

        assert "is not a scalar integer" in only_error(bag).message

    def test_an_array_mode_is_rejected(self, bag):
        checked(
            b.procedure("p", b.integer("mode", Intent.IN, "(n)", doc=self._mode_doc()), b.integer("n")),
            bag,
        )

        assert "is not a scalar integer" in only_error(bag).message

    def _split_mode_doc(self):
        return [
            "| Mode | Value | Procedure |",
            "|------|-------|-----------|",
            "| x | [[m(module):MODE_X(variable)]] | p_x |",
        ]

    def _scoped_optional(self, mode_doc):
        from codegen.ir.directives import Default, Directives, RequiredIfMode

        return b.procedure(
            "p",
            b.integer("mode", Intent.IN, doc=mode_doc),
            b.real(
                "tuning",
                Intent.IN,
                optional=True,
                doc="a tuning knob",
                directives=Directives(
                    default=Default("1.0_real64"),
                    required_if_mode=RequiredIfMode("mode", "m", "MODE_X"),
                ),
            ),
        )

    def test_a_default_with_a_required_if_mode_is_rejected_for_a_runtime_mode(self, bag):
        # the argument is always passed on, so "required in that mode" says nothing
        checked(self._scoped_optional(self._mode_doc()), bag)

        assert "both a default and a mode it is required in" in only_error(bag).message

    def test_a_default_with_a_required_if_mode_is_accepted_when_the_mode_splits(self, bag):
        # there the directive scopes the argument to its mode, and the default applies within it
        checked(self._scoped_optional(self._split_mode_doc()), bag)

        assert bag.errors == ()


class TestShapeArguments:
    def _shape_procedure(self, shape_argument, data=None):
        return b.procedure("p", data or b.real("data", Intent.IN, "(:)"), shape_argument)

    def test_a_well_formed_shape_argument_is_accepted(self, bag):
        checked(self._shape_procedure(b.integer("data_shape", Intent.IN, "(:)")), bag)

        assert bag.errors == ()

    def test_a_shape_argument_must_be_intent_in(self, bag):
        checked(self._shape_procedure(b.integer("data_shape", Intent.OUT, "(:)")), bag)

        assert "is not intent(in)" in only_error(bag).message

    def test_a_shape_argument_must_be_a_rank_1_integer(self, bag):
        checked(self._shape_procedure(b.integer("data_shape", Intent.IN)), bag)

        assert "is not a rank-1 integer array" in only_error(bag).message

    def test_a_real_shape_argument_is_rejected(self, bag):
        checked(self._shape_procedure(b.real("data_shape", Intent.IN, "(:)")), bag)

        assert "is not a rank-1 integer array" in only_error(bag).message

    def test_a_shape_argument_must_not_be_optional(self, bag):
        # the c_loc ordering depends on it being there
        checked(self._shape_procedure(b.integer("data_shape", Intent.IN, "(:)", optional=True)), bag)

        error = only_error(bag)
        assert "shape argument 'data_shape' is optional" in error.message
        assert "c_loc" in error.note

    def test_the_described_argument_must_be_flat(self, bag):
        checked(
            self._shape_procedure(
                b.integer("data_shape", Intent.IN, "(:)"),
                data=b.real("data", Intent.IN, "(:, :)"),
            ),
            bag,
        )

        assert "has a shape argument but is rank 2" in only_error(bag).message


class TestExtentArguments:
    def test_an_extent_must_not_be_optional(self, bag):
        checked(
            b.procedure("p", b.integer("n", Intent.IN, optional=True), b.real("v", Intent.IN, "(n)")),
            bag,
        )

        error = only_error(bag)
        assert "extent argument 'n' is optional, but 'v' needs it" in error.message
        assert "c_loc" in error.note

    def test_the_error_names_every_dependent_array(self, bag):
        checked(
            b.procedure(
                "p",
                b.integer("n", Intent.IN, optional=True),
                b.real("a", Intent.IN, "(n)"),
                b.real("c", Intent.IN, "(n)"),
            ),
            bag,
        )

        assert "'a', 'c'" in only_error(bag).message

    def test_a_non_extent_scalar_may_be_optional(self, bag):
        checked(b.procedure("p", b.integer("seed", Intent.IN, optional=True)), bag)

        assert bag.errors == ()


class TestErrorArgument:
    def test_a_well_formed_ierr_is_accepted(self, bag):
        checked(b.procedure("p", b.ierr()), bag)

        assert bag.errors == ()

    def test_an_intent_in_ierr_is_rejected(self, bag):
        checked(b.procedure("p", b.integer("ierr", Intent.IN, doc="Error code")), bag)

        assert "no error can be reported through it" in only_error(bag).message

    def test_a_non_integer_ierr_is_rejected(self, bag):
        checked(b.procedure("p", b.real("ierr", Intent.OUT, doc="Error code")), bag)

        assert "is not a scalar integer" in only_error(bag).message

    def test_an_array_ierr_is_rejected(self, bag):
        checked(b.procedure("p", b.integer("ierr", Intent.OUT, "(2)", doc="Error code")), bag)

        assert "is not a scalar integer" in only_error(bag).message


class TestWarnings:
    def test_a_procedure_without_a_summary_warns(self, bag):
        checked(b.procedure("p", meta=Meta(author="a", category="C-binding")), bag)

        assert bag.errors == ()
        assert "no summary meta tag" in bag.warnings[0].message

    def test_a_procedure_without_an_author_warns(self, bag):
        checked(b.procedure("p", meta=Meta(summary="s", category="C-binding")), bag)

        assert "no author meta tag" in bag.warnings[0].message

    def test_an_undocumented_argument_warns(self, bag):
        checked(b.procedure("p", b.integer("n", Intent.IN)), bag)

        assert bag.errors == ()
        assert "argument 'n' has no documentation" in bag.warnings[0].message

    def test_an_undocumented_result_does_not_warn(self, bag):
        checked(b.procedure("f", result=b.real("out", Intent.OUT, is_result=True)), bag)

        assert bag.warnings == () or all("result" not in w.message for w in bag.warnings)

    def test_an_undocumented_module_warns(self, bag):
        validate_module(b.module("m"), bag)

        assert "module has no documentation" in bag.warnings[0].message

    def test_warnings_never_stop_generation(self, bag):
        checked(b.procedure("p", b.integer("n", Intent.IN), meta=Meta(category="C-binding")), bag)

        bag.raise_if_errors()  # must not raise


class TestEverythingTogether:
    def test_a_clean_procedure_produces_nothing(self, bag):
        procedure = b.procedure(
            "normalize_unit_length",
            b.integer("n_dims", Intent.IN, doc="number of elements in `vector`"),
            b.real("vector", Intent.INOUT, "(n_dims)", doc="Vector to normalise"),
            b.ierr(),
        )

        checked(procedure, bag)

        assert bag.errors == ()
        assert bag.warnings == ()

    def test_every_problem_is_reported_in_one_run(self, bag):
        procedure = b.procedure(
            "p",
            b.character("s", Intent.IN, length=":", doc="a string"),
            b.real("tmp_work", Intent.IN, "(n)", doc="work"),
            b.integer("n", Intent.IN, optional=True, doc="extent"),
        )

        checked(procedure, bag)

        assert len(bag.errors) == 3
        assert any("deferred length" in m for m in messages(bag))
        assert any("is intent(in)" in m for m in messages(bag))
        assert any("is optional" in m for m in messages(bag))


class TestOptionalOutputs:
    """Refused: no binding can honour one, and every binding ignores it silently, so the
    Fortran declaration and the generated wrapper disagree about the same argument."""

    def test_an_optional_output_is_rejected(self, bag):
        checked(
            b.procedure(
                "p",
                b.real("values", Intent.IN, "(n)"),
                b.real("diagnostics", Intent.OUT, "(n)", optional=True),
                b.integer("n"),
            ),
            bag,
        )

        error = only_error(bag)
        assert "'diagnostics' is an optional output" in error.message
        assert "compute_influence" in error.note

    def test_an_optional_work_array_is_accepted(self, bag):
        # dropped from the allocating wrapper, so it never reaches a caller
        checked(
            b.procedure(
                "p",
                b.real("tmp_work", Intent.OUT, "(n)", optional=True),
                b.integer("n"),
            ),
            bag,
        )

        assert bag.errors == ()

    def test_an_optional_input_is_accepted(self, bag):
        checked(b.procedure("p", b.logical("compute_influence", Intent.IN, optional=True)), bag)

        assert bag.errors == ()

    def test_a_mandatory_output_is_accepted(self, bag):
        checked(
            b.procedure("p", b.real("results", Intent.OUT, "(n)"), b.integer("n")), bag
        )

        assert bag.errors == ()


class TestKernelModules:
    """The rules that hold for a kernel, which is read rather than wrapped and so is not
    reached by the exported-procedure rules above."""

    def kernel_module(self, *procedures):
        return b.module("tox_thing_kernel", *procedures, doc="a kernel module")

    def test_a_kernel_that_allocates_is_rejected(self, bag):
        module = self.kernel_module(
            b.procedure(
                "thing_kernel",
                b.real("x", Intent.IN, "(n)"),
                b.integer("n"),
                meta=Meta(summary="s", author="a"),
                allocatable_locals=("workspace",),
            )
        )

        validate_module(module, bag)

        error = only_error(bag)
        assert "'thing_kernel' allocates: 'workspace'" == error.message
        assert "tmp_" in error.note

    def test_a_helper_that_allocates_is_rejected_too(self, bag):
        # a kernel that allocates nothing but calls a helper that does is no better off
        module = self.kernel_module(
            b.procedure(
                "thing_helper",
                b.real("x", Intent.IN, "(n)"),
                b.integer("n"),
                meta=Meta(summary="s", author="a"),
                allocatable_locals=("scratch",),
            )
        )

        validate_module(module, bag)

        assert "'thing_helper' allocates" in only_error(bag).message

    def test_a_pointer_local_is_accepted(self, bag):
        # aliasing a buffer the caller handed in allocates nothing
        module = self.kernel_module(
            b.procedure(
                "thing_kernel",
                b.real("x", Intent.IN, "(n)"),
                b.integer("n"),
                meta=Meta(summary="s", author="a"),
            )
        )

        validate_module(module, bag)

        assert bag.errors == ()

    def test_an_exported_kernel_is_rejected(self, bag):
        module = self.kernel_module(
            b.procedure("thing_kernel", b.integer("n"))  # b.procedure exports by default
        )

        validate_module(module, bag)

        errors = messages(bag)
        assert "kernel 'thing_kernel' is exported" in errors
        assert "M_EXPORT_C" in bag.errors[0].note

    def test_an_exported_support_routine_is_accepted(self, bag):
        # a recommend routine has no wrapper, so DM_OUTPUT_FROM needs it exported
        module = self.kernel_module(b.procedure("thing_required_workspace", b.integer("n")))

        validate_module(module, bag)

        assert bag.errors == ()

    def test_a_kernel_named_for_the_alloc_wrapper_is_rejected(self, bag):
        module = self.kernel_module(
            b.procedure(
                "thing_alloc_kernel",
                b.integer("n"),
                meta=Meta(summary="s", author="a"),
            )
        )

        validate_module(module, bag)

        error = only_error(bag)
        assert "kernel 'thing_alloc_kernel' is named for the allocating wrapper" == error.message
        assert "tmp_" in error.note

    def test_an_ordinary_module_is_not_held_to_these_rules(self, bag):
        module = b.module(
            "tox_thing",
            b.procedure(
                "thing_alloc",
                b.integer("n"),
                allocatable_locals=("workspace",),
            ),
            doc="an ordinary module",
        )

        validate_module(module, bag)

        assert bag.errors == ()


class TestPrologue:
    """A DM_PROLOGUE must name something the wrapper can actually call.

    Every one of these used to pass silently: the emitter was the directive's only consumer,
    and an emitter has no author's line to point at.
    """

    def checked(self, guard_arguments=None, kernel_arguments=None, bag=None):
        from test_synthesize import prologue_kernel_module

        module = prologue_kernel_module(guard_arguments)
        if kernel_arguments is not None:
            crunch = module.procedure("crunch_kernel")
            crunch.arguments = tuple(kernel_arguments)
        project = b.project(module)
        validate_project(project, bag if bag is not None else self.bag)
        return self.bag

    @pytest.fixture(autouse=True)
    def _bag(self, bag):
        self.bag = bag

    def guard(self, *arguments):
        return tuple(arguments)

    def test_a_well_formed_prologue_is_accepted(self):
        self.checked()

        assert self.bag.errors == ()

    def test_a_prologue_that_does_not_exist_is_rejected(self):
        from codegen.ir.directives import Prologue
        from test_synthesize import prologue_kernel_module

        module = prologue_kernel_module()
        module.procedure("crunch_kernel").directives = replace(
            module.procedure("crunch_kernel").directives,
            prologue=Prologue("guardd", "tox_demo_kernel"),
        )
        validate_project(b.project(module), self.bag)

        error = only_error(self.bag)
        assert "prologue 'tox_demo_kernel:guardd' does not exist" == error.message
        assert "without a prologue at all" in error.note

    def test_a_prologue_without_handled_is_rejected(self):
        # the wrapper declares `handled` and branches on it regardless, so a prologue that
        # never sets it leaves that branch reading an undefined logical -- which compiles
        self.checked(guard_arguments=self.guard(
            b.real("values", Intent.IN, "(n)", doc="the data"),
            b.integer("n", Intent.IN, doc="length"),
            b.real("result", Intent.OUT, "(n)", doc="the answer"),
            b.ierr(),
        ))

        error = only_error(self.bag)
        assert "prologue 'guard' has no 'handled' argument" == error.message
        assert "undefined value" in error.note

    def test_a_handled_that_is_not_a_scalar_logical_is_rejected(self):
        self.checked(guard_arguments=self.guard(
            b.integer("n", Intent.IN, doc="length"),
            b.integer("handled", Intent.OUT, doc="not a logical"),
            b.ierr(),
        ))

        assert "'handled' is not a scalar logical" in only_error(self.bag).message

    def test_a_handled_the_prologue_reads_is_rejected(self):
        self.checked(guard_arguments=self.guard(
            b.integer("n", Intent.IN, doc="length"),
            b.logical("handled", Intent.IN, doc="read, not reported"),
            b.ierr(),
        ))

        assert "'handled' is not intent(out)" in only_error(self.bag).message

    def test_a_dummy_the_kernel_does_not_have_becomes_the_wrappers_own(self):
        # what the prologue derives *from* is the allocating tier's vocabulary, not the
        # kernel's: a threshold comes from a percentile, and the kernel takes the threshold
        self.checked(guard_arguments=self.guard(
            b.integer("n", Intent.IN, doc="length"),
            b.real("percentile", Intent.IN, doc="where to put the threshold"),
            b.real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            b.logical("handled", Intent.OUT, doc="dealt with"),
            b.ierr(),
        ))

        assert self.bag.errors == ()

    def test_a_near_miss_of_a_kernel_argument_is_refused_as_a_typo(self):
        # otherwise a misspelling reads as a new argument, and the prologue and the kernel
        # would work from different numbers with nothing to say so
        self.checked(guard_arguments=self.guard(
            b.integer("n", Intent.IN, doc="length"),
            b.real("valuess", Intent.IN, "(n)", doc="one letter off 'values'"),
            b.real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            b.logical("handled", Intent.OUT, doc="dealt with"),
            b.ierr(),
        ))

        error = only_error(self.bag)
        assert "'valuess' looks like a misspelling of 'values'" in error.message
        assert "different values" in error.note

    def test_a_work_array_is_not_a_dummy_the_kernel_does_not_have(self):
        # tmp_scratch is a kernel argument the allocating wrapper prepares, and a prologue
        # may take it -- that is the whole point of running below the allocations
        self.checked(guard_arguments=self.guard(
            b.integer("n", Intent.IN, doc="length"),
            b.real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            b.logical("handled", Intent.OUT, doc="dealt with"),
            b.ierr(),
        ))

        assert self.bag.errors == ()

    def test_what_the_prologue_writes_and_the_kernel_reads_becomes_a_local(self):
        # nobody outside the wrapper is involved, so the caller neither supplies nor
        # receives it -- and both intents already say so
        self.checked(
            guard_arguments=self.guard(
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.OUT, doc="written here"),
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
            ),
            kernel_arguments=[
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.IN, doc="only read here"),
                b.real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
                b.real("result", Intent.OUT, "(n)", doc="the answer"),
            ],
        )

        assert self.bag.errors == ()

    def test_it_is_refused_when_it_also_sizes_something_the_caller_gets(self):
        # then it cannot be dropped, so it stays an intent(in) dummy -- which Fortran will
        # not let the wrapper hand to something that writes it
        self.checked(
            guard_arguments=self.guard(
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.OUT, doc="written here"),
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
            ),
            kernel_arguments=[
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.IN, doc="read here, and sizes result"),
                b.real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
                b.real("result", Intent.OUT, "(room)", doc="the answer"),
            ],
        )

        error = only_error(self.bag)
        assert "cannot make a local" in error.message
        assert "sizes something the caller still" in error.note

    def test_two_alternative_producers_of_one_output_are_accepted(self):
        # kernel intent(out) means the prologue and the kernel each may produce it: the
        # prologue's value is what the caller gets on the handled path
        self.checked(
            guard_arguments=self.guard(
                b.integer("n", Intent.IN, doc="length"),
                b.real("result", Intent.OUT, "(n)", doc="answered here when degenerate"),
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
            ),
        )

        assert self.bag.errors == ()

    def test_a_prologue_may_not_produce_what_sizes_a_work_array(self):
        # tmp_scratch is allocated above the prologue, so its extent cannot be something the
        # prologue only fills afterwards -- the name resolves either way and computes rubbish
        self.checked(
            guard_arguments=self.guard(
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.OUT, doc="how much scratch is needed"),
                b.real("tmp_other", Intent.OUT, "(n)", doc="what makes it run late"),
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
            ),
            kernel_arguments=[
                b.real("values", Intent.IN, "(n)", doc="the data"),
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.INOUT, doc="how much scratch is needed"),
                b.real("tmp_scratch", Intent.OUT, "(room)", doc="scratch"),
                b.real("tmp_other", Intent.OUT, "(n)", doc="what makes it run late"),
                b.real("result", Intent.OUT, "(n)", doc="the answer"),
            ],
        )

        error = only_error(self.bag)
        assert "'tmp_scratch' is sized by 'room', which the prologue only fills afterwards" == error.message
        assert "runs below the work arrays it prepares" in error.note

    def test_a_prologue_may_not_rewrite_what_a_permutation_was_sorted_against(self):
        # the wrapper heapsorts values_perm against values above the prologue, so a prologue
        # that then rewrites values leaves an order describing data that no longer exists.
        # Nothing crashes -- the kernel is simply handed a wrong answer.
        self.checked(
            guard_arguments=self.guard(
                b.real("values", Intent.INOUT, "(n)", doc="rewritten here"),
                b.integer("n", Intent.IN, doc="length"),
                b.real("tmp_other", Intent.OUT, "(n)", doc="what makes it run late"),
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
            ),
            kernel_arguments=[
                b.real("values", Intent.INOUT, "(n)", doc="the data"),
                b.integer("n", Intent.IN, doc="length"),
                b.integer("values_perm", Intent.OUT, "(n)", doc="values, ordered"),
                b.real("tmp_other", Intent.OUT, "(n)", doc="what makes it run late"),
                b.real("result", Intent.OUT, "(n)", doc="the answer"),
            ],
        )

        error = only_error(self.bag)
        assert (
            "'values_perm' is sorted against 'values' before the prologue rewrites it"
            == error.message
        )

    def test_a_prologue_may_rewrite_the_base_of_a_permutation_it_builds_itself(self):
        # the counterpart of the test above: nothing sorts values_perm above the prologue,
        # because the prologue is what builds it, so rewriting values reads nothing early
        self.checked(
            guard_arguments=self.guard(
                b.real("values", Intent.INOUT, "(n)", doc="rewritten here"),
                b.integer("n", Intent.IN, doc="length"),
                b.integer("values_perm", Intent.OUT, "(n)", doc="built here too"),
                b.real("tmp_other", Intent.OUT, "(n)", doc="a work array"),
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
            ),
            kernel_arguments=[
                b.real("values", Intent.INOUT, "(n)", doc="the data"),
                b.integer("n", Intent.IN, doc="length"),
                b.integer("values_perm", Intent.OUT, "(n)", doc="values, ordered"),
                b.real("tmp_other", Intent.OUT, "(n)", doc="a work array"),
                b.real("result", Intent.OUT, "(n)", doc="the answer"),
            ],
        )

        assert self.bag.errors == ()

    def test_a_produced_value_nothing_above_reads_is_accepted(self):
        # `room` is the kernel's business alone -- no allocation, sort or recommend call
        # above the prologue names it, so filling it there is in time
        self.checked(
            guard_arguments=self.guard(
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.OUT, doc="something only the kernel reads"),
                b.real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
            ),
            kernel_arguments=[
                b.real("values", Intent.IN, "(n)", doc="the data"),
                b.integer("n", Intent.IN, doc="length"),
                b.integer("room", Intent.INOUT, doc="something only the kernel reads"),
                b.real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
                b.real("result", Intent.OUT, "(n)", doc="the answer"),
            ],
        )

        assert self.bag.errors == ()


class TestAPrologueFeedingARecommendRoutine:
    """The recommend calls sit above the prologue, so what they are passed must be there."""

    def test_a_producer_input_the_prologue_writes_is_rejected(self, bag):
        from codegen.ir.roles import analyse_project
        from codegen.synthesize import synthesize_wrappers
        from test_synthesize import prologue_feeding_a_producer_module

        project = synthesize_wrappers(
            b.project(prologue_feeding_a_producer_module()), diagnostics=DiagnosticBag()
        ).project
        analyse_project(project, DiagnosticBag())
        validate_project(project, bag)

        assert any(
            "'work_size' is passed 'budget', which the prologue only fills afterwards"
            == d.message
            for d in bag.errors
        ), messages(bag)


class TestPrologueAndTheModeSplit:
    """A prologue is declared once and runs in every wrapper, so it may only name what
    every wrapper has -- which is all of the kernel's arguments unless it splits per mode."""

    def checked(self, dummy, bag, split=True):
        from codegen.ir.directives import Prologue
        from codegen.ir.roles import analyse_project
        from test_synthesize import mode_split_kernel_module

        source = mode_split_kernel_module()
        if not split:
            # same kernel, but the table names no procedure per mode, so nothing splits
            table = source.procedure("detect_patterns_kernel").argument("pattern_mode").doc
            source = _without_the_procedure_column(source, table)
        kernel = source.procedure("detect_patterns_kernel")
        kernel.directives = replace(
            kernel.directives, prologue=Prologue("guard", "tox_demo_kernel")
        )
        source.procedures += (
            b.procedure(
                "guard",
                dummy,
                b.logical("handled", Intent.OUT, doc="dealt with"),
                b.ierr(),
                meta=Meta(summary="Guard", author="AUTHOR"),
            ),
        )
        source._adopt(source.procedures)
        # the real stage order: synthesis analyses each kernel (which is where a mode table
        # becomes roles.mode), then the project-wide pass, then validation
        from codegen.synthesize import synthesize_wrappers

        project = synthesize_wrappers(b.project(source), diagnostics=DiagnosticBag()).project
        analyse_project(project, DiagnosticBag())
        validate_project(project, bag)
        return bag

    def test_the_mode_argument_itself_is_rejected(self):
        # every wrapper fixes the mode in the kernel call and drops it from its signature,
        # so the prologue call would silently lose it
        bag = self.checked(
            b.integer("pattern_mode", Intent.IN, doc="which pattern"), DiagnosticBag()
        )

        assert any("'pattern_mode' is not on every wrapper" in m for m in messages(bag)), messages(bag)

    def test_a_mode_scoped_argument_is_rejected_when_the_mode_splits(self):
        bag = self.checked(
            b.real("threshold", Intent.IN, doc="dosage threshold"), DiagnosticBag()
        )

        assert any("'threshold' is not on every wrapper" in m for m in messages(bag)), messages(bag)

    def test_a_mode_scoped_argument_is_accepted_when_the_mode_does_not_split(self):
        # without a Procedure column there is one wrapper, which carries the argument as an
        # ordinary optional -- DM_REQUIRED_IF_MODE then only means "required at run time"
        bag = self.checked(
            b.real("threshold", Intent.IN, doc="dosage threshold"),
            DiagnosticBag(),
            split=False,
        )

        assert not any("is not on every wrapper" in m for m in messages(bag)), messages(bag)


def _without_the_procedure_column(source, doc):
    """The same kernel module with the mode table's Procedure column removed."""
    from codegen.ir.doc import Doc

    lines = [
        "Which pattern to detect.",
        "| Mode | Value |",
        "|------|-------|",
        "| dosage effect | [[tox_demo_kernel(module):MODE_DOSAGE(variable)]] |",
        "| subfunctionalisation | [[tox_demo_kernel(module):MODE_SUBFUNC(variable)]] |",
    ]
    source.procedure("detect_patterns_kernel").argument("pattern_mode").doc = Doc.parse(lines)
    return source


class TestPrologueThatCouldNeverRun:
    def test_a_prologue_on_a_kernel_with_no_work_arrays_is_rejected(self, bag):
        from test_synthesize import prologue_kernel_module

        source = prologue_kernel_module()
        crunch = source.procedure("crunch_kernel")
        # drop the work array, so no allocating wrapper is generated at all
        crunch.arguments = tuple(a for a in crunch.arguments if a.name != "tmp_scratch")
        validate_project(b.project(source), bag)

        error = only_error(bag)
        assert "on a kernel that generates no allocating wrapper" in error.message
        assert "would never be called" in error.note

    def test_it_is_accepted_when_the_kernel_does_generate_one(self, bag):
        from test_synthesize import prologue_kernel_module

        validate_project(b.project(prologue_kernel_module()), bag)

        assert bag.errors == ()
