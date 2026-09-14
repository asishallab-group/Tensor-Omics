"""Synthesising the validating wrapper from an implementation.

An implementation carries no `ierr` and no validation; the wrapper the generator builds from it
adds both and is the exported entry point. These tests hand-build an implementation and check the
injected wrapper, without going near Ford.
"""

from codegen.ir.directives import (
    Directives,
    Maximum,
    Minimum,
    OutputFrom,
    OutputFromMode,
    RequiredIfMode,
)
from codegen.ir.entities import Meta
from codegen.ir.types import Intent
from codegen.synthesize import synthesize_wrappers

from builders import ierr, integer, logical, module, procedure, project, real


def tmp_suffix_collision_impl_module():
    """tmp_ work arrays whose names also match the shape / mask suffix conventions.

    The tmp_ prefix must win: they are allocated as work arrays, not treated as a shape or
    a mask (which the binding derives rather than allocates).
    """
    return module(
        "tox_demo_impl",
        procedure(
            "work_impl",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            integer("tmp_shape", Intent.OUT, "(2)", doc="scratch that ends in _shape"),
            logical(
                "tmp_selection_mask", Intent.OUT, "(n)", doc="scratch that ends in _mask"
            ),
            real("result", Intent.OUT, "(n)", doc="output"),
            meta=Meta(summary="Work", author="AUTHOR"),
        ),
    )


def mask_count_impl_module():
    """An implementation with a selection mask and its `n_selected_` count (the mask convention)."""
    return module(
        "tox_demo_impl",
        procedure(
            "select_impl",
            real("values", Intent.IN, "(n_vecs)", doc="the data"),
            integer("n_vecs", Intent.IN, doc="number of vectors"),
            logical("vecs_selection_mask", Intent.IN, "(n_vecs)", doc="which vectors"),
            integer("n_selected_vecs", Intent.IN, doc="count of .true. in vecs_selection_mask"),
            real("result", Intent.OUT, "(n_selected_vecs)", doc="output"),
            meta=Meta(summary="Select", author="AUTHOR"),
        ),
    )


def distance_matrix_impl_module():
    """An implementation taking a square distance matrix (the distance-matrix convention)."""
    return module(
        "tox_demo_impl",
        procedure(
            "cluster_impl",
            real("distances", Intent.INOUT, "(n_points, n_points)", doc="distance matrix"),
            integer("n_points", Intent.IN, doc="number of points"),
            integer("out_labels", Intent.OUT, "(n_points)", doc="labels"),
            meta=Meta(summary="Cluster", author="AUTHOR"),
        ),
    )


def output_sizing_computed_impl_module():
    """A recommend-sized value that also gives the extent of a returned array.

    The allocating wrapper must NOT take it over: the caller (and the binding) needs it to
    size what comes back.
    """
    return module(
        "tox_demo_impl",
        procedure(
            "search_impl",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            integer(
                "n_subsets",
                Intent.IN,
                doc="how many subsets fit",
                directives=Directives(
                    output_from=OutputFrom(
                        "work_size", "subset_size", "tox_demo_impl", OutputFromMode.AUTO
                    )
                ),
            ),
            integer("subsets", Intent.OUT, "(n, n_subsets)", doc="the subsets found"),
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
            meta=Meta(summary="Search", author="AUTHOR"),
        ),
    )


def prologue_impl_module(guard_arguments=None, impl_arguments=None):
    """An implementation with a prologue, plus the prologue procedure it names.

    `guard_arguments` replaces the prologue's own dummies, for the cases where what it takes
    is the point: a work array, or a name the implementation does not have (which the validator
    refuses). `impl_arguments` replaces the implementation's, for the cases where what the wrapper
    prepares is the point -- a permutation, say.
    """
    from codegen.ir.directives import Prologue

    if guard_arguments is None:
        guard_arguments = (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            real("result", Intent.OUT, "(n)", doc="the answer"),
            logical("handled", Intent.OUT, doc="whether the call was dealt with"),
            ierr(),
        )

    return module(
        "tox_demo_impl",
        procedure(
            "guard",  # the prologue: may handle the call and skip the implementation
            *guard_arguments,
            meta=Meta(summary="Guard", author="AUTHOR"),
        ),
        procedure(
            "crunch_impl",
            *(impl_arguments if impl_arguments is not None else (
                real("values", Intent.IN, "(n)", doc="the data"),
                integer("n", Intent.IN, doc="length"),
                real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
                real("result", Intent.OUT, "(n)", doc="the answer"),
            )),
            meta=Meta(summary="Crunch", author="AUTHOR"),
            directives=Directives(prologue=Prologue("guard", "tox_demo_impl")),
        ),
    )


def mode_split_impl_module():
    """An implementation with a mode argument whose table names a procedure per mode."""
    mode_table = [
        "Which pattern to detect.",
        "| Mode | Value | Procedure |",
        "|------|-------|-----------|",
        "| dosage effect | [[tox_demo_impl(module):MODE_DOSAGE(variable)]] | detect_dosage_effect |",
        "| subfunctionalisation | [[tox_demo_impl(module):MODE_SUBFUNC(variable)]] | detect_subfunctionalization |",
    ]
    return module(
        "tox_demo_impl",
        procedure(
            "detect_patterns_impl",
            real("genes", Intent.IN, "(n)", doc="gene values"),
            integer("n", Intent.IN, doc="number of genes"),
            integer("pattern_mode", Intent.IN, doc=mode_table),
            real(
                "threshold",
                Intent.IN,
                optional=True,
                directives=Directives(
                    required_if_mode=RequiredIfMode(
                        "pattern_mode", "tox_demo_impl", "MODE_DOSAGE"
                    )
                ),
                doc="dosage threshold",
            ),
            real("result", Intent.OUT, "(n)", doc="pattern score per gene"),
            meta=Meta(summary="Detect paralog patterns", author="AUTHOR"),
        ),
    )

C_BINDING = Meta(summary="a summary", author="AUTHOR", category="C-binding")


def alloc_impl_module():
    """An implementation that needs work arrays: a scratch buffer, a permutation, and a
    recommend-sized buffer -- plus the recommend routine it is sized by."""
    return module(
        "tox_demo_impl",
        procedure(  # a recommend routine: exported, but not an implementation
            "work_size",
            integer("n", Intent.IN, doc="length"),
            integer("wsize", Intent.OUT, doc="recommended work size"),
            meta=C_BINDING,
        ),
        procedure(
            "crunch_impl",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length of `values`"),
            integer("values_perm", Intent.IN, "(n)", doc="ascending permutation of `values`"),
            real("tmp_scratch", Intent.OUT, "(n)", doc="scratch space"),
            real("tmp_work", Intent.OUT, "(wsize)", doc="work buffer"),
            integer(
                "wsize",
                Intent.IN,
                directives=Directives(
                    output_from=OutputFrom(
                        "wsize", "work_size", "tox_demo_impl", OutputFromMode.AUTO
                    )
                ),
                doc="work-buffer size",
            ),
            meta=Meta(summary="Crunch the data", author="AUTHOR"),  # an implementation: not exported
        ),
        path="src/tox/tox_demo_impl.F90",
    )


def optional_producer_input_impl_module():
    """A recommend routine that takes, mandatorily, what the implementation takes optionally.

    The allocating wrapper cannot forward an absent optional into a mandatory dummy, so it
    has to resolve the documented default into a local first.
    """
    from codegen.ir.directives import Default

    return module(
        "tox_demo_impl",
        procedure(  # a recommend routine whose `exact` is not optional
            "work_size",
            integer("n", Intent.IN, doc="length"),
            logical("exact", Intent.IN, doc="whether to size exactly"),
            integer("wsize", Intent.OUT, doc="recommended work size"),
            meta=C_BINDING,
        ),
        procedure(
            "crunch_impl",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length of `values`"),
            logical(
                "exact",
                Intent.IN,
                optional=True,
                directives=Directives(default=Default(".false.")),
                doc="whether to size exactly",
            ),
            real("tmp_work", Intent.OUT, "(wsize)", doc="work buffer"),
            integer(
                "wsize",
                Intent.IN,
                directives=Directives(
                    output_from=OutputFrom(
                        "wsize", "work_size", "tox_demo_impl", OutputFromMode.AUTO
                    )
                ),
                doc="work-buffer size",
            ),
            meta=Meta(summary="Crunch the data", author="AUTHOR"),
        ),
    )


def impl_module():
    """A one-impl module, as an author would write it: no ierr, no validation."""
    return module(
        "tox_demo_impl",
        procedure(
            "scale_vector_impl",
            real("vector", Intent.INOUT, "(n)", doc="the data to scale in place"),
            integer("n", Intent.IN, doc="length of `vector`"),
            real(
                "factor",
                Intent.IN,
                directives=Directives(minimum=Minimum("0.0_real64")),
                doc="scale factor",
            ),
            meta=Meta(summary="Scale a vector", author="AUTHOR"),  # an implementation: not exported
        ),
        path="src/tox/tox_demo_impl.F90",
    )


def tmp_perm_impl_module():
    """An implementation with a `tmp_`-prefixed permutation.

    Unlike a caller-facing `<base>_perm` (which the allocating wrapper seeds and sorts), a
    `tmp_..._perm` is a work buffer the implementation seeds and sorts itself, so the wrapper only
    allocates it -- the tox_normalization quantile-normalisation case.
    """
    return module(
        "tox_demo_impl",
        procedure(
            "rank_impl",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            real("tmp_column", Intent.OUT, "(n)", doc="scratch column"),
            integer("tmp_column_perm", Intent.OUT, "(n)", doc="scratch permutation"),
            real("result", Intent.OUT, "(n)", doc="ranked output"),
            meta=Meta(summary="Rank", author="AUTHOR"),
        ),
    )


def ierr_impl_module():
    """An implementation that itself declares `ierr` -- it propagates a sub-helper's failure."""
    return module(
        "tox_demo_impl",
        procedure(
            "risky_impl",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length of `values`"),
            real("result", Intent.OUT, "(n)", doc="the result"),
            ierr(),
            meta=Meta(summary="A fallible impl", author="AUTHOR"),
        ),
    )


def count_extent_impl_module():
    """An implementation whose extent is really a count: it carries a range that permits zero."""
    return module(
        "tox_demo_impl",
        procedure(
            "pick_impl",
            real("values", Intent.IN, "(n_values)", doc="the data"),
            integer("n_values", Intent.IN, doc="length of `values`"),
            integer(
                "n_selected",
                Intent.IN,
                directives=Directives(
                    minimum=Minimum("0_int32"), maximum=Maximum("n_values")
                ),
                doc="how many are selected (may be zero)",
            ),
            integer(
                "indices",
                Intent.IN,
                "(n_selected)",
                directives=Directives(
                    minimum=Minimum("1_int32"), maximum=Maximum("n_values")
                ),
                doc="selected indices",
            ),
            real("result", Intent.OUT, "(n_values)", doc="output"),
            meta=Meta(summary="Pick", author="AUTHOR"),
        ),
    )


class TestSynthesis:
    def test_a_generated_module_is_injected(self):
        result = synthesize_wrappers(project(impl_module()))

        assert result.project.module("tox_demo") is not None
        # the implementation module stays, inert
        assert result.project.module("tox_demo_impl") is not None

    def test_the_wrapper_drops_the_impl_suffix(self):
        result = synthesize_wrappers(project(impl_module()))

        wrapper = result.project.procedure("tox_demo", "scale_vector")
        assert wrapper is not None
        assert result.project.procedure("tox_demo", "scale_vector_impl") is None

    def test_the_wrapper_is_exported(self):
        result = synthesize_wrappers(project(impl_module()))

        wrapper = result.project.procedure("tox_demo", "scale_vector")
        assert wrapper.is_exported

    def test_the_wrapper_takes_the_impl_arguments_plus_ierr(self):
        result = synthesize_wrappers(project(impl_module()))

        wrapper = result.project.procedure("tox_demo", "scale_vector")
        names = [a.name for a in wrapper.arguments]
        assert names == ["vector", "n", "factor", "ierr"]
        assert wrapper.argument("ierr").intent is Intent.OUT

    def test_the_impl_keeps_no_ierr(self):
        # an implementation that already had an ierr would not get a second one; this one had none
        result = synthesize_wrappers(project(impl_module()))
        impl = result.project.procedure("tox_demo_impl", "scale_vector_impl")

        assert impl.argument("ierr") is None

    def test_a_impl_that_declares_ierr_gets_no_second_one(self):
        result = synthesize_wrappers(project(ierr_impl_module()))

        foo = result.project.procedure("tox_demo", "risky")
        names = [a.name for a in foo.arguments]
        assert names == ["values", "n", "result", "ierr"]

    def test_the_spec_links_wrapper_to_impl(self):
        result = synthesize_wrappers(project(impl_module()))

        (spec,) = result.specs
        assert spec.impl.name == "scale_vector_impl"
        assert spec.validating.name == "scale_vector"
        assert spec.module_name == "tox_demo"

    def test_no_allocating_wrapper_without_work_arrays(self):
        # scale_vector_impl takes no tmp_/perm/computed argument
        result = synthesize_wrappers(project(impl_module()))

        # the lone wrapper takes the plain name, so there is no `_expert` beside it
        assert result.project.procedure("tox_demo", "scale_vector_expert") is None
        assert result.specs[0].allocating is None

    def test_a_non_impl_procedure_is_not_wrapped(self):
        # a recommend routine sitting in the implementation module has no _impl suffix
        mod = module(
            "tox_demo_impl",
            procedure(
                "scale_vector_impl",
                real("vector", Intent.INOUT, "(n)"),
                integer("n", Intent.IN),
                meta=Meta(),
            ),
            procedure(
                "recommended_factor",
                integer("n", Intent.IN),
                real("factor", Intent.OUT),
                meta=Meta(),
            ),
        )
        result = synthesize_wrappers(project(mod))

        assert result.project.procedure("tox_demo", "scale_vector") is not None
        assert result.project.procedure("tox_demo", "recommended_factor") is None
        assert len(result.specs) == 1


class TestAllocatingSynthesis:
    def test_an_allocating_wrapper_is_generated(self):
        result = synthesize_wrappers(project(alloc_impl_module()))

        assert result.project.procedure("tox_demo", "crunch") is not None
        assert result.specs[0].allocating is not None

    def test_it_drops_the_work_arrays_permutation_and_computed_size(self):
        result = synthesize_wrappers(project(alloc_impl_module()))

        alloc = result.project.procedure("tox_demo", "crunch")
        names = [a.name for a in alloc.arguments]
        # values + n survive; values_perm, tmp_scratch, tmp_work, wsize are taken over
        assert names == ["values", "n", "ierr"]

    def test_a_permutation_whose_base_is_not_an_argument_is_kept(self):
        # `ranking_perm` orders something the implementation never receives, so the wrapper cannot
        # seed and sort it -- it is the caller's own data despite the name
        result = synthesize_wrappers(
            project(
                module(
                    "tox_demo_impl",
                    procedure(
                        "rank_impl",
                        real("values", Intent.IN, "(n)", doc="the data"),
                        integer("n", Intent.IN, doc="length"),
                        integer("ranking_perm", Intent.IN, "(n)", doc="a caller-supplied order"),
                        real("tmp_scratch", Intent.OUT, "(n)", doc="scratch"),
                        meta=Meta(summary="Rank", author="AUTHOR"),
                    ),
                )
            )
        )

        alloc = result.project.procedure("tox_demo", "rank")
        assert [a.name for a in alloc.arguments] == ["values", "n", "ranking_perm", "ierr"]

    def test_a_computed_size_that_sizes_a_returned_array_is_kept(self):
        # taking it over would leave neither the caller nor the binding able to size what
        # comes back, so only the tmp_ scratch is taken over here
        result = synthesize_wrappers(project(output_sizing_computed_impl_module()))

        alloc = result.project.procedure("tox_demo", "search")
        names = [a.name for a in alloc.arguments]
        assert names == ["values", "n", "n_subsets", "subsets", "ierr"]

    def test_the_validating_wrapper_keeps_everything(self):
        result = synthesize_wrappers(project(alloc_impl_module()))

        foo = result.project.procedure("tox_demo", "crunch_expert")
        names = [a.name for a in foo.arguments]
        assert names == [
            "values", "n", "values_perm", "tmp_scratch", "tmp_work", "wsize", "ierr"
        ]

    def test_both_variants_are_exported(self):
        result = synthesize_wrappers(project(alloc_impl_module()))

        assert result.project.procedure("tox_demo", "crunch_expert").is_exported
        assert result.project.procedure("tox_demo", "crunch").is_exported


class TestModeSplitSynthesis:
    def result(self):
        return synthesize_wrappers(project(mode_split_impl_module()))

    def test_one_procedure_per_mode_named_from_the_table(self):
        project_ = self.result().project

        assert project_.procedure("tox_demo", "detect_dosage_effect") is not None
        assert project_.procedure("tox_demo", "detect_subfunctionalization") is not None
        # the mode argument is gone, so there is no single runtime-mode procedure
        assert project_.procedure("tox_demo", "detect_patterns") is None

    def test_the_mode_argument_is_dropped(self):
        dosage = self.result().project.procedure("tox_demo", "detect_dosage_effect")

        assert "pattern_mode" not in [a.name for a in dosage.arguments]

    def test_a_required_in_this_mode_argument_is_mandatory(self):
        dosage = self.result().project.procedure("tox_demo", "detect_dosage_effect")

        threshold = dosage.argument("threshold")
        assert threshold is not None
        assert not threshold.optional
        # the now-unconditional requirement drops the directive
        assert threshold.directives.required_if_mode is None

    def test_a_mode_scoped_argument_with_a_default_stays_optional(self):
        # an argument the binding can default is never *required*: there the directive only
        # scopes it to its mode, so the per-mode wrapper keeps it optional
        from dataclasses import replace

        from codegen.ir.directives import Default

        source = mode_split_impl_module()
        impl = source.procedure("detect_patterns_impl")
        threshold = impl.argument("threshold")
        object.__setattr__(
            threshold,
            "directives",
            replace(threshold.directives, default=Default("0.5_real64")),
        )

        dosage = synthesize_wrappers(project(source)).project.procedure(
            "tox_demo", "detect_dosage_effect"
        )

        assert dosage.argument("threshold").optional

    def test_a_required_in_another_mode_argument_is_absent(self):
        subfunc = self.result().project.procedure(
            "tox_demo", "detect_subfunctionalization"
        )

        assert subfunc.argument("threshold") is None
        assert [a.name for a in subfunc.arguments] == ["genes", "n", "result", "ierr"]

    def test_the_mode_fix_is_recorded_on_the_spec(self):
        specs = self.result().specs
        dosage = next(s for s in specs if s.validating.name == "detect_dosage_effect")

        assert dosage.mode_fix.argument == "pattern_mode"
        assert dosage.mode_fix.parameter == "MODE_DOSAGE"
        assert dosage.mode_fix.module == "tox_demo_impl"
        assert dosage.impl.name == "detect_patterns_impl"

    def test_an_ordinary_impl_has_no_mode_fix(self):
        (spec,) = synthesize_wrappers(project(impl_module())).specs

        assert spec.mode_fix is None


class TestGeneratedModuleDocumentation:
    """The generated module publishes what the author wrote about the implementation module.

    It is the module Python imports, R builds its help pages from and Fortran `use`s, so its
    documentation has to be the family's documentation. It used to be a "do not edit" banner,
    which left every published module described by the same one sentence.
    """

    DOC = ["Percentiles and the empirical distribution function.", "", "Second paragraph."]

    def generated(self):
        impl = module(
            "tox_demo_impl",
            procedure(
                "work_impl",
                real("x", Intent.INOUT, doc="the value"),
                meta=Meta(summary="Work", author="AUTHOR"),
            ),
            doc=self.DOC,
            meta=Meta(summary="Descriptive statistics", author="MODULE_AUTHOR"),
        )
        return synthesize_wrappers(project(impl)).project.module("tox_demo")

    def test_the_documentation_is_the_implementation_module_s_own(self):
        assert [line.text for line in self.generated().doc.lines] == self.DOC

    def test_the_summary_is_the_author_s_first_line(self):
        # what Python puts in its module docstring, and Ford in every one-line listing
        assert self.generated().doc.summary == self.DOC[0]

    def test_the_meta_tags_come_across_too(self):
        assert self.generated().meta.summary == "Descriptive statistics"
        assert self.generated().meta.author == "MODULE_AUTHOR"

    def test_no_generated_note_is_added_to_the_ir(self):
        """It is per-file, so each emitter says it in its own words; carried here, Python
        would print it once from the doc and once from its own trailer."""
        assert not any(
            "do not edit" in line.text.lower() for line in self.generated().doc.lines
        )

    def test_a_re_export_parent_carries_its_parent_s_documentation(self):
        parent, first, second = split_family_modules()
        parent = module(
            "tox_demo_impl",
            uses=("tox_demo_left_impl", "tox_demo_right_impl"),
            doc=["Gathers the demo family."],
        )

        generated = synthesize_wrappers(
            project(parent, first, second)
        ).project.module("tox_demo")

        assert generated.doc.summary == "Gathers the demo family."


def split_family_modules():
    """A family split across two implementation modules, gathered by a parent that re-exports them.

    The parent holds no procedures -- only `use` lines -- which is what makes it a
    re-export rather than an implementation module of its own.
    """
    first = module(
        "tox_demo_left_impl",
        procedure(
            "left_impl",
            real("x", Intent.INOUT, doc="the value"),
            meta=Meta(summary="Left", author="AUTHOR"),
        ),
    )
    second = module(
        "tox_demo_right_impl",
        procedure(
            "right_impl",
            real("y", Intent.INOUT, doc="the value"),
            meta=Meta(summary="Right", author="AUTHOR"),
        ),
    )
    parent = module(
        "tox_demo_impl",
        uses=("tox_demo_left_impl", "tox_demo_right_impl"),
    )
    return [parent, first, second]


class TestReexportSynthesis:
    def result(self):
        return synthesize_wrappers(project(*split_family_modules()))

    def test_the_parent_is_generated_alongside_its_children(self):
        result = self.result()

        assert result.project.module("tox_demo") is not None
        assert result.project.module("tox_demo_left") is not None
        assert result.project.module("tox_demo_right") is not None

    def test_the_parent_re_exports_the_generated_children(self):
        parent = self.result().project.module("tox_demo")

        assert parent.uses == ("tox_demo_left", "tox_demo_right")
        assert parent.procedures == ()

    def test_the_parent_is_reported_as_a_re_export(self):
        result = self.result()

        assert result.reexports == ("tox_demo",)
        assert result.generated_names == {"tox_demo", "tox_demo_left", "tox_demo_right"}

    def test_a_child_that_generates_nothing_is_not_re_exported(self):
        """An implementation module of constants or recommend routines has no generated counterpart."""
        parent, first, _ = split_family_modules()
        parent = module(
            "tox_demo_impl",
            uses=("tox_demo_left_impl", "tox_demo_constants_impl"),
        )
        constants = module("tox_demo_constants_impl")

        result = synthesize_wrappers(project(parent, first, constants))

        assert result.project.module("tox_demo").uses == ("tox_demo_left",)
        assert result.project.module("tox_demo_constants") is None

    def test_a_parent_may_gather_another_parent(self):
        parent, first, second = split_family_modules()
        grandparent = module("tox_demo_all_impl", uses=("tox_demo_impl",))

        result = synthesize_wrappers(project(grandparent, parent, first, second))

        assert result.project.module("tox_demo_all").uses == ("tox_demo",)

    def test_a_impl_module_that_uses_another_is_not_a_re_export(self):
        """Only a module with no procedures of its own gathers; an implementation that `use`s a
        sibling (for a constant, say) still generates its own wrappers."""
        _, first, second = split_family_modules()
        first = module(
            "tox_demo_left_impl",
            procedure(
                "left_impl",
                real("x", Intent.INOUT, doc="the value"),
                meta=Meta(summary="Left", author="AUTHOR"),
            ),
            uses=("tox_demo_right_impl",),
        )

        result = synthesize_wrappers(project(first, second))

        assert result.reexports == ()
        assert result.project.procedure("tox_demo_left", "left") is not None


def mode_impl_module():
    """An implementation with a runtime mode argument: a table, but no per-mode procedure column.

    Without that third column the mode stays an argument of one wrapper, which is what the
    generated membership check is for.
    """
    # two columns: a third one is the per-mode-procedure opt-in, which is what this is not
    mode_table = [
        "which way to compute it",
        "| Mode | Value |",
        "|------|-------|",
        "| fast | [[tox_demo_impl(module):MODE_FAST(variable)]] |",
        "| exact | [[tox_demo_impl(module):MODE_EXACT(variable)]] |",
    ]
    return module(
        "tox_demo_impl",
        procedure(
            "compute_impl",
            real("values", Intent.INOUT, "(n)", doc="the data"),
            integer("mode", Intent.IN, doc=mode_table),
            integer("n", Intent.IN, doc="length of `values`"),
            meta=Meta(summary="Compute", author="AUTHOR"),
        ),
    )


def prologue_feeding_a_producer_module():
    """An implementation whose recommend routine is passed something the prologue writes.

    The recommend call is emitted above the prologue, so it would read the value before the
    prologue fills it -- and a name resolves the same either way, so it compiles.
    """
    from codegen.ir.directives import OutputFrom, OutputFromMode, Prologue

    return module(
        "tox_demo_impl",
        procedure(  # the recommend routine, sized from `budget`
            "work_size",
            integer("budget", Intent.IN, doc="how much room to plan for"),
            integer("wsize", Intent.OUT, doc="recommended work size"),
            meta=C_BINDING,
        ),
        procedure(
            "guard",
            integer("n", Intent.IN, doc="length"),
            integer("budget", Intent.OUT, doc="decided here -- too late for work_size"),
            real("tmp_work", Intent.OUT, "(wsize)", doc="work buffer"),
            logical("handled", Intent.OUT, doc="whether the call was dealt with"),
            ierr(),
            meta=Meta(summary="Guard", author="AUTHOR"),
        ),
        procedure(
            "crunch_impl",
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length of `values`"),
            integer("budget", Intent.INOUT, doc="how much room to plan for"),
            real("tmp_work", Intent.OUT, "(wsize)", doc="work buffer"),
            integer(
                "wsize",
                Intent.IN,
                directives=Directives(
                    output_from=OutputFrom(
                        "wsize", "work_size", "tox_demo_impl", OutputFromMode.AUTO
                    )
                ),
                doc="work-buffer size",
            ),
            real("result", Intent.OUT, "(n)", doc="the answer"),
            meta=Meta(summary="Crunch the data", author="AUTHOR"),
            directives=Directives(prologue=Prologue("guard", "tox_demo_impl")),
        ),
    )


class TestAnAllocatingWrapperWithoutWorkArrays:
    """`foo` exists when the two signatures would differ -- not only when something is
    taken over. A prologue that takes nothing over but asks for an argument of its own is the
    case that needs it: a degenerate-input guard with a tolerance to judge by."""

    def synthesised(self):
        from codegen.diagnostics import DiagnosticBag

        from builders import ierr, integer, logical, project, real

        # the prologue takes nothing over: `result` is intent(out) on both, so the two are
        # alternative producers of it rather than one feeding the other. `tolerance` is the
        # only thing that makes the signatures differ.
        guard = (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            real("tolerance", Intent.IN, doc="how degenerate is too degenerate"),
            real("result", Intent.OUT, "(n)", doc="answered here when degenerate"),
            logical("handled", Intent.OUT, doc="dealt with"),
            ierr(),
        )
        impl = (
            real("values", Intent.IN, "(n)", doc="the data"),
            integer("n", Intent.IN, doc="length"),
            real("result", Intent.OUT, "(n)", doc="the answer"),
        )
        return synthesize_wrappers(
            project(prologue_impl_module(guard, impl)), diagnostics=DiagnosticBag()
        )

    def test_nothing_is_taken_over_at_all(self):
        # the premise: without this the test would pass on the old rule too
        from codegen.config import CONVENTIONS
        from codegen.synthesize import taken_over_arguments

        spec = next(s for s in self.synthesised().specs if s.validating.name == "crunch_expert")
        assert taken_over_arguments(spec.impl.arguments, CONVENTIONS, spec.prologue) == []

    def test_the_allocating_wrapper_is_generated_anyway(self):
        module = self.synthesised().project.module("tox_demo")

        assert module.procedure("crunch") is not None
        assert [a.name for a in module.procedure("crunch").arguments] == [
            "values", "n", "result", "tolerance", "ierr",
        ]

    def test_the_expert_wrapper_does_not_take_it(self):
        module = self.synthesised().project.module("tox_demo")

        assert [a.name for a in module.procedure("crunch_expert").arguments] == [
            "values", "n", "result", "ierr",
        ]

    def test_it_counts_as_doing_more_than_allocating(self):
        # so both tiers reach Python and R: they genuinely differ
        spec = next(s for s in self.synthesised().specs if s.validating.name == "crunch_expert")

        assert spec.alloc_does_more
