"""The entity model: project, module, procedure, argument.

This is the keystone of the rework. Every entity is constructible directly:

    argument = Argument("vector", real64, Dimension(("n",)), Intent.INOUT)
    procedure = Procedure("normalize", arguments=[argument])

The previous model could not do this. Its `Procedure` took a `FortranSubroutine` from Ford
and its `Procedure_Argument` took a `FortranVariable`, so exercising any of the logic meant
parsing the whole project first -- which is why none of it was tested. Here the Ford
frontend is one more producer of these objects, and a test fixture is another.

Entities carry structure and documentation only. What an argument *means* -- that it is a
temporary, an extent, a mode -- is a semantic question answered in `roles.py`, and what its
C signature looks like is answered in `abi/`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar, Iterator, Sequence

from ..config import CONVENTIONS, Conventions
from ..diagnostics import SourceLocation
from .constants import ConstantError, ConstantEvaluator
from .directives import Directives
from .doc import Doc
from .types import Dimension, FortranType, Intent


class Entity:
    """Common shape of everything in the model: a name, a parent, a place, a doc.

    `parent` is filled in by the owning entity at construction, which is what lets a
    diagnostic report "argument 'x' in procedure 'p' in module 'm'" from any node.
    """

    entity_kind: ClassVar[str] = "entity"

    name: str
    parent: Entity | None
    location: SourceLocation
    doc: Doc

    def _adopt(self, children: Sequence[Entity]) -> None:
        for child in children:
            child.parent = self

    @property
    def qualified_name(self) -> str:
        names = [self.name]
        parent = self.parent
        while parent is not None and parent.name:
            names.append(parent.name)
            parent = parent.parent
        return ".".join(reversed(names))

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.name!r})"


class Argument(Entity):
    """A dummy argument, or a function result treated as one."""

    entity_kind: ClassVar[str] = "argument"

    def __init__(
        self,
        name: str,
        type: FortranType,
        dimension: Dimension = Dimension(),
        intent: Intent = Intent.INOUT,
        optional: bool = False,
        doc: Doc = Doc(),
        directives: Directives = Directives(),
        attributes: Sequence[str] = (),
        location: SourceLocation = SourceLocation(),
        is_result: bool = False,
    ):
        self.name = name
        self.type = type
        self.dimension = dimension
        self.intent = intent
        self.optional = optional
        self.doc = doc
        self.directives = directives
        self.attributes = tuple(a.lower() for a in attributes)
        self.location = location
        #: True for a function's result, which the C ABI turns into an intent(out) argument
        self.is_result = is_result
        self.parent: Procedure | None = None
        #: Filled in by `roles.analyse`; None until then
        self.roles = None

    @property
    def rank(self) -> int:
        return self.dimension.rank

    @property
    def is_scalar(self) -> bool:
        return self.dimension.is_scalar

    @property
    def is_array(self) -> bool:
        return not self.dimension.is_scalar

    @property
    def has_attribute(self):
        return self.attributes.__contains__

    @property
    def procedure(self) -> Procedure | None:
        return self.parent

    def with_name(self, name: str) -> Argument:
        """A copy under a different name, keeping everything else."""
        copy = Argument(
            name=name,
            type=self.type,
            dimension=self.dimension,
            intent=self.intent,
            optional=self.optional,
            doc=self.doc,
            directives=self.directives,
            attributes=self.attributes,
            location=self.location,
            is_result=self.is_result,
        )
        copy.parent = self.parent
        return copy


class ProcedureKind:
    SUBROUTINE = "subroutine"
    FUNCTION = "function"


@dataclass(frozen=True)
class Meta:
    """Ford meta tags on a procedure or module."""

    summary: str = ""
    author: str = ""
    category: str = ""


class Procedure(Entity):
    """A subroutine or function."""

    entity_kind: ClassVar[str] = "procedure"

    def __init__(
        self,
        name: str,
        arguments: Sequence[Argument] = (),
        result: Argument | None = None,
        doc: Doc = Doc(),
        directives: Directives = Directives(),
        meta: Meta = Meta(),
        location: SourceLocation = SourceLocation(),
        conventions: Conventions = CONVENTIONS,
        allocatable_locals: Sequence[str] = (),
    ):
        self.name = name
        self.arguments = tuple(arguments)
        self.result = result
        self.doc = doc
        self.directives = directives
        self.meta = meta
        self.location = location
        self.conventions = conventions
        #: names of the local variables declared `allocatable`. The body itself is never
        #: read, so this is how the generator sees that a procedure allocates: an
        #: `M_ALLOCATE` needs an allocatable to allocate into, and a kernel may not have one
        #: (`validate._check_kernel_allocates`).
        self.allocatable_locals = tuple(allocatable_locals)
        self.parent: Module | None = None

        self._adopt(self.arguments)
        if result is not None:
            result.parent = self

    @property
    def kind(self) -> str:
        return ProcedureKind.FUNCTION if self.is_function else ProcedureKind.SUBROUTINE

    @property
    def is_function(self) -> bool:
        return self.result is not None

    @property
    def is_exported(self) -> bool:
        """Whether the procedure carries the `category: C-binding` meta tag."""
        return self.meta.category.strip().lower() == (
            self.conventions.c_binding_category.lower()
        )

    @property
    def module(self) -> Module | None:
        return self.parent

    def argument(self, name: str) -> Argument | None:
        """Look an argument up by name, case-insensitively as Fortran resolves names."""
        lowered = name.lower()
        for argument in self.arguments:
            if argument.name.lower() == lowered:
                return argument
        return None

    @property
    def has_error_argument(self) -> bool:
        return self.argument(self.conventions.error_arg) is not None

    @property
    def is_alloc_variant(self) -> bool:
        """Whether this is the allocating half of an `_alloc` / expert pair."""
        return self.name.lower().endswith(self.conventions.alloc_suffix)

    @property
    def alloc_sibling(self) -> Procedure | None:
        """The `<name>_alloc` procedure in the *same module*, if there is one.

        The old generator asked `procedure.find_child(f"{name}_alloc")`, which searches a
        procedure's own children rather than its module, so it never found the sibling and
        the `_expert_c` naming rule never fired.
        """
        if self.module is None or self.is_alloc_variant:
            return None
        return self.module.procedure(f"{self.name}{self.conventions.alloc_suffix}")

    @property
    def has_alloc_sibling(self) -> bool:
        return self.alloc_sibling is not None

    def __iter__(self) -> Iterator[Argument]:
        return iter(self.arguments)


@dataclass
class Parameter:
    """A module-level named constant (`integer, parameter :: ERR_OK = 0`)."""

    entity_kind: ClassVar[str] = "parameter"

    name: str
    type: FortranType | None = None
    expression: str = ""
    doc: Doc = field(default_factory=Doc)
    location: SourceLocation = field(default_factory=SourceLocation)
    parent: Module | None = None


class Module(Entity):
    entity_kind: ClassVar[str] = "module"

    def __init__(
        self,
        name: str,
        procedures: Sequence[Procedure] = (),
        parameters: Sequence[Parameter] = (),
        doc: Doc = Doc(),
        meta: Meta = Meta(),
        location: SourceLocation = SourceLocation(),
        uses: Sequence[str] = (),
    ):
        self.name = name
        self.procedures = tuple(procedures)
        self.parameters = tuple(parameters)
        self.doc = doc
        self.meta = meta
        self.location = location
        #: the modules this one `use`s, sorted; a module that only re-exports its children
        #: has these and no procedures of its own
        self.uses = tuple(uses)
        self.parent: Project | None = None

        self._adopt(self.procedures)
        self._adopt(self.parameters)

    def procedure(self, name: str) -> Procedure | None:
        lowered = name.lower()
        for procedure in self.procedures:
            if procedure.name.lower() == lowered:
                return procedure
        return None

    def parameter(self, name: str) -> Parameter | None:
        lowered = name.lower()
        for parameter in self.parameters:
            if parameter.name.lower() == lowered:
                return parameter
        return None

    @property
    def exported_procedures(self) -> tuple[Procedure, ...]:
        return tuple(p for p in self.procedures if p.is_exported)

    @property
    def has_exports(self) -> bool:
        return bool(self.exported_procedures)

    @property
    def project(self) -> Project | None:
        return self.parent

    def __iter__(self) -> Iterator[Procedure]:
        return iter(self.procedures)


class Project(Entity):
    """Every parsed module. The root of the model."""

    entity_kind: ClassVar[str] = "project"

    def __init__(self, modules: Sequence[Module] = (), name: str = ""):
        self.name = name
        # sorted, so generated output does not depend on filesystem order
        self.modules = tuple(sorted(modules, key=lambda m: m.name.lower()))
        self.doc = Doc()
        self.location = SourceLocation()
        self.parent = None

        self._adopt(self.modules)

    def module(self, name: str) -> Module | None:
        lowered = name.lower()
        for module in self.modules:
            if module.name.lower() == lowered:
                return module
        return None

    def procedure(self, module_name: str, procedure_name: str) -> Procedure | None:
        module = self.module(module_name)
        return None if module is None else module.procedure(procedure_name)

    @property
    def modules_with_exports(self) -> tuple[Module, ...]:
        return tuple(m for m in self.modules if m.has_exports)

    @property
    def parameters(self) -> tuple[Parameter, ...]:
        return tuple(p for module in self.modules for p in module.parameters)

    def constant_values(self, evaluator: ConstantEvaluator | None = None) -> dict[str, object]:
        """Every module parameter that evaluates to a constant, for `DM_DEFAULT`.

        Parameters whose value is not a constant expression are skipped rather than
        reported: most of them are of no interest to any default, and a default that does
        refer to one fails with a message naming that default instead.
        """
        evaluator = evaluator or ConstantEvaluator()
        resolved: dict[str, object] = {}

        # Repeat until nothing new resolves, so a parameter defined in terms of another
        # is picked up regardless of the order modules were parsed in.
        remaining = {p.name: p.expression for p in self.parameters if p.expression}
        while remaining:
            progressed = False
            for name, expression in list(remaining.items()):
                try:
                    value = ConstantEvaluator(resolved, evaluator.intrinsics).evaluate(expression)
                except ConstantError:
                    continue
                resolved[name] = value
                del remaining[name]
                progressed = True
            if not progressed:
                break

        return resolved

    def __iter__(self) -> Iterator[Module]:
        return iter(self.modules)

    def __len__(self) -> int:
        return len(self.modules)
