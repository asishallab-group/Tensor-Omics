"""What Python and R share about calling a `DM_OUTPUT_FROM` producer's binding.

Each language writes the call its own way; which arguments get one, what the producer returns,
and how its error is mapped back to the consumer's arguments are the same facts in both, so
they are worked out once, here.
"""

from __future__ import annotations

from ..abi.model import CArgument, CWrapper
from ..ir.entities import Procedure
from ..ir.roles import OutputFromPlan


def auto_plans(wrapper: CWrapper) -> list[tuple[CArgument, OutputFromPlan]]:
    """Each argument the binding computes by calling a producer, with its plan."""
    found = []
    for argument in wrapper:
        roles = argument.source.roles if argument.source else None
        plan = roles.computed_from if roles is not None else None
        if plan is not None and plan.is_automatic:
            found.append((argument, plan))
    return found


def supplied_by(plan: OutputFromPlan, wrapper: CWrapper) -> list[tuple[str, int]]:
    """(producer argument, position of the consumer argument that supplied it), in call order.

    The producer reports a bad input by its own parameter's name, which the consumer's caller
    never wrote; this is what maps it to the consumer's argument. A constant supplies nothing a
    caller passed, so it is not among them. The position numbers the consumer's Fortran dummy
    list, as its error tables do.
    """
    positions = {a.name.lower(): i for i, a in enumerate(wrapper.procedure.arguments, start=1)}
    return [
        (supply.name, positions[supply.argument.lower()])
        for supply in plan.binding_inputs
        if supply.argument is not None and supply.argument.lower() in positions
    ]


def returns_single(producer: Procedure) -> bool:
    """Whether the producer's binding returns one value, rather than one per output."""
    outputs = [
        a for a in producer.arguments
        if a.intent and a.intent.is_output
        and a.name.lower() != "ierr"
        and not (a.roles and a.roles.is_temporary)
    ]
    if producer.result is not None:
        outputs.append(producer.result)
    return len(outputs) == 1
