"""What Python and R share about calling a `DM_OUTPUT_FROM` producer's binding.

Each language writes the call its own way; which arguments get one, what the producer returns,
and how its error is mapped back to the consumer's arguments are the same facts in both, so
they are worked out once, here.
"""

from __future__ import annotations

from ..abi.model import CArgument, CWrapper
from ..ir.entities import Procedure
from ..ir.roles import OutputFromPlan
from ..synthesize import binding_computes


def auto_plans(wrapper: CWrapper) -> list[tuple[CArgument, OutputFromPlan]]:
    """Each argument the binding computes by calling a producer, with its plan, in argument order.

    Only what this tier's binding computes (`synthesize.binding_computes`): in a generated
    pair's expert a data array is the caller's, and calling its producer anyway would overwrite
    what they passed. Python's and R's producer calls, Python's stashes and producer imports
    and R's error tables all follow from this one list. No producer is given another's output
    (`roles._refuses_chaining`), so the order is free.
    """
    return [
        (argument, argument.source.roles.computed_from)
        for argument in wrapper
        if binding_computes(argument, wrapper)
    ]


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
