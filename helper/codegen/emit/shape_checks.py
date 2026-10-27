"""Which axes of the arrays a caller passes the binding compares, and with what (#221 A-shape).

Fortran cannot check that an actual argument has the shape its explicit-shape dummy declares:
the wrapper receives an address, and a short buffer is read past. So the binding checks every
input axis it can, before the call, by one rule:

**every input axis that names an extent is compared with that extent's value, except the axis
the value was read from**; and **an axis the declaration fixes at a constant is compared with
that constant.**

An extent's value is read off the first input array that carries it (`extent_provider`), counted
from a mask (`n_selected_<x>`), computed by a `DM_OUTPUT_FROM` producer, or passed by the caller;
optional arrays are compared when given. An axis the declaration computes from other arguments
(`sum(reps_per_tissue)`) is not the binding's to check: C is given its actual extent and the C
wrapper compares it with the declaration (`CArgument.checked_against`), so here it names that
passed extent, which is read off this very axis.

This module decides; Python and R each write the comparisons out in their own language, so the
two cannot disagree about which axes are compared, and Python's derivation of an extent and its
check share one `extent_provider`.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..abi import extents
from ..abi.model import CArgument, Conversion, CWrapper, Origin


@dataclass(frozen=True)
class Position:
    """Where an input carries an extent: one of its axes, or all of its elements."""

    array: CArgument
    #: the axis as the binding sees the array -- a string array's length is not one -- or None
    #: for the element count of an array whose shape travels in a shape argument
    axis: int | None


@dataclass(frozen=True)
class Expected:
    """The value a compared axis must have, and where it comes from."""

    #: `"axis"` (read off another input), `"mask"` (its count of true values), `"variable"`
    #: (an argument the caller passes, or one a producer computes), `"constant"`
    origin: str
    #: the extent argument, for every origin but a constant
    extent: str | None = None
    #: where the value was read from, for `"axis"`
    position: Position | None = None
    #: the mask counted, for `"mask"`
    mask: str | None = None
    #: the value, for `"constant"`
    value: int | None = None


@dataclass(frozen=True)
class ShapeCheck:
    position: Position
    expected: Expected


def binding_axes(argument: CArgument) -> list[str]:
    """The extents of `argument` as the binding sees its axes: a string's length is not one."""
    extents_ = list(argument.dimension.extents)
    if argument.type.is_character and extents_:
        extents_ = extents_[1:]
    return extents_


def declared_shape(argument: CArgument) -> str:
    """The shape the author declared, axis by axis, for a message."""
    shown = list(argument.shown_extents)
    if argument.type.is_character and shown:
        shown = shown[1:]
    return ", ".join(shown)


def axes_of(extent: str, owner: CArgument) -> list[int]:
    """Which of `owner`'s axes carry `extent` -- every one, not the first.

    `distances(n_points, n_points)` claims to know `n_points` twice, and checking only the first
    axis lets a non-square array through to a view that reads, and for an `intent(inout)`
    argument writes, past the buffer.
    """
    return [index for index, name in enumerate(binding_axes(owner)) if name == extent]


def _carries(owner: CArgument | None) -> bool:
    """Whether `owner` is an array the caller passes, so its shape is there to read."""
    return (owner is not None and owner.is_array and owner.intent.is_input
            and not owner.is_temporary and owner.conversion is not Conversion.MODE
            and not (owner.type.is_character and owner.rank <= 1))


def extent_provider(extent: CArgument, wrapper: CWrapper) -> Position | None:
    """The input an extent's value is read off: the first one carrying it that must be passed.

    Not an optional, which may be absent, nor a work array, which the binding allocates *from*
    the extent. An array whose shape travels separately carries it as its element count.
    """
    roles = extent.source.roles if extent.source else None
    if roles is None:
        return None
    for owner in roles.extent_of:
        c_owner = wrapper.argument(owner.name)
        if not _carries(c_owner) or c_owner.optional:
            continue
        if c_owner.shape_arg is not None:
            return Position(c_owner, None)
        axes = axes_of(extent.name, c_owner)
        if axes:
            return Position(c_owner, axes[0])
    return None


def _expected(extent: CArgument, wrapper: CWrapper) -> Expected | None:
    roles = extent.source.roles if extent.source else None
    if roles is None or not roles.is_extent:
        return None
    if roles.is_computed:
        return Expected("variable", extent=extent.name)
    if roles.is_mask_count:
        return Expected("mask", extent=extent.name, mask=roles.mask_count_of.name)
    if roles.is_inferable_extent:
        position = extent_provider(extent, wrapper)
        if position is None:
            return None
        return Expected("axis", extent=extent.name, position=position)
    return Expected("variable", extent=extent.name)


def shape_checks(wrapper: CWrapper) -> list[ShapeCheck]:
    """Every comparison the binding makes before the call, in argument order.

    `wrapper` is the binding's model, every module constant in it already its value
    (`CWrapper.with_constant_sizes_evaluated`); `check_binding_extents` has refused a constant
    axis the grammar cannot evaluate.
    """
    checks = []
    for array in wrapper:
        if not _carries(array):
            continue
        if array.shape_arg is not None:
            # its extents travel in the shape argument; the extents it is declared with, if any,
            # are what its element count must equal
            position = Position(array, None)
            for name in (array.source.dimension.extents if array.source else ()):
                extent = wrapper.argument(name)
                if extent is None or extent.origin is Origin.EXTENT:
                    continue
                expected = _expected(extent, wrapper)
                if expected is not None and expected.position != position:
                    checks.append(ShapeCheck(position, expected))
            continue
        for axis, text in enumerate(binding_axes(array)):
            position = Position(array, axis)
            extent = wrapper.argument(text) if text.isidentifier() else None
            if extent is not None:
                if extent.origin is Origin.EXTENT:
                    # read off this very axis, for C: an assumed shape, or an axis the
                    # declaration computes, which the C wrapper compares (`checked_against`)
                    continue
                expected = _expected(extent, wrapper)
                if expected is not None and expected.position != position:
                    checks.append(ShapeCheck(position, expected))
                continue
            value = extents.evaluate(extents.parse(text))
            checks.append(ShapeCheck(position, Expected("constant", value=value)))
    return checks


def reason(expected: Expected) -> str:
    """Why an axis must have its value, up to the value itself: `'matrix' implies n_cols == `."""
    if expected.origin == "axis":
        return f"'{expected.position.array.name}' implies {expected.extent} == "
    if expected.origin == "mask":
        return f"'{expected.mask}' implies {expected.extent} == "
    if expected.origin == "variable":
        return f"{expected.extent} == "
    return "its declared shape fixes it at "
