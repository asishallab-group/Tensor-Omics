import re

import pytest

from codegen_reworked.diagnostics import DiagnosticBag
from codegen_reworked.ir.doc import Doc
from codegen_reworked.ir.errors import (
    ARG_POS_FACTOR,
    ErrorCatalogue,
    ErrorCode,
    ErrorGroup,
)

import builders as b
from conftest import REPO_ROOT

TOX_ERRORS = REPO_ROOT / "src/tox/tox_errors.F90"


@pytest.fixture(scope="module")
def source():
    return TOX_ERRORS.read_text()


@pytest.fixture
def bag():
    return DiagnosticBag()


def code(name, value, doc="something went wrong"):
    return b.parameter(name, str(value), doc=doc)


@pytest.fixture
def catalogue(bag):
    module = b.module(
        "tox_errors",
        parameters=(
            code("ERR_OK", 0, "no error, operation successful"),
            code("ERR_FILE_OPEN", 101, "could not open file"),
            code("ERR_INVALID_INPUT", 201, "invalid input arguments"),
            code("ERR_ALLOC_FAIL", 301, "memory allocation failed"),
            code("ERR_UNIT_NOT_CONNECTED", 5002, "unit not connected"),
            code("ERR_UNKNOWN", 9999, "unknown error"),
            code("STAT_CONVERGED", 1, "the algorithm converged"),
        ),
    )
    return ErrorCatalogue.from_module(module, bag)


class TestReadingTheCatalogue:
    def test_codes_are_read_with_their_values_and_messages(self, catalogue, bag):
        assert bag.errors == ()
        assert catalogue.code("ERR_INVALID_INPUT").value == 201
        assert catalogue.code("ERR_INVALID_INPUT").message == "invalid input arguments"

    def test_lookup_is_case_insensitive(self, catalogue):
        assert catalogue.code("err_invalid_input").value == 201

    def test_an_unknown_name_is_none(self, catalogue):
        assert catalogue.code("ERR_NOPE") is None

    def test_codes_are_ordered_by_value(self, catalogue):
        values = [c.value for c in catalogue]

        assert values == sorted(values)

    def test_parameters_that_are_not_codes_are_ignored(self, bag):
        module = b.module(
            "tox_errors", parameters=(code("ERR_OK", 0), b.parameter("MAX_STACK", "64"))
        )

        catalogue = ErrorCatalogue.from_module(module, bag)

        assert catalogue.code("MAX_STACK") is None
        assert len(catalogue) == 1

    def test_a_non_constant_code_is_reported(self, bag):
        module = b.module("tox_errors", parameters=(b.parameter("ERR_X", "f(x)"),))

        ErrorCatalogue.from_module(module, bag)

        assert "does not have a constant value" in bag.errors[0].message

    def test_a_non_integer_code_is_reported(self, bag):
        module = b.module("tox_errors", parameters=(b.parameter("ERR_X", "1.5"),))

        ErrorCatalogue.from_module(module, bag)

        assert "is not an integer" in bag.errors[0].message

    def test_an_undocumented_code_warns_and_falls_back_to_its_name(self, bag):
        module = b.module("tox_errors", parameters=(b.parameter("ERR_X", "1"),))

        catalogue = ErrorCatalogue.from_module(module, bag)

        assert bag.errors == ()
        assert "has no documentation to use as its message" in bag.warnings[0].message
        assert catalogue.code("ERR_X").message == "ERR_X"

    def test_two_codes_with_one_value_are_reported(self, bag):
        module = b.module(
            "tox_errors", parameters=(code("ERR_A", 201), code("ERR_B", 201))
        )

        ErrorCatalogue.from_module(module, bag)

        assert "share the value 201" in bag.errors[0].message

    def test_the_message_is_the_first_line_of_the_documentation(self, bag):
        parameter = b.parameter("ERR_X", "1")
        parameter.doc = Doc.parse(["could not open file", "", "check the path"])
        module = b.module("tox_errors", parameters=(parameter,))

        catalogue = ErrorCatalogue.from_module(module, bag)

        assert catalogue.code("ERR_X").message == "could not open file"


class TestErrorsVersusStatuses:
    def test_statuses_are_not_errors(self, catalogue):
        # a status code is an outcome and must never raise
        assert catalogue.code("STAT_CONVERGED").is_status
        assert not catalogue.code("ERR_INVALID_INPUT").is_status

    def test_errors_exclude_statuses_and_ok(self, catalogue):
        names = {c.name for c in catalogue.errors}

        assert "STAT_CONVERGED" not in names
        assert "ERR_OK" not in names
        assert "ERR_INVALID_INPUT" in names

    def test_statuses_lists_only_statuses(self, catalogue):
        assert [c.name for c in catalogue.statuses] == ["STAT_CONVERGED"]

    def test_ok_is_found(self, catalogue):
        assert catalogue.ok.value == 0
        assert catalogue.ok.is_ok

    def test_a_status_may_share_a_value_range_with_nothing(self, catalogue):
        # STAT_CONVERGED is 1: it must not be mistaken for an error of any group
        assert catalogue.code("STAT_CONVERGED").is_status


class TestGroups:
    @pytest.mark.parametrize(
        "value, group",
        [
            (101, ErrorGroup.IO),
            (123, ErrorGroup.IO),
            (200, ErrorGroup.INPUT),
            (210, ErrorGroup.INPUT),
            (301, ErrorGroup.MEMORY),
            (5002, ErrorGroup.RUNTIME),
            (9001, ErrorGroup.INTERNAL),
            (9999, ErrorGroup.INTERNAL),
            (1, ErrorGroup.OTHER),
            (4000, ErrorGroup.OTHER),
        ],
    )
    def test_the_documented_ranges_decide_the_group(self, value, group):
        assert ErrorGroup.of(value) is group

    def test_a_code_reports_its_group(self, catalogue):
        assert catalogue.code("ERR_FILE_OPEN").group is ErrorGroup.IO
        assert catalogue.code("ERR_ALLOC_FAIL").group is ErrorGroup.MEMORY

    def test_groups_lists_only_those_present_in_a_stable_order(self, catalogue):
        assert catalogue.groups() == (
            ErrorGroup.IO,
            ErrorGroup.INPUT,
            ErrorGroup.MEMORY,
            ErrorGroup.RUNTIME,
            ErrorGroup.INTERNAL,
        )


class TestDecoding:
    def test_a_plain_code_decodes_to_itself(self, catalogue):
        decoded = catalogue.decode(201)

        assert decoded.code == 201
        assert decoded.arg_pos == 0
        assert decoded.error.name == "ERR_INVALID_INPUT"
        assert decoded.is_error
        assert not decoded.is_argument_related

    def test_an_argument_position_is_unpacked(self, catalogue):
        # create_err_code packs it as 10000*arg_pos + error
        decoded = catalogue.decode(10000 * 3 + 201)

        assert decoded.code == 201
        assert decoded.arg_pos == 3
        assert decoded.error.name == "ERR_INVALID_INPUT"
        assert decoded.is_argument_related

    def test_ok_is_not_an_error(self, catalogue):
        assert not catalogue.decode(0).is_error

    def test_only_the_code_decides_whether_it_is_an_error(self, catalogue):
        # mirroring is_err, which compares get_err_code(ierr) against ERR_OK
        assert not catalogue.decode(10000 * 2).is_error

    def test_a_four_digit_code_still_decodes(self, catalogue):
        decoded = catalogue.decode(10000 * 2 + 9999)

        assert decoded.error.name == "ERR_UNKNOWN"
        assert decoded.arg_pos == 2

    def test_an_unmapped_code_decodes_without_an_error_object(self, catalogue):
        decoded = catalogue.decode(777)

        assert decoded.code == 777
        assert decoded.error is None
        assert decoded.is_error

    def test_encode_and_decode_round_trip(self, catalogue):
        ierr = catalogue.encode("ERR_ALLOC_FAIL", arg_pos=4)

        decoded = catalogue.decode(ierr)

        assert ierr == 40301
        assert decoded.error.name == "ERR_ALLOC_FAIL"
        assert decoded.arg_pos == 4

    def test_encoding_an_unknown_name_raises(self, catalogue):
        with pytest.raises(KeyError, match="no error code named"):
            catalogue.encode("ERR_NOPE")


class TestAgainstTheRealToxErrors:
    """Guards against the catalogue drifting from the Fortran it mirrors."""

    def test_the_arg_pos_factor_matches_create_err_code(self, source):
        # ARG_POS_FACTOR is duplicated from Fortran, so pin it to the actual expression
        match = re.search(r"ierr\s*=\s*(\d+)\s*\*\s*arg_pos\s*\+\s*error", source)

        assert match is not None, "create_err_code no longer packs arg_pos as N*arg_pos + error"
        assert int(match.group(1)) == ARG_POS_FACTOR

    def test_the_decoders_use_the_same_factor(self, source):
        assert re.search(rf"error\s*=\s*mod\(ierr,\s*{ARG_POS_FACTOR}\)", source) is not None
        assert re.search(rf"arg_pos\s*=\s*ierr/{ARG_POS_FACTOR}", source) is not None

    def test_every_real_error_code_fits_below_the_factor(self, source):
        # a code >= ARG_POS_FACTOR would be indistinguishable from an argument position
        values = [
            int(value)
            for value in re.findall(r"parameter\s*::\s*ERR_\w+\s*=\s*(\d+)", source)
        ]

        assert values, "no ERR_ parameters found"
        assert max(values) < ARG_POS_FACTOR

    def test_the_ok_code_is_named_as_configured(self, source):
        from codegen_reworked.config import CONVENTIONS

        assert f"parameter :: {CONVENTIONS.ok_code} = 0" in source

    def test_the_prefixes_are_the_ones_in_use(self, source):
        from codegen_reworked.config import CONVENTIONS

        assert f"{CONVENTIONS.error_code_prefix}OK" in source


class TestErrorCode:
    def test_a_code_without_documentation_falls_back_to_its_name(self):
        assert ErrorCode("ERR_X", 1).message == "ERR_X"

    def test_ok_is_recognised_by_value(self):
        assert ErrorCode("ERR_OK", 0).is_ok
        assert not ErrorCode("ERR_X", 1).is_ok
