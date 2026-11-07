import math

import pytest

from codegen.ir.constants import (
    ConstantError,
    ConstantEvaluator,
    Kind,
    tokenize,
)


@pytest.fixture
def evaluate():
    evaluator = ConstantEvaluator(constants={"PI": math.pi, "MAX_ITER": 300})
    return evaluator.evaluate


class TestIntegerLiterals:
    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("0", 0),
            ("123", 123),
            ("300_int32", 300),
            ("300_int64", 300),
            ("42_8", 42),
        ],
    )
    def test_integers(self, evaluate, expression, expected):
        result = evaluate(expression)

        assert result == expected
        assert isinstance(result, int)

    def test_an_unknown_kind_is_rejected(self, evaluate):
        with pytest.raises(ConstantError, match="unknown kind 'banana'"):
            evaluate("1_banana")


class TestRealLiterals:
    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("0.1_real64", 0.1),
            ("0.7", 0.7),
            ("1.", 1.0),
            (".5", 0.5),
            ("1.0e-3", 0.001),
            ("1.0E3", 1000.0),
            ("1.0d0", 1.0),
            ("1.0D3", 1000.0),
            ("1e3", 1000.0),
            ("2.5_real32", 2.5),
        ],
    )
    def test_reals(self, evaluate, expression, expected):
        result = evaluate(expression)

        assert result == pytest.approx(expected)
        assert isinstance(result, float)


class TestLogicalLiterals:
    @pytest.mark.parametrize(
        "expression, expected", [(".true.", True), (".false.", False), (".TRUE.", True)]
    )
    def test_logicals(self, evaluate, expression, expected):
        # the old evaluator replaced these with the empty string, so eval('') raised
        result = evaluate(expression)

        assert result is expected

    def test_a_logical_is_not_mistaken_for_a_real(self, evaluate):
        # '.5' and '.true.' both start with a dot
        assert evaluate(".5") == 0.5
        assert evaluate(".false.") is False


class TestKindedLogicalLiterals:
    """The guide asks for `.true._c_bool`, so a default may be written that way (G-R7)."""

    @pytest.mark.parametrize(
        "expression, expected",
        [(".true._c_bool", True), (".FALSE._C_BOOL", False), (".true._1", True)],
    )
    def test_the_kind_is_read_and_the_value_kept(self, evaluate, expression, expected):
        assert evaluate(expression) is expected

    def test_an_unknown_kind_is_rejected(self, evaluate):
        with pytest.raises(ConstantError, match="unknown kind 'banana'"):
            evaluate(".true._banana")


class TestCharacterLiterals:
    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("'x'", "x"),
            ('"x"', "x"),
            ("'hello world'", "hello world"),
            ("''''", "'"),
            ('""""', '"'),
            ("'it''s'", "it's"),
            ("''", ""),
        ],
    )
    def test_strings(self, evaluate, expression, expected):
        assert evaluate(expression) == expected


class TestNamedConstants:
    def test_a_module_parameter_resolves(self, evaluate):
        assert evaluate("PI") == pytest.approx(math.pi)

    def test_resolution_is_case_insensitive(self, evaluate):
        assert evaluate("pi") == pytest.approx(math.pi)
        assert evaluate("Pi") == pytest.approx(math.pi)

    def test_an_unknown_name_is_reported_not_guessed(self, evaluate):
        with pytest.raises(ConstantError, match="'NOPE' is not a known constant"):
            evaluate("NOPE")

    def test_a_name_containing_a_constant_name_is_not_substituted(self, evaluate):
        # the old evaluator replaced 'pi' inside any identifier, turning
        # 'spike_threshold' into 's3.141592653589793ke_threshold'
        with pytest.raises(ConstantError, match="'spike_threshold' is not a known constant"):
            evaluate("spike_threshold")

    def test_a_constant_named_like_a_prefix_is_matched_whole(self):
        evaluate = ConstantEvaluator(constants={"PI": 3.0, "PIVOT": 7}).evaluate

        assert evaluate("PIVOT") == 7


class TestIntrinsics:
    def test_achar(self, evaluate):
        assert evaluate("achar(9)") == "\t"

    def test_char(self, evaluate):
        assert evaluate("char(65)") == "A"

    def test_ichar_and_iachar(self, evaluate):
        assert evaluate("ichar('A')") == 65
        assert evaluate("iachar('A')") == 65

    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("abs(-3)", 3),
            ("max(1, 2)", 2),
            ("min(1, 2)", 1),
            ("nint(1.6)", 2),
            ("len('abc')", 3),
            ("int(3.9)", 3),
            ("real(3)", 3.0),
            ("sqrt(4.0)", 2.0),
        ],
    )
    def test_whitelisted_intrinsics(self, evaluate, expression, expected):
        assert evaluate(expression) == pytest.approx(expected)

    def test_mod_follows_fortran_truncation(self, evaluate):
        # Fortran's mod truncates toward zero, Python's % floors
        assert evaluate("mod(-7, 3)") == -1

    def test_mod_by_zero_is_reported(self, evaluate):
        with pytest.raises(ConstantError, match="mod by zero"):
            evaluate("mod(1, 0)")

    def test_a_kind_argument_is_accepted(self, evaluate):
        assert evaluate("int(3.9, int32)") == 3
        assert evaluate("real(3, real64)") == 3.0

    def test_nested_calls(self, evaluate):
        assert evaluate("achar(max(9, 8))") == "\t"

    def test_a_non_whitelisted_intrinsic_is_rejected(self, evaluate):
        # the elemental math intrinsics are whitelisted (acos et al.), but a special
        # function like the Bessel functions is not -- an unknown name is reported, not
        # silently turned into garbage the way the old str.replace evaluator did
        with pytest.raises(ConstantError, match="'bessel_j0' is not an intrinsic"):
            evaluate("bessel_j0(1.0_real64)")

    def test_the_error_lists_what_is_available(self, evaluate):
        with pytest.raises(ConstantError, match="achar"):
            evaluate("nope(1)")

    def test_an_intrinsic_call_that_fails_is_reported_with_the_expression(self, evaluate):
        with pytest.raises(ConstantError, match=r"cannot evaluate 'achar\('x'\)'"):
            evaluate("achar('x')")


class TestArithmetic:
    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("1 + 2", 3),
            ("5 - 2", 3),
            ("2 * 3", 6),
            ("1 + 2 * 3", 7),
            ("(1 + 2) * 3", 9),
            ("2.0 * PI", 2 * math.pi),
            ("-1", -1),
            ("- 1.5", -1.5),
            ("+1", 1),
            ("--1", 1),
            ("2 ** 3", 8),
            ("2 ** 3 ** 2", 512),
            ("MAX_ITER + 1", 301),
            ("1.0_real64 / 4.0_real64", 0.25),
        ],
    )
    def test_arithmetic(self, evaluate, expression, expected):
        assert evaluate(expression) == pytest.approx(expected)

    def test_integer_division_truncates_toward_zero(self, evaluate):
        assert evaluate("7 / 2") == 3
        assert evaluate("-7 / 2") == -3

    def test_real_division_does_not_truncate(self, evaluate):
        assert evaluate("7.0 / 2.0") == 3.5

    def test_division_by_zero_is_reported(self, evaluate):
        with pytest.raises(ConstantError, match="division by zero"):
            evaluate("1 / 0")

    def test_real_division_by_zero_is_reported(self, evaluate):
        with pytest.raises(ConstantError, match="division by zero"):
            evaluate("1.0 / 0.0")

    def test_a_string_cannot_be_negated(self, evaluate):
        with pytest.raises(ConstantError, match="cannot negate"):
            evaluate("-'x'")

    def test_a_logical_cannot_be_negated(self, evaluate):
        with pytest.raises(ConstantError, match="cannot negate"):
            evaluate("-.true.")

    def test_incompatible_operands_are_reported(self, evaluate):
        with pytest.raises(ConstantError, match="cannot evaluate"):
            evaluate("1 + 'x'")


class TestKinds:
    def test_a_bare_kind_is_not_a_value(self, evaluate):
        with pytest.raises(ConstantError, match="'int32' is a kind, not a value"):
            evaluate("int32")

    def test_a_kind_cannot_be_used_in_arithmetic(self, evaluate):
        with pytest.raises(ConstantError, match="cannot use a kind in arithmetic"):
            evaluate("int32 + 1")

    def test_kind_equality(self):
        assert Kind("int32") == Kind("int32")
        assert Kind("int32") != Kind("int64")
        assert Kind("int32") != "int32"


class TestMalformedExpressions:
    @pytest.mark.parametrize(
        "expression",
        ["", "   ", "(", "1 +", "1 2", "()", "1 @ 2", "'unterminated", "max(1,"],
    )
    def test_malformed_expressions_are_reported_not_crashed_on(self, evaluate, expression):
        with pytest.raises(ConstantError):
            evaluate(expression)

    def test_a_missing_bracket_is_named(self, evaluate):
        with pytest.raises(ConstantError, match="missing a closing bracket"):
            evaluate("(1 + 2")

    def test_trailing_junk_is_named(self, evaluate):
        with pytest.raises(ConstantError, match="unexpected"):
            evaluate("1 2")

    def test_an_unexpected_character_names_the_expression(self, evaluate):
        with pytest.raises(ConstantError, match=r"cannot parse '1 @ 2'"):
            evaluate("1 @ 2")


class TestCaseInsensitivity:
    """Fortran is case-insensitive, so every spelling must evaluate the same."""

    @pytest.mark.parametrize(
        "expression, expected",
        [
            (".TRUE.", True),
            (".False.", False),
            ("ACHAR(9)", "\t"),
            ("Achar(9)", "\t"),
            ("1.0D0", 1.0),
            ("1.0d0", 1.0),
            ("1.0E3", 1000.0),
            ("300_INT32", 300),
            ("300_Int32", 300),
            ("0.1_REAL64", 0.1),
            ("MAX(1, 2)", 2),
            ("Mod(7, 3)", 1),
        ],
    )
    def test_spellings_agree(self, evaluate, expression, expected):
        assert evaluate(expression) == pytest.approx(expected)

    def test_a_kind_argument_is_ignored_whatever_its_case(self, evaluate):
        assert evaluate("int(3.9, INT32)") == 3


class TestTokenize:
    def test_kinds_of_tokens(self):
        kinds = [token.kind for token in tokenize("achar(9) + PI * 2.0_real64")]

        assert kinds == ["name", "op", "number", "op", "op", "name", "op", "number"]

    def test_positions_are_recorded(self):
        tokens = tokenize("1 + 2")

        assert [token.position for token in tokens] == [0, 2, 4]

    def test_whitespace_is_skipped(self):
        assert len(tokenize("  1  +  2  ")) == 3

    def test_power_is_one_token(self):
        assert [token.kind for token in tokenize("2**3")] == ["number", "power", "number"]

    def test_a_logical_is_one_token(self):
        assert [token.kind for token in tokenize(".true.")] == ["logical"]

    def test_a_kinded_logical_is_one_token(self):
        assert [token.text for token in tokenize(".true._c_bool")] == [".true._c_bool"]


class TestRealWorldDefaults:
    """The DM_DEFAULT values actually used in the codebase."""

    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("0.1_real64", 0.1),
            ("300_int32", 300),
            ("PI", math.pi),
            ("achar(9)", "\t"),
        ],
    )
    def test_defaults_used_on_131_code_gen(self, evaluate, expression, expected):
        assert evaluate(expression) == pytest.approx(expected)


class TestMathIntrinsics:
    """The elemental numeric intrinsics (added 2026-07-29 so PI-style defaults evaluate)."""

    @pytest.mark.parametrize(
        "expression, expected",
        [
            # PI is a module parameter defined this way in f42_math_impl; evaluating its own
            # initialiser is what lets DM_DEFAULT(PI) resolve without injecting a value
            ("4.0_real64*atan(1.0_real64)", math.pi),
            ("atan(1.0_real64)", math.pi / 4),
            ("acos(-1.0_real64)", math.pi),
            ("2.0*asin(1.0)", math.pi),
            ("exp(0.0)", 1.0),
            ("log(1.0)", 0.0),
            ("sqrt(2.0)", math.sqrt(2)),
            ("hypot(3.0, 4.0)", 5.0),
            ("sin(0.0)", 0.0),
            ("cos(0.0)", 1.0),
            ("tanh(0.0)", 0.0),
        ],
    )
    def test_elemental_math(self, evaluate, expression, expected):
        assert evaluate(expression) == pytest.approx(expected)

    def test_a_type_model_inquiry_is_still_not_evaluable(self, evaluate):
        # huge/tiny/epsilon depend on the kind of the argument, which the evaluator does
        # not track, so they are deliberately left out and must be reported, not guessed
        with pytest.raises(ConstantError, match="not an intrinsic"):
            evaluate("huge(1_int32)")


class TestFortranArithmetic:
    """The value is the one a Fortran caller gets, not Python's (issue #221, B9)."""

    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("nint(2.5)", 3),
            ("nint(-2.5)", -3),
            ("nint(2.4999)", 2),
            # floor(x + 0.5) would round this up: the addition itself rounds to 1.0
            ("nint(0.49999999999999994)", 0),
            ("nint(-0.49999999999999994)", 0),
        ],
    )
    def test_nint_rounds_half_away_from_zero(self, evaluate, expression, expected):
        result = evaluate(expression)

        assert result == expected
        assert isinstance(result, int)

    def test_anint_rounds_half_away_from_zero_and_stays_real(self, evaluate):
        assert evaluate("anint(2.5)") == 3.0
        assert evaluate("anint(-2.5)") == -3.0
        assert isinstance(evaluate("anint(2.5)"), float)

    @pytest.mark.parametrize(
        "expression, expected",
        [("2**-1", 0), ("(-2)**-1", 0), ("1**-2", 1), ("(-1)**-1", -1), ("(-1)**-2", 1)],
    )
    def test_an_integer_to_a_negative_power_is_an_integer(self, evaluate, expression, expected):
        result = evaluate(expression)

        assert result == expected
        assert isinstance(result, int)

    def test_zero_to_a_negative_power_is_refused(self, evaluate):
        with pytest.raises(ConstantError, match="zero raised to a negative power"):
            evaluate("0**-1")

    def test_a_real_to_a_negative_power_stays_real(self, evaluate):
        assert evaluate("2.0**-1") == 0.5

    def test_integer_division_is_exact(self, evaluate):
        # int(a / b) goes through a float, which has no 9007199254740993
        assert evaluate("9007199254740993_int64 / 1") == 9007199254740993
        assert evaluate("-9007199254740993_int64 / 2") == -4503599627370496
        assert evaluate("mod(9007199254740993_int64, 2)") == 1

    def test_a_real_mod_is_exact(self, evaluate):
        assert evaluate("mod(7.5, 2.0)") == 1.5
        assert evaluate("mod(-7.5, 2.0)") == -1.5

    @pytest.mark.parametrize(
        "expression",
        ["2147483648", "2147483647 + 1", "-2147483647 - 2", "2**31", "2**63", "46341 * 46341",
         "int(3.0e10)", "nint(1.0e10)", "abs(-2147483647 - 1)"],
    )
    def test_overflowing_the_default_integer_kind_is_refused(self, evaluate, expression):
        with pytest.raises(ConstantError, match="overflows an integer of kind 4"):
            evaluate(expression)

    def test_the_kind_s_own_range_is_allowed(self, evaluate):
        assert evaluate("-2147483647 - 1") == -2147483648
        assert evaluate("2147483647") == 2147483647

    @pytest.mark.parametrize(
        "expression, expected",
        [("2_int64**62", 2**62), ("2147483647_int64 + 1", 2**31),
         ("2147483647 + 1_int64", 2**31), ("int(3.0e10, int64)", 30000000000),
         ("nint(1.0e10, int64)", 10000000000)],
    )
    def test_a_wider_kind_holds_more(self, evaluate, expression, expected):
        assert evaluate(expression) == expected

    def test_overflowing_int64_is_refused(self, evaluate):
        with pytest.raises(ConstantError, match="overflows an integer of kind 8"):
            evaluate("2_int64**63")

    def test_a_huge_power_is_refused_without_being_computed(self, evaluate):
        with pytest.raises(ConstantError, match="overflows every integer kind"):
            evaluate("2_int64**100000000000_int64")

    def test_a_parameter_s_kind_carries_into_the_expression(self):
        from codegen.ir.constants import Integer

        evaluate = ConstantEvaluator(constants={"BIG": Integer(2**40, 8)}).evaluate

        assert evaluate("BIG * 2") == 2**41

    def test_a_negative_real_to_a_real_power_is_refused_not_complex(self, evaluate):
        with pytest.raises(ConstantError, match="negative real raised to a real power"):
            evaluate("(-8.0)**(1.0/3.0)")

    def test_a_negative_real_to_an_integer_power_is_fine(self, evaluate):
        assert evaluate("(-2.0)**3") == -8.0

    @pytest.mark.parametrize("expression", ["1.0e400", "1.0e308 * 10.0"])
    def test_a_real_overflow_is_refused(self, evaluate, expression):
        with pytest.raises(ConstantError, match="overflows a real"):
            evaluate(expression)

    def test_the_result_is_a_plain_int(self, evaluate):
        # the kind is the evaluator's business; what it hands out is an ordinary int
        assert type(evaluate("2 + 3")) is int


class TestSignPrecedence:
    """A sign binds below `**` and `*`, as in Fortran (review M1 of #221's WP1).

    Every expected value is what gfortran 16.1 and ifx 2026.1 give for the same expression
    as a `parameter`. `(-a)**b` used to be computed, so `DM_DEFAULT(-2.0_real64**2)` reached
    Python and R as +4.
    """

    @pytest.mark.parametrize(
        "expression, expected",
        [
            ("-2**2", -4),
            ("- 3**2 + 1", -8),
            ("-1**2", -1),
            ("-2**3**2", -512),
            ("-2*3", -6),
            ("+2**2", 4),
            ("-2.0_real64**2", -4.0),
            # the gfortran/ifx extension: a sign straight after an operator signs the
            # operand that follows, its `**` included
            ("2*-3**2", -18),
            ("1 - -2**2", 5),
            ("2.0_real64**-2**2", 0.0625),
            ("2**-1", 0),
        ],
    )
    def test_a_sign_applies_to_the_whole_power(self, evaluate, expression, expected):
        result = evaluate(expression)

        assert result == pytest.approx(expected)
        assert type(result) is type(expected)

    def test_a_negated_real_power_is_not_a_negative_base(self, evaluate):
        # -(2.0**0.5), which Fortran computes; only (-2.0)**0.5 has no real value
        assert evaluate("-2.0_real64**0.5") == pytest.approx(-math.sqrt(2.0))
        with pytest.raises(ConstantError, match="negative real raised to a real power"):
            evaluate("(-2.0_real64)**0.5")
