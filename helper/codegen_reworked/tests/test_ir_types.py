import pytest

from codegen_reworked.ir.types import (
    BaseType,
    CharacterLength,
    Dimension,
    FortranType,
    Intent,
    LengthKind,
    UnsupportedTypeError,
    split_top_level,
)


class TestIntent:
    @pytest.mark.parametrize(
        "intent, is_input, is_output",
        [
            (Intent.IN, True, False),
            (Intent.OUT, False, True),
            (Intent.INOUT, True, True),
        ],
    )
    def test_direction(self, intent, is_input, is_output):
        assert intent.is_input is is_input
        assert intent.is_output is is_output


class TestBaseType:
    @pytest.mark.parametrize(
        "text, expected",
        [
            ("integer", BaseType.INTEGER),
            ("REAL", BaseType.REAL),
            ("  Character  ", BaseType.CHARACTER),
            ("logical", BaseType.LOGICAL),
            ("complex", BaseType.COMPLEX),
            ("type", BaseType.DERIVED),
        ],
    )
    def test_parse_is_case_and_space_insensitive(self, text, expected):
        assert BaseType.parse(text) is expected

    def test_parse_rejects_unknown_types(self):
        with pytest.raises(UnsupportedTypeError, match="unsupported type 'quaternion'"):
            BaseType.parse("quaternion")


class TestCharacterLength:
    @pytest.mark.parametrize(
        "text, kind, expr",
        [
            ("*", LengthKind.ASSUMED, None),
            (":", LengthKind.DEFERRED, None),
            ("n", LengthKind.EXPLICIT, "n"),
            ("123", LengthKind.EXPLICIT, "123"),
            ("  n_chars  ", LengthKind.EXPLICIT, "n_chars"),
        ],
    )
    def test_parse(self, text, kind, expr):
        length = CharacterLength.parse(text)

        assert length.kind is kind
        assert length.expr == expr

    @pytest.mark.parametrize(
        "text, is_constant",
        [("123", True), ("n", False), ("*", False), (":", False)],
    )
    def test_is_constant_distinguishes_a_literal_from_a_reference(self, text, is_constant):
        assert CharacterLength.parse(text).is_constant is is_constant

    def test_assumed_and_deferred_flags(self):
        assert CharacterLength.parse("*").is_assumed
        assert CharacterLength.parse(":").is_deferred
        assert not CharacterLength.parse("n").is_assumed
        assert not CharacterLength.parse("n").is_deferred

    @pytest.mark.parametrize("text", ["*", ":", "n", "123"])
    def test_str_round_trips(self, text):
        assert str(CharacterLength.parse(text)) == text

    def test_explicit_length_requires_an_expression(self):
        with pytest.raises(ValueError, match="explicit length needs an expression"):
            CharacterLength(LengthKind.EXPLICIT)

    def test_non_explicit_length_rejects_an_expression(self):
        with pytest.raises(ValueError, match="explicit length needs an expression"):
            CharacterLength(LengthKind.ASSUMED, "n")


class TestFortranType:
    def test_numeric_types_require_a_kind(self):
        # a default kind is processor-dependent, so it has no defensible C mapping
        for base in (BaseType.INTEGER, BaseType.REAL, BaseType.COMPLEX):
            with pytest.raises(ValueError, match="kind is required"):
                FortranType(base)

    def test_logical_needs_no_kind(self):
        assert FortranType(BaseType.LOGICAL).kind is None

    def test_character_requires_a_length(self):
        with pytest.raises(ValueError, match="character type needs a length"):
            FortranType(BaseType.CHARACTER)

    def test_non_character_rejects_a_length(self):
        with pytest.raises(ValueError, match="length is meaningless"):
            FortranType(BaseType.INTEGER, kind="int32", length=CharacterLength.parse("n"))

    def test_derived_type_requires_a_name(self):
        with pytest.raises(ValueError, match="derived type needs a name"):
            FortranType(BaseType.DERIVED)

    def test_only_derived_type_carries_a_name(self):
        with pytest.raises(ValueError, match="derived type needs a name"):
            FortranType(BaseType.LOGICAL, derived_name="point")

    def test_character_rejects_a_non_c_kind(self):
        with pytest.raises(ValueError, match="unsupported character kind"):
            FortranType(BaseType.CHARACTER, kind="utf8", length=CharacterLength.parse("n"))

    def test_types_are_immutable(self):
        # the old model had to mutate a type after construction to fix up intent
        type_ = FortranType(BaseType.INTEGER, kind="int32")

        with pytest.raises(Exception):
            type_.kind = "int64"

    def test_types_compare_by_value(self):
        assert FortranType(BaseType.REAL, kind="real64") == FortranType(
            BaseType.REAL, kind="real64"
        )


class TestConversion:
    def test_character_always_needs_conversion(self):
        assert FortranType(BaseType.CHARACTER, length=CharacterLength.parse("n")).needs_conversion

    def test_default_logical_needs_conversion(self):
        assert FortranType(BaseType.LOGICAL).needs_conversion

    def test_c_bool_logical_needs_no_conversion(self):
        # the one case where the copy can be skipped
        type_ = FortranType(BaseType.LOGICAL, kind="c_bool")

        assert type_.is_c_bool
        assert not type_.needs_conversion

    def test_c_bool_detection_is_case_insensitive(self):
        assert FortranType(BaseType.LOGICAL, kind="C_BOOL").is_c_bool

    def test_a_kinded_but_non_c_bool_logical_still_converts(self):
        type_ = FortranType(BaseType.LOGICAL, kind="int8")

        assert not type_.is_c_bool
        assert type_.needs_conversion

    @pytest.mark.parametrize(
        "type_",
        [
            FortranType(BaseType.INTEGER, kind="int32"),
            FortranType(BaseType.REAL, kind="real64"),
            FortranType(BaseType.COMPLEX, kind="real64"),
        ],
    )
    def test_numeric_types_need_no_conversion(self, type_):
        assert not type_.needs_conversion


class TestFortranTypeStr:
    @pytest.mark.parametrize(
        "type_, expected",
        [
            (FortranType(BaseType.INTEGER, kind="int32"), "integer(int32)"),
            (FortranType(BaseType.REAL, kind="real64"), "real(real64)"),
            (FortranType(BaseType.LOGICAL), "logical"),
            (FortranType(BaseType.LOGICAL, kind="c_bool"), "logical(c_bool)"),
            (
                FortranType(BaseType.CHARACTER, length=CharacterLength.parse("*")),
                "character(len=*)",
            ),
            (
                FortranType(
                    BaseType.CHARACTER, kind="c_char", length=CharacterLength.parse("1")
                ),
                "character(len=1, kind=c_char)",
            ),
            (FortranType(BaseType.DERIVED, derived_name="point"), "type(point)"),
        ],
    )
    def test_str(self, type_, expected):
        assert str(type_) == expected


class TestDimension:
    def test_scalar_has_no_extents(self):
        dimension = Dimension()

        assert dimension.rank == 0
        assert dimension.is_scalar
        assert not dimension

    @pytest.mark.parametrize(
        "text, extents",
        [
            ("(n, m)", ("n", "m")),
            ("n, m", ("n", "m")),
            ("(:)", (":",)),
            ("(:, :)", (":", ":")),
            ("(*)", ("*",)),
            ("(  n_genes ,  n_tissues  )", ("n_genes", "n_tissues")),
            ("()", ()),
            ("", ()),
            ("(1:n)", ("1:n",)),
        ],
    )
    def test_parse(self, text, extents):
        assert Dimension.parse(text).extents == extents

    def test_parse_respects_nesting(self):
        # the old generator split on every comma and carried a TODO about this
        dimension = Dimension.parse("(int(n, int32), m)")

        assert dimension.extents == ("int(n, int32)", "m")
        assert dimension.rank == 2

    def test_extents_are_stripped_on_construction(self):
        assert Dimension((" n ", "m ")).extents == ("n", "m")

    @pytest.mark.parametrize(
        "text, assumed_shape, assumed_size, explicit",
        [
            ("(n)", False, False, True),
            ("(n, m)", False, False, True),
            ("(:)", True, False, False),
            ("(:, :)", True, False, False),
            ("(*)", False, True, False),
            ("(n, *)", False, True, False),
            ("()", False, False, False),
        ],
    )
    def test_shape_classification(self, text, assumed_shape, assumed_size, explicit):
        dimension = Dimension.parse(text)

        assert dimension.is_assumed_shape is assumed_shape
        assert dimension.is_assumed_size is assumed_size
        assert dimension.is_explicit_shape is explicit

    def test_index_of_finds_an_extent(self):
        dimension = Dimension.parse("(n_genes, n_tissues)")

        assert dimension.index_of("n_tissues") == 1
        assert dimension.index_of("n_absent") is None

    def test_slicing_yields_a_dimension(self):
        dimension = Dimension.parse("(strlen, n, m)")

        assert dimension[1:] == Dimension(("n", "m"))
        assert dimension[0] == "strlen"

    def test_iteration_and_length(self):
        dimension = Dimension.parse("(n, m)")

        assert list(dimension) == ["n", "m"]
        assert len(dimension) == 2

    def test_with_extents_returns_a_new_dimension(self):
        original = Dimension.parse("(:)")
        replaced = original.with_extents(["n"])

        assert replaced == Dimension(("n",))
        assert original == Dimension((":",))

    @pytest.mark.parametrize("text, expected", [("(n, m)", "(n, m)"), ("()", "")])
    def test_str(self, text, expected):
        assert str(Dimension.parse(text)) == expected

    def test_dimensions_compare_by_value(self):
        assert Dimension.parse("(n)") == Dimension(("n",))


class TestSplitTopLevel:
    @pytest.mark.parametrize(
        "text, expected",
        [
            ("a, b", ["a", "b"]),
            ("a", ["a"]),
            ("int(n, int32), m", ["int(n, int32)", "m"]),
            ("f(g(a, b), c), d", ["f(g(a, b), c)", "d"]),
            ("a(1), b[2, 3]", ["a(1)", "b[2, 3]"]),
            ("'x, y', z", ["'x, y'", "z"]),
            ('"x, y", z', ['"x, y"', "z"]),
            ("  a  ,  b  ", ["a", "b"]),
        ],
    )
    def test_split(self, text, expected):
        assert split_top_level(text) == expected

    def test_custom_separator(self):
        assert split_top_level("a:b", separator=":") == ["a", "b"]

    @pytest.mark.parametrize("text", ["f(a", "a)", "f(a))"])
    def test_unbalanced_brackets_are_rejected(self, text):
        with pytest.raises(ValueError, match="unbalanced brackets"):
            split_top_level(text)

    def test_unterminated_string_is_rejected(self):
        with pytest.raises(ValueError, match="unterminated string"):
            split_top_level("'a, b")
