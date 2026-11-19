import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    bitmask_n_words, bitmask_test, bitmask_from_logical, bitmask_to_logical,
    bitmasks_from_logical_2D, bitmasks_to_logical_2D,
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT

# The values are tested in Fortran. These tests check that each procedure can be called, its
# return type and shape, and its documented errors, plus one hand-derived value per procedure:
# it proves the binding passes the layout through unchanged. Flag i is bit mod(i-1, 32) of word
# (i-1)//32 + 1, least significant first, so the flags [T, T, F, T] are the word 0b1011 = 11.
FLAGS_1011 = np.array([True, True, False, True])
WORD_1011 = 11


def decode_with_numpy(words, n_bits):
    """The published layout, decoded the way the module documentation says numpy does it."""
    words = np.ascontiguousarray(words, dtype=np.int32)
    return np.unpackbits(words.view(np.uint8), bitorder="little")[:n_bits].astype(bool)


def test_bitmask_n_words():
    n_words = bitmask_n_words(4)
    assert isinstance(n_words, int), f"expected an int, got {type(n_words).__name__}"
    assert n_words == 1, f"4 flags take 1 word, got {n_words}"
    assert bitmask_n_words(32) == 1, "32 flags fill exactly 1 word"
    assert bitmask_n_words(33) == 2, "33 flags take 2 words"

    # an empty mask takes no words, as the paralog code's mask_chunk_count did
    assert bitmask_n_words(0) == 0, "0 flags take 0 words"

    assert_error(lambda: bitmask_n_words(-1), "n_bits = -1 must be rejected", ERR_INVALID_INPUT)


def test_bitmask_from_logical():
    bitmask = bitmask_from_logical(1, FLAGS_1011)
    assert bitmask.dtype == np.int32, f"expected int32 words, got {bitmask.dtype}"
    assert bitmask.shape == (1,), f"expected shape (1,), got {bitmask.shape}"
    assert not bitmask.flags.writeable, "a result is a value and must be read-only"
    assert bitmask[0] == WORD_1011, f"[T, T, F, T] must pack to {WORD_1011}, got {bitmask[0]}"

    # flag 33 is bit 0 of the second word
    flags = np.zeros(33, dtype=bool)
    flags[32] = True
    assert (bitmask_from_logical(2, flags) == [0, 1]).all(), "flag 33 must be bit 0 of word 2"

    assert_error(lambda: bitmask_from_logical(2, FLAGS_1011),
                 "n_words = 2 for 4 flags must be rejected", ERR_INVALID_INPUT)
    assert_error(lambda: bitmask_from_logical(0, FLAGS_1011),
                 "n_words = 0 for 4 flags must be rejected", ERR_INVALID_INPUT)
    assert_error(lambda: bitmask_from_logical(0, np.array([], dtype=bool)),
                 "no flags (n_bits = 0) must be rejected", ERR_INVALID_INPUT)


def test_bitmask_to_logical():
    flags = bitmask_to_logical(4, np.array([WORD_1011], dtype=np.int32))
    assert flags.dtype == np.bool_, f"expected bool flags, got {flags.dtype}"
    assert flags.shape == (4,), f"expected shape (4,), got {flags.shape}"
    assert not flags.flags.writeable, "a result is a value and must be read-only"
    assert (flags == FLAGS_1011).all(), f"{WORD_1011} must unpack to [T, T, F, T], got {flags}"

    assert_error(lambda: bitmask_to_logical(0, np.array([0], dtype=np.int32)),
                 "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
    assert_error(lambda: bitmask_to_logical(4, np.array([WORD_1011, 0], dtype=np.int32)),
                 "2 words for 4 flags must be rejected", ERR_INVALID_INPUT)
    assert_error(lambda: bitmask_to_logical(33, np.array([WORD_1011], dtype=np.int32)),
                 "1 word for 33 flags must be rejected", ERR_INVALID_INPUT)


def test_bitmask_test():
    bitmask = np.array([WORD_1011], dtype=np.int32)
    is_set = bitmask_test(4, bitmask, 4)
    assert isinstance(is_set, bool), f"expected a bool, got {type(is_set).__name__}"
    assert is_set, "flag 4 of 11 is set"
    assert not bitmask_test(4, bitmask, 3), "flag 3 of 11 is clear"

    assert_error(lambda: bitmask_test(4, bitmask, 0), "i_bit = 0 must be rejected", ERR_INVALID_INPUT)
    assert_error(lambda: bitmask_test(4, bitmask, 5), "i_bit = n_bits + 1 must be rejected",
                 ERR_INVALID_INPUT)
    assert_error(lambda: bitmask_test(0, bitmask, 1), "n_bits = 0 must be rejected", ERR_INVALID_INPUT)
    assert_error(lambda: bitmask_test(4, np.array([WORD_1011, 0], dtype=np.int32), 1),
                 "2 words for 4 flags must be rejected", ERR_INVALID_INPUT)


def test_bitmasks_from_logical_2D():
    # 33 flags per mask take 2 words, and every column starts on a word of its own: column 2's
    # flag 1 is bit 0 of its own first word, not bit 1 of column 1's second word
    flags = np.zeros((33, 2), dtype=bool, order="F")
    flags[:, 0] = True
    flags[0, 1] = True
    flags[32, 1] = True
    bitmasks = bitmasks_from_logical_2D(2, flags)
    assert bitmasks.dtype == np.int32, f"expected int32 words, got {bitmasks.dtype}"
    assert bitmasks.shape == (2, 2), f"expected shape (2, 2), got {bitmasks.shape}"
    assert bitmasks.flags.f_contiguous, "the masks must be column-major"
    assert not bitmasks.flags.writeable, "a result is a value and must be read-only"
    # column 1: 32 set flags are all bits of word 1 (-1 as int32), then flag 33 alone
    assert (bitmasks[:, 0] == [-1, 1]).all(), f"column 1 must be [-1, 1], got {bitmasks[:, 0]}"
    assert (bitmasks[:, 1] == [1, 1]).all(), f"column 2 must be [1, 1], got {bitmasks[:, 1]}"

    assert_error(lambda: bitmasks_from_logical_2D(1, flags),
                 "n_words = 1 for 33 flags must be rejected", ERR_INVALID_INPUT)
    assert_error(lambda: bitmasks_from_logical_2D(1, np.zeros((0, 2), dtype=bool, order="F")),
                 "no flags (n_bits = 0) must be rejected", ERR_INVALID_INPUT)


def test_bitmasks_to_logical_2D():
    # two masks of 4 flags: one word each
    bitmasks = np.array([[WORD_1011, 0]], dtype=np.int32, order="F")
    flags = bitmasks_to_logical_2D(4, bitmasks)
    assert flags.dtype == np.bool_, f"expected bool flags, got {flags.dtype}"
    assert flags.shape == (4, 2), f"expected shape (4, 2), got {flags.shape}"
    assert flags.flags.f_contiguous, "the flags must be column-major"
    assert not flags.flags.writeable, "a result is a value and must be read-only"
    assert (flags[:, 0] == FLAGS_1011).all(), f"column 1 must unpack to [T, T, F, T], got {flags[:, 0]}"
    assert not flags[:, 1].any(), "column 2 is the empty mask"

    # the round trip through the padded layout
    padded_flags = np.zeros((33, 2), dtype=bool, order="F")
    padded_flags[::3, 0] = True
    padded_flags[32, 1] = True
    round_trip = bitmasks_to_logical_2D(33, bitmasks_from_logical_2D(2, padded_flags))
    assert (round_trip == padded_flags).all(), "unpacking must undo packing"

    assert_error(lambda: bitmasks_to_logical_2D(0, bitmasks), "n_bits = 0 must be rejected",
                 ERR_INVALID_INPUT)
    assert_error(lambda: bitmasks_to_logical_2D(33, bitmasks),
                 "1 word per mask for 33 flags must be rejected", ERR_INVALID_INPUT)


def test_layout_matches_numpy_unpackbits():
    # the module documentation promises numpy decodes the words directly
    n_bits = 70
    flags = np.zeros(n_bits, dtype=bool)
    flags[[0, 5, 31, 32, 40, 63, 64, 69]] = True
    bitmask = bitmask_from_logical(bitmask_n_words(n_bits), flags)
    assert (decode_with_numpy(bitmask, n_bits) == bitmask_to_logical(n_bits, bitmask)).all(), \
        "np.unpackbits(bitorder='little') must decode what bitmask_to_logical decodes"
    assert (decode_with_numpy(bitmask, n_bits) == flags).all(), "np.unpackbits must recover the flags"


def test_word_with_only_bit_31_set():
    # flag 32 alone is the word 0x80000000, the most negative int32. In R this word is NA_integer_
    # (issue #206); in Python it is an ordinary value.
    flags = np.zeros(32, dtype=bool)
    flags[31] = True
    bitmask = bitmask_from_logical(1, flags)
    assert bitmask[0] == np.int32(-2**31), f"flag 32 alone must pack to -2**31, got {bitmask[0]}"
    assert bitmask_test(32, bitmask, 32), "flag 32 must test as set"
    assert not bitmask_test(32, bitmask, 31), "flag 31 must test as clear"
    assert (bitmask_to_logical(32, bitmask) == flags).all(), "the word must unpack to flag 32 alone"


if __name__ == "__main__":
    run_all_tests(globals().values())
