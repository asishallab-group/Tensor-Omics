from sys import argv as args, exit
import numpy as np
from tensor_omics import error_handling as _errors
from tensor_omics.error_handling import ARG_POS_FACTOR, ToxError

ANSI_START = "\033["
FG256_START = f"{ANSI_START}38;5;"

COLORS = {
    "green": 154,
    "copper": 214,
    "dark_copper": 208,
    "red": 196,
    "light_gray": 252,
    "yellow": 226,
    "cream": 255,
    "error": 222
}


def fg256(color_name=""):
    if color_name in COLORS:
        return f"{FG256_START}{COLORS[color_name]}m"
    else:
        return f"{ANSI_START}0m"


def cprint(s, **kwargs):
    colored = s
    for substr in s.split("@"):
        color_name = substr.split(".", 1)[0]
        colored = colored.replace(f"@{color_name}.", fg256(color_name))
    print(colored + fg256(), **kwargs)


def run_all_tests(functions, script_file_name=args[0], test_only=True):
    """Run all functions in `functions`. If test_only is True, run only those with name starting with `test_`, like `test_func`"""
    if test_only:
        all_tests = [func for func in functions if callable(func) and func.__name__.startswith("test_")]
    else:
        all_tests = list(functions)

    script_file_name = f"@light_gray.{script_file_name}"

    passed = 0
    failed = 0
    skipped = 0
    cprint(f"@cream.Running tests of '{script_file_name}@cream.'...")
    for test_func in all_tests:
        test_name = test_func.__name__
        try:
            test_func()
            cprint(f"@green.✓ @copper.{test_name} @green.passed@cream..")
            passed += 1
        except (AssertionError, ToxError) as e:
            cprint(f"@red.✗ @dark_copper.{test_name} @red.FAILED@cream.: @error.{e}")
            failed += 1
        # except Exception as e:
        #     # Some tests are expected to raise exceptions in certain cases
        #     if "Note:" in str(e) or "acceptable" in str(e).lower():
        #         print(f"~ {test_name} skipped (expected behavior): {e}")
        #         skipped += 1
        #     else:
        #         raise e
        #         tb = e.__traceback__.tb_next
        #         filename = tb.tb_frame.f_code.co_filename
        #         line = tb.tb_lineno
        #         print(f"✗ {test_name} FAILED with unexpected error: {e}, line {line}, file '{filename}'")
        #         failed += 1

    cprint(f"@cream.\nSummary: @green.{passed} passed@cream., @red.{failed} failed@cream., @yellow.{skipped} skipped")
    if (failed):
        exit(1)

    cprint(f"@cream.All tests in '{script_file_name}@cream.' passed successfully.")


#: every ERR_* the generated module exports, keyed by code, so a failure can name what arrived
_CODE_NAMES = {value: name for name, value in vars(_errors).items()
               if name.startswith("ERR_") and isinstance(value, int)}


def _code_label(code):
    name = _CODE_NAMES.get(code)
    return f"{name} ({code})" if name else f"code {code}"


def assert_error(func, msg, code=None):
    """Fail unless calling `func` raises the expected error.

    `code` is a bare tox_errors code -- ERR_INVALID_INPUT, not the `ierr` the Fortran packs
    the argument position into. check_err_code strips the position before it raises, so
    ToxError.code is already bare, and a packed value here is rejected as a mistake.

    Omit `code` to demand an error that carries no tox code at all: the binding's own
    ValueError or TypeError, raised before the library is ever called.
    """
    if code is not None and code >= ARG_POS_FACTOR:
        raise AssertionError(
            f"assert_error takes a bare error code, not a packed ierr: {code} means "
            f"{_code_label(code % ARG_POS_FACTOR)} at argument {code // ARG_POS_FACTOR}")
    try:
        func()
    # ToxError subclasses RuntimeError, so this clause has to come first or the code check
    # below is dead
    except ToxError as error:
        if code is None:
            raise AssertionError(
                f"{msg}: expected an error that is not a tensor-omics error, but got "
                f"{_code_label(error.code)}: {error}") from None
        if error.code != code:
            raise AssertionError(
                f"{msg}: expected {_code_label(code)}, but got "
                f"{_code_label(error.code)}: {error}") from None
    except (RuntimeError, ValueError, TypeError, FileNotFoundError) as error:
        if code is not None:
            raise AssertionError(
                f"{msg}: expected {_code_label(code)}, but got "
                f"{type(error).__name__}: {error}") from None
    else:
        raise AssertionError(f"{msg}: nothing was raised")
