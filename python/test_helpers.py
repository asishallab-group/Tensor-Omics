from sys import argv as args, exit
import numpy as np

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
            cprint(f"@green.✓ @copper.{test_name} @green.passed.")
            passed += 1
        except AssertionError as e:
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


def assert_error(func, msg):
    try:
        func()
        assert False, msg
    except AssertionError as e:
        raise e
    except (RuntimeError, ValueError, FileNotFoundError) as e:
        pass
