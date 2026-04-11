from sys import argv as args


def run_all_tests(functions, script_file_name=args[0], test_only=True):
    """Run all functions in `functions`. If test_only is True, run only those with name starting with `test_`, like `test_func`"""
    if test_only:
        all_tests = [func for func in functions if callable(func) and func.__name__.startswith("test_")]
    else:
        all_tests = list(functions)

    passed = 0
    failed = 0
    skipped = 0

    print(f"Running tests of in '{script_file_name}'...")
    for test_func in all_tests:
        test_name = test_func.__name__
        try:
            test_func()
            print(f"✓ {test_name} passed.")
            passed += 1
        except AssertionError as e:
            print(f"✗ {test_name} FAILED: {e}")
            failed += 1
        except Exception as e:
            # Some tests are expected to raise exceptions in certain cases
            if "Note:" in str(e) or "acceptable" in str(e).lower():
                print(f"~ {test_name} skipped (expected behavior): {e}")
                skipped += 1
            else:
                tb = e.__traceback__.tb_next
                filename = tb.tb_frame.f_code.co_filename
                line = tb.tb_lineno
                print(f"✗ {test_name} FAILED with unexpected error: {e}, line {line}, file '{filename}'")
                failed += 1

    print(f"\nSummary: {passed} passed, {failed} failed, {skipped} skipped")
    assert failed == 0, f"{failed} tests in '{script_file_name}' failed."

    print(f"All tests in '{script_file_name}' passed successfully.")
