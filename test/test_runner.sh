#!/bin/bash
# filepath: test/test_runner.sh

# This script assumes you have fpm installed and configured
# The project root is assumed to be the directory containing fpm.toml

# Navigate to the project root (if not already there)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# Ensure the test executable is built
echo "Building test executable..."
# Corrected fpm build command: target the specific test executable by name
# We specify --test followed by the name from fpm.toml's [test] section
fpm build --profile test --test test_runner_suite

# Check if build was successful
if [ $? -ne 0 ]; then
    echo "FPM build failed. Exiting."
    exit 1
fi

# Determine the path to the compiled test executable
# fpm typically puts test executables under build/fpm_test/<test_name>
TEST_EXEC="./build/fpm_test/test_runner_suite" # Name from fpm.toml [test] name

# Double-check if the executable exists
if [ ! -f "$TEST_EXEC" ]; then
    echo "Error: Compiled test executable not found at '$TEST_EXEC'."
    echo "Please check your fpm.toml [test] configuration and fpm build output."
    exit 1
fi

# Pass all arguments from the shell script to the compiled Fortran executable
echo "Running tests with arguments: $@"
"$TEST_EXEC" "$@"

# Check the exit status of the Fortran executable
if [ $? -ne 0 ]; then
    echo "Tests FAILED."
    exit 1
else
    echo "All tests PASSED."
    exit 0
fi