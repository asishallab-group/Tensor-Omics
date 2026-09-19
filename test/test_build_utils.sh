#!/bin/bash
# Tests of build_utils.sh, run by run_all_tests.sh from the repository root. Exits non-zero if a
# case fails.
#
# remove_stale_libraries gets what `fpm build --list` prints -- every library, then every object
# -- and removes a library that is older than one of *its own* objects. The cases build a fake
# tree in a temporary directory and set its file times with `touch -d`: a root package "root"
# whose sources lie in src/ (objects src_*.o) and a dependency "dep" in the directory dep/
# (objects dep_src_*.o), as fpm names them for tensor_omics and its test_framework.

source build_utils.sh

scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT
failed=0

# build_tree <root library> <dependency library> <root object> <dependency object>: the modify
# times, as seconds past a fixed minute
build_tree() {
  rm -rf "$scratch/build"
  mkdir -p "$scratch/build"
  touch -d "2026-01-01 00:00:$1" "$scratch/build/libroot.so"
  touch -d "2026-01-01 00:00:$2" "$scratch/build/libdep.so"
  touch -d "2026-01-01 00:00:$3" "$scratch/build/src_a.F90.o"
  touch -d "2026-01-01 00:00:$4" "$scratch/build/dep_src_b.F90.o"
  fpm_list=" $scratch/build/libroot.so
 $scratch/build/libdep.so
 $scratch/build/src_a.F90.o
 $scratch/build/dep_src_b.F90.o"
}

# check <case> <library> <expected: kept|removed>
check() {
  declare state=kept
  [[ -e $scratch/build/$2 ]] || state=removed
  if [[ $state == "$3" ]]; then
    cecho "$1: $2 $state: ${COLOR_GREEN}success"
  else
    cecho "$1: $2 $state, expected $3: ${COLOR_RED}failed"
    failed=1
  fi
}

# an earlier build compiled root sources and failed before linking: the root library is stale
build_tree 10 10 20 05
remove_stale_libraries "$fpm_list" root 2>/dev/null
check "root library older than its own object" libroot.so removed
check "root library older than its own object" libdep.so kept

# the dependency is only older than a root object, not than its own: nothing of its is stale
build_tree 30 10 20 05
remove_stale_libraries "$fpm_list" root 2>/dev/null
check "dependency library older than a root object only" libroot.so kept
check "dependency library older than a root object only" libdep.so kept

# the dependency's own object is newer than it: the dependency library is stale
build_tree 30 10 20 15
remove_stale_libraries "$fpm_list" root 2>/dev/null
check "dependency library older than its own object" libroot.so kept
check "dependency library older than its own object" libdep.so removed

# every library newer than its objects: nothing to do
build_tree 30 30 20 20
remove_stale_libraries "$fpm_list" root 2>/dev/null
check "all libraries up to date" libroot.so kept
check "all libraries up to date" libdep.so kept

# the root package's name, as build.sh passes it, comes from fpm.toml
if [[ $(fpm_package_name) == tensor_omics ]]; then
  cecho "fpm_package_name reads fpm.toml: ${COLOR_GREEN}success"
else
  cecho "fpm_package_name reads fpm.toml: got '$(fpm_package_name 2>/dev/null)': ${COLOR_RED}failed"
  failed=1
fi

exit $failed
