#!/bin/bash

gfortran -g -O0 -fcheck=all -fbacktrace -lxxhash \
    src/tox_errors.F90 \
    src/config.F90 \
    src/tox_conversions.F90 \
    src/tox_array_read_write/array_utils.F90 \
    src/safeguard.F90 \
    src/f42_utils.F90 \
    src/k_d_tree.F90 \
    src/tox_data/xxh3_hashmap_module.F90 \
    src/tox_data/tox_data_tools.F90 \
    src/tox_loess.F90 \
    test/debug_loess_gdb.f90 \
    -o debug_loess_gdb