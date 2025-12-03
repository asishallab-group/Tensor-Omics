#!/bin/bash

gfortran -g -O0 -fcheck=all -fbacktrace \
    src/tox_errors.F90 \
    src/f42_utils.F90 \
    src/k_d_tree.F90 \
    src/tox_loess.F90 \
    test/debug_loess_gdb.f90 \
    -o debug_loess_gdb