#!/bin/bash
gfortran -c -fPIC TensorOmics_Interface.f90
gfortran -c -fPIC TensorOmics_mod.f90
gfortran -shared -o libtensoromics.so TensorOmics_Interface.o TensorOmics_mod.o



