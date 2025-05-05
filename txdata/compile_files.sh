#!/bin/bash
gfortran -shared -fPIC -o build/libtensoromics.so src/TensorOmics_Interface.f90 src/TensorOmics_mod.f90
mv tensoromics_interface.mod src/tensoromics_interface.mod


