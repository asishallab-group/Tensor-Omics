#!/bin/bash
mkdir -p build
# This command uses the -shared flag to create libF42CsvReader.so
gfortran -O3 -march=native -ffast-math -funroll-loops -fPIC \
  -shared src/*.f90 -o build/libF42CsvReader.so
echo "Shared library build/libF42CsvReader.so created successfully."