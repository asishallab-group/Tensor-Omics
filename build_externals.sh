#!/bin/bash
# build_externals.sh
source build_utils.sh
COMPILER=$(get_compiler)
mkdir -p external/lib

FFLAGS_F77="-O2 -std=legacy -ffixed-line-length-none -fallow-argument-mismatch -fPIC -w"
SRC_F="external/loess/loessf.f external/loess/supp.f external/loess/d1mach.f external/loess/dqrsl.f external/loess/dsvdc.f external/loess/daxpy.f external/loess/ddot.f external/loess/dnrm2.f external/loess/dscal.f external/loess/dswap.f external/loess/dcopy.f external/loess/drot.f"

echo "Compiling LOESS kernels..."
for f in $SRC_F; do
    $COMPILER $FFLAGS_F77 -c $f -o external/lib/$(basename "$f" .f).o
done

$COMPILER $(get_flags) -c external/loess/drotg.f90 -o external/lib/drotg.o

ar rcs external/lib/libloess_netlib.a external/lib/*.o