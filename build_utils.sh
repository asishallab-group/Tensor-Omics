function get_alignment() {
  ALIGN=32
  # Detect capabilities in order of descending priority:
  if lscpu | grep -q amx; then
  ALIGN=128
  elif lscpu | grep -q avx512; then
  ALIGN=64
  elif lscpu | grep -q avx2; then
  ALIGN=32
  elif lscpu | grep -q sse2; then
  ALIGN=16
  fi
  echo $ALIGN
}

function get_compiler() {
  # Detect compiler and choose appropriate profile:
  if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
    echo ifx
  elif [[ "$FC" == "nvfortran" ]]; then
    echo nvfortran
  else
    echo gfortran
  fi
}

function get_flags() {
  # Detect compiler and choose appropriate profile:
  if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
    echo "-O2 -fopenmp-target-do-concurrent -warn all -diag-enable=all -qopenmp -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
  elif [[ "$FC" == "nvfortran" ]]; then
    echo "-O2 -Mconcur -Mstack_arrays -fPIC -fopenmp -stdpar=multicore"
  else
    echo "-O2 -march=native -mtune=native -fopenmp -ffast-math -funroll-loops -ftree-vectorize -fassociative-math -fPIC"
  fi
}

function get_module_flag() {
  if [[ $1 ]]; then
    if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
      echo "-module $1"
    elif [[ "$FC" == "nvfortran" ]]; then
      echo "-module $1"
    else
      echo "-J$1"
    fi
  fi
}

function stderr() {
  echo "$@" >&2
}

function check_build() {
  if [[ "$@" ]]; then
    missing=()
    for c in "$@"; do
      if [[ ! $(find build -name "$c") ]]; then
        missing+=(" '$c'")
      fi
    done
    if [[ "$missing" ]]; then
      stderr "Missing files:$(IFS=', ';echo "${missing[*]}")"
      return 1
    fi
  fi

  mod_count=$(find build -name "*.mod" | wc -l)
  obj_count=$(find build -name "*.o" | wc -l) 
  so_count=$(find build -name "*.so" | wc -l)

  if [ $mod_count -eq 0 ] && [ ! $obj_count -eq 0 ] && [ ! $so_count -eq 0 ]; then
    stderr "Missing .mod files"
    return 2
  elif [ $mod_count -eq 0 ] && [ $obj_count -eq 0 ] && [ ! $so_count -eq 0 ]; then
    stderr "Missing .mod and .o files"
    return 3
  elif [ $mod_count -eq 0 ] && [ $obj_count -eq 0 ] && [ $so_count -eq 0 ]; then
    stderr "Missing .mod and .o and .so files"
    return 4
  elif [ $mod_count -eq 0 ] && [ ! $obj_count -eq 0 ] && [ $so_count -eq 0 ]; then
    stderr "Missing .mod and .so files"
    return 5
  elif [ ! $mod_count -eq 0 ] && [ $obj_count -eq 0 ] && [ ! $so_count -eq 0 ]; then
    stderr "Missing .o files"
    return 6
  elif [ ! $mod_count -eq 0 ] && [ $obj_count -eq 0 ] && [ $so_count -eq 0 ]; then
    stderr "Missing .o and .so files"
    return 7
  elif [ ! $mod_count -eq 0 ] && [ ! $obj_count -eq 0 ] && [ $so_count -eq 0 ]; then
    stderr "Missing .so files"
    return 8
  fi
}

function check_exit_code() {
  code=$?
  if [[ ! $code -eq 0 ]]; then
    if [[ "$@" ]]; then
      echo "$@ - Exit code: $code"
    fi
    exit $code
  fi
}