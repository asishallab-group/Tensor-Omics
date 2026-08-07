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
  # Libraries
  echo -en "-lzip -lxxhash "

  # Detect compiler and choose appropriate profile:
  if [[ "$FC" == "ifx" || "$FC" == "ifort" ]]; then
    # echo "-O0 -g -traceback -check all -warn all -diag-enable=all -fPIC"
    echo "-O2 -fopenmp-target-do-concurrent -warn all -diag-enable=all -qopenmp -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
  elif [[ "$FC" == "nvfortran" ]]; then
    echo "-O2 -Mconcur -fPIC -fopenmp -stdpar=multicore"
  else
    echo "-O2 -march=native -mtune=native -fopenmp -funroll-loops -ftree-vectorize -fPIC"
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

function handle_args() {
  MAX_PERF_FLAG=""
  ARGS=""
  DIRECTIVES=""
  for arg in "$@"; do
    if [[ "$arg" == "--max-performance" ]]; then
      MAX_PERF_FLAG="-DMAX_PERFORMANCE"
    elif [[ "$arg" == -D* ]]; then
      DIRECTIVES="$DIRECTIVES $arg"
    # genericly handle optional flags
    elif [[ "$arg" == --* ]]; then
      declare undashed=${arg:2}
      declare key=${undashed%=*}
      declare val=${undashed##$key}
      if [[ ! $val ]];then val=1;fi
      declare -g "$(echo "$key" | sed 's/\W/_/g; s/\w/\U&/g')=$val"
    else
      ARGS="$ARGS $arg"
    fi
  done

  # if extra directives are added, a clean build is necessary. Otherwise fpm doesn't recompile
  if [[ $DIRECTIVES ]]; then
    CLEAN_BUILD=1
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
      echo -en "$@ - "
    fi
    echo "Exit code: $code"
    exit $code
  fi
}

function generate_fpm_toml() {
  extra_libs=   # space, tab or comma separated list, like: "lib1, lib2"
  if [[ "$2" == "ifx" ]]; then
    extra_libs="iomp5"
  fi

  awk -v extra_libs="$extra_libs" '
{
  line = $0

  # match category, like "build" from [build] or "test.dependencies" from [test.dependencies]
  match($0, /^[ \t]*\[[ \t]*([a-z\.]+)[ \t]*\]/, arr)

  if (arr[1]) {
    category = arr[1]
  }

  if (category == "build") {
    # match: link = [ "lib_1.0" , "lib_2.0" ]
    # and extract the array elements
    match($0, /^[ \t]*link[ \t]*=[ \t]*\[([ \ta-z,",_0-9\.]+)\]/, arr)

    if (arr[1]) {
      # unify separators, trim start and wrap each lib in: ",\"<lib>\""
      gsub(/[\t,]/," ",extra_libs)
      sub(/^ +/,"",extra_libs)
      gsub(/[^ ]+/,",\"&\"",extra_libs)

      line = sprintf("link = [%s %s]", arr[1], extra_libs)
    }
  }

  print line
}
' $1
}