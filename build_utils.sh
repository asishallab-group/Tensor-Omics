function init() {
  ALIGN=$(get_alignment)

  # handle_args overwrites:
  # ALIGN if --align=<align> specified
  # FC if --fc=<compiler> specified
  handle_args "$@"

  # --compiler beats global $COMPILER beats --fc beats global $FC
  COMPILER=$(get_compiler)
  if [[ -z $(which $COMPILER) ]]; then
    stderr "$COMPILER not installed"
    exit 1
  fi
  FLAGS=$(get_flags)
}

function utils_fpm() {
  declare prefix="fpm build"
  declare libpath="$LD_LIBRARY_PATH"
  if [[ "$1" == "test" ]]; then
    prefix="fpm test --target ${2:-run_tests}"
    libpath=build:"$libpath"
  elif [[ "$1" == "list" ]]; then
    prefix="fpm build --list"
  fi
  LD_LIBRARY_PATH="$libpath" $prefix --compiler $COMPILER --flag "$FLAGS $DIRECTIVES" --flag "-DDEFAULT_ALIGNMENT=$ALIGN" --flag "-I." -- $ARGS
}

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

# gets compiler from context, it uses
# 1. $COMPILER if defined (by --compiler)
# 2. else $FC
# falls back to gfortran if the set compiler is not known
function get_compiler() {
  declare compiler=${COMPILER:-$FC}
  declare default=gfortran

  # Detect compiler and choose appropriate profile:
  if [[ "$compiler" == "ifx" || "$compiler" == "ifort" ]]; then
    echo ifx
  elif [[ "$compiler" == "nvfortran" ]]; then
    echo nvfortran
  else
    if [[ $compiler ]]; then
      if [[ $compiler != "$default" ]]; then
        stderr "Compiler '$compiler' not officially supported by Tensor Omics, trying '$default' instead"
      fi
    else
      stderr "No compiler specified, using '$default'. To specify the compiler, use --compiler=<compiler> or the env variables \$FC, \$COMPILER"
    fi
    echo $default
  fi
}

function get_flags() {
  # Libraries: greps the libraries from .fpm.toml and translates them from '"<lib>"' to '-l<lib>'
  printf "%s" "$(grep -oP 'link = \[\K.*\]' .fpm.toml)" | sed 's/ //g; s/"/-l/g; s/-l,/ /g; s/-l]/ /g;'

  if [[ $MAX_PERFORMANCE ]]; then
    echo -en "-DMAX_PERFORMANCE "
  fi

  # Detect compiler and choose appropriate profile:
  if [[ "$COMPILER" == "ifx" || "$COMPILER" == "ifort" ]]; then
    # echo "-O2 -parallel -warn all -diag-enable=all -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
    echo "-O0 -g -traceback -check all -warn all -diag-enable=all -fPIC"
  elif [[ "$COMPILER" == "nvfortran" ]]; then
    echo "-O2 -Mconcur -fPIC -fopenmp -stdpar=multicore"
  else
    echo "-O0 -g -fbacktrace -fcheck=all -fPIC"
  fi
}

function handle_args() {
  ARGS=""
  DIRECTIVES=""
  
  for arg in "$@"; do
    if [[ "$arg" == -D* ]]; then
      DIRECTIVES="$DIRECTIVES $arg"
    # genericly handle optional flags
    elif [[ "$arg" == --* ]]; then
      declare undashed=${arg:2}
      declare key=${undashed%%=*}
      # extract value after first '=' if present, else set to 1
      declare val=$(echo "$undashed" | sed 's/^'$key'\(=\(.*\)\?\)\?/\2/g')
      : ${val:=1}
      declare -g "$(echo "$key" | sed 's/\W/_/g; s/\w/\U&/g')=$val"
    else
      ARGS="$ARGS $arg"
    fi
  done
}

function stderr() {
  echo "$@" >&2
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