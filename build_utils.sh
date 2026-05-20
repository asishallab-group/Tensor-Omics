function init() {
  # handle_args overwrites:
  # FC if --fc=<compiler> specified
  handle_args "$@"

  # --compiler beats global $COMPILER beats --fc beats global $FC
  get_compiler

  if [[ -z $(command -v $COMPILER) ]]; then
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
  LD_LIBRARY_PATH="$libpath" $prefix --features "$COMPILER_FEATURE" --compiler $COMPILER --flag "$FLAGS $DIRECTIVES" --flag "-I." -- $ARGS
}

# gets compiler from context, it uses
# 1. $COMPILER if defined (by --compiler)
# 2. else $FC
# falls back to gfortran if the set compiler is not known
function get_compiler() {
  declare compiler=${COMPILER:-$FC}
  declare default=gfortran
  COMPILER_FEATURE=

  # Detect compiler and choose appropriate profile:
  if [[ "$compiler" == "ifx" ]]; then
    COMPILER_FEATURE=ifx
    COMPILER=ifx
  elif [[ "$compiler" == "nvfortran" ]]; then
    COMPILER_FEATURE=nvfortran
    COMPILER=nvfortran
  else
    if [[ $compiler ]]; then
      if [[ $compiler != "$default" ]]; then
        if [[ $I_WANT_TO_USE_THIS_COMPILER ]]; then
          COMPILER_FEATURE=unknown-compiler
          COMPILER="$compiler"
          return
        else
          stderr "Compiler '$compiler' not officially supported by Tensor Omics, trying '$default' instead. Use '--i-want-to-use-this-compiler' to run with '$compiler' anyway."
        fi
      fi
    else
      stderr "No compiler specified, using '$default'. To specify the compiler, use --compiler=<compiler> or the env variables \$FC, \$COMPILER"
    fi
    COMPILER_FEATURE=$default
    COMPILER=$default
  fi
}

function get_flags() {
  if [[ "$OVERRIDE_FLAGS" ]]; then
    echo "$OVERRIDE_FLAGS"
    return
  fi

  if [[ $MAX_PERFORMANCE ]]; then
    echo -en "-DMAX_PERFORMANCE "
    # Detect compiler and choose appropriate profile:
    if [[ "$COMPILER" == "ifx" ]]; then
      # echo "-O0 -g -traceback -check all -warn all -diag-enable=all -fPIC"
      echo "-O3 -warn all -diag-enable=all -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
    elif [[ "$COMPILER" == "nvfortran" ]]; then
      echo "-O3 -Mconcur -fPIC -fopenmp -stdpar=multicore"
    else
      echo "-O3 -march=native -mtune=native -fopenmp -funroll-loops -ftree-vectorize -fPIC"
    fi
  else
    echo "-fPIC"
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
      val="${undashed#"$key"}"      # strip leading $key
      val="${val#=}"                # strip leading '=' if present
      : ${val:=1}

      declare varname="$key"

      # Replace non-alphanumeric with _
      varname="${varname//[^a-zA-Z0-9]/_}"

      # Uppercase everything
      varname="${varname^^}"

      declare -g "$varname=$val"
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
