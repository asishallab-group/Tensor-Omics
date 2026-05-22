COLOR_GREEN="\033[38;5;154m"
COLOR_COPPER="\033[38;5;214m"
COLOR_DARK_COPPER="\033[38;5;208m"
COLOR_RED="\033[38;5;196m"
COLOR_LIGHT_GRAY="\033[38;5;252m"
COLOR_YELLOW="\033[38;5;226m"
COLOR_CREAM="\033[38;5;255m"
COLOR_ERROR="\033[38;5;222m"
COLOR_RESET="\033[0m"

function init() {
  handle_args "$@"

  # --compiler beats global $TOX_COMPILER beats global $FC
  get_compiler

  if [[ -z $(command -v $COMPILER) ]]; then
    error "'$(echo_compiler $COMPILER)' not installed or accessible in current scope"
    exit 1
  fi
  FLAGS=$(get_flags)
}

function utils_fpm() {
  cecho "${COLOR_CREAM}Using compiler: $(echo_compiler $COMPILER)"
  declare prefix="fpm build"
  declare libpath="$LD_LIBRARY_PATH"
  if [[ "$1" == "test" ]]; then
    prefix="fpm test --target ${2:-run_tests}"
    libpath=build:"$libpath"
  elif [[ "$1" == "list" ]]; then
    prefix="fpm build --list"
  fi
  LD_LIBRARY_PATH="$libpath" $prefix --features "$COMPILER_FEATURE" --compiler "$COMPILER" --flag "$FLAGS $DIRECTIVES" --link-flag "-Lexternal" --flag "-I." -- $ARGS
  rm -f build/cache.toml  # can cause issues (when switching branches and external libs are missing), but doesn't affect compilation when missing
}

# gets compiler from context, it uses
# 1. $TOX_COMPILER if defined (by --compiler)
# 2. else $FC
# falls back to gfortran if the set compiler is not known
function get_compiler() {
  declare compiler=${TOX_COMPILER:-$FC}
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
        if [[ $TOX_I_WANT_TO_USE_THIS_COMPILER ]]; then
          COMPILER_FEATURE=unknown-compiler
          COMPILER="$compiler"
          return
        else
          warning "Compiler '$(echo_compiler $compiler)' not officially supported by Tensor Omics, trying '$(echo_compiler $default)' instead.
Use '$COLOR_LIGHT_GRAY--i-want-to-use-this-compiler$COLOR_CREAM' to run with '$(echo_compiler $compiler)' anyway.
Use '$COLOR_LIGHT_GRAY--override-flags$COLOR_CREAM' to define additional compiler-related flags like '$COLOR_LIGHT_GRAY--override-flags=\"-O3 -fPIC\"$COLOR_RESET'
"
        fi
      fi
    else
      warning "No compiler specified, using '$(echo_compiler $default)'. To specify the compiler, use $COLOR_LIGHT_GRAY--compiler=<$(echo_compiler compiler)$COLOR_LIGHT_GRAY>$COLOR_CREAM or the env variables $COLOR_LIGHT_GRAY\$$(echo_compiler FC)$COLOR_CREAM, $COLOR_LIGHT_GRAY\$$(echo_compiler TOX_COMPILER)"
    fi
    COMPILER_FEATURE=$default
    COMPILER=$default
  fi
}

function get_flags() {
  if [[ "$TOX_OVERRIDE_FLAGS" ]]; then
    echo "$TOX_OVERRIDE_FLAGS"
    return
  fi

  if [[ $TOX_MAX_PERFORMANCE ]]; then
    echo -en "-DMAX_PERFORMANCE "
    # Detect compiler and choose appropriate profile:
    if [[ "$COMPILER" == "ifx" ]]; then
      # echo "-O0 -g -traceback -check all -warn all -diag-enable=all -fPIC"
      echo "-O3 -warn all -diag-enable=all -xHost -align array64byte -qopt-zmm-usage=high -qopt-prefetch=3 -qopt-matmul -fPIC"
    elif [[ "$COMPILER" == "nvfortran" ]]; then
      echo "-O3 -Mconcur -fPIC -fopenmp -stdpar=multicore"
    elif [[ "$COMPILER" == "gfortran" ]]; then
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

      declare -g "TOX_${varname}=$val"
    else
      ARGS="$ARGS $arg"
    fi
  done
}

function find_and_mv_libs() {
  while IFS= read -r line; do
    if [[ $line == *.so || $line == *.a ]]; then
      # remove leading whitespaces
      lib="${line#"${line%%[![:space:]]*}"}"
      cp "${lib}" "$2" 2>/dev/null
    fi
  done <<< "$1"
}

function echo_compiler() {
  echo "$COLOR_COPPER$1$COLOR_CREAM"
}

function cecho() {
  echo -e "$COLOR_CREAM$@$COLOR_RESET"
}

function error() {
  stderr "${COLOR_RED}Error$COLOR_CREAM: $@"
  exit 1
}

function warning() {
  stderr "${COLOR_DARK_COPPER}Warning$COLOR_CREAM: $@"
}

function stderr() {
  cecho "$@" >&2
}

function check_exit_code() {
  code=$?
  if [[ ! $code -eq 0 ]]; then
    if [[ "$@" ]]; then
      error "$@ - ${COLOR_RED}Exit code${COLOR_CREAM}: $code"
    else
      error "${COLOR_RED}Exit code${COLOR_CREAM}: $code"
    fi
  fi
}
