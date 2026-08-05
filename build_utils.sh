DIRECTIVES=

# define colors if output is not being piped
if [[ -t 1 && -t 2 ]]; then
  COLOR_GREEN="\033[38;5;154m"
  COLOR_COPPER="\033[38;5;214m"
  COLOR_DARK_COPPER="\033[38;5;208m"
  COLOR_RED="\033[38;5;196m"
  COLOR_LIGHT_GRAY="\033[38;5;252m"
  COLOR_YELLOW="\033[38;5;226m"
  COLOR_CREAM="\033[38;5;255m"
  COLOR_ERROR="\033[38;5;222m"
  COLOR_RESET="\033[0m"
else
  DIRECTIVES="-DNO_COLORS"
fi

function init() {
  handle_args "$@"
  # --compiler beats global $TOX_COMPILER beats global $FC
  get_compiler

  if [[ -z $(command -v $COMPILER) ]]; then
    error "'$(echo_compiler $COMPILER)' not installed or accessible in current scope"
    exit 1
  fi
  get_flags_and_features
  get_c_flags
}

# Flags for the C sources (the R .Call shims in src/generated/bindings/r). fpm's --flag is Fortran
# only, so C needs its own. The shims are guarded by NO_R_BINDING / NO_C_BINDING, so when
# either is set they compile to empty objects that need no R headers; otherwise they need R's
# include path. If the R layer is wanted but R is not installed, drop it with a warning.
function get_c_flags() {
  C_FLAGS="-fPIC $DIRECTIVES"
  if [[ "$DIRECTIVES" == *NO_R_BINDING* || "$DIRECTIVES" == *NO_C_BINDING* ]]; then
    return
  fi
  if [[ -z $(command -v R) ]]; then
    warning "'$(echo_compiler R)' not found -- building without the R binding.
Install R to include it, or pass '$COLOR_LIGHT_GRAY--directive=NO_R_BINDING$COLOR_CREAM' to silence this."
    DIRECTIVES="$DIRECTIVES -DNO_R_BINDING"
    C_FLAGS="$C_FLAGS -DNO_R_BINDING"
    TOX_CLEAN_BUILD=1
    return
  fi
  C_FLAGS="$C_FLAGS $(R CMD config --cppflags)"
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
  LD_LIBRARY_PATH="$libpath" $prefix --features "$FEATURES" --compiler "$COMPILER" --flag "$FLAGS $DIRECTIVES" --c-flag "$C_FLAGS" --link-flag "-Lexternal" --flag "-I." -- $ARGS
  exit_code=$?
  rm -f build/cache.toml  # can cause issues (when switching branches and external libs are missing), but doesn't affect compilation when missing
  (exit $exit_code)
}

# gets compiler from context, it uses
# 1. $TOX_COMPILER if defined (by --compiler)
# 2. else $FC
# falls back to gfortran if the set compiler is not known
function get_compiler() {
  declare compiler=${TOX_COMPILER:-$FC}
  declare default=gfortran
  FEATURES=

  # Detect compiler and choose appropriate profile:
  if [[ "$compiler" == "ifx" ]]; then
    FEATURES=ifx
    COMPILER=ifx
  elif [[ "$compiler" == "nvfortran" ]]; then
    FEATURES=nvfortran
    COMPILER=nvfortran
  else
    if [[ $compiler ]]; then
      if [[ $compiler != "$default" ]]; then
        if [[ $TOX_I_WANT_TO_USE_THIS_COMPILER ]]; then
          FEATURES=unknown-compiler
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
    FEATURES=$default
    COMPILER=$default
  fi
}

function get_flags_and_features() {
  if [[ "$TOX_OVERRIDE_FLAGS" ]]; then
    FLAGS="$TOX_OVERRIDE_FLAGS"
    FEATURES=
    return
  fi

  if [[ $TOX_MAX_PERFORMANCE ]]; then
    FEATURES="$FEATURES,optimization"
    FLAGS="-DMAX_PERFORMANCE -O3"
  elif [[ $TOX_DIAGNOSTICS ]]; then
    FLAGS="-O0"
  fi
  if [[ $TOX_DIAGNOSTICS ]]; then
    FEATURES="$FEATURES,diagnostics"
  fi
  FEATURES="$FEATURES,default"
}

function handle_args() {
  ARGS=""
  
  for arg in "$@"; do
    if [[ "$arg" == --* ]]; then
      declare undashed=${arg:2}
      declare key=${undashed%%=*}

      # extract value after first '=' if present, else set to 1
      val="${undashed#"$key"}"      # strip leading $key
      if [[ "$val" == *=* ]]; then
        val="${val#=}"                # strip leading '=' if present
      else
        val=1
      fi

      declare varname="$key"

      # Replace non-alphanumeric with _
      varname="${varname//[^a-zA-Z0-9]/_}"

      # Uppercase everything
      varname="TOX_${varname^^}"

      if [[ "$varname" == "TOX_DIRECTIVE" ]]; then
        DIRECTIVES="$DIRECTIVES -D${val}"
        TOX_CLEAN_BUILD=1
      else
        declare -g "${varname}=$val"
      fi
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
