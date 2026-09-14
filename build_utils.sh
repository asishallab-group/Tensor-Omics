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
  # C_ONLY_FLAGS is what the C sources get and the Fortran ones do not. $DIRECTIVES goes to
  # both, and fpm keys its build directories by the Fortran flags, so a directive never needs
  # a clean build. A C-only flag does: the library sits in the Fortran-keyed directory, so
  # switching one back reuses the objects without relinking -- check_build_state hashes these.
  # That is also why the R fallback adds NO_R_BINDING to $DIRECTIVES, not to the C flags alone.
  C_ONLY_FLAGS="-fPIC"
  C_COMPILER=  # empty: fpm picks the C compiler that goes with $COMPILER
  if [[ "$DIRECTIVES" == *NO_R_BINDING* || "$DIRECTIVES" == *NO_C_BINDING* ]]; then
    :
  elif [[ -z $(command -v R) ]]; then
    warning "'$(echo_compiler R)' not found -- building without the R binding.
Install R to include it, or pass '$COLOR_LIGHT_GRAY--directive=NO_R_BINDING$COLOR_CREAM' to silence this."
    DIRECTIVES="$DIRECTIVES -DNO_R_BINDING"
  else
    C_ONLY_FLAGS="$C_ONLY_FLAGS $(R CMD config --cppflags)"
    if [[ "$COMPILER" == "nvfortran" ]]; then
      # R configures its headers for the C compiler R was built with -- Rconfig.h decides, for
      # one, whether R_ext/Boolean.h declares an enum with a fixed base type, which nvc cannot
      # parse. So the R shims, the only C sources, are compiled with R's own compiler.
      declare -a r_cc=($(R CMD config CC))
      C_COMPILER="${r_cc[0]}"
      C_ONLY_FLAGS="${r_cc[*]:1} $C_ONLY_FLAGS"
    fi
  fi
  C_FLAGS="$C_ONLY_FLAGS $DIRECTIVES"
}

# fpm decides what to recompile from each source's own content, and keys its build directories
# by the compiler's name and the Fortran flags. Whatever else changes what the library should
# contain is invisible to it, so it is hashed here and answered with a clean build:
#   - the hand-written headers: fpm never hashes what a source #includes (fortran-lang/fpm#358).
#     tox_marshal.h is not among them -- the generator stamps its hash into every shim instead.
#   - fpm.toml and the link flags: a changed link library alone does not relink the library.
#   - the C-only flags: see get_c_flags.
#   - the Fortran compiler's version: an upgrade keeps the name, and the old .mod files with it.
# This replaces the clean build on every branch switch, which was a stand-in for the same
# thing, and it also catches these changes when no branch changed -- a pull, a stash pop, an
# edit to macros.h. One marker per compiler, as a clean build deletes only that compiler's.
function check_build_state() {
  declare state
  state=$(
    {
      "$COMPILER" --version 2>/dev/null | head -1
      echo "C_ONLY_FLAGS=$C_ONLY_FLAGS"
      echo "C_COMPILER=$C_COMPILER"
      echo "LINK_FLAGS=$LINK_FLAGS"
      cat fpm.toml
      # git lists the headers without walking the tree by hand (the repo may sit under /mnt/c)
      { git ls-files --cached --others --exclude-standard -- '*.h' ':!src/generated/**' 2>/dev/null \
          || find src -name '*.h' -not -path 'src/generated/*'; } | LC_ALL=C sort -u | xargs -r -d '\n' sha256sum
    } | sha256sum | cut -c1-16
  )
  declare marker="build/.${COMPILER}.${state}.buildstate"
  if [[ ! -f "$marker" ]]; then
    TOX_CLEAN_BUILD=1
    rm -f "build/.${COMPILER}."*.buildstate "build/.${COMPILER}."*.branch build/.branch
    : > "$marker"
  fi
}

# Regenerates the bindings and the generated Fortran wrappers from `src/` before compiling, so a
# source change and its generated layers can never drift apart in a build. Everything under
# `src/generated` is committed as well, so a machine without the generator's dependencies still
# builds -- it just builds what is in the tree, which is what the warning below says.
# `--skip-code-generation` skips the stage outright.
function generate_code() {
  if [[ "$TOX_SKIP_CODE_GENERATION" ]]; then
    return
  fi

  declare python=$(command -v python3 || command -v python)
  declare hint="pass '$COLOR_LIGHT_GRAY--skip-code-generation$COLOR_CREAM' to silence this."
  if [[ -z "$python" ]]; then
    warning "'$(echo_compiler python3)' not found -- building the generated sources as they stand in the tree.
Install Python to regenerate them, or $hint"
    return
  fi
  if ! "$python" -c "import ford" >/dev/null 2>&1; then
    warning "'$(echo_compiler ford)' not found -- building the generated sources as they stand in the tree.
Run '$COLOR_LIGHT_GRAY$(basename $python) -m pip install ford$COLOR_CREAM' to regenerate them, or $hint"
    return
  fi

  cecho "${COLOR_CREAM}Generating the bindings from $(echo_compiler src/)"
  "$python" helper/generate_code.py
  check_exit_code "Code generation failed"
}

function utils_fpm() {
  cecho "${COLOR_CREAM}Using compiler: $(echo_compiler $COMPILER)"
  # --tests: fpm rebuilds a target's dependents only within the run that rebuilds the target,
  # and records nothing for the next run. A plain `fpm build` leaves the tests out of the model,
  # so a changed module (a parameter value, an interface) recompiled the library but left the
  # test objects compiled against the old .mod -- and the `fpm test` that followed found both
  # up to date. Building them here puts them in the same run, and only the stale ones rebuild.
  declare -a prefix=(fpm build --tests)
  declare libpath="$LD_LIBRARY_PATH"
  if [[ "$1" == "test" ]]; then
    prefix=(fpm test --target "${2:-run_tests}")
    libpath=build:"$libpath"
    if [[ $TOX_DEBUG ]]; then
      prefix+=(--runner "gdb --args")
    fi
  elif [[ "$1" == "list" ]]; then
    prefix=(fpm build --list)
  fi
  LD_LIBRARY_PATH="$libpath" "${prefix[@]}" --features "$FEATURES" --compiler "$COMPILER" --flag "$FLAGS $DIRECTIVES" --c-flag "$C_FLAGS" ${C_COMPILER:+--c-compiler "$C_COMPILER"} --link-flag "$LINK_FLAGS" -- $ARGS
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
  # -Lexternal is where build.sh puts the loess archives it builds, so no override removes it
  LINK_FLAGS="-Lexternal"
  if [[ "$COMPILER" == "nvfortran" ]]; then
    # nvfortran links through whichever `ld` comes first on PATH, and a linker searches only its
    # own default directories. When that is not the system linker nvfortran was configured
    # against -- a Homebrew binutils, say -- it finds neither libc nor libgcc_s. The directories
    # are read from nvfortran's own configuration, so this is right on any machine, and a no-op
    # where the system linker is first anyway. (NVCOMPILER_LINKER is NVIDIA's own switch for the
    # linker itself, for anyone who prefers to fix it in the environment.)
    LINK_FLAGS="$LINK_FLAGS$(nvfortran -show 2>/dev/null |
      awk -F= '/^(DEFSTDOBJDIR|GCCDIR) /{gsub(/[ \t]+$/, "", $2); printf " -L%s", $2}')"
  fi
  if [[ "$TOX_OVERRIDE_LINK_FLAGS" ]]; then
    # The compiler's own feature -- what get_compiler put in $FEATURES -- carries nothing but
    # its link libraries, so replacing those means leaving it out. They cannot be passed via
    # --override-flags instead: fpm never hands --flag to the link of the shared library.
    FEATURES=
    LINK_FLAGS="$LINK_FLAGS $TOX_OVERRIDE_LINK_FLAGS"
  fi
  if [[ "$TOX_OVERRIDE_FLAGS" ]]; then
    FLAGS="$TOX_OVERRIDE_FLAGS"
    FEATURES=
    return
  fi

  if [[ $TOX_MAX_PERFORMANCE ]]; then
    FEATURES="$FEATURES,optimization"
    FLAGS="-O3"
    if [[ $TOX_DEBUG ]]; then
      declare ans=y
      if [[ -z $TOX_YES ]]; then
        read -n 1 -p "--debug doesn't work well with --max-performance, continue anyway? (y/n) " ans
        echo
      fi
      if [[ $ans == y ]]; then
        cecho "${COLOR_DARK_COPPER}Continuing"
        export TOX_YES=1  # test_runner.sh re-invokes build.sh as a subprocess (possibly several times, e.g. the kinds test); avoid re-prompting there
      else
        cecho "${COLOR_RED}Cancelled"
        exit 1
      fi
    fi
  else
    # Stated for the default build too: once `--features` is passed, fpm adds none of its own
    # profile flags, so the level used to be whatever each compiler assumes without one --
    # -O0 for gfortran, -O2 for ifx -- and the two default builds were quietly different.
    FLAGS="-O0"
  fi
  if [[ $TOX_DIAGNOSTICS || $TOX_DEBUG ]]; then
    FEATURES="$FEATURES,diagnostics"
  fi
  FEATURES="$FEATURES,default"
  FEATURES="${FEATURES#,}"  # the compiler's feature is absent under --override-link-flags
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
