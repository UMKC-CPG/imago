#!/usr/bin/env bash
# envs.sh -- switch the active Imago toolchain between build "flavors".
#
# This file is meant to be *sourced* (not executed), normally from your
# shell startup right after sourcing the Imago rc file, e.g.:
#
#     source "$IMAGO_DIR/.imago/imagorc"
#     source "$IMAGO_BIN/envs.sh"
#
# WHAT A FLAVOR IS
# ---------------
# The production install lives at $IMAGO_DIR/bin (the normal
# `make install` target).  A *flavor* is an alternative build of the
# Fortran engine -- e.g. an AddressSanitizer build for leak hunting, or
# an audit build with runtime checks and traps -- kept in its own
# directory under $IMAGO_DIR/envs/<flavor>/bin.  These exist to support
# the bug-squashing campaign (see dev/DEBUG.md) and, later, profiling
# and parallel-debug builds.
#
# Only the compiled Fortran executables differ between flavors; the
# Python helper scripts and the share/ database are compiler-agnostic
# and identical.  So a flavor's bin is just a set of symlinks to the
# production bin with the engine executables (imago, imagoG, ...)
# overlaid by the flavor build.  Activating a flavor therefore only
# repoints IMAGO_BIN (and PATH / PYTHONPATH) at that bin; the shared
# share/ and .imago/ are reused as-is.  No data is duplicated, and the
# production install is never modified.
#
# USAGE
# -----
#     imago_env [flavor]      activate a flavor; no arg (or "release"/
#                             "prod") restores the production install
#     imago_env --list        list the flavors that are built
#     imago_env --build NAME   build flavor NAME from preset gfortran-NAME
#                              and (re)assemble its bin
#     imago_env --help        show this usage
#
# Examples:
#     imago_env asan          # run the engine under AddressSanitizer
#     imago_env audit         # run with checks + FP traps + warnings
#     imago_env               # back to the optimized production build


# Internal helper: echo $1 (a colon-separated PATH-like list) with every
#   occurrence of the directory $2 removed.  Used to strip a previously
#   activated flavor bin before prepending the new one, so repeated
#   switches never accumulate stale entries.
_imago_path_drop() {
    local list="$1" want="$2" out="" entry
    local IFS=:
    for entry in $list; do
        [ "$entry" = "$want" ] && continue
        out="${out:+$out:}$entry"
    done
    printf '%s' "$out"
}


imago_env() {
    if [ -z "$IMAGO_DIR" ]; then
        echo "imago_env: IMAGO_DIR is not set; source the Imago rc" \
             "file first." >&2
        return 1
    fi

    case "$1" in
        --help|-h)
            # Print the usage block from the header of this file.
            sed -n '/^# USAGE/,/^$/p' \
                "${BASH_SOURCE[0]:-$IMAGO_BIN/envs.sh}" \
                | sed 's/^# \{0,1\}//'
            return 0
            ;;
        --list|-l)
            echo "Production: $IMAGO_DIR/bin"
            if [ -d "$IMAGO_DIR/envs" ]; then
                local d
                for d in "$IMAGO_DIR"/envs/*/bin; do
                    [ -d "$d" ] && echo "Flavor:    $d"
                done
            fi
            return 0
            ;;
        --build|-b)
            _imago_env_build "$2"
            return $?
            ;;
    esac

    # Resolve the requested flavor to a bin directory.  An empty flavor
    #   or "release"/"prod"/"production" means the production install.
    local flavor="${1:-release}"
    local target_bin
    case "$flavor" in
        release|prod|production)
            flavor="release"
            target_bin="$IMAGO_DIR/bin"
            ;;
        *)
            target_bin="$IMAGO_DIR/envs/$flavor/bin"
            ;;
    esac

    if [ ! -d "$target_bin" ]; then
        echo "imago_env: flavor '$flavor' is not built ($target_bin)." >&2
        echo "           build it with:  imago_env --build $flavor" >&2
        return 1
    fi

    # Strip the previously activated Imago bin (production or a prior
    #   flavor) from PATH / PYTHONPATH, then prepend the new one.
    #   IMAGO_BIN is the single source of truth imago.py uses to find
    #   the engine and the helper scripts.
    local previous="${_IMAGO_ACTIVE_BIN:-$IMAGO_DIR/bin}"
    PATH="$(_imago_path_drop "$PATH" "$previous")"
    PYTHONPATH="$(_imago_path_drop "$PYTHONPATH" "$previous")"

    export IMAGO_BIN="$target_bin"
    export PATH="$target_bin${PATH:+:$PATH}"
    export PYTHONPATH="$target_bin${PYTHONPATH:+:$PYTHONPATH}"
    export _IMAGO_ACTIVE_BIN="$target_bin"

    # AddressSanitizer flavors want LeakSanitizer at exit; set a sane
    #   default the user can still override beforehand.
    if [ "$flavor" != "${flavor%asan*}" ] && [ -z "$ASAN_OPTIONS" ]; then
        export ASAN_OPTIONS="detect_leaks=1:abort_on_error=1"
    fi

    echo "imago_env: active flavor = $flavor"
    echo "           IMAGO_BIN     = $IMAGO_BIN"
}


# Build flavor $1 from the matching CMake preset (gfortran-$1) and
#   (re)assemble its bin as symlinks to the production toolchain with
#   the freshly built engine executables overlaid.  A full preset build
#   is run so every engine executable (not just imago/imagoG) is the
#   instrumented one.
_imago_env_build() {
    local flavor="$1"
    if [ -z "$flavor" ]; then
        echo "imago_env --build: a flavor name is required" \
             "(e.g. asan, audit)." >&2
        return 1
    fi
    local preset="gfortran-$flavor"
    local tree="$IMAGO_DIR/build/$preset"

    echo "imago_env: building preset '$preset' ..."
    cmake --preset "$preset" || return 1
    cmake --build --preset "$preset" -j || return 1

    local envbin="$IMAGO_DIR/envs/$flavor/bin"
    echo "imago_env: assembling $envbin ..."
    rm -rf "$envbin"
    mkdir -p "$envbin"
    # Symlink the whole production toolchain (scripts + helper exes) so
    #   every helper resolves and tracks production automatically.
    ln -s "$IMAGO_DIR"/bin/* "$envbin"/ 2>/dev/null

    # Overlay every freshly built executable whose name matches a
    #   production binary (the Fortran engine and its auxiliaries).
    local exe name
    while IFS= read -r exe; do
        name="$(basename "$exe")"
        if [ -e "$IMAGO_DIR/bin/$name" ]; then
            ln -sf "$exe" "$envbin/$name"
        fi
    done < <(find "$tree" -type f -executable \
                  ! -name '*.so' ! -name '*.so.*' ! -name '*.a')

    echo "imago_env: flavor '$flavor' ready.  Activate with:" \
         " imago_env $flavor"
}
