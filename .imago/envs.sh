#!/usr/bin/env bash
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

# envs.sh -- shell helpers for the Imago toolchain.  It defines two
#   functions: imago_env, which switches the active toolchain between
#   build "flavors", and imago_scripts, which refreshes the production
#   scripts in bin/ without recompiling the Fortran engine.
#
# This file is meant to be *sourced* (not executed), normally from your
# shell startup right after sourcing the Imago rc file, e.g.:
#
#     source "$IMAGO_DIR/.imago/imagorc"
#     source "$IMAGO_DIR/.imago/envs.sh"
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
#     imago_env [flavor]      activate a flavor; no arg (or "prod" /
#                             "production") restores the production build
#     imago_env --list        list the flavors that are built
#     imago_env --add NAME    build preset gfortran-NAME and assemble
#                             the flavor bin, so `imago_env NAME` works
#     imago_env --help        show this usage
#
# Examples:
#     imago_env mpi           # the parallel build (MPI + ELPA)
#     imago_env serial        # the serial fallback build
#     imago_env asan          # run the engine under AddressSanitizer
#     imago_env audit         # run with checks + FP traps + warnings
#     imago_env               # back to the installed production build


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
                "${BASH_SOURCE[0]:-$IMAGO_DIR/.imago/envs.sh}" \
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
        --add|-a)
            _imago_env_add "$2"
            return $?
            ;;
    esac

    # Resolve the requested flavor to a bin directory.  With no
    #   argument, or "prod"/"production", the target is the installed
    #   production build -- whatever it happens to be (serial today,
    #   parallel once the default flips; ARCHITECTURE 4.2).  Every other
    #   name, including "serial" and "mpi", is a flavor under
    #   envs/<name>/bin.  The switcher deliberately makes no assumption
    #   about which build production is, so a build name is never
    #   conflated with "the installed default".
    local flavor="${1:-prod}"
    local target_bin
    case "$flavor" in
        prod|production)
            flavor="prod"
            target_bin="$IMAGO_DIR/bin"
            ;;
        *)
            target_bin="$IMAGO_DIR/envs/$flavor/bin"
            ;;
    esac

    if [ ! -d "$target_bin" ]; then
        echo "imago_env: flavor '$flavor' is not built ($target_bin)." >&2
        echo "           add it with:  imago_env --add $flavor" >&2
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


# Add flavor $1: build it from the matching CMake preset (gfortran-$1)
#   and (re)assemble its bin as symlinks to the production toolchain
#   with the freshly built engine executables overlaid.  A full preset
#   build is run so every engine executable (not just imago/imagoG) is
#   the instrumented one.
_imago_env_add() {
    local flavor="$1"
    if [ -z "$flavor" ]; then
        echo "imago_env --add: a flavor name is required" \
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


# imago_scripts -- refresh the production toolchain (the Python/bash
#   scripts, the kaleidoscope package, and the .imago/ resource-control
#   files) into $IMAGO_DIR/bin and $IMAGO_DIR/.imago WITHOUT recompiling
#   the Fortran engine.
#
# The engine is compiled per *flavor* (see imago_env); rebuilding it has
#   nothing to do with updating a script, so this path never invokes the
#   compiler.  It drives CMake's install step for only the script-related
#   install components -- "scripts", "kaleidoscope" and "rc", tagged in
#   src/scripts/CMakeLists.txt -- and "cmake --install" only runs install
#   rules, it never builds a target.  Routing through CMake keeps a single
#   source of truth for "what is a script" (the SCRIPTS list in that
#   CMakeLists), so there is no parallel copy list here to drift out of
#   step with it.
#
# Options:
#     --no-rc            leave the .imago/ resource-control files alone
#     --no-kaleidoscope  skip the kaleidoscope flight-runner package
#
# The first call configures the production tree build/release if it does
#   not exist yet.  That is a *configure*, not a build: CMake still probes
#   the Fortran compiler when it configures the project, so the compiler
#   wrapper must be on PATH for that one step (activate the cpg
#   environment).  Every later call skips straight to the install and is
#   instant.
imago_scripts() {
    if [ -z "$IMAGO_DIR" ]; then
        echo "imago_scripts: IMAGO_DIR is not set; source the Imago" \
             "rc file first." >&2
        return 1
    fi

    # Parse the two opt-out switches.  Both default to "on" so that a
    #   bare "imago_scripts" refreshes the whole toolchain.
    local install_rc=1 install_kaleidoscope=1 arg
    for arg in "$@"; do
        case "$arg" in
            --no-rc)           install_rc=0 ;;
            --no-kaleidoscope) install_kaleidoscope=0 ;;
            *)
                echo "imago_scripts: unknown option '$arg'" >&2
                echo "  usage: imago_scripts [--no-rc]" \
                     "[--no-kaleidoscope]" >&2
                return 1
                ;;
        esac
    done

    # The production install lives in the hand-driven build/release tree,
    #   whose install prefix is $IMAGO_DIR (BUILD.md).  Configure it once
    #   if it is missing; configure builds nothing, so the engine is never
    #   compiled on this path.
    local tree="$IMAGO_DIR/build/release"
    if [ ! -f "$tree/CMakeCache.txt" ]; then
        echo "imago_scripts: configuring $tree (one-time) ..."
        mkdir -p "$tree" || return 1
        if ! ( cd "$tree" && cmake ../.. ); then
            echo "imago_scripts: configure failed -- is the compiler" \
                 "wrapper (h5fc) on PATH?  Activate cpg and retry." >&2
            return 1
        fi
    fi

    # Install only the requested components.  "--prefix $IMAGO_DIR" makes
    #   the destination explicit and independent of how the tree happened
    #   to be configured, so scripts always land in the production bin.
    echo "imago_scripts: installing scripts into $IMAGO_DIR/bin ..."
    cmake --install "$tree" --prefix "$IMAGO_DIR" \
          --component scripts || return 1

    if [ "$install_kaleidoscope" -eq 1 ]; then
        echo "imago_scripts: installing the kaleidoscope package ..."
        cmake --install "$tree" --prefix "$IMAGO_DIR" \
              --component kaleidoscope || return 1
    fi

    if [ "$install_rc" -eq 1 ]; then
        echo "imago_scripts: installing rc files into $IMAGO_DIR" \
             "/.imago ..."
        cmake --install "$tree" --prefix "$IMAGO_DIR" \
              --component rc || return 1
    fi

    echo "imago_scripts: done."
}
