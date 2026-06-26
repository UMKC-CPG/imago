# Top level CMake file for Imago

# Set the minimum required version, the project name, and enable Fortran.
cmake_minimum_required(VERSION 3.1.0)
project(Imago)
enable_language (Fortran)

# Make sure that the default build type is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Establish the installation location on the basis of the $IMAGO_DIR
#   environment variable that is defined in the .imagorc file. Make sure that
#   you source this file (source .imagorc) before you try to build.
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set (CMAKE_INSTALL_PREFIX $ENV{IMAGO_DIR})
   set_property (CACHE CMAKE_INSTALL_PREFIX PROPERTY VALUE $ENV{IMAGO_DIR})
endif (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
message("CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")

# Get the full path of the compiler from the FC environment variable and
#   extract the compiler name from that. Then report the name and the ID (which
#   is found by cmake itself). The idea is that FC may be equal to h5fc but
#   the underlying compiler that h5fc wraps may be gfortran, intel, pgf90,
#   etc. We need all of those pieces of information. (Note that FC should be
#   defined before attempting to build.)
# Another tricky note is that the Fortran compiler that was used to compile
#   the HDF5 library must also be used as the compiler for Imago. Otherwise the
#   module files will not be properly recognized. (At least it is the case that
#   gfortran and intel cannot read eachother's .mod files. I.e. the HDF5.mod
#   file that is created when HDF5 is compiled must be readable by the compiler
#   that you are using to compile Imago.)
get_filename_component (CMAKE_Fortran_COMPILER $ENV{FC} PROGRAM)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("Fortran_COMPILER_NAME= ${Fortran_COMPILER_NAME}")
message("CMAKE_Fortran_COMPLIER= ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_COMPILER_ID= ${CMAKE_Fortran_COMPILER_ID}")

# Define the compiler and linker flags that are needed for the different
#   compilers. So far only gfortran and ifort are considered. You can just
#   change the appropriate line or add a new elseif block as you see fit if
#   you have a particular compiler requirement. I think that this might be
#   easier than attempting to read the user's $FCFLAGS and $LDFLAGS variables
#   because the user may change those to compile other programs and may not
#   retain a record of what they were when he/she compiles Imago. If the user
#   stores it here then he/she will have a record. I guess just be careful to
#   retain that information if you get a new version from github. Or something.
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   set (CMAKE_Fortran_FLAGS_DEBUG "-Og -g -fcheck=all -fimplicit-none -Wall -fbacktrace -DDEBUG")
   #   if (DEFINED ENV{MKLROOT})
   #   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fimplicit-none -Wall -fdefault-integer-8 -m64 -I$ENV{MKLROOT}/include")
   #   set (CMAKE_EXE_LINKER_FLAGS " -L$ENV{MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
   #else()
      set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fimplicit-none -Wall")
      set (CMAKE_EXE_LINKER_FLAGS "-L/usr/lib -I/usr/include -llapack -lblas") # -L,-I for nanoHUB
      #set (CMAKE_EXE_LINKER_FLAGS "-L/usr/lib -I/usr/include") # -L,-I for nanoHUB
      #endif (DEFINED ENV{MKLROOT})
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -heap-arrays")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -no-prec-div -fp-model fast=2 -heap-arrays")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -no-prec-div -fp-model fast=2 -heap-arrays -xHost")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-warn all -pad -ipo -O3 -no-prec-div -fp-model fast=2 -heap-arrays -xHost")
   #   set (CMAKE_Fortran_FLAGS_RELEASE "-warn all -pad -ipo -O3 -no-prec-div -fp-model fast=2 -heap-arrays -xHost -qopt-report=5 -qopt-report-annotate=text -qopt-report-annotate-position=both")
#   set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -heap-arrays")
   set (CMAKE_Fortran_FLAGS_DEBUG "-g -check all -traceback -fp-stack-check -heap-arrays -DDEBUG")
   set (CMAKE_EXE_LINKER_FLAGS "-mkl=sequential")
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
message("CMAKE_Fortran_FLAGS_RELEASE= ${CMAKE_Fortran_FLAGS_RELEASE}")
message("CMAKE_Fortran_FLAGS_DEBUG= ${CMAKE_Fortran_FLAGS_DEBUG}")
message("CMAKE_EXE_LINKER_FLAGS= ${CMAKE_EXE_LINKER_FLAGS}")

# ----------------------------------------------------------------------
# Phase 0a -- opt-in debug instrumentation (default OFF).
#
# These options ADD instrumentation on top of the normal build-type
#   flags set above.  They exist to support the bug-squashing campaign
#   documented in dev/DEBUG.md: they produce instrumented binaries for
#   runtime leak / uninitialized-value / floating-point hunting without
#   disturbing the known-good default build.
#
# Crucially, with every option left at its default (OFF / "none") the
#   compile and link flags are byte-for-byte identical to the historical
#   build.  The instrumentation is only ever appended, never substituted,
#   and the append itself is skipped when nothing is enabled.  Turn the
#   options on only for a throwaway / instrumented build tree, e.g.:
#     cmake -DCMAKE_BUILD_TYPE=Debug -DIMAGO_SANITIZE=address ../..
#
# Each option maps to one debugging capability:
#   IMAGO_CHECKS     runtime checks: array bounds, etc.
#   IMAGO_FPE_TRAP   trap floating-point exceptions (invalid, divide-by-
#                      zero, overflow) so a bad value halts at its origin
#   IMAGO_INIT_SNAN  initialize reals to a signaling NaN so that any use
#                      of an uninitialized value traps instead of silently
#                      producing a wrong answer
#   IMAGO_WARN_EXTRA the extra, audit-grade compiler warning set
#   IMAGO_SANITIZE   select a sanitizer: none, address, undefined, leak
# ----------------------------------------------------------------------
option (IMAGO_CHECKS     "Runtime checks (array bounds, etc.)"     OFF)
option (IMAGO_FPE_TRAP   "Trap FP exceptions (inv/zero/overflow)"  OFF)
option (IMAGO_INIT_SNAN  "Init reals to signaling NaN"             OFF)
option (IMAGO_WARN_EXTRA "Extra audit-grade warning set"           OFF)
set (IMAGO_SANITIZE "none" CACHE STRING
   "Sanitizer: none, address, undefined, or leak")

# Accumulate the requested instrumentation for the detected compiler
#   into a single string, then append it once so that it composes with
#   whatever build type (Release / Debug) was chosen.  Every fragment
#   carries its own leading space, and the final append is guarded so
#   that an empty set leaves CMAKE_Fortran_FLAGS exactly as it was.
set (IMAGO_INSTRUMENT_FLAGS "")
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   if (IMAGO_CHECKS)
      set (IMAGO_INSTRUMENT_FLAGS "${IMAGO_INSTRUMENT_FLAGS} -fcheck=all")
   endif ()
   if (IMAGO_FPE_TRAP)
      set (IMAGO_INSTRUMENT_FLAGS
         "${IMAGO_INSTRUMENT_FLAGS} -ffpe-trap=invalid,zero,overflow")
      set (IMAGO_INSTRUMENT_FLAGS "${IMAGO_INSTRUMENT_FLAGS} -fbacktrace")
   endif ()
   if (IMAGO_INIT_SNAN)
      set (IMAGO_INSTRUMENT_FLAGS
         "${IMAGO_INSTRUMENT_FLAGS} -finit-real=snan -finit-integer=-99999")
   endif ()
   if (IMAGO_WARN_EXTRA)
      set (IMAGO_INSTRUMENT_FLAGS
         "${IMAGO_INSTRUMENT_FLAGS} -Wextra -Wuninitialized")
      set (IMAGO_INSTRUMENT_FLAGS
         "${IMAGO_INSTRUMENT_FLAGS} -Wrealloc-lhs-all -Warray-temporaries")
      set (IMAGO_INSTRUMENT_FLAGS
         "${IMAGO_INSTRUMENT_FLAGS} -Wimplicit-interface")
   endif ()
   if (NOT IMAGO_SANITIZE STREQUAL "none")
      set (IMAGO_INSTRUMENT_FLAGS
         "${IMAGO_INSTRUMENT_FLAGS} -fsanitize=${IMAGO_SANITIZE}")
      set (CMAKE_EXE_LINKER_FLAGS
         "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=${IMAGO_SANITIZE}")
   endif ()
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
# Intel-compiler equivalents are added in Phase 0c.

if (NOT IMAGO_INSTRUMENT_FLAGS STREQUAL "")
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}${IMAGO_INSTRUMENT_FLAGS}")
endif ()
message("IMAGO_INSTRUMENT_FLAGS=${IMAGO_INSTRUMENT_FLAGS}")

# Proceed to the source directory and its subdirectories.
add_subdirectory(src)
