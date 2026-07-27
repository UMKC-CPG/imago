## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

# CMake file for Imago application: makeKPoints

# Document the current location
message("At the makeKPoints application.")

add_executable(makekpoints ../kinds.f90 ../constants.f90 makekpoints.F90)

install(TARGETS makekpoints DESTINATION bin)
