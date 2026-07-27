## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

# CMake file for Imago application: gaussFit

# Document the current location
message("At the gaussFit application.")

add_executable(gaussFit ../kinds.f90 ../readDataSerial.f90 dataModules.f90
   gaussFit.f90)

install(TARGETS gaussFit DESTINATION bin)
