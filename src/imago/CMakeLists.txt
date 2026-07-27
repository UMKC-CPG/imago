## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

# imago level CMake file for the Imago electronic-structure engine.

# Document that we are at the imago level.
message("At the imago level.")

# List all subdirectories to be made.
add_subdirectory(real)
add_subdirectory(complex)
add_subdirectory(auxiliary)
