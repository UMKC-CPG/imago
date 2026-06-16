# Source level CMake file for Imago

# Document that we are at the source level.
message("At the source level.")

# List all subdirectories to be made.  The test suite (src/tests) is
# deliberately NOT a subdirectory here: tests are not a deployable
# artifact, so they are not installed.  Run them in place from
# src/tests instead (e.g. `python3 -m pytest` from that directory).
add_subdirectory(applySpaceGroup)
add_subdirectory(atomSCF)
add_subdirectory(contract)
add_subdirectory(gaussFit)
add_subdirectory(makeKPoints)
add_subdirectory(scripts)
add_subdirectory(data)
add_subdirectory(imago)
