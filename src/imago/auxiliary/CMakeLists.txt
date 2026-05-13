# CMake file for Imago application: auxiliary tools.

# Document the current location
message("At the auxiliary imago application level.")

add_executable(imagoKKc ../../kinds.f90 ../../constants.f90
   imagoKKc.f90)
add_executable(imagoRhoV ../../kinds.f90 ../../constants.f90
   imagoRhoV.f90)
add_executable(gaussIntgTest ../../kinds.f90 ../../constants.f90
   gaussIntgTest.f90)

install(TARGETS imagoKKc imagoRhoV gaussIntgTest DESTINATION bin)
