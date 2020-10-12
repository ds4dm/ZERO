message("Loading PATH")
set(ENV{PATH_LICENSE_STRING}
    "2617827524&Courtesy&&&USR&64785&11_12_2017&1000&PATH&GEN&31_12_2020&0_0_0&5000&0_0"
)
find_library(
  PATH_LIBRARY
  NAMES path path47 path50
  HINTS ${PATH_DIR}/lib
  PATH_SUFFIXES lib)

find_library(
  LUSOL_LIBRARY
  NAMES lusol
  HINTS ${PATH_DIR}/lib
  PATH_SUFFIXES lib)
message("PathLIB is located at ${PATH_LIBRARY}")
set(PATH_INCLUDE_DIR "${PATH_DIR}/include")

message("LUSOL is located at ${LUSOL_LIBRARY}")
