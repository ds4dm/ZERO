cmake_minimum_required(VERSION 3.14)

set(ENV{PATH_LICENSE_STRING}
    "2617827524&Courtesy&&&USR&64785&11_12_2017&1000&PATH&GEN&31_12_2020&0_0_0&5000&0_0"
)
# ##############################################################################
# ZERO - PATH WRAPPER
# ##############################################################################

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

add_library(PATHWrapper STATIC ${PATH_DIR}/src/PathLCP.cpp)
message("PathLIB is located at ${PATH_LIBRARY}")

target_include_directories(
  PATHWrapper
  PUBLIC ${PATH_INCLUDE_DIR}
  PRIVATE ${PATH_DIR}/src)

target_link_libraries(PATHWrapper PUBLIC ${LUSOL_LIBRARY} ${PATH_LIBRARY} m)

set_target_properties(
  PATHWrapper
  PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)
