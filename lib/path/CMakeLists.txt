cmake_minimum_required(VERSION 3.14)

set(ENV{PATH_LICENSE_STRING} "2617827524&Courtesy&&&USR&64785&11_12_2017&1000&PATH&GEN&31_12_2020&0_0_0&5000&0_0")
set(PATH_INCLUDES "${PATH_DIR}/path_include")
set(PATHWRAPPER_INCLUDE_DIR "${PATH_DIR}/include")
# ##############################################################################
# ZERO - PATH WRAPPER
# ##############################################################################
add_library(
        PATHWrapper
        SHARED ${PATH_DIR}/PathLCP.c
)
message("PathLIB is located at ${PATH_LIB}")
set_target_properties(
        PATHWrapper PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(
        PATHWrapper
        PRIVATE m
        PUBLIC ${PATH_LIB})
target_include_directories(
        PATHWrapper
        PUBLIC ${PATH_INCLUDES} ${PATHWRAPPER_INCLUDE_DIR}
        PRIVATE ${PATH_DIR})



