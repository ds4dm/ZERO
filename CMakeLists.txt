# ##############################################################################
# ############### CONFIGURATION #############
# Here you can edit with your paths #####
# ##############################################################################

if(ROSSOBIANCO)
  set(BOOST_ROOT "/home/x86_64-unknown-linux_ol7-gnu/boost-1.70.0")
  set(ARMA_PATH "/home/draggabr/.local/armadillo")
  set(ARMA_LIB "${ARMA_PATH}/lib64/libarmadillo.so")
else()
  set(BOOST_ROOT "/usr/local/boost-1.70.0/")
  set(ARMA_PATH "/usr/local")
  set(ARMA_LIB "/usr/local/lib/libarmadillo.dylib")
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ##############################################################################
# ############### DO NOT MODIFY #############
# ##############################################################################
cmake_minimum_required(VERSION 3.14)
option(CXX "enable C++ compilation" ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(EXECUTABLE_OUTPUT_PATH "bin")
project(
  ZERO
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION "ZERO is a multi-purpose game solver")

include_directories(${ARMA_PATH}/include)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(GUROBI_DIR "/home/gurobi/9.0.0/linux64" "/Library/gurobi900/mac64"
               "C:\\dev\\gurobi900\\win64")
find_package(GUROBI REQUIRED)

set(Boost_USE_MULTITHREADED TRUE)
find_package(Boost 1.70.0 REQUIRED COMPONENTS program_options
                                              unit_test_framework log log_setup)

configure_file(src/support/CMakeConfig.h.in include/support/version.h)

# ZERO - Main library
add_library(
  ZERO
  src/support/utils.cpp
  src/support/codes.cpp
  src/games/games.cpp
  src/lcp/lcp.cpp
  src/lcp/poly_lcp.cpp
  src/lcp/outer_lcp.cpp
  src/games/qpmp.cpp
  src/games/nash.cpp
  src/games/epec.cpp
  # src/games/ipg.cpp
  src/interfaces/epec_models.cpp
  src/algorithms/EPEC/epec_polybase.cpp
  src/algorithms/EPEC/epec_innerapproximation.cpp
  src/algorithms/EPEC/epec_combinatorialpne.cpp
  src/algorithms/EPEC/epec_fullenumeration.cpp
  src/algorithms/EPEC/epec_outerapproximation.cpp
  # src/algorithms/IPG/ipg_oracle.cpp
)
set_target_properties(
  ZERO PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
                  VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(
  ZERO
  PRIVATE blas lapack pthread m
  PUBLIC Boost::log Boost::log_setup ${ARMA_LIB} ${GUROBI_LIBRARY}
         ${GUROBI_CXX_LIBRARY})
target_include_directories(
  ZERO
  PUBLIC ${GUROBI_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include
         ${CMAKE_CURRENT_SOURCE_DIR}/lib ${Boost_INCLUDE_DIR}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Examples add_executable(EPECExample examples/example.cpp)
# target_link_libraries(EPECExample EPECSolve)

# EPEC Command line interface
add_executable(EPECShell src/shells/epec.cpp)
target_link_libraries(EPECShell Boost::program_options ZERO)
set_target_properties(
  EPECShell
  PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)

# EPEC Tests
add_executable(EPECTests test/EPEC/EPEC_test.cpp)
target_include_directories(EPECTests PRIVATE test)
target_link_libraries(EPECTests ZERO Boost::unit_test_framework)

message("\n\nZERO cMake configuration terminated:")
message("\tWorking directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("\tArmadillo (lib) : ${ARMA_PATH}/include")
message("\tArmadillo (include) : ${ARMA_LIB}")
message("\tGurobi (lib) : ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY}")
message("\tGurobi (include) : ${GUROBI_INCLUDE_DIRS}")
message("\tBoost (include) : ${Boost_INCLUDE_DIR} ")
message("\tBoost (lib) : ${Boost_LIBRARIES}")
# ##############################################################################
