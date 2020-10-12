cmake_minimum_required(VERSION 3.14)

# Project setting
project(
  ZERO
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION "ZERO is a multi-purpose game solver")
set(ZERO_AUTHORS "Gabriele Dragotto, Sriram Sankaranarayanan")

# CXX Options
option(CXX "enable C++ compilation" ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# Additional options
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(ZERO_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release ... FORCE)
endif()

# ##############################################################################
# CUSTOM HINTS AND CUSTOM CONFIGURATIONS
# ##############################################################################
set(ortools_dir "/Users/gabrieledragotto/Documents/or-tools")
set(GUROBI_DIR "/home/gurobi/9.0.0/linux64" "/Library/gurobi900/mac64"
               "C:\\dev\\gurobi900\\win64")

set(PATH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/path")

# ##############################################################################
# ##############################################################################

# CMake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# Explicit declaration for custom CMake modules
include(${CMAKE_MODULE_PATH}/FindGUROBI.cmake)
include(${CMAKE_MODULE_PATH}/FindSphinx.cmake)
include(${CMAKE_MODULE_PATH}/FindConan.cmake)
include(${CMAKE_MODULE_PATH}/FindPATH.cmake)
# find_package(ortools CONFIG REQUIRED)
# include(${CMAKE_MODULE_PATH}/FindORTools.cmake)

conan_cmake_run(
  CONANFILE
  conanfile.txt
  BASIC_SETUP
  CMAKE_TARGETS
  NO_OUTPUT_DIRS
  SKIP_STD
  KEEP_RPATHS
  BUILD
  missing
  ENV
  OUTPUT_QUIET)

configure_file(src/support/CMakeConfig.h.in include/support/version.h)

# ##############################################################################
# ZERO - MAIN LIBRARY
# ##############################################################################
add_library(
  ZERO
  src/support/utils.cpp
  src/support/codes.cpp
  src/mathopt/mathopt.cpp
  src/games/games.cpp
  src/mathopt/lcp/lcp.cpp
  src/mathopt/lcp/poly_lcp.cpp
  src/mathopt/lcp/outer_lcp.cpp
  src/mathopt/mp_param/mp_param.cpp
  src/mathopt/mp_param/qp_param.cpp
  src/mathopt/mp_param/ip_param.cpp
  src/solvers/PathSolver.cpp
  src/games/nash.cpp
  src/games/epec.cpp
  src/games/ipg.cpp
  src/games/algorithms/EPEC/epec_polybase.cpp
  src/games/algorithms/EPEC/epec_innerapproximation.cpp
  src/games/algorithms/EPEC/epec_combinatorialpne.cpp
  src/games/algorithms/EPEC/epec_fullenumeration.cpp
  src/games/algorithms/EPEC/epec_outerapproximation.cpp
  src/games/algorithms/IPG/ipg_oracle.cpp
  src/interfaces/epec_models.cpp
  src/interfaces/ipg_models.cpp)
set_target_properties(
  ZERO PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
                  VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(
  ZERO
  PRIVATE blas lapack pthread m ${LUSOL_LIBRARY} ${PATH_LIBRARY}
  PUBLIC CONAN_PKG::boost CONAN_PKG::armadillo CONAN_PKG::openblas
         ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY})
target_include_directories(
  ZERO
  PUBLIC ${GUROBI_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include
         ${PATH_INCLUDE_DIR} ${CONAN_INCLUDE_DIRS}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
# ##############################################################################
# ##############################################################################

# ##############################################################################
# ZERO - DOCS Compiled only for local dev
# ##############################################################################
if(NOT DEFINED ROSSOBIANCO)
  add_subdirectory("docs")
endif()
# ##############################################################################
# ##############################################################################

# ##############################################################################
# ZERO - EPEC EXAMPLES
# ##############################################################################
# Examples add_executable(EPECExample examples/example.cpp)
# target_link_libraries(EPECExample EPECSolve)
# ##############################################################################
# ##############################################################################

# ##############################################################################
# ZERO - EPEC COMMAND LINE
# ##############################################################################
add_executable(EPECShell src/shells/epec.cpp)
target_link_libraries(EPECShell CONAN_PKG::boost ZERO)
set_target_properties(
  EPECShell
  PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)
# ##############################################################################
# ##############################################################################

# ##############################################################################
# ZERO - IPG COMMAND LINE
# ##############################################################################
add_executable(IPGShell src/shells/ipg.cpp)
target_link_libraries(IPGShell CONAN_PKG::boost ZERO)
set_target_properties(
  IPGShell PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
                      VISIBILITY_INLINES_HIDDEN ON)
# ##############################################################################
# ##############################################################################

# ##############################################################################
# ZERO - EPEC TESTS
# ##############################################################################
add_executable(EPECTests test/EPEC/epec.cpp)
target_include_directories(EPECTests PRIVATE test)
target_link_libraries(EPECTests ZERO CONAN_PKG::boost)
set_target_properties(
  EPECTests
  PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
             VISIBILITY_INLINES_HIDDEN ON)
# ##############################################################################
# ##############################################################################

# ##############################################################################
# ZERO - IPG TESTS
# ##############################################################################
add_executable(IPGTests test/IPG/RandomKP.cpp)
target_include_directories(IPGTests PRIVATE test)
target_link_libraries(IPGTests ZERO CONAN_PKG::boost)
set_target_properties(
  IPGTests PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
                      VISIBILITY_INLINES_HIDDEN ON)
# ##############################################################################
# ##############################################################################

message("\n\nZERO Configuration:")
message("\tVersion:\t\t\t${PROJECT_NAME}-${PROJECT_VERSION}")
message("\tCopyright:\t\t\t${ZERO_AUTHORS}\n\n")
message("\tWorking directory:\t${CMAKE_CURRENT_SOURCE_DIR}")
message("\tConan (libs):\t\t${CONAN_LIBS}")
message("\tConan (includes):\t${CONAN_INCLUDE_DIRS}")
message("\tGurobi (lib):\t\t${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY}")
message("\tGurobi (include):\t${GUROBI_INCLUDE_DIRS}")
message("\tPATH (include):\t${PATH_INCLUDE_DIR}")
message("\tPATH (library):\t${PATH_LIBRARY}")
