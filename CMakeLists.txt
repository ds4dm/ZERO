############################################
################ CONFIGURATION #############
#### Here you can edit with your paths #####
############################################

if (ROSSOBIANCO)
    set(BOOST_ROOT "/home/x86_64-unknown-linux_ol7-gnu/boost-1.70.0")
    set(ARMA_PATH "/home/draggabr/.local/armadillo")
    set(ARMA_LIB "${ARMA_PATH}/lib64/libarmadillo.so")
ELSE ()
    set(BOOST_ROOT "/usr/local/boost-1.70.0/")
    set(ARMA_PATH "/usr/local")
    set(ARMA_LIB "/usr/local/lib/libarmadillo.dylib")
    set(Boost_USE_MULTITHREADED TRUE)
ENDIF ()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

############################################
################ DO NOT MODIFY #############
############################################
cmake_minimum_required(VERSION 3.8)
option(CXX "enable C++ compilation" ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
project(
        ZERO
        VERSION 1.0.0
        LANGUAGES CXX
        DESCRIPTION "ZERO is a multi-purpose game solver"
)

include_directories(${ARMA_PATH}/include)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(GUROBI_DIR "/home/gurobi/9.0.0/linux64" "/Library/gurobi900/mac64" "C:\\dev\\gurobi900\\win64")

find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.70.0 REQUIRED COMPONENTS program_options unit_test_framework log log_setup)
include_directories(${BOOST_ROOT}/include)

include_directories("src" "include" "lib")
configure_file(src/CMakeConfig.h.in ${PROJECT_SOURCE_DIR}/include/version.h)


#EPECSOLVE - Main library
add_library(ZERO
            src/games/games.cpp
            src/games/epec.cpp
            src/games/qpmp.cpp
            src/games/ipg.cpp
            src/games/nash.cpp
            src/lcp/LCP.cpp
            src/lcp/PolyLCP.cpp
            src/Utils.cpp
            src/interfaces/EPECModels.cpp
            src/algorithms/EPEC/PolyBase.cpp
            src/algorithms/EPEC/InnerApproximation.cpp
            src/algorithms/EPEC/CombinatorialPNE.cpp
            src/algorithms/EPEC/FullEnumeration.cpp
            src/algorithms/EPEC/OuterApproximation.cpp
            src/lcp/OuterLCP.cpp)
target_include_directories(ZERO PRIVATE src)
set_target_properties(ZERO PROPERTIES
                      POSITION_INDEPENDENT_CODE ON
                      )
target_link_libraries(ZERO
                      Boost::log Boost::log_setup
                      -O2 blas lapack pthread m
                      ${ARMA_LIB}
                      ${GUROBI_LIBRARY}
                      optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})


#Examples
#add_executable(EPECExample
#               examples/example.cpp)
#target_link_libraries(EPECExample EPECSolve)

#EPEC Command line interface
add_executable(EPECShell
               src/shells/EPEC.cpp)
target_link_libraries(EPECShell
                      Boost::program_options
                      Boost::log Boost::log_setup
                      ZERO)

#EPEC Tests
add_executable(EPECTests
               test/EPEC/EPEC.cpp)
target_link_libraries(EPECTests ZERO Boost::unit_test_framework)


message("\n\nZERO cMake configuration terminated:")
message("\tWorking directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("\tArmadillo (lib) : ${ARMA_PATH}/include")
message("\tArmadillo (include) : ${ARMA_LIB}")
message("\tGurobi (lib) : ${GUROBI_LIBRARY}")
message("\tGurobi (include) : ${GUROBI_INCLUDE_DIRS}")
message("\tBoost (include) : ${BOOST_ROOT}/include ")
############################################
