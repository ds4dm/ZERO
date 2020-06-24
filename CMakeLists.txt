############################################
################ CONFIGURATION #############
#### Here you can edit with your paths #####
############################################

if (ROSSOBIANCO)
    set(BOOST_ROOT "/home/x86_64-unknown-linux_ol7-gnu/boost-1.70.0")
    set(ARMA_PATH "/home/draggabr/.local/armadillo")
    set(ARMA_LIB "${ARMA_PATH}/lib64/libarmadillo.so")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gdwarf-3")
ELSE ()
    set(BOOST_ROOT "/usr/local/boost-1.70.0/")
    set(ARMA_PATH "/usr/local")
    set(ARMA_LIB "/usr/local/lib/libarmadillo.dylib")
    set(Boost_USE_MULTITHREADED TRUE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++11-narrowing -Wdefaulted-function-deleted")
ENDIF ()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

############################################
################ DO NOT MODIFY #############
############################################
cmake_minimum_required(VERSION 3.8)
option(CXX "enable C++ compilation" ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frtti -fvisibility-inlines-hidden -fvisibility=hidden")
#set(CMAKE_SHARED_LINKER_FLAGS "lm -O2 -lopenblas -llapack -lpthread")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
project(
        EPECSolve_OG
        VERSION 3.0.0
        LANGUAGES CXX
        DESCRIPTION "Equilibrium Problems with Equilibrium Constraints solver"
)
add_compile_options(-std=c++11)

set(GUROBI_DIR "/home/gurobi/9.0.0/linux64" "/Library/gurobi900/mac64" "C:\\dev\\gurobi900\\win64")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(Boost_USE_STATIC_LIBS ON)
if (ROSSOBIANCO)
    find_package(Boost 1.70.0 REQUIRED COMPONENTS program_options unit_test_framework log log_setup)
else ()
    find_package(Boost 1.72.0 REQUIRED COMPONENTS program_options unit_test_framework log log_setup)
endif ()

find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})
include_directories(${BOOST_ROOT}/include)
include_directories(${ARMA_PATH}/include)

include_directories("src" "include" "lib")
configure_file(src/CMakeConfig.h.in ${PROJECT_SOURCE_DIR}/include/version.h)


#EPECSOLVE - Main library
add_library(EPECSolve
            src/games/games.cpp
            src/games/epec.cpp
            src/games/qpmp.cpp
            src/games/ipg.cpp
            src/games/nash.cpp
            src/lcp/LCP.cpp
            src/lcp/PolyLCP.cpp
            src/Utils.cpp
            src/Models.cpp
            src/algorithms/PolyBase.cpp
            src/algorithms/InnerApproximation.cpp
            src/algorithms/CombinatorialPNE.cpp
            src/algorithms/FullEnumeration.cpp
            src/algorithms/OuterApproximation.cpp
            src/lcp/OuterLCP.cpp)
target_include_directories(EPECSolve PRIVATE src)
set_target_properties(EPECSolve PROPERTIES
                      POSITION_INDEPENDENT_CODE ON
                      )
target_link_libraries(EPECSolve
                      Boost::program_options
                      Boost::log Boost::log_setup
                      ${GUROBI_LIBRARY}
                      ${ARMA_LIB}
                      optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})


#Examples
add_executable(EPECExample
               examples/example.cpp)
target_link_libraries(EPECExample EPECSolve)

#Command line interface
add_executable(EPEC
               src/CommandLineInterface.cpp)
target_link_libraries(EPEC EPECSolve)

#Tests
add_executable(EPECtest
               test/EPEC.cpp)
target_link_libraries(EPECtest EPECSolve Boost::unit_test_framework)


message("\n\nEPEC Cmake configuration terminated:")
message("\tWorking directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("\tArmadillo (lib) : ${ARMA_PATH}/include")
message("\tArmadillo (include) : ${ARMA_LIB}")
message("\tGurobi (lib) : ${GUROBI_LIBRARY}")
message("\tGurobi (include) : ${GUROBI_INCLUDE_DIRS}")
message("\tBoost (include) : ${BOOST_ROOT}/include ")
############################################
