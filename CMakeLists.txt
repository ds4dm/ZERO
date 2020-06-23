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
set(GUROBI_PATH_LINUX "/home/gurobi/9.0.0/linux64")
set(GUROBI_PATH_MAC "/Library/gurobi900/mac64")
set(GUROBI_PATH_WINDOWS "C:\\dev\\gurobi900\\win64")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


############################################
################ DO NOT MODIFY #############
############################################
cmake_minimum_required(VERSION 3.10)
project(
        EPECSolve_OG
        VERSION 3.0.0
        LANGUAGES CXX
        DESCRIPTION "Equilibrium Problems with Equilibrium Constraints solver"
)
add_compile_options(-std=c++11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_SHARED_LINKER_FLAGS "-lgurobi_c++ -lgurobi90 -lm -O2 -lopenblas -llapack -lpthread")

#From https://github.com/LiangliangNan/PolyFit/blob/master/cmake/FindGUROBI.cmake
if (NOT GUROBI_FOUND)
    # Hardcoded search paths
    set(SEARCH_PATHS_FOR_HEADERS
        "$ENV{GUROBI_HOME}/include"
        "${GUROBI_PATH_MAC}/include"
        "${GUROBI_PATH_LINUX}/include"
        "${GUROBI_PATH_WINDOWS}/include"
        )

    set(SEARCH_PATHS_FOR_LIBRARIES
        "$ENV{GUROBI_HOME}/lib"
        "/opt/gurobi900/linux64/lib"
        "${GUROBI_PATH_MAC}/lib"
        "${GUROBI_PATH_LINUX}/lib"
        "${GUROBI_PATH_WINDOWS}/lib"
        )

    find_path(GUROBI_INCLUDE_DIR gurobi_c++.h
              PATHS ${SEARCH_PATHS_FOR_HEADERS}
              )


    find_library(GUROBI_C_LIBRARY
                 NAMES gurobi90
                 PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
                 )

    find_library(GUROBI_CXX_LIBRARY_DEBUG
                 NAMES gurobi_c++ gurobi_c++mdd2017
                 PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
                 )

    find_library(GUROBI_CXX_LIBRARY_RELEASE
                 NAMES gurobi_c++ gurobi_c++md2017
                 PATHS ${SEARCH_PATHS_FOR_LIBRARIES}
                 )

    # setup header file directories
    set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})

    # setup libraries files
    set(GUROBI_LIBRARIES ${GUROBI_C_LIBRARY}
        debug ${GUROBI_CXX_LIBRARY_DEBUG}
        optimized ${GUROBI_CXX_LIBRARY_RELEASE})

endif ()
include(FindPackageHandleStandardArgs)
set(Boost_USE_STATIC_LIBS ON)
if (ROSSOBIANCO)
    find_package(Boost 1.70.0 REQUIRED COMPONENTS program_options unit_test_framework log log_setup)
else ()
    find_package(Boost 1.72.0 REQUIRED COMPONENTS program_options unit_test_framework log log_setup)
endif ()
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_INCLUDE_DIRS)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARIES)
mark_as_advanced(GUROBI_LIBRARIES GUROBI_INCLUDE_DIRS GUROBI_INCLUDE_DIR)
include_directories(${GUROBI_INCLUDE_DIR})
include_directories(${BOOST_ROOT}/include)
include_directories(${ARMA_PATH}/include)

include_directories("src" "include")
configure_file(src/CMakeConfig.h.in ${PROJECT_SOURCE_DIR}/include/version.h)

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
                      #CXX_VISIBILITY_PRESET hidden
                      )
target_link_libraries(EPECSolve
                      Boost::program_options
                      Boost::log Boost::log_setup
                      ${GUROBI_LIBRARIES}
                      ${ARMA_LIB})
add_executable(EPEC
               src/CommandLineInterface.cpp)
#set_target_properties(EPEC PROPERTIES CXX_VISIBILITY_PRESET hidden)
add_executable(EPECExample
               examples/example.cpp)
add_executable(EPECtest
               test/EPEC.cpp)
target_link_libraries(EPECExample EPECSolve)
target_link_libraries(EPEC EPECSolve)
target_link_libraries(EPECtest EPECSolve Boost::unit_test_framework)
message("\n\nEPEC Cmake configuration terminated:")
message("\tWorking directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("\tArmadillo (lib) : ${ARMA_PATH}/include")
message("\tArmadillo (include) : ${ARMA_LIB}")
message("\tGurobi (lib) : ${GUROBI_LIBRARIES}")
message("\tGurobi (include) : ${GUROBI_INCLUDE_DIRS}")
message("\tBoost (include) : ${BOOST_ROOT}/include ")
############################################
