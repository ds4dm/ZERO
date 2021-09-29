cmake_minimum_required(VERSION 3.16)

#######################################
# PROJECT SETTING AND VERSIONING
#######################################
include(cmake/ZERO.cmake)
project(
        ZERO
        VERSION ${ZERO_VERSION_MAJOR}.${ZERO_VERSION_MINOR}.${ZERO_VERSION_PATCH}
        LANGUAGES CXX
        DESCRIPTION ${ZERO_DESC})

message("Configuring ZEROVersion")
configure_file(src/support/CMakeConfig.h.in include/support/version.h)

#######################################
# PROJECT OPTIONS
#######################################
option(ZERO_INSTALL_DEPS "Install ZERO dependencies through conda" OFF)


#######################################
# NON-STANDARD FLAGS
#######################################
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release ... FORCE)
endif ()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -frtti -fvisibility-inlines-hidden -fvisibility=hidden")
#DEBUG AND PROFILING ON APPLE SILICON
set(CMAKE_OSX_ARCHITECTURES "x86_64")
if (CMAKE_BUILD_TYPE EQUAL "DEBUG" AND NOT ROSSOBIANCO)
    message("Using MacOs debug setting")
    set(CMAKE_CXX_COMPILER "clang++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNEDEBUG -fno-inline-functions-called-once -fno-inline-functions -fno-optimize-sibling-calls -stdlib=libc++ -lc++abi")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -gdwarf-3")
endif ()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(ZERO_DEP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

#######################################
# CUSTOM HINTS FOR DEPENDENCIES
#######################################
SET(ZERO_GUROBIHINT "" CACHE STRING "Hint for the Gurobi Path")
set(GUROBI_DIR "/home/gurobi/9.1.1/linux64/" "/Library/gurobi910/mac64" "C:\\dev\\gurobi900\\win64" ${ZERO_GUROBIHINT})
set(PATH_DIR "${ZERO_DEP_DIR}/path")
set(CONDA_ENV "${ZERO_DEP_DIR}/conda")
set(SCIP_ROOT_DIR "${CONDA_ENV}" )
#set(COIN_HINT "/Users/gabrieledragotto/dist" "/local_workspace/draggabr/coin")
set(COIN_HINT "${CONDA_ENV}")

# Get Conda
find_package(Conda REQUIRED)
if (NOT ${FOUND_CONDA})
    message("Warning. No Conda found.")
endif ()
#Install deps through Conda?
if (ZERO_INSTALL_DEPS)
    execute_process(COMMAND conda create --yes --prefix ${CONDA_ENV})
    MESSAGE("COMMAND conda install -c conda-forge --prefix ${CONDA_ENV} armadillo coin-or-cgl boost rapidjson coin-or-osi -y")
    execute_process(COMMAND conda install -c conda-forge --prefix ${CONDA_ENV} armadillo coin-or-cgl boost rapidjson coin-or-osi -y)
endif ()
set(CONDA_INCLUDE_DIR "${CONDA_ENV}/include")

# Boost
set(Boost_NO_SYSTEM_PATHS ON)
set(BOOST_ROOT "${CONDA_ENV}")
set(Boost_NO_CMAKE ON)
find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework HINTS ${CONDA_ENV} PATHS ${CONDA_ENV} NO_DEFAULT_PATH)
# Gurobi
find_package(GUROBI REQUIRED)
# Path
find_package(PATH REQUIRED)
# Armadillo
find_package(Armadillo REQUIRED HINTS ${CONDA_ENV} PATHS ${CONDA_ENV} NO_DEFAULT_PATH)
# CoinUtils and CGL
find_package(Coin REQUIRED)
# Scip
#find_package(SCIP REQUIRED)
# CustomCoin
include(cmake/InstallCustomCoin.cmake)

# Update flags on Conda, loguru, and Armadillo
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{DEBUG_CXXFLAGS} $ENV{DEBUG_CPPFLAGS} -DARMA_NO_DEBUG")
add_definitions(-DLOGURU_WITH_STREAMS=1)

#######################################
# ZERO - MAIN LIBRARY
#######################################

add_library(
        ZERO
        include/zero.h
        lib/loguru/loguru.cpp
        src/support/utils.cpp
        src/support/codes.cpp
        src/mathopt/mathopt.cpp
        src/mathopt/lcp/lcp.cpp
        src/mathopt/lcp/poly_lcp.cpp
        src/mathopt/mp_param/mp_param.cpp
        src/mathopt/mp_param/qp_param.cpp
        src/mathopt/mp_param/ip_param.cpp
        src/solvers/PathSolver.cpp
        src/games/nash.cpp
        src/games/epec.cpp
        src/games/ipg.cpp
        src/games/algorithms/EPEC/epec_polybase.cpp
        src/games/algorithms/EPEC/epec_innerapp.cpp
        src/games/algorithms/EPEC/epec_combPNE.cpp
        src/games/algorithms/EPEC/epec_fullenum.cpp
        src/games/algorithms/EPEC/epec_cutandplay.cpp
        src/games/algorithms/IPG/ipg_cutandplay.cpp
        src/games/algorithms/IPG/ipg_zeroregrets.cpp
        src/interfaces/epec_models.cpp
        src/interfaces/ipg_models.cpp)
set_target_properties(
        ZERO PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(
        ZERO
        PRIVATE blas lapack pthread m ${LUSOL_LIBRARY} ${PATH_LIBRARY} ${CMAKE_DL_LIBS}
        PUBLIC ${ARMADILLO_LIBRARIES} ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY}
        ${COIN_UTILS} ${COIN_CGL}  ${COIN_OSI} libCustomCoin)
target_include_directories(
        ZERO
        PUBLIC ${GUROBI_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${PATH_INCLUDE_DIR} ${CONDA_INCLUDE_DIR} ${COIN_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/loguru/
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)


#######################################
# ZERO - EPEC COMMAND LINE
#######################################
add_executable(EPECShell src/shells/epec.cpp)
target_link_libraries(EPECShell Boost::program_options Boost::unit_test_framework ${CMAKE_DL_LIBS} ZERO)
set_target_properties(
        EPECShell
        PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)

#######################################
# ZERO - IPG COMMAND LINE
#######################################
add_executable(IPGShell src/shells/ipg.cpp)
target_link_libraries(IPGShell Boost::program_options ${CMAKE_DL_LIBS} ZERO)
set_target_properties(
        IPGShell PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)

#######################################
# ZERO - TESTS
#######################################
add_executable(ZEROTests
               test/tests.cpp
               test/core_tests.cpp
               test/epec_test.cpp
               test/ipg_tests.cpp)
target_compile_definitions(ZEROTests PRIVATE BOOST_ALL_DYN_LINK)
target_include_directories(ZEROTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ZEROTests ZERO ${CMAKE_DL_LIBS} Boost::unit_test_framework)
set_target_properties(
        ZEROTests
        PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)

#######################################
# ZERO - Compiling an example
#######################################
add_executable(IPGExample examples/IPG.cpp)
target_include_directories(IPGExample PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(IPGExample ZERO)
set_target_properties(
        IPGExample PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)


#######################################
# OUTPUT
#######################################
message("\n\nZERO Configuration:")
message("\tVersion:\t\t\t${PROJECT_NAME}-${PROJECT_VERSION}")
message("\tCopyright:\t\t\t${ZERO_AUTHORS}\n\n")
message("\tWorking directory:\t${CMAKE_CURRENT_SOURCE_DIR}")
message("\tConda (includes):\t${CONDA_INCLUDE_DIR}")
message("\tBoost (include):\t${Boost_INCLUDE_DIRS}")
message("\tBoost (lib):\t\t${Boost_LIBRARIES}")
message("\tGurobi (include):\t${GUROBI_INCLUDE_DIRS}")
message("\tGurobi (lib):\t\t${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY}")
message("\tArmadillo (include):\t${ARMADILLO_INCLUDE_DIRS}")
message("\tArmadillo (lib):\t\t${ARMADILLO_LIBRARIES}")
message("\tCoin (lib):\t\t${COIN_UTILS} ${COIN_CGL} ${COIN_OSI}")
message("\tCoin (include):\t\t${COIN_INCLUDE_DIR}")
message("\tPATH (include):\t${PATH_INCLUDE_DIR}")
message("\tPATH (library):\t${PATH_LIBRARY}")
#message("\tSCIP (include):\t${SCIP_INCLUDE_DIRS}")
#message("\tSCIP (lib):\t\t${SCIP_LIBRARIES}")