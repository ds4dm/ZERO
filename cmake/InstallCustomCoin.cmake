message("Installing custom COIN-OR Package")
set(CUSTCOIN_DIR "${ZERO_DEP_DIR}/CustomCoin")

include_directories(${COIN_INCLUDE_DIR})
include_directories(${CUSTCOIN_DIR})
include_directories(${GUROBI_INCLUDE_DIRS})


add_library(libCustomCoin STATIC ${CUSTCOIN_DIR}/OsiGrbSolverInterface.cpp
            #${CUSTCOIN_DIR}/CglKnapsackCoverZERO.cpp
            )

if (NOT WIN32 OR MINGW)
    set_target_properties(libCustomCoin PROPERTIES PREFIX "")
endif ()

set_target_properties(libCustomCoin PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${CONDA_ENV}/lib"
                      ARCHIVE_OUTPUT_DIRECTORY "${CONDA_ENV}/lib")
configure_file(${CUSTCOIN_DIR}/OsiGrbSolverInterface.hpp ${CONDA_ENV}/include/coin/OsiGrbSolverInterface.hpp COPYONLY)