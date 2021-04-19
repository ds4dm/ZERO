message(Installing OsiGrb)
set(OSIGRB_DIR "${ZERO_DEP_DIR}/OsiGrb")

include_directories(${CONDA_INCLUDE_DIR}/coin)
include_directories(${OSIGRB_DIR})
include_directories(${GUROBI_INCLUDE_DIRS})


add_library(libOsiGrb STATIC ${OSIGRB_DIR}/OsiGrbSolverInterface.cpp)

if (NOT WIN32 OR MINGW)
    set_target_properties(libOsiGrb PROPERTIES PREFIX "")
endif ()

set_target_properties(libOsiGrb PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${CONDA_ENV}/lib"
                      ARCHIVE_OUTPUT_DIRECTORY "${CONDA_ENV}/lib")
configure_file(${OSIGRB_DIR}/OsiGrbSolverInterface.hpp ${CONDA_ENV}/include/coin/OsiGrbSolverInterface.hpp COPYONLY)