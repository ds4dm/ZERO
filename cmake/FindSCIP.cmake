find_path(SCIP_INCLUDE_DIRS scip/scip.h HINTS ${SCIP_ROOT_DIR}/include PATH_SUFFIXES src)
find_library(
        SCIP_LIBRARIES
        NAMES scip libscip
        HINTS ${SCIP_ROOT_DIR}/lib
        PATH_SUFFIXES lib lib/static)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCIP DEFAULT_MSG SCIP_LIBRARIES SCIP_INCLUDE_DIRS)