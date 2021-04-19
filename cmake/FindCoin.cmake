message("Searching for Coin")
unset(COIN_CGL)
unset(COIN_UTILS)
find_library(
        COIN_CGL
        NAMES Cgl
        HINTS ${COIN_HINT} ${CONDA_ENV}
        PATH_SUFFIXES lib)

find_library(
        COIN_UTILS
        NAMES CoinUtils
        HINTS ${COIN_HINT} ${CONDA_ENV}
        PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Coin DEFAULT_MSG COIN_CGL)
find_package_handle_standard_args(Coin DEFAULT_MSG COIN_UTILS)
