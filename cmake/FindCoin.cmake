message("Searching for Coin")
unset(COIN_CGL)
unset(COIN_UTILS)
find_library(
        COIN_CGL
        NAMES Cgl
        HINTS ${COIN_HINT}
        PATH_SUFFIXES lib)

find_library(
        COIN_OSI
        NAMES Osi
        HINTS ${COIN_HINT}
        PATH_SUFFIXES lib)

find_library(
        COIN_UTILS
        NAMES CoinUtils
        HINTS ${COIN_HINT}
        PATH_SUFFIXES lib)
#find_path(COIN_INCLUDE_DIR CoinUtility.hpp HINTS ${COIN_HINT}/include PATH_SUFFIXES coin-or)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Coin DEFAULT_MSG COIN_CGL)
find_package_handle_standard_args(Coin DEFAULT_MSG COIN_OSI)
find_package_handle_standard_args(Coin DEFAULT_MSG COIN_UTILS)
#find_package_handle_standard_args(Coin DEFAULT_MSG COIN_INCLUDE_DIR)
