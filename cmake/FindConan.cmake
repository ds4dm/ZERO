# ##############################################################################
# CONAN
# ##############################################################################
# TODO: temporary workaround
set(CONAN_GIT
    "https://raw.githubusercontent.com/dforsten/cmake-conan/be876ddc9a45401aa842ad45e3de75e597c6d263/conan.cmake"
)
set(CONAN_SHA1 "f87634521cc40038bed27198bf8f627da6e4c9ed")
if(NOT (EXISTS "${CMAKE_BINARY_DIR}/conan.cmake"))
  set(CONAN_DOWNLOAD True)
else()
  file(SHA1 "${CMAKE_BINARY_DIR}/conan.cmake" LOCAL_SHA1)
  if(NOT (${LOCAL_SHA1} STREQUAL ${CONAN_SHA1}))
    set(CONAN_DOWNLOAD True)
    message("Updating Conan (SH1 mismatch)")
  endif()
endif()
if(CONAN_DOWNLOAD)
  message(STATUS "Downloading conan.cmake")
  file(
    DOWNLOAD ${CONAN_GIT} "${CMAKE_BINARY_DIR}/conan.cmake"
    TLS_VERIFY ON
    EXPECTED_HASH SHA1=${CONAN_SHA1})
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)
