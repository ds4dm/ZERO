# ##############################################################################
# CONAN
# ##############################################################################
# TODO: temporary workaround to fix 0.15 bugs on mac
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

find_program(CONAN_CMD conan)
if(NOT CONAN_CMD)
  message(FATAL_ERROR "Conan executable not found!")
endif()
execute_process(COMMAND ${CONAN_CMD} --version
                OUTPUT_VARIABLE CONAN_VERSION_OUTPUT)
string(REGEX MATCH ".*Conan version ([0-9]+\.[0-9]+\.[0-9]+)" FOO
             "${CONAN_VERSION_OUTPUT}")
set(CONAN_VERSION_REQUIRED 1.0.0)
if(${CMAKE_MATCH_1} VERSION_LESS ${CONAN_VERSION_REQUIRED})
  message(
    FATAL_ERROR
      "Conan outdated. Installed: ${CONAN_VERSION}, \
        required: ${CONAN_VERSION_REQUIRED}. Consider updating via 'pip \
        install conan --upgrade'.")
endif()

execute_process(COMMAND ${CONAN_CMD} remote list OUTPUT_VARIABLE CONAN_REMOTES)

if(NOT ("${CONAN_REMOTES}" MATCHES ".*darcamo-bintray:.*"))
  message(STATUS "FindConan.cmake: adding the repository dracamo-bintray")
  execute_process(COMMAND ${CONAN_CMD} remote add -i 1 darcamo-bintray
                          https://api.bintray.com/conan/darcamo/cppsim)
endif()
