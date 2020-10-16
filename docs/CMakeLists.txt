cmake_minimum_required(VERSION 3.14)

include(../cmake/ZERO.cmake)
set(CMAKE_CXX_STANDARD 14)
project(
  ZERODoc
  VERSION ${ZERO_V}
  LANGUAGES CXX
  DESCRIPTION ${ZERO_DESC})

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release ... FORCE)
endif()
option(CXX "enable C++ compilation" ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH ../cmake)
include(${CMAKE_MODULE_PATH}/FindSphinx.cmake)

find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_SOURCE_DIR}/html/)
set(SPHINX_INDEX_FILE ${CMAKE_CURRENT_SOURCE_DIR}/index.html)
file(REMOVE_RECURSE SPHINX_BUILD ${CMAKE_CURRENT_SOURCE_DIR}/api/)
file(REMOVE_RECURSE SPHINX_BUILD ${CMAKE_CURRENT_SOURCE_DIR}/html/)

# Custom resources
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/html/)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/html/html)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/support_files/res
     DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/html/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/support_files/index.html
     DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/html/html/)
# Nojekyll for gh-pages
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/html/.nojekyll " ")

file(GLOB_RECURSE ZERO_SOURCE_LIST "../src/*.cpp" "../include/*.h")

set(ZERO_SOURCE_LIST_ESCAPED "")
foreach(FILE ${ZERO_SOURCE_LIST})
  set(ZERO_SOURCE_LIST_ESCAPED "${ZERO_SOURCE_LIST_ESCAPED} ${FILE}")
endforeach()

set(SPHINX_IN ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in)
set(SPHINX_OUT ${CMAKE_CURRENT_SOURCE_DIR}/conf.py)

configure_file(${SPHINX_IN} ${SPHINX_OUT} @ONLY)

add_custom_command(
  OUTPUT ${SPHINX_INDEX_FILE}
  COMMAND ${SPHINX_EXECUTABLE} -b html ${SPHINX_SOURCE} ${SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
  MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
  COMMENT "Generating documentation with Sphinx")

# Nice named target so we can run the job easily
add_custom_target(ZERODocs ALL DEPENDS ${SPHINX_INDEX_FILE})

# Add an install target to install the docs
include(GNUInstallDirs)
install(DIRECTORY ${SPHINX_BUILD} DESTINATION ${CMAKE_INSTALL_DOCDIR})

message("\n\nZERO Docs loaded")
