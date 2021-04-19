message("Loading PATH")
set(ENV{PATH_LICENSE_STRING}
    "2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0"
)

set(PATH_URL "http://pages.cs.wisc.edu/~ferris/path/julia")

if(WIN32)
  set(PATH_FILES "lusol.dll" "path50.dll")
endif()
if(UNIX)
  set(PATH_FILES "liblusol.so" "libpath50.so")
  if(APPLE)
    set(PATH_FILES ${PATH_FILES} "liblusol.dylib" "libpath50.dylib")
  endif()
endif()

foreach(lib ${PATH_FILES})
  if(NOT (EXISTS "${PATH_DIR}/lib/${lib}"))
    file(DOWNLOAD ${PATH_URL}/${lib} "${PATH_DIR}/lib/${lib}" TLS_VERIFY ON)
  endif()
endforeach()

find_library(
  PATH_LIBRARY
  NAMES path path47 path50
  HINTS ${PATH_DIR}/lib
  PATH_SUFFIXES lib)

find_library(
  LUSOL_LIBRARY
  NAMES lusol
  HINTS ${PATH_DIR}/lib
  PATH_SUFFIXES lib)
message("PathLIB is located at ${PATH_LIBRARY}")
set(PATH_INCLUDE_DIR "${PATH_DIR}/include")

find_package_handle_standard_args(PATH DEFAULT_MSG LUSOL_LIBRARY)
find_package_handle_standard_args(PATH DEFAULT_MSG PATH_LIBRARY)
