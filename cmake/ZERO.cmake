set(ZERO_AUTHORS "Gabriele Dragotto, Sriram Sankaranarayanan")
set(ZERO_V 1.0.0)
set(ZERO_DESC "ZERO is a multi-purpose game solver")

# CXX Options
option(CXX "enable C++ compilation" ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release ... FORCE)
endif()
