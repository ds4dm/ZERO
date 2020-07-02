message("Loading Google OR Tools")
include(FetchContent)
FetchContent_Declare(
  or-tools
  GIT_REPOSITORY https://github.com/google/or-tools.git
  GIT_TAG master)

FetchContent_MakeAvailable(or-tools)
