execute_process(COMMAND "conda --version" RESULT_VARIABLE CONDA_RET)
string(REPLACE "conda " "" CONDA_RET EQUAL CONDA_RET EQUAL)
string(REPLACE ". " "" CONDA_RET EQUAL CONDA_RET EQUAL)
if (CONDA_RET_EQUAL GREATER_EQUAL 3)
    message("Found Conda")
    set(CONDA_FOUND TRUE)
elseif ()
    message("Conda not found.")
endif ()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Conda DEFAULT_MSG)