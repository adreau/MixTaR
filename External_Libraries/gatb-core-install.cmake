# download the project from the git repository

execute_process(
    COMMAND git clone --branch v1.1.0  "https://github.com/GATB/gatb-core/"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
