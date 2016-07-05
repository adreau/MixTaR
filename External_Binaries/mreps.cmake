# mreps build

set(mreps_PREFIX "${CMAKE_SOURCE_DIR}/mreps")
set(mreps_INSTALL_DIR "${CMAKE_SOURCE_DIR}")

message("mreps_PREFIX='${mreps_PREFIX}'")
message("mreps_INSTALL_DIR='${mreps_INSTALL_DIR}'")

ExternalProject_Add(mreps
  PREFIX ${mreps_PREFIX}
  GIT_REPOSITORY "https://github.com/gregorykucherov/mreps"
  INSTALL_DIR ${mreps_INSTALL_DIR}

#  CONFIGURE_COMMAND ${mreps_PREFIX}/src/mreps/configure --prefix=${mreps_PREFIX/src/mreps}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND cd ${mreps_PREFIX}/src/mreps && make 
  INSTALL_COMMAND ""
)

#set(mreps_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/mreps/include")
#set(mreps_LIBRARIES "${CMAKE_SHARED_LIBRARY_PREFIX}mreps${CMAKE_SHARED_LIBRARY_SUFFIX}")
#include_directories(${mreps_INCLUDE_DIRS})

