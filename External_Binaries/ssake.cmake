# ssake build

set(ssake_PREFIX "${CMAKE_SOURCE_DIR}/ssake")
set(ssake_INSTALL_DIR "${CMAKE_SOURCE_DIR}/ssake")

message("ssake_PREFIX='${ssake_PREFIX}'")
message("ssake_INSTALL_DIR='${ssake_INSTALL_DIR}'")

ExternalProject_Add(ssake
  URL "http://www.bcgsc.ca/platform/bioinfo/software/ssake/releases/3.8.2/ssake_v3-8-2.tar.gz"
  DOWNLOAD_DIR ${ssake_PREFIX}
  SOURCE_DIR ${ssake_PREFIX}
  BINARY_DIR ${ssake_PREFIX}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

#set(ssake_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/ssake/include")
#set(ssake_LIBRARIES "${CMAKE_SHARED_LIBRARY_PREFIX}ssake${CMAKE_SHARED_LIBRARY_SUFFIX}")
#include_directories(${ssake_INCLUDE_DIRS})

