# seqAn build

set(seqAn_PREFIX "${CMAKE_SOURCE_DIR}/seqAn")
set(seqAn_INSTALL_DIR "${CMAKE_SOURCE_DIR}")

message("seqAn_PREFIX='${seqAn_PREFIX}'")
message("seqAn_INSTALL_DIR='${seqAn_INSTALL_DIR}'")

find_package(OpenMP)
find_library(BOOST_LIBRARY spl HINTS ${LIB_DIR}/boost_1_54_0/lib)

ExternalProject_Add(seqAn
  PREFIX ${CMAKE_SOURCE_DIR}/seqAn
  GIT_REPOSITORY "https://github.com/seqan/seqan"
  GIT_TAG "seqan-v2.0.0"
  INSTALL_DIR ${seqAn_INSTALL_DIR}

#  CONFIGURE_COMMAND ${seqAn_PREFIX}/src/seqAn/configure --prefix=${seqAn_PREFIX/src/seqAn}
  BUILD_COMMAND make
  INSTALL_COMMAND ""
)

#set(seqAn_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/seqAn/include")
#set(seqAn_LIBRARIES "${CMAKE_SHARED_LIBRARY_PREFIX}seqAn${CMAKE_SHARED_LIBRARY_SUFFIX}")
#include_directories(${seqAn_INCLUDE_DIRS})

