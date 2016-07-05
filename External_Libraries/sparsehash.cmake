# sparsehash build

set(sparsehash_PREFIX "${CMAKE_SOURCE_DIR}/sparsehash")
set(sparsehash_INSTALL_DIR "${CMAKE_SOURCE_DIR}")

message("sparsehash_PREFIX='${sparsehash_PREFIX}'")
message("sparsehash_INSTALL_DIR='${sparsehash_INSTALL_DIR}'")

ExternalProject_Add(sparsehash
  PREFIX ${CMAKE_SOURCE_DIR}/sparsehash
  GIT_REPOSITORY "https://github.com/sparsehash/sparsehash/"
  GIT_TAG "sparsehash-2.0.2"
  INSTALL_DIR ${sparsehash_INSTALL_DIR}

  CONFIGURE_COMMAND ${sparsehash_PREFIX}/src/sparsehash/configure --prefix=${sparsehash_PREFIX}/src/sparsehash
  BUILD_COMMAND make
  INSTALL_COMMAND make install
)

#set(sparsehash_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/sparsehash/include")
#set(sparsehash_LIBRARIES "${CMAKE_SHARED_LIBRARY_PREFIX}sparsehash${CMAKE_SHARED_LIBRARY_SUFFIX}")
#include_directories(${sparsehash_INCLUDE_DIRS})

