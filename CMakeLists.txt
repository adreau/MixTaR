cmake_minimum_required (VERSION 2.8)
project (MixTaR)

set(CMAKE_CXX_FLAGS "-O3")
#set(BOOST_ROOT External_Libraries/boost_1_54_0/)
include_directories(
  External_Libraries/gatb-core/gatb-core/src/
  External_Libraries/gatb-core/src/gatb-core-build/include
  External_Libraries/sparsehash/src/sparsehash/include
  External_Libraries/gatb-core/gatb-core/thirdparty

  External_Libraries/seqAn/src/seqAn/include/
  )
find_library(GATB_LIBRARY gatbcore HINTS ${LIB_DIR}/build/probt-prefix/src/probt/lib)

set(MIXTAR_SRC 
  src/main.cpp 
  src/smallAssemblyFind.cpp 
  src/smallAssemblyOverlap.cpp)

# define the binary name to compile
add_executable(mixtar ${MIXTAR_SRC})
# link to specified library
target_link_libraries(mixtar gatbcore)
# install binary in specified directory
install(TARGETS mixtar RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
