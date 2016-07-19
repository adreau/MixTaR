cmake_minimum_required (VERSION 2.8)
project (MixTaR)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

include_directories(
  External_Libraries/gatb-core/gatb-core/src/
  External_Libraries/gatb-core/src/gatb-core-build/include
  External_Libraries/sparsehash/src/sparsehash/include
  External_Libraries/gatb-core/gatb-core/thirdparty
  External_Libraries/seqAn/src/seqAn/include/
  )

set(MIXTAR_SRC 
  src/main.cpp)
#  src/smallAssemblyFind.cpp
#  src/smallAssemblyOverlap.cpp)

# define library directories
link_directories(../External_Libraries/gatb-core/src/gatb-core-build/lib/)

# define the binary name to compile
add_executable(mixtar ${MIXTAR_SRC})

# link to specified library
target_link_libraries(mixtar gatbcore hdf5 dl z pthread)

# install binary in specified directory
install(TARGETS mixtar RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})
