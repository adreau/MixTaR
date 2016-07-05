cmake_minimum_required (VERSION 2.8)
project (MixTaR_External_Libraries)

set(LIB_DIR "${PROJECT_SOURCE_DIR}/External_Libraries")

include(ExternalProject)
include(sparsehash.cmake)
include(gatb-core-install.cmake)
include(gatb-core.cmake)
include(seqAn.cmake)
