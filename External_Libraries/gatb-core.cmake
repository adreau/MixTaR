# gatb-core build
set(gatb_prefix ${CMAKE_HOME_DIRECTORY}/gatb-core)
set(CMAKE_INSTALL_PREFIX ${CMAKE_HOME_DIRECTORY}/)

# Add the external project.
ExternalProject_Add(gatb-core
  PREFIX ${gatb_prefix}
  SOURCE_DIR ${gatb_prefix}/gatb-core
  INSTALL_DIR ${CMAKE_HOME_DIRECTORY}
  INSTALL_COMMAND ""
  )