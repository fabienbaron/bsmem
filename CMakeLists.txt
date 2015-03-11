cmake_minimum_required(VERSION 2.8)
project(bsmem C Fortran)

# Set this to the directory in which you want BSMEM installed if the default
# location isn't acceptable.
#SET(CMAKE_INSTALL_PREFIX ~/local)

#=============================================================================
# All done editing
SET(BSMEM_VERSION "1.6")

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

add_subdirectory(src)
