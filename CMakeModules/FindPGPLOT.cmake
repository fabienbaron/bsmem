# - Try to find the PGPLOT
# Once done this will define
#  PGPLOT_FOUND - System has PGPLOT
#  PGPLOT_INCLUDE_DIRS - The PGPLOT include directories
#  PGPLOT_LIBRARIES - The libraries needed to use PGPLOT (none)

FIND_LIBRARY(PGPLOT_LIBRARY
    NAMES pgplot
    HINTS ${PGPLOT_HINT}/lib ${CMAKE_INSTALL_PREFIX}/lib ${PGPLOT_HINT}/lib64 ${CMAKE_INSTALL_PREFIX}/lib64 /usr/lib/lib64
    DOC "PGPLOT library.")

FIND_LIBRARY(CPGPLOT_LIBRARY 
    NAMES cpgplot
    HINTS ${CPGPLOT_HINT}/lib ${CMAKE_INSTALL_PREFIX}/lib ${CPGPLOT_HINT}/lib64 ${CMAKE_INSTALL_PREFIX}/lib64 /usr/lib/lib64
    DOC "CPGPLOT library.")

FIND_PATH(PGPLOT_INCLUDE_DIR 
    NAMES cpgplot.h
    HINTS ${PGPLOT_HINT}/include ${CMAKE_INSTALL_PREFIX}/include
    DOC "PGPLOT include directory.")

# concatenate the two libraries into one variable:
set(PGPLOT_LIBRARIES ${PGPLOT_LIBRARY} ${CPGPLOT_LIBRARY})
set(PGPLOT_INCLUDE_DIRS ${PGPLOT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PGPLOT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PGPLOT DEFAULT_MSG
                                  PGPLOT_LIBRARIES PGPLOT_INCLUDE_DIR)

mark_as_advanced(PGPLOT_INCLUDE_DIR PGPLOT_LIBRARY)
