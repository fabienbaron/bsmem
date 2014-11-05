# - Try to find the NFFT
# Once done this will define
#  NFFT_FOUND - System has NFFT
#  NFFT_INCLUDE_DIRS - The NFFT include directories
#  NFFT_LIBRARIES - The libraries needed to use NFFT (none)

FIND_LIBRARY(NFFT_LIBRARY 
    NAMES nfft3
    HINTS ${NFFT_HINT}/lib
    DOC "NFFT library.")

FIND_PATH(NFFT_INCLUDE_DIR 
    NAMES nfft3.h
    HINTS ${NFFT_HINT}/include
    DOC "NFFT include directory.")

set(NFFT_LIBRARIES ${NFFT_LIBRARY})
set(NFFT_INCLUDE_DIRS ${NFFT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NFFT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NFFT DEFAULT_MSG
                                  NFFT_LIBRARY NFFT_INCLUDE_DIR)

mark_as_advanced(NFFT_INCLUDE_DIR NFFT_LIBRARY)
