cmake_minimum_required(VERSION 2.8) 

#=============================================================================
# Copyright 2012 Brian Kloppenborg
#
# This code is licensed under the MIT License.  
#
# The MIT License
#
# License for the specific language governing rights and limitations under
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#=============================================================================

# Find required packages
FIND_PACKAGE(CFITSIO REQUIRED)
FIND_PACKAGE(FFTW3 REQUIRED)
FIND_PACKAGE(PGPLOT REQUIRED)
FIND_PACKAGE(NFFT REQUIRED)
FIND_PACKAGE(X11 REQUIRED)
# Set includes
INCLUDE_DIRECTORIES(${CFITSIO_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${PGPLOT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${NFFT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${X11_INCLUDE_DIRS})

# Enable optimization.
#LIST(APPEND CMAKE_Fortran_FLAGS "-O5")
#LIST(APPEND CMAKE_C_FLAGS "-march=native -O2 -pipe")

LIST(APPEND CMAKE_C_FLAGS "-ggdb -fomit-frame-pointer -pipe")
# Collect the source:
FILE(GLOB SOURCE *.f *.c)

add_executable(bsmem ${SOURCE})
SET_TARGET_PROPERTIES(bsmem PROPERTIES LINKER_LANGUAGE C)
target_link_libraries(bsmem ${CFITSIO_LIBRARIES} ${FFTW3_LIBRARIES} ${PGPLOT_LIBRARIES} ${NFFT_LIBRARIES} ${X11_LIBRARIES})


# Set install directory
set_property(TARGET bsmem PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
INSTALL(TARGETS bsmem DESTINATION bin)

