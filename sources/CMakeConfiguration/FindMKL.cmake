# Find Intel Math Karnel Library
# 
# This sets the following variables:
#
#   MKL_INCLUDES   - Directory containing MKL include files
#   MKL_LIBRARIES  - Path to MKL libraries
#
# Users may wish to set:
#  MKL_DIR         - Location to start searching for MKL libraries
#

SET(MKL_DIR "${MKL_DIR}" CACHE PATH "MKL root diretory")

FIND_PATH(MKL_INCLUDES NAMES mkl.h
    HINTS ${MKL_DIR} PATH_SUFFIXES include
)

# list up MKL library names
IF (${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
	#64 bit
	set(MKL_LIB_NAMES mkl_intel_lp64 mkl_core)
	set(MKL_PATH_SUFFIXES "lib/intel64" "lib/em64t")
ELSE()
	#32 bit
	set(MKL_LIB_NAMES mkl_core mkl_intel mkl_solver)
	set(MKL_PATH_SUFFIXES "lib/32")
ENDIF()

IF (PARALLEL_USE_OPENMP)
	IF(CMAKE_C_COMPILER EQUAL "icc")
        list(APPEND MKL_LIB_NAMES mkl_intel_thread iomp5 pthread)
	ELSE()
        list(APPEND MKL_LIB_NAMES mkl_gnu_thread)
	ENDIF()
ELSE()
        list(APPEND MKL_LIB_NAMES mkl_sequential)
ENDIF()
#MESSAGE (STATUS "MKL_LIB_NAMES ${MKL_LIB_NAMES}")

# find the libraries
foreach (lib_name ${MKL_LIB_NAMES})
	find_library(${lib_name}_LIBRARY ${lib_name} HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
    IF(${lib_name}_LIBRARY)
        LIST(APPEND MKL_LIBRARIES ${${lib_name}_LIBRARY})
    ENDIF()
    #MESSAGE (STATUS "${lib_name} - ${${lib_name}_LIBRARY}")
endforeach()
MESSAGE (STATUS "MKL libraries found: ${MKL_LIBRARIES}")

# 
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDES MKL_LIBRARIES)
IF(MKL_FOUND)
	mark_as_advanced(MKL_DIR MKL_INCLUDES MKL_LIBRARIES)
ENDIF()

