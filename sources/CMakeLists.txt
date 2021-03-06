#### OGS-5 Project ####

# Specify minimum CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3)

# Project name
project(OGS)

###############
### Options ###
###############
# Configurations
option(OGS_FEM "Default configuration" OFF)
option(OGS_FEM_SP "Sparse matrix configuration" OFF)
option(OGS_FEM_JFNK "Jacobi free Newton-Krylov method for H2M problem" OFF)
option(OGS_FEM_PETSC "Use PETSc parallel solver" OFF)
option(OGS_FEM_PETSC_GEMS "GEMS with PETSc parallel solver" OFF)
option(OGS_FEM_GEMS "GEMS configuration" OFF)
option(OGS_FEM_IPQC "IPhreeqc configuration" OFF)
option(OGS_FEM_PQC "Phreeqc configuration" OFF)
option(OGS_FEM_BRNS "BRNS configuration" OFF)
option(OGS_FEM_CHEMAPP "CHEMAPP configuration" OFF)
option(OGS_FEM_LIS "Library of Iterative Solvers for Linear Systems configuration" OFF)
option(OGS_FEM_CAP "CAP configuration" OFF)
option(OGS_FEM_MKL "Math kernel library configuration" OFF)
option(BLUE_G "Blue/G optimization" OFF)

# supported on Linux and Windows
option(OGS_FEM_MPI "Message passing interface configuration" OFF)
# needs to be set of MPI and GEMS is used
option(PARALLEL_USE_MPI "Message passing interface configuration also for GEMS coupling" OFF)
# needs to be set for parallel mode of KinReact module
OPTION(OGS_FEM_MPI_KRC "KinReact mpi config." OFF) 
# needs to be set for SUNDIALS CVODE solver library
OPTION(OGS_KRC_CVODE "SUNDIALS CVODE solver library configuration" OFF)

### CMake setup ###
include(scripts/cmake/CMakeSetup.cmake)

### Compiler setup ###
include(scripts/cmake/CompilerSetup.cmake)

### Find directories and libraries ###
include (scripts/cmake/Find.cmake)

### General project setup ###
include(scripts/cmake/GeneralProjectSetup.cmake)

# Set default configuration when no other config is given
if (NOT (OGS_FEM OR OGS_FEM_SP OR OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS OR OGS_FEM_PQC OR OGS_FEM_IPQC OR OGS_FEM_BRNS OR OGS_FEM_CHEMAPP OR OGS_FEM_LIS OR OGS_FEM_MKL OR OGS_FEM_MPI OR OGS_FEM_PETSC OR OGS_FEM_CAP))
	message (STATUS "No configuration specified. Default confuguration is used.")
	set (OGS_FEM ON)
endif ()

# Check if a valid OGS configuration is given (see Macros.cmake)
CHECK_CONFIG()

# Additional options
option(OGS_PACKAGING "Creating installers / packages" OFF)
OPTION_REQUIRES(OGS_PACKAGING_ZIP "Do you want to package as zip?" OGS_PACKAGING)
option(OGS_PYTHON "Enable python scripting interface" OFF)
option(OGS_BUILD_UTILITIES "Build additional utility programs?" OFF)
option(OGS_CMAKE_DEBUG "Show additional cmake debug information" OFF)
option(OGS_ONE_BIG_GTEST "Should all gtests be collected to one ctest" ON)
option(OGS_BUILD_INFO "Should build information be generated" ON)
option(CMAKE_CMD_ARGS "Build ogs5 storing the CMake command line args (hint: use cmake.ogs.sh" OFF )
if(GCC)
	if(GPROF_PATH)
		option(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
		if(OGS_PROFILE AND DOT_TOOL_PATH)
			option(OGS_OUTPUT_PROFILE "Enables pdf output of profiling results." OFF)
		endif() # OGS_PROFILE AND DOT_TOOL_PATH
	endif() # GPROF_PATH

	OPTION_REQUIRES(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OGS_BUILD_TESTS)
endif() # GCC
option (OGS_BUILD_TESTS "Enables building of tests." OFF)
option (OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies."
OFF)
option (OGS_DELETE_EDGES_AFTER_INIT "Delete mesh edges after initialization if possible" OFF)
option (OGS_OUTPUT_TIME_ALLPCS "Output clock time for all processes" OFF)

mark_as_advanced(FORCE OGS_PYTHON OGS_CMAKE_DEBUG OGS_ONE_BIG_GTEST OGS_BUILD_INFO CMAKE_CMD_ARGS OGS_DELETE_EDGES_AFTER_INIT)

###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
set (OGS_VERSION_MAJOR 0)
set (OGS_VERSION_MINOR 2)
set (OGS_VERSION_PATCH 5)
#set (OGS_RELEASE_PERSONS "")
set (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}")
#set (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}(${OGS_RELEASE_PERSONS})")
# should we use ISO dates: 2010-12-31 ?
set (OGS_DATE "21.04.2020")

# Check for number of processors
include(ProcessorCount)
ProcessorCount(PROCESSOR_COUNT)
if(PROCESSOR_COUNT EQUAL 0)
	message(WARNING "Processor count could not be detected. Setting to one processor.")
	set(PROCESSOR_COUNT 1)
else()
	message(STATUS "Number of processors: ${PROCESSOR_COUNT}")
endif() # PROCESSOR_COUNT

###################
### Fem configs ###
###################
if(OGS_FEM)
	message (STATUS	"Configuring for FEM command line" )
endif()

if(OGS_FEM_SP)
	message (STATUS	"Configuring for FEM command line with new sparse matrix" )
	set(OGS_USE_NEW_EQS ON BOOL)
endif()

if(OGS_FEM_JFNK)
	message (STATUS "Configuring for Jacobi free Newton-Krylov method for H2M problem" )
	set(OGS_USE_NEW_EQS ON BOOL)
	add_definitions(-DJFNK_H2M)
endif()

if(OGS_FEM_PETSC)
	message (STATUS "Configuring for PETSc" )
	set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake/cmake/findPETSC)

	find_package(PETSc REQUIRED QUIET)
	if(PETSC_VERSION VERSION_GREATER 3.3)
		add_definitions(-DUSEPETSC34 -DUSE_PETSC)
		message(STATUS "found version greater 3.3, version is ${PETSC_VERSION}")
	else()
		add_definitions(-DUSE_PETSC)
	endif()
	include_directories( ${PETSC_INCLUDES} )
endif()

if(OGS_FEM_PETSC_GEMS)
	message (STATUS "Configuring for GEMS with PETSc" )
	include_directories( ${Boost_INCLUDE_DIR} )
	link_directories( ${Boost_LIBRARY_DIR} )
	set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake/cmake/findPETSC)

	find_package(PETSc REQUIRED QUIET)
	if(PETSC_VERSION VERSION_GREATER 3.3)
		add_definitions(-DUSEPETSC34 -DUSE_PETSC -DNO_ERROR_CONTROL -DGEM_REACT -DGEMlib -DIPMGEMPLUGIN)
		message(STATUS "found version greater 3.3, version is ${PETSC_VERSION}")
	else()
		add_definitions(-DUSE_PETSC -DNO_ERROR_CONTROL -DGEM_REACT -DGEMlib -DIPMGEMPLUGIN)
	endif()
	include_directories( ${PETSC_INCLUDES} )

	if (MSVC)
		set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
		set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	endif ()
endif()

if(OGS_FEM_MPI)
	message (STATUS "Configuring for FEM command line with MPI" )
	set(PARALLEL_USE_MPI ON BOOL "Use Message Passing Interface")
	set(OGS_USE_NEW_EQS ON BOOL)
endif()

if(OGS_FEM_MPI_KRC)
	MESSAGE (STATUS "Configuring for FEM command line with MPI_KRC" )
	add_definitions(-DUSE_MPI_KRC)
	#set(PARALLEL_USE_MPI ON BOOL "Use Message Passing Interface")
endif()

if(OGS_FEM_GEMS)
	message (STATUS	"Configuring for FEM command line with GEMS" )
	set(OGS_USE_NEW_EQS ON BOOL)
	add_definitions(-DNO_ERROR_CONTROL -DGEM_REACT -DGEMlib -DIPMGEMPLUGIN )
	if (MSVC)
		set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
		set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	endif ()

	include_directories( ${Boost_INCLUDE_DIR} )
	link_directories( ${Boost_LIBRARY_DIR} )

	if(PARALLEL_USE_MPI)
		add_definitions(-DUSE_MPI_GEMS)
	endif()
endif()

if(OGS_FEM_IPQC)
	message (STATUS "Configuring for FEM command line with IPQC" )
	add_definitions( -DOGS_FEM_IPQC )
	if(PARALLEL_USE_MPI)
		set(OGS_USE_NEW_EQS ON BOOL)
	endif()
endif()

if(OGS_FEM_PQC)
	message (STATUS "Configuring for FEM command line with PQC" )
	message (STATUS "Configuration: g++ and LIBPHREEQC")
	set(OGS_USE_NEW_EQS ON BOOL)
	add_definitions(-DLIBPHREEQC -D_POSIX_SOURCE -DNO_ERROR_CONTROL)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		add_definitions ( "-DMDL_DEBUG" )
	endif() # Debug
endif()

if(OGS_FEM_BRNS)
	message (STATUS	"Configuring for FEM command line with BRNS" )
	if (UNIX)
		find_package(BLAS)
		find_package(LAPACK)
	endif()
	set(OGS_USE_NEW_EQS ON BOOL)
	add_definitions(-DBRNS -DNO_ERROR_CONTROL -D_POSIX_SOURCE)
	if(PARALLEL_USE_MPI)
		add_definitions(-DUSE_MPI_BRNS)
	endif()
endif()

if(OGS_FEM_LIS)
	message (STATUS "Configuring for FEM command line with LIS" )
	set(OGS_USE_NEW_EQS ON BOOL)
	add_definitions(-DIPMGEMPLUGIN -DLIS)
	set(PARALLEL_USE_OPENMP ON BOOL "Use Library of Iterative Solvers")
	set(LIS ON BOOL "Use Library of Iterative Solvers")
	mark_as_advanced(PARALLEL_USE_OPENMP LIS)
endif()

if(OGS_FEM_MKL)
	message (STATUS	"Configuring for FEM command line with MKL" )
	set(OGS_USE_NEW_EQS ON BOOL)
	add_definitions(-DNO_ERROR_CONTROL -DMKL -DLIS)
	set(PARALLEL_USE_OPENMP ON BOOL "Use Library of Iterative Solvers")
	set(LIS ON BOOL "Use Library of Iterative Solvers")
	set(MKL ON BOOL "Use Math Kernel Library")
	mark_as_advanced(PARALLEL_USE_OPENMP LIS MKL)
endif()

if(OGS_FEM_CHEMAPP)
	message (STATUS	"Configuring for FEM command line with CHEMAPP" )
	add_definitions(-DCHEMAPP)
endif()

if(OGS_FEM_CAP)
	message (STATUS  "Configuring for FEM command line with CAP" )
	set( PRJ_EXT "-CAP" )
	add_definitions(-DOGS_FEM_CAP)
	if (UNIX) # Only supported on Linux
		include_directories(/usr/lib64)
		link_directories(/usr/lib64)
	endif ()
endif()

### End fem configs ###

IF(OGS_KRC_CVODE)
	MESSAGE (STATUS	"Configuring with SUNDIALS CVODE solver lib for KinReact " )
	ADD_DEFINITIONS(-DOGS_KRC_CVODE)
ENDIF(OGS_KRC_CVODE)

if(OGS_USE_NEW_EQS)
	message(STATUS "Use New EQS class")
	add_definitions(-DNEW_EQS)
endif()

if(PARALLEL_USE_MPI)
	message(STATUS "Enabling MPI")
	find_package(MPI REQUIRED)
	add_definitions(-DUSE_MPI)
	
    #IF(OGS_FEM_MPI_KRC)
   	#	ADD_DEFINITIONS(-DUSE_MPI_KRC)
    #ENDIF(OGS_FEM_MPI_KRC)
		
endif()

if(OGS_FEM_MPI_KRC)
	message(STATUS "Enabling MPI for KRC")
	find_package(MPI REQUIRED)
#	add_definitions(-DUSE_MPI)
 	add_definitions(-DUSE_MPI_KRC)
endif()


if(MSVC)
	add_definitions(-DMSVC)
endif() # MSVC

if(OGS_DELETE_EDGES_AFTER_INIT)
	add_definitions(-DOGS_DELETE_EDGES_AFTER_INIT)
endif()

if(OGS_OUTPUT_TIME_ALLPCS)
	add_definitions(-DOUTPUT_TIME_FOR_ALLPCS)
endif()

# Add subdirectories with the projects
add_subdirectory( ThirdParty )
include_directories (SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty)

add_subdirectory( Base )
add_subdirectory( MathLib )
add_subdirectory( GEO )
add_subdirectory( MSH )
add_subdirectory( MSHGEOTOOLS )
add_subdirectory( FEM )
add_subdirectory( GCC )
add_subdirectory( TOOLS )

if(OGS_FEM_CHEMAPP)
	add_subdirectory( EQL )
	link_directories( ${CMAKE_SOURCE_DIR}/EQL )
endif()
if (OGS_FEM_CAP)
	add_subdirectory( ThirdParty/CMP )
endif ()
add_subdirectory( FileIO )
add_subdirectory( OGS )


# Create Configure.h
configure_file (Base/Configure.h.in ${PROJECT_BINARY_DIR}/Base/Configure.h)
configure_file (Base/BuildInfo.h.in ${PROJECT_BINARY_DIR}/Base/BuildInfo.h)

if (BENCHMARK_DIR_FOUND OR GTEST_FOUND)
	enable_testing()
endif ()

if (EXAMPLEDATA_DIR_FOUND)
	install (DIRECTORY ${EXAMPLEDATA_DIR_FOUND} DESTINATION . PATTERN .svn EXCLUDE)
endif ()

if (OGS_BUILD_TESTS)
	add_subdirectory(tests)
endif ()

if (OGS_PACKAGING)
	# Add packaging configuration
	message (STATUS "Enabling packaging")
	include( scripts/cmake/Pack.cmake)
endif ()

if (OGS_CMAKE_DEBUG)
	include(ListAllCMakeVariableValues)
	list_all_cmake_variable_values()
endif ()

if (BENCHMARK_DIR_FOUND)
	set(BENCHMARK_TIMEOUT 1800) # in s, 30 minutes timeout on normal benchmarks
	set(EXCEEDING_BENCHMARK_TIMEOUT 86400) # 1 day timeout on exceeding benchmarks
	string (COMPARE EQUAL "$ENV{HOSTNAME}" "dev2.intern.ufz.de" HOST_IS_DEV2)
	if (${HOST_IS_DEV2})
		message (STATUS "On dev2: file comparing enabled")
	endif ()

	file (GLOB BENCHMARK_CONFIGS "${PROJECT_SOURCE_DIR}/../benchmarks/*.cmake")
	foreach (BENCHMARK_CONFIG ${BENCHMARK_CONFIGS})
		include ("${BENCHMARK_CONFIG}")
	endforeach()

	unset (COPY_BENCHMARKS_TO_REF CACHE)
endif ()

if (OGS_BUILD_UTILITIES)
	add_subdirectory (UTL/MSHGEOTOOLS/)
	add_subdirectory (UTL/GIS2FEM/)
endif ()

## Documentation ##
if(DOXYGEN_FOUND)
	option(DOCS_GENERATE_DIAGRAMS "Use the DOT tool to generate class diagrams." OFF)
	option(DOCS_GENERATE_CALL_GRAPHS "Generate call dependency graphs." OFF)
	option(DOCS_GENERATE_COLLABORATION_GRAPHS "Generate collaboration graphs." OFF)
	if(DOT_TOOL_PATH AND DOCS_GENERATE_DIAGRAMS)
		set(DOT_FOUND "YES" CACHE INTERNAL "")
	endif() # DOT_TOOL_PATH AND DOCS_GENERATE_DIAGRAMS
	if(DOCS_GENERATE_CALL_GRAPHS)
		set(DOCS_GENERATE_CALL_GRAPHS_STRING "YES" CACHE INTERNAL "")
	endif() # DOCS_GENERATE_CALL_GRAPHS
	if(DOCS_GENERATE_COLLABORATION_GRAPHS)
		set(DOCS_GENERATE_COLLABORATION_GRAPHS_STRING "YES" CACHE INTERNAL "")
	endif() # DOCS_GENERATE_COLLABORATION_GRAPHS
	get_filename_component(DOT_TOOL_PATH_ONLY ${DOT_TOOL_PATH} PATH)
	configure_file(scripts/doxygen/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile)
	add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating source code documentation with Doxygen." VERBATIM)
endif() # DOXYGEN_FOUND
