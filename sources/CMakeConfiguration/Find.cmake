############################
### Find OGS directories ###
############################

IF(DEFINED BENCHMARK_DIR)
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py ${BENCHMARK_DIR})
ELSE()
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py ${PROJECT_SOURCE_DIR}/../benchmarks)
ENDIF()

IF(DEFINED EXAMPLEDATA_DIR)
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${EXAMPLEDATA_DIR})
ELSE()
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${PROJECT_SOURCE_DIR}/../ExampleData)
ENDIF()

IF(DEFINED TESTDATA_DIR)
	FIND_PATH(TESTDATA_DIR_FOUND testdata.dummy ${TESTDATA_DIR})
ELSE()
	FIND_PATH(TESTDATA_DIR_FOUND testdata.dummy ${PROJECT_SOURCE_DIR}/../testdata)
ENDIF()

######################
### Find tools     ###
######################

# Find Python interpreter
FIND_PACKAGE (PythonInterp)

# Find Subversion
FIND_PACKAGE(Subversion)

# Find Git
FIND_PACKAGE(Git)

# msysGit on Windows
IF(WIN32 AND GIT_FOUND)
	FIND_PACKAGE(MsysGit)
ENDIF() # WIN32 AND GIT_FOUND

# Find dot tool from graphviz
FIND_PROGRAM(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
FIND_PACKAGE(Doxygen)

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof")

FIND_PACKAGE(cppcheck)

# Find Exuberant ctags or BBEdit for code completion
FIND_PROGRAM(CTAGS_TOOL_PATH ctags DOC "Exuberant ctags")
FIND_PROGRAM(BBEDIT_TOOL_PATH bbedit DOC "BBEdit Editor")
IF(BBEDIT_TOOL_PATH)
	ADD_CUSTOM_TARGET(ctags
		bbedit --maketags
		WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
		COMMENT "Creating tags..." VERBATIM
	)
	ADD_CUSTOM_COMMAND(TARGET ctags POST_BUILD
		COMMAND mv -f tags ../tags
		WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
		COMMENT "Moving tags..." VERBATIM
	)
ELSE()
	IF(CTAGS_TOOL_PATH)
		ADD_CUSTOM_TARGET(ctags
			ctags -R --fields=+iamS -f ${CMAKE_SOURCES_DIR}/../tags
			WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
			COMMENT "Creating tags..." VERBATIM
		)
	ENDIF()
ENDIF()

## Unix tools ##
# Date
FIND_PROGRAM(DATE_TOOL_PATH date PATHS ${MSYSGIT_BIN_DIR})
# Grep
FIND_PROGRAM(GREP_TOOL_PATH grep PATHS ${MSYSGIT_BIN_DIR})
# Unzip
FIND_PROGRAM(UNZIP_TOOL_PATH unzip PATHS ${MSYSGIT_BIN_DIR})

# Hide these variables for the CMake user
MARK_AS_ADVANCED(DOT_TOOL_PATH GPROF_PATH CTAGS_TOOL_PATH BBEDIT_TOOL_PATH
	UNZIP_TOOL_PATH
)
########################
### Find other stuff ###
########################

# Check if on Jenkins
IF(NOT $ENV{JENKINS_URL} STREQUAL "")
	SET(JENKINS_URL $ENV{JENKINS_URL})
	SET(JENKINS_JOB_NAME $ENV{JOB_NAME})
ENDIF()


######################
### Find libraries ###
######################
IF(OGS_FEM_PETSC OR OGS_NO_EXTERNAL_LIBS)
	RETURN()
ENDIF()

FIND_PATH (OGS_LIBS_DIR_FOUND geotiff.lib ${PROJECT_SOURCE_DIR}/../Libs/libgeotiff)

# Find precompiled libraries (for BRNS GEMS LIS)
FIND_PATH (OGS_PRECOMPILED_LIBS_DIR_FOUND BrnsDll.lib ${PROJECT_SOURCE_DIR}/../Libs/precompiled)
IF (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
	LINK_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
ELSE (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	IF (WIN32)
		IF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
			MESSAGE (FATAL_ERROR "Precompiled libraries not found! Make sure to also check out the trunk/Libs directory beneath your sources directory.")
		ENDIF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
	ELSE (WIN32)
		IF (OGS_FEM_LIS)
			MESSAGE (FATAL_ERROR "Precompiled libraries not found! Make sure to also check out the trunk/Libs directory beneath your sources directory.")
		ENDIF (OGS_FEM_LIS)
	ENDIF (WIN32)
ENDIF (OGS_PRECOMPILED_LIBS_DIR_FOUND)


FIND_PACKAGE( Shapelib )
IF(Shapelib_FOUND)
	ADD_DEFINITIONS(-DShapelib_FOUND)
ENDIF() # Shapelib_FOUND

## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD ON CACHE BOOL "" )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT AND NOT HAVE_PTHREADS)
	SET (HAVE_PTHREADS TRUE CACHE BOOL "Is PThreads found.")
	MESSAGE (STATUS "pthread library found." )
ENDIF ()
IF(HAVE_PTHREADS)
  ADD_DEFINITIONS(-DHAVE_PTHREADS)
ENDIF()
MARK_AS_ADVANCED(CMAKE_THREAD_PREFER_PTHREAD)

## boost (see FindBoost.cmake for more options) ##
##kg44 this configuration works for boost and petsc on a cray
OPTION(Boost_USE_STATIC_LIBS "" ON)
OPTION(Boost_USE_MULTITHREADED "" ON)
OPTION(Boost_USE_STATIC_RUNTIME "" ON)
MARK_AS_ADVANCED(Boost_USE_STATIC_LIBS Boost_USE_MULTITHREADED Boost_USE_STATIC_RUNTIME)

IF(NOT OGS_FEM_GEMS AND NOT OGS_FEM_PETSC_GEMS)
	IF(NOT OGS_DONT_USE_BOOST)
		FIND_PACKAGE( Boost 1.50.0 COMPONENTS filesystem system regex)
	ENDIF()
ELSE()
	# Boost with threads is required for GEMS
	FIND_PACKAGE(Boost 1.50.0 COMPONENTS system thread REQUIRED)
        MESSAGE(STATUS "** Boost root: ${BOOST_ROOT}")
        MESSAGE(STATUS "** Boost include: ${Boost_INCLUDE_DIR}")
        MESSAGE(STATUS "** Boost libraries: ${Boost_LIBRARY_DIRS}")
        MESSAGE(STATUS "** Boost libraries: ${Boost_LIBRARIES}")
ENDIF()

## VTK ##
IF (OGS_LIBS_DIR_FOUND AND NOT VTK_DIR)
	SET (VTK_DIR ${OGS_LIBS_DIR}/VTK/build)
ENDIF () # OGS_LIBS_DIR_FOUND
IF(NOT OGS_DONT_USE_VTK)
	FIND_PACKAGE( VTK )
ENDIF()
IF(VTK_FOUND)
	ADD_DEFINITIONS(-DVTK_FOUND)
	FIND_PACKAGE(QVTK)
	IF(NOT QVTK_FOUND AND OGS_USE_QT)
		MESSAGE(FATAL_ERROR "QVTK was not found but is required for OGS_BUILD_GUI! On Ubuntu it can be installed via 'sudo apt-get install libvtk5-qt4-dev'")
	ENDIF()
ENDIF()

## NetCDF ##
IF("${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_PATCH_VERSION}" VERSION_GREATER 5.6)
	FIND_PATH(VTK_NETCDF_FOUND netcdf.h
		PATHS ${VTK_INCLUDE_DIRS}/vtknetcdf ${VTK_SOURCE_DIR}/Utilities/vtknetcdf
		PATH_SUFFIXES include
		NO_DEFAULT_PATH)
ENDIF()

IF(VTK_NETCDF_FOUND)
	ADD_DEFINITIONS(-DVTK_NETCDF_FOUND)
	INCLUDE_DIRECTORIES(
		${VTK_NETCDF_FOUND} ${VTK_DIR}/Utilities ${VTK_NETCDF_FOUND}/..
		${VTK_NETCDF_FOUND}/../.. ${VTK_NETCDF_FOUND}/../cxx)
ELSE()
	SET(NETCDF_CXX TRUE)
	FIND_PACKAGE(NetCDF)
	IF(NOT NETCDF_FOUND AND OGS_USE_QT)
		MESSAGE(FATAL_ERROR "NetCDF was not found but is required for OGS_BUILD_GUI!")
	ENDIF()
ENDIF()

## geotiff ##
IF(NOT MSVC)
	FIND_PACKAGE( LibTiff )
ENDIF() # NOT MSVC
FIND_PACKAGE( LibGeoTiff )
IF(libgeotiff_FOUND)
	ADD_DEFINITIONS(-Dlibgeotiff_FOUND)
ENDIF() # libgeotiff_FOUND

IF (OGS_PYTHON)
	FIND_PACKAGE (PythonLibs 2.5 REQUIRED)
ENDIF (OGS_PYTHON)

## Qt4 library ##
IF(NOT OGS_DONT_USE_QT)
	FIND_PACKAGE( Qt4 4.5)
ENDIF(NOT OGS_DONT_USE_QT)

IF ( QT4_FOUND )
	# OPTION(OGS_GUI OFF )
	# this is needed to correctly link the qt libraries through target_link_libraries
	# By default only QtCore and QtGui modules are enabled
	# other modules must be enabled like this:
	SET( QT_USE_QTOPENGL TRUE )
	SET( QT_USE_QTSQL TRUE )
	SET( QT_USE_QTTEST TRUE )
	SET( QT_USE_QTXML TRUE )
	IF (QT_QTXMLPATTERNS_FOUND)
		set( QT_USE_QTXMLPATTERNS TRUE )
	ENDIF (QT_QTXMLPATTERNS_FOUND)
	INCLUDE( ${QT_USE_FILE} )
	ADD_DEFINITIONS(${QT_DEFINITIONS})
ENDIF (QT4_FOUND )

## VRPN ##
FIND_PACKAGE( VRPN )

## VRED ##
FIND_PATH (VRED_DIR_FOUND vrNodePtr.h ${VRED_DIR}/include/vred)

## OpenSG ##
FIND_PACKAGE( OpenSG COMPONENTS OSGBase OSGSystem)

IF(OGS_FEM_MKL)
	# Find MKLlib
	FIND_PACKAGE( MKL REQUIRED )
	INCLUDE_DIRECTORIES (${MKL_INCLUDE_DIR})
ENDIF()

IF(OGS_FEM_LIS OR OGS_FEM_MKL)
	# Find LISlib
	FIND_PACKAGE( LIS REQUIRED )
	set (NEW_EQS ON)
	add_definitions(
		-o3
		-DIPMGEMPLUGIN
	)
ENDIF()

# Find OpenMP
IF(PARALLEL_USE_OPENMP)
	FIND_PACKAGE( OpenMP REQUIRED )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
ENDIF(PARALLEL_USE_OPENMP)

IF(PARALLEL_USE_MPI)
	IF (WIN32)
#		MESSAGE (FATAL_ERROR "Aborting: MPI is only supported under UNIX/LINUX!")
		#ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
		FIND_PACKAGE(MPI REQUIRED)
	ENDIF(WIN32)
	IF(UNIX)

# If there is an mpi compiler find it and interogate (farther below) it for the include
# and lib dirs otherwise we will continue to search from ${_MPI_BASE_DIR}.

		IF( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
			find_program(MPI_COMPILER
				NAMES mpic++ mpicxx mpiCC mpicc
				HINTS "${_MPI_BASE_DIR}"
				PATH_SUFFIXES bin
				DOC "MPI compiler. Used only to detect MPI compilation flags.")
			IF(MPI_COMPILER)

			MESSAGE (STATUS  "CMake version is less than 2.8, MPI compiler is set directly" )
			mark_as_advanced(MPI_COMPILER)
				SET(CMAKE_C_COMPILER ${MPI_COMPILER})
				SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
			ENDIF(MPI_COMPILER)
		ELSE( ${CMAKE_MAJOR_VERSION}  EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
			FIND_PACKAGE(MPI REQUIRED)
		ENDIF( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
	ENDIF(UNIX)
ENDIF(PARALLEL_USE_MPI)

