############################
### Find OGS directories ###
############################

if(DEFINED BENCHMARK_DIR)
	find_path (BENCHMARK_DIR_FOUND copy.py ${BENCHMARK_DIR})
else()
	find_path (BENCHMARK_DIR_FOUND copy.py ${PROJECT_SOURCE_DIR}/../benchmarks)
endif()

if(DEFINED EXAMPLEDATA_DIR)
	find_path (EXAMPLEDATA_DIR_FOUND points.gli ${EXAMPLEDATA_DIR})
else()
	find_path (EXAMPLEDATA_DIR_FOUND points.gli ${PROJECT_SOURCE_DIR}/../ExampleData)
endif()

if(DEFINED TESTDATA_DIR)
	find_path(TESTDATA_DIR_FOUND testdata.dummy ${TESTDATA_DIR})
else()
	find_path(TESTDATA_DIR_FOUND testdata.dummy ${PROJECT_SOURCE_DIR}/../testdata)
endif()

######################
### Find tools     ###
######################

# Find Python interpreter
find_package (PythonInterp)

# Find Subversion
find_package(Subversion)

# Find Git
find_package(Git)

# msysGit on Windows
if(WIN32 AND GIT_FOUND)
	find_package(MsysGit)
endif() # WIN32 AND GIT_FOUND

# Find dot tool from graphviz
find_program(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
find_package(Doxygen)

# Find gnu profiler gprof
find_program(GPROF_PATH gprof DOC "GNU profiler gprof")

find_package(cppcheck)

# Find Exuberant ctags or BBEdit for code completion
find_program(CTAGS_TOOL_PATH ctags DOC "Exuberant ctags")
find_program(BBEDIT_TOOL_PATH bbedit DOC "BBEdit Editor")
if(BBEDIT_TOOL_PATH)
	add_custom_target(ctags
		bbedit --maketags
		WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
		COMMENT "Creating tags..." VERBATIM
	)
	add_custom_command(TARGET ctags POST_BUILD
		COMMAND mv -f tags ../tags
		WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
		COMMENT "Moving tags..." VERBATIM
	)
else()
	if(CTAGS_TOOL_PATH)
		add_custom_target(ctags
			ctags -R --fields=+iamS -f ${CMAKE_SOURCES_DIR}/../tags
			WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
			COMMENT "Creating tags..." VERBATIM
		)
	endif()
endif()

## Unix tools ##
# Date
find_program(DATE_TOOL_PATH date PATHS ${MSYSGIT_BIN_DIR})
# Grep
find_program(GREP_TOOL_PATH grep PATHS ${MSYSGIT_BIN_DIR})
# Unzip
find_program(UNZIP_TOOL_PATH unzip PATHS ${MSYSGIT_BIN_DIR})

# Hide these variables for the CMake user
mark_as_advanced(DOT_TOOL_PATH GPROF_PATH CTAGS_TOOL_PATH BBEDIT_TOOL_PATH
	UNZIP_TOOL_PATH
)
########################
### Find other stuff ###
########################

# Check if on Jenkins
if(NOT $ENV{JENKINS_URL} STREQUAL "")
	set(JENKINS_URL $ENV{JENKINS_URL})
	set(JENKINS_JOB_NAME $ENV{JOB_NAME})
endif()


######################
### Find libraries ###
######################
if(OGS_FEM_PETSC OR OGS_NO_EXTERNAL_LIBS)
	return()
endif()

find_path (OGS_LIBS_DIR_FOUND geotiff.lib ${PROJECT_SOURCE_DIR}/../Libs/libgeotiff)

# Find precompiled libraries (for BRNS GEMS LIS)
find_path (OGS_PRECOMPILED_LIBS_DIR_FOUND BrnsDll.lib ${PROJECT_SOURCE_DIR}/../Libs/precompiled)
if (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	include_directories (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
	link_directories (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
else ()
	if (WIN32)
		if (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
			message (FATAL_ERROR "Precompiled libraries not found! Make sure to also check out the trunk/Libs directory beneath your sources directory.")
		endif ()
	else ()
		if (OGS_FEM_LIS)
			message (FATAL_ERROR "Precompiled libraries not found! Make sure to also check out the trunk/Libs directory beneath your sources directory.")
		endif ()
	endif ()
endif ()


## pthread ##
set ( CMAKE_THREAD_PREFER_PTHREAD ON CACHE BOOL "" )
find_package( Threads )
if ( CMAKE_USE_PTHREADS_INIT AND NOT HAVE_PTHREADS)
	set (HAVE_PTHREADS TRUE CACHE BOOL "Is PThreads found.")
	message (STATUS "pthread library found." )
endif ()
if(HAVE_PTHREADS)
  add_definitions(-DHAVE_PTHREADS)
endif()
mark_as_advanced(CMAKE_THREAD_PREFER_PTHREAD)

## boost (see FindBoost.cmake for more options) ##
##kg44 this configuration works for boost and petsc on a cray
option(Boost_USE_STATIC_LIBS "" ON)
option(Boost_USE_MULTITHREADED "" ON)
option(Boost_USE_STATIC_RUNTIME "" ON)
mark_as_advanced(Boost_USE_STATIC_LIBS Boost_USE_MULTITHREADED Boost_USE_STATIC_RUNTIME)

if(NOT OGS_FEM_GEMS AND NOT OGS_FEM_PETSC_GEMS)
	if(NOT OGS_DONT_USE_BOOST)
		find_package( Boost 1.50.0 COMPONENTS filesystem system regex)
	endif()
else()
	# Boost with threads is required for GEMS
	find_package(Boost 1.50.0 COMPONENTS system thread REQUIRED)
        message(STATUS "** Boost root: ${BOOST_ROOT}")
        message(STATUS "** Boost include: ${Boost_INCLUDE_DIR}")
        message(STATUS "** Boost libraries: ${Boost_LIBRARY_DIRS}")
        message(STATUS "** Boost libraries: ${Boost_LIBRARIES}")
endif()

if(OGS_FEM_MKL)
	# Find MKLlib
	find_package( MKL REQUIRED )
	include_directories (${MKL_INCLUDE_DIR})
endif()

if(OGS_FEM_LIS OR OGS_FEM_MKL)
	# Find LISlib
	find_package( LIS REQUIRED )
	set (NEW_EQS ON)
	add_definitions(
		-o3
		-DIPMGEMPLUGIN
	)
endif()

# Find OpenMP
if(PARALLEL_USE_OPENMP)
	find_package( OpenMP REQUIRED )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
endif()
