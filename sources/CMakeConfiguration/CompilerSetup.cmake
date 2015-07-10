INCLUDE(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
INCLUDE(SetDefaultBuildType)
INCLUDE(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Release)
INCLUDE(MSVCMultipleProcessCompile) # /MP Switch for VS

IF (WIN32)
	## For Visual Studio compiler
	IF (MSVC)
		ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS)
		# Sets warning level 3 and ignores some warnings
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267")
		SET(GCC OFF)

		DisableCompilerFlag(DEBUG /RTC1)

    # Set $PATH to Visual Studio bin directory. Needed for finding dumpbin.exe
    IF (MSVC80)
      SET(ENV{PATH} "$ENV{PATH};$ENV{VS80COMNTOOLS}..\\..\\VC\\bin")
    ENDIF ()
    IF (MSVC90)
      SET(ENV{PATH} "$ENV{PATH};$ENV{VS90COMNTOOLS}..\\..\\VC\\bin")
    ENDIF ()
    IF (MSVC10)
      SET(ENV{PATH} "$ENV{PATH};$ENV{VS100COMNTOOLS}..\\..\\VC\\bin")
    ENDIF ()

	ELSE (MSVC)
#FOR CYGWIN.  25.02.2010. WW
		MESSAGE (STATUS "Might be GCC under cygwin.")
		SET(GCC ON)
#		MESSAGE (FATAL_ERROR "Aborting: On Windows only the Visual Studio compiler is supported!")
	ENDIF (MSVC)
ENDIF (WIN32)

### For GNU C/CXX. WW
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	SET(GCC ON)
	IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
		IF (OGS_FEM_PETSC_GEMS)
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
		ELSE(OGS_FEM_PETSC_GEMS)
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
		ENDIF(OGS_FEM_PETSC_GEMS)
	ENDIF(NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
	# -g
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra -fno-nonansi-builtins")

	EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	STRING(REPLACE "\n" "" GCC_VERSION ${GCC_VERSION})
	MESSAGE(STATUS "GCC_VERSION: ${GCC_VERSION}")
	IF (NOT (GCC_VERSION VERSION_LESS 4.8) ) 
	  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs") # suppress warnings in Eigen
	ENDIF()
  
	# would be cool: -Woverloaded-virtual, would be overkill: -Weffc++
	ADD_DEFINITIONS(-DGCC)

	IF (OGS_PROFILE)
		IF( NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
			MESSAGE(STATUS "When using profiling you should set CMAKE_BUILD_TYPE to Release.")
		ENDIF()
		SET(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
	ENDIF (OGS_PROFILE)
ENDIF() # CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC

IF(OGS_COVERAGE)
  INCLUDE(CodeCoverage)
ENDIF()
