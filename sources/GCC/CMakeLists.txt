# Source files 
set( HEADERS
#	Activity.h
	Density.h	
	DiffAL.h
	DiffOH.h
	dllheader.h
	eqbrm.h
	export.h
	Fluid.h
	HKF.h
	IAPWS-IF97.h
	IO.h
	logk.h
	NR.h
	PDE.h
	PITZcalc.h
	PITZdata.h
	slop98.h
	species.h
	VLE.h
	WATEQdata.h
)

set ( SOURCES
#	Activity.cpp
	Density.cpp	
	DiffAL.cpp
	DiffOH.cpp
	eqbrm.cpp
	export.cpp
	Fluid.cpp
	HKF.cpp
	IAPWS-IF97.cpp
	IO.cpp
	logk.cpp
	NR.cpp
	PDE.cpp
	PITZcalc.cpp
	PITZdata.cpp
	VLE.cpp
)

# Create the library
add_library( GCC STATIC 
	     ${SOURCES}
             ${HEADERS} 
)
include_directories(
        ${CMAKE_SOURCE_DIR}/GCC
)


##SET_TARGET_PROPERTIES(GCC PROPERTIES LINKER_LANGUAGE CXX)

IF (WIN32)
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN")
ENDIF (WIN32)

IF (UNIX)
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN -D__unix")
ENDIF (UNIX)
