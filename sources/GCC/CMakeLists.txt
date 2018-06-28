set( HEADERS
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
add_library( GCC STATIC ${HEADERS} ${SOURCES} )
include_directories( ${CMAKE_SOURCE_DIR}/GCC 
${CMAKE_SOURCE_DIR}/../Libs/WellDoubletControl

)

if (WIN32)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN")
endif ()

if (UNIX)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN -D__unix")
endif ()
