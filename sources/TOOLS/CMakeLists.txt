set( HEADERS
	OGS_WDC.h
	OGS_contraflow.h
	fancyTimer.h
)

set( SOURCES
	OGS_contraflow.cpp
	OGS_WDC.cpp
	logger.cpp
)

include_directories(
        ${CMAKE_SOURCE_DIR}/Base
        ${CMAKE_SOURCE_DIR}/FEM
        ${CMAKE_SOURCE_DIR}/FileIO
        ${CMAKE_SOURCE_DIR}/GEO
        ${CMAKE_SOURCE_DIR}/LIB
        ${CMAKE_SOURCE_DIR}/MathLib
        ${CMAKE_SOURCE_DIR}/MSH
        ${CMAKE_SOURCE_DIR}/TOOLS
        ${CMAKE_SOURCE_DIR}/../Libs/WellDoubletControl
        ${CMAKE_SOURCE_DIR}/../Libs/Contraflow
)


# Create the library
add_library( TOOLS STATIC ${HEADERS} ${SOURCES} )


target_link_libraries(TOOLS
	FEM
	Base
	FileIO
	GCC
	GEO
	MSH
)
