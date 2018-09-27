set( HEADERS
	OGS_WDC.h
	fancyTimer.h
)

set( SOURCES
	OGS_WDC.cpp
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
