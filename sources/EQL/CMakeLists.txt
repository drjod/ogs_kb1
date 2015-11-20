# Source files
set ( SOURCES 
	cacint.h
	cacint.cpp
	eqlink.h
	eqlink.cpp
)

# Create the library
add_library ( EQL STATIC ${SOURCES} )

include_directories (
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/Mathlib
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/FEM
)
