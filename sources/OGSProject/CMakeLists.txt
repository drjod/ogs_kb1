# Source files
set( SOURCES
	gs_project.cpp
)


# Header files
set( HEADERS
	gs_project.h
)

# and include additional header
include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/FEM
)


# Create the library
add_library( OGSProject
	${SOURCES}
	${HEADERS}
)

target_link_libraries( OGSProject
	FEM
)
