# Source files
SET( SOURCES
	Writer.cpp
	GMSInterface.cpp
	GocadInterface.cpp
	OGSIOVer4.cpp
	PetrelInterface.cpp	
	StationIO.cpp
	readNonBlankLineFromInputStream.cpp
	#FEMIO.cpp
	FEMIO/GeoIO.cpp
	FEMIO/ProcessIO.cpp
	FEMIO/BoundaryConditionIO.cpp
	GeoIO/Gmsh2GeoIO.cpp
	MeshIO/GMSHInterface.cpp
	MeshIO/GMSHFixedMeshDensity.cpp
	MeshIO/GMSHAdaptiveMeshDensity.cpp
	MeshIO/GMSHPoint.cpp
	MeshIO/GMSHLine.cpp
	MeshIO/GMSHLineLoop.cpp
	MeshIO/GMSHPolygonTree.cpp
	MeshIO/LegacyVtkInterface.cpp
	MeshIO/OGSMeshIO.cpp
	MeshIO/TetGenInterface.cpp
	XmlIO/RapidXMLInterface.cpp
)

# Header files
SET( HEADERS
	Writer.h
	GMSInterface.h
	GocadInterface.h
	OGSIOVer4.h
	PetrelInterface.h
	StationIO.h
	readNonBlankLineFromInputStream.h
	#FEMIO.h
	FEMIO/GeoIO.h
	FEMIO/ProcessIO.h
	FEMIO/BoundaryConditionIO.h
	GeoIO/Gmsh2GeoIO.h
	MathIO/CRSIO.h
	MeshIO/GMSHInterface.h
	MeshIO/GMSHMeshDensityStrategy.h
	MeshIO/GMSHNoMeshDensity.h
	MeshIO/GMSHFixedMeshDensity.h
	MeshIO/GMSHAdaptiveMeshDensity.h
	MeshIO/GMSHPoint.h
	MeshIO/GMSHLine.h
	MeshIO/GMSHLineLoop.h
	MeshIO/GMSHPolygonTree.h
	MeshIO/LegacyVtkInterface.h
	MeshIO/OGSMeshIO.h
	MeshIO/TetGenInterface.h
	XmlIO/RapidXMLInterface.h
)

IF (OGS_USE_QT)
	IF (VTK_FOUND)
		INCLUDE( ${VTK_USE_FILE} )
	ENDIF ()
	IF (NETCDF_FOUND OR VTK_NETCDF_FOUND)
		SET (HEADERS NetCDFInterface.h ${HEADERS})
		SET (SOURCES NetCDFInterface.cpp ${SOURCES})
	ENDIF()
ENDIF ()

IF (Shapelib_FOUND)
	SET( HEADERS ${HEADERS} SHPInterface.h )
	SET (SOURCES ${SOURCES} SHPInterface.cpp )
	INCLUDE_DIRECTORIES (${Shapelib_INCLUDE_DIR})
ENDIF (Shapelib_FOUND)

IF (QT4_FOUND)
	SET( HEADERS ${HEADERS} 
		XmlIO/XMLInterface.h 
		XmlIO/XmlGmlInterface.h
		XmlIO/XmlStnInterface.h
		XmlIO/XmlCndInterface.h
		XmlIO/XmlGspInterface.h
		XmlIO/XmlLutReader.h
		FEFLOWInterface.h
	)
	SET (SOURCES ${SOURCES} 
		XmlIO/XMLInterface.cpp
		XmlIO/XmlGmlInterface.cpp
		XmlIO/XmlStnInterface.cpp
		XmlIO/XmlCndInterface.cpp
		XmlIO/XmlGspInterface.cpp
		FEFLOWInterface.cpp 
	)
ENDIF (QT4_FOUND)

IF(OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS)
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

# Create the library
ADD_LIBRARY( FileIO STATIC
	${SOURCES}
	${HEADERS}
)

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/FileIO/MathIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/FileIO/GeoIO
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/GEM
	${CMAKE_SOURCE_DIR}/Qt/Base
)

IF (QT4_FOUND)
	INCLUDE_DIRECTORIES(
		../Qt/VtkVis
	)
ENDIF()

TARGET_LINK_LIBRARIES( FileIO
	Base
	#	MathLib
	GEO
	FEM
)

IF (Shapelib_FOUND)
	TARGET_LINK_LIBRARIES( FileIO ${Shapelib_LIBRARIES} )
ENDIF (Shapelib_FOUND)

IF(OGS_USE_QT AND VTK_FOUND)
	TARGET_LINK_LIBRARIES( FileIO vtkCommon )
ENDIF()

IF(OGS_PACKAGING AND OGS_USE_QT)
	# Install xml schema files
	FILE (GLOB XML_SCHEMA_FILES . "*.xs?")
	INSTALL (FILES ${XML_SCHEMA_FILES} DESTINATION bin)
ENDIF()
