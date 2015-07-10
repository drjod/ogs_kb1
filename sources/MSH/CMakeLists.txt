# Source files

set( SOURCES_FEM
	admin.h
	GridAdapter.h
	msh_core.h
	msh_edge.h
	msh_elem.h
	msh_faces.h
	msh_lib.h
	msh_mesh.h
	msh_node.h
	MshEditor.h
	MSHEnums.h
	MeshNodesAlongPolyline.h
	
	GridAdapter.cpp
	msh_core.cpp
	msh_edge.cpp
	msh_elem.cpp
	msh_faces.cpp
	#msh_gen.cpp
	msh_lib.cpp
	#msh_mat.cpp
	msh_mesh.cpp
	msh_node.cpp
	MshEditor.cpp
	MSHEnums.cpp

	MeshNodesAlongPolyline.cpp
	#msh_pnt.cpp		# references removed elements.h nodes.h
	
	
	#NW/delaunay2D.h
	#NW/delaunay3D.h
	#NW/delaunay2D.cpp
	#NW/delaunay3D.cpp
	
	#TM/dtm_calculate.h
	#TM/dtm_crowd.h
	#TM/dtm_element.h
	#TM/dtm_error.h
	#TM/dtm_figure.h
	#TM/dtm_fixed.h
	#TM/dtm_gpnode.h
	#TM/dtm_gptetra.h
	#TM/dtm_gptriangle.h
	#TM/dtm_group.h
	#TM/dtm_laplas.h
	#TM/dtm_neighbor.h
	#TM/dtm_node.h
	#TM/dtm_plane.h
	#TM/dtm_point.h
	#TM/dtm_stdio.h
	#TM/dtm_surface.h
	#TM/dtm_tetgen.h
	#TM/dtm_tetra.h
	#TM/dtm_timer.h
	#TM/dtm_triangle.h
	#TM/dtmesh.h
	#TM/dtm_calculate.cpp
	#TM/dtm_crowd.cpp
	#TM/dtm_error.cpp
	#TM/dtm_figure.cpp
	#TM/dtm_gpnode.cpp
	#TM/dtm_gptetra.cpp
	#TM/dtm_gptriangle.cpp
	#TM/dtm_laplas.cpp
	#TM/dtm_neighbor.cpp
	#TM/dtm_node.cpp
	#TM/dtm_plane.cpp
	#TM/dtm_point.cpp
	#TM/dtm_stdio.cpp
	#TM/dtm_surface.cpp
	#TM/dtm_tetgen.cpp
	#TM/dtm_tetra.cpp
	#TM/dtm_triangle.cpp
	#TM/dtmesh.cpp
	#TM/main.cpp
)

if(OGS_USE_QT) #By WW
   set( SOURCES 
        ${SOURCES_FEM}
		MeshQualityChecker.h
		MeshQualityShortestLongestRatio.h
		MeshQualityEquiAngleSkew.h
		MeshQualityChecker.cpp
		MeshQualityShortestLongestRatio.cpp
		MeshQualityArea.cpp
		MeshQualityVolume.cpp
		MeshQualityEquiAngleSkew.cpp
        
		Delaunay/DT_Matrix.h
		Delaunay/tricommon.h
		Delaunay/tri_algebra.cpp
		Delaunay/tri_geometry.cpp
		Delaunay/tri_main.cpp
		Delaunay/tri_triangulation.cpp
		Delaunay/tri_utility.cpp        
        )
else (OGS_USE_QT)
		set( SOURCES ${SOURCES_FEM})
endif (OGS_USE_QT)

IF(OGS_FEM_PETSC OR OGS_FEM_PETSC_GEMS) #By WW
   SET( SOURCES 
        ${SOURCES_FEM}
       msh_lib_ext.cpp  
      )
ENDIF(OGS_FEM_PETSC OR OGS_FEM_PETSC_GEMS) 



include_directories(
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/LIB
	${CMAKE_SOURCE_DIR}/OGSProject
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

# Create the library
add_library( MSH STATIC ${SOURCES} )

target_link_libraries (
	MSH
	FEM
	FileIO
	GEO
	${CMAKE_THREAD_LIBS_INIT}
)
