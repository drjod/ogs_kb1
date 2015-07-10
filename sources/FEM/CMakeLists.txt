# Source files

set( SOURCES
	Stiff_Bulirsch-Stoer.h
    StepperBulischStoer.h
	BoundaryCondition.h
	break.h
	Constrained.h
	conversion_rate.h
	display.h
	eos.h
	fem_ele.h
	fem_ele_std.h
	fem_ele_vec.h
	FEMCondition.h
	FEMEnums.h
	femlib.h
	files0.h
	InitialCondition.h
	makros.h
	mathlib.h
	matrix_class.h
	SparseMatrixDOK.h
	memory.h
	pcs_dm.h
	problem.h
	prototyp.h
	rf_bc_new.h
	rf_fct.h
	rf_fluid_momentum.h
	rf_ic_new.h
	rf_kinreact.h
	rf_mfp_new.h
	rf_mmp_new.h
	rf_msp_new.h
	rf_node.h
	rf_num_new.h
	rf_out_new.h
	rf_pcs.h
	rf_random_walk.h
	rf_react.h
	rf_react_int.h
	rf_st_new.h
	rf_tim_new.h
	rfmat_cp.h
	solver.h
	SourceTerm.h
	timer.h
	tools.h
	vtk.h
	Eclipse.h
	DUMUX.h
	GeoInfo.h
	ProcessInfo.h
	DistributionInfo.h
	LinearFunctionData.h
	Output.h
	GeoInfo.cpp
	ProcessInfo.cpp
	DistributionInfo.cpp
	LinearFunctionData.cpp
	Stiff_Bulirsch-Stoer.cpp
	BoundaryCondition.cpp
	break.cpp
	conversion_rate.cpp
	display.cpp
	eos.cpp
	fem_ele.cpp
	fem_ele_std.cpp
	fem_ele_std1.cpp
	fem_ele_vec.cpp
	femlib.cpp
	FEMCondition.cpp
	FEMEnums.cpp
	files0.cpp
	InitialCondition.cpp
	mathlib.cpp
	matrix_class.cpp
    SparseMatrixDOK.cpp
	memory.cpp
	pcs_dm.cpp
	problem.cpp
	rf_bc_new.cpp
	rf_fct.cpp
	rf_fluid_momentum.cpp
	rf_ic_new.cpp
	rf_kinreact.cpp
	rf_mfp_new.cpp
	rf_mmp_new.cpp
	rf_msp_new.cpp
	rf_node.cpp
	rf_num_new.cpp
	rf_out_new.cpp
	rf_pcs.cpp
	rf_random_walk.cpp
	rf_react.cpp
	rf_react_int.cpp
	rf_st_new.cpp
	rf_tim_new.cpp
	rfmat_cp.cpp	
	SourceTerm.cpp
	timer.cpp
	tools.cpp
	vtk.cpp
	Eclipse.cpp
	DUMUX.cpp
	Output.cpp
    CAP_IO.h 
    CAP_IO.cpp
    rf_react_cap.h 
    rf_react_cap.cpp 	
)

IF(OGS_FEM_PETSC OR OGS_FEM_PETSC_GEMS)
 SET( SOURCES ${SOURCES} rf_pcs1.cpp fct_mpi.h fct_mpi.cpp)
ELSE(OGS_FEM_PETSC OR OGS_FEM_PETSC_GEMS)
  SET( SOURCES ${SOURCES} par_ddc.h par_ddc.cpp)
## New sparse class. WW
IF(OGS_FEM_SP OR OGS_FEM_MPI OR OGS_FEM_GEMS OR OGS_FEM_BRNS OR OGS_FEM_PQC OR OGS_FEM_MKL OR OGS_FEM_LIS OR OGS_FEM_JFNK )
	SET( SOURCES ${SOURCES} equation_class.h equation_class.cpp )
    IF(OGS_FEM_JFNK)
       SET( SOURCES ${SOURCES} rf_pcs1.cpp ) 
    ENDIF(OGS_FEM_JFNK)	
ELSE(OGS_FEM_SP OR OGS_FEM_MPI OR OGS_FEM_GEMS OR OGS_FEM_BRNS OR OGS_FEM_PQC OR OGS_FEM_MKL OR OGS_FEM_LIS OR OGS_FEM_JFNK)
	SET( SOURCES ${SOURCES} solver.h solver.cpp matrix_routines.h matrix_routines.cpp)
ENDIF(OGS_FEM_SP OR OGS_FEM_MPI OR OGS_FEM_GEMS OR OGS_FEM_BRNS OR OGS_FEM_PQC OR OGS_FEM_MKL OR OGS_FEM_LIS OR OGS_FEM_JFNK)
ENDIF()

IF(OGS_FEM_BRNS)
	SET( SOURCES ${SOURCES} rf_REACT_BRNS.h rf_REACT_BRNS.cpp ) 
ENDIF(OGS_FEM_BRNS)

IF(OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS)
	SET( SOURCES ${SOURCES} rf_REACT_GEM.h rf_REACT_GEM.cpp ) 
ENDIF(OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS)


include_directories(
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/GEM
	${CMAKE_SOURCE_DIR}/OGSProject
	${CMAKE_SOURCE_DIR}/GCC
)


IF((NOT OGS_FEM_GEMS) AND (NOT OGS_FEM_PETSC_GEMS))
        include_directories(${CMAKE_SOURCE_DIR}/FEM/Eigen)
ENDIF((NOT OGS_FEM_GEMS) AND (NOT OGS_FEM_PETSC_GEMS))

IF (LIS)
  INCLUDE_DIRECTORIES(${LIS_INCLUDE_DIR})
ENDIF (LIS)

IF (MKL)
  INCLUDE_DIRECTORIES(${MKL_INCLUDES})
ENDIF (MKL)

IF(OGS_FEM_PQC)
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/PQC)
ENDIF (OGS_FEM_PQC)

IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

# Create the library
add_library( FEM STATIC ${SOURCES} )

target_link_libraries( 
	FEM
	Base
	FileIO
 	GEO
 	MSH
	OGSProject
	GCC
)

IF(OGS_FEM_PQC)
   target_link_libraries(FEM PQC)
ENDIF(OGS_FEM_PQC)

if (WIN32)
	if (OGS_FEM_GEMS)
		target_link_libraries( FEM GEM ) 
	endif (OGS_FEM_GEMS)
	if (LIS)
		target_link_libraries( FEM ${LIS_LIBRARIES} ) 
		set_target_properties( FEM PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libcmt" )
	endif (LIS)
	if (MKL)
		target_link_libraries( FEM ${MKL_LIBRARIES} )
	endif (MKL)
	if (OGS_FEM_BRNS)
		target_link_libraries( FEM BrnsDll )
	endif (OGS_FEM_BRNS)
	if (OGS_FEM_CHEMAPP)
		target_link_libraries( FEM ca_vc_e )
	endif (OGS_FEM_CHEMAPP)
	if (PARALLEL_USE_MPI)
		SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} )
		target_link_libraries( FEM ${MPI_LIBRARIES} )
	endif (PARALLEL_USE_MPI)
endif (WIN32)

if (UNIX)
	if (OGS_FEM_BRNS)
		target_link_libraries( FEM dl ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
	endif (OGS_FEM_BRNS)

	if (LIS)
		target_link_libraries( FEM ${LIS_LIBRARIES} )
	endif (LIS)

	if (MKL)
		target_link_libraries( FEM ${MKL_LIBRARIES} )
	endif (MKL)
	if (PARALLEL_USE_MPI)
#		SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} )
		target_link_libraries( FEM ${MPI_LIBRARIES} )
	endif (PARALLEL_USE_MPI)
endif (UNIX)

IF (OGS_FEM_CHEMAPP)
	INCLUDE_DIRECTORIES ( ${CMAKE_SOURCE_DIR}/EQL )
	TARGET_LINK_LIBRARIES ( FEM EQL )
ENDIF (OGS_FEM_CHEMAPP)


IF (OGS_FEM_CAP)
	INCLUDE_DIRECTORIES (  ${CMAKE_SOURCE_DIR}/CMP )
	TARGET_LINK_LIBRARIES ( 	
	 FEM
 	 GEO
 	 MSH
 	 CMP 
 	)
    IF (UNIX)
        TARGET_LINK_LIBRARIES ( 	
          FEM
          ChemAppCS_64
          EChemAppS_64
      )
    ENDIF (UNIX)
ENDIF (OGS_FEM_CAP)


