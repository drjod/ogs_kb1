set( HEADERS
	AnalyticalGeometry.h
	EarClippingTriangulation.h
	InterpolationAlgorithms/CubicSpline.h
	InterpolationAlgorithms/InverseDistanceInterpolation.h
	InterpolationAlgorithms/LinearIntervalInterpolation.h
	InterpolationAlgorithms/PiecewiseLinearInterpolation.h
	LinAlg/DenseDirectLinearSolver.h
	LinAlg/DirectLinearSolver.h
	LinAlg/GaussAlgorithm.h
	LinAlg/IterativeLinearSolver.h
	LinAlg/LinearSolver.h
	LinAlg/TriangularSolve.h
	LinAlg/VectorNorms.h
	LinkedTriangle.h
	MathTools.h
	Matrix.h
	max.h
	Vector3.h
)

set( SOURCES
	AnalyticalGeometry.cpp
	EarClippingTriangulation.cpp
	InterpolationAlgorithms/CubicSpline.cpp
	InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp
	LinAlg/TriangularSolve.cpp
	LinkedTriangle.cpp
	MathTools.cpp
)

if(OGS_FEM_PETSC OR OGS_FEM_PETSC_GEMS)
	set( SOURCES ${SOURCES} PETSC/PETScLinearSolver.h PETSC/PETScLinearSolver.cpp)
endif()

add_library( MathLib STATIC ${HEADERS} ${SOURCES} )

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MathLib/LinAlg
	${CMAKE_SOURCE_DIR}/MathLib/PETSC
# ${CMAKE_SOURCE_DIR}/../Libs/WellDoubletControl

)

target_link_libraries (
	MathLib
	Base
)

