# Source files
set( HEADERS
	Vector3.h
	Matrix.h
	MathTools.h
	AnalyticalGeometry.h
	max.h
	LinkedTriangle.h
	EarClippingTriangulation.h
	# directory InterpolationAlgorithms
	InterpolationAlgorithms/CubicSpline.h
	InterpolationAlgorithms/InverseDistanceInterpolation.h
	InterpolationAlgorithms/LinearIntervalInterpolation.h
	InterpolationAlgorithms/PiecewiseLinearInterpolation.h
	# directory LinAlg
	LinAlg/VectorNorms.h
	LinAlg/LinearSolver.h
	LinAlg/DirectLinearSolver.h
	LinAlg/DenseDirectLinearSolver.h
	LinAlg/GaussAlgorithm.h
#	LinAlg/HouseholderAlgorithm.h
	LinAlg/TriangularSolve.h
	LinAlg/IterativeLinearSolver.h
)

set( SOURCES
	MathTools.cpp
	AnalyticalGeometry.cpp
	LinkedTriangle.cpp
	EarClippingTriangulation.cpp
	# directory InterpolationAlgorithms
	InterpolationAlgorithms/CubicSpline.cpp
	InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp
	# directory LinAlg
#	LinAlg/GaussAlgorithm.cpp
	LinAlg/TriangularSolve.cpp
#	LinAlg/HouseholderAlgorithm.cpp
)

IF(OGS_FEM_PETSC OR OGS_FEM_PETSC_GEMS)
   set( HEADERS  ${HEADERS} PETSC/PETScLinearSolver.h)
   set( SOURCES  ${SOURCES} PETSC/PETScLinearSolver.cpp)

ENDIF()
# Create the library
add_library( MathLib STATIC ${HEADERS} ${SOURCES} )

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MathLib/LinAlg
	${CMAKE_SOURCE_DIR}/MathLib/PETSC
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
)

target_link_libraries (
	MathLib
	Base
)

