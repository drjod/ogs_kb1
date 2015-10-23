
if(OGS_FEM_IPQC)
	add_subdirectory(IPQC)
endif()

if(OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS)
	add_subdirectory(GEM)
endif()

if(OGS_FEM_PQC)
	add_subdirectory ( PQC )
endif()
