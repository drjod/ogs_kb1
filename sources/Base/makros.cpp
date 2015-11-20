
#include "makros.h"

std::string FileName;
std::string FilePath;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || \
        defined(USE_MPI_GEMS) || defined(USE_MPI_BRNS) || defined(USE_MPI_KRC) || defined(USE_PETSC)
int mysize;
int myrank;
#endif
