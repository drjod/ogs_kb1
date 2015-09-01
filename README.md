## Kiel Branch One of OpenGeoSys http://www.opengeosys.org/




### Cheatsheet 

for OGS on RZ cluster Kiel

#### PBS
qstat
qstat –u $LOGIN
qstat –q angus
qdel –f $JOBID
mpi: #PBS -l select=1:ncpus=2:mpiprocs=2:mem=1gb:place=scatter
omp: #PBS -l select=1:ncpus=1:ompthreads=4:mem=1gb:place=group=host

#### Compilation

Script compileInKiel.sh for RZ (and NEC) cluster Kiel in repository tUNIX

##### Compiler
icc, icpc: /cluster/Software/intel14/composer_xe_2013_sp1/bin
     /cluster/Software/intel14/composer_xe_2013_sp1/bin/compilervars.sh intel64
mpiicc, mpiicpc: /cluster/Software/intel1502/impi/5.0.3.048/intel64/bin
     /cluster/Software/intel1502/impi/5.0.3.048/intel64/bin/mpivars.sh

Composer: /cluster/Software/intel14/composer_xe_2013_sp1.0.080
export PATH=$PATH:/cluster/Software/intel14/composer_xe_2013_sp1.0.080/compiler/lib/intel64

#### Libs

##### MKL
/cluster/Software/intel14/composer_xe_2013_sp1.0.080/mkl
/cluster/Software/intel14/composer_xe_2013_sp1.0.080/mkl/bin/mklvars.sh intel64

#####PETSC 
/work_j/SoftwareSL/Dpetsc/Dintel14/petsc-3.3-p4
module load petsc-3.3-p4-intel	

#### Mesh partitioning 

Script partition.sh for domain decomposition with METIS in repository icbc/remote/easyPeasy (calls partmesh from  ufz/mesh_partition)

For OGE_FEM_MPI: ./partition.sh numberOfPartitions –e –asci path (or no path)
For OGS_FEM_PETSC: ./partition.sh numberOfPartitions –n –bin path (or no path)

numberOfPartitions = 2, 3, 4,...

#### Numerics

##### Linear solver options

1 Gauss
2 BiCGStab OGS_FEM_MPI
3 BiCG OGS_FEM_MPI
4 QMRGG Stab
5 CG OGS_FEM_MPI
6 CGNR
7 CGS OGS_FEM_MPI
8 Richardson
9 JOR
10 SOR
11 AMG1R5
12 UMF
13 GMRes
805 Pardiso MKL

##### Preconditioners

0 No
1 Jacobi (for OGS_FEM_MPI)
100 ILU

##### Coupled Simulation

Example *.num file:
```
$OVERALL_COUPLING
2 10    ; minNumberOfIterations maxNumberOfIterations 
#NUMERICS
$PCS_TYPE
LIQUID_FLOW
$ELE_MASS_LUMPING
1
$LINEAR_SOLVER
; method error_tolerance max_iterations theta precond storage
2 6 1.e-014 5000 1.0 1 2
$NON_LINEAR_ITERATIONS
;type -- error_method -- max_iterations -- relaxation -- tolerance(s)
PICARD LMAX 2 0.0 1.
$ELE_GAUSS_POINTS
2
$COUPLING_CONTROL
LMAX 1.  ; Error tolerance in coupling loop

#NUMERICS
$PCS_TYPE
MASS_TRANSPORT
$ELE_GAUSS_POINTS
2
$LINEAR_SOLVER
2 6 1.e-014 5000 1.0 1 2
$NON_LINEAR_ITERATIONS
;type -- error_method -- max_iterations -- relaxation -- tolerance(s)
PICARD LMAX .1 0.0 1.0e-7
$COUPLING_CONTROL
LMAX 1e-1  ; Error tolerance in coupling loop

#STOP
```

##### PETSC

Solver bcgs, gmres, etc.
Preconditioner jacobi, bjacobi, sor, asm, mg
```
$LINEAR_SOLVER
petsc solver preconditioner errorTolerance maxInterations Theta
```

##### Convergence

0 r<e
1 r/b<e
2 r_n1/r_n<e
3 if r_n1>1 then r_n1/r_n < e else r<e
4 r/x<e
5 r_n1/max(x,b,r_n)<e
6 ???

##### Matrix storage

1) vollbesetzte Matrix
( param1 = Dimension )
(2) nur A[i,j]!=0 werden gespeichert (Sparse)
( param1 = Dimension )
(3) symmetrische sparse Matrix
fuer Preconditioner "incomplete Cholesky"
nur Aik !=0 mit k>=i werden gespeichert
( param1 = Dimension )
(4) unsymmetrische sparse Matrix
fuer Preconditioner "incomplete LDU-Zerlegung"
nur Aik !=0 werden gespeichert
( param1 = Dimension )

##### Flux corrected transport (FCT)

In *num:
```
$FEM_FCT
method prelimiterType const_alpha
```
method 1: linearized FCT
prelimiterType 0: just cancel, 1: minmod, 2: superbee
fct_const_apha -1: off (is default), [0.0, 1.0] 0: upwind, 1: galerkin

#### Adaptive time stepping

Example *tim file:
```
#TIME_STEPPING
$PCS_TYPE
LIQUID_FLOW
$TIME_END
1e10
$TIME_CONTROL
SELF_ADAPTIVE
4 1.3 ; multiply time step sitze by 1.3 if number of iterations <= 4 (if condition satisfied, this is selected, the following ignored)
7 1.0 ; multiply time step sitze by 1. if number of iterations <= 7
10 0.7 ; else multiply time step sitze by 0.7 ( the 10 in last line does nothing, can give more than three lines)
MAX_TIME_STEP
1e10
MIN_TIME_STEP
1e3
INITIAL_STEP_SIZE
1e3
ITERATIVE_TYPE
COUPLED ; for outer coupling loop (take keyword LINEAR for solver, NONLINEAR for Picard)
#STOP
```

#### Large pore volume cells (LPVC)

See [wiki page](https://github.com/drjod/ogs_kb1/wiki/LPVC)
```
#MEDIUM PROPERTIES
…
$ELEMENT_VOLUME_MULTIPLYER
1 8 ; 2D case x, y (z) - 1 parameter if 1D - 3 parameter if 3D
```

#### BC / IC GRADIENT

```
$DIS_TYPE
GRADIENT ref_depth ref_depth_value ref_depth_gradient
```

##### Output

##### Fluxes

Darcy flux
1st phase at nodes: VELOCITY_X1 VELOCITY_Y1 VELOCITY_Z1
1st phase at gauss point nr 0: VELOCITY1_X VELOCITY1_Y VELOCITY1_Z

To get fick flux,add 
```
$PCS_TYPE
  MASS_TRANSPORT
```  
and to get fourier flux, add  
```
$PCS_TYPE
  HEAT_TRANSPORT
```  
to output instance. 

##### Change to initial state

In *ic: set flag with keyword 
```
$STORE_VALUES 
```
in *.out: Use prefix DELTA_, e.g.
```
$DAT_TYPE
DELTA_PRESSURE1
DELTA_TEMPERATURE1 
DELTA_CONCENTRATION1
```

##### Mass balancing toolkit

See [wiki page](https://github.com/drjod/ogs_kb1/wiki/Balancing-toolkit)

Example:
```
#OUTPUT
$PCS_TYPE
LIQUID_FLOW
$DAT_TYPE
CONTENT -1 ; or mmp index
$TIM_TYPE
STEPS 1
#OUTPUT
$PCS_TYPE
LIQUID_FLOW
$DAT_TYPE
TOTAL_FLUX
$TIM_TYPE
STEPS 1
```

##### Nom-neighbor node connections

See [wiki page] (https://github.com/drjod/ogs_kb1/wiki/NNNC)

Example:
```
#SOURCE_TERM
$PCS_TYPE
HEAT_TRANSPORT
$PRIMARY_VARIABLE
TEMPERATURE1
$GEO_TYPE
SURFACE SURF_BOTTOM_LEFT
$DIS_TYPE
CONSTANT_NEUMANN 1.0
$CONNECTED_GEOMETRY
SURFACE SURF_BOTTOM_RIGHT
$CONNECT_PARAMETERS
1e6 1 ; exchange coefficient verbosity level = 0, 1
$CONNECT_MODE
2 226 0 0 1 1.e-10
; mode ; 0: symmetric, 1: non-symmetric (downwind fixed), 2 variable (dependent on velocity in reference element)
; 2 ref_element_number n_ref_x, n_ref_y, n_ref_z minimum_velocity_abs
```




  
  
