## Kiel Branch One of OpenGeoSys http://www.opengeosys.org/




### Cheat sheet 

for OGS on RZ cluster Kiel

#### PBS
qstat <br>
qstat –u $LOGIN <br>
qstat –q angus <br>
qdel –f $JOBID <br>
mpi: #PBS -l select=1:ncpus=2:mpiprocs=2:mem=1gb:place=scatter <br>
omp: #PBS -l select=1:ncpus=1:ompthreads=4:mem=1gb:place=group=host <br>

An example bash script for an MPI-Job with 8 cores on the angus queue (rz cluster):
```
#!/bin/bash
#PBS -o screenout.txt
#PBS -j oe
#PBS -r n
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=8:mem=3gb
#PBS -l place=scatter
#PBS -q angus
#PBS -N test

cd $PBS_O_WORKDIR

. /usr/share/Modules/init/bash

. /cluster/Software/intel1502/composer_xe_2015.2.164/bin/compilervars.sh  intel64
. /cluster/Software/intel1502/composer_xe_2015.2.164/mkl/bin/intel64/mklvars_intel64.sh
. /cluster/Software/intel1502/impi/5.0.3.048/intel64/bin/mpivars.sh

time mpirun -r rsh -machinefile $PBS_NODEFILE -n 8 ogs_OGS_FEM_MPI testCase

qstat -f $PBS_JOBID
exit
```
In this example, name of ogs is 'ogs_OGS_FEM_MPI' and of input files 'testCase'. Put this script with ogs and input files into one folder, cd there and execute the script. Output will be written into screenout.txt. <br>
To change the number of cores, you modify the number (8) in two locations of the script. These are in the PBS script command ncpus=8 and the mpirun command $PBS_NODEFILE -n 8. Important: Each node hosts 16 cores. You can select more than 16 cores by taking more cpus. For instance, you get with the PBS- command   
```
#PBS -l select=3:ncpus=8:mem=64gb
```
24 cores (set also $PBS_NODEFILE -n 24). <br>
The example script sets maximum wall time as 2 hours,  3G of memory are available for the simulation (each cluster node has 64 (128)GB available in total). 


#### Compilation

Script compileInKiel.sh for RZ cluster, NEC Cluster, Lokstedt server Kiel in repository tUNIX <br>

#### Mesh partitioning 

Script partition.sh for domain decomposition with METIS in repository icbc/remote/easyPeasy (calls partmesh from  ufz/mesh_partition) <br>

For OGE_FEM_MPI: ./partition.sh numberOfPartitions –e –asci path (or no path) <br>
For OGS_FEM_PETSC: ./partition.sh numberOfPartitions –n –bin path (or no path) <br>

numberOfPartitions = 2, 3, 4,... <br>

#### Numerics

##### Linear solver options

With OGS_FEM, OGS_FEM_SP: <br>
1 Gauss <br>
2 BiCGStab (OGS_FEM_MPI) <br>
3 BiCG (OGS_FEM_MPI) <br>
4 QMRGG Stab <br>
5 CG (OGS_FEM_MPI) <br>
6 CGNR <br>
7 CGS (OGS_FEM_MPI) <br>
8 Richardson <br>
9 JOR <br>
10 SOR <br>
11 AMG1R5 <br>
12 UMF <br>
13 GMRes <br>

With OGS_FEM_MKL
805 Pardiso <br>

##### Preconditioners

0 No <br>
1 Jacobi (for OGS_FEM_MPI) <br>
100 ILU <br>

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

Solver bcgs, gmres, etc. <br>
Preconditioner jacobi, bjacobi, sor, asm, mg
```
$LINEAR_SOLVER
petsc solver preconditioner errorTolerance maxInterations Theta
```

##### Convergence
```
0 r<e 
1 r/b<e
2 r_n1/r_n<e
3 if r_n1>1 then r_n1/r_n < e else r<e
4 r/x<e
5 r_n1/max(x,b,r_n)<e 
6 ??? 
```
##### Matrix storage

1) vollbesetzte Matrix <br>
( param1 = Dimension ) <br>
(2) nur A[i,j]!=0 werden gespeichert (Sparse) <br>
( param1 = Dimension ) <br>
(3) symmetrische sparse Matrix <br>
fuer Preconditioner "incomplete Cholesky" <br>
nur Aik !=0 mit k>=i werden gespeichert <br>
( param1 = Dimension ) <br>
(4) unsymmetrische sparse Matrix <br>
fuer Preconditioner "incomplete LDU-Zerlegung" <br>
nur Aik !=0 werden gespeichert <br>
( param1 = Dimension ) <br>

##### Flux corrected transport (FCT)

In *num:
```
$FEM_FCT
method prelimiterType const_alpha
```
method 1: linearized FCT <br>
prelimiterType 0: just cancel, 1: minmod, 2: superbee <br>
fct_const_apha -1: off (is default), [0.0, 1.0] 0: upwind, 1: galerkin <br>

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

See [wiki page](https://github.com/drjod/ogs_kb1/wiki/LPVC) <br>
Example:
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

Darcy flux for LIQUID_FLOW, GROUNDWATER_FLOW, RICHARDS_FLOW<br>
at nodes: VELOCITY_X1 VELOCITY_Y1 VELOCITY_Z1 <br>
at gauss point nr 0: VELOCITY1_X VELOCITY1_Y VELOCITY1_Z <br>

To get fick flux, add 
```
$PCS_TYPE
  MASS_TRANSPORT
```  
and to get fourier flux, add  
```
$PCS_TYPE
  HEAT_TRANSPORT
```  
to output instance.  <br>

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

See [wiki page](https://github.com/drjod/ogs_kb1/wiki/Balancing-toolkit) <br>

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

See [wiki page] (https://github.com/drjod/ogs_kb1/wiki/NNNC) <br>

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




  
  
