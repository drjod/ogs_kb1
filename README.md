## Kiel Branch One of OpenGeoSys http://www.opengeosys.org/

### Cheat sheet 

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

#### Reload results

You can store primary variable values with the DAT_TYPE PRIMARY_VARIABLES and reload them as initial conditions with the DIS_TYPE DIRECT into your next simulation 

Example:

- Generate output file with pressure from LIQUID_FLOW at time 1000:
```
#OUTPUT
$PCS_TYPE
LIQUID_FLOW
$NOD_VALUES
PRESSURE1
$GEO_TYPE
DOMAIN
$DAT_TYPE
PRIMARY_VARIABLES
$TIM_TYPE
1000
```
- Read the file: 
```
#INITIAL_CONDITION
$PCS_TYPE
LIQUID_FLOW
$PRIMARY_VARIABLE
PRESSURE1
$DIS_TYPE
DIRECT LIQUID_FLOW_domain_primary_variables.txt
```
  
#### Large pore volume cells (LPVC)

See [wiki page](https://github.com/drjod/ogs_kb1/wiki/LPVC) <br>
Example:
```
#MEDIUM PROPERTIES
…
$ELEMENT_VOLUME_MULTIPLYER
1 1 8 ; 2D case x, y (z) - first number is mode, then 1 parameter if 1D, 2 parameter if 2D, 3 parameter if 3D
```
modes <br> 
0: increase storage only <br>
1: increase storage and adapt Laplace, advection, velocity


#### BC / IC GRADIENT

```
$DIS_TYPE
GRADIENT ref_depth ref_depth_value ref_depth_gradient
```

#### Output

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

##### Averaging

Implemented for surfaces, e.g. to get an average surface temperature
```
$NOD_VALUES
TEMPERATURE1
$GEO_TYPE
 SURFACE sfc_1
$DIS_TYPE
 AVERAGE
```

#### Mass balancing toolkit

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

#### Cutting and plumbing

##### Copy Gauss point velocities
Copy GP velocities from one element to a number of other elements:
```
#PROCESS
 $PCS_TYPE
  LIQUID_FLOW
 $NUM_TYPE
  NEW
 $COPY_GP_VELOCITIES  ; Give total number of indices, then index to copy from, the indices to copy to.
   2  44  46
  ```

##### Non-neighbor node connections (NNNC)

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
<<<<<<< HEAD

###### Process coupling 

Flux for RICHARDS_FLOW to couple to OVERLAND_FLOW:
```
#SOURCE_TERM
 $PCS_TYPE
  RICHARDS_FLOW
 $PRIMARY_VARIABLE
  PRESSURE1
 $GEO_TYPE
  POLYLINE INTERFACE
 $DIS_TYPE_CONDITION
  CONSTANT_NEUMANN 1
  PCS OVERLAND_FLOW
  HEAD
  2.83333333e-2 1.e-3 0.0 0.0 ; leakage, rill height, given value, residual permeability 
```
Couple HEAT_TRANSPORT to a given temperature of 10:
```
#SOURCE_TERM
 $PCS_TYPE
  HEAT_TRANSPORT
 $PRIMARY_VARIABLE
  TEMPERATURE1
 $GEO_TYPE
  POLYLINE INTERFACE
 $DIS_TYPE_CONDITION
  CONSTANT_NEUMANN 1
  HEAT_TRANSPORT
  GIVEN_TEMPERATURE
  0.1 0.0 10.0 0.0 ; leakage, rill height, given value, residual permeability 
```


###### Deactivated subdomain
=======
##### Deactivated subdomain
>>>>>>> develop

```
#PROCESS
 $PCS_TYPE
  LIQUID_FLOW
 $DEACTIVATED_SUBDOMAIN
  2  ; number of deactivated subdomains
  1 2   ; patch indices (mmp groups)
```
##### Material density
Different fluid density models for different patch indices (mmp groups). Supported are density models 1-8 and 14, e.g.:
```
#FLUID_PROPERTIES
 $FLUID_TYPE
  LIQUID
 $MATERIAL_DENSITY
  0 1 1000    ; density model 1 for patch index 0
  1 8 0       ; density model 8 for patch index 1
  ```
##### Material dependent fluid
Different fluid properties for different patch indices (mmp groups).

Example:

Generate two liquids

```
#FLUID_PROPERTIES
 $FLUID_TYPE
   LIQUID_WELL
   
#FLUID_PROPERTIES
 $FLUID_TYPE
   LIQUID_AQ  
```

and select liquids as:

```
#MEDIUM_PROPERTIES
 $DEPENDENT_FLUID
   WELL
   
#MEDIUM_PROPERTIES
 $DEPENDENT_FLUID
   AQ
  ```


  
#### Fluid properties

##### Density models

0: Curve rho(x) <br>
1: rho = const <br>
2: rho(p) = rho_0*(1+beta_p*(p-p_0)) <br>
3: rho(C) = rho_0*(1+beta_C*(C-C_0)) <br>
4: rho(T) = rho_0*(1+beta_T*(T-T_0)) <br>
5: rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0)) <br>
6: rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0)) <br>
7: Pefect gas <br>
8: Density output AB-model <br>
10: density from fct-file with temperature-pressure values <br>
11: Redlich-Kwong EOS for different fluids <br>
12: Peng-Robinson EOS for different fluids <br>
13: Helmholtz free Energy <br>
14: Exponential law <br>
15: Amagat's law for mixture <br>
18: Density at nodes from the phase transition model <br>
19: Density from GEMS <br>
20: rho(p,T, C) for water, range p < 100 MPa, 0 <= T <= 350 °C  <br>
21: rho(p,C,T) = rho_0*(1+beta_p*(p-p_0)+beta_C*(C-C_0)+beta_T*(T-T_0)) <br>
23: Density depends on concentration of salt and dissolved CO2 
26: Dalton's law + ideal gas for use with TNEQ/TES <br><br>

Example specifications in input file *mfp:<br>
1. Model 0 - rho = rho(x): <br>
```
$DENSITY 
 0 1 
```
Parameter "1" refers to the first data table in a *.rfd ascii input file, which you will have to provide in this case:<br>
```
;Curve 1  Temp Density
#CURVES
 273.15  999.8675792
 274.15  999.9265054
 ...
#STOP
```
<br>
2.  Model 4 - rho = rho(T): <br>
```
$DENSITY
 4 1.000000e+003	0	0.2	; C0 is 0. beta_C (drho_dC) is 0.2
```
<br>

More information: <br>
CFluidProperties::Read(std::ifstream*)<br>
CFluidProperties::Density(double*) [separate calculation for solving PDEs and density output]<br>

###### Viscosity models
0: Curve my(x) <br>
1: my = const <br>
2: my(p) = my_0*(1+gamma_p*(p-p_0)) <br>
3: my^l(T) acoording toYaws et al. (1976) <br>
4: my^g(T) acoording to Marsily (1986) <br>
5: my^g(p,T) acoording toReichenberg (1971) <br>
6: my(C,T) <br>
8: my(p,C,T) <br>
9: my(rho,T) <br>
15: VTPR-EoS <br>
18: ViscositY at nodes from the phase transition model <br>
19: Viscosity for GEMS  <br>
26: Wilke (see Poling, B. E.; Prausnitz, J. M.; John Paul, O. & Reid, R. C. The properties of gases and liquids McGraw-Hill New York, 2001, 5: page 9.21) <br> <br>

More information: <br>
CFluidProperties::Read(std::ifstream*)<br>
CFluidProperties::Viscosity(double* variables)

#### To run OGS on RZ cluster Kiel

##### PBS
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
To change the number of cores, you modify the number (8) in two locations of the script. These are in the PBS script command ncpus=8 and the mpirun command $PBS_NODEFILE -n 8. Important: Each node hosts 16 cores. You can select more than 16 cores by taking more nodes. For instance, you get with the PBS- command   
```
#PBS -l select=3:ncpus=8:mem=64gb
```
24 cores (set also $PBS_NODEFILE -n 24). <br>
The example script sets maximum wall time as 2 hours,  3G of memory are available for the simulation (each cluster node has 64 (128)GB available in total). 


#### Compilation

You find the newest script compileInKiel.sh for RZ cluster, NEC Cluster, Lokstedt server Kiel in repository tUNIX. Instructions are in the script.  <br>

#### Mesh partitioning 

Script partition.sh for domain decomposition with METIS in repository icbc/remote/easyPeasy (calls partmesh from  ufz/mesh_partition) <br>

For OGE_FEM_MPI: ./partition.sh numberOfPartitions –e –asci path (or no path) <br>
For OGS_FEM_PETSC: ./partition.sh numberOfPartitions –n –bin path (or no path) <br>

numberOfPartitions = 2, 3, 4,... <br>
  
  
